#!/usr/bin/env bash

#SBATCH -J Nested
#SBATCH -e Launching_Alignment_jobs.err
#SBATCH -o Launching_Alignment_jobs.out

#SBATCH --nodes=1                            # Number of nodes
#SBATCH --mem=1G                             # RAM per node
#SBATCH --cpus-per-task=1                    # CPUs per node
#SBATCH --time=0-00:00:30

RAN_ON=$(date +%Y_%m_%d_%Hhr%Mmin)
NUMBER=0

# Directories where all the data sets are located
INVESTIGATOR=
PROJECT=
PROCESS=
PROCESS2=

ADAPTERS=/bbtools/bbmap/resources
LIBRARY=/genomes
SCRATCH=/scratch/myaccount
CUTNRUN=$SCRATCH/Documents/cutnrun_data
DATA=$CUTNRUN/${INVESTIGATOR}_data/$PROJECT
SCRIPTS=$DATA/scripts
FASTQ=$DATA/fastq_files
TRIMMED=$DATA/trimmed_fastq_files
SORTED=$DATA/bam_files_sorted
INDEX=$LIBRARY/Bowtie2_indices/human_GRCh38_v104_genome
NOMT=$DATA/bam_files_wo_mtDNA
PICARD=/picard/picard/build/libs/picard.jar
MYTMP=$SCRATCH/tmp
NODUPS=$DATA/bam_files_wo_duplicates
EXCLUSION=$LIBRARY/exclusion_lists
EXCLUDED=$DATA/bam_files_w_keeping_duplicates_after_using_exclusion_list
FLAGSTAT=$DATA/flagstat_reports
BIGWIGS1=$DATA/bigwig_files_normalized_w_RPGC
BIGWIGS2=$DATA/bigwig_files_not_normalized

export TMP=/scratch/$USER
export TEMP=/scratch/$USER
export TMPDIR=/scratch/$USER
export TEMPDIR=/scratch/$USER
export TMP_DIR=/scratch/$USER
export MYTMP=/scratch/$USER

mkdir -p $TRIMMED
mkdir -p $SORTED
mkdir -p $NOMT
mkdir -p $NODUPS
mkdir -p $MYTMP
mkdir -p $EXCLUDED
mkdir -p $FLAGSTAT
mkdir -p $BIGWIGS1
mkdir -p $BIGWIGS2

cd $FASTQ
for FILE in `ls *R1.fastq.gz | cut -d  "." -f 1`;
do
NUMBER=$((${NUMBER}+1))
cat << EOF > $SCRIPTS/${RAN_ON}_${FILE}.sh
#!/usr/bin/bash

#SBATCH -J CUTN_$NUMBER
#SBATCH -e $SCRIPTS/${RAN_ON}_${FILE}_Trimming_alignment_filtering_and_bigwig_creation.err
#SBATCH -o $SCRIPTS/${RAN_ON}_${FILE}_Trimming_alignment_filtering_and_bigwig_creation.out

#SBATCH --nodes=1                            # Number of nodes
#SBATCH --mem=40G                            # RAM per node
#SBATCH --cpus-per-task=12                   # CPUs per node
#SBATCH --time=0-05:00:00

module load bbtools/39.01
module load bowtie2/2.5.0
module load samtools/1.16.1
module load bedtools/2.29.1
module load jdk/18.0.1.1
module load picard/2.27.5
module load python
module load anaconda

export TMP=/scratch/$USER
export TEMP=/scratch/$USER
export TMPDIR=/scratch/$USER
export TEMPDIR=/scratch/$USER
export TMP_DIR=/scratch/$USER
export MYTMP=/scratch/$USER

ADAPTERS=$ADAPTERS
DATA=$DATA
CONCAT=$CONCAT
TRIMMED=$TRIMMED
SORTED=$SORTED
INDEX=$INDEX
NOMT=$NOMT
PICARD=$PICARD
MYTMP=$MYTMP
NODUPS=$NODUPS
EXCLUSION=$EXCLUSION
EXCLUDED=$EXCLUDED
FLAGSTAT=$FLAGSTAT
BIGWIGS1=$BIGWIGS1
BIGWIGS2=$BIGWIGS2

# Trim adapter sequences using bbduk from bbtools/BBMap
bbduk.sh -Xmx38g \
    threads=10 \
    in1=$FASTQ/${FILE}.R1.fastq.gz \
    in2=$FASTQ/${FILE}.R2.fastq.gz \
    out1=$TRIMMED/${FILE}.R1.fastq.gz \
    out2=$TRIMMED/${FILE}.R2.fastq.gz \
    ref=$ADAPTERS/truseq.fa.gz \
    ftm=5 \
    ktrim=r k=23 mink=11 hdist=1 \
    minlen=30 qtrim=rl trimq=15 \
    tpe tbo
   
# Align paired-end reads with bowtie2
bowtie2 -p 10 -x $INDEX/GRCh38_v104 \
    -1 $TRIMMED/${FILE}.R1.fastq.gz \
    -2 $TRIMMED/${FILE}.R2.fastq.gz | \
    samtools view -q 30 -b -@ 10 | \
    samtools sort -@ 10 -o $SORTED/${FILE}.bam

##Index the resorted files
samtools index -@ 10 $SORTED/${FILE}.bam

## Remove the reads matching the mitochondrial chromosome
samtools idxstats -@ 10 $SORTED/${FILE}.bam | cut -f 1 | grep -v "MT" | \
     xargs samtools view -b $SORTED/${FILE}.bam > $NOMT/${FILE}.bam

# "chrM" for murine genome
# "MT" for human genome

## Remove duplicates with picard
java -Xmx38g -Djava.io.tmpdir=$MYTMP -jar $PICARD MarkDuplicates \
    -I $NOMT/${FILE}.bam \
    -O $NODUPS/${FILE}.bam \
    -M $NODUPS/${FILE}_picard_metrics.txt \
    -REMOVE_DUPLICATES true \
    --TMP_DIR $MYTMP
# Remove reads in ENCODE exclusion list
bedtools intersect -v -abam $NOMT/${FILE}.bam -b $EXCLUSION/GRCh38_unified_blacklist.bed \
    > $EXCLUDED/${FILE}.bam

## Count the number of mapped reads if extraction is needed or for general information
echo "$FILE" >> $FLAGSTAT/${RAN_ON}_${FILE}_report.txt
samtools flagstat -@ 10 $EXCLUDED/${FILE}.bam >> $FLAGSTAT/${RAN_ON}_${FILE}_report.txt
echo "" >> $FLAGSTAT/${RAN_ON}_${FILE}_report.txt 

# Index bam files
samtools index -@ 10 $EXCLUDED/${FILE}.bam

## Convert .bam files into .bw (bigwig) files
conda activate for_CUTNRUN

bamCoverage --bam $EXCLUDED/${FILE}.bam \
    -o $BIGWIGS1/${FILE}.bw \
    --numberOfProcessors 10 \
    --normalizeUsing RPGC \
    --binSize 10 \
    --effectiveGenomeSize 2862010428  # For 150 bp - Number to use for human reads when multimapped are removed 

bamCoverage --bam $EXCLUDED/${FILE}.bam \
   -o $BIGWIGS2/${FILE}.bw \
   --numberOfProcessors 10 \
   --normalizeUsing None \
   --binSize 10 \
   --effectiveGenomeSize 2862010428

conda deactivate

EOF

sbatch $SCRIPTS/${RAN_ON}_${FILE}.sh

rm -f $SCRIPTS/${RAN_ON}_${FILE}.sh

done


