#!/usr/bin/bash

#SBATCH -J Nested
#SBATCH -e Launching_read_extraction_jobs.err
#SBATCH -o Launching_read_extraction_jobs.out

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
PROCESS2=normalized_w_human_reads
NORMFILE=Normalization_file_w_human_reads.csv

SCRATCH=/scratch/myaccount
CUTNRUN=$SCRATCH/Documents/cutnrun_data
DATA=$CUTNRUN/${INVESTIGATOR}_data/$PROJECT
SCRIPTS=$DATA/scripts
EXCLUDED=$DATA/bam_files_wo_duplicates_after_using_exclusion_list$PROCESS
NORM=$DATA/normalization_files
SORTED2=$DATA/bam_files_${PROCESS2}_and_sorted
BEDGRAPHS=$DATA/bedgraph_files_$PROCESS2
BIGWIGS=$DATA/bigwig_files_$PROCESS2

# Create the directories necessary to save all the different types of files
mkdir -p $SORTED2
mkdir -p $BEDGRAPHS
mkdir -p $BIGWIGS

# Create a for loop to create a script going through each .bam file also using a heredoc
cd $EXCLUDED

for FILE in `ls *.bam | cut -d "." -f 1`;
do
NUMBER=$((${NUMBER}+1))

cat $NORM/$NORMFILE | while read LINE;
do read SAMPLES RATIOS READS;
        if [[ $FILE = $SAMPLES ]];
        then 
        FILE=$SAMPLES
        RATIO=$RATIOS
        READ=$READS
        fi     

cat << EOF > $SCRIPTS/${RAN_ON}_${FILE}.sh
#!/usr/bin/bash

#SBATCH -J CUTN_$NUMBER
#SBATCH -e $SCRIPTS/${RAN_ON}_${FILE}_QC_before_trimming.err
#SBATCH -o $SCRIPTS/${RAN_ON}_${FILE}_QC_before_trimming.out

#SBATCH --nodes=1                            # Number of nodes
#SBATCH --mem=20G                            # RAM per node
#SBATCH --cpus-per-task=12                   # CPUs per node
#SBATCH --time=0-03:00:00

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

EXCLUDED=$EXCLUDED
NORM=$NORM
SORTED2=$SORTED2
BEDGRAPHS=$BEDGRAPHS
BIGWIGS=$BIGWIGS

# Extract the same number of reads for all samples generated with the same antibody
# Directly followed by resorting the  files
samtools view -@ 10 -hs $RATIO $EXCLUDED/${FILE}.bam | \
    samtools sort -@ 10 -o $SORTED2/extract_${READ}M_${FILE}.bam

# Index the resorted files
samtools index -@ 10 $SORTED2/extract_${READ}M_${FILE}.bam

# Create .bedgraph files
bedtools genomecov -ibam $SORTED2/extract_${READ}M_${FILE}.bam \
    -bg \
    > $BEDGRAPHS/extract_${READ}M_${FILE}.bedgraph

# Create .bigwig files
conda activate for_CUTNRUN

bamCoverage --bam $SORTED2/extract_${READ}M_${FILE}.bam \
     -o $BIGWIGS/extract_${READ}M_${FILE}.bw \
    --numberOfProcessors 10 \
    --binSize 10 \
    --normalizeUsing None
    
    #RPGC \
    #--effectiveGenomeSize 2805636331

# 2701495761 # When multimapped reads are removed and for 50 bp reads - human hg38 genome
# 2747877777 # When multimapped reads are removed and for 75 bp reads - human hg38 genome
# 2805636331 # hen multimapped reads are removed and for 100 bp reads - human hg38 genome
# 2862010578 # When multimapped reads are removed and for 150 bp reads - human hg38 genome  

# 2913022398 # When multimapped areads are kept and for 150 bp reads - human hg38 genome
##    2652783500 # This is for mm10 (UCSC) or GRCm38 (NCBI)
##   --extendReads
## More genome size here:
## https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html

conda deactivate

EOF

done

sbatch $SCRIPTS/${RAN_ON}_${FILE}.sh

rm -f $SCRIPTS/${RAN_ON}_${FILE}.sh

done

