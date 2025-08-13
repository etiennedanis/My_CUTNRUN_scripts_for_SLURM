#!/usr/bin/env bash

RAN_ON=$(date +%Y_%m_%d_%H_%M_%S)
NUMBER=0

INVESTIGATOR=
PROJECT=

FASTQSCREEN=/references/genomes/FastQ_Screen_Genomes
SCRATCH=/scratch/myaccount
CUTNRUN=$SCRATCH/Documents/cutnrun_data
DATA=$CUTNRUN/${INVESTIGATOR}_data/$PROJECT
SCRIPTS=$DATA/scripts
TRIMMED=$DATA/trimmed_fastq_files
SCREENO=$DATA/fastq_screen_files_w_mycoplasma

mkdir -p $SCREENO
cd $TRIMMED

#for FILE in `ls *.fastq.gz | cut -d "." -f 1`; # for single-end sequencing reads
for FILE in `ls *.fastq.gz | cut -d "." -f 1,2`; # for paired-end sequencing reads
do
NUMBER=$((${NUMBER}+1)) 
cat << EOF > $SCRIPTS/${RAN_ON}_${FILE}_script.sh
#!/usr/bin/bash

#SBATCH -J SCREEN_$NUMBER
#SBATCH -e $SCRIPTS/${RAN_ON}_${FILE}_FastQ_Screen_w_mycoplasma.err
#SBATCH -o $SCRIPTS/${RAN_ON}_${FILE}_FastQ_Screen_w_mycoplasma.out

#SBATCH --nodes=1                            # Number of nodes
#SBATCH --mem=10G                            # RAM per node
#SBATCH --cpus-per-task=8                   # CPUs per node
#SBATCH --time=1:00:00

module load bowtie2/2.5.0
module load perl
module load anaconda

FASTQSCREEN=$FASTQSCREEN
TRIMMED=$TRIMMED
SCREENO=$SCREENO

conda activate for_Fastq_Screen

/projects/software/anaconda/envs/for_Fastq_Screen/bin/fastq_screen --conf ${FASTQSCREEN}/w_mycoplasma_fastq_screen.conf \
    --threads 6 \
    --nohits \
    --outdir $SCREENO \
    --top 200000 \
    --aligner bowtie2 \
    $TRIMMED/${FILE}.fastq.gz

conda deactivate

EOF

sbatch $SCRIPTS/${RAN_ON}_${FILE}_script.sh

rm -f $SCRIPTS/${RAN_ON}_${FILE}_script.sh

done

# --threads: To specify the number of threads                                                                                                                      
# --nohits: If you want the sequences that do not map to anything to be written into a new file                                                                   
# --outdir: To specify the output directory                                                                                                                       
# --top 200000: To ony evaluate the top 200,000 sequences                  

# --conf: For the location of the configuration, which includes the location of the genomes and of bowtie2     



