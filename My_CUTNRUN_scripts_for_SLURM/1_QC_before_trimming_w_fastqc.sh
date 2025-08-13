#!/usr/bin/env bash

#SBATCH -J Nested
#SBATCH -e QC_before_trimming_jobs.err
#SBATCH -o QC_before_trimming_jobs.out

#SBATCH --nodes=1                            # Number of nodes
#SBATCH --mem=1G                             # RAM per node
#SBATCH --cpus-per-task=1                    # CPUs per node
#SBATCH --time=0-00:00:30


RAN_ON=$(date +%Y_%m_%d_%H_%M_%S)
NUMBER=0

INVESTIGATOR=
PROJECT=

SCRATCH=/scratch/myaccount
CUTNRUN=$SCRATCH/Documents/cutnrun_data
DATA=$CUTNRUN/${INVESTIGATOR}_data/$PROJECT
SCRIPTS=$DATA/scripts
FASTQ=$DATA/fastq_files
FASTQC=$DATA/qc_files_before_trimming

mkdir -p $FASTQC
cd $FASTQ

for FILE in *.gz; do mv "$file" "${file//-/_}"; done

rename "_R1_001" ".R1" *.gz
rename "_R2_001" ".R2" *.gz

#for FILE in `ls *.fastq.gz | cut -d "." -f 1`; # for single-end sequencing reads
for FILE in `ls *.fastq.gz | cut -d "." -f 1,2`; # for paired-end sequencing reads
do
echo $FILE
NUMBER=$((${NUMBER}+1)) 
cat << EOF > $SCRIPTS/${RAN_ON}_${FILE}_script.sh
#!/usr/bin/bash

#SBATCH -J QC_$NUMBER
#SBATCH -e $SCRIPTS/${RAN_ON}_${FILE}_QC_before_trimming.err
#SBATCH -o $SCRIPTS/${RAN_ON}_${FILE}_QC_before_trimming.out

#SBATCH --nodes=1                            # Number of nodes
#SBATCH --mem=4G                             # RAM per node
#SBATCH --cpus-per-task=5                    # CPUs per node
#SBATCH --time=01:00:00

module load jdk/1.8
module load fastqc/0.11.9

FASTQ=$FASTQ
FASTQC=$FASTQC

cd $FASTQ

fastqc -t 4 -o $FASTQC ${FILE}.fastq.gz

EOF

sbatch $SCRIPTS/${RAN_ON}_${FILE}_script.sh

rm -f $SCRIPTS/${RAN_ON}_${FILE}_script.sh

done


