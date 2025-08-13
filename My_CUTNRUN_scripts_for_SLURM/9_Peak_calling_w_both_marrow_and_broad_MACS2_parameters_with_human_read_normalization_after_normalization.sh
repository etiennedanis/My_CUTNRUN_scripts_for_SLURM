#!/usr/bin/env bash

#SBATCH -J Nested
#SBATCH -e Peak_calling_launch_after_normalization.err
#SBATCH -o Peak_calling_launch_after_normalization.out

#SBATCH --nodes=1                            # Number of nodes
#SBATCH --mem=1G                             # RAM per node
#SBATCH --cpus-per-task=1                    # CPUs per node
#SBATCH --time=0-00:00:30

RAN_ON=$(date +%Y_%m_%d_%Hhr%Mmin)
NUMBER=0

INVESTIGATOR=
PROJECT=

PROCESS0=
PROCESS1=normalized
PROCESS2=and_w_broad_peak_detection
PROCESS3=and_w_narrow_peak_detection
TYPE=macs

SCRATCH=/scratch/myaccount
CUTNRUN=$SCRATCH/Documents/cutnrun_data
DATA=$CUTNRUN/${INVESTIGATOR}_data/$PROJECT
SORTED2=$DATA/bam_files_normalized_w_human_reads_and_sorted
MACS2FILES=$DATA/macs2_files_${PROCESS1}_${PROCESS2}$PROCESS0
MACS2FILES2=$DATA/macs2_files_${PROCESS1}_${PROCESS3}$PROCESS0
SCRIPTS=$DATA/scripts

export TMP=/scratch/$USER
export TEMP=/scratch/$USER
export TMPDIR=/scratch/$USER
export TEMPDIR=/scratch/$USER

mkdir -p $MACS2FILES
mkdir -p $MACS2FILES2

cd $SORTED2

for FILE in `ls *.bam | cut -d "." -f 1`;
do
NUMBER=$((${NUMBER}+1))
cat << EOF > $SCRIPTS/${RAN_ON}_${FILE}.sh
#!/usr/bin/bash

#SBATCH -J MACS2_$NUMBER
#SBATCH -e $SCRIPTS/${RAN_ON}_${FILE}_Peak_calling_after_normalization.err
#SBATCH -o $SCRIPTS/${RAN_ON}_${FILE}_Peak_calling_after_normalization.out

#SBATCH --nodes=1                            # Number of nodes
#SBATCH --mem=10G                            # RAM per node
#SBATCH --cpus-per-task=8                    # CPUs per node
#SBATCH --time=0-01:00:00

module load anaconda
conda activate for_CUTNRUN

SORTED2=$SORTED2
MACS2FILES=$MACS2FILES

export TMP=/scratch/$USER
export TEMP=/scratch/$USER
export TMPDIR=/scratch/$USER
export TEMPDIR=/scratch/$USER

# Call broad peaks with macs2 not using the input as a control (no control)
macs2 callpeak -t $SORTED2/${FILE}.bam \
    -f BAMPE \
    --gsize hs --broad --broad-cutoff 0.1 \
    --name ${FILE} \
    --outdir $MACS2FILES

# Call narrow peaks with macs2 not using the input as a control (no control)
macs2 callpeak -t $SORTED2/${FILE}.bam \
    -f BAMPE \
    --gsize hs --qvalue 0.01 \
    --name ${FILE} \
    --outdir $MACS2FILES2

conda deactivate

EOF

sbatch $SCRIPTS/${RAN_ON}_${FILE}.sh

rm -f $SCRIPTS/${RAN_ON}_${FILE}.sh

done

