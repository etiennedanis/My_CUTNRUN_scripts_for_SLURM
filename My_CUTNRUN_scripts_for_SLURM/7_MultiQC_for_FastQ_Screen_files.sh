#!/usr/bin/bash

#SBATCH -J Multiqc
#SBATCH -e Multiqc_on_FastQ_Screen_output.err
#SBATCH -o Multiqc_on_FastQ_Screen_output.out

#SBATCH --nodes=1                            # Number of nodes
#SBATCH --mem=4G                             # RAM per node
#SBATCH --cpus-per-task=2                    # CPUs per node

#SBATCH --time=0-00:05:00

INVESTIGATOR=
PROJECT=

CUTNRUN=/scratch/myaccount/Documents/cutnrun_data
DATA=$CUTNRUN/${INVESTIGATOR}_data/$PROJECT
FASTQC=$DATA/fastq_screen_files_w_mycoplasma

module load multiqc/1.14

cd $FASTQC

multiqc .


