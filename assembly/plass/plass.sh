#!/bin/bash
#SBATCH --mem-per-cpu=640
#SBATCH --time=48:00:00
#SBATCH --job-name=plass
#SBATCH --mail-user=anelisse.dominicci@upr.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Load project
source /home/humberto/adominicci/miniconda3/etc/profile.d/conda.sh
conda activate roble


plass assemble SRR*_[12].qc.fq.gz assembly.fas tmp

