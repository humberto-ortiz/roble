#!/bin/bash
#SBATCH --mem-per-cpu=16000
#SBATCH --time=48:00:00
#SBATCH --job-name=metaspades
#SBATCH --mail-user=anelisse.dominicci@upr.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Load project
source /home/humberto/adominicci/miniconda3/etc/profile.d/conda.sh
conda activate roble

metaspades.py -o metaspades_out --pe1-1 left.fq.gz --pe1-2 right.fq.gz
