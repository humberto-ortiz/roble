#!/bin/bash
#SBATCH --mem-per-cpu=512
#SBATCH --time=48:00:00
#SBATCH --job-name=fastqc-qc
#SBATCH --mail-user=humberto.ortiz@upr.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Load project
source /home/humberto/humberto/miniconda3/etc/profile.d/conda.sh
conda activate roble

fastqc --threads=8 -outdir=recommended *.qc.fq.gz
