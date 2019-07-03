#!/bin/bash
#SBATCH --mem-per-cpu=15000
#SBATCH --time=100:00:00
#SBATCH --job-name=normalize
#SBATCH --mail-user=anelisse.dominicci@upr.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Load project
source /home/humberto/adominicci/miniconda3/etc/profile.d/conda.sh
conda activate roble

normalize-by-median.py -k 20 -C 20 --max-memory-usage 64G -p --savegraph normC20k20.ct --report report.txt --report-frequency 1000000 *.pe.qc.fq
