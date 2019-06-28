#!/bin/bash
#SBATCH --mem-per-cpu=1024
#SBATCH --time=48:00:00
#SBATCH --job-name=interleave
#SBATCH --mail-user=anelisse.dominicci@upr.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Load project
source /home/humberto/adominicci/miniconda3/etc/profile.d/conda.sh
conda activate roble

for filename in *_1.qc.fq.gz
do
  # first, make the base by removing fastq.gz
  base=$(basename $filename _1.qc.fq.gz)
  echo $base

  # now, construct the R2 filename by replacing R1 with R2
  baseR2=${base}_2
  baseR1=${base}_1
  echo $baseR2

 # finally, run interleave-reads.py
  interleave-reads.py --no-reformat --gzip -o ${base}.pe.qc.fq.gz ${baseR1}.qc.fq.gz ${baseR2}.qc.fq.gz
done
