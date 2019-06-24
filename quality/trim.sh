#!/bin/bash
#SBATCH --mem-per-cpu=1024
#SBATCH --time=48:00:00
#SBATCH --job-name=Trimmomatic
#SBATCH --mail-user=humberto.ortiz@upr.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Load project
source /home/humberto/humberto/miniconda3/etc/profile.d/conda.sh
conda activate roble

# Delete any old orphan reads
rm -f orphans.qc.fq.gz

for filename in *_1.fastq.gz
do
  # first, make the base by removing fastq.gz
  base=$(basename $filename .fastq.gz)
  echo $base

  # now, construct the R2 filename by replacing R1 with R2
  baseR2=${base/_1/_2}
  echo $baseR2

 # finally, run Trimmomatic
  trimmomatic PE -threads 8 ${base}.fastq.gz ${baseR2}.fastq.gz \
    ${base}.qc.fq.gz s1_se \
    ${baseR2}.qc.fq.gz s2_se \
    ILLUMINACLIP:adapters.txt:2:40:15 \
    LEADING:2 TRAILING:2 \
    SLIDINGWINDOW:4:2 \
    MINLEN:25

  # save the orphans
  gzip -9c s1_se s2_se >> orphans.qc.fq.gz
  rm -f s1_se s2_se
done
