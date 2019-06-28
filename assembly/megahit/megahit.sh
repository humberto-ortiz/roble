#!/bin/bash
#SBATCH --mem-per-cpu=15000
#SBATCH --time=48:00:00
#SBATCH --job-name=megahit
#SBATCH --mail-user=anelisse.dominicci@upr.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Load project
source /home/humberto/adominicci/miniconda3/etc/profile.d/conda.sh
conda activate roble

#echo ~/miniconda3/pkgs/megahit-1.0.3-0/megahit -1 $(ls *_1.qc.fq.gz | tr '\n' ,) -2 $( ls *_2.qc.fq.gz | tr '\n' ,) orphans.qc.fq.gz -o megahit_out -t 8

/home/humberto/adominicci/miniconda3/pkgs/megahit-1.0.3-0/megahit \
  -1 SRR5256888_1.qc.fq.gz,SRR5256985_1.qc.fq.gz,SRR5256987_1.qc.fq.gz,SRR5256988_1.qc.fq.gz \
  -2 SRR5256888_2.qc.fq.gz,SRR5256985_2.qc.fq.gz,SRR5256987_2.qc.fq.gz,SRR5256988_2.qc.fq.gz orphans.qc.fq.gz \
  -o megahit_out -t 8
