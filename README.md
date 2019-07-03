# Introduction
The project consists of analyzing those 4-shotgun datasets of the microbial community from the rhizospheres of the *Tabebuia heterophylla*. The main objective is to find microbial differences between the 4 soil samples of the “roble” *Tabebuia heterophylla*. The samples were collected in 2012, in 4 towns in Puerto Rico, and in 3 different soil types. The towns are Guayama; in low altitude volcanic soil, Guánica; in low altitude limestone soil, Maricao; in high altitude serpentine soil and Cabo Rojo; in low altitude serpentine soil. They collected a sample of each one, except for the serpentine we have from a high and low altitude zone. We want to analyze and establish the differences in microbial composition, of bacteria, fungi, prokaryotes and microeukaryotes between different sites and geology. Once we can identify those groups, we will establish which genes are present in each one, and what do they code. Also, we want to identify the metabolism of these genomes, see if there are groups present with nitrogen fixing capacity, also to be able to identify microbes that convert plant biomass to sugars, identify genes of Crispr/cas and genes resistant to antibiotics. In order to that we are doing the following:

# What are we doing?
In order to carry out the research, we went to GOLD, which is where the samples is available, since “Genomes Online Database”, which is a part of the World Wide Web, it is a public database that gives us access to the genome and metagenome sequencing projects uploaded there (Mukherjee, Stamatis et al.).  However, the data was extracted from SRA, which is a bioinformatics database that functions as a public repository for DNA sequencing data which provides shorts reads (Sequence Read Archive). In order to understand the first part of the experiment and be able to continue with the metagenomic part of the project that is what we are going to do this summer. For doing that we have several steps. The first step is that I will read the protocols of the programs they used to carry out the genome assembly, which in this case was SPAdes an assemble that was designed for single cell and multi-cells sequencing, providing information about genomes (Bankevich et al.). In addition, those used for gene calling methods: FragGeneScan version 1.16, which is a gene predictor, functions as predicting fragmented genes and genes with frameshifts in addition to the complex genes (Rho et al.), Prokaryotic GeneMark.hmm version 2.8, also a gene predictor, designed to improve the gene prediction quality in terms of finding exact gene boundaries (Lukashin and Borodovsky), Metagene Annotator version 1.0, which predicts all kinds of prokaryotic genes froms a single or a set of genomic sequences (Noguchi, Tanigu) (Noguchi, Park, et al.) and Prodigal version 2.5, that can identify genes in short coding sequences with a high degree of accuracy  (Hyatt et al.).

When I finally end reading all the papers, to acknowledge what and how they did the sequencing, we start working in the terminal.

In the terminal I got access to the mega-computer called “Boqueron” and in that way we then can downloaded the data from SRA to “Boqueron". Before having access to it, we downloaded miniconda3, which is a system that handles packages that contain many programs that perform different functions. 
When we saw the data, we saw that it was very large, the files had approximate 60GB, as I mentioned there are 4 samples of 60GB each, they were too big. We spent a whole week trying to download the files, because the connection between the computer “NCBI” and “Boqueron” was failing. Therefore, Humberto had to download them to a supercomputer called "Hulk" and then transfer them to "Boquerón". Then we run a command called “fastq-dump” to compress the data. We try doing “fasterq-dump” but it didn’t work. We achieve to compress each file approximate to 20GB. We open our environment called “Roble” and there we created a directory called "quality" that had 2 files one called “raw” (where we deposited the raw data” and the other one “trimmed” (where we are going to deposit the data after the quality control.
## Quality Control
Therefore, Humberto wrote a script to be able to see a report of the samples before doing the quality control, to see how they were and compare them with the results once they went through the quality control. This script ran the program called "FASTQC":

```
#!/bin/bash
#SBATCH --mem-per-cpu=512
#SBATCH --time=48:00:00
#SBATCH --job-name=fastqc-qc
#SBATCH --mail-user=anelisse.dominicci@upr.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Load project
source /home/humberto/adominicci/miniconda3/etc/profile.d/conda.sh
conda activate roble

fastqc --threads=8 -outdir=trim *.qc.fq.gz

```

After about 2 hours, we obtained the reports of each sample that they told us were in poor condition. Then we ran "MULTIQC" which is another program that gives us a general report of all the samples together. We obtained the same results, but in a general way.
We realized that there were many repeated sequences and that they were due to the adapters that are put to the extracted DNA when it is processed by ILUMINA, which is what produces the sequencing. Then, in order to obtain the cleanest possible sequencing and the most compressed we decided to do quality control, which is recommended to do, even once. We ran the script Humberto wrote to trim the data and eliminating the adapters. Is important to notice that we ran the trim script with a quality of 5, so that way we would eliminate all sequences less that a quality of 5.
```

#!/bin/bash
#SBATCH --mem-per-cpu=1024
#SBATCH --time=48:00:00
#SBATCH --job-name=Trimmomatic
#SBATCH --mail-user=anelisse.dominicci@upr.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Load project
source /home/humberto/adominicci/miniconda3/etc/profile.d/conda.sh
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
    ILLUMINACLIP:Combinados:2:40:15 \
    LEADING:5 TRAILING:5 \
    SLIDINGWINDOW:4:5 \
    MINLEN:25

  # save the orphans
  gzip -9c s1_se s2_se >> orphans.qc.fq.gz
  rm -f s1_se s2_se
done


Adapters:

>PrefixIndexR1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>PrefixIndexR2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>PrefixIndexR1rev
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>PrefixIndexR2rev
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
```


We waited for a full day and got the trim, then we ran the FASTQC and MULTIQC scripts again with the data trim, to see how the samples had been given and decide if we were going to run the trim again.
Seeing the report, we decided that it was all good, the samples were in good condition, so we moved the data to the directory “quality/trimmed”.

## Assemblies
We decided that it was time to create the assemblies, referring to aligning and merging fragments from a longer DNA sequence in order to reconstruct the original sequence. We had to install the three assemblers that we would use with Conda, creating a folder for each one. We will use Megahit, an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph, and we run it with this script:
```
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

#echo ~/miniconda3/pkgs/megahit-1.0.3-0/megahit -1 $(ls *_1.qc.fq.gz | tr '\n' $

/home/humberto/adominicci/miniconda3/pkgs/megahit-1.0.3-0/megahit \
  -1 SRR5256888_1.qc.fq.gz,SRR5256985_1.qc.fq.gz,SRR5256987_1.qc.fq.gz,SRR52569$
  -2 SRR5256888_2.qc.fq.gz,SRR5256985_2.qc.fq.gz,SRR5256987_2.qc.fq.gz,SRR52569$
  -o megahit_out -t 8
```

We had several problems creating the script, but with a little of work we solve in a way that the program could read it.
Then we did the same with PLASS, which is another assembler, but it translates the sequences into protein sequences.

```
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

```

Finally, we did the same with MetaSpades, which is an assembler too, which is part of Spades. We had to join all the _1 sequence and put it to getter in a folder called “left”. And then the sequences that are the reverse _2 in a folder called “right”. Since this assembler read the command/script, for 2 folders to create the assembly.
```
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
```
The 3-assembler failed, because they needed more memory than the one, we give to each one. Therefore, we decided to do a "Digital Normalization" that eliminates possible errors and reduce the sequences if there are too many or too few. This process is also part of the quality control. 
```
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

# Delete any old orphan reads
rm -f orphans.qc.fq.gz

for filename in *_1.qc.fq.gz
do
  # first, make the base by removing fastq.gz
base=$(basename $filename .qc.fq.gz)
  echo $base

  # now, construct the R2 filename by replacing R1 with R2
  baseR2=${base/_1/_2}
  echo $baseR2

 # finally, run interleave-reads.py
  interleave-reads.py -o ${base}.pe.qc.fq.gz ${base}.qc.fq.gz ${baseR2}.qc.fq.gz
done

```


The “Digital Normalization” (interleave-reads.py) was done, but we had to eliminate the results. First because the 3rd sample gives us an error. And secondly, because in the first script and we forgot to put --gzip and --no-reformat (so when we runned again, it doesn't change the names of the sequences, nor check if they are paired sequences. So here the new script with the corrections:

```
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

# finally, run interleave-reads.py
  interleave-reads.py --no-reformat --gzip -o ${base}.pe.qc.fq.gz ${baseR1}.qc.$
done

```
# Note
Also, we found in the web page of IMG under the name of “Tabebuia heterophylla” the assemblies of the project made by other person.  (https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=FindGenomes&page=displayTaxonList&searchFilter=all&searchTerm=Tabebuia%20heterophylla&file=all78716&allDataFiltersFile=allGenomeDataFilters78716)
