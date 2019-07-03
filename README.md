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
We realized that there were many repeated sequences and that they were due to the adapters that are put to the extracted DNA when it is processed by ILUMINA, which is what produces the sequencing. Then, in order to obtain the cleanest possible sequencing and the most compressed we decided to do quality control, which is recommended to do, even once. To explore what different adapters were present in the samples we run the script:

```
front=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
end=ATCTCGTATGCCGTCTTCTGCTTG

#zgrep -o -e "${front}......${end}" $*

zcat $* | sed -n 's/.*GATCGGAAGAGCACACGTCTGAACTCCAGTCAC\(......\)ATCTCGTATGCCGT$$TGCCGTCTTCTGCTTG.*/\1/p'
```
Then knowing that we ran the script Humberto wrote to trim the data and eliminating the adapters. Is important to notice that we ran the trim script with a quality of 5, so that way we would eliminate all sequences less that a quality of 5.

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
```
File called “Combinados” that is run inside the script trim.
```

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
## Beginning Digital Normalization

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
## Note
Also, we found in the web page of IMG under the name of “Tabebuia heterophylla” the assemblies of the project made by other person.  (https://img.jgi.doe.gov/cgi-bin/m/main.cgi?section=FindGenomes&page=displayTaxonList&searchFilter=all&searchTerm=Tabebuia%20heterophylla&file=all78716&allDataFiltersFile=allGenomeDataFilters78716)

## Interleave

Well, the command that we gave with the interleave.sh script did not work. Because when we ran the script that we wrote to do normalize, it told us that the samples were not paired.
Therefore, we had to download a copy of the khmer tool script and modified, so that when we ran our interleave.sh script and got the results, we could then run normalize.sh script.

To explain it better, I will use the same samples as visual. All this happened since we ran the quality control with the fastq-dump script and fasterq-dump. We ran the first 2 samples (SRR5256888 and SRR5256985) with fastq-dump and they came out this way:

Example of sample #1 SRR5256888_1.qc.fq.gz 

```
@SRR5256888.1 1 length=150
CGAGACTCTTGCGCGCGGCGGCGTCGGACA
+SRR5256888.1 1 length=150
CCCFFFFFHHHHHJJJJJJJJGD@BDDDDD
@SRR5256888.2 2 length=150
GTCGCCTCTTCCGATCCAATTGCACAAGCA
+SRR5256888.2 2 length=150
CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJ
````

But to explore whether quality control could be simplified and much better, we used the fasterq-dump, that is the update of fastq-dump. It came up this way:

Sample #3 SRR5256987_1.qc.fq.gz 
```
@SRR5256987.12.1 12 length=150
AGCTCCTCCCACCAGGTCGCAGCGCGACCGTCGAGCTCTCCGTAGACGCCGGTCGCCTCCTTCCGCTGATCGTCGACGAGCACGACCCCGGCGAACCCAACCCGCAGGAGCGAGGTCACCTCCCGAACGAGCGGCCGCGCCACG
+SRR5256987.12.1 12 length=150
1=DDFFFHFHHHJJJJIGIJJIGHIJIJJJJJJHHFFFFDEDDDDDDDDBDDBDDDDDBDDDDDDBDDDDDDDDBDDDDDBDDDDDDDDDDDDDDBBBDDDDB@>B@<<<89>@B>CCCDCBDD>BB.99<<@DBBB@<BDBBD
```

The difference between the samples, is the name of each sample. Since each program put it in a different way. As we can see the first sample is called:
```
@ SRR5256888.1 1 length = 150
``` 

and the third sample:
```
@ SRR5256987.12.1 12 length = 150
```

This one has an extra period and an extra number compared to the names of the other samples. We decided to run the fourth sample with fastq-dump. 

The problem is that as mentioned previously the 3rd sample gives us an error and for this small problem we had to modify the scripts.

Therefore, we order conda to install Khmer, in order to modify the script of Khmer tools, so interleave could process the data and fix the files names.  The script is called: interleave-reads.py

```
#!/home/humberto/adominicci/miniconda3/envs/roble/bin/python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
# 	 derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=invalid-name,missing-docstring
""" 
Interleave left and right reads.

Take two files containing left & right reads from a paired-end sequencing run,
and interleave them.

% scripts/interleave-reads.py <R1> <R2> [ -o <outputfile> ]

By default, output is sent to stdou
"""

import screed
import sys
import textwrap
from khmer import __version__
from khmer.kfile import check_input_files, check_space
from khmer.khmer_args import sanitize_help, KhmerArgumentParser
from khmer.khmer_args import FileType as khFileType
from khmer.kfile import (add_output_compression_type, get_file_writer,
from khmer.utils import (write_record_pair, check_is_left, check_is_right,
                         check_is_pair)

try:
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest


def get_parser():
    epilog = """\
    The output is an interleaved set of reads, with each read in <R1> paired
    with a read in <R2>. By default, the output goes to stdout unless
    :option:`-o`/:option:`--output` is specified.

    As a "bonus", this file ensures that if read names are not already
    formatted properly, they are reformatted consistently, such that
    This reformatting can be switched off with the
    :option:`--no-reformat` flag.

    Example::

        interleave-reads.py tests/test-data/paired.fq.1 \\
                tests/test-data/paired.fq.2 -o paired.fq"""
    parser = KhmerArgumentParser(
  description='Produce interleaved files from R1/R2 paired files',
        epilog=textwrap.dedent(epilog))

    parser.add_argument('left')
    parser.add_argument('right')
    parser.add_argument('-o', '--output', metavar="filename",
                        type=khFileType('wb'),
                        default=sys.stdout)
    parser.add_argument('--no-reformat', default=False, action='store_true',
                        help='Do not reformat read names or enforce\
                              consistency')
    parser.add_argument('-f', '--force', default=False, action='store_true',
                        help='Overwrite output file if it exists')
    add_output_compression_type(parser)
    return parser


def main():
    args = sanitize_help(get_parser()).parse_args()

    check_input_files(args.left, args.force)
    check_input_files(args.right, args.force)
    check_space([args.left, args.right], args.force)

    s1_file = args.left
    s2_file = args.right

    print("Interleaving:\n\t%s\n\t%s" % (s1_file, s2_file), file=sys.stderr)

    outfp = get_file_writer(args.output, args.gzip, args.bzip)

    counter = 0
    screed_iter_1 = screed.open(s1_file)
    screed_iter_2 = screed.open(s2_file)
    for read1, read2 in zip_longest(screed_iter_1, screed_iter_2):
        if read1 is None or read2 is None:

            print(("ERROR: Input files contain different number"
                   " of records."), file=sys.stderr)
            sys.exit(1)

        if counter % 100000 == 0:
            print('...', counter, 'pairs', file=sys.stderr)
        counter += 1

        name1 = read1.name.split(" ")[0]
        name2 = read2.name.split(" ")[0]
        fields = name1.split(".")
        if (3 == len(fields)) :
           name1= ".".join(fields [:2])
        fields = name2.split(".")
        if (3 == len(fields)) :
           name2= ".".join(fields [:2])

        if not args.no_reformat:
            if not check_is_left(name1):
                name1 += '/1'
            if not check_is_right(name2):
                name2 += '/2'

            read1.name = name1
            read2.name = name2

            if not check_is_pair(read1, read2):
                print("ERROR: This doesn't look like paired data! "
                      "%s %s" % (read1.name, read2.name), file=sys.stderr)
                sys.exit(1)

        write_record_pair(read1, read2, outfp)

    print('final: interleaved %d pairs' % counter, file=sys.stderr)
    print('output written to', describe_file_handle(outfp), file=sys.stderr)
if __name__ == '__main__':
    main()
```

We specifically modified this part:
```

name1 = read1.name.split(" ")[0]
        name2 = read2.name.split(" ")[0]
        fields = name1.split(".")
        if (3 == len(fields)) :
           name1= ".".join(fields [:2])
        fields = name2.split(".")
        if (3 == len(fields)) :
           name2= ".".join(fields [:2])
```

That is a command to modify the name of the sample, in order to process and that at the end t all the names of each sample have the same format. We commanded that the program split the name as follow @ SRRxxx.yy -yy- lengh = 150 and eliminate the part after the second period. Then we include the modified script called "interleave-reads.py" in our script called: "interleve.sh":
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
  ./interleave-reads.py -o ${base}.pe.qc.fq.gz ${baseR1}.qc.fq.gz ${baseR2}.qc.$
Done
```

What then turned out to be:
```
SRR5256888.pe.qc.fq

@SRR5256888.1/1
CGAGACTCTTGCGCGCGGCGGCGTCGGACA
+
CCCFFFFFHHHHHJJJJJJJJGD@BDDDDD
@SRR5256888.1/2
AAGGCCCGTGACAGCCACGAGCTTATCGTCGATGGCTAGTTCGCCCGTCACCTTGAAGAAGGGCTGGCTGTAGTAATACGAGGCCTGCTCGCGCTCGGACTTCTTGCTGAAGCCGCCTTCGCCCTGCAGCACCAGCGGACGATC
+
CCCFFFFFHHHHHJJJJJJJJIJJJJJJIJJIIJIJJIIGHIIIJJGHHFBDEFACEEEEDDDDDDDDDDBCDCCDCCCDDDDDDDDDDCDDDDDDDD>BDDDDDCCD@CCACDDDD@@C?BDD>ABCCCDDB88A>9@<>9<B
```

## Normalize Script
Therefore, we decided to run the normalize script "normalize-by-median.py"

normalize.sh
```

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

normalize-by-median.py -k 20 -C 20 --max-memory-usage 64G -p --savegraph normC2$normC20k20.ct --report report.txt --report-frequency 1000000 *.pe.qc.fq
```

For the script we modify the memory and the time, because we know that is going to take a lot of memory and time.

## Next Step

Right now, we are waiting for those results, and we thinking that our next step should be to do a partition (which is what we will probably do) or run the assemblies again (Megahit, PLASS and MetaSPAdes).

