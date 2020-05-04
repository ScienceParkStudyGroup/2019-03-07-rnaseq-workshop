
# Alignment to a reference genome

<img src="../images/RNAseqWorkflow.png" height="400" >

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to
choose from and, while there is no gold standard, there are some tools that are better suited for particular NGS analyses. We will be
using the [hisat2](http://ccb.jhu.edu/software/hisat2/index.shtml).

The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome


# Setting up



### Index the reference genome
Our first step is to index the reference genome for use by hisat2. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. Indexing the reference only has to be run once. The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment.

~~~
$ cd ~/RNAseq070319/general

$ mkdir genomeIndex

$ STAR --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles AtChromosome1.fa --runThreadN 8

~~~


While the index is created, you will see output something like this:

~~~
Apr 29 16:55:14 ..... Started STAR run
Apr 29 16:55:14 ... Starting to generate Genome files
Apr 29 16:55:16 ... starting to sort  Suffix Array. This may take a long time...
Apr 29 16:55:16 ... sorting Suffix Array chunks and saving them to disk...
Apr 29 16:56:13 ... loading chunks from disk, packing SA...
Apr 29 16:56:26 ... writing Suffix Array to disk ...
Apr 29 16:56:27 ... Finished generating suffix array
Apr 29 16:56:27 ... starting to generate Suffix Array index...
Apr 29 16:56:48 ... writing SAindex to disk
Apr 29 16:57:00 ..... Finished successfully

~~~



The indexing should have produced 8 star index files. Use the following command to see if they're really there.

~~~
$ ls -l genomeIndex/
~~~


result should be:
~~~
-rw-r--r-- 1 tbliek genseq-local          9 Apr 29 16:55 chrLength.txt
-rw-r--r-- 1 tbliek genseq-local         14 Apr 29 16:55 chrNameLength.txt
-rw-r--r-- 1 tbliek genseq-local          5 Apr 29 16:55 chrName.txt
-rw-r--r-- 1 tbliek genseq-local         11 Apr 29 16:55 chrStart.txt
-rw-r--r-- 1 tbliek genseq-local   30670848 Apr 29 16:55 Genome
-rw-r--r-- 1 tbliek genseq-local        290 Apr 29 16:55 genomeParameters.txt
-rw-r--r-- 1 tbliek genseq-local  249672325 Apr 29 16:56 SA
-rw-r--r-- 1 tbliek genseq-local 1565873616 Apr 29 16:56 SAindex

~~~




### Align reads to reference genome

creating the aligtnment (bam-files) is done in two steps. first the aligning it self with the use of hisat2. After that the alignment file will be filtered to only contain the reads that actualy map to the genome.This is done with [sam flags](https://broadinstitute.github.io/picard/explain-flags.html) in samtools view (with the '-F 4' all the unmapped reads will be removed).   

First of course we will need to create a directory to output the alignment files

~~~
$ cd ~/RNAseq070319/

$ mkdir mapped
~~~



Running STAR to align ( or map ) the reads and optionaly filter and sort them.

In contrast to most apps or programms, STAR does not have a help function.
running STAR -h or STAR --help will result in a n error. For information on what arguments to use you need to 
have a look at the manual.
[STAR](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).

~~~
$ cd ~/RNAseq070319/trimmed/

$  hisat2  -p 2 --dta -x ../general/ath -U sub06_qc.fq | samtools view -Sb -F 4 -o ../mapped/sub06_qc.bam
~~~



Next we want to make a loop to do all the files

It's good again to first start with a 'dry' run with the use of echo

~~~
$ cd ~/RNAseq070319/trimmed/

$for infile in *.fq
 do
   outfile="$(basename $infile .fq)”.bam
   echo "hisat2 -p 2 --dta -x ../general/ath -U $infile | samtools view -Sb -F 4 -o ../mapped/$outfile"
 done
~~~


If the commands look good rerun but this time without the echo.

~~~
$for infile in *.fq
 do
    outfile="$(basename $infile .fq)”.bam
    hisat2 -p 2 --dta -x ../general/ath -U $infile | samtools view -Sb -F 4 -o ../mapped/$outfile
 done
~~~


When running the hisat2 | samtools view, you will see output something like this:

~~~
899979 reads; of these:
  899979 (100.00%) were unpaired; of these:
    160982 (17.89%) aligned 0 times
    581749 (64.64%) aligned exactly 1 time
    157248 (17.47%) aligned >1 times
82.11% overall alignment rate
~~~



#### SAM/BAM format
The [SAM file](https://github.com/adamfreedman/knowyourdata-genomics/blob/gh-pages/lessons/01-know_your_data.md#aligned-reads-sam),
is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we do not
have time to go in detail of the features of the SAM format, the paper by
[Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

**The compressed binary version of SAM is called a BAM file.** We use this version to reduce size and to allow for *indexing*, which enables efficient random access of the data contained within the file.

The file begins with a **header**, which is optional. The header is used to describe source of data, reference sequence, method of
alignment, etc., this will change depending on the aligner being used. Following the header is the **alignment section**. Each line
that follows corresponds to alignment information for a single read. Each alignment line has **11 mandatory fields** for essential
mapping information and a variable number of other fields for aligner specific information. An example entry from a SAM file is
displayed below with the different fields highlighted.

![sam_bam_1](../images/sam_bam_1.png)

![sam_bam2](../images/sam_bam2.png)


### Creating the counts file

For downstream application for each of the samples the number of reads that maps within a gene has to be determent.
Featurecounts from the subread package can do this.


~~~
$ cd ~/RNAseqWorkshop/mapped

$ featureCounts -O -t mRNA -g ID -a ../general/annotation.all_transcripts.exon_features.ath.gff3 -o counts.txt *.bam
~~~



The file produced by featureCounts is a tab-delimited file.

For the ideal shiny that we will be using for the downstream analysis we need to make some adjustments.

the first line of the file containing the command used to run featureCounts, needs to be removed.(sill in the mapped directory)

~~~
$ awk '{if (NR!=1){print}} counts.txt > countsNew.txt
~~~



last we need to remove some columns that contain info on the gene (chromosome, startingposition, finalposition, length)

~~~
$ cut -f2,3,4,5,6 --complement countsNew.txt > countsFinal.txt
~~~



Now you should have a counts file (tab-delimited)

~~~
$ head countsNew.txt
~~~




| Geneid  sub06_qc.bam    sub07_qc.bam    sub08_qc.bam    sub21_qc.bam    sub23_qc.bam    sub24_qc.bam |
|------------------------------------------------------------------------------------------------------|
| AT1G01010       0       0       6       6       3       10                                           |
| AT1G01020       2       4       3       1       2       1                                            |
| AT1G03987       0       0       0       0       0       0                                            |
| AT1G01030       0       0       0       0       2       1                                            |
| AT1G03993       0       0       0       0       0       0                                            |
| AT1G01040       17      17      22      3       3       4                                            |
| AT1G01046       0       0       0       0       0       0                                            |
| ath-miR838      0       0       0       0       0       0                                            |
| AT1G01050       16      23      15      21      22      30                                           |
| AT1G03997       0       0       0       0       0       0                                            |
| AT1G01060       0       0       1       0       0       2                                            |
| AT1G01070       0       1       2       1       2       4                                            |
| AT1G04003       0       0       0       0       0       0                                            |
| AT1G01080       51      28      26      25      25      21                                           |
| AT1G01090       108     90      103     30      46      32                                           |
