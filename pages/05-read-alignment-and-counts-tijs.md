
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

$ hisat2-build -p 5 ath.fas ath --quiet

~~~


While the index is created, you will see output something like this:

~~~
Settings:
  Output files: "ath.*.ht2"
  Line rate: 6 (line is 64 bytes)
  Lines per side: 1 (side is 64 bytes)
  Offset rate: 4 (one in 16)
  FTable chars: 10
 …
  numLines: 622300
  gbwtTotLen: 39827200
  gbwtTotSz: 39827200
  reverse: 0
  linearFM: Yes
Total time for call to driver() for forward index: 00:01:09

~~~



The indexing should have produced 8 hisat2 index files (.ht2). Use the following command to see if they're really there.

~~~
$ ls *.ht2
~~~


result should be:
~~~
ath.1.ht2  ath.2.ht2  ath.3.ht2  ath.4.ht2  ath.5.ht2  ath.6.ht2  ath.7.ht2  ath.8.ht2
~~~




### Align reads to reference genome

creating the aligtnment (bam-files) is done in two steps. first the aligning it self with the use of hisat2. After that the alignment file will be filtered to only contain the reads that actualy map to the genome.This is done with [sam flags](https://broadinstitute.github.io/picard/explain-flags.html) in samtools view (with the '-F 4' all the unmapped reads will be removed).   

First of course we will need to create a directory to output the alignment files

~~~
$ cd ~/RNAseq070319/

$ mkdir mapped
~~~



Running hisat2 to align ( or map ) the reads and pipe the result through samtools view to remove the non-mapping reads.

~~~
$ cd ~/RNAseq070319/trimmed/

$  hisat2  -p 2 -x ../general/ath -U sub06qc.fq | samtools view -Sb -F 4 -o ../mapped/sub06_qc.bam
~~~



Next we want to make a loop to do all the files

It's good again to first start with a 'dry' run with the use of echo

~~~
$ cd ~/RNAseq070319/trimmed/

$For filename in *.fq
 do
   outfile="$(basename "$infile" .fq)”.bam
   echo hisat2 -p 2 —dta -x ../general/ath -U filename | samtools view -Sb -F 4 -o ../mapped/$outfile
 done
~~~


If the commands look good rerun but this time without the echo.

~~~
$For filename in *.fq
 do
    hisat2  -p 2 —dta -x ../general/ath -U filename | samtools view -Sb -F 4 -o ../mapped/"$(basename "$infile" .fq)”.bam
 done
~~~


When running the hisat2 | samtools view, you will see output something like this:

~~~
output of hisat2





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

![sam_bam1](../images/sam_bam1.png)

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
