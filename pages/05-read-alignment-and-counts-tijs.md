
# Alignment to a reference genome

<img src="../images/RNAseqWorkflow.png" height="400" >

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to
choose from and, while there is no gold standard, there are some tools that are better suited for particular NGS analyses. In this tutorial we will be using [STAR](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) but also 
a tool like [hisat2](http://ccb.jhu.edu/software/hisat2/index.shtml) does the job.

The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome


# Setting up



### Index the reference genome
Our first step is to index the reference genome for use by STAR. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. Indexing the reference only has to be run once. The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment (index files are not exchangeable between tools).

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



The indexing should have produced 8 star index files. Use the following command to see if they're really there. Take note that depending on the genome size these files can be pretty big. Be sure to have enough disk space.

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

In some tools like hisat2 creating the sequence alignment files (bam-files) is done in two steps. first the aligning it self. After that the alignment file will be filtered for instance to only contain the reads that actualy map to the genome. This is done with [sam flags](https://broadinstitute.github.io/picard/explain-flags.html) in samtools view (with the '-F 4' all the unmapped reads will be removed). STAR on the other hand has a build in filter and also a sort function. So the output is ready to use for downstream tools.  



First of course we will need to create a directory to output the alignment files

~~~
$ cd ~/RNAseq070319/

$ mkdir mapped
~~~


Running STAR to align ( or map ) the reads and optionaly filter and sort them.

In contrast to most apps or programms, STAR does not have a help function.
running STAR -h or STAR --help will result in a n error. For information on what arguments to use you can or need to 
use have a look at the 
[STAR manual.](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).


Here are some examples of comman used arguments.

| argument   | meaning |
| ------- | ---------- |
| `--runThreads` | number of threads |
| `--genomeDir` | /path/to/genomeDir |
| `--readFilesIn` | /path/to/read1 [/path/to/read2] |
| `--readFilesCommand zcat` | when make use of gzipped fastq files |
| `--outFileNamePrefix` | /path/to/output file name |
| `--outSAMtype` | BAM/SAM or None  [optional: SortedByCoordinate] |
| `--outReadsUnmapped` | [default: None] Fastx ; output in separate fasta/fastq file |
| `--outFilterMultimapNmax` | [default: 10] max number of alignments accepted |
| `--outFilterMismatchNmax` | [default: 10] max number of mismatches accepted |
| `--outFilterMismatchNoverLmax`  | [default: 0.3] max fraction of mismatches mapped length |
| `--outFilterMismatchNoverReadLmax` | [default: 1.0] max fraction of mismatches read length |
| `--alignEndsType` | EndToEnd force end-to-end alignment, don't soft-clip |



For now we will be using STAR with the following arguments
~~~

$  STAR --genomeDir genomeindex --runThreadN 2 --readFilesIn {i} --readFilesCommand zcat --outFileNamePrefix --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --outFilterMultimapNmax 1

~~~



Next we want to make a loop to do all the files

It's good again to first start with a 'dry' run with the use of echo

~~~

$ for infile in trimmed/*.fq
 do
   outfile="$(basename $infile .fq)”
   echo "STAR --genomeDir genomeIndex --runThreadN 2 --readFilesIn trimmed/$infile --readFilesCommand zcat --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --alignEndsType EndToEnd --outFilterMultimapNmax 1"
 done
~~~


If the commands look good, rerun but this time without the echo.

~~~
$for infile in *.fq
 do
   outfile="$(basename $infile .fq)”
   STAR --genomeDir genomeIndex --runThreadN 2 --readFilesIn trimmed/$infile --readFilesCommand zcat --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --alignEndsType EndToEnd --outFilterMultimapNmax 1
 done
~~~


When running the STAR command, you will see output something like this:

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


