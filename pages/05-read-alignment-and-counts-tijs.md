
# Alignment to a reference genome

![workflow_align](../img/variant_calling_workflow_align.png)

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
$ cd ~/RNAseqWorkshop/general

$ hisat2-build -p 5 ath.fas ath --quiet

~~~
{: .bash}

While the index is created, you will see output something like this:

~~~
[bwa_index] Pack FASTA... 0.04 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 1.05 seconds elapse.
[bwa_index] Update BWT... 0.03 sec
[bwa_index] Pack forward-only FASTA... 0.02 sec
[bwa_index] Construct SA from BWT and Occ... 0.57 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index data/ref_genome/ecoli_rel606.fasta
[main] Real time: 1.765 sec; CPU: 1.715 sec
~~~
{: .output}

### Align reads to reference genome

creating the aligtnment (bam-files) is done in two steps. first the aligning it self with the use of hisat2. After that the alignment file will be filtered to only contain the reads that actualy map to the genome.This is done with [sam flags](https://broadinstitute.github.io/picard/explain-flags.html) in samtools view (with the '-F 4' all the unmapped reads will be removed).   

First of course we will need to create a directory to output the alignment files 

~~~
$ cd ~/RNAseqWorkshop/

$ mkdir mapped
~~~
{: .bash}

Next we have to enter the 'trimmed' directory, and loop through the fq files and get them aliugned/mapped to the indexed genome.

It's good again to first start with a 'dry' run with the use of echo

~~~
$ cd ~/RNAseqWorkshop/trimmed/

$For filename in *.fq
 do
    echo hisat2  -p 2 —dta -x ../general/ath -U filename | samtools view -Sb -F 4 -o ../mapped/"$(basename "$infile" .fq)”.bam
 done
~~~
{: .bash}

If the commands look good rerun but this time without the echo.

~~~
$For filename in *.fq
 do
    hisat2  -p 2 —dta -x ../general/ath -U filename | samtools view -Sb -F 4 -o ../mapped/"$(basename "$infile" .fq)”.bam
 done
~~~




~~~
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 77446 sequences (10000033 bp)...
[M::process] read 77296 sequences (10000182 bp)...
[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (48, 36728, 21, 61)
[M::mem_pestat] analyzing insert size distribution for orientation FF...
[M::mem_pestat] (25, 50, 75) percentile: (420, 660, 1774)
[M::mem_pestat] low and high boundaries for computing mean and std.dev: (1, 4482)
[M::mem_pestat] mean and std.dev: (784.68, 700.87)
[M::mem_pestat] low and high boundaries for proper pairs: (1, 5836)
[M::mem_pestat] analyzing insert size distribution for orientation FR...
~~~
{: .output}


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

![sam_bam1](../img/sam_bam.png)


![sam_bam2](../img/sam_bam3.png)
