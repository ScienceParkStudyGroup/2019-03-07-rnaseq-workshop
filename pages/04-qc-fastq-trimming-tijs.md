# Quality control and trimming



## Running FastQC  

We will now assess the quality of the reads that we downloaded. First, we need to make an output directory for the fastqc results to be stored. This we want to do in the 'RNAseq070319' directory.

~~~
$ cd ~/RNAseq070319

$ mkdir fastqc
~~~
{: .bash}

Next we need to get to the directory thay actually contains the the fastq files.

~~~
$ cd ~/RNAseq070319/rawReads
~~~
{: .bash}


Running fastqc uses the following command

~~~
fastqc -o ../fastqc $filename
~~~
{: .bash}

Of course we don't want to do y=this for all the samples seperately so we can loop through the list of samples and run them all sequentially

With the use of echo you can start off with a "dry run"

~~~
$ for filename in *.fastq
  do
    echo fastqc -o ../fastqc $filename
  done
~~~
{: .bash}

The echo command only prints the commands to the screen, and doesn't really run it.

~~~
fastqc -o ../fastqc sub06.fastq
fastqc -o ../fastqc sub07.fastq
fastqc -o ../fastqc sub08.fastq
fastqc -o ../fastqc sub21.fastq
fastqc -o ../fastqc sub23.fastq
fastqc -o ../fastqc sub24.fastq
~~~


If it looks good remove the echo and go for it.

~~~
$ for filename in *.fastq
  do
    fastqc -o ../fastqc $filename
  done
~~~
{: .bash}

You will see an automatically updating output message telling you the
progress of the analysis. It shoul look something like this:

~~~
Started analysis of sub06.fastq
Approx 5% complete for sub06.fastq
Approx 10% complete for sub06.fastq
Approx 15% complete for sub06.fastq
Approx 20% complete for sub06.fastq
…
Approx 85% complete for sub24.fastq
Approx 90% complete for sub24.fastq
Approx 95% complete for sub24.fastq
Approx 100% complete for sub24.fastq
Analysis complete for sub24.fastq
~~~


In total, it should take about five minutes for FastQC to run on all
six of our FASTQ files.

If the command doesn't run or you want more information on fastqc, run the following to get the help page.

~~~
$ fastqc -h
~~~

But if all went right, the FastQC program will have created several new files within our
`~/RNAseq070319/fastqc` directory.

~~~
$ cd ~/RNAseq070319/fastqc
$ ls
~~~
{: .bash}

~~~
sub06_fastqc.html  sub07_fastqc.zip   sub21_fastqc.html  sub23_fastqc.zip
sub06_fastqc.zip   sub08_fastqc.html  sub21_fastqc.zip   sub24_fastqc.html
sub07_fastqc.html  sub08_fastqc.zip   sub23_fastqc.html  sub24_fastqc.zip
~~~



## Viewing the FastQC results

If we were working on our local computers, we'd be able to display each of these
HTML files as a webpage:

~~~
$ open sub06_fastqc.html
~~~
{: .bash}

However, if you try this on our genseq instance, you'll get an error:

~~~
Couldn't get a file descriptor referring to the console
~~~


This is because the genseq instance we're using doesn't have any web
browsers installed on it, so the remote computer doesn't know how to
open the file. We want to look at the webpage summary reports, so
let's transfer them to our local computers (i.e. your laptop).

To transfer a file from a remote server to our own machines, we will
use `scp`.

First we
will make a new directory on our computer to store the HTML files
we're transfering. Let's put it on our desktop for now. Open a new
tab in your terminal program (you can use the pull down menu at the
top of your screen or the Cmd+t keyboard shortcut) and type:

~~~
$ mkdir -p ~/Desktop/fastqc_html
~~~
{: .bash}

Now we can transfer our HTML files to our local computer using `scp`.

~~~
$ scp tbliek@genseq-cn02.science.uva.nl:~/RNAseq070319/fastqc/*.html ~/Desktop/fastqc_html
~~~
{: .bash}

As a reminder, the first part
of the command `tbliek@genseq-cn02.science.uva.nl` is
the address for your remote computer. Make sure you replace everything
after `dcuser@` with your instance number (the one you used to log in).

The second part starts with a `:` and then gives the absolute path
of the files you want to transfer from your remote computer. Don't
forget the `:`. We used a wildcard (`*.html`) to indicate that we want all of
the HTML files.

The third part of the command gives the absolute path of the location
you want to put the files. This is on your local computer and is the
directory we just created `~/Desktop/fastqc_html`.

You should see a status output like this:

~~~
SRR2584863_1_fastqc.html                      100%  249KB 152.3KB/s   00:01    
SRR2584863_2_fastqc.html                      100%  254KB 219.8KB/s   00:01    
SRR2584866_1_fastqc.html                      100%  254KB 271.8KB/s   00:00    
SRR2584866_2_fastqc.html                      100%  251KB 252.8KB/s   00:00    
SRR2589044_1_fastqc.html                      100%  249KB 370.1KB/s   00:00    
SRR2589044_2_fastqc.html                      100%  251KB 592.2KB/s   00:00  
~~~


Now we can go to our new directory and open the HTML files.

~~~
$ cd ~/Desktop/fastqc_html/
$ open *.html
~~~
{: .bash}

Your computer will open each of the HTML files in your default web
browser. Depending on your settings, this might be as six separate
tabs in a single window or six separate browser windows.

> ## Exercise
>
> Discuss your results with a neighbor. Which sample(s) looks the best
> in terms of per base sequence quality? Which sample(s) look the
> worst?
>
>> ## Solution
>> All of the reads contain usable data, but the quality decreases toward
>> the end of the reads.
> {: .solution}
{: .challenge}

## Decoding the other FastQC outputs
We've now looked at quite a few "Per base sequence quality" FastQC graphs, but there are nine other graphs that we haven't talked about! Below we have provided a brief overview of interpretations for each of these plots. It's important to keep in mind

+ **Per tile sequence quality**: the machines that perform sequencing are divided into tiles. This plot displays patterns in base quality along these tiles. Consistently low scores are often found around the edges, but hot spots can also occur in the middle if an air bubble was introduced at some point during the run.
+ **Per sequence quality scores**: a density plot of quality for all reads at all positions. This plot shows what quality scores are most common.
+ **Per base sequence content**: plots the proportion of each base position over all of the reads. Typically, we expect to see each base roughly 25% of the time at each position, but this often fails at the beginning or end of the read due to quality or adapter content.
+ **Per sequence GC content**: a density plot of average GC content in each of the reads.  
+ **Per base N content**: the percent of times that 'N' occurs at a position in all reads. If there is an increase at a particular position, this might indicate that something went wrong during sequencing.  
+ **Sequence Length Distribution**: the distribution of sequence lengths of all reads in the file. If the data is raw, there is often on sharp peak, however if the reads have been trimmed, there may be a distribution of shorter lengths.
+ **Sequence Duplication Levels**: A distribution of duplicated sequences. In sequencing, we expect most reads to only occur once. If some sequences are occurring more than once, it might indicate enrichment bias (e.g. from PCR). If the samples are high coverage (or RNA-seq or amplicon), this might not be true.  
+ **Overrepresented sequences**: A list of sequences that occur more frequently than would be expected by chance.
+ **Adapter Content**: a graph indicating where adapater sequences occur in the reads.

## Working with the FastQC text output

Now that we've looked at our HTML reports to get a feel for the data,
let's look more closely at the other output files. Go back to the tab
in your terminal program that is connected to your AWS instance
(the tab lab will start with `dcuser@ip`) and make sure you're in
our results subdirectory.   

~~~
$ cd ~/RNAseq070319/fastqc/
$ ls
~~~
{: .bash}

~~~
sub06_fastqc.html  sub07_fastqc.zip   sub21_fastqc.html  sub23_fastqc.zip
sub06_fastqc.zip   sub08_fastqc.html  sub21_fastqc.zip   sub24_fastqc.html
sub07_fastqc.html  sub08_fastqc.zip   sub23_fastqc.html  sub24_fastqc.zip
~~~


Our `.zip` files are compressed files. They each contain multiple
different types of output files for a single input FASTQ file. To
view the contents of a `.zip` file, we can use the program `unzip`
to decompress these files. Let's try doing them all at once using a
wildcard.

~~~
$ unzip *.zip
~~~
{: .bash}

~~~
Archive:  SRR2584863_1_fastqc.zip
caution: filename not matched:  sub07_fastqc.zip
caution: filename not matched:  sub08_fastqc.zip
caution: filename not matched:  sub21_fastqc.zip
caution: filename not matched:  sub22_fastqc.zip
caution: filename not matched:  sub24_fastqc.zip
~~~


This didn't work. We unzipped the first file and then got a warning
message for each of the other `.zip` files. This is because `unzip`
expects to get only one zip file as input. We could go through and
unzip each file one at a time, but this is very time consuming and
error-prone. Someday you may have 500 files to unzip!

A more efficient way is to use a `for` loop like we learned in the Shell Genomics lesson to iterate through all of
our `.zip` files. Let's see what that looks like and then we'll
discuss what we're doing with each line of our loop.

~~~
$ for filename in *.zip
> do
> unzip $filename
> done
~~~
{: .bash}

In this example, the input is six filenames (one filename for each of our `.zip` files).
Each time the loop iterates, it will assign a file name to the variable `filename`
and run the `unzip` command.
The first time through the loop,
`$filename` is `sub06_fastqc.zip`.
The interpreter runs the command `unzip` on `sub06_fastqc.zip`.
For the second iteration, `$filename` becomes
`Ssub07_fastqc.zip`. This time, the shell runs `unzip` on `sun07_fastqc.zip`.
It then repeats this process for the four other `.zip` files in our directory.


When we run our `for` loop, you will see output that starts like this:

~~~
Archive:  sub06_fastqc.zip
   creating: sub06_fastqc/
   creating: sub06_fastqc/Icons/
   creating: sub06_2_fastqc/Images/
  inflating: sub06_2_fastqc/Icons/fastqc_icon.png  
  inflating: sub06_2_fastqc/Icons/warning.png  
  inflating: sub06_2_fastqc/Icons/error.png  
  inflating: sub06_2_fastqc/Icons/tick.png  
  inflating: sub06_2_fastqc/summary.txt  
  inflating: sub06_2_fastqc/Images/per_base_quality.png  
  inflating: sub06_2_fastqc/Images/per_tile_quality.png  
  inflating: sub06_2_fastqc/Images/per_sequence_quality.png  
  inflating: sub06_2_fastqc/Images/per_base_sequence_content.png  
  inflating: sub06_2_fastqc/Images/per_sequence_gc_content.png  
  inflating: sub06_fastqc/Images/per_base_n_content.png  
  inflating: sub06_fastqc/Images/sequence_length_distribution.png  
  inflating: sub06_fastqc/Images/duplication_levels.png  
  inflating: sub06_fastqc/Images/adapter_content.png  
  inflating: sub06_fastqc/fastqc_report.html  
  inflating: sub06_fastqc/fastqc_data.txt  
  inflating: sub06_fastqc/fastqc.fo  
~~~


The `unzip` program is decompressing the `.zip` files and creating
a new directory (with subdirectories) for each of our samples, to
store all of the different output that is produced by FastQC. There
are a lot of files here. The one we're going to focus on is the
`summary.txt` file.

If you list the files in our directory now you will see:

~~~
SRR2584863_1_fastqc       SRR2584866_1_fastqc       SRR2589044_1_fastqc
SRR2584863_1_fastqc.html  SRR2584866_1_fastqc.html  SRR2589044_1_fastqc.html
SRR2584863_1_fastqc.zip   SRR2584866_1_fastqc.zip   SRR2589044_1_fastqc.zip
SRR2584863_2_fastqc       SRR2584866_2_fastqc       SRR2589044_2_fastqc
SRR2584863_2_fastqc.html  SRR2584866_2_fastqc.html  SRR2589044_2_fastqc.html
SRR2584863_2_fastqc.zip   SRR2584866_2_fastqc.zip   SRR2589044_2_fastqc.zip
~~~
{:. output}

The `.html` files and the uncompressed `.zip` files are still present,
but now we also have a new directory for each of our samples. We can
see for sure that it's a directory if we use the `-F` flag for `ls`.

~~~
$ ls -F
~~~
{: .bash}

~~~
SRR2584863_1_fastqc/      SRR2584866_1_fastqc/      SRR2589044_1_fastqc/
SRR2584863_1_fastqc.html  SRR2584866_1_fastqc.html  SRR2589044_1_fastqc.html
SRR2584863_1_fastqc.zip   SRR2584866_1_fastqc.zip   SRR2589044_1_fastqc.zip
SRR2584863_2_fastqc/      SRR2584866_2_fastqc/      SRR2589044_2_fastqc/
SRR2584863_2_fastqc.html  SRR2584866_2_fastqc.html  SRR2589044_2_fastqc.html
SRR2584863_2_fastqc.zip   SRR2584866_2_fastqc.zip   SRR2589044_2_fastqc.zip
~~~


Let's see what files are present within one of these output directories.

~~~
$ ls -F SRR2584863_1_fastqc/
~~~
{: .bash}

~~~
fastqc_data.txt  fastqc.fo  fastqc_report.html	Icons/	Images/  summary.txt
~~~


Use `less` to preview the `summary.txt` file for this sample.

~~~
$ less SRR2584863_1_fastqc/summary.txt
~~~
{: .bash}

~~~
PASS    Basic Statistics        SRR2584863_1.fastq
PASS    Per base sequence quality       SRR2584863_1.fastq
PASS    Per tile sequence quality       SRR2584863_1.fastq
PASS    Per sequence quality scores     SRR2584863_1.fastq
WARN    Per base sequence content       SRR2584863_1.fastq
WARN    Per sequence GC content SRR2584863_1.fastq
PASS    Per base N content      SRR2584863_1.fastq
PASS    Sequence Length Distribution    SRR2584863_1.fastq
PASS    Sequence Duplication Levels     SRR2584863_1.fastq
PASS    Overrepresented sequences       SRR2584863_1.fastq
WARN    Adapter Content SRR2584863_1.fastq
~~~


The summary file gives us a list of tests that FastQC ran, and tells
us whether this sample passed, failed, or is borderline (`WARN`). Remember to quit from `less` you enter `q`.

## Documenting Our Work

We can make a record of the results we obtained for all our samples
by concatenating all of our `summary.txt` files into a single file
using the `cat` command. We'll call this `full_report.txt` and move
it to `~/dc_workshop/docs`.

~~~
$ cat */summary.txt > ~/dc_workshop/docs/fastqc_summaries.txt
~~~
{: .bash}

> ## Exercise
>
> Which samples failed at least one of FastQC's quality tests? What
> test(s) did those samples fail?
>
>> ## Solution
>>
>> We can get the list of all failed tests using `grep`.
>>
>> ~~~
>> $ cd ~/dc_workshop/docs
>> $ grep FAIL fastqc_summaries.txt
>> ~~~
>> {: .bash}
>>
>> ~~~
>> FAIL    Per base sequence quality       SRR2584863_2.fastq.gz
>> FAIL    Per tile sequence quality       SRR2584863_2.fastq.gz
>> FAIL    Per base sequence content       SRR2584863_2.fastq.gz
>> FAIL    Per base sequence quality       SRR2584866_1.fastq.gz
>> FAIL    Per base sequence content       SRR2584866_1.fastq.gz
>> FAIL    Adapter Content SRR2584866_1.fastq.gz
>> FAIL    Adapter Content SRR2584866_2.fastq.gz
>> FAIL    Adapter Content SRR2589044_1.fastq.gz
>> FAIL    Per base sequence quality       SRR2589044_2.fastq.gz
>> FAIL    Per tile sequence quality       SRR2589044_2.fastq.gz
>> FAIL    Per base sequence content       SRR2589044_2.fastq.gz
>> FAIL    Adapter Content SRR2589044_2.fastq.gz
>> ~~~
>>
>>
> {: .solution}
{: .challenge}


# Other notes  -- Optional

> ## Quality Encodings Vary
>
> Although we've used a particular quality encoding system to demonstrate interpretation of
> read quality, different sequencing machines use different encoding systems. This means that,
> depending on which sequencer you use to generate your data, a `#` may not be an indicator of
> a poor quality base call.
>
> This mainly relates to older Solexa/Illumina data,
> but it's essential that you know which sequencing platform was
> used to generate your data, so that you can tell your quality control program which encoding
> to use. If you choose the wrong encoding, you run the risk of throwing away good reads or
> (even worse) not throwing away bad reads!
{: .callout}


> ## Same Symbols, Different Meanings
>
> Here we see `>` being used a shell prompt, whereas `>` is also
> used to redirect output.
> Similarly, `$` is used as a shell prompt, but, as we saw earlier,
> it is also used to ask the shell to get the value of a variable.
>
> If the *shell* prints `>` or `$` then it expects you to type something,
> and the symbol is a prompt.
>
> If *you* type `>` or `$` yourself, it is an instruction from you that
> the shell to redirect output or get the value of a variable.
{: .callout}








# Trimming and filtering

Before we will do the alignment we need to remove sequences of low quality and sequences that are to short (below 25 bases).
Also in this case we will trim down long sequences to 100 bases, quality of the Ion-torrent reads drops the further it gets.
When making use of illumina reads this is not as much of a problem and 3'-trimming would then be a waste of data.

To start off make a directory trimmed for the output and then back to the rawReads directory.

~~~
$ cd ~/RNAseq070319/
$ mkdir trimmed
$ cd ~/RNAseq070319/rawReads/
~~~
{: .bash}



The trimming and quality filtering will be done with trimmomatic.
In the programm the following arguments can be used.

| step   | meaning |
| ------- | ---------- |
| `SE` or `PE` | Reads are single end or paired end. |
| `ILLUMINACLIP` | Perform adapter removal |
| `SLIDINGWINDOW` | Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. |
| `LEADING`  | Cut bases off the start of a read, if below a threshold quality.  |
|  `TRAILING` |  Cut bases off the end of a read, if below a threshold quality. |
| `CROP`  |  Cut the read to a specified length. |
|  `HEADCROP` |  Cut the specified number of bases from the start of the read. |
| `MINLEN`  |  Drop an entire read if it is below a specified length. |
|  `TOPHRED33` | Convert quality scores to Phred-33.  |
|  `TOPHRED64` |  Convert quality scores to Phred-64. |


To run this on a single sample it looks something like this
~~~
trimmomatic SE -phred33 -threads 2 sub06.fastq ../trimmed/sub06_qc.fq ILLUMINACLIP:../adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
~~~
{: .bash}


Of cource we don't want to do this for all the reads seperately so lets create a loop through all the fastq files.

When doing the fastqc only input files needed to be specified. In this case both the input and a matching output filenames need to be given.
this can be done with the help of 'basename'


~~~
$ for fastq in *.fastq
do
 echo inputfile $fastq
 base="$(basename $fastq .fastq)"_qc.fq
 echo outputfile $base
 echo
done
~~~
{: .bash}

This be be producing the following list

~~~
inputfile sub06.fastq
outputfile sub06_qc.fq

inputfile sub07.fastq
outputfile sub07_qc.fq

inputfile sub08.fastq
outputfile sub08_qc.fq

inputfile sub21.fastq
outputfile sub21_qc.fq

inputfile sub23.fastq
outputfile sub23_qc.fq

inputfile sub24.fastq
outputfile sub24_qc.fq
~~~


Next we can start writing the trimmomatic loop.
Again starting with a dry run with echo.

~~~
$ for fastq in *.fastq
do
  outputFile="$(basename $fastq .fastq)"_qc.fq
  echo trimmomatic SE -phred33 -threads 2 $fastq ../trimmed/$outputFile ILLUMINACLIP:../adapters.fasta:2:30:10 LEADING:3 TRAILING:3     SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
done
~~~
{: .bash}

should be producing something like this

~~~
trimmomatic SE -phred33 -threads 2 sub06.fastq ../trimmed/sub06_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
trimmomatic SE -phred33 -threads 2 sub07.fastq ../trimmed/sub07_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
trimmomatic SE -phred33 -threads 2 sub08.fastq ../trimmed/sub08_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
trimmomatic SE -phred33 -threads 2 sub21.fastq ../trimmed/sub21_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
trimmomatic SE -phred33 -threads 2 sub23.fastq ../trimmed/sub23_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
trimmomatic SE -phred33 -threads 2 sub24.fastq ../trimmed/sub24_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
~~~


If it al seems ok rerun with out 'echo'

~~~
$ for fastq in *.fastq
do
    trimmomatic SE -phred33 -threads 2 $fastq ../trimmed/"$(basename "$fastq" .fastq)"_qc.fq ILLUMINACLIP:../adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
done
~~~
{: .bash}

The following should appear:

~~~
TrimmomaticSE: Started with arguments:
 -phred33 -threads 2 sub06.fastq ../trimmed/sub06_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
…
Input Reads: 1000000 Surviving: 887553 (88.76%) Dropped: 112447 (11.24%)
TrimmomaticSE: Completed successfully
~~~


It's possible to scroll up to check if the percentage of surviving & dropped is within the same range in all of the samples.
