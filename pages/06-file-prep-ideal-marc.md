---
title: "File preparation for online analysis with ideal"
teaching: 15
exercises: 0
questions:
- "How can I analyse the counts obtained from the DESeq2 analysis?"
objectives:
- "Understand what is the content of the counts and design files."
- "Be able to transform the counts file into a desired format using the Shell."
keypoints:
- "sed, awk and cut commands"
- "counts and design files for ideal."
---

# Outline
1. The __ideal__ application for rnaseq analysis.
2. Description of the `counts.tsv` file.
3. Parsing (= transformation) of the `counts.tsv` file.
4. Description of the `design.tsv` file.

# What is ideal?
__ideal__ stands for: "Interactive Differential Expression AnaLysis". It is a Shiny online application that runs R code without having you to type a single line of code.       
![](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/ideal_presentation.png)  

We will use it as it will help us to explore the output of our trimming+alignment+count steps.  
We will need two files for ideal: a `counts.tsv` file and a `design.tsv` file.

# Input files for ideal
## Description of the counts file
The __featureCounts__ program used read alignment and annotation information to compute counting values for each annotated gene. That information is stored in the `counts.tsv` file.    
This file is a tabulated file (columns are separated by a tabulation).     

#### Removing the first line
We can see that the first line is a comment line:  
` # Program:featureCounts v1.6.3; Command:"featureCounts" "-t" "exon" "-a" "/home/tbliek/RNAseqWorkshop/genome/annotation.all_transcripts.exon_features.ath.gff3" "-o" "counts.txt" "sub06_qc.bam" "sub07_qc.bam" "sub08_qc.bam" "sub21_qc.bam" "sub23_qc.bam" "sub24_qc.bam"`   
We need to remove that line because the `ideal` application won't be able to process it.  
For that, we are going to use the __awk__ language that can be run from the Shell directly.

__Removing the first line with awk:__ in the Shell, type the following command:    
 `awk '{if (NR!=1){print}}' counts.tsv > counts.parsed.tsv`

 Check that it did the trick: in the Shell, type:  
 `head counts.parsed.tsv`  
 The comment line should be removed.

### Removing unnecessary columns (chr,strand,etc.)
Now our file looks like:   

| Geneid    | Chr                           | Start                         | End                           | Strand      | Length | sub06_qc.bam | sub07_qc.bam | sub08_qc.bam | sub21_qc.bam | sub23_qc.bam | sub24_qc.bam |
|-----------|-------------------------------|-------------------------------|-------------------------------|-------------|--------|--------------|--------------|--------------|--------------|--------------|--------------|
| AT1G01010 | Chr1 | 3631 |5899 | +| 1688   | 0            | 0            | 6            | 6            | 3            | 10           |
| AT1G01020 |  Chr1       |                               |                               |             |        | 2            | 4            | 3            | 1            | 2            | 1            |
| AT1G03987 | ...                           |                               |                               |             |        |              |              |              |              |              |              |
| ...       | ...                           | ...                           | ...                           | ...         | ...    | ...          | ...          | ...          | ...          | ...          | ...          |


We need to remove a few columns so that the file looks like:  

| Geneid    | sub06_qc.bam | sub07_qc.bam | sub08_qc.bam | sub21_qc.bam | sub23_qc.bam | sub24_qc.bam |
|-----------|--------------|--------------|--------------|--------------|--------------|--------------|
| AT1G01010 | 0            | 0            | 6            | 6            | 3            | 10           |
| AT1G01020 | 2            | 4            | 3            | 1            | 2            | 1            |
| AT1G03987 |              |              |              |              |              |              |
| ...       | ...          | ...          | ...          | ...          | ...          | ...          |

__Removing a few columns with cut:__ in the Shell, type `cut -f2,3,4,5,6 --complement counts.parsed.tsv > counts.parsed.cut.tsv`

| command | description |
| -------- | -------- |
| cut | program name for cutting sections from each line (cut==keep in output)|
| -f | specify the fields (columns) you want to cut |
|--complement | since we want to remove columns 2 to 6, we use the complement option |  

### Renaming samples to remove the "\_qc.bam"
To have nice and short sample names ("sub06"), we will search and replace the "\_qc.bam" inside the entire file.    
One nice simple tool to do that is _sed_ and the command is:  
`sed 's/_qc.bam//g' counts.parsed.cut.tsv > counts.parsed.cut.renamed.tsv`  

| command | description |
| -------- | -------- |
| sed | program name that stands for __s__tream __ed__itor__ |
| s | s for substitute |
| _qc_bam | the string we want to substitute |
| // | this means that we want to replace it with an empty string |
| g | this means that we want to do it globally for the whole file |

A complete but understandable guide for sed is available [here](http://www.grymoire.com/Unix/Sed.html).

Check that it worked: `head counts.parsed.cut.renamed.tsv`

### Rename the file
Finally, we are going to rename that file because the name is too long. In the Shell, type:  
`mv counts.parsed.renamed.tsv counts.final.tsv`

## Design file
The design file tells the relationship between samples and experimental conditions.

| SampleName | condition |
|------------|-----------|
| sub06 | control |
| sub07 | control |
| sub08 | control |
| sub21 | drought |
| sub23 | drought |
| sub24 | drought |  

This file is in the right format for __ideal__. No need to do anything on this file.

# Navigate to a public instance of the ideal application

## Public online instances
:de: at the Institute of Medical Biostatistics, Epidemiology and Informatics (IMBEI) in Mainz: http://shiny.imbei.uni-mainz.de:3838/ideal/

:computer: using the installed instance of __ideal__ on the genseq cluster: http://genseq-h0.science.uva.nl/shiny/ideal/

__Warning__:
<p style="border:3px; border-style:solid; border-color:#FF0000; padding: 1em;">Warning: it might be that the analysis takes too much time on both public instances due to everybody trying to access it at the same time.<br>If this is the case, you can try to start RStudio and install it locally.</p>

## Local installation in RStudio

### Installation
Run these commands in RStudio:  
```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ideal")
```

### Load the library and load files
- First, load the library by typing `library(ideal)` in the RStudio interactive console.   
- Then, start ideal by typing `ideal` in the RStudio console.

You can then follow the lesson as the rest of the class.
