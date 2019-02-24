---
title: "File preparation for online analysis with ideal"
teaching: 15
exercises: 0
questions:
- "How can I analyse the counts obtained from the DESeq2 analysis?"
objectives:
- "Understand what is a Shiny application."
- "Load the counts and experimental design files into ideal."
- "Understand the dds class from DESeq."
keypoints:
- "Use the `Data Setup` panel to load the `counts` and `design` files."
---
## Outline
1. Description of the `counts.tsv` and the `design.tsv` files
2. Parsing (= clean-up) of the `counts.txt` and of the `design.tsv` files.
2. Generation a `DESeqDataSet` data object in the __ideal__ online Shiny application.

### Description of the counts file
The __featureCounts__ program used read alignment and annotation information to compute counting values for each of the annotation gene in the genome. That information is stored in the `counts.tsv` file.  
This file is a tabulated file (columns are separated by a tabulation).     

#### Removing the first line
We can see that the first line is a comment line:  
` # Program:featureCounts v1.6.3; Command:"featureCounts" "-t" "exon" "-a" "/home/tbliek/RNAseqWorkshop/genome/annotation.all_transcripts.exon_features.ath.gff3" "-o" "counts.txt" "sub06_qc.bam" "sub07_qc.bam" "sub08_qc.bam" "sub21_qc.bam" "sub23_qc.bam" "sub24_qc.bam"`   
We need to remove that line because the `ideal` application won't be able to process it.  
For that, we are going to use the __awk__ language that can be run from the Shell directly.

__Removing the first line with awk:__ in the Shell, type the following command  
 `awk '{if (NR!=1){print}} counts.txt > counts.parsed.txt`  
 Check that it did the trick `head counts.parsed.txt`. The comment line should be removed.

#### Removing unnecessary columns (chr,strand,etc.)
Now our file looks like:   

| Geneid    | Chr                           | Start                         | End                           | Strand      | Length | sub06_qc.bam | sub07_qc.bam | sub08_qc.bam | sub21_qc.bam | sub23_qc.bam | sub24_qc.bam |
|-----------|-------------------------------|-------------------------------|-------------------------------|-------------|--------|--------------|--------------|--------------|--------------|--------------|--------------|
| AT1G01010 | Chr1;Chr1;Chr1;Chr1;Chr1;Chr1 | 3631;3996;4486;4706;5174;5439 | 3913;4276;4605;5095;5326;5899 | +;+;+;+;+;+ | 1688   | 0            | 0            | 6            | 6            | 3            | 10           |
| AT1G01020 |  Chr1;Chr1;.... n times       |                               |                               |             |        | 2            | 4            | 3            | 1            | 2            | 1            |
| AT1G03987 | ...                           |                               |                               |             |        |              |              |              |              |              |              |
| ...       | ...                           | ...                           | ...                           | ...         | ...    | ...          | ...          | ...          | ...          | ...          | ...          |


We need to remove a few columns so that the file looks like:  

| Geneid    | sub06_qc.bam | sub07_qc.bam | sub08_qc.bam | sub21_qc.bam | sub23_qc.bam | sub24_qc.bam |
|-----------|--------------|--------------|--------------|--------------|--------------|--------------|
| AT1G01010 | 0            | 0            | 6            | 6            | 3            | 10           |
| AT1G01020 | 2            | 4            | 3            | 1            | 2            | 1            |
| AT1G03987 |              |              |              |              |              |              |
| ...       | ...          | ...          | ...          | ...          | ...          | ...          |

__Removing a few columns with cut:__ in the Shell, type `cut -f2,3,4,5,6 --complement counts.parsed.txt > counts.final.txt`

### Design file
