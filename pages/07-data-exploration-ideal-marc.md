---
title: "Exploring RNA-Seq results interactively with ideal"
teaching: 45
exercises: 0
questions:
- "How can I analyze the counts obtained from the DESeq2 analysis?"
objectives:
- "Understand the model used by DESeq2 to compute differential expression"
- "Understand the different units for counts (TPM etc.)
- "Correlate samples based on counts."
- "Understand the experimental factor from which to compute differential expression"
- ""
keypoints:
- "fi
tted "
- "Raw counts, normalized counts, TPM "
---

# Data exploration and analysis using the ideal
For each panel, you can select the following button to get help: <br/>
![Help button](../images/help-button.png)  

## Data setup
We will first create a dds object ("DESeq object") that will contain all necessary information for further processing.   
 Select the __Data Setup__ panel: ![Data setup](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/data_setup.png)   

 Import the `counts.final.tsv` and the `design.tsv` files. Upon upload, you should have:  
 ![import](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/import-tables.png)     

You can take a look at how __ideal__ imported your counts.
![preview](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/count-design-preview.png)      


## Counts overview
<- add picture of Counts overview panel ->  
### Data scale
__Data scale in the table__: here you can change the units of the counts. Try different ones:
- Counts (raw)
- Counts (normalized)
- Log10 (pseudocount of 1 added)

__Question:__ what is the unit difference between raw and normalized counts? Any idea why?

### Basic summary
In the __Basic summary__ subsection, you can see a summary of the count data.   
Try to change the __"Threshold on the row sums of the counts"__ value from 0 to some other bigger number.

### Sample correlations
Here you can visualise the correlation between your samples.

## Extract results
<- add picture of Extract Results panel ->

### Select the experimental factor
Since we have only one experimental factor that is quite easy: select the "condition" factor that contains the _drought_ and _control_ levels.

__Question:__ what should you   
Then select the numerator (_drought_) and denominator (control).

### Result table

Explanation of the results:
> The first column, baseMean, is a just the average of the normalized count values, divided by the size factors, taken over all samples in the DESeqDataSet. The remaining four columns refer to a specific contrast, namely the comparison of the trt level over the untrt level for the factor variable dex. We will find out below how to obtain other contrasts.

### Diagnostic plots
- p-value histogram
- stratified p-value histogram
- log2 fold change histogram

## Summary Plots
- MA plot: mean of normalized counts _versus_ log fold change
- Zoom in a section of the MA plot + select a gene to display the boxplot
- Volcano plot: log2 Fold Change _versus_ -log10(pvalue)  
- Heatmap: built from the region you have selected in the MA plot

__Exercise__: in the MA plot, zoom on a region with nonDE genes (black). Look at the heatmap. Redo this with a region with DE genes (red). What can you say about the clustered heatmap?

## Gene of interest
Here you can look for your favorite gene of interest.      
