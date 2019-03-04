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
We will first create a _dds_ object ("DESeqDataSet object") that will contain all necessary information for further processing.   
 Select the __Data Setup__ panel: ![Data setup](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/data_setup.png)   

### Step 1: import data tables
 Import the `counts.final.tsv` and the `design.tsv` files. Upon upload, you should have:  
 ![import](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/import-tables.png)       

You can take a look at how __ideal__ imported your counts.
![preview](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/count-design-preview.png)      

### Step 2: generate the dds object
From the authors of DESeq2:   
> The object class used by the DESeq2 package to store the read counts and the intermediate estimated quantities during statistical analysis is the DESeqDataSet, which will usually be represented in the code here as an object dds.

Generate the _dds_ object by selecting the _condition_ for the experimental design. This will indicate that we want to compare the levels (_control_ and _drought_) of the _condition_ column from the `design.tsv`.
![dds](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/generate-dds.png)    

Once this is done, you'll see a short summary of the _DESeqDataSet_ object that you've generated.  Skip the two optional steps and move on to Step 3.

### Step 3: run DESeq2
Finally, we compute a standard differential expression analysis by running the `DESeq` function. If you want to learn more about the theory behind it, have a look at the manual of DESeq2 [here](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory).   

Click on `Run DESeq`:  
![deseq](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/run-deseq.png)    
You can see a short overview of the analysis. In particular, how many genes have non-zero values (= are detected).   
Moreover, you can see the number of genes that have a log fold change (LFC) superior or inferior to zero (up-regulated or down-regulated).


## Counts overview
Let's move on and have a quick look at the counts. Select the __Counts Overview__ panel: ![counts-overview](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/counts-overview.png)    

### Data scale
__Data scale in the table__: here you can change the units of the counts.
![counts](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/counts-overview.png)      

Try different ones:
- Counts (raw)
- Counts (normalized)
- Log10 (pseudocount of 1 added)

__Question:__ what is the unit difference between raw and normalized counts? Any idea why?

### Basic summary
In the __Basic summary__ subsection, you can see a summary of the count data.     
![basic](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/basic-summary.png)   

 Try to change the __"Threshold on the row sums of the counts"__ value from 0 to some other bigger number.

 __Question__: select a few values from 10 to 10,000. What happens to the number of detected genes? What does it tell you about the underlying distribution of genes in terms of expression values?  

 You can also filter on the column sums (on samples and not on genes).
 The _dds_ object can then be updated to keep filtered genes.

### Sample correlations
Here you can visualise the correlation between your samples.

## Extract results
Next, we want to extract the differential expression analysis results. We can do so by going to the __Extract Results__ panel: ![extract](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/extract-results.png)   

### Select the experimental factor
Since we have only one experimental factor that is quite easy: select the "condition" factor that contains the _drought_ and _control_ levels.

![extract](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/extract-results.png)   

To calculate a log fold change, select the numerator (_drought_) and denominator (_control_).

__Question__: what would mean a positive log fold change?

### Select the False Discovery Rate
For each gene, we are testing the _null hypothesis_, that is: "the gene of interest is not differentially regulated between my experimental conditions (i.e. drought and control)". Since we are testing this null hypothesis thousands of times, we will overestimate the number of false positives (type I error) and thus overestimate the number of differentially expressed genes if we don't make any corrections.  
For instance, you have 10,000 genes being tested and you choose a p-value cutoff of 0.05 (5%). That means that you can say that a gene is differentially regulated with a type I error of 5% (in other words, a confidence of 95%).  
But since you have 10,000 genes, that means that 10,000 genes x 5% = 500 genes will be false positives (called differential while they are not).   
To reduce the number of false positives, we are going to control (but not eliminate) the number of false positives using a False Discovery Rate (corrected p-value).

You can change the FDR value here: ![fdr](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/fdr.png)   


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
