---
title: "Exploring RNA-Seq results interactively with ideal"
teaching: 45
exercises: 0
questions:
- "How can I analyze the counts obtained from the DESeq2 analysis?"
objectives:
- "Understand the model used by DESeq2 to compute differential expression"
- "Understand the different units for counts (TPM etc.)"
- "Correlate samples based on counts."
- "Understand the experimental factor from which to compute differential expression"
keypoints:
- "Raw counts, normalized counts, TPM"
- "p-value and false discovery rate"
- "log fold change"
- "MA plot, Volcanot plot"
---
# Download files
If you haven't managed to get the `counts.final.tsv` file, no worry! Just download it here: https://doi.org/10.5281/zenodo.2583083    

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
![counts](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/data-scale.png)      

Try different units:
- Counts (raw)
- Counts (normalized)
- Log10 (pseudocount of 1 added)

__Question:__ why is the difference between raw and normalized counts?

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
To calculate a log fold change, select the numerator (_drought_) and denominator (_control_).    
![xp](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/experimental-factor.png)   


__Question__: what would mean a positive log fold change?

### Select the False Discovery Rate
For each gene, we are testing the _null hypothesis_, that is: "the gene of interest is not differentially regulated between my experimental conditions (i.e. drought and control)". Since we are testing this null hypothesis thousands of times, we will overestimate the number of false positives (type I error) and thus overestimate the number of differentially expressed genes if we don't make any corrections.  
For instance, you have 10,000 genes being tested and you choose a p-value cutoff of 0.05 (5%). That means that you can say that a gene is differentially regulated with a type I error of 5% (in other words, a confidence of 95%).  
But since you have 10,000 genes, that means that 10,000 genes x 5% = 500 genes will be false positives (called differential while they are not).   
To reduce the number of false positives, we are going to control (but not eliminate) the number of false positives using a False Discovery Rate (corrected p-value).

To make a long story short, you can remember this:
> if you repeat a test enough times, you’re going to find an effect…but that effect may not actually exist.  

You can change the FDR value here: ![fdr](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/fdr.png)   

Click on "Extract the results!" ![extract the results](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/extract-results-button.png) when ready.

### Result table
Explanation of the result table:

| column | description | example |
|--------|-------------|---------|
| unnamed | gene name | ATCG00470|
| baseMean | average of the normalized count values divided by the size factors | 3038.11149196245 |
| log2FoldChange | log2 scaled fold change between drought and control (log2FC = log2(drought) - log2(control)) |3.66236695846517	|
| lfcSE | Standard error estimate for the log2 fold change estimate. |0.290498413254083	|
| stat | Statistical value associated (origin?) |-10.9229847880748	|
| pvalue | the associated p-value |8.95051429777092e-28	|
| padj | the p-value corrected using the false discovery rate |1.64403046621456e-24|
| symbol | gene symbol (if an annotation is available) |-|


### Diagnostic plots
The first interesting plot is the p-value histogram. It looks really good because we have a high peak of genes with a highly significant p-value (left side). And the non-significant differential genes have p-values scattered between 0 and 1. For a deeper explanation see [this blog post](http://varianceexplained.org/statistics/interpreting-pvalue-histogram/).    
![pvalue distribution](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/pvalues_distribution.png)     

The second interesting diagnostic plot is the stratified p-value histogram. You can see that expression level has a consequence for calling differential genes.
![pvalue stratified distribution](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/pvalues_distribution_stratified.png)     

The last but not least diagnostic plot that is interesting is the log2 fold change histogram. Since most genes should not be differential between your conditions, the log2 fold change should be centered around zero and looks like a Guaussian law.
![pvalue stratified distribution](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/log2fc_histogram.png)     

## Summary Plots
Let's now look at the __Summary Plots__ section: ![summary plots](https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/summary-plots.png)     
We'll explain it as we go through:  

### MA plot: mean of normalized counts _versus_ log fold change
From the [DESeq2 documentation](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#ma-plot):
> An MA-plot (Dudoit et al. 2002) provides a useful overview for the distribution of the estimated coefficients in the model, e.g. the comparisons of interest, across all genes. On the y-axis, the “M” stands for “minus” – subtraction of log values is equivalent to the log of the ratio – and on the x-axis, the “A” stands for “average”. You may hear this plot also referred to as a mean-difference plot, or a Bland-Altman plot.

__To do list__:
- Zoom in a section of the MA plot + select a gene to display the boxplot
- Volcano plot: log2 Fold Change _versus_ -log10(pvalue)  
- Heatmap: built from the region you have selected in the MA plot

__Exercise__: in the MA plot, zoom on a region with nonDE genes (black). Look at the heatmap. Redo this with a region with DE genes (red). What can you say about the clustered heatmap?
