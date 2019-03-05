########################
### library installation
########################
if (!"BiocManager" %in% installed.packages()){
  install.packages("BiocManager")
  library("BiocManager",quietly = F)
} else {
    library("BiocManager",quietly = F)
}

if (!"org.At.tair.db" %in% installed.packages()){
  BiocManager::install('org.At.tair.db',version = "3.8")
} else {
  library("org.At.tair.db",quietly = F)
}

if (!"GenomeInfoDb" %in% installed.packages()){
  BiocManager::install("GenomeInfoDb",version = "3.8")
  library("GenomeInfoDb",quietly = F)
} else {
  library("GenomeInfoDb",quietly = F)
}


if (!"ideal" %in% installed.packages()){
  BiocManager::install("ideal",version = "3.8")
  library("ideal",quietly = F)
} else {
  library("ideal",quietly = F)
}


########################
### data import
########################
mytable <- read.delim("https://raw.githubusercontent.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/gh-pages/files/counts.final.tsv",sep="\t")
mydesign <- read.delim("https://raw.githubusercontent.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/gh-pages/files/design.tsv",sep="\t")
mydesign

rownames(mydesign) <- mydesign$SampleName
rownames(mytable) <- mytable$Geneid
head(mytable)
mytable <- mytable[,-1]
head(mytable)
dds <- DESeq2::DESeqDataSetFromMatrix(mytable,mydesign,~1)
myanno <- pcaExplorer::get_annotation_orgdb(dds, "org.At.tair.db","TAIR")
ideal(dds,annotation_obj = myanno)
