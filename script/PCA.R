## Usage:
# Made customized PCA plot functions by revising the DESeq2's buildin PCA function, see "PCA_functions.R".
# Here, run the functions to geneate PCA plots.


## setup
library(tidyverse)
library(DESeq2)
source("PCA_functions.R")

## provided data, meta to create DESeqDataSet object
# create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~sampletype)

# data transformation 
# vst: do both size factor normalization+transformation
vsd <- vst(dds)

## use the PCA functions
# function input:
# object=vsd
# intgroup="sampletype"
# PC_1=1
# PC_2=2
# ntop=500
# PCA_label <- meta$PCA_label

# check the plots
p1 <- PCAdotlab(object=vsd, intgroup="sampletype", PC_1=1, PC_2=2, ntop=500, PCA_label=PCA_label)
p2 <- PCAdot(object=vsd, intgroup="sampletype", PC_1=1, PC_2=2, ntop=500)
p3 <- PCAlab(object=vsd, intgroup="sampletype", PC_1=1, PC_2=2, ntop=500, PCA_label=PCA_label)

print(p1)
print(p2)
print(p3)

# save files
file_dir <- "../report/PCA/"
ggsave(p1,
       filename = paste0(file_dir, "PCA", "_12_", intgroup, ".pdf"),
       height = 7, width = 7,units = "in")
