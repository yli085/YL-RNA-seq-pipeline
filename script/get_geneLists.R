## script to get gene lists of all DE comparisons for the metascape analysis
# input: limma output "gene_ls", and cutoff conditions
# output: excel file with gene lists, one tab per comparison.

library(tidyverse)
library(writexl)

# makePaddedDataFrame() function ----------------------------------------------
# the function for making padded dataframe
na.pad <- function(x,len){
  x[1:len]
}

makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}


# input  --------------------------------------------------------------
## selection 2 cutoff conditions
gene_ls <- readRDS("../report/report_wo_aging/DE_analysis/DE_sum/gene_ls.RDS")
comparisons <- names(gene_ls)
top_ls <- list()

## two cutoff conditions for selecting DEGs or top genes for the DE comparisons
# condition 1 for C-WT_AG, comparison 3
padj.cutoff=0.05
lfc.cutoff=0
# condition 2 for others
rpval.cutoff = 0.05


# make list for functional analysis ---------------------------------------

for (i in seq_along(comparisons)) {
  
  # this fl is used for new list naming
  fl = comparisons[i]
  gene_df <- gene_ls[[i]]
  
  if (i==3) {
    # use condition 1
    up_df <- gene_df %>% dplyr::filter(logFC > lfc.cutoff, adj.P.Val < padj.cutoff)
    down_df <- gene_df %>% dplyr::filter(logFC < -lfc.cutoff, adj.P.Val < padj.cutoff)
    
  } else {
    # use condition 2
    up_df <- gene_df %>% dplyr::filter(logFC > 0, P.Value < rpval.cutoff)
    down_df <- gene_ls[[fl]] %>% dplyr::filter(logFC < 0, P.Value < rpval.cutoff)
  }
  
  # make the gene lists a dataframe, though uneven row numbers
  up <- up_df[[1]]
  down <- down_df[[1]]
  
  top_ls[[fl]] <- makePaddedDataFrame(list(up=up,down=down))
  
}

# save files
file_dir <- "../report/report_wo_aging/functional/"
filenm <- paste0(file_dir, "metacape_genelist.xlsx")

write_xlsx(top_ls, path=filenm)
# saveRDS(top_ls, paste0(file_dir, "top_ls.RDS"))
