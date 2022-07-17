## script to generate a heatmap for the top 20 most significantly increased and decreased genes (a total of 40)
# also save those genes in excel files and RDS

library(tidyverse)
library(writexl)
library(gplots)
library(ggplot2)

## get the top gene from one comparison
gene_ls <- readRDS("../report/report_wo_aging/DE_analysis/DE_sum/gene_ls.RDS")
comparisons <- names(gene_ls)

# gene name annotation
gene_id_df <- readRDS("/Users/YL/Documents/Nutstore/Nutstore/TB_NGS/RNAseq_method/genome_annotations/ensembl_ids_conversion.RDS")

# get the de result table
df_ls <- vector("list", length(comparisons))
names(df_ls) <- comparisons

igenes <-vector("list", length(comparisons))
names(igenes) <- comparisons

## get the de table for the igenes for each comparison group
for (i in seq_along(comparisons)) {
  
  de_df <- gene_ls[[i]]
  # select top 20 most significantly changed genes and combine them.
  de_up_df <- de_df %>% filter(logFC > 0) %>% head(n=20) %>% 
    mutate(Regulation = 'Up')
  de_down_df <- de_df %>% filter(logFC < 0) %>% head(n=20)%>% 
    mutate(Regulation = 'Down')
  
  top_df <- bind_rows(de_up_df, de_down_df) %>% 
    left_join(gene_id_df, by=c('gene'='gene_id')) %>% 
    select('gene', 'gene_name', 'entrez_description', 'Regulation', 'AveExpr',
           'logFC', 'P.Value', 'adj.P.Val')
  
  # rename the first few columns
  names(top_df)[1:3] <- c("Gene_ID", "Gene_symbol", "Gene_Description")
  
  df_ls[[i]] <- top_df
  
  # get the gene names
  igenes[[i]] <- as.character(top_df$Gene_symbol)

}

# save files
heatmap_dir <- "../report/report_wo_aging/DE_analysis/heatmap/"
# dir.create(heatmap_dir)
write_xlsx(df_ls, path=paste0(heatmap_dir, 'top_degs.xlsx'), col_names = TRUE, format_headers = TRUE)
saveRDS(igenes, paste0(heatmap_dir, 'top_degs.RDS'))

# heatmap -----------------------------------------------------------------
# read in logCPM
lcpm_df <- read.table("../report/report_wo_aging/gene_expression/logCPM.txt", header=T, stringsAsFactors = F)

# add gene name and description to the expression matrix
lcpm_df <- lcpm_df %>% left_join(gene_id_df, by=c("gene"="gene_id")) %>% 
  select(-c("gene", "ensembl_id", "mgi_id", "entrez_id", "entrez_description")) %>% 
  relocate(gene_name) %>% 
  rename(Gene=gene_name)

# make and save heatmaps
for (i in seq_along(comparisons)) {
  
  igenes_lcpm_df <- lcpm_df[match(igenes[[i]], lcpm_df$Gene), ]
  row.names(igenes_lcpm_df) <- NULL
  
  igenes_lcpm_mat <- igenes_lcpm_df %>% 
    column_to_rownames('Gene') %>% 
    as.matrix()
  
  ## heatmap
  filenm <- paste0(heatmap_dir, "top_genes_heatmap_", comparisons[i], ".pdf")
  pdf(file=filenm, height = 7,width = 10)
  heatmap.2(x=igenes_lcpm_mat,
            Rowv = F,
            Colv =F, 
            dendrogram = "none",
            scale="row",
            col="bluered",
            trace="none",
            # the side bar color need to customized each time
            ColSideColors=c(rep("green", 5), rep("orange", 4), rep("purple", 3), 
                            rep("grey", 7)),
            key = TRUE,
            keysize = 1,
            main="Heatmap of the most significantly changed genes")
  
  dev.off()
  
}
