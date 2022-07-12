## script to get the numbers of DEGs or top genes using a set of cutoffs
# input is the limma output "gene_ls"
# output is the excel table of gene numbers, one for 4 padj&lfc cutoffs, the other is for 2 raw p-value cutoffs
# add a plot later...

library(tidyverse)
library(writexl)


# DEG numbers-------------------------------------------------------------
# input
gene_ls <- readRDS('../report/report_wo_aging/DE_analysis/DE_sum/gene_ls.RDS')
# comparison groups
comparisons <- names(gene_ls)

# condition groups
padj_vec <- c(0.05, 0.05, 0.01, 0.01)
lfc_vec <- c(0, 0.58, 0.58, 1)
conditions <- c()
for (i in seq_along(padj_vec)) {
  conditions[i] <-  paste0("padj", padj_vec[i], "&lfc", lfc_vec[i])
}

# the number matrix
tgene_num <- matrix(NA, nrow=length(comparisons), ncol=length(conditions))
rownames(tgene_num) <- comparisons
colnames(tgene_num) <- conditions


# The number of genes passing filters from all the conditions
for (i in seq_along(comparisons)) {
  for (j in seq_along(conditions)) {
    df <- gene_ls[[i]] %>% dplyr::filter(abs(logFC) > lfc_vec[j], adj.P.Val < padj_vec[j])
    tgene_num[i, j] <- dim(df)[1]
  }
}

tgene_num_df <- as.data.frame(tgene_num) %>% rownames_to_column("Comparisons")

# save files
write_xlsx(tgene_num_df, "../report/report_wo_aging/DE_analysis/DE_sum/DEG_numbers.xlsx")


# Top genes using raw P-value ---------------------------------------------------------
# condition groups
rawP_vec <- c(0.01, 0.05)
conditions <- paste0("rawP", rawP_vec)

# the number matrix
tgene_num <- matrix(NA, nrow=length(comparisons), ncol=length(conditions))
rownames(tgene_num) <- comparisons
colnames(tgene_num) <- conditions


# The number of genes passing filters from all the conditions
for (i in seq_along(comparisons)) {
  for (j in seq_along(conditions)) {
    df <- gene_ls[[i]] %>% dplyr::filter(P.Value < rawP_vec[j])
    tgene_num[i, j] <- dim(df)[1]
  }
}

tgene_num_df <- as.data.frame(tgene_num) %>% rownames_to_column("Comparisons")

# save files
write_xlsx(tgene_num_df, "../report/report_wo_aging/DE_analysis/DE_sum/top_gene_numbers.xlsx")

