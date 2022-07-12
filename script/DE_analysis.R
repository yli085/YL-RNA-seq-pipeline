## script to perform DE analysis, as well as selecting DEGs
# usage: 
# input are count matrix and meta data
# output include annotating DEGs, volcano plots, logCPM tables.
# general steps: first run the first several sections to get the overal DE calculate results and logCPM,
# and then manually adjust DEG cutoffs to extract DEGs.

# data folder structure
# report: DE_analysis + gene_expression
# within DE_analysis: DE_sum + DEG_lfcxx_padjxx
# within DEG_lfcxx_padjxx: DEG_output + volcano_plots (+ annotated DEG files)

# Setup -------------------------------------------------------------------
setwd("/Users/YL/Documents/Nutstore/Nutstore/TB_NGS/Projects/Neuro-007-All/Neuro_007_brain_75PE/scripts")

library(tidyverse)
library(readxl)
library(writexl)
library(limma)
library(edgeR)
library(ggrepel) 


# data input --------------------------------------------------------------
#@ data input: data, meta with the group info
data0 <- read_delim("../data/counts_Rmatrix.txt", delim="\t", col_types=list(Geneid="c", .default="i"))
meta0 <- read_excel("../data/meta/meta_Neuro_007_Brain.xlsx") # sheet 1

# data: make the Geneid into rowname of the count dataframe
data <- data0 %>% column_to_rownames("Geneid")
# select data
data <- data %>% select(1:19) 

# meta: rename the samples in the meta table
meta <- meta0 %>% column_to_rownames("Sample_rename")
# select data
meta <- meta[1:19,]

# check data and meta's sample names
all(names(data) == row.names(meta))
# rename data's column names if FALSE
# names(data) <- row.names(meta)

# meta: factorize groups etc. The levels will be shown at the legend of PCA plot
meta$sampletype <- factor(meta$Group2, levels=c("A", "B", "C", "WT_AG"))
meta$Gender <- factor(meta$Gender, levels=c("Male", "Female"))
str(meta)

## set the DE comparisons and output directory
# comparisons, use a dash to separate groups.
comparisons <- c("A-C", "B-C", "C-WT_AG")

# output dir for save results
out_dir <- "../report/"
# create sub folders
de_dir <- paste0(out_dir, "DE_analysis/")
expr_dir <- paste0(out_dir, "gene_expression/")
dir.create(de_dir)
dir.create(expr_dir)


# DE analysis: calculate statistics -------------------------------------------------------------------
## read data into a DGElist object
group <- as_factor(meta$sampletype)

# create a DGElist object
x <- DGEList(counts=data, group=group)

## remove lowly expressed genes
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
# dim(x)

## data normalization
x <- calcNormFactors(x, method = "TMM")
# x$samples$norm.factors

## create a design matrix and contrasts for limma
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))

contr.matrix <- makeContrasts(
  contrasts= comparisons,
  levels = levels(meta$sampletype))

## run limma calculation
par(mfrow=c(1,2))
v <- voom(x, design, plot=T)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

## calculate statistics
# create an empty list
gene_ls <- list()

for (i in seq_along(comparisons)) {
  gene_ls[[i]] <- topTable(efit, sort.by = "P", coef=i, number=Inf)
  gene_ls[[i]] <- gene_ls[[i]] %>% rownames_to_column("gene")
}

names(gene_ls) <- comparisons

# create the folder and save files
de_sum_dir <- paste0(de_dir, "DE_sum/")
dir.create(de_sum_dir)
saveRDS(gene_ls, file = paste0(de_sum_dir, "gene_ls.RDS"))


# All gene results table  ------------------------------------------------
# gene id annotation file
gene_id_df <- readRDS("/Users/YL/Documents/Nutstore/Nutstore/TB_NGS/RNAseq_method/genome_annotations/ensembl_ids_conversion.RDS")

# create an empty list
gene_ls2 <- list()

for (i in seq_along(comparisons)) {
  # select useful columns
  gene_ls2[[i]] <- gene_ls[[i]] %>%
    select(gene, AveExpr, logFC, P.Value, adj.P.Val)
  # add comparison names to the column
  names(gene_ls2[[i]])[3:5] <- paste0(names(gene_ls2[[i]])[3:5], " (", comparisons[i], ")")
}

# combine all tables together
gene_df <- gene_ls2 %>% purrr::reduce(left_join, by=c("gene", "AveExpr"))

# add gene annotations
gene_df <- gene_df %>% left_join(gene_id_df, by=c("gene"="gene_id")) %>%
  select(-c("entrez_id", "ensembl_id", "mgi_id")) %>%
  relocate(any_of(c("gene_name", "entrez_description")),  .after = gene)

# rename the first few columns
names(gene_df)[1:3] <- c("Gene_ID", "Gene_symbol", "Gene_Description")

# save files
write_xlsx(gene_df, path=paste0(de_sum_dir, "All_genes_DE_results.xlsx"))
# saveRDS(gene_df, file=paste0(de_sum_dir, "All_genes_DE_results.RDS"))


# log CPM data-------------------------------------------------------------
## keep ALL genes and TMM normalization on the counts
group <- as_factor(meta$sampletype)
x <- DGEList(counts=data, group=group)

# calculate and save the normalization factor in x
x <- calcNormFactors(x, method = "TMM")
# calculate normalized cpm
norm_cpm <- cpm(x, normalized.lib.sizes=T, log = F)
# log
lcpm <- log2(norm_cpm+0.1)
lcpm_df <- as.data.frame(lcpm) %>% rownames_to_column("gene")

# save files into "gene_expression" folder
filenm <- paste0(expr_dir, "logCPM.txt")
write.table(lcpm_df, filenm, sep="\t", row.names = F, quote = F, na="")


# DEG selection cutoffs and the summary table-------------------------------------------------------------
padj.cutoff=0.05
lfc.cutoff=0

# save the following results in this folder
deg_dir <- paste0(de_dir, "DEG_lfc", lfc.cutoff, "_padj", padj.cutoff, "/")
# create this folder
dir.create(deg_dir)

## get the DEG number sum table using the limma data
sumTable<- summary(decideTests(efit, lfc=lfc.cutoff,p.value=padj.cutoff))
sumDf <- as.data.frame(unclass(t(sumTable)))
sumDf <- sumDf %>% 
  rownames_to_column("Comparison") %>% 
  mutate("All DEG" = Down + Up) %>% 
  relocate(any_of(c("All DEG", "Up", "Down")), .after=Comparison)
sumDf

# save files
filenm = paste0(deg_dir, "DEG_summary.xlsx")
write_xlsx(sumDf, path=filenm)


# extract DEGs ------------------------------------------------------------
## extract DEGs based on cutoffs
# gene_ls <- readRDS(paste0(out_dir, "DE_analysis/DE_sum/gene_ls.RDS"))

# create the empty list
deg_ls <- list()

# create the output folder
deg_raw_dir <- paste0(deg_dir, "DEG_output/")
dir.create(deg_raw_dir)

for (i in seq_along(comparisons)) {
  deg_ls[[i]] <- gene_ls[[i]] %>% dplyr::filter(abs(logFC) > lfc.cutoff, adj.P.Val < padj.cutoff)
  # save files
  filenm <- paste0(deg_raw_dir, comparisons[i], "_DEG.txt")
  write.table(deg_ls[[i]], filenm, sep="\t", row.names = F, quote = F, na="")
}


# annotate DEGs --------------------------------------------------------
gene_id_df <- readRDS("/Users/YL/Documents/Nutstore/Nutstore/TB_NGS/RNAseq_method/genome_annotations/ensembl_ids_conversion.RDS")

for (i in seq_along(comparisons)) {
  
  degAnno <- deg_ls[[i]] %>% 
    select(gene, AveExpr, logFC, P.Value, adj.P.Val)
  
  # add gene annotations
  degAnno <- degAnno %>% left_join(gene_id_df, by=c("gene"="gene_id")) %>% 
    select(-c("entrez_id", "ensembl_id", "mgi_id")) 
  
  # add Regulation direction column and reorder columns
  degAnno <- degAnno %>% 
    mutate(Regulation = ifelse(degAnno$logFC > 0, "Up", "Down")) %>% 
    relocate(any_of(c("gene_name", "entrez_description", "Regulation")),  .after = gene)
  
  # rename the first few columns
  names(degAnno)[1:3] <- c("Gene_ID", "Gene_symbol", "Gene_Description")
  
  # save files
  filenm <- paste0(deg_dir, comparisons[i], "_DEG.xlsx")
  write_xlsx(degAnno, path=filenm) 
}

# volcano plot ------------------------------------------------------------

# create output directory
deg_vol_dir <- paste0(deg_dir, "volcano_plots/")
dir.create(deg_vol_dir)


for (i in seq_along(comparisons)){
  
  df <- gene_ls[[i]]
  df <- df %>% left_join(gene_id_df, by=c("gene"="gene_id"))
  df <- df %>% dplyr::select(ensembl_id, logFC, adj.P.Val, gene_name)
  
  df <- df %>% mutate(threshold = adj.P.Val < padj.cutoff & abs(logFC) > lfc.cutoff)
  df <- df %>% arrange(-threshold, adj.P.Val) %>% mutate(genelabels = "")
  df$genelabels[1:10] <- df$gene_name[1:10]
  
  df <- df %>% mutate(dir = 
                        case_when(df$logFC > lfc.cutoff & df$adj.P.Val < padj.cutoff ~ "Up-regulated",
                                  df$logFC < -lfc.cutoff & df$adj.P.Val < padj.cutoff ~ "Down-regulated",
                                  TRUE ~ "Nonsignificant"))
  
  # volcano plot
  # plot tile
  vs_groups <- unlist(strsplit(comparisons[i], split="-"))
  
  # save files
  filenm <- paste0(deg_vol_dir, comparisons[i], "_volcano_plot.pdf")
  
  pdf(file=filenm, width=7, height=5)
  p1 <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = dir), size=1, alpha = 0.8) +
    scale_color_manual(name = "Directionan",
                       values = c("Up-regulated" = "#F8766D", "Down-regulated" = "#7CAE00", "Nonsignficant" = "grey70"))+
    geom_text_repel(aes(label = genelabels),
                    size= rel(3),
                    box.padding = 0.2,
                    show.legend=FALSE,
                    # max.overlaps= Inf
    ) +
    ggtitle(label = paste0(vs_groups[1], " VS ", vs_groups[2]))+
    geom_hline(yintercept = -log10(padj.cutoff), linetype='dashed', size=0.5, color="grey40") +
    geom_vline(xintercept= c(-lfc.cutoff, lfc.cutoff), linetype='dashed', size=0.5, color="grey40")+
    labs(x = "log2 fold change", 
         y = "-log10 adjusted p-value") + 
    theme_bw() +
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5, face = 'bold'),
          axis.title = element_text(size = rel(1.25), face = 'bold'))+
    theme(legend.text = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'),
          legend.position = 'right')
  print(p1)
  dev.off()
  
}
