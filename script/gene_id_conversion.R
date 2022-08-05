# read gtf and extract genes----------

library(rtracklayer)

# gtf is gff2, version = 2
gtf0 <- readGFF('/Users/YL/Documents/Nutstore/Nutstore/TB_NGS/RNAseq_method/genome_annotations/mouse.gencode.vM27.primary_assembly.annotation.protein_coding_only.gtf.gz', version=2L)
# head(gtf0)
# dim(gtf0)

gtf <- gtf0[with(gtf0, type=='gene'),]
# head(gtf)
# dim(gtf)

gtf2 <- gtf %>% separate(gene_id, into = "ensembl_id", sep="\\.", remove=FALSE, extra="drop")

# check duplicates
which(duplicated(gtf2$ensembl_id)) # None
which(duplicated(gtf2$gene_name)) # 25

# filter <- which(duplicated(genes)) # empty 
# gtf <- gtf[-filter,]
# rownames(gtf) <- genes[-filter]

saveRDS(gtf2, file='genome_annotations/GRCm39.vM27.gencode.protein.coding.genes.gtf.rds')


# make gene_id conversion table --------------------------------

gtf <- readRDS("/Users/YL/Documents/Nutstore/Nutstore/TB_NGS/RNAseq_method/genome_annotations/GRCm39.vM27.gencode.protein.coding.genes.gtf.rds")
id_df <- gtf[c("gene_id", "ensembl_id", "gene_name", "mgi_id")]


## biomart: ensembl to entrez table
df <- read.delim("/Users/YL/Documents/Nutstore/Nutstore/TB_NGS/RNAseq_method/genome_annotations/ensembl_to_entrez.txt", stringsAsFactors = F, header=T, col.names= c("ensembl_id", "ensembl_id_v", "entrez_id", "entrez_description"))
sum(!is.na(df$entrez_id)) # 27749

# compare the gencode gtf and biomart table, remove duplicated match
id_df <- id_df %>% left_join(df, by="ensembl_id")
id_df <- id_df %>% filter(!duplicated(id_df$ensembl_id))

# remove emsembl version column which are consistent between two sources
# all(id_df$gene_id == id_df$ensembl_id_v) # TRUE
id_df <- id_df %>% dplyr::select(-ensembl_id_v) 

dim(id_df)[1] # 21885
sum(!is.na(id_df$gene_name)) # 21885 All with gene name
sum(!is.na(id_df$entrez_id)) # 21038

saveRDS(id_df, file = "/Users/YL/Documents/Nutstore/Nutstore/TB_NGS/RNAseq_method/genome_annotations/ensembl_ids_conversion.RDS")
