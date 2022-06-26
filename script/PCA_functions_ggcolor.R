# script including three PCA function: 
# PCAdotlab: use both dot to represent sample with label annotation.
# PCAdot: only use dot to represent sample.
# PCAlab: use label to represent sample.
# input
# object: data 
# PC_1 and PC_2: the numbers  of two PCs
# intgroup="condition"
# ntop: use n number of top genes, default=500.
# PCA_label: the sample labels used in the PCA plot, i.e. meta$PCA_label. Not required in PCAdot.
# Note: color is the default ggcolor

# library
library(ggrepel) 
library(RColorBrewer)


# PCAdotlab --------------------


# each sample has a text label and a dot
PCAdotlab <- function(object=vsd, intgroup="sampletype", PC_1=1, PC_2=2, ntop=500, PCA_label){
  
  ## get the PC data
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,PC_1], PC2=pca$x[,PC_2], group=group, intgroup.df, name=colnames(object))
  
  ## generate a plot
  # data point labels are the sample ID
  d$label = PCA_label
  
  # ggplot: text as labels
  ggplot(data=d, aes_string(x="PC1", y="PC2", color = "group")) + 
    coord_fixed() + 
    theme_bw() +
    geom_point(aes(color=group),alpha =0.8, size=3) +
    geom_text_repel(data=d,
                    aes(label=label),
                    box.padding = 0.25,
                    show.legend=FALSE
    )+
    labs(title = "PCA plot") +
    xlab(paste0("PC",PC_1,": ",round(percentVar[PC_1] * 100),"% variance")) +
    ylab(paste0("PC",PC_2,": ",round(percentVar[PC_2] * 100),"% variance")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(color = "black", size=10, face="bold"),
          axis.text.y = element_text(color = "black", size=10, face="bold"),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.border = element_rect(color = "black", fill=NA, size=1))
  
}


# PCAdot ---------------------------------------------


PCAdot <- function(object=vsd, intgroup="sampletype", PC_1=1, PC_2=2, ntop=500){
  
  ## get the PC data
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,PC_1], PC2=pca$x[,PC_2], group=group, intgroup.df, name=colnames(object))
  
  ## generate a plot
  # ggplot: text as labels
  ggplot(data=d, aes_string(x="PC1", y="PC2", color = "group")) + 
    coord_fixed() + 
    theme_bw() +
    geom_point(aes(color=group),alpha =0.8, size=3) +
    labs(title = "PCA plot") +
    xlab(paste0("PC",PC_1,": ",round(percentVar[PC_1] * 100),"% variance")) +
    ylab(paste0("PC",PC_2,": ",round(percentVar[PC_2] * 100),"% variance")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(color = "black", size=10, face="bold"),
          axis.text.y = element_text(color = "black", size=10, face="bold"),
          # legend.position = "none",
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.border = element_rect(color = "black", fill=NA, size=1))
  
}


# PCAlab ---------------------------------------------

PCAlab <- function(object=vsd, intgroup="sampletype", PC_1=1, PC_2=2, ntop=500, PCA_label) {
  
  ## get the PC data
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,PC_1], PC2=pca$x[,PC_2], group=group, intgroup.df, name=colnames(object))
  
  ## generate a plot
  # data point labels are the sample ID
  d$label = PCA_label
  
  # ggplot: text as labels
  ggplot(data=d, aes_string(x="PC1", y="PC2", color = "group")) + 
    coord_fixed() + 
    theme_bw() +
    # geom_point(aes(color=group), size=3) +
    geom_text(data=d,
              aes(label=label),
              # box.padding = 0.5,
              show.legend=T, size = 4
    )+
    labs(title = "PCA plot") +
    xlab(paste0("PC",PC_1,": ",round(percentVar[PC_1] * 100),"% variance")) +
    ylab(paste0("PC",PC_2,": ",round(percentVar[PC_2] * 100),"% variance")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(color = "black", size=10, face="bold"),
          axis.text.y = element_text(color = "black", size=10, face="bold"),
          # legend.position = "none",
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.border = element_rect(color = "black", fill=NA, size=1))
  
  
}
