################################################################################
################################################################################
################################################################################
################################################################################
############################ BULK GCB PCA & DESEQ2 #############################

## Plan:
## 1. PCA of bulk GCB data 
## 2. Identift TEs upregulated in all portions of the GCB
## 3. Identify TEs specific to LZ, DZ, PB, MB

#################################### SETUP #####################################

library(knitr)
library(tidyverse)
library(matrixStats)
library(data.table)
library(PCAtools)
library(DESeq2)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(edgeR)
library(ashr)
library(cowplot)
library(wesanderson)
library(UpSetR)
library(EnhancedVolcano)

################################### LOAD DATA ##################################

load("r_outputs/02-GCB_Bulk_filt_counts.Rdata")
load("r_outputs/01-metadata.Rdata")
load("r_outputs/01-refs.Rdata")

remove(BL_metadata, DLBCL_metadata, FL_metadata, all_metadata)

############################# FUNCTION SCREE PLOT ##############################

pca_standard <- function(tform, metadata, var) {
  
  removeVar <- var
  pca.obj <- PCAtools::pca(assay(tform), 
                           metadata=metadata, 
                           removeVar=removeVar)
  
  cat(sprintf('Removed %d pct low variance variables, %d retained\n', 
              removeVar*100, length(pca.obj$xvars)))
  
  varline <- 50
  varline.x <- min(which(cumsum(pca.obj$variance) >= varline))
  
  horn <- PCAtools::parallelPCA(assay(tform), removeVar = removeVar)
  elbow <- PCAtools::findElbowPoint(pca.obj$variance)
  
  screeplot <-PCAtools::screeplot(pca.obj,
                                  axisLabSize = 6,
                                  components = getComponents(pca.obj, 1:30),
                                  title=paste("Retrotranscriptome SCREE",
                                              metadata$cancer_type[1],
                                              sep=" "),
                                  hline=varline, vline=c(varline.x, horn$n, elbow)
  ) +
    geom_label(aes(x=varline.x+1, y=50, 
                   label = paste0(varline, '% var'), vjust = -1)) +
    geom_label(aes(x = horn$n + 1, y = 50,
                   label = 'Horn\'s', vjust = -1)) +
    geom_label(aes(x = elbow + 1, y = 50,
                   label = 'Elbow method', vjust = -1))
  
  
  cat(sprintf('%d PCs for Elbow method\n', elbow))
  cat(sprintf('%d PCs for Horn method\n', horn$n))
  cat(sprintf('%d PCs needed to explain %d percent of variation\n', 
              varline.x, varline))
  
  return(pca.obj)
}

################################ GCB BULK DESEQ ################################

### DESeq2 (HERVs Only)

bulk.countdat <- GCB_Bulk.filt.herv
cat(sprintf('%d variables\n', nrow(bulk.countdat)))

stopifnot(all(colnames(bulk.countdat) == rownames(bulk_metadata)))

GCB.dds <- DESeq2::DESeqDataSetFromMatrix(countData = bulk.countdat,
                                            colData = bulk_metadata,
                                            design = ~ Cell_type + 0)

GCB.dds <- DESeq2::DESeq(GCB.dds, parallel=T)
GCB.tform <- DESeq2::varianceStabilizingTransformation(GCB.dds, blind=FALSE)

## PCA
GCB.herv.pca.obj <-
  pca_standard(tform = GCB.tform, 
               metadata = bulk_metadata, 
               var = 0.1)

# 3 PCs for Elbow method
# 2 PCs for Horn method
# 1 PCs needed to explain 50 percent of variation

######################## GCB BULK DESEQ (HERVS + GENES) ########################

### DESeq2 (HERVs + Genes)

bulk.g.countdat <- GCB_Bulk.filt.comb
cat(sprintf('%d variables\n', nrow(bulk.g.countdat)))

stopifnot(all(colnames(bulk.g.countdat) == rownames(bulk_metadata)))

GCB.g.dds <- DESeq2::DESeqDataSetFromMatrix(countData = bulk.g.countdat,
                                            colData = bulk_metadata,
                                            design = ~ Cell_type + 0)

GCB.g.dds <- DESeq2::DESeq(GCB.g.dds, parallel=T)
GCB.g.tform <- DESeq2::varianceStabilizingTransformation(GCB.g.dds, blind=FALSE)

## PCA
GCB.g.pca.obj <-
  pca_standard(tform = GCB.g.tform, 
               metadata = bulk_metadata, 
               var = 0.1)

# 3 PCs for Elbow method
# 4 PCs for Horn method
# 1 PCs needed to explain 50 percent of variation

###################### BULK GCB BIPLOTS HERVs & GENES ##########################

## Biplot with projects (HERVs and genes)

biplot(GCB.g.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "Cell_type",
       colkey = c("Naive B" = pal_jco("default", alpha = 0.7)(5)[1],
                  "Germinal Center B" = pal_jco("default", alpha = 0.7)(5)[2],
                  "Memory B" = pal_jco("default", alpha = 0.7)(5)[3],
                  "Dark Zone Germinal Center B" = pal_jco("default", alpha = 0.7)(5)[4],
                  "Light Zone Germinal Center B" = pal_jco("default", alpha = 0.7)(5)[5]),
       legendPosition = "right")  +
  theme_cowplot()

## Biplot with projects (HERVs only)
biplot(GCB.herv.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "Cell_type",
       colkey = c("Naive B" = pal_jco("default", alpha = 0.7)(5)[1],
                  "Germinal Center B" = pal_jco("default", alpha = 0.7)(5)[2],
                  "Memory B" = pal_jco("default", alpha = 0.7)(5)[3],
                  "Dark Zone Germinal Center B" = pal_jco("default", alpha = 0.7)(5)[4],
                  "Light Zone Germinal Center B" = pal_jco("default", alpha = 0.7)(5)[5]),
       legendPosition = "right")  +
  theme_cowplot()

################################ SET THRESHOLDS ################################

fc=2 # fold change
l2fc=log2(fc) # log 2 fold change
lfc.cutoff <- 1.5
pval=0.001 # p value threshold

############################## TOP GENES & HERVS ###############################

gcb_res <- list(
  "DZ" = DESeq2::results(GCB.g.dds, contrast=c(+1, -1/4, -1/4, -1/4, -1/4), alpha=pval),
  "GCB" = DESeq2::results(GCB.g.dds, contrast=c(-1/4, +1, -1/4, -1/4, -1/4), alpha=pval),
  "LZ" = DESeq2::results(GCB.g.dds, contrast=c(-1/4, -1/4, +1, -1/4, -1/4), alpha=pval),
  "MB" = DESeq2::results(GCB.g.dds, contrast=c(-1/4, -1/4, -1/4, +1, -1/4), alpha=pval),
  "NB" = DESeq2::results(GCB.g.dds, contrast=c(-1/4, -1/4, -1/4, -1/4, 1), alpha=pval),
  "DZvLZ" = DESeq2::results(GCB.g.dds, contrast=c(+1, 0, +1, 0, 0), alpha=pval),
  "MBvsNB" = DESeq2::results(GCB.g.dds, contrast=c(0, 0, 0, +1, +1), alpha=pval)
)

gcb_res <- lapply(gcb_res, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig <- lapply(gcb_res, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig[[n]])
}

############################### TOP HERVs ONLY #################################

gcb_res_herv <- list(
  "DZ" = DESeq2::results(GCB.dds, contrast=c(+1, -1/4, -1/4, -1/4, -1/4), alpha=pval),
  "GCB" = DESeq2::results(GCB.dds, contrast=c(-1/4, +1, -1/4, -1/4, -1/4), alpha=pval),
  "LZ" = DESeq2::results(GCB.dds, contrast=c(-1/4, -1/4, +1, -1/4, -1/4), alpha=pval),
  "MB" = DESeq2::results(GCB.dds, contrast=c(-1/4, -1/4, -1/4, +1, -1/4), alpha=pval),
  "NB" = DESeq2::results(GCB.dds, contrast=c(-1/4, -1/4, -1/4, -1/4, 1), alpha=pval),
  "DZvLZ" = DESeq2::results(GCB.dds, contrast=c(+1, 0, +1, 0, 0), alpha=pval),
  "MBvsNB" = DESeq2::results(GCB.dds, contrast=c(0, 0, 0, +1, +1), alpha=pval)
)

gcb_res_herv <- lapply(gcb_res_herv, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})


sig_herv <- lapply(gcb_res_herv, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})


for (n in names(sig_herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig_herv[[n]])
}

############################ UPSET GENES & HERVS ###############################

upvars <- lapply(sig[1:5], function(r) rownames(subset(r, log2FoldChange>0)))
downvars <- lapply(sig[1:5], function(r) rownames(subset(r, log2FoldChange<0)))

upset(fromList(upvars), sets=c("DZ", "GCB", "LZ", "MB", "NB"),  
      keep.order = T, order.by='degree', decreasing=F)

upset(fromList(downvars), sets=c("DZ", "GCB", "LZ", "MB", "NB"),  
      keep.order = T, order.by='degree', decreasing=F)

up.binmat <- fromList(upvars)
rn <- do.call(c, upvars)
rn <- rn[!duplicated(rn)]
rownames(up.binmat) <- rn
rm(rn)

dn.binmat <- fromList(downvars)
rn <- do.call(c, downvars)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat) <- rn
rm(rn)

############################# UPSET HERVs ONLY #################################

upvars_hervs <- lapply(sig_herv[1:5], function(r) rownames(subset(r, log2FoldChange>0)))
downvars_hervs <- lapply(sig_herv[1:5], function(r) rownames(subset(r, log2FoldChange<0)))

upset(fromList(upvars_hervs), sets=c("DZ", "GCB", "LZ", "MB", "NB"),  
      keep.order = T, order.by='degree', decreasing=F)

upset(fromList(downvars_hervs), sets=c("DZ", "GCB", "LZ", "MB", "NB"),  
      keep.order = T, order.by='degree', decreasing=F)

up.binmat.hervs <- fromList(upvars_hervs)
rn <- do.call(c, upvars_hervs)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.hervs) <- rn
rm(rn)

dn.binmat.hervs <- fromList(downvars_hervs)
rn <- do.call(c, downvars_hervs)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.hervs) <- rn
rm(rn)

################################# HEATMAPS #####################################

################################### SETUP ######################################

# Colors for cells
cols <- rgb_gsea(palette = c("default"), n = 14, alpha = 0.7, reverse = FALSE)

# Annotation column
df <- as.data.frame(colData(GCB.dds)[,c("BioSample","Cell_type")])
df <- subset(df, select = -c(1))

# Create colors for each group
annoCol <-  pal_jco("default", alpha = 0.7)(5)
names(annoCol) <- unique(df$Cell_type)
annoCol <- list(Cell_type = annoCol)

######################### UPREGULATED IN ALL GROUPS ############################

top.genes.hervs <- rownames(up.binmat)
top.hervs <- rownames(up.binmat.hervs)

pheatmap(assay(GCB.tform)[top.hervs,], 
         main="Upregulated HERVs, all clusters",
         cluster_rows=TRUE,
         show_rownames=FALSE,
         show_colnames = FALSE,
         color = cols,
         scale="row",
         breaks=seq(-3,3,length.out=14),
         labels_row = gene_table[top.hervs,]$display,
         cluster_cols=TRUE, 
         treeheight_row=0,
         annotation_col=df,
         annotation_colors = annoCol)

pheatmap(assay(GCB.g.tform)[top.genes.hervs,],
         main="Upregulated Genes and HERVs, all clusters",
         cluster_rows=TRUE,
         show_rownames=FALSE,
         show_colnames = FALSE,
         color = cols,
         scale="row",
         breaks=seq(-3,3,length.out=14),
         labels_row = gene_table[top.genes.hervs,]$display,
         cluster_cols=TRUE, 
         treeheight_row=0,
         annotation_col=df,
         annotation_colors = annoCol)


########################### DE HERVs PER GC SITE ###############################


makeheatmap <- function(topgenes, ...) {
  args <- list(...)
  mat <- assay(GCB.tform)[unique(topgenes), ]
  rowdist <- as.dist((1 - cor(t(mat), method='spearman'))/2)
  
  annotation_col <- df
  cols <- rgb_gsea(palette = c("default"), n = 14, alpha = 0.7, reverse = FALSE)
  
  annotation_row <- data.frame(
    row.names = rownames(mat),
    Group=retro.annot[rownames(mat),]$family,
    Chrom=retro.annot[rownames(mat),]$chrom
  )
  
  pheatmap(mat,
           color=cols,
           scale="row", breaks=seq(-3,3,length.out=14),
           clustering_distance_rows = rowdist,
           clustering_method="average",
           annotation_col = df,
           annotation_colors = annoCol,
           show_colnames=F, 
           show_rownames = T,
           fontsize = 8, fontsize_row = 6, border_color = NA,
           legend=T,
           treeheight_col=10,
           treeheight_row=10, 
           cutree_cols=4,
           ...
  )
}


for(clust in c("DZ", "GCB", "LZ", "MB", "NB")) {
  tg <- rownames(sig_herv[[clust]][1:75,])
  p <- makeheatmap(tg, main=paste0('DE in cluster ', clust))
  print(p)
}

############################# VOLCANO DZ VS LZ #################################

EnhancedVolcano(sig_herv$DZ,
                lab = rownames(sig_herv$DZ),
                x = 'log2FoldChange',
                y = 'pvalue')

EnhancedVolcano(sig_herv$LZ,
                lab = rownames(sig_herv$LZ),
                x = 'log2FoldChange',
                y = 'pvalue')

EnhancedVolcano(sig_herv$NB,
                lab = rownames(sig_herv$NB),
                x = 'log2FoldChange',
                y = 'pvalue')

EnhancedVolcano(sig_herv$MB,
                lab = rownames(sig_herv$MB),
                x = 'log2FoldChange',
                y = 'pvalue')

EnhancedVolcano(sig_herv$DZvLZ,
                lab = rownames(sig_herv$DZvLZ),
                x = 'log2FoldChange',
                y = 'pvalue')


########################## PLOT INDIVIDUAL HERVs ###############################

plot.counts <- function(df, gene) {
  
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="Cell_type", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=Cell_type, y=count, fill=Cell_type))  +
    geom_boxplot() +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("Cell Type") +
    ylab("Counts") +
    scale_x_discrete(labels=c("Dark Zone Germinal Center B" = "DZ",
                              "Germinal Center B" = "GCB",
                              "Light Zone Germinal Center B" = "LZ",
                              "Memory B" = "MB",
                              "Naive B" = "NB")) +
    scale_fill_manual(values = c("Naive B" = pal_jco("default", alpha = 0.7)(5)[1],
                                 "Germinal Center B" = pal_jco("default", alpha = 0.7)(5)[2],
                                 "Memory B" = pal_jco("default", alpha = 0.7)(5)[3],
                                 "Dark Zone Germinal Center B" = pal_jco("default", alpha = 0.7)(5)[4],
                                 "Light Zone Germinal Center B" = pal_jco("default", alpha = 0.7)(5)[5])) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    ggtitle(gene) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)
}


plot.counts(GCB.g.dds, "HML2_3q12.3")
plot.counts(GCB.g.dds, "ERV316A3_6p21.33c")
plot.counts(GCB.g.dds, "MER4_17q21.2d")
plot.counts(GCB.g.dds, "HML3_5p15.33d")
plot.counts(GCB.g.dds, "HERVEA_5q22.2")
plot.counts(GCB.g.dds, "HERVH_7q11.21")
plot.counts(GCB.g.dds, "HERVH_1q32.1b")
plot.counts(GCB.g.dds, "HERVW_1q32.1")
plot.counts(GCB.g.dds, "HML5_1q22")



                       
                       