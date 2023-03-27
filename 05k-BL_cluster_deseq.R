################################################################################
################################################################################
################################################################################
################################################################################
###################### BURKITT LYMPHOMA HERV CLUSTER DESEQ #####################

## Plan:
## 1.

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
library(pheatmap)

################################### LOAD DATA ##################################

load("r_outputs/02-BL_filt_counts.Rdata")
load("r_outputs/05i-BL_pca_ccp_clusters_metadata.Rdata")
load("r_outputs/01-refs.Rdata")

################################# METADATA SETUP ###############################

BL_metadata$clust.retro.k3 <- clust.df$clust.retro.k3

#################################### DESEQ2 ####################################

lfc.cutoff <- 1.5
pval=0.001 # p value threshold

stopifnot(all(colnames(BL.filt.comb) == rownames(BL_metadata)))

BL.hc.dds <- DESeq2::DESeqDataSetFromMatrix(BL.filt.comb, 
                                            BL_metadata, 
                                            ~ clust.retro.k3 + 0)


BL.hc.dds <- DESeq2::DESeq(BL.hc.dds, parallel=T)
BL.hc.tform <- DESeq2::varianceStabilizingTransformation(BL.hc.dds, blind=FALSE)

############################## TOP GENES & HERVS ###############################

## > resultsNames(BL.hc.dds)
## [1] "clust.retro.k3C1" "clust.retro.k3C2" "clust.retro.k3C3"

BL_res <- list(
  "C1" = DESeq2::results(BL.hc.dds, contrast=c(+1, -1/2, -1/2), alpha=pval),
  "C2" = DESeq2::results(BL.hc.dds, contrast=c(-1/2, +1, -1/2), alpha=pval),
  "C3" = DESeq2::results(BL.hc.dds, contrast=c(-1/2, -1/2, +1), alpha=pval),
  "C1vC2" = DESeq2::results(BL.hc.dds, contrast=c("clust.retro.k3", "C1", "C2"), alpha=pval),
  "C1vC3" = DESeq2::results(BL.hc.dds, contrast=c("clust.retro.k3", "C1", "C2"), alpha=pval),
  "C2vC3" = DESeq2::results(BL.hc.dds, contrast=c("clust.retro.k3", "C2", "C3"), alpha=pval)
  )

BL_res <- lapply(BL_res, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig <- lapply(BL_res, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig[[n]])
}

############################ UPSET GENES & HERVS ###############################

upvars_BL <- lapply(sig[1:3], function(r) rownames(subset(r, log2FoldChange>0)))
downvars_BL <- lapply(sig[1:3], function(r) rownames(subset(r, log2FoldChange<0)))

upset(fromList(upvars_BL), sets=c("C1", "C2", "C3"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))

upset(fromList(downvars_BL), sets=c("C1", "C2", "C3"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))

up.binmat.BL <- fromList(upvars_BL)
rn <- do.call(c, upvars_BL)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.BL) <- rn
rm(rn)

dn.binmat.BL <- fromList(downvars_BL)
rn <- do.call(c, downvars_BL)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.BL) <- rn
rm(rn)

################################# HEATMAPS #####################################

################################### SETUP ######################################

# Colors for cells
cols <- rgb_gsea(palette = c("default"), n = 14, alpha = 0.7, reverse = FALSE)

# Annotation column
df <- as.data.frame(BL_metadata[c("clust.retro.k3", "clinical_variant",
                                  "ebv_status", "subgroup", "gender")])

# Create colors for each group
annoCol <-  c(wes_palette("Chevalier1")[1], 
              wes_palette("Chevalier1")[2], 
              wes_palette("Chevalier1")[3])
names(annoCol) <- c("C1", "C2", "C3")
annoCol <- list(clust.retro.k3 = annoCol)

######################### UPREGULATED IN ALL GROUPS ############################

top.genes.hervs <- rownames(up.binmat.BL)

pdf("plots/05k-BL_top_hervs_genes_upregulated_all.pdf", height=15, width=15)
pheatmap(assay(BL.hc.tform)[top.genes.hervs,], 
         main="Upregulated Genes and HERvs, all clusters",
         cluster_rows=TRUE,
         show_rownames=FALSE,
         show_colnames = FALSE,
         color = cols,
         scale="row",
         breaks=seq(-3,3,length.out=14),
         cluster_cols=TRUE, 
         treeheight_row=0,
         annotation_col=df,
         annotation_colors = annoCol)
dev.off()

########################### DE GENES PER CLUSTER ###############################

makeheatmap <- function(topgenes, ...) {
  args <- list(...)
  mat <- assay(BL.hc.tform)[unique(topgenes), ]
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
           labels_row = gene_table[topgenes,]$display,
           fontsize = 8, fontsize_row = 6, border_color = NA,
           legend=T,
           treeheight_col=10,
           treeheight_row=10, 
           cutree_cols=4,
           ...
  )
}


for(clust in c("C1", "C2", "C3")) {
  tg <- rownames(sig[[clust]][1:75,])
  p <- makeheatmap(tg, main=paste0('DE in cluster ', clust))
  pdf(paste0("plots/05k-BL_top_de_hervs_", clust, ".pdf"), height=7, width=7)
  print(p)
  dev.off()
}

############################# VOLCANO DZ VS LZ #################################

pdf("plots/05k-BL_volcano_C1_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig$C1,
                lab = gene_table[rownames(sig$C1),]$display,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs C1')
dev.off()

pdf("plots/05k-BL_volcano_C2_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig$C2,
                lab = gene_table[rownames(sig$C2),]$display,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs C2')
dev.off()

pdf("plots/05k-BL_volcano_C3_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig$C3,
                lab = gene_table[rownames(sig$C3),]$display,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs C3')
dev.off()

pdf("plots/05k-BL_volcano_C1.pdf", height=8, width=8)
EnhancedVolcano(sig$C3,
                lab = rownames(sig$C3),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs C3')
dev.off()

########################## PLOT INDIVIDUAL HERVs ###############################

plot.counts <- function(df, gene) {
  
  title <- gene_table[gene,]$display
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="clust.retro.k3", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=clust.retro.k3, y=count, fill=clust.retro.k3))  +
    geom_boxplot() +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("Cluster") +
    ylab("Counts") +
    scale_fill_manual(values = c("C1" = wes_palette("Chevalier1")[1], 
                                 "C2" = wes_palette("Chevalier1")[2],
                                 "C3" = wes_palette("Chevalier1")[3])) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)
}

p1 <- plot.counts(BL.hc.dds, "ENSG00000104972.16")
p2 <- plot.counts(BL.hc.dds, "ENSG00000156508.19")
p3 <-  plot.counts(BL.hc.dds, "ENSG00000214182.5")
p4 <- plot.counts(BL.hc.dds, "ENSG00000103512.15")
p5 <- plot.counts(BL.hc.dds, "ENSG00000102898.12")
p6 <- plot.counts(BL.hc.dds, "ENSG00000196418.13")
p7 <- plot.counts(BL.hc.dds, "ENSG00000226803.9")
p8 <- plot.counts(BL.hc.dds, "ENSG00000213742.7")
p9 <- plot.counts(BL.hc.dds, "ENSG00000186063.13")

pdf("plots/05k-BL_clusters_individual_hervs.pdf", height=12, width=12)
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9,
          nrow = 3, 
          ncol = 3,
          labels = "AUTO")
dev.off()

pdf("plots/05k-BL_clusters_SNORD3A.pdf", height=3, width=3)
plot.counts(BL.hc.dds, "ENSG00000263934.5")
dev.off()

pdf("plots/05k-BL_clusters_YWHAG.pdf", height=3, width=3)
plot.counts(BL.hc.dds, "ENSG00000170027.7")
dev.off()

plot_grid(plot.counts(BL.hc.dds, "ENSG00000171401.15"), 
          plot.counts(BL.hc.dds, "ENSG00000185479.6"), 
          plot.counts(BL.hc.dds, "ENSG00000170465.10"), 
          plot.counts(BL.hc.dds, "ENSG00000205420.11"),
          nrow = 2, 
          ncol = 2,
          labels = "AUTO")


pdf("plots/05k-BL_clusters_HARLEQUIN_1q32.1.pdf", height=3, width=3)
plot.counts(BL.hc.dds, "HARLEQUIN_1q32.1")
dev.off()

pdf("plots/05k-BL_clusters_MER101_16p12.2a.pdf", height=3, width=3)
plot.counts(BL.hc.dds, "MER101_16p12.2a")
dev.off()

pdf("plots/05k-BL_clusters_MER101_16p12.2a.pdf", height=3, width=3)
plot.counts(BL.hc.dds, "MER101_16p12.2a")
dev.off()