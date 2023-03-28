################################################################################
################################################################################
################################################################################
################################################################################
######################## ALL LYMPHOMA HERV CLUSTER DESEQ #######################

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

load("r_outputs/05b-all_lymphoma_pca_dds.Rdata")
load("r_outputs/01-metadata.Rdata")
load("r_outputs/01-refs.Rdata")

################################ SET THRESHOLDS ################################

fc=2 # fold change
l2fc=log2(fc) # log 2 fold change
lfc.cutoff <- 1.5
pval=0.001 # p value threshold

############################## TOP GENES & HERVS ###############################

resultsNames(all.g.dds)
# [1] "cancer_typeBL"    "cancer_typeDLBCL" "cancer_typeFL"                       

all_res <- list(
  "BL" = DESeq2::results(all.g.dds, contrast=c(+1, -1/2, -1/2), alpha=pval),
  "DLBCL" = DESeq2::results(all.g.dds, contrast=c(-1/2, +1, -1/2), alpha=pval),
  "FL" = DESeq2::results(all.g.dds, contrast=c(-1/2, -1/2, +1), alpha=pval),
  "BLvDLBCL" = DESeq2::results(all.g.dds, contrast=c("cancer_type", "BL", "DLBCL"), alpha=pval),
  "BLvFL" = DESeq2::results(all.g.dds, contrast=c("cancer_type", "BL", "FL"), alpha=pval),
  "DLBCLvFL" = DESeq2::results(all.g.dds, contrast=c("cancer_type", "DLBCL", "FL"), alpha=pval)
)

all_res <- lapply(all_res, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig <- lapply(all_res, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig[[n]])
}

############################### TOP HERVs ONLY #################################

resultsNames(all.dds)
# [1] "cancer_typeBL"    "cancer_typeDLBCL" "cancer_typeFL"                       

all_res_herv <- list(
  "BL" = DESeq2::results(all.dds, contrast=c(+1, -1/2, -1/2), alpha=pval),
  "DLBCL" = DESeq2::results(all.dds, contrast=c(-1/2, +1, -1/2), alpha=pval),
  "FL" = DESeq2::results(all.dds, contrast=c(-1/2, -1/2, +1), alpha=pval),
  "BLvDLBCL" = DESeq2::results(all.dds, contrast=c("cancer_type", "BL", "DLBCL"), alpha=pval),
  "BLvFL" = DESeq2::results(all.dds, contrast=c("cancer_type", "BL", "FL"), alpha=pval),
  "DLBCLvFL" = DESeq2::results(all.dds, contrast=c("cancer_type", "DLBCL", "FL"), alpha=pval)
)

all_res_herv <- lapply(all_res_herv, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})


sig_herv <- lapply(all_res_herv, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

sink(file = "r_outputs/05l-all_lymp_herv_deseq.txt")
for (n in names(sig_herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig_herv[[n]])
}
sink(file = NULL)

sink(file = "r_outputs/05l-all_lymp_hervv_deseq_top5.txt")
for (n in names(sig_herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  print(head(sig_herv[[n]]))
}
sink(file = NULL)

############################ UPSET GENES & HERVS ###############################

upvars <- lapply(sig[1:6], function(r) rownames(subset(r, log2FoldChange>0)))
downvars <- lapply(sig[1:6], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05l-all_lymph_upset_upvars.pdf", height=5, width=7)
upset(fromList(upvars), sets=c("DLBCL", "BL", "FL"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

pdf("plots/05l-all_lymph_upset_dnwars.pdf", height=5, width=7)
upset(fromList(downvars), sets=c("DLBCL", "BL", "FL"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

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

pdf("plots/05l-all_lymph_upset_upvars_hervs.pdf", height=5, width=7)
upset(fromList(upvars_hervs), sets=c("DLBCL", "BL", "FL"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 2, 1.5, 1.5, 1.5, 1.5))
dev.off()

pdf("plots/05l-all_lymph_upset_dnvars_hervs.pdf", height=5, width=7)
upset(fromList(downvars_hervs), sets=c("DLBCL", "BL", "FL"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 2, 1.5, 1.5, 1.5, 1.5))
dev.off()

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
df <- as.data.frame(colData(all.dds)[,c("cancer_type","subtype")])

# Create colors for each group
annoCol <-  pal_jco("default", alpha = 0.7)(3)
names(annoCol) <- unique(df$cancer_type)
annoCol2 <- c("#DD8D29", "#798E87", "#F4B5BD", "#9C964A", "#46ACC8", "#E2D200",
              "#CDC08C", "#FAD77B", "#E58601", "#B40F20", "#85D4E3")
names(annoCol2) <- unique(df$subtype)
annoCol <- list(cancer_type = annoCol, subtype = annoCol2)

######################### UPREGULATED IN ALL GROUPS ############################

top.genes.hervs <- rownames(up.binmat)
top.hervs <- rownames(up.binmat.hervs)

pdf("plots/05l-all_lymphoma_top_hervs_upregulated_all.pdf", height=10, width=10)
pheatmap(assay(all.dds)[top.hervs,], 
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
dev.off()


########################### DE HERVs PER GC SITE ###############################

makeheatmap <- function(topgenes, ...) {
  args <- list(...)
  mat <- assay(all.tform)[unique(topgenes), ]
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

for(clust in c("BL", "DLBCL", "FL")) {
  tg <- rownames(sig_herv[[clust]][1:75,])
  p <- makeheatmap(tg, main=paste0('DE in cluster ', clust))
  pdf(paste0("plots/05l-all_lymphoma_top_de_hervs_", clust, ".pdf"), height=7, width=7)
  print(p)
  dev.off()
}

######################### VOLCANO BETWEEN LYMPHOMAS ############################

pdf("plots/05l-all_lymphoma_volcano_DLBCL_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$DLBCL,
                lab = rownames(sig_herv$DLBCL),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs DLBCL')
dev.off()

pdf("plots/05l-all_lymphoma_BL_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$BL,
                lab = rownames(sig_herv$BL),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs BL')
dev.off()


pdf("plots/05l-all_lymphoma_volcano_DLBCL_v_BL.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$BLvDLBCL,
                lab = rownames(sig_herv$BLvDLBCL),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'DLBCL vs BL')
dev.off()


pdf("plots/05l-all_lymphoma_volcano_FL_v_BL.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$BLvFL,
                lab = rownames(sig_herv$BLvFL),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'FL vs BL')
dev.off()

pdf("plots/05l-all_lymphoma_volcano_FL_v_DLBCL.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$DLBCLvFL,
                lab = rownames(sig_herv$DLBCLvFL),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'FL vs DLBCL')
dev.off()

################################### SAVE DATA ##################################

upvars_all_l <- upvars
upvars_all_l_hervs <- upvars_hervs
downvars_all_l <- downvars
downvars_all_l_hervs <- downvars_hervs
sig_herv_all_l <- sig_herv
sig_all_l <- sig

save(upvars_all_l, upvars_all_l_hervs, downvars_all_l, downvars_all_l_hervs,
     sig_herv_all_l, sig_all_l, file = "r_outputs/05l-all_lymphoma_vars.Rdata")   

################################# CLUSTER SIZES ################################

load("r_outputs/02-all_lymphoma_filt_counts.Rdata")
reads.left <- as.data.frame(colSums(all.counts.filt.comb))
colnames(reads.left) <- c("reads")
reads.left$sample <- rownames(reads.left)

reads.left$clust <- all_metadata$cancer_type

reads.left %>%
  group_by(clust) %>%
  summarise_at(vars(reads), list(name = mean))

# A tibble: 3 × 2
# clust      name
# <chr>     <dbl>
#   1 BL    46223018.
# 2 DLBCL 46417908.
# 3 FL    35957779.

