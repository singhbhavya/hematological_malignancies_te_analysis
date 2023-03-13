################################################################################
################################################################################
################################################################################
################################################################################
######################## BULK AGIRRE GCB PCA & DESEQ2 ##########################

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

load("r_outputs/02-GCB_Agirre_filt_counts.Rdata")
load("r_outputs/01-metadata.Rdata")
load("r_outputs/01-refs.Rdata")

remove(BL_metadata, DLBCL_metadata, FL_metadata, all_metadata, bulk_metadata)

################################# METADATA SETUP ###############################
 
agirre_metadata$Cell <- agirre_metadata$source_name
agirre_metadata$Cell[agirre_metadata$Cell == "Naive"] <- "Naive B"
agirre_metadata$Cell[agirre_metadata$Cell == "Centroblast"] <- "Dark Zone Germinal Center B"
agirre_metadata$Cell[agirre_metadata$Cell == "Centrocyte"] <- "Light Zone Germinal Center B"
agirre_metadata$Cell[agirre_metadata$Cell == "Memory"] <- "Memory B"
agirre_metadata$Cell[agirre_metadata$Cell == "Tonsilar plasma cell"] <- "Plasmablasts"
agirre_metadata$Cell[agirre_metadata$Cell == "Bone Marrow plasma cell"] <- "Bone Marrow plasma cell"

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

############################### GCB AGIRRE DESEQ ###############################

### DESeq2 (HERVs Only)

bulk.countdat <- GCB_Agirre.filt.herv
cat(sprintf('%d variables\n', nrow(bulk.countdat)))

stopifnot(all(colnames(bulk.countdat) == rownames(agirre_metadata)))

Agirre.dds <- DESeq2::DESeqDataSetFromMatrix(countData = bulk.countdat,
                                          colData = agirre_metadata,
                                          design = ~ Cell + 0)

Agirre.dds <- DESeq2::DESeq(Agirre.dds, parallel=T)
Agirre.tform <- DESeq2::varianceStabilizingTransformation(Agirre.dds, blind=FALSE)

## PCA
Agirre.herv.pca.obj <-
  pca_standard(tform = Agirre.tform, 
               metadata = agirre_metadata, 
               var = 0.1)

# 4 PCs for Elbow method
# 6 PCs for Horn method
# 3 PCs needed to explain 50 percent of variation


######################## GCB BULK DESEQ (HERVS + GENES) ########################

### DESeq2 (HERVs + Genes)

bulk.g.countdat <- GCB_Agirre.filt.comb
cat(sprintf('%d variables\n', nrow(bulk.g.countdat)))

stopifnot(all(colnames(bulk.g.countdat) == rownames(agirre_metadata)))

Agirre.g.dds <- DESeq2::DESeqDataSetFromMatrix(countData = bulk.g.countdat,
                                            colData = agirre_metadata,
                                            design = ~ Cell + 0)

Agirre.g.dds <- DESeq2::DESeq(Agirre.g.dds, parallel=T)
Agirre.g.tform <- DESeq2::varianceStabilizingTransformation(Agirre.g.dds, blind=FALSE)

## PCA
Agirre.g.pca.obj <-
  pca_standard(tform = Agirre.g.tform, 
               metadata = agirre_metadata, 
               var = 0.1)

# 4 PCs for Elbow method
# 6 PCs for Horn method
# 2 PCs needed to explain 50 percent of variation

##################### AGIRRE GCB BIPLOTS HERVs & GENES #########################

## Biplot with projects (HERVs and genes)

pdf("plots/05f-gcb_agirre_pca_genes_hervs.pdf", height=5, width=6)
biplot(Agirre.g.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "Cell",
       colkey = c("Naive B" = pal_jco("default", alpha = 0.8)(7)[1],
                  "Memory B" = pal_jco("default", alpha = 0.8)(7)[3],
                  "Dark Zone Germinal Center B" = pal_jco("default", alpha = 0.8)(57)[4],
                  "Light Zone Germinal Center B" = pal_jco("default", alpha = 0.8)(7)[5],
                  "Plasmablasts" = pal_jco("default", alpha = 0.8)(7)[6],
                  "Bone Marrow plasma cell" = pal_jco("default", alpha = 0.8)(7)[7]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)
dev.off()

pdf("plots/05f-gcb_agirre_pca_hervs.pdf", height=5, width=6)
## Biplot with projects (HERVs only)
biplot(Agirre.herv.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "Cell",
       colkey = c("Naive B" = pal_jco("default", alpha = 0.8)(7)[1],
                  "Memory B" = pal_jco("default", alpha = 0.8)(7)[3],
                  "Dark Zone Germinal Center B" = pal_jco("default", alpha = 0.8)(57)[4],
                  "Light Zone Germinal Center B" = pal_jco("default", alpha = 0.8)(7)[5],
                  "Plasmablasts" = pal_jco("default", alpha = 0.8)(7)[6],
                  "Bone Marrow plasma cell" = pal_jco("default", alpha = 0.8)(7)[7]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)
dev.off()

################################ SET THRESHOLDS ################################

fc=2 # fold change
l2fc=log2(fc) # log 2 fold change
lfc.cutoff <- 1.5
pval=0.001 # p value threshold

############################## TOP GENES & HERVS ###############################

# resultsNames(Agirre.g.dds)
# [1] "CellBone.Marrow.plasma.cell" "CellDark.Zone.Germinal.Center.B"  "CellLight.Zone.Germinal.Center.B"
# [4] "CellMemory.B""CellNaive.B" "CellPlasmablasts"                              

agirre_res <- list(
  "BMPC" = DESeq2::results(Agirre.g.dds, contrast=c(+1, -1/5, -1/5, -1/5, -1/5, -1/5), alpha=pval),
  "DZ" = DESeq2::results(Agirre.g.dds, contrast=c(-1/5, +1, -1/5, -1/5, -1/5, -1/5), alpha=pval),
  "LZ" = DESeq2::results(Agirre.g.dds, contrast=c(-1/5, -1/5, +1, -1/5, -1/5, -1/5), alpha=pval),
  "MB" = DESeq2::results(Agirre.g.dds, contrast=c(-1/5, -1/5, -1/5, +1, -1/5, -1/5), alpha=pval),
  "NB" = DESeq2::results(Agirre.g.dds, contrast=c(-1/5, -1/5, -1/5, -1/5, +1, -1/5), alpha=pval),
  "PB" = DESeq2::results(Agirre.g.dds, contrast=c(-1/5, -1/5, -1/5, -1/5, -1/5, +1), alpha=pval),
  "DZvLZ" = DESeq2::results(Agirre.g.dds, contrast=c("Cell", "Dark Zone Germinal Center B", 
                                                  "Light Zone Germinal Center B"), alpha=pval),
  "MBvsNB" = DESeq2::results(Agirre.g.dds, contrast=c("Cell", "Memory B", "Naive B"), alpha=pval)
)

agirre_res <- lapply(agirre_res, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig <- lapply(agirre_res, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig[[n]])
}

############################### TOP HERVs ONLY #################################

# resultsNames(Agirre.dds)
# [1] "CellBone.Marrow.plasma.cell" "CellDark.Zone.Germinal.Center.B"  "CellLight.Zone.Germinal.Center.B"
# [4] "CellMemory.B""CellNaive.B" "CellPlasmablasts"                              

agirre_res_herv <- list(
  "BMPC" = DESeq2::results(Agirre.dds, contrast=c(+1, -1/5, -1/5, -1/5, -1/5, -1/5), alpha=pval),
  "DZ" = DESeq2::results(Agirre.dds, contrast=c(-1/5, +1, -1/5, -1/5, -1/5, -1/5), alpha=pval),
  "LZ" = DESeq2::results(Agirre.dds, contrast=c(-1/5, -1/5, +1, -1/5, -1/5, -1/5), alpha=pval),
  "MB" = DESeq2::results(Agirre.dds, contrast=c(-1/5, -1/5, -1/5, +1, -1/5, -1/5), alpha=pval),
  "NB" = DESeq2::results(Agirre.dds, contrast=c(-1/5, -1/5, -1/5, -1/5, +1, -1/5), alpha=pval),
  "PB" = DESeq2::results(Agirre.dds, contrast=c(-1/5, -1/5, -1/5, -1/5, -1/5, +1), alpha=pval),
  "DZvLZ" = DESeq2::results(Agirre.dds, contrast=c("Cell", "Dark Zone Germinal Center B", 
                                                     "Light Zone Germinal Center B"), alpha=pval),
  "MBvsNB" = DESeq2::results(Agirre.dds, contrast=c("Cell", "Memory B", "Naive B"), alpha=pval)
)

agirre_res_herv <- lapply(agirre_res_herv, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig_herv <- lapply(agirre_res_herv, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

sink(file = "r_outputs/05f-agirre_herv_deseq.txt")
for (n in names(sig_herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig_herv[[n]])
}
sink(file = NULL)

sink(file = "r_outputs/05f-agirre_herv_deseq_top5.txt")
for (n in names(sig_herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  print(head(sig_herv[[n]]))
}
sink(file = NULL)

############################ UPSET GENES & HERVS ###############################

upvars_agirre <- lapply(sig[1:6], function(r) rownames(subset(r, log2FoldChange>0)))
downvars_agirre <- lapply(sig[1:6], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05f-gcb_agirre_upvars.pdf", height=5, width=7)
upset(fromList(upvars_agirre), sets=c("BMPC", "DZ", "LZ", "MB", "NB", "PB"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

pdf("plots/05f-gcb_agirre_upset_dnwars.pdf", height=5, width=7)
upset(fromList(downvars_agirre), sets=c("BMPC", "DZ", "LZ", "MB", "NB", "PB"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

up.binmat.agirre <- fromList(upvars_agirre)
rn <- do.call(c, upvars_agirre)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.agirre) <- rn
rm(rn)

dn.binmat.agirre <- fromList(downvars_agirre)
rn <- do.call(c, downvars_agirre)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.agirre) <- rn
rm(rn)

############################# UPSET HERVs ONLY #################################

upvars_agirre_hervs <- lapply(sig_herv[1:6], function(r) rownames(subset(r, log2FoldChange>0)))
downvars_agirre_hervs <- lapply(sig_herv[1:6], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05f-gcb_agirre_upset_upvars_hervs.pdf", height=5, width=7)
upset(fromList(upvars_agirre_hervs), sets=c("BMPC", "DZ", "LZ", "MB", "NB", "PB"),  
      keep.order = T, order.by='degree', decreasing=F)
dev.off()

pdf("plots/05f-gcb_agirre_upset_dnvars_hervs.pdf", height=5, width=7)
upset(fromList(downvars_agirre_hervs), sets=c("BMPC", "DZ", "LZ", "MB", "NB", "PB"),  
      keep.order = T, order.by='degree', decreasing=F)
dev.off()

up.binmat.hervs.agirre <- fromList(upvars_agirre_hervs)
rn <- do.call(c, upvars_agirre_hervs)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.hervs.agirre) <- rn
rm(rn)

dn.binmat.hervs.agirre <- fromList(downvars_agirre_hervs)
rn <- do.call(c, downvars_agirre_hervs)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.hervs.agirre) <- rn
rm(rn)

################################# HEATMAPS #####################################

################################### SETUP ######################################

# Colors for cells
cols <- rgb_gsea(palette = c("default"), n = 14, alpha = 0.7, reverse = FALSE)

# Annotation column
df <- as.data.frame(colData(Agirre.dds)[,c("BioSample","Cell")])
df <- subset(df, select = -c(1))

# Create colors for each group
annoCol <-  c("#0073C2CC", "#868686CC", "#CD534CCC", "#7AA6DCCC",
                 "#003C67CC", "#8F7700CC")
names(annoCol) <- c("Naive B", "Memory B", "Dark Zone Germinal Center B",
                    "Light Zone Germinal Center B", "Plasmablasts",
                    "Bone Marrow plasma cell")
annoCol <- list(Cell = annoCol)

######################### UPREGULATED IN ALL GROUPS ############################

top.genes.hervs <- rownames(up.binmat.agirre)
top.hervs <- rownames(up.binmat.hervs.agirre)

pdf("plots/05f-gcb_agirre_top_hervs_upregulated_all.pdf", height=10, width=10)
pheatmap(assay(Agirre.tform)[top.hervs,], 
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

pdf("plots/05f-gcb_agirre_top_hervs_genes_upregulated_all.pdf", height=10, width=10)
pheatmap(assay(Agirre.g.tform)[top.genes.hervs,],
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
dev.off()

########################### DE HERVs PER GC SITE ###############################

makeheatmap <- function(topgenes, ...) {
  args <- list(...)
  mat <- assay(Agirre.tform)[unique(topgenes), ]
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


for(clust in c("BMPC", "DZ", "LZ", "NB", "PB")) {
  tg <- rownames(sig_herv[[clust]][1:75,])
  p <- makeheatmap(tg, main=paste0('DE in cluster ', clust))
  pdf(paste0("plots/05f-gcb_agirre_top_de_hervs_", clust, ".pdf"), height=7, width=7)
  print(p)
  dev.off()
}

for(clust in c("MB")) {
  tg <- rownames(sig_herv[[clust]][1:20,])
  p <- makeheatmap(tg, main=paste0('DE in cluster ', clust))
  pdf(paste0("plots/05f-gcb_agirre_top_de_hervs_", clust, ".pdf"), height=5, width=7)
  print(p)
  dev.off()
}

############################# VOLCANO DZ VS LZ #################################

pdf("plots/05f-gcb_agirre_volcano_DZ_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$DZ,
                lab = rownames(sig_herv$DZ),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs DZ')
dev.off()

pdf("plots/05f-gcb_agirre_volcano_LZ_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$LZ,
                lab = rownames(sig_herv$LZ),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs LZ')
dev.off()

pdf("plots/05f-gcb_agirre_volcano_NB_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$NB,
                lab = rownames(sig_herv$NB),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs NB')
dev.off()

pdf("plots/05f-gcb_agirre_volcano_MB_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$MB,
                lab = rownames(sig_herv$MB),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs MB')
dev.off()

pdf("plots/05f-gcb_agirre_volcano_PB_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$PB,
                lab = rownames(sig_herv$PB),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs PB')
dev.off()

pdf("plots/05f-gcb_agirre_volcano_BMPC_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$BMPC,
                lab = rownames(sig_herv$BMPC),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs BMPC')
dev.off()

pdf("plots/05f-gcb_agirre_volcano_DZ_v_LZ.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$DZvLZ,
                lab = rownames(sig_herv$DZvLZ),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'LZ vs DZ')
dev.off()

pdf("plots/05f-gcb_agirre_volcano_NB_v_MB.pdf", height=8, width=8)
EnhancedVolcano(sig_herv$MBvsNB,
                lab = rownames(sig_herv$MBvsNB),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'NB vs MB')
dev.off()

########################## PLOT INDIVIDUAL HERVs ###############################

plot.counts <- function(df, gene) {
  
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="Cell", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=Cell, y=count, fill=Cell))  +
    geom_boxplot() +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("Cell Type") +
    ylab("Counts") +
    scale_x_discrete(labels=c("Dark Zone Germinal Center B" = "DZ",
                              "Bone Marrow plasma cell" = "BMPC",
                              "Light Zone Germinal Center B" = "LZ",
                              "Memory B" = "MB",
                              "Naive B" = "NB",
                              "Plasmablasts" = "PB")) +
    scale_fill_manual(values = c("Naive B" = pal_jco("default", alpha = 0.8)(7)[1],
                                 "Memory B" = pal_jco("default", alpha = 0.8)(7)[3],
                                 "Dark Zone Germinal Center B" = pal_jco("default", alpha = 0.8)(7)[4],
                                 "Light Zone Germinal Center B" = pal_jco("default", alpha = 0.8)(7)[5],
                                 "Plasmablasts" = pal_jco("default", alpha = 0.8)(7)[6],
                                 "Bone Marrow plasma cell" = pal_jco("default", alpha = 0.8)(7)[7])) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    ggtitle(gene) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)
}


p1 <- plot.counts(Agirre.g.dds, "HML2_3q12.3")
p2 <- plot.counts(Agirre.g.dds, "ERV316A3_6p21.33c")
p3 <- plot.counts(Agirre.g.dds, "MER4_17q21.2d")
p4 <- plot.counts(Agirre.g.dds, "HML3_5p15.33d")
p5 <- plot.counts(Agirre.g.dds, "HERVEA_5q22.2")
p6 <- plot.counts(Agirre.g.dds, "HERVH_7q11.21")
p7 <- plot.counts(Agirre.g.dds, "HML2_1q22")
p8 <- plot.counts(Agirre.g.dds, "HERVIP10F_2q21.2")
p9 <- plot.counts(Agirre.g.dds, "MER61_12q15a")

pdf("plots/05f-gcb_agirre_individual_hervs.pdf", height=12, width=12)
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9,
          nrow = 3, 
          ncol = 3,
          labels = "AUTO")
dev.off()

pdf("plots/05f-agirre_gcb_healthy_HARLEQUIN_17q21.31.pdf", height=3, width=3)
plot.counts(Agirre.g.dds, "HARLEQUIN_17q21.31")
dev.off()

pdf("plots/05f-agirre_gcb_healthy_HARLEQUIN_1q32.1.pdf", height=3, width=3)
plot.counts(Agirre.g.dds, "HARLEQUIN_1q32.1")
dev.off()


pdf("plots/05f-agirre_gcb_healthy_LILRB1.pdf", height=3, width=3)
plot.counts(Agirre.g.dds, "ENSG00000104972.16")
dev.off()

pdf("plots/05f-agirre_gcb_healthy_LILRB1.pdf", height=3, width=3)
plot.counts(Agirre.g.dds, "ENSG00000104972.16")
dev.off()

pdf("plots/05f-agirre_gcb_healthy_MER101_16p12.2a.pdf", height=3, width=3)
plot.counts(Agirre.g.dds, "MER101_16p12.2a")
dev.off()

################################### SAVE DATA ##################################


sig_herv_agirre <- sig_herv
sig_agirre <- sig

save(upvars_agirre, upvars_agirre_hervs, downvars_agirre, downvars_agirre_hervs,
     sig_herv_agirre, sig_agirre, file = "r_outputs/05f-agirre_vars.Rdata")

