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
library(ComplexUpset)
library(EnhancedVolcano)
library(pheatmap)
library(fgsea)
library(viridis)

################################### LOAD DATA ##################################

load("r_outputs/02-DLBCL_filt_counts.Rdata")
load("r_outputs/05o-DLBCL_pca_ccp_clusters_metadata.Rdata")
load("r_outputs/01-refs.Rdata")

# Add clusters to metadata
DLBCL_metadata$clust.retro.k2 <- clust.df$clust.retro.k2
DLBCL_metadata$clust.retro.k3 <- clust.df$clust.retro.k3
DLBCL_metadata$clust.retro.k4 <- clust.df$clust.retro.k4
DLBCL_metadata$clust.retro.k5 <- clust.df$clust.retro.k5
DLBCL_metadata$clust.retro.k7 <- clust.df$clust.retro.k7
DLBCL_metadata$clust.retro.k9 <- clust.df$clust.retro.k9

################################################################################
################################################################################
#################################### CLUST 7 ###################################
################################################################################
################################################################################


############################# FUNCTION SCREE PLOT ##############################

pdf("plots/05o-dlbcl_subtypes_biplot_clust.k7.pdf", height=5, width=6)
biplot(DLBCL.herv.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       pointSize = 2, 
       ellipse = FALSE,
       colby = "clust.retro.k7",
       shape = "project", 
       shapekey = c("NCICCR-DLBCL" = 15, "TCGA-DLBC" = 8),
       colkey = ggsci::pal_npg(palette = c("nrc"))(9),
       legendPosition = "right") +
  theme_cowplot() +
  theme(aspect.ratio = 1)

dev.off()


#################################### DESEQ2 ####################################

# p value threshold 
lfc.cutoff <- 1.5
pval=0.001 

# Sanity check
stopifnot(all(colnames(DLBCL.filt.comb) == rownames(DLBCL_metadata)))


# Run Deseq2
DLBCL.k7.dds <- DESeq2::DESeqDataSetFromMatrix(DLBCL.filt.comb, 
                                               DLBCL_metadata, 
                                            ~ clust.retro.k7 + 0)


DLBCL.k7.dds <- DESeq2::DESeq(DLBCL.k7.dds, parallel=T)
DLBCL.k7.tform <- DESeq2::varianceStabilizingTransformation(DLBCL.k7.dds, blind=FALSE)

################################# DESEQ2 HERVs ################################# 

stopifnot(all(colnames(DLBCL.filt.herv) == rownames(DLBCL_metadata)))

DLBCL.k7.herv.dds <- DESeq2::DESeqDataSetFromMatrix(DLBCL.filt.herv, 
                                                    DLBCL_metadata, 
                                                 ~ clust.retro.k7 + 0)


DLBCL.k7.herv.dds <- DESeq2::DESeq(DLBCL.k7.herv.dds, parallel=T)
DLBCL.k7.herv.tform <- DESeq2::varianceStabilizingTransformation(DLBCL.k7.herv.dds, 
                                                              blind=FALSE)


############################## TOP GENES AND HERVS ############################# 

res.k7 <- list(
  "C1" = DESeq2::results(DLBCL.k7.dds, contrast=c(+1, -1/6, -1/6, -1/6, -1/6, -1/6, -1/6), alpha=pval),
  "C2" = DESeq2::results(DLBCL.k7.dds, contrast=c(-1/6, +1, -1/6, -1/6, -1/6, -1/6, -1/6), alpha=pval),
  "C3" = DESeq2::results(DLBCL.k7.dds, contrast=c(-1/6, -1/6, +1, -1/6, -1/6, -1/6, -1/6), alpha=pval),
  "C4" = DESeq2::results(DLBCL.k7.dds, contrast=c(-1/6, -1/6, -1/6, +1, -1/6, -1/6, -1/6), alpha=pval),
  "C5" = DESeq2::results(DLBCL.k7.dds, contrast=c(-1/6, -1/6, -1/6, -1/6, +1, -1/6, -1/6), alpha=pval),
  "C6" = DESeq2::results(DLBCL.k7.dds, contrast=c(-1/6, -1/6, -1/6, -1/6, -1/6, +1, -1/6), alpha=pval),
  "C7" = DESeq2::results(DLBCL.k7.dds, contrast=c(-1/6, -1/6, -1/6, -1/6, -1/6, -1/6, +1), alpha=pval)
)

res.k7 <- lapply(res.k7, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})


sig.k7 <- lapply(res.k7, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig.k7)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig.k7[[n]])
}

for (n in names(sig.k7)) {
  cat("\n#--- Contrast", n, "---#\n")
  print(head(sig.k7[[n]]))
}

################################ TOP HERVS ONLY ################################ 

res.k7.herv <- list(
  "C1" = DESeq2::results(DLBCL.k7.herv.dds, contrast=c(+1, -1/6, -1/6, -1/6, -1/6, -1/6, -1/6), alpha=pval),
  "C2" = DESeq2::results(DLBCL.k7.herv.dds, contrast=c(-1/6, +1, -1/6, -1/6, -1/6, -1/6, -1/6), alpha=pval),
  "C3" = DESeq2::results(DLBCL.k7.herv.dds, contrast=c(-1/6, -1/6, +1, -1/6, -1/6, -1/6, -1/6), alpha=pval),
  "C4" = DESeq2::results(DLBCL.k7.herv.dds, contrast=c(-1/6, -1/6, -1/6, +1, -1/6, -1/6, -1/6), alpha=pval),
  "C5" = DESeq2::results(DLBCL.k7.herv.dds, contrast=c(-1/6, -1/6, -1/6, -1/6, +1, -1/6, -1/6), alpha=pval),
  "C6" = DESeq2::results(DLBCL.k7.herv.dds, contrast=c(-1/6, -1/6, -1/6, -1/6, -1/6, +1, -1/6), alpha=pval),
  "C7" = DESeq2::results(DLBCL.k7.herv.dds, contrast=c(-1/6, -1/6, -1/6, -1/6, -1/6, -1/6, +1), alpha=pval)
)

res.k7.herv <- lapply(res.k7.herv, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})


sig.k7.herv <- lapply(res.k7.herv, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig.k7.herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig.k7.herv[[n]])
}

for (n in names(sig.k7.herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  print(head(sig.k7.herv[[n]]))
}


############################ UPSET GENES & HERVS ###############################

upvars.DLBCL.k7 <- lapply(sig.k7[1:7], function(r) rownames(subset(r, log2FoldChange>0)))
downvars.DLBCL.k7 <- lapply(sig.k7[1:7], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05q-DLBCL_k7_upset_upvars.pdf", height=5, width=7)

ComplexUpset::upset(fromList(upvars.DLBCL.k7), 
                    intersect = c(names(upvars.DLBCL.k7)),
                    intersections = list( 
                      c("C1"), 
                      c("C2"), 
                      c("C3"),
                      c("C4"), 
                      c("C5"),
                      c("C6"),
                      c("C7")), 
                    queries = list(
                      upset_query(set=c("C1"), 
                                  color = "#E64B35FF", 
                                  fill = "#E64B35FF"),
                      upset_query(set=c("C2"), 
                                  color = "#4DBBD5FF", 
                                  fill = "#4DBBD5FF"),
                      upset_query(set=c("C3"), 
                                  color = "#00A087FF", 
                                  fill = "#00A087FF"),
                      upset_query(set=c("C4"), 
                                  color = "#3C5488FF", 
                                  fill = "#3C5488FF"),
                      upset_query(set=c("C5"), 
                                  color = "#F39B7FFF", 
                                  fill = "#F39B7FFF"),
                      upset_query(set=c("C6"), 
                                  color = "#8491B4FF", 
                                  fill = "#8491B4FF"),
                      upset_query(set=c("C7"), 
                                  color = "#91D1C2FF", 
                                  fill = "#91D1C2FF")
                    ))

dev.off()

pdf("plots/05q-DLBCL_k7_upset_downvars.pdf", height=5, width=7)
ComplexUpset::upset(fromList(downvars.DLBCL.k7), 
                    intersect = c(names(upvars.DLBCL.k7)),
                    intersections = list( 
                      c("C1"), 
                      c("C2"), 
                      c("C3"),
                      c("C4"), 
                      c("C5"),
                      c("C6"),
                      c("C7")), 
                    queries = list(
                      upset_query(set=c("C1"), 
                                  color = "#E64B35FF", 
                                  fill = "#E64B35FF"),
                      upset_query(set=c("C2"), 
                                  color = "#4DBBD5FF", 
                                  fill = "#4DBBD5FF"),
                      upset_query(set=c("C3"), 
                                  color = "#00A087FF", 
                                  fill = "#00A087FF"),
                      upset_query(set=c("C4"), 
                                  color = "#3C5488FF", 
                                  fill = "#3C5488FF"),
                      upset_query(set=c("C5"), 
                                  color = "#F39B7FFF", 
                                  fill = "#F39B7FFF"),
                      upset_query(set=c("C6"), 
                                  color = "#8491B4FF", 
                                  fill = "#8491B4FF"),
                      upset_query(set=c("C7"), 
                                  color = "#91D1C2FF", 
                                  fill = "#91D1C2FF")
                    ))
dev.off()

up.binmat.DLBCL.k7 <- fromList(upvars.DLBCL.k7)
rn <- do.call(c, upvars.DLBCL.k7)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.DLBCL.k7) <- rn
rm(rn)

dn.binmat.DLBCL.k7 <- fromList(downvars.DLBCL.k7)
rn <- do.call(c, downvars.DLBCL.k7)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.DLBCL.k7) <- rn
rm(rn)

############################## UPSET HERVS ONLY ################################

upvars.DLBCL.k7.herv <- lapply(sig.k7.herv[1:7], function(r) rownames(subset(r, log2FoldChange>0)))
downvars.DLBCL.k7.herv <- lapply(sig.k7.herv[1:7], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05q-DLBCL_k7_upset_upvars_hervs.pdf", height=5, width=7)

ComplexUpset::upset(fromList(upvars.DLBCL.k7.herv), 
                    intersect = c(names(upvars.DLBCL.k7.herv)),
                    intersections = list( 
                      c("C1"), 
                      c("C2"), 
                      c("C3"),
                      c("C4"), 
                      c("C5"),
                      c("C6"),
                      c("C7")), 
                    queries = list(
                      upset_query(set=c("C1"), 
                                  color = "#E64B35FF", 
                                  fill = "#E64B35FF"),
                      upset_query(set=c("C2"), 
                                  color = "#4DBBD5FF", 
                                  fill = "#4DBBD5FF"),
                      upset_query(set=c("C3"), 
                                  color = "#00A087FF", 
                                  fill = "#00A087FF"),
                      upset_query(set=c("C4"), 
                                  color = "#3C5488FF", 
                                  fill = "#3C5488FF"),
                      upset_query(set=c("C5"), 
                                  color = "#F39B7FFF", 
                                  fill = "#F39B7FFF"),
                      upset_query(set=c("C6"), 
                                  color = "#8491B4FF", 
                                  fill = "#8491B4FF"),
                      upset_query(set=c("C7"), 
                                  color = "#91D1C2FF", 
                                  fill = "#91D1C2FF")
                    ))

dev.off()

pdf("plots/05q-DLBCL_k7_upset_downvars_hervs.pdf", height=5, width=7)
ComplexUpset::upset(fromList(downvars.DLBCL.k7.herv), 
                    intersect = c(names(upvars.DLBCL.k7)),
                    intersections = list( 
                      c("C1"), 
                      c("C2"), 
                      c("C3"),
                      c("C4"), 
                      c("C5"),
                      c("C6"),
                      c("C7")), 
                    queries = list(
                      upset_query(set=c("C1"), 
                                  color = "#E64B35FF", 
                                  fill = "#E64B35FF"),
                      upset_query(set=c("C2"), 
                                  color = "#4DBBD5FF", 
                                  fill = "#4DBBD5FF"),
                      upset_query(set=c("C3"), 
                                  color = "#00A087FF", 
                                  fill = "#00A087FF"),
                      upset_query(set=c("C4"), 
                                  color = "#3C5488FF", 
                                  fill = "#3C5488FF"),
                      upset_query(set=c("C5"), 
                                  color = "#F39B7FFF", 
                                  fill = "#F39B7FFF"),
                      upset_query(set=c("C6"), 
                                  color = "#8491B4FF", 
                                  fill = "#8491B4FF"),
                      upset_query(set=c("C7"), 
                                  color = "#91D1C2FF", 
                                  fill = "#91D1C2FF")
                    ))
dev.off()

up.binmat.DLBCL.k7.hervs <- fromList(upvars.DLBCL.k7.herv)
rn <- do.call(c, upvars.DLBCL.k7.herv)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.DLBCL.k7.hervs) <- rn
rm(rn)

dn.binmat.DLBCL.k7.hervs <- fromList(downvars.DLBCL.k7.herv)
rn <- do.call(c, downvars.DLBCL.k7.herv)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.DLBCL.k7.hervs) <- rn
rm(rn)

################################# HEATMAPS #####################################
################################### SETUP ######################################

# Colors for cells
cols <- rgb_gsea(palette = c("default"), n = 14, alpha = 0.7, reverse = FALSE)

# Annotation column
df <- as.data.frame(colData(DLBCL.k7.dds)[,c("COO_class", "LymphGen_call",
                                             "EcoTyper_call", "clust.retro.k7")])

# Create colors for each group
annoCol <- list(COO_class = c("GCB" = "royalblue", 
                                "ABC" = "red3", 
                                "Unclass" = "lightblue", 
                                "Missing" = "grey"),
                EcoTyper_call = c("S1" = "#bcc779",
                                     "S2" = "#60aaec",
                                     "S3" = "#Ecaf60",
                                     "S4" = "#B789EE",
                                     "S5" = "#E0AF71",
                                     "Unassigned" = "grey",
                                     "Missing" = "grey"),
                LymphGen_call = c("A53" = "#e78bf0",
                                     "BN2" = "#e28743",
                                     "BN2/A53" = "#AF6C3A",
                                     "BN2/EZB" = "#c79875",
                                     "BN2/MCD" = "#7a4c29",
                                     "BN2/ST2" = "#58361d",
                                     "EZB" = "#3d3aaf",
                                     "EZB/A53" = "#7775c7",
                                     "EZB/MCD" = "#8b89cf",
                                     "EZB/N1/ST2/A53" = "#b1b0df",
                                     "EZB/ST2/A53" = "#d8d8ef",
                                     "MCD" = "#4d8347",
                                     "MCD/A53" = "#82a87e",
                                     "MCD/ST2" = "#b8cdb5",
                                     "N1" = "#e84646",
                                     "N1/ST2" = "#f19090",
                                     "Other" = "#e9f190",
                                     "ST2" = "#38baaa",
                                     "NA" = "grey"),
                clust.retro.k7 = c("C1" = "#E64B35FF",
                                    "C2" = "#4DBBD5FF",
                                    "C3" = "#00A087FF",
                                    "C4" = "#3C5488FF",
                                    "C5" = "#F39B7FFF",
                                    "C6" = "#8491B4FF",
                                    "C7" = "#91D1C2FF")
)

######################### UPREGULATED IN ALL GROUPS ############################
top.genes <- rownames(up.binmat.DLBCL.k7)

pdf("plots/05q-DLBCL_k7_geneshervs_upregulated_all.pdf", height=10, width=10)
pheatmap(assay(DLBCL.k7.dds)[top.genes,], 
         main="Upregulated Genes and HERVs, all clusters",
         cluster_rows=TRUE,
         show_rownames=FALSE,
         show_colnames = FALSE,
         color = cols,
         scale="row",
         breaks=seq(-3,3,length.out=14),
         labels_row = gene_table[top.genes,]$display,
         cluster_cols=TRUE, 
         treeheight_row=0,
         annotation_col=df,
         annotation_colors = annoCol)
dev.off()

top.hervs <- rownames(up.binmat.DLBCL.k7.hervs)

pdf("plots/05q-DLBCL_k7_hervs_upregulated_all.pdf", height=10, width=10)
pheatmap(assay(DLBCL.k7.herv.dds)[top.hervs,], 
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
  mat <- assay(DLBCL.k7.tform)[unique(topgenes), ]
  rowdist <- as.dist((1 - cor(t(mat), method='spearman'))/2)
  
  annotation_col <- df
  cols <- rgb_gsea(palette = c("default"), n = 14, alpha = 0.7, reverse = FALSE)
  
  if (topgenes[1] %in% retro.annot.v2$Locus) {
    annoRow <- as.data.frame(retro.annot.v2[,c("TE_type", "Locus")])
    annoRow <- annoRow[topgenes,]
    annoRow <- subset(annoRow, select = -c(2))
    annotation_row <- annoRow
  } else {
    annotation_row <- NULL
  }
  
  pheatmap(mat,
           color=cols,
           scale="row", breaks=seq(-3,3,length.out=14),
           clustering_distance_rows = rowdist,
           clustering_method="average",
           annotation_col = df,
           annotation_colors = annoCol,
           labels_row = gene_table[topgenes,]$display,
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

for(clust in names(sig.k7.herv)) {
  tg <- rownames(sig.k7.herv[[clust]][1:75,])
  p <- makeheatmap(tg, main=paste0('DE in cluster ', clust))
  pdf(paste0("plots/05p-all_lymphoma_subtype_top_de_hervs_", clust, ".pdf"), height=7, width=7)
  print(p)
  dev.off()
}

for(clust in names(sig.k7)) {
  tg <- rownames(sig.k7[[clust]][1:75,])
  p <- makeheatmap(tg, main=paste0('DE in cluster ', clust))
  pdf(paste0("plots/05p-all_lymphoma_subtype_top_de_geeneshervs_", clust, ".pdf"), height=7, width=7)
  print(p)
  dev.off()
}

################################### PATHWAYS ###################################

pathways.hallmark <- gmtPathways("gsea/h.all.v2023.1.Hs.symbols.gmt")
pathways.immune <- gmtPathways("gsea/c7.immunesigdb.v2023.1.Hs.symbols.gmt")
pathways.bp <- gmtPathways("gsea/c5.go.bp.v2023.1.Hs.symbols.gmt")
pathways.kegg <- gmtPathways("gsea/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt")
pathways.biocarta <- gmtPathways("gsea/c2.cp.biocarta.v2023.1.Hs.symbols.gmt.txt")

################################ HALLMARK FGSEA ################################

make.fsgsea <- function(pathway, fgsea.res, clust_name, pathway_name) {
  
  fgsea.res$SYMBOL <- gene_table[rownames(fgsea.res),]$display
  
  fgsea.res <- as.data.frame(fgsea.res) %>% 
    dplyr::select(SYMBOL, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(stat))
  
  fgsea.ranks <- deframe(fgsea.res)
  
  fgsea.out <- fgsea(pathways=pathway, 
                     stats=fgsea.ranks, 
                     nPermSimple = 10000,
                     eps=0)
  
  k4.fgseaResTidy <- k4.fgsea %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  return(fgsea.out)
  
  assign(paste0(clust_name, ".", pathway_name, ".fgsea.out"), fgsea.out, envir = .GlobalEnv )
  }


fsgsea.hallmarks.k7 <- list(
  "C1.hallmark" = make.fsgsea(pathways.hallmark, res.k7$C1, "C1", "hallmark"),
  "C2.hallmark" = make.fsgsea(pathways.hallmark, res.k7$C2, "C2", "hallmark"),
  "C3.hallmark" = make.fsgsea(pathways.hallmark, res.k7$C3, "C3", "hallmark"),
  "C4.hallmark" = make.fsgsea(pathways.hallmark, res.k7$C4, "C4", "hallmark"),
  "C5.hallmark" = make.fsgsea(pathways.hallmark, res.k7$C5, "C5", "hallmark"),
  "C6.hallmark" = make.fsgsea(pathways.hallmark, res.k7$C6, "C6", "hallmark"),
  "C7.hallmark" = make.fsgsea(pathways.hallmark, res.k7$C7, "C7", "hallmark")
)

pdf("plots/05q-DLBCL_k7_all_c_hallmarks.pdf", height=10, width=7)
for(clust in names(fsgsea.hallmarks.k7)) {
  fgseaResTidy <- fsgsea.hallmarks.k7[[clust]] %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  rownames(fgseaResTidy) <- fgseaResTidy$pathway
  
  p <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=clust) + 
    theme_minimal()
  
  print(p)
}
dev.off()

fsgsea.hallmarks.k7.summ <- as.data.frame(do.call(cbind, 
                                                  lapply(fsgsea.hallmarks.k7, 
                                           function(x) x[, c("NES")])))

rownames(fsgsea.hallmarks.k7.summ) <- fsgsea.hallmarks.k7$C1.hallmark$pathway
rownames(fsgsea.hallmarks.k7.summ) <- gsub("HALLMARK_","",rownames(fsgsea.hallmarks.k7.summ))
colnames(fsgsea.hallmarks.k7.summ) <- gsub(".hallmark.NES","",colnames(fsgsea.hallmarks.k7.summ))

pdf("plots/05q-DLBCL_k7_all_hallmarks_heatmap.pdf", height=10, width=5)
pheatmap(fsgsea.hallmarks.k7.summ, 
         cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = TRUE,
         color = viridis_pal()(10),
         cluster_cols=TRUE, 
         treeheight_row=0)
dev.off()
  
################################# IMMUNE FGSEA #################################

fsgsea.immune.k7 <- list(
  "C1" = make.fsgsea(pathways.immune, res.k7$C1, "C1", "immune"),
  "C2" = make.fsgsea(pathways.immune, res.k7$C2, "C2", "immune"),
  "C3" = make.fsgsea(pathways.immune, res.k7$C3, "C3", "immune"),
  "C4" = make.fsgsea(pathways.immune, res.k7$C4, "C4", "immune"),
  "C5" = make.fsgsea(pathways.immune, res.k7$C5, "C5", "immune"),
  "C6" = make.fsgsea(pathways.immune, res.k7$C6, "C6", "immune"),
  "C7" = make.fsgsea(pathways.immune, res.k7$C7, "C7", "immune")
)

pdf("plots/05q-DLBCL_k7_all_c_immune.pdf", height=10, width=7)
for(clust in names(fsgsea.immune.k7)) {
  fgseaResTidy.immune <- fsgsea.immune.k7[[clust]] %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  rownames(fgseaResTidy.immune) <- fgseaResTidy.immune$pathway
  fgseaResTidy.immune <- fgseaResTidy.immune[fgseaResTidy.immune$pathway %like% "GC", ]
  
  p <- ggplot(fgseaResTidy.immune, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=clust) + 
    theme_minimal()
  
  print(p)
}
dev.off()

fgseaResTidy.immune.summ <- as.data.frame(do.call(cbind, 
                                                  lapply(fsgsea.immune.k7, 
                                                         function(x) x[, c("NES")])))
rownames(fgseaResTidy.immune.summ) <- fsgsea.immune.k7$C1$pathway

fgseaResTidy.immune.summ <- fgseaResTidy.immune.summ[rownames(fgseaResTidy.immune.summ) %like% "GC", ]

colnames(fgseaResTidy.immune.summ) <- gsub(".NES","",colnames(fgseaResTidy.immune.summ))

pdf("plots/05q-DLBCL_k7_gc_immune_heatmap.pdf", height=10, width=5)
pheatmap(fgseaResTidy.immune.summ, 
         cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = TRUE,
         color = viridis_pal()(10),
         cluster_cols=TRUE, 
         treeheight_row=0)
dev.off()

################################## KEGG FGSEA ##################################

fsgsea.kegg.k7 <- list(
  "C1" = make.fsgsea(pathways.kegg, res.k7$C1, "C1", "KEGG"),
  "C2" = make.fsgsea(pathways.kegg, res.k7$C2, "C2", "KEGG"),
  "C3" = make.fsgsea(pathways.kegg, res.k7$C3, "C3", "KEGG"),
  "C4" = make.fsgsea(pathways.kegg, res.k7$C4, "C4", "KEGG"),
  "C5" = make.fsgsea(pathways.kegg, res.k7$C5, "C5", "KEGG"),
  "C6" = make.fsgsea(pathways.kegg, res.k7$C6, "C6", "KEGG"),
  "C7" = make.fsgsea(pathways.kegg, res.k7$C7, "C7", "KEGG")
)

pdf("plots/05q-DLBCL_k7_all_c_kegg.pdf", height=20, width=7)
for(clust in names(fsgsea.kegg.k7)) {
  fgseaResTidy.kegg <- fsgsea.kegg.k7[[clust]] %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  rownames(fsgsea.kegg.k7[[clust]]) <- fsgsea.kegg.k7[[clust]]$pathway
  
  rownames(fgseaResTidy.kegg) <- fgseaResTidy.kegg$pathway
  
  p <- ggplot(fgseaResTidy.kegg, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=clust) + 
    theme_minimal()
  
  print(p)
  
  }
dev.off()


fgseaResTidy.kegg.summ <- as.data.frame(do.call(cbind, 
                                                  lapply(fsgsea.kegg.k7, 
                                                         function(x) x[, c("NES")])))
rownames(fgseaResTidy.kegg.summ) <- fsgsea.kegg.k7$C1$pathway

colnames(fgseaResTidy.kegg.summ) <- gsub(".NES","",colnames(fgseaResTidy.kegg.summ))

pdf("plots/05q-DLBCL_k7_kegg_heatmap.pdf", height=35, width=10)
pheatmap(fgseaResTidy.kegg.summ, 
         cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = TRUE,
         color = viridis_pal()(10),
         cluster_cols=TRUE, 
         treeheight_row=0)
dev.off()

################################ BIOCARTA FGSEA ################################

fsgsea.biocarta.k7 <- list(
  "C1" = make.fsgsea(pathways.biocarta, res.k7$C1, "C1", "Biocarta"),
  "C2" = make.fsgsea(pathways.biocarta, res.k7$C2, "C2", "Biocarta"),
  "C3" = make.fsgsea(pathways.biocarta, res.k7$C3, "C3", "Biocarta"),
  "C4" = make.fsgsea(pathways.biocarta, res.k7$C4, "C4", "Biocarta"),
  "C5" = make.fsgsea(pathways.biocarta, res.k7$C5, "C5", "Biocarta"),
  "C6" = make.fsgsea(pathways.biocarta, res.k7$C6, "C6", "Biocarta"),
  "C7" = make.fsgsea(pathways.biocarta, res.k7$C7, "C7", "Biocarta")
)

pdf("plots/05q-DLBCL_k7_all_c_biocarta.pdf", height=30, width=7)
for(clust in names(fsgsea.biocarta.k7)) {
  fgseaResTidy.biocarta <- fsgsea.biocarta.k7[[clust]] %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  rownames(fsgsea.biocarta.k7[[clust]]) <- fsgsea.biocarta.k7[[clust]]$pathway
  
  rownames(fgseaResTidy.biocarta) <- fgseaResTidy.biocarta$pathway
  
  p <- ggplot(fgseaResTidy.biocarta, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=clust) + 
    theme_minimal()
  
  print(p)
  
}
dev.off()

for(clust in names(fsgsea.biocarta.k7)) {
  
  fsgsea.biocarta.k7[[clust]]$NES[fsgsea.biocarta.k7[[clust]]$padj >= 0.05] <- NA 
}

fgseaResTidy.biocarta.summ <- as.data.frame(do.call(cbind, 
                                                lapply(fsgsea.biocarta.k7, 
                                                       function(x) x[, c("NES")])))
rownames(fgseaResTidy.biocarta.summ) <- fsgsea.biocarta.k7$C1$pathway

colnames(fgseaResTidy.biocarta.summ) <- gsub(".NES","",colnames(fgseaResTidy.biocarta.summ))
fgseaResTidy.biocarta.summ <- 
  fgseaResTidy.biocarta.summ[rowSums(is.na(fgseaResTidy.biocarta.summ)) != ncol(fgseaResTidy.biocarta.summ), ]

pdf("plots/05q-DLBCL_k7_biocarta_heatmap.pdf", height=40, width=10)
pheatmap(data.matrix(fgseaResTidy.biocarta.summ), 
         cluster_rows=FALSE,
         show_rownames=TRUE,
         show_colnames = TRUE,
         color = viridis_pal()(10),
         cluster_cols=FALSE, 
         treeheight_row=0)
dev.off()

