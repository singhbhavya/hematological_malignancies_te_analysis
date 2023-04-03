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
library(fgsea)

################################### LOAD DATA ##################################

load("r_outputs/02-BL_filt_counts.Rdata")
load("r_outputs/05i-BL_pca_ccp_clusters_metadata.Rdata")
load("r_outputs/01-refs.Rdata")
load("r_outputs/01-metadata.Rdata")

################################################################################
################################################################################
#################################### CLUST 2 ###################################
################################################################################
################################################################################

#################################### DESEQ2 ####################################

# p value threshold 
lfc.cutoff <- 1.5
pval=0.001 

# Sanity check
stopifnot(all(colnames(BL.filt.comb) == rownames(BL_metadata)))

# Add metadata
BL_metadata$clust.retro.k2 <- clust.df$clust.retro.k2

# Run Deseq2
BL.k2.dds <- DESeq2::DESeqDataSetFromMatrix(BL.filt.comb, 
                                            BL_metadata, 
                                            ~ clust.retro.k2 + 0)


BL.k2.dds <- DESeq2::DESeq(BL.k2.dds, parallel=T)
BL.k2.tform <- DESeq2::varianceStabilizingTransformation(BL.k2.dds, blind=FALSE)

################################# DESEQ2 HERVs ################################# 

hervs.to.keep <- intersect(rownames(BL.filt.herv), 
                           retro.annot$locus[retro.annot$chrom != "chrY"])

BL.filt.herv.y <- BL.filt.herv[hervs.to.keep,] 
stopifnot(all(colnames(BL.filt.herv.y) == rownames(BL_metadata)))

BL.k2.herv.dds <- DESeq2::DESeqDataSetFromMatrix(BL.filt.herv.y, 
                                                  BL_metadata, 
                                                  ~ clust.retro.k2 + 0)


BL.k2.herv.dds <- DESeq2::DESeq(BL.k2.herv.dds, parallel=T)
BL.k2.herv.tform <- DESeq2::varianceStabilizingTransformation(BL.k2.herv.dds, 
                                                               blind=FALSE)


############################## TOP GENES & HERVS ###############################

# > resultsNames(BL.k2.dds)
# > "clust.retro.k2C1" "clust.retro.k2C2"

BL.res.k2 <- list(
  "C1" = DESeq2::results(BL.k2.dds, contrast=c(+1, -1), alpha=pval),
  "C2" = DESeq2::results(BL.k2.dds, contrast=c(-1, +1), alpha=pval))

BL.res.k2 <- lapply(BL.res.k2, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig.k2 <- lapply(BL.res.k2, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig.k2)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig.k2[[n]])
}

############################ UPSET GENES & HERVS ###############################

upvars.BL.k2 <- lapply(sig.k2[1:2], function(r) rownames(subset(r, log2FoldChange>0)))
downvars.BL.k2 <- lapply(sig.k2[1:2], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05k-BL_k2_upset_upvars.pdf", height=5, width=7)
upset(fromList(upvars.BL.k2), sets=c("C1", "C2"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

up.binmat.BL.k2 <- fromList(upvars.BL.k2)
rn <- do.call(c, upvars.BL.k2)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.BL.k2) <- rn
rm(rn)

dn.binmat.BL.k2 <- fromList(downvars.BL.k2)
rn <- do.call(c, downvars.BL.k2)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.BL.k2) <- rn
rm(rn)

################################### TOP HERVs ################################## 

# > resultsNames(BL.k2.dds)
# > "clust.retro.k2C1" "clust.retro.k2C2"

BL.res.k2.herv <- list(
  "C1" = DESeq2::results(BL.k2.herv.dds, contrast=c(+1, -1), alpha=pval),
  "C2" = DESeq2::results(BL.k2.herv.dds, contrast=c(-1, +1), alpha=pval))

BL.res.k2.herv <- lapply(BL.res.k2.herv, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig.k2.herv <- lapply(BL.res.k2.herv, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig.k2.herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig.k2.herv[[n]])
}

################################## UPSET HERVs ################################# 

upvars.BL.k2.hervs <- lapply(sig.k2.herv[1:2], function(r) rownames(subset(r, log2FoldChange>0)))
downvars.BL.k2.hervs <- lapply(sig.k2.herv[1:2], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05k-BL_k2_upset_hervs_upvars.pdf", height=5, width=7)
upset(fromList(upvars.BL.k2.hervs), sets=c("C1", "C2"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.

up.binmat.BL.k2.hervs <- fromList(upvars.BL.k2.hervs)
rn <- do.call(c, upvars.BL.k2.hervs)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.BL.k2.hervs) <- rn
rm(rn)

dn.binmat.BL.k2.hervs <- fromList(downvars.BL.k2.hervs)
rn <- do.call(c, downvars.BL.k2.hervs)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.BL.k2.hervs) <- rn
rm(rn)

################################# HEATMAPS #####################################

################################### SETUP ######################################

# Colors for cells
cols <- rgb_gsea(palette = c("default"), n = 14, alpha = 0.7, reverse = FALSE)

# Annotation column
df <- as.data.frame(BL_metadata[c("clust.retro.k2", "clinical_variant",
                                  "ebv_status", "subgroup", "gender",
                                  "MYC_SV_partner")])

# Create colors for each group
annoCol.k2 <-  c(wes_palette("Chevalier1")[1], 
              wes_palette("Chevalier1")[2])
names(annoCol.k2) <- c("C1", "C2")
annoCol.gender <- c(wes_palette("GrandBudapest2")[1:2])
names(annoCol.gender) <- c("Male", "Female")
annoCol.ebv <- c(wes_palette("Zissou1")[1],
                 wes_palette("Zissou1")[4])
names(annoCol.ebv) <- c("EBV-negative", "EBV-positive")
annoCol.subgroup <- wes_palette("Darjeeling2")
names(annoCol.subgroup) <- c("DGG-BL", "IC-BL", "DLBCL-B", "DLBCL-C", "Q53-BL")
annoCol.clinvar <- c(wes_palette("Darjeeling1")[1:2])
names(annoCol.clinvar) <- c("Endemic BL", "Sporadic BL")
annoCol.te.type <- c(wes_palette("IsleofDogs1")[1:3])
names(annoCol.te.type) <- c("INTERGENIC","EXONIC","INTRONIC"  )
annoCol.myc.sv <- c(wes_palette("IsleofDogs2")[1:4])
names(annoCol.myc.sv) <- c("BCL6", "IGH", "IGK", "IGL")

annoCol <- list(clust.retro.k2 = annoCol.k2,
                gender = annoCol.gender,
                ebv_status = annoCol.ebv,
                clinical_variant = annoCol.clinvar,
                subgroup = annoCol.subgroup,
                TE_type = annoCol.te.type,
                MYC_SV_partner = annoCol.myc.sv)

######################### UPREGULATED IN ALL GROUPS ############################

top.genes.hervs <- rownames(up.binmat.BL.k2)
top.hervs <- rownames(up.binmat.BL.k2.hervs)

pdf("plots/05k-BL_k2_no_y_top_hervs_genes_upregulated_all.pdf", height=12, width=12)
pheatmap(assay(BL.k2.tform)[top.genes.hervs,], 
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

pdf("plots/05k-BL_k2_no_y_top_hervs_upregulated_all.pdf", height=12, width=12)
pheatmap(assay(BL.k2.tform)[top.hervs,], 
         main="Upregulated HERvs, all clusters",
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
  mat <- assay(BL.k2.tform)[unique(topgenes), ]
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
           annotation_row = annotation_row,
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


for(clust in c("C1", "C2")) {
  tg <- rownames(sig.k2[[clust]][1:75,])
  p <- makeheatmap(tg, main=paste0('Genes/HERVs DE in cluster ', clust))
  pdf(paste0("plots/05k-BL_top_de_genes_hervs_", clust, ".pdf"), height=7, width=7)
  print(p)
  dev.off()
}

for(clust in c("C1", "C2")) {
  tg <- rownames(sig.k2.herv[[clust]][1:75,])
  p <- makeheatmap(tg, main=paste0('HERVs DE in cluster ', clust))
  pdf(paste0("plots/05k-BL_top_de_hervs_", clust, ".pdf"), height=7, width=7)
  print(p)
  dev.off()
}


############################## VOLCANO C1 v C2 #################################

pdf("plots/05k-BL_k2_volcano_C1_v_C2.pdf", height=8, width=8)
EnhancedVolcano(sig.k2$C1,
                lab = gene_table[rownames(sig.k2$C1),]$display,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'C2 v C1')
dev.off()

########################## PLOT INDIVIDUAL HERVs ###############################

plot.counts <- function(df, gene) {
  
  title <- gene_table[gene,]$display
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="clust.retro.k2", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=clust.retro.k2, y=count, fill=clust.retro.k2))  +
    geom_boxplot() +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("Cluster") +
    ylab("Counts") +
    scale_fill_manual(values = c("C1" = wes_palette("Chevalier1")[1], 
                                 "C2" = wes_palette("Chevalier1")[2])) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)
}


plot.counts(BL.k2.dds, "HARLEQUIN_7q32.1b")
plot.counts(BL.k2.dds, "ERVLB4_3p24.1b")


################################### PATHWAYS ###################################

pathways.hallmark <- gmtPathways("gsea/h.all.v2023.1.Hs.symbols.gmt")
pathways.immune <- gmtPathways("gsea/c7.immunesigdb.v2023.1.Hs.symbols.gmt")
pathways.bp <- gmtPathways("gsea/c5.go.bp.v2023.1.Hs.symbols.gmt")

################################ HALLMARK FGSEA ################################

k2.fgsea.res <- DESeq2::results(BL.k2.dds, contrast=c("clust.retro.k2", "C1", "C2"), 
                                alpha=pval)

k2.fgsea.res$SYMBOL <- gene_table[rownames(k2.fgsea.res),]$display

k2.fgsea.res <- as.data.frame(k2.fgsea.res) %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))

k2.ranks <- deframe(k2.fgsea.res)

k2.fgsea <- fgsea(pathways=pathways.hallmark, stats=k2.ranks, eps=0)

k2.fgseaResTidy <- k2.fgsea %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
k2.fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(k2.fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

################################# IMMUNE FGSEA #################################

k2.immune.fgsea <- fgsea(pathways=pathways.immune, stats=k2.ranks, eps=0)

k2.immune.fgseaResTidy <- k2.immune.fgsea %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
k2.immune.fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

k2.immune.pathways.subset <- k2.immune.fgseaResTidy[k2.immune.fgseaResTidy$pathway %like% "GC", ]

ggplot(k2.immune.fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

ggplot(k2.immune.pathways.subset, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

################################### BP FGSEA ###################################

k2.bp.fgsea <- fgsea(pathways=pathways.bp, stats=k2.ranks, eps=0, 
                     nPermSimple = 10000)

k2.bp.fgseaResTidy <- k2.bp.fgsea %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
k2.bp.fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(k2.bp.fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.001)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

################################################################################
################################################################################
################################### CLIN VAR ###################################
################################################################################
################################################################################

################################## DESEQ2 ALL ################################## 

stopifnot(all(colnames(BL.filt.comb) == rownames(BL_metadata)))

BL.ebvclin.dds <- DESeq2::DESeqDataSetFromMatrix(BL.filt.comb, 
                                            BL_metadata, 
                                            ~ 0 + subtype)


BL.ebvclin.dds <- DESeq2::DESeq(BL.ebvclin.dds, parallel=T)
BL.ebvclin.tform <- DESeq2::varianceStabilizingTransformation(BL.ebvclin.dds, blind=FALSE)

################################# DESEQ2 HERVs ################################# 

BL.ebvclin.herv.dds <- DESeq2::DESeqDataSetFromMatrix(BL.filt.herv.y, 
                                                  BL_metadata, 
                                                  ~ subtype + 0)


BL.ebvclin.herv.dds <- DESeq2::DESeq(BL.ebvclin.herv.dds, parallel=T)
BL.ebvclin.herv.tform <- DESeq2::varianceStabilizingTransformation(BL.ebvclin.herv.dds, 
                                                               blind=FALSE)

############################## TOP GENES & HERVS ###############################

## > resultsNames(BL.hc.dds)
## [1] "subtypeEndemic.BL.EBV.negative"  "subtypeEndemic.BL.EBV.positive"  "subtypeSporadic.BL.EBV.negative"
## [4] "subtypeSporadic.BL.EBV.positive"

BL.subtype.res <- list(
  "Endemic EBV Negative" = DESeq2::results(BL.ebvclin.dds, contrast=c(+1, -1/3, -1/3, -1/3), alpha=pval),
  "Endemic EBV Positive" = DESeq2::results(BL.ebvclin.dds, contrast=c(-1/3, +1, -1/3, -1/3), alpha=pval),
  "Sporadic EBV Negative" = DESeq2::results(BL.ebvclin.dds, contrast=c(-1/3, -1/3, +1, -1/3), alpha=pval),
  "Sporadic EBV Positive" = DESeq2::results(BL.ebvclin.dds, contrast=c(-1/3, -1/3, -1/3, +1), alpha=pval),
  "EBV Positive" = DESeq2::results(BL.ebvclin.dds, contrast=c(-1/2, +1/2, -1/2, +1/2), alpha=pval),
  "EBV Negative" = DESeq2::results(BL.ebvclin.dds, contrast=c(+1/2, -1/2, +1/2, -1/2), alpha=pval),
  "Endemic" = DESeq2::results(BL.ebvclin.dds, contrast=c(+1/2, +1/2, -1/2, -1/2), alpha=pval),
  "Sporadic" = DESeq2::results(BL.ebvclin.dds, contrast=c(-1/2, -1/2, +1/2, +1/2), alpha=pval)
)

BL.subtype.res <- lapply(BL.subtype.res, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig.subtype <- lapply(BL.subtype.res, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig.subtype)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig.subtype[[n]])
}

################################## TOP HERVS ###################################

BL.subtype.herv.res <- list(
  "Endemic EBV Negative" = DESeq2::results(BL.ebvclin.herv.dds, contrast=c(+1, -1/3, -1/3, -1/3), alpha=pval),
  "Endemic EBV Positive" = DESeq2::results(BL.ebvclin.herv.dds, contrast=c(-1/3, +1, -1/3, -1/3), alpha=pval),
  "Sporadic EBV Negative" = DESeq2::results(BL.ebvclin.herv.dds, contrast=c(-1/3, -1/3, +1, -1/3), alpha=pval),
  "Sporadic EBV Positive" = DESeq2::results(BL.ebvclin.herv.dds, contrast=c(-1/3, -1/3, -1/3, +1), alpha=pval),
  "EBV Positive" = DESeq2::results(BL.ebvclin.herv.dds, contrast=c(-1/2, +1/2, -1/2, +1/2), alpha=pval),
  "EBV Negative" = DESeq2::results(BL.ebvclin.herv.dds, contrast=c(+1/2, -1/2, +1/2, -1/2), alpha=pval),
  "Endemic" = DESeq2::results(BL.ebvclin.herv.dds, contrast=c(+1/2, +1/2, -1/2, -1/2), alpha=pval),
  "Sporadic" = DESeq2::results(BL.ebvclin.herv.dds, contrast=c(-1/2, -1/2, +1/2, +1/2), alpha=pval)
)

BL.subtype.herv.res <- lapply(BL.subtype.herv.res, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig.herv.subtype <- lapply(BL.subtype.herv.res, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig.herv.subtype)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig.herv.subtype[[n]])
}

############################ UPSET GENES & HERVS ###############################

upvars.BL.subtype <- lapply(sig.subtype[1:8], function(r) rownames(subset(r, log2FoldChange>0)))
downvars.BL.subtype <- lapply(sig.subtype[1:8], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05k-BL_subtype_upvars.pdf", height=5, width=7)
upset(fromList(upvars.BL.subtype), sets=c("Endemic EBV Negative",
                                      "Endemic EBV Positive",
                                      "Sporadic EBV Negative", 
                                      "Sporadic EBV Positive"),
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

upset(fromList(upvars.BL.subtype), sets=c("EBV Positive",
                                          "EBV Negative",
                                          "Endemic", 
                                          "Sporadic"),
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))

pdf("plots/05k-BL_subtype_downvars.pdf", height=5, width=7)
upset(fromList(downvars.BL.subtype), sets=c("Endemic EBV Negative",
                                        "Endemic EBV Positive",
                                        "Sporadic EBV Negative", 
                                        "Sporadic EBV Positive"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

upset(fromList(upvars.BL.subtype), sets=c("EBV Positive",
                                          "EBV Negative"),
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))

upset(fromList(upvars.BL.subtype), sets=c("Endemic",
                                          "Sporadic"),
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))

up.binmat.BL.subtype <- fromList(upvars.BL.subtype)
rn <- do.call(c, upvars.BL.subtype)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.BL.subtype) <- rn
rm(rn)

dn.binmat.BL.subtype <- fromList(downvars.BL.subtype)
rn <- do.call(c, downvars.BL.subtype)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.BL.subtype) <- rn
rm(rn)

############################### UPSET HERVs ONLY ############################### 

upvars.BL.herv.subtype <- lapply(sig.herv.subtype[1:8], function(r) rownames(subset(r, log2FoldChange>0)))
downvars.BL.herv.subtype <- lapply(sig.herv.subtype[1:8], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05k-BL_subtype_herv_upvars.pdf", height=5, width=7)
upset(fromList(upvars.BL.herv.subtype), sets=c("Endemic EBV Negative",
                                          "Endemic EBV Positive",
                                          "Sporadic EBV Negative", 
                                          "Sporadic EBV Positive"),
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

upset(fromList(upvars.BL.herv.subtype), sets=c("EBV Positive",
                                          "EBV Negative",
                                          "Endemic", 
                                          "Sporadic"),
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))


pdf("plots/05k-BL_subtype_herv_downvars.pdf", height=5, width=7)
upset(fromList(downvars.BL.herv.subtype), sets=c("Endemic EBV Negative",
                                            "Endemic EBV Positive",
                                            "Sporadic EBV Negative", 
                                            "Sporadic EBV Positive"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

upset(fromList(upvars.BL.herv.subtype), sets=c("EBV Positive",
                                          "EBV Negative"),
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))

upset(fromList(upvars.BL.herv.subtype), sets=c("Endemic",
                                          "Sporadic"),
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))

up.binmat.BL.subtype.herv <- fromList(upvars.BL.herv.subtype)
rn <- do.call(c, upvars.BL.herv.subtype)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.BL.subtype.herv) <- rn
rm(rn)

dn.binmat.BL.subtype.herv <- fromList(downvars.BL.herv.subtype)
rn <- do.call(c, downvars.BL.herv.subtype)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.BL.subtype.herv) <- rn
rm(rn)


################################# HEATMAPS #####################################

################################### SETUP ######################################

# Colors for cells
cols <- rgb_gsea(palette = c("default"), n = 14, alpha = 0.7, reverse = FALSE)

# Annotation column
df <- as.data.frame(colData(BL.ebvclin.dds)[,c("clinical_variant","ebv_status", "subtype",
                                           "subgroup", "MYC_SV_partner", "gender")])


######################### UPREGULATED IN ALL GROUPS ############################

top.genes.hervs <- rownames(up.binmat.BL.subtype)
top.hervs <- rownames(up.binmat.BL.subtype.herv)

annoRow <- as.data.frame(retro.annot.v2[,c("TE_type", "Locus")])
annoRow <- annoRow[top.hervs,]
annoRow <- subset(annoRow, select = -c(2))

pdf("plots/05k-BL_top_hervs_upregulated_all.pdf", height=10, width=10)
pheatmap(assay(BL.ebvclin.herv.tform)[top.hervs,], 
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
         annotation_col = df,
         annotation_colors = annoCol,
         annotation_row=annoRow)
dev.off()

pdf("plots/05f-gcb_agirre_top_hervs_genes_upregulated_all.pdf", height=10, width=10)
pheatmap(assay(BL.ebvclin.tform)[top.genes.hervs,],
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
         annotation_col = df,
         annotation_colors = annoCol)
dev.off()


for(clust in c("Endemic EBV Positive",
               "Sporadic EBV Negative")) {
  tg <- rownames(sig.subtype[[clust]][1:75,])
  p <- makeheatmap(tg, main=paste0('Genes/HERVs DE in cluster ', clust))
  pdf(paste0("plots/05k-BL_top_de_genes_hervs_", clust, ".pdf"), height=7, width=7)
  print(p)
  dev.off()
}

