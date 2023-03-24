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
retro.annot.v2 <- read.csv("/efs/projects/hematological_malignancies_te_analysis/refs/TE_annotation.v2.0.tsv",
                           sep = "\t")
rownames(retro.annot.v2) <- retro.annot.v2$Locus

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

################################################################################
################################################################################
################################################################################
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

########################## GCB BULK DESEQ (GENES ONLY) #########################

### DESeq2 (HERVs + Genes)

bulk.gonly.countdat <- GCB_Agirre.filt.tx
cat(sprintf('%d variables\n', nrow(bulk.gonly.countdat)))

stopifnot(all(colnames(bulk.gonly.countdat) == rownames(agirre_metadata)))

Agirre.gonly.dds <- DESeq2::DESeqDataSetFromMatrix(countData = bulk.gonly.countdat,
                                               colData = agirre_metadata,
                                               design = ~ Cell + 0)

Agirre.gonly.dds <- DESeq2::DESeq(Agirre.gonly.dds, parallel=T)
Agirre.gonly.tform <- DESeq2::varianceStabilizingTransformation(Agirre.gonly.dds, blind=FALSE)

## PCA
Agirre.gonly.pca.obj <-
  pca_standard(tform = Agirre.gonly.tform, 
               metadata = agirre_metadata, 
               var = 0.1)

# 4 PCs for Elbow method
# 5 PCs for Horn method
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

pdf("plots/05f-gcb_agirre_pca_genes.pdf", height=5, width=6)
## Biplot with projects (HERVs only)
biplot(Agirre.gonly.pca.obj, 
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
                  "Dark Zone Germinal Center B" = pal_jco("default", alpha = 0.8)(7)[4],
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

################################################################################
################################################################################
################################################################################
####################### ANALYSIS 1: EACH CELL TYPE VS ALL ######################
################################################################################
################################################################################
################################################################################

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

########################### TOP GENES & HERVS NO PB ############################

# resultsNames(Agirre.g.dds)
# [1] "CellBone.Marrow.plasma.cell" "CellDark.Zone.Germinal.Center.B"  "CellLight.Zone.Germinal.Center.B"
# [4] "CellMemory.B""CellNaive.B" "CellPlasmablasts"                              

agirre_res_no_pb <- list(
  "DZ" = DESeq2::results(Agirre.g.dds, contrast=c(0, +1, -1/3, -1/3, -1/3, 0), alpha=pval),
  "LZ" = DESeq2::results(Agirre.g.dds, contrast=c(0, -1/3, +1, -1/3, -1/3, 0), alpha=pval),
  "MB" = DESeq2::results(Agirre.g.dds, contrast=c(0, -1/3, -1/3, +1, -1/3, 0), alpha=pval),
  "NB" = DESeq2::results(Agirre.g.dds, contrast=c(0, -1/3, -1/3, -1/3, +1, 0), alpha=pval),
  "DZvLZ" = DESeq2::results(Agirre.g.dds, contrast=c("Cell", "Dark Zone Germinal Center B", 
                                                     "Light Zone Germinal Center B"), alpha=pval),
  "MBvsNB" = DESeq2::results(Agirre.g.dds, contrast=c("Cell", "Memory B", "Naive B"), alpha=pval)
)

agirre_res_no_pb <- lapply(agirre_res_no_pb, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig_no_pb <- lapply(agirre_res_no_pb, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig_no_pb)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig_no_pb[[n]])
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

############################ TOP HERVs ONLY NO PB ##############################

agirre_res_herv_no_pb <- list(
  "DZ" = DESeq2::results(Agirre.dds, contrast=c(0, +1, -1/3, -1/3, -1/3, 0), alpha=pval),
  "LZ" = DESeq2::results(Agirre.dds, contrast=c(0, -1/3, +1, -1/3, -1/3, 0), alpha=pval),
  "MB" = DESeq2::results(Agirre.dds, contrast=c(0, -1/3, -1/3, +1, -1/3, 0), alpha=pval),
  "NB" = DESeq2::results(Agirre.dds, contrast=c(0, -1/3, -1/3, -1/3, +1, 0), alpha=pval),
  "DZvLZ" = DESeq2::results(Agirre.dds, contrast=c("Cell", "Dark Zone Germinal Center B", 
                                                   "Light Zone Germinal Center B"), alpha=pval),
  "MBvsNB" = DESeq2::results(Agirre.dds, contrast=c("Cell", "Memory B", "Naive B"), alpha=pval)
)

agirre_res_herv_no_pb <- lapply(agirre_res_herv_no_pb, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig_herv_no_pb <- lapply(agirre_res_herv_no_pb, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

sink(file = "r_outputs/05f-agirre_herv_no_pb_deseq.txt")
for (n in names(sig_herv_no_pb)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig_herv_no_pb[[n]])
}
sink(file = NULL)

sink(file = "r_outputs/05f-agirre_herv_no_pb_deseq_top5.txt")
for (n in names(sig_herv_no_pb)) {
  cat("\n#--- Contrast", n, "---#\n")
  print(head(sig_herv_no_pb[[n]]))
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

######################### UPSET GENES & HERVS NO PB ############################

upvars_agirre_no_pb <- lapply(sig_no_pb[1:4], function(r) rownames(subset(r, log2FoldChange>0)))
downvars_agirre_no_pb <- lapply(sig_no_pb[1:4], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05f-gcb_agirre_no_pb_upvars.pdf", height=5, width=7)
upset(fromList(upvars_agirre_no_pb), sets=c("DZ", "LZ", "MB", "NB"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

pdf("plots/05f-gcb_agirre_no_pb_upset_dnwars.pdf", height=5, width=7)
upset(fromList(downvars_agirre_no_pb), sets=c("DZ", "LZ", "MB", "NB"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

up.binmat.agirre.no.pb <- fromList(upvars_agirre_no_pb)
rn <- do.call(c, upvars_agirre_no_pb)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.agirre.no.pb) <- rn
rm(rn)

dn.binmat.agirre.no.pb <- fromList(downvars_agirre_no_pb)
rn <- do.call(c, downvars_agirre_no_pb)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.agirre.no.pb) <- rn
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

########################## UPSET HERVs ONLY NO PB ##############################

upvars_agirre_hervs_no_pb <- lapply(sig_herv_no_pb[1:4], function(r) rownames(subset(r, log2FoldChange>0)))
downvars_agirre_hervs_no_pb <- lapply(sig_herv_no_pb[1:4], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05f-gcb_agirre_upset_upvars_hervs_no_pb.pdf", height=5, width=7)
upset(fromList(upvars_agirre_hervs_no_pb), sets=c("DZ", "LZ", "MB", "NB"),  
      keep.order = T, order.by='degree', decreasing=F)
dev.off()

pdf("plots/05f-gcb_agirre_upset_dnvars_hervs_no_pb.pdf", height=5, width=7)
upset(fromList(downvars_agirre_hervs_no_pb), sets=c("DZ", "LZ", "MB", "NB"),  
      keep.order = T, order.by='degree', decreasing=F)
dev.off()

up.binmat.hervs.agirre.no.pb <- fromList(upvars_agirre_hervs_no_pb)
rn <- do.call(c, upvars_agirre_hervs_no_pb)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.hervs.agirre.no.pb) <- rn
rm(rn)

dn.binmat.hervs.agirre.no.pb <- fromList(downvars_agirre_hervs_no_pb)
rn <- do.call(c, downvars_agirre_hervs_no_pb)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.hervs.agirre.no.pb) <- rn
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
annoCol2 <- pal_aaas("default", alpha=0.7)(3)
names(annoCol2) <- unique(retro.annot.v2$TE_type)
annoCol <- list(Cell = annoCol, TE_type = annoCol2)

######################### UPREGULATED IN ALL GROUPS ############################

top.genes.hervs <- rownames(up.binmat.agirre)
top.hervs <- rownames(up.binmat.hervs.agirre)

annoRow <- as.data.frame(retro.annot.v2[,c("TE_type", "Locus")])
annoRow <- annoRow[top.hervs,]
annoRow <- subset(annoRow, select = -c(2))

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
         annotation_row=annoRow,
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
  
  annoRow <- as.data.frame(retro.annot.v2[,c("TE_type", "Locus")])
  annoRow <- annoRow[topgenes,]
  annoRow <- subset(annoRow, select = -c(2))
  
  pheatmap(mat,
           color=cols,
           scale="row", breaks=seq(-3,3,length.out=14),
           clustering_distance_rows = rowdist,
           clustering_method="average",
           annotation_col = df,
           annotation_colors = annoCol,
           annotation_row = annoRow,
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

################################ FAMILY LEVEL UP ############################### 

upreg_agirre_hervs_df <- do.call(rbind, lapply(upvars_agirre_hervs, data.frame))
colnames(upreg_agirre_hervs_df) <- c("herv")
upreg_agirre_hervs_df$cell_type <- rownames(upreg_agirre_hervs_df)
upreg_agirre_hervs_df$cell_type <- gsub("\\..*","",upreg_agirre_hervs_df$cell_type)
upreg_agirre_hervs_df$family <- retro.annot$family[match(upreg_agirre_hervs_df$herv, 
                                                        retro.annot$locus)]

upreg_agirre_families <-
  upreg_agirre_hervs_df %>% dplyr::count(family, cell_type, sort = TRUE) 

upreg_agirre_families<-upreg_agirre_families[!(upreg_agirre_families$cell_type=="DZvLZ"),]

pdf("plots/05f-gcb_agirre_upreg_families_all.pdf", height=8, width=6)
ggplot(upreg_agirre_families, aes(fill=reorder(family, -n), y=cell_type, x=n)) + 
  geom_bar(position="fill", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7)),
                    breaks = unique(retro.annot$family),
                    labels = unique(retro.annot$family)) + 
  
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="TE Family")) +
  ylab("Upregulated HERVs by Family") +
  xlab("Proportion of HERVs") + 
  theme(legend.position = c("right")) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2))
dev.off()

############################### FAMILY LEVEL DOWN ##############################

down_agirre_hervs_df <- do.call(rbind, lapply(downvars_agirre_hervs, data.frame))
colnames(down_agirre_hervs_df) <- c("herv")
down_agirre_hervs_df$cell_type <- rownames(down_agirre_hervs_df)
down_agirre_hervs_df$cell_type <- gsub("\\..*","",down_agirre_hervs_df$cell_type)
down_agirre_hervs_df$family <- retro.annot$family[match(down_agirre_hervs_df$herv, 
                                                       retro.annot$locus)]

downreg_agirre_families <-
  down_agirre_hervs_df %>% dplyr::count(family, cell_type, sort = TRUE) 

downreg_agirre_families<-downreg_agirre_families[!(downreg_agirre_families$cell_type=="DZvLZ"),]

pdf("plots/05f-gcb_agirre_downreg_families_all.pdf", height=8, width=6)
ggplot(downreg_agirre_families, aes(fill=reorder(family, -n), y=cell_type, x=n)) + 
  geom_bar(position="fill", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7)),
                    breaks = unique(retro.annot$family),
                    labels = unique(retro.annot$family)) + 
  
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="TE Family")) +
  ylab("Downregulated HERVs by Family") +
  xlab("Proportion of HERVs") + 
  theme(legend.position = c("right")) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2))
dev.off()

############################# MERGED FAMILY UP/DOWN ############################


updown_family <- do.call(rbind, (list(Upregulated = upreg_agirre_families, 
                                      Dowregulated = downreg_agirre_families)))
updown_family$expression <- rownames(updown_family)
updown_family$expression <- gsub("\\..*","",updown_family$expression)
rownames(updown_family)<-NULL

pdf("plots/05f-gcb_agirre_updown_families_all.pdf", height=8, width=10)
ggplot(updown_family, aes(fill=reorder(family, -n), y=cell_type, x=n)) + 
  geom_bar(position="fill", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7)),
                    breaks = unique(retro.annot$family),
                    labels = unique(retro.annot$family)) + 
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="TE Family")) +
  ylab("Downregulated HERVs by Family") +
  xlab("Proportion of HERVs") + 
  theme(legend.position = c("right")) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2)) +
  facet_wrap(~ expression, ncol = 2)
dev.off()

#################################### TE TYPE ###################################

rownames(upreg_agirre_hervs_df) <- NULL
upreg_agirre_hervs_df$type <- retro.annot.v2$TE_type[
  match(upreg_agirre_hervs_df$herv, retro.annot.v2$Locus)]

upreg_agirre_te_type <-
  upreg_agirre_hervs_df %>% dplyr::count(type, cell_type, sort = TRUE) 

pdf("plots/05f-gcb_agirre_hervs_type.pdf", height=4, width=4)
ggplot(upreg_agirre_te_type, aes(fill=reorder(type, -n), y=cell_type, x=n)) + 
  geom_bar(position="stack", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_aaas("default", alpha=0.7)(3))) + 
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="HERV Type")) +
  ylab("Upregulated HERVs by Type") +
  xlab("Number of HERVs") + 
  theme(legend.position = c("right")) +
  theme(aspect.ratio = 1)
dev.off()

rownames(down_agirre_hervs_df) <- NULL
down_agirre_hervs_df$type <- retro.annot.v2$TE_type[
  match(down_agirre_hervs_df$herv, retro.annot.v2$Locus)]

downreg_agirre_te_type <-
  down_agirre_hervs_df %>% dplyr::count(type, cell_type, sort = TRUE) 

pdf("plots/05f-gcb_agirre_downreg_hervs_type.pdf", height=4, width=4)
ggplot(downreg_agirre_te_type, aes(fill=reorder(type, -n), y=cell_type, x=n)) + 
  geom_bar(position="stack", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_aaas("default", alpha=0.7)(3))) + 
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="HERV Type")) +
  ylab("Downregulated HERVs by Type") +
  xlab("Number of HERVs") + 
  theme(legend.position = c("right")) +
  theme(aspect.ratio = 1)
dev.off()


########################### FAMILY LEVEL UP NO PB ##############################

upreg_hervs_no_pb_df <- do.call(rbind, lapply(upvars_agirre_hervs_no_pb, data.frame))
colnames(upreg_hervs_no_pb_df) <- c("herv")
upreg_hervs_no_pb_df$cell_type <- rownames(upreg_hervs_no_pb_df)
upreg_hervs_no_pb_df$cell_type <- gsub("\\..*","",upreg_hervs_no_pb_df$cell_type)
upreg_hervs_no_pb_df$family <- retro.annot$family[match(upreg_hervs_no_pb_df$herv, 
                                                         retro.annot$locus)]

upreg_families_no_pb <-
  upreg_hervs_no_pb_df %>% dplyr::count(family, cell_type, sort = TRUE) 

upreg_families_no_pb<-upreg_families_no_pb[!(upreg_families_no_pb$cell_type=="DZvLZ"),]

pdf("plots/05f-gcb_agirre_upreg_families.pdf", height=8, width=6)
ggplot(upreg_families_no_pb, aes(fill=reorder(family, -n), y=cell_type, x=n)) + 
  geom_bar(position="fill", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7)),
                    breaks = unique(retro.annot$family),
                    labels = unique(retro.annot$family)) + 
  
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="TE Family")) +
  ylab("Upregulated HERVs by Family") +
  xlab("Proportion of HERVs") + 
  theme(legend.position = c("right")) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2))
dev.off()

############################ FAMILY LEVEL DOWN NO PB ###########################

down_hervs_no_pb_df <- do.call(rbind, lapply(downvars_agirre_hervs_no_pb, data.frame))
colnames(down_hervs_no_pb_df) <- c("herv")
down_hervs_no_pb_df$cell_type <- rownames(down_hervs_no_pb_df)
down_hervs_no_pb_df$cell_type <- gsub("\\..*","",down_hervs_no_pb_df$cell_type)
down_hervs_no_pb_df$family <- retro.annot$family[match(down_hervs_no_pb_df$herv, 
                                                        retro.annot$locus)]

downreg_families_no_pb <-
  down_hervs_no_pb_df %>% dplyr::count(family, cell_type, sort = TRUE) 

downreg_families_no_pb<-downreg_families_no_pb[!(downreg_families_no_pb$cell_type=="DZvLZ"),]

pdf("plots/05f-gcb_agirre_families.pdf", height=8, width=6)
ggplot(downreg_families_no_pb, aes(fill=reorder(family, -n), y=cell_type, x=n)) + 
  geom_bar(position="fill", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7)),
                    breaks = unique(retro.annot$family),
                    labels = unique(retro.annot$family)) + 
  
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="TE Family")) +
  ylab("Downregulated HERVs by Family") +
  xlab("Proportion of HERVs") + 
  theme(legend.position = c("right")) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2))
dev.off()


########################## MERGED FAMILY UP/DOWN NO PB #########################


updown_family <- do.call(rbind, (list(Upregulated = upreg_families_no_pb, 
                                      Dowregulated = downreg_families_no_pb)))
updown_family$expression <- rownames(updown_family)
updown_family$expression <- gsub("\\..*","",updown_family$expression)
rownames(updown_family)<-NULL

pdf("plots/05f-gcb_agirre_updown_families.pdf", height=8, width=10)
ggplot(updown_family, aes(fill=reorder(family, -n), y=cell_type, x=n)) + 
  geom_bar(position="fill", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7)),
                    breaks = unique(retro.annot$family),
                    labels = unique(retro.annot$family)) + 
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="TE Family")) +
  ylab("Downregulated HERVs by Family") +
  xlab("Proportion of HERVs") + 
  theme(legend.position = c("right")) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2)) +
  facet_wrap(~ expression, ncol = 2)
dev.off()


############################ FAMILY LEVEL WITH PB ##############################

upreg_hervs_pb <- do.call(rbind, lapply(upvars_agirre_hervs, data.frame))
colnames(upreg_hervs_pb) <- c("herv")
upreg_hervs_pb$cell_type <- rownames(upreg_hervs_pb)
upreg_hervs_pb$cell_type <- gsub("\\..*","",upreg_hervs_pb$cell_type)
upreg_hervs_pb$family <- retro.annot$family[match(upreg_hervs_no_pb_df$herv, 
                                                        retro.annot$locus)]

upreg_families_pb <-
  upreg_hervs_pb %>% dplyr::count(family, cell_type, sort = TRUE) 

upreg_families_pb<-upreg_families_pb[!(upreg_families_pb$cell_type=="DZvLZ"),]

down_hervs_pb <- do.call(rbind, lapply(downvars_agirre_hervs, data.frame))
colnames(down_hervs_pb) <- c("herv")
down_hervs_pb$cell_type <- rownames(down_hervs_pb)
down_hervs_pb$cell_type <- gsub("\\..*","",down_hervs_pb$cell_type)
down_hervs_pb$family <- retro.annot$family[match(down_hervs_pb$herv, 
                                                       retro.annot$locus)]

downreg_families_pb <-
  down_hervs_pb %>% dplyr::count(family, cell_type, sort = TRUE) 

downreg_families_pb<-downreg_families_pb[!(downreg_families_pb$cell_type=="DZvLZ"),]

updown_family <- do.call(rbind, (list(Upregulated = upreg_families_pb, 
                                      Dowregulated = downreg_families_pb)))
updown_family$expression <- rownames(updown_family)
updown_family$expression <- gsub("\\..*","",updown_family$expression)
rownames(updown_family)<-NULL

pdf("plots/05f-gcb_agirre_updown_families_with_pb.pdf", height=8, width=10)
ggplot(updown_family, aes(fill=reorder(family, -n), y=cell_type, x=n)) + 
  geom_bar(position="fill", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7)),
                    breaks = unique(retro.annot$family),
                    labels = unique(retro.annot$family)) + 
  
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="TE Family")) +
  ylab("Downregulated HERVs by Family") +
  xlab("Proportion of HERVs") + 
  theme(legend.position = c("right")) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2)) +
  facet_wrap(~ expression, ncol = 2)

dev.off()

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


################################################################################
################################################################################
################################################################################
######################## ANALYSIS 2: NAIVE AS BASELINE #########################
################################################################################
################################################################################
################################################################################

agirre_res_naive <- list(
  "BMPC" = DESeq2::results(Agirre.g.dds, contrast=c("Cell", "Bone Marrow plasma cell", 
                                                    "Naive B"), alpha=pval),
  "DZ" = DESeq2::results(Agirre.g.dds, contrast=c("Cell", "Dark Zone Germinal Center B", 
                                                  "Naive B"), alpha=pval),
  "LZ" = DESeq2::results(Agirre.g.dds, contrast=c("Cell", "Light Zone Germinal Center B", 
                                                  "Naive B"), alpha=pval),
  "MB" = DESeq2::results(Agirre.g.dds, contrast=c("Cell", "Memory B", 
                                                  "Naive B"), alpha=pval),
  "NB" = DESeq2::results(Agirre.g.dds, contrast=c(-1/5, -1/5, -1/5, -1/5, +1, -1/5), alpha=pval),
  "PB" = DESeq2::results(Agirre.g.dds, contrast=c("Cell", "Plasmablasts", 
                                                  "Naive B"))
)

agirre_res_naive <- lapply(agirre_res_naive, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r
})

sig_naive <- lapply(agirre_res_naive, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

for (n in names(sig_naive)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig_naive[[n]])
}

############################### TOP HERVs ONLY #################################

agirre_res_naive_herv <- list(
  "BMPC" = DESeq2::results(Agirre.dds, contrast=c("Cell", "Bone Marrow plasma cell", 
                                                    "Naive B"), alpha=pval),
  "DZ" = DESeq2::results(Agirre.dds, contrast=c("Cell", "Dark Zone Germinal Center B", 
                                                  "Naive B"), alpha=pval),
  "LZ" = DESeq2::results(Agirre.dds, contrast=c("Cell", "Light Zone Germinal Center B", 
                                                  "Naive B"), alpha=pval),
  "MB" = DESeq2::results(Agirre.dds, contrast=c("Cell", "Memory B", "Naive B"), alpha=pval),
  "NB" = DESeq2::results(Agirre.dds, contrast=c(-1/5, -1/5, -1/5, -1/5, +1, -1/5), alpha=pval),
  "PB" = DESeq2::results(Agirre.dds, contrast=c("Cell", "Plasmablasts", 
                                                  "Naive B"))
)

agirre_res_naive_herv <- lapply(agirre_res_naive_herv, function(r) {
  r$display <- gene_table[rownames(r),]$display
  r$class <- gene_table[rownames(r),]$gene_type
  r$TE_type <- retro.annot.v2[rownames(r),]$TE_type
  r
})


sig_naive_herv <- lapply(agirre_res_naive_herv, function(r) {
  s <- subset(r, padj < pval & abs(log2FoldChange) > lfc.cutoff)
  s[order(s$padj),]
})

sink(file = "r_outputs/05f-agirre_naive_herv_deseq.txt")
for (n in names(sig_naive_herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig_naive_herv[[n]])
}
sink(file = NULL)

sink(file = "r_outputs/05f-agirre_naive_herv_deseq_top5.txt")
for (n in names(sig_naive_herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  print(head(sig_naive_herv[[n]]))
}
sink(file = NULL)


############################ UPSET GENES & HERVS ###############################

upvars_agirre_naive <- lapply(sig_naive[1:6], function(r) rownames(subset(r, log2FoldChange>0)))
downvars_agirre_naive <- lapply(sig_naive[1:6], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05f-gcb_agirre_naive_upvars.pdf", height=5, width=7)
upset(fromList(upvars_agirre_naive), sets=c("BMPC", "DZ", "LZ", "MB", "NB", "PB"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

pdf("plots/05f-gcb_agirre_naive_upset_dnwars.pdf", height=5, width=7)
upset(fromList(downvars_agirre_naive), sets=c("BMPC", "DZ", "LZ", "MB", "NB", "PB"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

up.binmat.agirre.naive <- fromList(upvars_agirre_naive)
rn <- do.call(c, upvars_agirre_naive)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.agirre.naive) <- rn
rm(rn)

dn.binmat.agirre.naive <- fromList(downvars_agirre_naive)
rn <- do.call(c, downvars_agirre_naive)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.agirre.naive) <- rn
rm(rn)

############################### UPSET HERVS ONLY ###############################

upvars_agirre_naive_herv <- lapply(sig_naive_herv[1:6], function(r) rownames(subset(r, log2FoldChange>0)))
downvars_agirre_naive_herv <- lapply(sig_naive_herv[1:6], function(r) rownames(subset(r, log2FoldChange<0)))

pdf("plots/05f-gcb_agirre_naive_upvars_hervs.pdf", height=5, width=7)
upset(fromList(upvars_agirre_naive_herv), sets=c("DZ", "LZ", "MB", "NB", "PB", "BMPC"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

pdf("plots/05f-gcb_agirre_naive_upset_dnwars_hervs.pdf", height=5, width=7)
upset(fromList(downvars_agirre_naive_herv), sets=c("DZ", "LZ", "MB", "NB", "PB", "BMPC"),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1))
dev.off()

up.binmat.agirre.naive.hervs <- fromList(upvars_agirre_naive_herv)
rn <- do.call(c, upvars_agirre_naive_herv)
rn <- rn[!duplicated(rn)]
rownames(up.binmat.agirre.naive.hervs) <- rn
rm(rn)

dn.binmat.agirre.naive.hervs <- fromList(downvars_agirre_naive_herv)
rn <- do.call(c, downvars_agirre_naive_herv)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat.agirre.naive.hervs) <- rn
rm(rn)


######################### UPREGULATED IN ALL GROUPS ############################

top.genes.hervs <- rownames(up.binmat.agirre.naive)
top.hervs <- rownames(up.binmat.agirre.naive.hervs)

annoRow <- as.data.frame(retro.annot.v2[,c("TE_type", "Locus")])
annoRow <- annoRow[top.hervs,]
annoRow <- subset(annoRow, select = -c(2))

pdf("plots/05f-gcb_agirre_naive_top_hervs_upregulated_all.pdf", height=10, width=10)
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
         annotation_row=annoRow,
         annotation_colors = annoCol)
dev.off()

pdf("plots/05f-gcb_agirre_naive_top_hervs_genes_upregulated_all.pdf", height=10, width=10)
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


for(clust in c("BMPC", "DZ", "LZ", "NB", "PB")) {
  tg <- rownames(sig_naive_herv[[clust]][1:75,])
  p <- makeheatmap(tg, main=paste0('DE in cluster ', clust))
  pdf(paste0("plots/05f-gcb_agirre_naive_top_de_hervs_", clust, ".pdf"), height=7, width=7)
  print(p)
  dev.off()
}

for(clust in c("MB")) {
  tg <- rownames(sig_naive_herv[[clust]][1:39,])
  p <- makeheatmap(tg, main=paste0('DE in cluster ', clust))
  pdf(paste0("plots/05f-gcb_agirre_naive_top_de_hervs_", clust, ".pdf"), height=5, width=7)
  print(p)
  dev.off()
}

################################# VOLCANO PLOTS ################################

pdf("plots/05f-gcb_agirre_naive_volcano_DZ_v_NB.pdf", height=8, width=8)
EnhancedVolcano(sig_naive_herv$DZ,
                lab = rownames(sig_naive_herv$DZ),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'NB vs DZ')
dev.off()

pdf("plots/05f-gcb_agirre_naive_volcano_LZ_v_NB.pdf", height=8, width=8)
EnhancedVolcano(sig_naive_herv$LZ,
                lab = rownames(sig_naive_herv$LZ),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'MB vs LZ')
dev.off()

pdf("plots/05f-gcb_agirre_naive_volcano_NB_v_all.pdf", height=8, width=8)
EnhancedVolcano(sig_naive_herv$NB,
                lab = rownames(sig_naive_herv$NB),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'All vs NB')
dev.off()

pdf("plots/05f-gcb_agirre_naive_volcano_MB_v_NB.pdf", height=8, width=8)
EnhancedVolcano(sig_naive_herv$MB,
                lab = rownames(sig_naive_herv$MB),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'NB vs MB')
dev.off()

pdf("plots/05f-gcb_agirre_naive_volcano_PB_v_NB.pdf", height=8, width=8)
EnhancedVolcano(sig_naive_herv$PB,
                lab = rownames(sig_naive_herv$PB),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'NB vs PB')
dev.off()

pdf("plots/05f-gcb_agirre_naive_volcano_BMPC_v_NBl.pdf", height=8, width=8)
EnhancedVolcano(sig_naive_herv$BMPC,
                lab = rownames(sig_naive_herv$BMPC),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'NB vs BMPC')
dev.off()

############################### FAMILY LEVEL UP ################################

upreg_agirre_hervs_nb_df <- do.call(rbind, lapply(upvars_agirre_naive_herv, data.frame))
colnames(upreg_agirre_hervs_nb_df) <- c("herv")
upreg_agirre_hervs_nb_df$cell_type <- rownames(upreg_agirre_hervs_nb_df)
upreg_agirre_hervs_nb_df$cell_type <- gsub("\\..*","",upreg_agirre_hervs_nb_df$cell_type)
upreg_agirre_hervs_nb_df$family <- retro.annot$family[match(upreg_agirre_hervs_nb_df$herv, 
                                                     retro.annot$locus)]

upreg_agirre_families_naive <-
  upreg_agirre_hervs_nb_df %>% dplyr::count(family, cell_type, sort = TRUE) 

pdf("plots/05f-gcb_agirre_naive_upreg_families.pdf", height=8, width=6)
ggplot(upreg_agirre_families_naive, aes(fill=reorder(family, -n), y=cell_type, x=n)) + 
  geom_bar(position="fill", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7)),
                    breaks = unique(retro.annot$family),
                    labels = unique(retro.annot$family)) + 
  
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="TE Family")) +
  ylab("Upregulated HERVs by Family") +
  xlab("Proportion of HERVs") + 
  theme(legend.position = c("right")) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2))
dev.off()

############################### FAMILY LEVEL DOWN ##############################

down_agirre_hervs_nb_df <- do.call(rbind, lapply(downvars_agirre_naive_herv, data.frame))
colnames(down_agirre_hervs_nb_df) <- c("herv")
down_agirre_hervs_nb_df$cell_type <- rownames(down_agirre_hervs_nb_df)
down_agirre_hervs_nb_df$cell_type <- gsub("\\..*","",down_agirre_hervs_nb_df$cell_type)
down_agirre_hervs_nb_df$family <- retro.annot$family[match(down_agirre_hervs_nb_df$herv, 
                                                       retro.annot$locus)]

downreg_agirre_families_naive <-
  down_agirre_hervs_nb_df %>% dplyr::count(family, cell_type, sort = TRUE) 

pdf("plots/05f-gcb_agirre_naive_downreg_families.pdf", height=8, width=6)
ggplot(downreg_agirre_families_naive, aes(fill=reorder(family, -n), y=cell_type, x=n)) + 
  geom_bar(position="fill", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7)),
                    breaks = unique(retro.annot$family),
                    labels = unique(retro.annot$family)) + 
  
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="TE Family")) +
  ylab("Downregulated HERVs by Family") +
  xlab("Proportion of HERVs") + 
  theme(legend.position = c("right")) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2))
dev.off()

############################# MERGED FAMILY UP/DOWN ############################


updown_family <- do.call(rbind, (list(Upregulated = upreg_agirre_families_naive, 
                                      Dowregulated = downreg_agirre_families_naive)))
updown_family$expression <- rownames(updown_family)
updown_family$expression <- gsub("\\..*","",updown_family$expression)
rownames(updown_family)<-NULL

pdf("plots/05f-gcb_agirre_naive_updown_families.pdf", height=8, width=10)
ggplot(updown_family, aes(fill=reorder(family, -n), y=cell_type, x=n)) + 
  geom_bar(position="fill", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_futurama("planetexpress")(12), 
                               pal_npg("nrc", alpha = 0.7)(10),
                               pal_jco("default", alpha=0.7)(10),
                               pal_nejm("default", alpha=0.7)(8),
                               pal_tron("legacy", alpha=0.7)(7),
                               pal_lancet("lanonc", alpha=0.7)(9),
                               pal_startrek("uniform", alpha=0.7)(7)),
                    breaks = unique(retro.annot$family),
                    labels = unique(retro.annot$family)) + 
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="TE Family")) +
  ylab("Downregulated HERVs by Family") +
  xlab("Proportion of HERVs") + 
  theme(legend.position = c("right")) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2)) +
  facet_wrap(~ expression, ncol = 2)
dev.off()

#################################### TE TYPE ###################################

rownames(upreg_agirre_hervs_nb_df) <- NULL
upreg_agirre_hervs_nb_df$type <- retro.annot.v2$TE_type[
  match(upreg_agirre_hervs_nb_df$herv, retro.annot.v2$Locus)]

upreg_agirre_naive_te_type <-
  upreg_agirre_hervs_nb_df %>% dplyr::count(type, cell_type, sort = TRUE) 

pdf("plots/05f-gcb_agirre_naive_upreg_hervs_type.pdf", height=4, width=4)
ggplot(upreg_agirre_naive_te_type, aes(fill=reorder(type, -n), y=cell_type, x=n)) + 
  geom_bar(position="stack", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_aaas("default", alpha=0.7)(3))) + 
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="HERV Type")) +
  ylab("Upregulated HERVs by Type") +
  xlab("Number of HERVs") + 
  theme(legend.position = c("right")) +
  theme(aspect.ratio = 1)
dev.off()

rownames(down_agirre_hervs_nb_df) <- NULL
down_agirre_hervs_nb_df$type <- retro.annot.v2$TE_type[
  match(down_agirre_hervs_nb_df$herv, retro.annot.v2$Locus)]

downreg_agirre_naive_te_type <-
  down_agirre_hervs_nb_df %>% dplyr::count(type, cell_type, sort = TRUE) 

pdf("plots/05f-gcb_agorre_naive_downreg_hervs_type.pdf", height=4, width=4)
ggplot(downreg_agirre_naive_te_type, aes(fill=reorder(type, -n), y=cell_type, x=n)) + 
  geom_bar(position="stack", stat="identity", colour="black", size=0.3) + 
  scale_fill_manual(values = c(pal_aaas("default", alpha=0.7)(3))) + 
  coord_flip() +
  theme_pubclean() +  
  guides(fill=guide_legend(title="HERV Type")) +
  ylab("Downregulated HERVs by Type") +
  xlab("Number of HERVs") + 
  theme(legend.position = c("right")) +
  theme(aspect.ratio = 1)
dev.off()



################################### SAVE DATA ##################################


sig_herv_agirre <- sig_herv
sig_agirre <- sig
sig_agirre_no_pb <- sig_no_pb
sig_herv_agirre_no_pb <- sig_herv_no_pb
sig_agirre_naive <- sig_naive
sig_agirre_naive_herv <- sig_naive_herv

save(upvars_agirre, upvars_agirre_hervs, downvars_agirre, downvars_agirre_hervs,
     sig_herv_agirre, sig_agirre, file = "r_outputs/05f-agirre_vars.Rdata")

save(upvars_agirre_hervs_no_pb, upvars_agirre_hervs_no_pb, downvars_agirre_no_pb, downvars_agirre_hervs_no_pb,
     sig_herv_agirre_no_pb, sig_agirre_no_pb, file = "r_outputs/05f-agirre_vars_no_pb.Rdata")

save(upvars_agirre_naive, upvars_agirre_naive_herv, downvars_agirre_naive, downvars_agirre_naive_herv,
     sig_agirre_naive, sig_agirre_naive_herv, file = "r_outputs/05f-agirre_vars_naive.Rdata")

load("r_outputs/05f-agirre_vars.Rdata")
load("r_outputs/05f-agirre_vars_no_pb.Rdata")