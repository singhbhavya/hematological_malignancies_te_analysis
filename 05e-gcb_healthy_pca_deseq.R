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
library(edgeR)
library(ashr)
library(cowplot)
library(wesanderson)
library(UpSetR)

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
       legendPosition = "right")  +
  theme_cowplot()

## Biplot with projects (HERVs only)
biplot(GCB.herv.pca.obj, 
       lab = NULL,
       showLoadings = TRUE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "Cell_type",
       legendPosition = "right")  +
  theme_cowplot()

################################ SET THRESHOLDS ################################

fc=2 # fold change
l2fc=log2(fc) # log 2 fold change
lfc.cutoff <- 1.5
pval=0.05 # p value threshold

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



