################################################################################
################################################################################
################################################################################
################################################################################
######################### LYMPHOMA / B-CELL REGRESSION ######################### 

## Plan:
##

#################################### SETUP #####################################

load("r_outputs/05f-agirre_vars_naive.Rdata")
load("r_outputs/05l-all_lymphoma_vars.Rdata")  
load("r_outputs/02-all_lymphoma_filt_counts.Rdata")
load("r_outputs/01-refs.Rdata")

################################ SET THRESHOLDS ################################

fc=2 # fold change
l2fc=log2(fc) # log 2 fold change
lfc.cutoff <- 1.5
pval=0.001 # p value threshold

################################ SUBTYPE DESEQ  ################################

all_metadata <- all_metadata %>% replace(is.na(.), "Unclass")
all_metadata$subtype[all_metadata$subtype == "Missing"] <- "Unclass"
all_metadata$subtype[all_metadata$subtype == "FOLLICULAR GRADE 1"] <- "FL"
all_metadata$subtype[all_metadata$subtype == "FOLLICULAR GRADE 3A"] <- "FL"
all_metadata$subtype[all_metadata$subtype == "FOLLICULAR GRADE 2"] <- "FL"

hervs.to.keep <- intersect(rownames(all.counts.filt.herv), 
                           retro.annot$locus[retro.annot$chrom != "chrY"])

all.counts.filt.herv.y <- all.counts.filt.herv[hervs.to.keep,] 

all.countdat <- all.counts.filt.herv.y
cat(sprintf('%d variables\n', nrow(all.countdat)))

stopifnot(all(colnames(all.countdat) == rownames(all_metadata)))

all.dds <- DESeq2::DESeqDataSetFromMatrix(countData = all.countdat,
                                          colData = all_metadata,
                                          design = ~ subtype + 0)

all.dds <- DESeq2::DESeq(all.dds, parallel=T)
all.tform <- DESeq2::varianceStabilizingTransformation(all.dds, blind=FALSE)

############################### TOP HERVs ONLY #################################

all_res_herv <- list(
  "DLBCL ABC" = DESeq2::results(all.dds, contrast=c(+1, -1/7, -1/7, -1/7, -1/7, -1/7, 
                                                    -1/7, -1/7), alpha=pval),
  "BL Endemic EBV Negative" = DESeq2::results(all.dds, contrast=c(-1/7, +1, -1/7, -1/7, -1/7, -1/7, 
                                                                  -1/7, -1/7), alpha=pval),
  "BL Endemic EBV Positive" = DESeq2::results(all.dds, contrast=c(-1/7, -1/7, +1, -1/7, -1/7, -1/7, 
                                                                  -1/7, -1/7), alpha=pval),
  "FL" = DESeq2::results(all.dds, contrast=c(-1/7, -1/7, -1/7, +1, -1/7, -1/7, 
                                             -1/7, -1/7), alpha=pval),
  "DLBCL GCB" = DESeq2::results(all.dds, contrast=c(-1/7, -1/7, -1/7, -1/7, +1, -1/7, 
                                             -1/7, -1/7), alpha=pval),
  "BL Sporadic EBV Negative" = DESeq2::results(all.dds, contrast=c(-1/7, -1/7, -1/7, -1/7, -1/7, +1, 
                                                    -1/7, -1/7), alpha=pval),
  "BL Sporadic EBV Positive" = DESeq2::results(all.dds, contrast=c(-1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 
                                                                   +1, -1/7), alpha=pval),
  "DLBCL Unclass" = DESeq2::results(all.dds, contrast=c(-1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 
                                                                   -1/7, +1), alpha=pval)
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

for (n in names(sig_herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  summary(sig_herv[[n]])
}

for (n in names(sig_herv)) {
  cat("\n#--- Contrast", n, "---#\n")
  print(head(sig_herv[[n]]))
}


############################ UPSET GENES & HERVS ###############################

upvars <- lapply(sig_herv[1:8], function(r) rownames(subset(r, log2FoldChange>0)))
downvars <- lapply(sig_herv[1:8], function(r) rownames(subset(r, log2FoldChange<0)))
allvars <- lapply(sig_herv[1:8], function(r) rownames(subset(r, log2FoldChange>0 | log2FoldChange<0)))

pdf("plots/05p-all_lymph_subtype_upset_upvars.pdf", height=7, width=11)
upset(fromList(upvars), sets=c(names(sig_herv)),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1.7))
dev.off()

pdf("plots/05p-all_lymph_subtype_upset_dnwars.pdf", height=7, width=11)
upset(fromList(downvars), sets=c(names(sig_herv)),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1.7))
dev.off()

up.binmat <- fromList(upvars)
rn <- do.call(c, upvars)
rn <- rn[!duplicated(rn)]
rownames(up.binmat) <- rn
rm(rn)

up.binmat.unique <- subset(up.binmat, rowSums(up.binmat, na.rm = TRUE) == 1)

dn.binmat <- fromList(downvars)
rn <- do.call(c, downvars)
rn <- rn[!duplicated(rn)]
rownames(dn.binmat) <- rn
rm(rn)

dn.binmat.unique <- subset(dn.binmat, rowSums(dn.binmat, na.rm = TRUE) == 1)

all.binmat <- fromList(allvars)
rn <- do.call(c, allvars)
rn <- rn[!duplicated(rn)]
rownames(all.binmat) <- rn
rm(rn)

all.binmat.unique <- subset(all.binmat, rowSums(all.binmat, na.rm = TRUE) == 1)

intersect(rownames(all.binmat.unique), sig_herv$`DLBCL ABC`$display)

cor.test(samp_tes[,2], samp_tes[,3], 
         method="spearman", exact = FALSE)


################################# HEATMAPS #####################################
################################### SETUP ######################################

# Colors for cells
cols <- rgb_gsea(palette = c("default"), n = 14, alpha = 0.7, reverse = FALSE)

# Annotation column
df <- as.data.frame(colData(all.dds)[,c("cancer_type","subtype")])

# Create colors for each group
annoCol <-  pal_jco("default", alpha = 0.7)(3)
names(annoCol) <- unique(df$cancer_type)
annoCol2 <- c("red3", "lightblue", "royalblue", wes_palette("Darjeeling1")[2], 
              "#ae5c00", wes_palette("Darjeeling1")[4], "#006053", pal_jco("default", alpha = 0.7)(3)[3])
names(annoCol2) <- unique(df$subtype)
annoCol <- list(cancer_type = annoCol, subtype = annoCol2)

######################### UPREGULATED IN ALL GROUPS ############################
top.hervs <- rownames(up.binmat)

pdf("plots/05p-all_lymphoma_subtype_top_hervs_upregulated_all.pdf", height=10, width=10)
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

for(clust in names(sig_herv)) {
  tg <- rownames(sig_herv[[clust]][1:75,])
  p <- makeheatmap(tg, main=paste0('DE in cluster ', clust))
  pdf(paste0("plots/05p-all_lymphoma_subtype_top_de_hervs_", clust, ".pdf"), height=7, width=7)
  print(p)
  dev.off()
}

