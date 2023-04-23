################################################################################
################################################################################
################################################################################
################################################################################
######################### LYMPHOMA / B-CELL REGRESSION ######################### 

## Plan:
##
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

#################################### SETUP #####################################

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

############################# SUBTYPE DESEQ GENES ##############################

all.countdat <- all.counts.filt.comb
cat(sprintf('%d variables\n', nrow(all.countdat)))

stopifnot(all(colnames(all.countdat) == rownames(all_metadata)))

all.g.dds <- DESeq2::DESeqDataSetFromMatrix(countData = all.countdat,
                                          colData = all_metadata,
                                          design = ~ subtype + 0)

all.g.dds <- DESeq2::DESeq(all.g.dds, parallel=T)
all.g.tform <- DESeq2::varianceStabilizingTransformation(all.g.dds, blind=FALSE)


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

############################# TOP GENES AND HERVS ############################## 

all_res <- list(
  "DLBCL ABC" = DESeq2::results(all.g.dds, contrast=c(+1, -1/7, -1/7, -1/7, -1/7, -1/7, 
                                                    -1/7, -1/7), alpha=pval),
  "BL Endemic EBV Negative" = DESeq2::results(all.g.dds, contrast=c(-1/7, +1, -1/7, -1/7, -1/7, -1/7, 
                                                                  -1/7, -1/7), alpha=pval),
  "BL Endemic EBV Positive" = DESeq2::results(all.g.dds, contrast=c(-1/7, -1/7, +1, -1/7, -1/7, -1/7, 
                                                                  -1/7, -1/7), alpha=pval),
  "FL" = DESeq2::results(all.g.dds, contrast=c(-1/7, -1/7, -1/7, +1, -1/7, -1/7, 
                                             -1/7, -1/7), alpha=pval),
  "DLBCL GCB" = DESeq2::results(all.g.dds, contrast=c(-1/7, -1/7, -1/7, -1/7, +1, -1/7, 
                                                    -1/7, -1/7), alpha=pval),
  "BL Sporadic EBV Negative" = DESeq2::results(all.g.dds, contrast=c(-1/7, -1/7, -1/7, -1/7, -1/7, +1, 
                                                                   -1/7, -1/7), alpha=pval),
  "BL Sporadic EBV Positive" = DESeq2::results(all.g.dds, contrast=c(-1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 
                                                                   +1, -1/7), alpha=pval),
  "DLBCL Unclass" = DESeq2::results(all.g.dds, contrast=c(-1/7, -1/7, -1/7, -1/7, -1/7, -1/7, 
                                                        -1/7, +1), alpha=pval)
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

for (n in names(sig)) {
  cat("\n#--- Contrast", n, "---#\n")
  print(head(sig[[n]]))
}

upvars.g <- lapply(sig[1:8], function(r) rownames(subset(r, log2FoldChange>0)))
downvars.g <- lapply(sig[1:8], function(r) rownames(subset(r, log2FoldChange<0)))
allvars.g <- lapply(sig[1:8], function(r) rownames(subset(r, log2FoldChange>0 | log2FoldChange<0)))


############################ UPSET GENES & HERVS ###############################

upvars <- lapply(sig_herv[1:8], function(r) rownames(subset(r, log2FoldChange>0)))
downvars <- lapply(sig_herv[1:8], function(r) rownames(subset(r, log2FoldChange<0)))
allvars <- lapply(sig_herv[1:8], function(r) rownames(subset(r, log2FoldChange>0 | log2FoldChange<0)))

pdf("plots/05p-all_lymph_subtype_upset_upvars.pdf", height=7, width=11)
UpSetR::upset(fromList(upvars), sets=c(names(sig_herv)),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1.7))
dev.off()

pdf("plots/05p-all_lymph_subtype_compexupset_upvars.pdf", height=4, width=6)
ComplexUpset::upset(fromList(upvars), 
                    intersect = c(names(upvars)),
                    intersections = list( 
                      c("DLBCL ABC"), 
                      c("BL Endemic EBV Negative"), 
                      c("BL Endemic EBV Positive"),
                      c("FL"), 
                      c("DLBCL GCB"),
                      c("BL Sporadic EBV Negative"),
                      c("BL Sporadic EBV Positive"),
                      c("DLBCL Unclass")), 
                    queries = list(
                      upset_query(set=c("DLBCL ABC"), 
                                  color = "red3", 
                                  fill = "red3"),
                      upset_query(set=c("BL Endemic EBV Negative"), 
                                  color = "#ae5c00", 
                                  fill = "#ae5c00"),
                      upset_query(set=c("BL Endemic EBV Positive"), 
                                  color = wes_palette("Darjeeling1")[2], 
                                  fill = wes_palette("Darjeeling1")[2]),
                      upset_query(set=c("FL"), 
                                  color = "#868686B2", 
                                  fill = "#868686B2"),
                      upset_query(set=c("DLBCL GCB"), 
                                  color = "royalblue", 
                                  fill = "royalblue"),
                      upset_query(set=c("BL Sporadic EBV Negative"), 
                                  color = wes_palette("Darjeeling1")[4], 
                                  fill = wes_palette("Darjeeling1")[4]),
                      upset_query(set=c("BL Sporadic EBV Positive"), 
                                  color = "#006053", 
                                  fill = "#006053"),
                      upset_query(set=c("DLBCL Unclass"), 
                                  color = "lightblue", 
                                  fill = "lightblue")
                    ),
                    set_sizes=(upset_set_size())
                      + theme(axis.text.x=element_text(angle=45)))

dev.off()

pdf("plots/05p-all_lymph_subtype_upset_dnwars.pdf", height=7, width=11)
UpSetR::upset(fromList(downvars), sets=c(names(sig_herv)),  
      keep.order = T, order.by='degree', decreasing=F,
      text.scale = c(1.5, 1.5, 1, 1, 1.7, 1.7))
dev.off()

pdf("plots/05p-all_lymph_subtype_complexupset_dnwars.pdf", height=7, width=11)
ComplexUpset::upset(fromList(downvars), 
                    intersect = c(names(downvars)),
                    intersections = list( 
                      c("DLBCL ABC"), 
                      c("BL Endemic EBV Negative"), 
                      c("BL Endemic EBV Positive"),
                      c("FL"), 
                      c("DLBCL GCB"),
                      c("BL Sporadic EBV Negative"),
                      c("BL Sporadic EBV Positive"),
                      c("DLBCL Unclass")), 
                    queries = list(
                      upset_query(set=c("DLBCL ABC"), 
                                  color = "red3", 
                                  fill = "red3"),
                      upset_query(set=c("BL Endemic EBV Negative"), 
                                  color = "#ae5c00", 
                                  fill = "#ae5c00"),
                      upset_query(set=c("BL Endemic EBV Positive"), 
                                  color = wes_palette("Darjeeling1")[2], 
                                  fill = wes_palette("Darjeeling1")[2]),
                      upset_query(set=c("FL"), 
                                  color = "#868686B2", 
                                  fill = "#868686B2"),
                      upset_query(set=c("DLBCL GCB"), 
                                  color = "royalblue", 
                                  fill = "royalblue"),
                      upset_query(set=c("BL Sporadic EBV Negative"), 
                                  color = wes_palette("Darjeeling1")[4], 
                                  fill = wes_palette("Darjeeling1")[4]),
                      upset_query(set=c("BL Sporadic EBV Positive"), 
                                  color = "#006053", 
                                  fill = "#006053"),
                      upset_query(set=c("DLBCL Unclass"), 
                                  color = "lightblue", 
                                  fill = "lightblue")
                    ))

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

################################### PATHWAYS ###################################

pathways.hallmark <- gmtPathways("gsea/h.all.v2023.1.Hs.symbols.gmt")
pathways.immune <- gmtPathways("gsea/c7.immunesigdb.v2023.1.Hs.symbols.gmt")
pathways.bp <- gmtPathways("gsea/c5.go.bp.v2023.1.Hs.symbols.gmt")
pathways.kegg <- gmtPathways("gsea/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt")
pathways.biocarta <- gmtPathways("gsea/c2.cp.biocarta.v2023.1.Hs.symbols.gmt.txt")

################################ FGSEA FUNCTION ################################

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
  
  return(fgsea.out)
  
  assign(paste0(clust_name, ".", pathway_name, ".fgsea.out"), fgsea.out, envir = .GlobalEnv )
}

################################# B CELL FGSEA #################################

load("r_outputs/05f-agirre_vars.Rdata")

gene.sets.b.cell <- lapply(upvars_agirre[1:6], function(r) gene_table[r[1:150],]$display)
herv.sets.b.cell <- lapply(upvars_agirre_hervs[1:6], function(r) r[1:25])

gene.herv.sets.b.cell <- list("DZ" = append(gene.sets.b.cell$DZ, herv.sets.b.cell$DZ),
                              "LZ" = append(gene.sets.b.cell$LZ, herv.sets.b.cell$LZ),
                              "PB" = append(gene.sets.b.cell$PB, herv.sets.b.cell$PB),
                              "BMPC" = append(gene.sets.b.cell$BMPC, herv.sets.b.cell$BMPC),
                              "NB" = append(gene.sets.b.cell$NB, herv.sets.b.cell$NB),
                              "MB" = append(gene.sets.b.cell$MB, herv.sets.b.cell$MB)
)

b.cell.all.lymph <- list(
  "DLBCL ABC" = make.fsgsea(gene.herv.sets.b.cell, all_res$`DLBCL ABC`, "DLBCL ABC", "b.cell"),
  "BL Endemic EBV Negative" = make.fsgsea(gene.herv.sets.b.cell, all_res$`BL Endemic EBV Negative`, "BL Endemic EBV Negative", "b.cell"),
  "BL Endemic EBV Positive" = make.fsgsea(gene.herv.sets.b.cell, all_res$`BL Endemic EBV Positive`, "BL Endemic EBV Positive", "b.cell"),
  "FL" = make.fsgsea(gene.herv.sets.b.cell, all_res$FL, "FL", "b.cell"),
  "DLBCL GCB" = make.fsgsea(gene.herv.sets.b.cell, all_res$`DLBCL GCB`, "DLBCL GCB", "b.cell"),
  "BL Sporadic EBV Negative" = make.fsgsea(gene.herv.sets.b.cell, all_res$`BL Sporadic EBV Negative`, "BL Sporadic EBV Negative", "b.cell"),
  "BL Sporadic EBV Positive" = make.fsgsea(gene.herv.sets.b.cell, all_res$`BL Sporadic EBV Positive`, "BL Sporadic EBV Positive", "b.cell"),
  "DLBCL Unclass" = make.fsgsea(gene.herv.sets.b.cell, all_res$`DLBCL Unclass`, "DLBCL Unclass", "b.cell")
)

pdf("plots/05p-all_lymph_bcell.pdf", height=10, width=7)
for(clust in names(b.cell.all.lymph)) {
  fgseaResTidy <- b.cell.all.lymph[[clust]] %>%
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

# Get longform df
b.cell.all.lymph.summary <- rbindlist(b.cell.all.lymph, idcol = "index")


pdf("plots/05p-all_lymph_bcell_bubble.pdf", height=5, width=6)
ggplot(b.cell.all.lymph.summary, aes(x = index, 
                                     y = pathway, 
                                     size = -log(padj), 
                                     color = NES)) +
  geom_point() +
  scale_size(name = "-log (P value)", range = c(1, 10)) + 
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))+
  scale_colour_gradientn(colors = viridis_pal()(10)) +
  xlab("Lymphoma Subtype") +
  ylab("B Cell Signature") +
  theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"))
dev.off()


################################ FAMILY LEVEL UP ############################### 

upreg.hervs.df <- do.call(rbind, lapply(upvars, data.frame))
colnames(upreg.hervs.df) <- c("herv")
upreg.hervs.df$cancer_type <- rownames(upreg.hervs.df)
upreg.hervs.df$cancer_type <- gsub("\\..*","",upreg.hervs.df$cancer_type)
upreg.hervs.df$family <- retro.annot$family[match(upreg.hervs.df$herv, 
                                                         retro.annot$locus)]

upreg.families <-
  upreg.hervs.df %>% dplyr::count(family, cancer_type, sort = TRUE) 

down.hervs.df <- do.call(rbind, lapply(downvars, data.frame))
colnames(down.hervs.df) <- c("herv")
down.hervs.df$cancer_type <- rownames(down.hervs.df)
down.hervs.df$cancer_type <- gsub("\\..*","",down.hervs.df$cancer_type)
down.hervs.df$family <- retro.annot$family[match(down.hervs.df$herv, 
                                                        retro.annot$locus)]

downreg.families <-
  down.hervs.df %>% dplyr::count(family, cancer_type, sort = TRUE) 

updown.family <- do.call(rbind, (list(Upregulated = upreg.families, 
                                      Dowregulated = downreg.families)))
updown.family$expression <- rownames(updown.family)
updown.family$expression <- gsub("\\..*","",updown.family$expression)
rownames(updown.family)<-NULL

pdf("plots/05p-lymphoma_updown_families_all.pdf", height=8, width=8)
ggplot(updown.family, aes(fill=reorder(family, -n), y=cancer_type, x=n)) + 
  geom_bar(position="stack", stat="identity", colour="black", size=0.3) + 
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
  theme_cowplot() +  
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  guides(fill=guide_legend(title="TE Family")) +
  ylab(NULL) +
  xlab("Number of HERV Loci") + 
  theme(legend.position = c("right"),
        plot.margin = margin(10, 10, 10, 40),
        axis.line=element_blank()) + 
  guides(fill = guide_legend(title = "HERV family", ncol = 2)) +
  facet_wrap(~ expression, ncol = 2)
dev.off()


################################# PLOT COUNTS ################################## 


plot.counts <- function(df, gene, title) {
  
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="subtype", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=subtype, y=count, fill=subtype))  +
    geom_boxplot(notch = TRUE) +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("subtype") +
    ylab("Counts") +
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1) +
    scale_fill_manual(values = c("Endemic BL EBV-positive" = wes_palette("Darjeeling1")[2], 
                                 "Sporadic BL EBV-positive" = "#006053",
                                 "Sporadic BL EBV-negative" = wes_palette("Darjeeling1")[4],
                                 "Endemic BL EBV-negative" = "#ae5c00",
                                 "GCB" = "royalblue", 
                                 "ABC" = "red3", 
                                 "Unclass" = "lightblue", 
                                 "Missing" = "grey",
                                 "FOLLICULAR GRADE 1" = "#F1BB7B",
                                 "FOLLICULAR GRADE 2" = "#FD6467",
                                 "FOLLICULAR GRADE 3A" = "#5B1A18")) + 
    scale_y_log10(labels = label_comma()) 
}

plot.counts(all.g.dds, "ENSG00000143379.12", "SETDB1") +
  stat_compare_means(comparisons = list(c("ABC", "Endemic BL EBV-positive")))
#SETDB1
plot.counts(all.g.dds, "ENSG00000182578.14", "CSF1R") +
  stat_compare_means(comparisons = c("ABC", "Endemic BL EBV-positive"))#CSF1R


