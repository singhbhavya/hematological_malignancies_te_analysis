################################################################################
################################################################################
################################################################################
################################################################################
############################ GCB CELL TYPE MARKERS #############################

## Plan
## - Integration
## - Cell type annotation
## - Find markers
## - Create clusters
## - Azimuth 

#################################### SETUP #####################################

library(tidyverse)
library(Matrix)
library(scater)
library(Seurat)
library(SeuratData)
library(edgeR)
library(SeuratObject)
library(rtracklayer)
library(Azimuth)
library(patchwork)
library(cowplot)
library(data.table)

################################## LOAD DATA ###################################

load("r_outputs/03-gcb_seurat_objects_individual.Rdata")
load("r_outputs/03-gcb_seurat_objects_merged.Rdata")

############################# DETERMINE DIMENSIONS ##############################

determine_dim <- function(sobj_norm){
  sobj_norm <- FindVariableFeatures(object = sobj_norm, 
                                    selection.method = "vst", 
                                    nfeatures = 5000)
  sobj_norm <- ScaleData(object = sobj_norm)
  sobj_norm <- RunPCA(object = sobj_norm)
  sobj_norm <- JackStraw(sobj_norm, num.replicate = 100)
  sobj_norm <- ScoreJackStraw(sobj_norm, dims = 1:20)
  return(JackStrawPlot(sobj_norm, dims = 1:15))
}

SAMN13191512.seurat.p <- determine_dim(SAMN13191512.seurat.norm)+ggtitle("GCB 3B")
# All 15 PCs are significant
SAMN13191513.seurat.p <- determine_dim(SAMN13191513.seurat.norm)+ggtitle("GCB 3A")
# All 15 PCs are significant
SAMN14979745.seurat.p <- determine_dim(SAMN14979745.seurat.norm)+ggtitle("GCB 1")
# All 15 PCs are significant 
SAMN14979762.seurat.p <- determine_dim(SAMN14979762.seurat.norm)+ggtitle("GCB 2")
# All 15 PCa are significant 

############################### STANDARD SEURAT ################################

standard_seurat <- function(sobj_norm){
  sobj_norm <- FindVariableFeatures(object = sobj_norm, 
                                    selection.method = "vst", 
                                    nfeatures = 5000)
  sobj_norm <- ScaleData(object = sobj_norm)
  sobj_norm <- RunPCA(object = sobj_norm)
  sobj_norm <- FindNeighbors(object = sobj_norm)
  sobj_norm <- FindClusters(object = sobj_norm)
  sobj_norm <- RunTSNE(object = sobj_norm)
  sobj_norm <- RunUMAP(object = sobj_norm, dims=1:15)
  return(DimPlot(object = sobj_norm, reduction = 'umap'))
}

SAMN13191512.seurat.dim <- standard_seurat(SAMN13191512.seurat.norm)+ggtitle("GCB 3A")
SAMN13191513.seurat.dim <- standard_seurat(SAMN13191513.seurat.norm)+ggtitle("GCB 3B")
SAMN14979745.seurat.dim <- standard_seurat(SAMN14979745.seurat.norm)+ggtitle("GCB 1")
SAMN14979762.seurat.dim <- standard_seurat(SAMN14979762.seurat.norm)+ggtitle("GCB 2")

pdf("plots/04-individual_gcb_umaps.pdf", height=10, width=10)
plot_grid(SAMN14979762.seurat.dim, SAMN14979745.seurat.dim, 
          SAMN13191512.seurat.dim, SAMN13191513.seurat.dim)
dev.off()

GCB.norm.merged.dim <- standard_seurat(GCB.norm.merged)+ggtitle("All GCB Normalized and Merged")
GCB.merged.norm.dim <-  standard_seurat(GCB.merged.norm)+ggtitle("All GCB Merged and Normalized")

pdf("plots/04-combined_gcb_norm_merged_umap.pdf", height=5, width=5)
GCB.norm.merged.dim
dev.off()

pdf("plots/04-combined_gcb_merged_norm_umap.pdf", height=5, width=5)
GCB.merged.norm.dim
dev.off()

############################### AZIMUTH MAPPING ################################

mapped <- Azimuth::RunAzimuth(
  GCB.norm.merged,
  reference = "tonsilref",
  do.adt = TRUE,
)

saveRDS(mapped, file = "r_outputs/04-azimuth_mapping.Rds")

results <- list()
if('impADT' %in% Assays(mapped)) {
  results$impADT <- mapped[['impADT']]
}
if('ref.umap' %in% Reductions(mapped)) {
  results$umap <- mapped[['ref.umap']]
}

results$pred.df <- mapped@meta.data %>%
  tibble::rownames_to_column('cell') %>%
  dplyr::select(
    cell,
    dplyr::matches('predicted.celltype.l\\d$'),
    dplyr::matches('predicted.celltype.l\\d.score$'),
    mapping.score
  ) %>% as.data.frame

saveRDS(results, file = "r_outputs/04-azimuth_results.Rds")


########################### ADD AZIMUTH TO METADATA ############################

GCB.norm.merged <- Seurat::AddAzimuthResults(GCB.norm.merged, 
                                             file="r_outputs/04-azimuth_results.Rds")

# Examine metadata
GCB.norm.merged@meta.data %>% head

# the predicted celltype, score, and mapping.score columns are added to the
# metadata but are all NA. Add it manually here:
GCB.norm.merged@meta.data$predicted.celltype.l1 <- results$pred.df$predicted.celltype.l1
GCB.norm.merged@meta.data$predicted.celltype.l2 <- results$pred.df$predicted.celltype.l2
GCB.norm.merged@meta.data$predicted.celltype.l1.score <- results$pred.df$predicted.celltype.l1.score
GCB.norm.merged@meta.data$predicted.celltype.l2.score <- results$pred.df$predicted.celltype.l2.score
GCB.norm.merged@meta.data$mapping.score <- results$pred.df$mapping.score

# See that values are added
GCB.norm.merged@meta.data %>% head

saveRDS(GCB.norm.merged, file = 'r_outputs/04-GCB.norm.merged.mapped.Rds')

############################### FIND L1 MARKERS ################################

#

Idents(object=GCB.norm.merged) <- "predicted.celltype.l1"

l1.markers <- FindAllMarkers(
  GCB.norm.merged,
  test.use = 'wilcox',    # default = 'wilcox'
  only.pos = FALSE,        # default = FALSE
  min.pct = 0.1,         # default = 0.1
  logfc.threshold = 0.25, # default = 0.25
  return.thresh = 0.01    # default = 0.01
)

# Add feature metadata
fmeta <- GCB.norm.merged[['RNA']]@meta.features %>%
  tibble::rownames_to_column() %>%
  select(rowname, symbol, feattype, te_class, te_family)

stopifnot(all(l1.markers$gene %in% fmeta$rowname))
orig.rownames <- rownames(l1.markers)
l1.markers <- dplyr::left_join(l1.markers, fmeta, by=c('gene' = 'rowname'))
rownames(l1.markers) <- orig.rownames
l1.markers[is.na(l1.markers)] <- ''
rm(orig.rownames, fmeta)

saveRDS(l1.markers, file = 'r_outputs/04-GCB.norm.merged.l1.markers.rds')

############################### FIND L2 MARKERS ################################

Idents(object=GCB.norm.merged) <- "predicted.celltype.l2"

l2.markers <- FindAllMarkers(
  GCB.norm.merged,
  test.use = 'wilcox',    # default = 'wilcox'
  only.pos = FALSE,        # default = FALSE
  min.pct = 0.1,         # default = 0.1
  logfc.threshold = 0.25, # default = 0.25
  return.thresh = 0.01    # default = 0.01
)

# Add feature metadata
fmeta <- GCB.norm.merged[['RNA']]@meta.features %>%
  tibble::rownames_to_column() %>%
  select(rowname, symbol, feattype, te_class, te_family)

stopifnot(all(l2.markers$gene %in% fmeta$rowname))
orig.rownames <- rownames(l2.markers)
l2.markers <- dplyr::left_join(l2.markers, fmeta, by=c('gene' = 'rowname'))
rownames(l2.markers) <- orig.rownames
l2.markers[is.na(l2.markers)] <- ''
rm(orig.rownames, fmeta)

saveRDS(l2.markers, file = 'r_outputs/04-GCB.norm.merged.l2.markers.rds')


################################# PLOT UMAPs ###################################

umap1_lim <- range(GCB.norm.merged[['umap.proj']]@cell.embeddings[,1])
umap2_lim <- range(GCB.norm.merged[['umap.proj']]@cell.embeddings[,2])

pdf('plots/04-GCB.norm.merged.mapped.l1.pdf', width=10, height=10)
nclust <- length(unique(GCB.norm.merged@meta.data$predicted.celltype.l1))
p <- DimPlot(GCB.norm.merged,
             repel=TRUE,
             group.by = "predicted.celltype.l1",
             label = TRUE,
             label.size = 9,
             cols=Seurat::DiscretePalette(nclust, 'polychrome'),
             raster = FALSE,
             raster.dpi = c(1200, 1200)
)
p + xlim(umap1_lim) + ylim(umap2_lim) + NoLegend() + theme(plot.title=element_blank())
dev.off()

nclust <- length(unique(GCB.norm.merged@meta.data$predicted.celltype.l2))
p_l2 <- DimPlot(GCB.norm.merged,
        repel = TRUE,
        label.size = 4,
        pt.size = 1,
        group.by = "predicted.celltype.l2",
        label = TRUE,
        label.box = FALSE,
        cols=c(Seurat::DiscretePalette(32, c('glasbey')),
               Seurat::DiscretePalette(36, c('polychrome'))),
        raster = FALSE) 

pdf('plots/04-GCB.norm.merged.mapped.l2_nolab.pdf', width=12, height=8)
p_l2 + theme_cowplot() +
  xlim(umap1_lim) + ylim(umap2_lim) + NoLegend() +  theme(plot.title=element_blank()) 
dev.off()

pdf('plots/04-GCB.norm.merged.mapped.l2_lab.pdf', width=15, height=8) 
p_l2 + xlim(umap1_lim) + ylim(umap2_lim) + theme_cowplot() + 
  theme(plot.title=element_blank()) +
  guides(color = guide_legend(override.aes = list(size=5), ncol=2) ) 
dev.off()

############################## L1 FEATURE PLOTS ################################
# condition_pct.1: percentage of cells where the gene is detected in the cluster 
# for condition
# condition_pct.2: percentage of cells where the gene is detected on average in 
# the other clusters for condition


herv_markers_l1 <- l1.markers[l1.markers$te_class=='LTR',]

nrow(herv_markers_l1) # There are 57 HERV markers.
length(unique(herv_markers_l1$gene)) # There are 26 unique HERVs.

## CD4 T: MER41-15q21.2a
## Natural killer: HERVH-12p13.31d, HARLEQUIN-17q12, MER4-22q12.3, 
##                 HUERSP1-11p15.4b, ERV316A3-2q22.2b
## Naive B: ERV316A3-2q22.2b
## Monocyte/macrophage: MER34B-1q23.3a, MER34B-4q31.3, HARLEQUIN-17q25.3a, 
##                      ERVLB4-1p34.2a, MER4B-7p14.3
## Naive CD4 T: ERV316A3-2q22.2b
## Activated Naive B: HARLEQUIN-1q32.1
## Memory B: ERV316A3-2q22.2b
## Plasmacytoid dendritic: MER4-10q23.31b, HERV3-Yq11.223b, HARLEQUIN-1q32.1
## Plasma: HUERSP2-19q13.2, HARLEQUIN-1q32.1
## Cycling B: ERV316A3-8q13.3a
## Germinal center B: HARLEQUIN-1q32.1

herv_markers_l2 <- l2.markers[l2.markers$te_class=='LTR',]

nrow(herv_markers_l2) # There are 144 HERV markers.
length(unique(herv_markers_l2$gene)) # There are 61 unique HERVs.

############################## L2 FEATURE PLOTS ################################

fpcols <- c('#eeeeeeFF', viridis::viridis(6))

## Plasmablast 
PB_all <- herv_markers_l2[herv_markers_l2$cluster == 'PB', 'gene']
lapply(PB_all, function(x) herv_markers_l2[herv_markers_l2$gene == x,])
b_x <- c('HML6-19q13.43b', 'ERV316A3-8q13.3a')
pdf('plots/04-gcb_combined_l2_herv_markers_PB.pdf', width=6, height=3)
p <- FeaturePlot(GCB.norm.merged, b_x, cols=fpcols, ncol=2, raster=FALSE)
p & xlim(umap1_lim) & ylim(umap2_lim) & NoAxes() & NoLegend() & theme(plot.title=element_text(size=8)) 
dev.off()

## HARLEQUIN-1q32.1
pdf('plots/04-gcb_combined_l2_HARLEQUIN-1q32.1.pdf', width=3, height=3)
p <- FeaturePlot(GCB.norm.merged, "HARLEQUIN-1q32.1", cols=fpcols, raster=FALSE)
p & xlim(umap1_lim) & ylim(umap2_lim) & NoAxes() & NoLegend() & theme(plot.title=element_text(size=8))
dev.off()

## PDC (Plasmacytoid dendritic cells)
PDC_all <- herv_markers_l2[herv_markers_l2$cluster == 'PDC', 'gene']
lapply(PDC_all, function(x) herv_markers_l2[herv_markers_l2$gene == x,])
b_x <- c("ERV316A3-2q22.2b", "ERVLE-4q24e", "HARLEQUIN-1q32.1", 
         "HARLEQUIN-10q23.1", "HML1-1q32.1", "HERVEA-5q22.2", "ERV316A3-8q13.3a" )
pdf('plots/04-gcb_combined_l2_herv_markers_PDC.pdf', width=9, height=9)
p <- FeaturePlot(GCB.norm.merged, b_x, cols=fpcols, ncol=3, raster=FALSE)
p & xlim(umap1_lim) & ylim(umap2_lim) & NoAxes() & NoLegend() & theme(plot.title=element_text(size=8))
dev.off()

## All memory B cell clusters
MB_all <- herv_markers_l2[herv_markers_l2$cluster  %like%  'MB', 'gene']
lapply(MB_all, function(x) herv_markers_l2[herv_markers_l2$gene == x,])
b_x <- c("ERV316A3-2q22.2b", "ERVLE-4q24e", "ERV316A3-8q13.3a", 
         "HARLEQUIN-1q32.1", "HERVS71-19q13.12a", "MER101-12p13.31b")
pdf('plots/04-gcb_combined_l2_herv_markers_MB.pdf', width=9, height=6)
p <- FeaturePlot(GCB.norm.merged, b_x, cols=fpcols, ncol=3, raster=FALSE)
p & xlim(umap1_lim) & ylim(umap2_lim) & NoAxes() & NoLegend() & theme(plot.title=element_text(size=8))
dev.off()

## HARLEQUIN-17q25
pdf('plots/04-gcb_combined_l2_HARLEQUIN-17q25.pdf', width=3, height=3)
p <- FeaturePlot(GCB.norm.merged, c("HARLEQUIN-17q25.3a", "HARLEQUIN-17q25.3b"), 
                 cols=fpcols, ncol=2, raster=FALSE)
p & xlim(umap1_lim) & ylim(umap2_lim) & NoAxes() & theme(plot.title=element_text(size=8))
dev.off()

################################# SAVE FILES ###################################

save(GCB.norm.merged, mapped, results,
     l1.markers, l2.markers, 
     herv_markers_l1, herv_markers_l2, GCB_metadata,
     file="r_outputs/04-gcb.norm.merged.azimuth.Rdata")

