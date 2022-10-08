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

################################ FIND MARKERS ##################################

Idents(object=GCB.norm.merged) <- "predicted.celltype.l1"

l1.markers <- FindAllMarkers(
  GCB.norm.merged,
  test.use = 'wilcox',    # default = 'wilcox'
  only.pos = TRUE,        # default = FALSE
  min.pct = 0.25,         # default = 0.1
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


