################################################################################
################################################################################
################################################################################
################################################################################
########################## BURKITT LYMPHOMA CLUSTERS ###########################

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
library(glmnet)
library(Boruta)
library(UpSetR)
library(c060)
library(randomForest)
library(rpart)
library(rpart.plot)
library(ConsensusClusterPlus)

################################### LOAD DATA ##################################

load("r_outputs/01-refs.Rdata")
load("r_outputs/05b-BL_pca_dds.Rdata")
load("r_outputs/01-metadata.Rdata")

##################################### PCA ###################################### 

# Just re-plotting here to remind myself what the PCA plots look like
biplot(BL.herv.pca.obj, 
       lab = NULL,
       showLoadings = TRUE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "ebv_status",
       shape = "clinical_variant", 
       shapekey = c("Endemic BL" = 15, "Sporadic BL" = 8),
       colkey = c("EBV-positive" = wes_palette("Darjeeling1")[2], 
                  "EBV-negative" = wes_palette("Darjeeling1")[4]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

biplot(BL.g.pca.obj, 
       lab = NULL,
       showLoadings = TRUE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "ebv_status",
       shape = "clinical_variant", 
       shapekey = c("Endemic BL" = 15, "Sporadic BL" = 8),
       colkey = c("EBV-positive" = wes_palette("Darjeeling1")[2], 
                  "EBV-negative" = wes_palette("Darjeeling1")[4]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

################################## CLUSTERING ##################################

tDat <- assay(BL.tform)[BL.herv.pca.obj$xvars, ]
cDat <- sweep(tDat, 1, apply(tDat, 1, median, na.rm=T))

BL.ccp.obj <- ConsensusClusterPlus(cDat, maxK=3, reps=1000, pItem=0.8, pFeature=0.8,
                                title="plots", clusterAlg="km", distance='euclidean',
                                seed=12345, plot="pdf")
icl <- calcICL(BL.ccp.obj, title="plots", plot='pdf')

stopifnot(all(rownames(BL.herv.pca.obj$metadata) == names(BL.ccp.obj[[2]]$consensusClass)))

clust.df <- lapply(3, function(k) {
  factor(paste0('C', BL.ccp.obj[[k]]$consensusClass), levels=paste0('C', 1:k))
})%>% bind_cols() %>% data.frame
colnames(clust.df) <- "clust.retro.k3"
rownames(clust.df) <- names(BL.ccp.obj[[2]]$consensusClass)

BL.herv.pca.obj$metadata <- cbind(BL.herv.pca.obj$metadata, clust.df)
BL.herv.pca.obj$metadata <- cbind(BL.herv.pca.obj$metadata, BL_metadata$MYC_SV)
BL.herv.pca.obj$metadata <- cbind(BL.herv.pca.obj$metadata, BL_metadata$MYC_SV_partner)
BL.herv.pca.obj$metadata <- cbind(BL.herv.pca.obj$metadata, BL_metadata$Total_N_SSM)
BL.herv.pca.obj$metadata <- cbind(BL.herv.pca.obj$metadata, BL_metadata$subgroup)
colnames(BL.herv.pca.obj$metadata) <- c("case", "project_id", "submitter_id",
                                        "sample_type", "tissue_type", "tumor_descriptor",
                                        "cohort", "clinical_variant", "ebv_status",
                                        "ebv_genome_type", "gender", "age_at_diagnosis",
                                        "anatomic_site_classification", "tissue_source_site",
                                        "cancer_type", "subtype", "clust.retro.k3", 
                                        "MYC_SV", "MYC_SV_partner", "Total_N_SSM",
                                        "subgroup")

wrap <- function(nclu) {
  columnname <- paste0('clust.retro.k', nclu)
  
  gp <- PCAtools::biplot(BL.herv.pca.obj,
                         colby = columnname,
                         colkey = c("C1" = wes_palette("Chevalier1")[1], 
                                    "C2" = wes_palette("Chevalier1")[2],
                                    "C3" = wes_palette("Chevalier1")[3]),
                         ellipse = TRUE,
                         ellipseConf = 0.9,
                         shape='ebv_status',
                         hline = 0,
                         vline = 0,
                         legendPosition = 'right',
                         lab=NULL,
                         title=paste0("Retrotranscriptome CC PCA")
  ) + theme_cowplot() +
    theme(aspect.ratio=1)
  print(gp)
  print(table(BL.herv.pca.obj$metadata[ ,columnname], 
              BL.herv.pca.obj$metadata$ebv_status))
}

pdf("plots/05i-BL_ccp_PCA.pdf", height=7, width=7)
wrap(3)
dev.off()

cplots <- lapply(3, function(nclu) {
  p1 <- factoextra::fviz_cluster(
    list(data=t(cDat), cluster=BL.ccp.obj[[nclu]]$consensusClass),
    ellipse.type="norm", ellipse.level=0.9,
    geom="point",
    palette = "npg"
  ) + theme_cowplot()
  sil <- cluster::silhouette(BL.ccp.obj[[nclu]]$consensusClass, dist(t(cDat), method = "euclidean"))
  p2 <- factoextra::fviz_silhouette(
    cluster::silhouette(BL.ccp.obj[[nclu]]$consensusClass, dist(t(cDat), method = "euclidean")),
    palette = "npg"
  )
  return(list(p1,p2))
})
cplots <- unlist(cplots, recursive = FALSE)

pdf("plots/05i-BL_clusters.pdf", height=5, width=10)
gridExtra::grid.arrange(grobs = cplots, ncol=2)
dev.off()

sil.score <- sapply(3, function(k) {
  mean(cluster::silhouette(BL.ccp.obj[[k]]$consensusClass, 
                           dist(t(cDat),
                                method = "euclidean"))[,3])
})

pdf("plots/05i-BL_eigen.pdf", height=8, width=10)
PCAtools::eigencorplot(BL.herv.pca.obj,
                       metavars = c('ebv_status','clinical_variant','subgroup','gender',
                                    'tissue_source_site','MYC_SV_partner',
                                    'subtype', 'Total_N_SSM', 'anatomic_site_classification',
                                    "clust.retro.k3"),
                       col=wes_palette("Zissou1", 12, type = "continuous"),
                       signifSymbols = c('****', '***', '**', '*', ''),
                       signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))
dev.off()

# Add clusters to metadata
BL_metadata$clust.retro.k3 <- clust.df

################################## SAVE FILES ##################################

save(BL.herv.pca.obj, BL.ccp.obj, tDat, cDat, clust.df, BL_metadata,
     file="r_outputs/05i-BL_pca_ccp_clusters_metadata.Rdata")

