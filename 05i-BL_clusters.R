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
load("r_outputs/02-BL_filt_counts.Rdata")
load("r_outputs/01-metadata.Rdata")

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
################################### ALL HERVs ##################################
################################################################################
################################################################################

################################ BL DESEQ HERVs ################################

### DESeq2 (HERVs Only)

BL.countdat <- BL.filt.herv
cat(sprintf('%d variables\n', nrow(BL.countdat)))

stopifnot(all(colnames(BL.countdat) == rownames(BL_metadata)))

BL.dds <- DESeq2::DESeqDataSetFromMatrix(countData = BL.countdat,
                                            colData = BL_metadata,
                                            design = ~ 1)

BL.dds <- DESeq2::DESeq(BL.dds, parallel=T)
BL.tform <- DESeq2::varianceStabilizingTransformation(BL.dds, blind=FALSE)

## PCA
BL.herv.pca.obj <-
  pca_standard(tform = BL.tform, 
               metadata = BL_metadata, 
               var = 0.7)
# 4 PCs for Elbow method
# 15 PCs for Horn method
# 16 PCs needed to explain 50 percent of variation

################################## CLUSTERING ##################################

maxK <- 9
tDat <- assay(BL.tform)[BL.herv.pca.obj$xvars, ]
cDat <- sweep(tDat, 1, apply(tDat, 1, median, na.rm=T))

ccp.obj <- ConsensusClusterPlus(cDat, maxK=maxK, reps=1000, pItem=0.8, pFeature=0.8,
                                title="plots/bl_icl", clusterAlg="km", distance='euclidean',
                                seed=12345, plot="pdf")
icl <- calcICL(ccp.obj, title="plots/bl_icl", plot='pdf')

stopifnot(all(rownames(BL.herv.pca.obj$metadata) == names(ccp.obj[[2]]$consensusClass)))

clust.df <- lapply(2:maxK, function(k) {
  factor(paste0('C', ccp.obj[[k]]$consensusClass), levels=paste0('C', 1:k))
})%>% bind_cols() %>% data.frame
colnames(clust.df) <- paste0('clust.retro.k', 2:maxK)
rownames(clust.df) <- names(ccp.obj[[2]]$consensusClass)

BL.herv.pca.obj$metadata <- cbind(
  BL.herv.pca.obj$metadata, clust.df)

wrap <- function(nclu) {
  columnname <- paste0('clust.retro.k', nclu)
  colkey <- ccp.obj[[nclu]]$clrs[[3]]
  names(colkey) <- paste0('C', 1:nclu)
  
  gp <- PCAtools::biplot(BL.herv.pca.obj,
                         colby = columnname,
                         colkey = colkey,
                         ellipse = FALSE,
                         ellipseConf = 0.9,
                         shape = "clinical_variant", 
                         shapekey = c("Endemic BL" = 15, "Sporadic BL" = 8),
                         hline = 0,
                         vline = 0,
                         legendPosition = 'right',
                         lab=NULL,
                         title=paste0("Retrotranscriptome PCA (", columnname, ")")
  ) + 
    xlim(-45, 45) + 
    ylim(-35, 35) + 
    theme(aspect.ratio=1)
  print(gp)
  print(table(BL.herv.pca.obj$metadata[ ,columnname], 
              BL.herv.pca.obj$metadata$ebv_status))
  print(table(BL.herv.pca.obj$metadata[ ,columnname], 
              BL.herv.pca.obj$metadata$clinical_variant))
  print(table(BL.herv.pca.obj$metadata[ ,columnname], 
              BL.herv.pca.obj$metadata$gender))
  print(table(BL.herv.pca.obj$metadata[ ,columnname], 
              BL.herv.pca.obj$metadata$tissue_source_site))
  print(table(BL.herv.pca.obj$metadata[ ,columnname], 
              BL.herv.pca.obj$metadata$MYC_SV))
}

pdf("plots/05i-bl_cluster_biplot.pdf")
for(nclu in 2:maxK) {
  wrap(nclu)
}
dev.off()

cplots <- lapply(2:maxK, function(nclu) {
  p1 <- factoextra::fviz_cluster(
    list(data=t(cDat), cluster=ccp.obj[[nclu]]$consensusClass),
    # ellipse.type="norm", 
    # ellipse.level=0.9,
    geom="point",
    palette = "npg"
  ) + theme_cowplot()
  sil <- cluster::silhouette(ccp.obj[[nclu]]$consensusClass, dist(t(cDat), method = "euclidean"))
  p2 <- factoextra::fviz_silhouette(
    cluster::silhouette(ccp.obj[[nclu]]$consensusClass, dist(t(cDat), method = "euclidean")),
    palette = "npg"
  )
  return(list(p1,p2))
})
cplots <- unlist(cplots, recursive = FALSE)

pdf("plots/05i-bl_clusters_analysis.pdf", height=11.5, width=11.5)
gridExtra::marrangeGrob(grobs = cplots, layout_matrix=matrix(1:8,4,2,T))
dev.off()

sil.score <- sapply(2:maxK, function(k) {
  mean(cluster::silhouette(ccp.obj[[k]]$consensusClass, dist(t(cDat), method = "euclidean"))[,3])
})
ggplot(data.frame(x=2:maxK, y=sil.score)) + geom_point(aes(x,y)) + geom_line(aes(x,y)) + 
  xlab("k (num clusters)") + ylab('average silhouette score') + theme_cowplot()
ggsave("plots/05i_bl_clusters_sil_score.pdf")


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

save(BL.herv.pca.obj, ccp.obj, tDat, cDat, clust.df, BL_metadata,
     file="r_outputs/05i-BL_pca_ccp_clusters_metadata.Rdata")

