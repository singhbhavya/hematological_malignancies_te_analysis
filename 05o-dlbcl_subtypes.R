################################################################################
################################################################################
################################################################################
################################################################################
########################## DLBCL SUBTYPES AND FEATURES #########################

## Plan:
##

#################################### SETUP #####################################

library(tidyverse)
library(DESeq2)
library(matrixStats)
library(glmnet)
library(Boruta)
library(UpSetR)
library(c060)
library(randomForest)
library(rpart)
library(rpart.plot)
library(wesanderson)
library(ConsensusClusterPlus)
library(cowplot)
library(PCAtools)

################################### LOAD DATA ##################################

load("r_outputs/02-DLBCL_filt_counts.Rdata")
load("r_outputs/01-refs.Rdata")
remove(DLBCL.counts.mfilt.comb,DLBCL.counts.mfilt.tx, DLBCL.counts.mfilt.rtx,
       DLBCL.counts.mfilt.herv, DLBCL.counts.mfilt.l1)
retro.annot.v2 <- read.csv("/efs/projects/hematological_malignancies_te_analysis/refs/TE_annotation.v2.0.tsv",
                           sep = "\t")
rownames(retro.annot.v2) <- retro.annot.v2$Locus

################################### LOAD DATA ##################################

intergenic.exonic <- intersect(rownames(DLBCL.filt.herv), 
                        rownames(retro.annot.v2[retro.annot.v2$TE_type == "INTERGENIC" |
                                                  retro.annot.v2$TE_type == "EXONIC",]))

intergenic <- intersect(rownames(DLBCL.filt.herv), 
                        rownames(retro.annot.v2[retro.annot.v2$TE_type == "INTERGENIC",]))

DLBCL.filt.int.ex <- DLBCL.filt.herv[intergenic.exonic,]
DLBCL.filt.int <- DLBCL.filt.herv[intergenic,]

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

############################## DLBCL DESEQ HERVs ###############################

### DESeq2 (HERVs Only)

DLBCL_metadata$COO_class <- DLBCL_metadata$COO_class %>% 
  replace(is.na(.), "Unclass")

DLBCL.countdat <- DLBCL.filt.herv
cat(sprintf('%d variables\n', nrow(DLBCL.countdat)))

stopifnot(all(colnames(DLBCL.countdat) == rownames(DLBCL_metadata)))

DLBCL.dds <- DESeq2::DESeqDataSetFromMatrix(countData = DLBCL.countdat,
                                            colData = DLBCL_metadata,
                                            design = ~ project)

DLBCL.dds <- DESeq2::DESeq(DLBCL.dds, parallel=T)
DLBCL.tform <- DESeq2::varianceStabilizingTransformation(DLBCL.dds, blind=FALSE)

## PCA
DLBCL.herv.pca.obj <-
  pca_standard(tform = DLBCL.tform, 
               metadata = DLBCL_metadata, 
               var = 0.9)
# 4 PCs for Elbow method
# 25 PCs for Horn method
# 22 PCs needed to explain 50 percent of variation

############################## DLBCL DESEQ GENES ###############################


### DESeq2 (Genes Only)

DLBCL_metadata$COO_class <- DLBCL_metadata$COO_class %>% 
  replace(is.na(.), "Unclass")

DLBCL.countdat <- DLBCL.filt.tx
cat(sprintf('%d variables\n', nrow(DLBCL.countdat)))

stopifnot(all(colnames(DLBCL.countdat) == rownames(DLBCL_metadata)))

DLBCL.g.dds <- DESeq2::DESeqDataSetFromMatrix(countData = DLBCL.countdat,
                                            colData = DLBCL_metadata,
                                            design = ~ project)

DLBCL.g.dds <- DESeq2::DESeq(DLBCL.g.dds, parallel=T)
DLBCL.g.tform <- DESeq2::varianceStabilizingTransformation(DLBCL.g.dds, blind=FALSE)

## PCA
DLBCL.g.pca.obj <-
  pca_standard(tform = DLBCL.g.tform, 
               metadata = DLBCL_metadata, 
               var = 0.9)
# 3 PCs for Elbow method
# 44 PCs for Horn method
# 20 PCs needed to explain 50 percent of variation

############################### DLBCL BIPLOTS ##################################

############################ Color keys and function ########################### 

colkey_coo <- c("GCB" = "royalblue", 
                "ABC" = "red3", 
                "Unclass" = "lightblue", 
                "Missing" = "grey")

colkey_ecotyper <- c("S1" = "#bcc779",
                     "S2" = "#60aaec",
                     "S3" = "#Ecaf60",
                     "S4" = "#B789EE",
                     "S5" = "#E0AF71",
                     "Unassigned" = "grey",
                     "Missing" = "grey")

colkey_lymphgen <- c("A53" = "#e78bf0",
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
                     "NA" = "grey")

dlbcl_biplots <- function(pca.obj, col_by, colkey, ellipse_arg) {
  biplot(pca.obj, 
         lab = NULL,
         showLoadings = FALSE,
         pointSize = 2, 
         ellipse = ellipse_arg,
         colby = col_by,
         shape = "project", 
         shapekey = c("NCICCR-DLBCL" = 15, "TCGA-DLBC" = 8),
         colkey = colkey,
         legendPosition = "right") +
    theme_cowplot() 
}

################################### Lymphgen ###################################

lymphgen.herv <-
  dlbcl_biplots(DLBCL.herv.pca.obj, "LymphGen_call", colkey_lymphgen, FALSE)

lymphgen.g <-
  dlbcl_biplots(DLBCL.g.pca.obj, "LymphGen_call", colkey_lymphgen, FALSE) + 
  scale_x_reverse() + scale_y_reverse()

lymphgen.l <- get_legend(
  lymphgen.g + 
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "right")
)


pdf("plots/05o_dlbcl_lymphgen.pdf", height=5.5, width=10)
plot_grid(plot_grid(lymphgen.herv + theme(legend.position = c("none")),
                    lymphgen.g + theme(legend.position = c("none"))), 
          lymphgen.l,
          rel_widths = c(3, .6))
dev.off()

remove(lymphgen.herv, lymphgen.g, lymphgen.l)

##################################### COO ######################################

coo.herv <-
  dlbcl_biplots(DLBCL.herv.pca.obj, "COO_class", colkey_coo, FALSE)

coo.g <-
  dlbcl_biplots(DLBCL.g.pca.obj, "COO_class", colkey_coo, FALSE) + 
  scale_x_reverse() + scale_y_reverse()

coo.l <- get_legend(
  coo.g + 
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "right")
)

pdf("plots/05o_dlbcl_coo.pdf", height=5, width=10)
plot_grid(plot_grid(coo.herv + theme(legend.position = c("none")),
                    coo.g + theme(legend.position = c("none"))), 
          coo.l,
          rel_widths = c(3, .4))
dev.off()

remove(coo.herv, coo.g, coo.l)

################################### EcoTyper ################################### 

eco.herv <-
  dlbcl_biplots(DLBCL.herv.pca.obj, "EcoTyper_call", colkey_ecotyper, FALSE)

eco.g <-
  dlbcl_biplots(DLBCL.g.pca.obj, "EcoTyper_call", colkey_ecotyper, FALSE) + 
  scale_x_reverse() + scale_y_reverse()

eco.l <- get_legend(
  eco.g + 
    guides(color = guide_legend(ncol = 1)) +
    theme(legend.position = "right")
)

pdf("plots/05o_dlbcl_ecotyper.pdf", height=5, width=10)
plot_grid(plot_grid(eco.herv + theme(legend.position = c("none")),
                    eco.g + theme(legend.position = c("none"))), 
          eco.l,
          rel_widths = c(3, .5))
dev.off()

remove(eco.herv, eco.g, eco.l)

# Remove supervised dds and pca objects
remove(DLBCL.herv.pca.obj, DLBCL.g.pca.obj, DLBCL.dds, DLBCL.g.dds)

################################# UNSUP CLUST ################################## 

### DESeq2 (HERVs Only)

DLBCL.countdat <- DLBCL.filt.herv
cat(sprintf('%d variables\n', nrow(DLBCL.countdat)))
stopifnot(all(colnames(DLBCL.countdat) == rownames(DLBCL_metadata)))

DLBCL.dds <- DESeq2::DESeqDataSetFromMatrix(countData = DLBCL.countdat,
                                            colData = DLBCL_metadata,
                                            design = ~ 1)

DLBCL.dds <- DESeq2::DESeq(DLBCL.dds, parallel=T)
DLBCL.tform <- DESeq2::varianceStabilizingTransformation(DLBCL.dds, blind=FALSE)

## PCA
DLBCL.herv.pca.obj <-
  pca_standard(tform = DLBCL.tform, 
               metadata = DLBCL_metadata, 
               var = 0.9)


################################# UNSUP CLUST ################################## 

maxK <- 9
tDat <- assay(DLBCL.tform)[DLBCL.herv.pca.obj$xvars, ]
cDat <- sweep(tDat, 1, apply(tDat, 1, median, na.rm=T))

ccp.obj <- ConsensusClusterPlus(cDat, maxK=maxK, reps=1000, pItem=0.8, pFeature=0.8,
                                title="plots/dlbcl_icl", clusterAlg="km", distance='euclidean',
                                seed=12345, plot="pdf")
icl <- calcICL(ccp.obj, title="plots/dlbcl_icl", plot='pdf')

stopifnot(all(rownames(DLBCL.herv.pca.obj$metadata) == names(ccp.obj[[2]]$consensusClass)))

clust.df <- lapply(2:maxK, function(k) {
  factor(paste0('C', ccp.obj[[k]]$consensusClass), levels=paste0('C', 1:k))
})%>% bind_cols() %>% data.frame
colnames(clust.df) <- paste0('clust.retro.k', 2:maxK)
rownames(clust.df) <- names(ccp.obj[[2]]$consensusClass)

DLBCL.herv.pca.obj$metadata <- cbind(
  DLBCL.herv.pca.obj$metadata, clust.df)

wrap <- function(nclu) {
  columnname <- paste0('clust.retro.k', nclu)
  colkey <- ccp.obj[[nclu]]$clrs[[3]]
  names(colkey) <- paste0('C', 1:nclu)
  
  gp <- PCAtools::biplot(DLBCL.herv.pca.obj,
                         colby = columnname,
                         colkey = colkey,
                         ellipse = FALSE,
                         ellipseConf = 0.9,
                         shape = "project", 
                         shapekey = c("NCICCR-DLBCL" = 15, "TCGA-DLBC" = 8),
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
  print(table(DLBCL.herv.pca.obj$metadata[ ,columnname], 
              DLBCL.herv.pca.obj$metadata$COO_class))
  print(table(DLBCL.herv.pca.obj$metadata[ ,columnname], 
              DLBCL.herv.pca.obj$metadata$LymphGen_call))
  print(table(DLBCL.herv.pca.obj$metadata[ ,columnname], 
              DLBCL.herv.pca.obj$metadata$EcoTyper_call))
  print(table(DLBCL.herv.pca.obj$metadata[ ,columnname], 
              DLBCL.herv.pca.obj$metadata$DblHit_call))
  print(table(DLBCL.herv.pca.obj$metadata[ ,columnname], 
              DLBCL.herv.pca.obj$metadata$RCHOP_like_chemo))
}

pdf("plots/05o-dlbcl_cluster_biplot.pdf")
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

pdf("plots/05o-dlbcl_clusters_analysis.pdf", height=11.5, width=11.5)
gridExtra::marrangeGrob(grobs = cplots, layout_matrix=matrix(1:8,4,2,T))
dev.off()

sil.score <- sapply(2:maxK, function(k) {
  mean(cluster::silhouette(ccp.obj[[k]]$consensusClass, dist(t(cDat), method = "euclidean"))[,3])
})
ggplot(data.frame(x=2:maxK, y=sil.score)) + geom_point(aes(x,y)) + geom_line(aes(x,y)) + 
  xlab("k (num clusters)") + ylab('average silhouette score') + theme_cowplot()
ggsave("plots/05o_dlbcl_clusters_sil_score.pdf")

pdf("plots/05o_dlbcl_clusters_eigen.pdf", height=10, width=10)
PCAtools::eigencorplot(DLBCL.herv.pca.obj,
                       metavars = c('project', 'gender', 'COO_class', 'LymphGen_call',
                                    'EcoTyper_call', 
                                    'clust.retro.k2','clust.retro.k3','clust.retro.k4','clust.retro.k5',
                                    'clust.retro.k6', 'clust.retro.k7', 'clust.retro.k8',
                                    'clust.retro.k9'),
                       col=viridis::viridis_pal()(12),
                       signifSymbols = c('****', '***', '**', '*', ''),
                       signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
)
dev.off()

sink(file = "r_outputs/05o-dlbcl_clusters.txt")
print(table(DLBCL.herv.pca.obj$metadata$clust.retro.k2, DLBCL.herv.pca.obj$metadata$COO_class))
print(table(DLBCL.herv.pca.obj$metadata$clust.retro.k3, DLBCL.herv.pca.obj$metadata$COO_class))
print(table(DLBCL.herv.pca.obj$metadata$clust.retro.k4, DLBCL.herv.pca.obj$metadata$COO_class))
print(table(DLBCL.herv.pca.obj$metadata$clust.retro.k5, DLBCL.herv.pca.obj$metadata$COO_class))
print(table(DLBCL.herv.pca.obj$metadata$clust.retro.k6, DLBCL.herv.pca.obj$metadata$COO_class))
print(table(DLBCL.herv.pca.obj$metadata$clust.retro.k7, DLBCL.herv.pca.obj$metadata$COO_class))
print(table(DLBCL.herv.pca.obj$metadata$clust.retro.k8, DLBCL.herv.pca.obj$metadata$COO_class))
print(table(DLBCL.herv.pca.obj$metadata$clust.retro.k9, DLBCL.herv.pca.obj$metadata$COO_class))
sink(file = NULL)


# Add clusters to metadata
DLBCL_metadata$clust.retro.k2 <- clust.df$clust.retro.k2
DLBCL_metadata$clust.retro.k3 <- clust.df$clust.retro.k3
DLBCL_metadata$clust.retro.k4 <- clust.df$clust.retro.k4
DLBCL_metadata$clust.retro.k5 <- clust.df$clust.retro.k5
DLBCL_metadata$clust.retro.k7 <- clust.df$clust.retro.k7
DLBCL_metadata$clust.retro.k9 <- clust.df$clust.retro.k9

# reformat
clust.df$clust.retro.k7 <- paste0("H", paste0("H", DLBCL_metadata$clust.retro.k7))
DLBCL_metadata$clust.retro.k7 <- paste0("H", DLBCL_metadata$clust.retro.k7)
DLBCL.herv.pca.obj$metadata$clust.retro.k7 <- paste0("H", 
                                                     DLBCL.herv.pca.obj$metadata$clust.retro.k7)



################################## SAVE FILES ##################################

save(DLBCL.herv.pca.obj, ccp.obj, tDat, cDat, clust.df, DLBCL_metadata,
     file="r_outputs/05o-DLBCL_pca_ccp_clusters_metadata.Rdata")

################################################################################
################################################################################
############################## #INTERGENIC HERVS ###############################
################################################################################
################################################################################


############################## DLBCL DESEQ HERVs ###############################

### DESeq2 (HERVs Only)

DLBCL.countdat <- DLBCL.filt.int.ex
cat(sprintf('%d variables\n', nrow(DLBCL.countdat)))

stopifnot(all(colnames(DLBCL.countdat) == rownames(DLBCL_metadata)))

DLBCL.int.ex.dds <- DESeq2::DESeqDataSetFromMatrix(countData = DLBCL.countdat,
                                            colData = DLBCL_metadata,
                                            design = ~ project)

DLBCL.int.ex.dds <- DESeq2::DESeq(DLBCL.int.ex.dds, parallel=T)
DLBCL.int.ex.tform <- DESeq2::varianceStabilizingTransformation(DLBCL.int.ex.dds, blind=FALSE)

## PCA
DLBCL.int.ex.pca.obj <-
  pca_standard(tform = DLBCL.int.ex.tform, 
               metadata = DLBCL_metadata, 
               var = 0.9)
# 4 PCs for Elbow method
# 21 PCs for Horn method
# 20 PCs needed to explain 50 percent of variation


################################### BIPLOTS ####################################

dlbcl_biplots(DLBCL.int.ex.pca.obj, "COO_class", colkey_coo, FALSE)
dlbcl_biplots(DLBCL.herv.pca.obj, "LymphGen_call", colkey_lymphgen, FALSE)
dlbcl_biplots(DLBCL.herv.pca.obj, "EcoTyper_call", colkey_ecotyper, FALSE)

################################# UNSUP CLUST ################################## 

maxK <- 9
tDat.int.ex <- assay(DLBCL.int.ex.tform)[DLBCL.int.ex.pca.obj$xvars, ]
cDat.int.ex <- sweep(tDat.int.ex, 1, apply(tDat.int.ex, 1, median, na.rm=T))

ccp.int.ex.obj <- ConsensusClusterPlus(cDat.int.ex, 
                                       maxK=maxK, reps=1000, pItem=0.8, pFeature=0.8,
                                title="plots/dlbcl_icl.int.ex", clusterAlg="km", distance='euclidean',
                                seed=12345, plot="pdf")
icl <- calcICL(ccp.int.ex.obj, title="plots/dlbcl_icl.int.ex", plot='pdf')

stopifnot(all(rownames(DLBCL.int.ex.pca.obj$metadata) ==
                names(ccp.int.ex.obj[[2]]$consensusClass)))

clust.int.ex.df <- lapply(2:maxK, function(k) {
  factor(paste0('C', ccp.int.ex.obj[[k]]$consensusClass), levels=paste0('C', 1:k))
})%>% bind_cols() %>% data.frame
colnames(clust.int.ex.df) <- paste0('clust.retro.k', 2:maxK)
rownames(clust.int.ex.df) <- names(ccp.int.ex.obj[[2]]$consensusClass)

DLBCL.int.ex.pca.obj$metadata <- cbind(
  DLBCL.int.ex.pca.obj$metadata, clust.int.ex.df)

wrap <- function(nclu) {
  columnname <- paste0('clust.retro.k', nclu)
  colkey <- ccp.int.ex.obj[[nclu]]$clrs[[3]]
  names(colkey) <- paste0('C', 1:nclu)
  
  gp <- PCAtools::biplot(DLBCL.int.ex.pca.obj,
                         colby = columnname,
                         colkey = colkey,
                         ellipse = FALSE,
                         ellipseConf = 0.9,
                         shape = "project", 
                         shapekey = c("NCICCR-DLBCL" = 15, "TCGA-DLBC" = 8),
                         hline = 0,
                         vline = 0,
                         legendPosition = 'right',
                         lab=NULL,
                         title=paste0("Intergenic/Exonic HERV PCA (", columnname, ")")
  ) + 
    xlim(-45, 45) + 
    ylim(-35, 35) + 
    theme(aspect.ratio=1)
  print(gp)
  print(table(DLBCL.int.ex.pca.obj$metadata[ ,columnname], 
              DLBCL.int.ex.pca.obj$metadata$COO_class))
  print(table(DLBCL.int.ex.pca.obj$metadata[ ,columnname], 
              DLBCL.int.ex.pca.obj$metadata$LymphGen_call))
  print(table(DLBCL.int.ex.pca.obj$metadata[ ,columnname], 
              DLBCL.int.ex.pca.obj$metadata$EcoTyper_call))
  print(table(DLBCL.int.ex.pca.obj$metadata[ ,columnname], 
              DLBCL.int.ex.pca.obj$metadata$DblHit_call))
  print(table(DLBCL.int.ex.pca.obj$metadata[ ,columnname], 
              DLBCL.int.ex.pca.obj$metadata$RCHOP_like_chemo))
}

pdf("plots/05o-dlbcl_int_ex_cluster_biplot.pdf")
for(nclu in 2:maxK) {
  wrap(nclu)
}
dev.off()

cplots <- lapply(2:maxK, function(nclu) {
  p1 <- factoextra::fviz_cluster(
    list(data=t(cDat.int.ex), 
         cluster=ccp.int.ex.obj[[nclu]]$consensusClass),
    # ellipse.type="norm", 
    # ellipse.level=0.9,
    geom="point",
    palette = "npg"
  ) + theme_cowplot()
  sil <- cluster::silhouette(ccp.int.ex.obj[[nclu]]$consensusClass, 
                             dist(t(cDat.int.ex), 
                                  method = "euclidean"))
  p2 <- factoextra::fviz_silhouette(
    cluster::silhouette(ccp.int.ex.obj[[nclu]]$consensusClass, 
                        dist(t(cDat.int.ex), 
                             method = "euclidean")),
    palette = "npg"
  )
  return(list(p1,p2))
})
cplots <- unlist(cplots, recursive = FALSE)

pdf("plots/05o-dlbcl_int_ex_clusters_analysis.pdf", height=11.5, width=11.5)
gridExtra::marrangeGrob(grobs = cplots, layout_matrix=matrix(1:8,4,2,T))
dev.off()

sil.score <- sapply(2:maxK, function(k) {
  mean(cluster::silhouette(ccp.int.ex.obj[[k]]$consensusClass, 
                           dist(t(cDat.int.ex), method = "euclidean"))[,3])
})

ggplot(data.frame(x=2:maxK, y=sil.score)) + 
  geom_point(aes(x,y)) + 
  geom_line(aes(x,y)) + 
  xlab("k (num clusters)") + 
  ylab('average silhouette score') + theme_cowplot()
ggsave("plots/05o_dlbcl_int_ex_clusters_sil_score.pdf")

pdf("plots/05o_dlbcl_clusters_int_ex_eigen.pdf", height=10, width=10)
PCAtools::eigencorplot(DLBCL.int.ex.pca.obj,
                       metavars = c('project', 'gender', 'COO_class', 'LymphGen_call',
                                    'EcoTyper_call', 
                                    'clust.retro.k2','clust.retro.k3','clust.retro.k4','clust.retro.k5',
                                    'clust.retro.k6', 'clust.retro.k7', 'clust.retro.k8',
                                    'clust.retro.k9'),
                       col=viridis::viridis_pal()(12),
                       signifSymbols = c('****', '***', '**', '*', ''),
                       signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
)
dev.off()

sink(file = "r_outputs/05o-dlbcl_int_ex_clusters.txt")
print(table(DLBCL.int.ex.pca.obj$metadata$clust.retro.k2, DLBCL.int.ex.pca.obj$metadata$COO_class))
print(table(DLBCL.int.ex.pca.obj$metadata$clust.retro.k3, DLBCL.int.ex.pca.obj$metadata$COO_class))
print(table(DLBCL.int.ex.pca.obj$metadata$clust.retro.k4, DLBCL.int.ex.pca.obj$metadata$COO_class))
print(table(DLBCL.int.ex.pca.obj$metadata$clust.retro.k5, DLBCL.int.ex.pca.obj$metadata$COO_class))
print(table(DLBCL.int.ex.pca.obj$metadata$clust.retro.k6, DLBCL.int.ex.pca.obj$metadata$COO_class))
print(table(DLBCL.int.ex.pca.obj$metadata$clust.retro.k7, DLBCL.int.ex.pca.obj$metadata$COO_class))
print(table(DLBCL.int.ex.pca.obj$metadata$clust.retro.k8, DLBCL.int.ex.pca.obj$metadata$COO_class))
print(table(DLBCL.int.ex.pca.obj$metadata$clust.retro.k9, DLBCL.int.ex.pca.obj$metadata$COO_class))
sink(file = NULL)




