################################################################################
################################################################################
################################################################################
################################################################################
################################ LYMPHOMA PCAs #################################

## Plan:
## DLBCL individual DESEQ2 and PCA
## FL individual DESEQ2 and PCA
## BL individual DESEQ2 and PCA
## All together DESEQ2 and PCA (Maybe not?)

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
library(sva)

################################### LOAD DATA ##################################

load("r_outputs/02-DLBCL_filt_counts.Rdata")
load("r_outputs/02-BL_filt_counts.Rdata")
load("r_outputs/02-FL_filt_counts.Rdata")
load("r_outputs/01-metadata.Rdata")
load("r_outputs/01-refs.Rdata")
load("r_outputs/02-all_lymphoma_filt_counts.Rdata")

################################# FUNCTION SCREE PLOT #################################

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
  
  print(screeplot)
  
  return(pca.obj)
}

############################## DLBCL DESEQ HERVs ###############################

### DESeq2 (HERVs Only)

DLBCL_metadata$COO_class <- DLBCL_metadata$COO_class %>% 
  replace(is.na(.), "Unclass")


DLBCL.countdat <- DLBCL.filt.herv
cat(sprintf('%d variables\n', nrow(DLBCL.countdat)))

stopifnot(all(colnames(DLBCL.countdat) == rownames(DLBCL_metadata)))

DLBCL.dds <- DESeq2::DESeqDataSetFromMatrix(countData = DLBCL.countdat,
                                      colData = DLBCL_metadata,
                                      design = ~ project + COO_class)

DLBCL.dds <- DESeq2::DESeq(DLBCL.dds, parallel=T)
DLBCL.tform <- DESeq2::varianceStabilizingTransformation(DLBCL.dds, blind=FALSE)

## PCA
DLBCL.herv.pca.obj <-
  pca_standard(tform = DLBCL.tform, 
               metadata = DLBCL_metadata, 
               var = 0.9)
# 4 PCs for Elbow method
# 22 PCs for Horn method
# 22 PCs needed to explain 50 percent of variation

############################# DLBCL BIPLOTS HERVs ###############################

## Biplot with projects (only HERVs)

biplot(DLBCL.herv.pca.obj, 
       lab = NULL,
       showLoadings = TRUE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "project",
       colkey = c("NCICCR-DLBCL" = "orange", 
                  "TCGA-DLBC" = "lightblue"),
       legendPosition = "right")  +
  theme_cowplot()

ggsave("plots/05b-dlbcl_hervs_biplot_pc1_pc2_project.pdf", height = 6, width = 8)

# Biplot with COO call (only HERVs)

biplot(DLBCL.herv.pca.obj, 
       lab = NULL,
       showLoadings = TRUE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       ellipse = TRUE,
       colby = "COO_class",
       shape = "project", shapekey = c("NCICCR-DLBCL" = 15, "TCGA-DLBC" = 8),
       colkey = c("GCB" = "royalblue", 
                  "ABC" = "red3", 
                  "Unclass" = "lightblue", 
                  "Missing" = "grey"),
       legendPosition = "right") +
  theme_cowplot()

ggsave("plots/05b-dlbcl_hervs_biplot_pc1_pc2_coo.pdf", height = 6, width = 8)

############################ DLBCL LOADINGS HERVs ##############################

rangeRetain <- 0.01
PCAtools::plotloadings(DLBCL.herv.pca.obj,
                       title=paste0("DLBCL HERV loadings"),
                       rangeRetain = rangeRetain,
                       caption = paste0('Top ', rangeRetain * 100, '% variables'),
                       subtitle = 'PC1-PC5',
                       shapeSizeRange = c(3,3),    
                       labSize = 3.0,    
                       shape = 24,
                       col = c('white', 'pink'),
                       drawConnectors = TRUE
)

ggsave("plots/05b-dlbcl_hervs_loadings_plot.pdf", height = 10, width = 10)

########################## DLBCL DESEQ HERVS & GENES ###########################

DLBCL.g.countDat <- DLBCL.filt.comb
cat(sprintf('%d variables\n', nrow(DLBCL.g.countDat)))

stopifnot(all(colnames(DLBCL.g.countDat) == rownames(DLBCL_metadata)))

DLBCL.g.dds <- DESeq2::DESeqDataSetFromMatrix(countData = DLBCL.g.countDat,
                                      colData = DLBCL_metadata,
                                      design = ~ project + COO_class)

DLBCL.g.dds <- DESeq2::DESeq(DLBCL.g.dds, parallel=T)
DLBCL.g.tform <- DESeq2::varianceStabilizingTransformation(DLBCL.g.dds, 
                                                           blind=FALSE)

## scree plot
DLBCL.g.pca.obj <-
  pca_standard(tform = DLBCL.g.tform, 
               metadata = DLBCL_metadata, 
               var = 0.9)

# 3 PCs for Elbow method
# 45 PCs for Horn method
# 21 PCs needed to explain 50 percent of variation

all(rownames(DLBCL.g.pca.obj$loadings) %in% rownames(gene_table))
rownames(DLBCL.g.pca.obj$loadings) <- 
  gene_table[rownames(DLBCL.g.pca.obj$loadings), 'display']


######################## DLBCL BIPLOTS HERVs & GENES ###########################

## Biplot with projects (HERVs and genes)

biplot(DLBCL.g.pca.obj, 
       lab = NULL,
       showLoadings = TRUE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "project",
       colkey = c("NCICCR-DLBCL" = "orange", 
                  "TCGA-DLBC" = "lightblue"),
       legendPosition = "right")  +
  theme_cowplot()

ggsave("plots/05b-dlbcl_hervsgenes_biplot_pc1_pc2_project.pdf", height = 6, width = 8)

# Biplot with COO call (only HERVs)

biplot(DLBCL.g.pca.obj, 
       lab = NULL,
       showLoadings = TRUE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       ellipse = TRUE,
       drawConnectors = TRUE,
       colby = "COO_class",
       shape = "project", shapekey = c("NCICCR-DLBCL" = 15, "TCGA-DLBC" = 8),
       colkey = c("GCB" = "royalblue", 
                  "ABC" = "red3", 
                  "Unclass" = "lightblue", 
                  "Missing" = "grey"),
       legendPosition = "right") +
  theme_cowplot()

ggsave("plots/05b-dlbcl_hervsgenes_biplot_pc1_pc2_coo.pdf", height = 6, width = 8)

######################## DLBCL LOADINGS HERVs & GENES ##########################

rangeRetain <- 0.01
PCAtools::plotloadings(DLBCL.g.pca.obj,
                       title=paste0("DLBCL HERV & Gene loadings"),
                       rangeRetain = rangeRetain,
                       caption = paste0('Top ', rangeRetain * 100, '% variables'),
                       subtitle = 'PC1-PC5',
                       shapeSizeRange = c(3,3),    
                       labSize = 3.0,    
                       shape = 24,
                       col = c('white', 'pink'),
                       drawConnectors = TRUE
)

ggsave("plots/05b-dlbcl_hervsgenes_loadings_plot.pdf", height = 10, width = 10)



################################ BL DESEQ HERVs ################################

### DESeq2 (HERVs Only)

BL.countdat <- BL.filt.herv
cat(sprintf('%d variables\n', nrow(BL.countdat)))

stopifnot(all(colnames(BL.countdat) == rownames(BL_metadata)))

BL.dds <- DESeq2::DESeqDataSetFromMatrix(countData = BL.countdat,
                                            colData = BL_metadata,
                                            design = ~ clinical_variant +
                                              ebv_status)

BL.dds <- DESeq2::DESeq(BL.dds, parallel=T)
BL.tform <- DESeq2::varianceStabilizingTransformation(BL.dds, blind=FALSE)

## PCA
BL.herv.pca.obj <-
  pca_standard(tform = BL.tform, 
               metadata = BL_metadata, 
               var = 0.1)

# var = 0.1, 4, 19, 20
# var = 0.3, 4, 18, 19
# var = 0.5, 4, 18, 18
# var = 0.7, 4, 15, 16
# var = 0.8, 5, 15, 15
# Final var 0.7:
# 4 PCs for Elbow method
# 15 PCs for Horn method
# 16 PCs needed to explain 50 percent of variation


############################## BL BIPLOTS HERVs ################################

## Biplot with projects (only HERVs)

biplot(BL.herv.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       ellipse = FALSE,
       colby = "clinical_variant",
       colkey = c("Endemic BL" = wes_palette("Zissou1")[2], 
                  "Sporadic BL" = wes_palette("Zissou1")[4]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("plots/05b-BL_hervs_biplot_pc1_pc2_clinvariant.pdf", height = 6, width = 6)

# Biplot with COO call (only HERVs)

biplot(BL.herv.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "ebv_status",
       ellipse = FALSE,
       shape = "clinical_variant", 
       shapekey = c("Endemic BL" = 15, "Sporadic BL" = 8),
       colkey = c("EBV-positive" = wes_palette("Darjeeling1")[2], 
                  "EBV-negative" = wes_palette("Darjeeling1")[4]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("plots/05b-BL_hervs_biplot_pc1_pc2_ebvclinvar.pdf", height = 6, width = 6)

biplot(BL.herv.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "gender",
       ellipse = FALSE,
       shape = "ebv_status",
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("plots/05b-BL_hervs_biplot_pc1_pc2_gender.pdf", height = 6, width = 6)

eigencorplot(BL.herv.pca.obj,
             metavars = c('ebv_status','clinical_variant','subgroup','gender',
                          'tissue_source_site','tumor_biopsy','MYC_SV_partner',
                          'subtype', 'Total_N_SSM', 'anatomic_site_classification'))

pairsplot(BL.herv.pca.obj)

############################## BL LOADINGS HERVs ###############################

rangeRetain <- 0.01
PCAtools::plotloadings(BL.herv.pca.obj,
                       title=paste0("BL HERV loadings"),
                       rangeRetain = rangeRetain,
                       caption = paste0('Top ', rangeRetain * 100, '% variables'),
                       subtitle = 'PC1-PC5',
                       shapeSizeRange = c(3,3),    
                       labSize = 3.0,    
                       shape = 24,
                       col = c('white', 'pink'),
                       drawConnectors = TRUE
)

ggsave("plots/05b-BL_hervs_loadings_plot.pdf", height = 10, width = 10)

############################# BL DESEQ GENES ONLY ##############################

BL.gonly.countdat <- BL.filt.tx
cat(sprintf('%d variables\n', nrow(BL.gonly.countdat)))

stopifnot(all(colnames(BL.gonly.countdat) == rownames(BL_metadata)))

BL.gonly.dds <- DESeq2::DESeqDataSetFromMatrix(countData = BL.gonly.countdat,
                                           colData = BL_metadata,
                                           design = ~ clinical_variant +
                                             ebv_status)

BL.gonly.dds <- DESeq2::DESeq(BL.gonly.dds, parallel=T)
BL.gonly.tform <- DESeq2::varianceStabilizingTransformation(BL.gonly.dds, blind=FALSE)

## PCA
BL.gonly.pca.obj <-
  pca_standard(tform = BL.gonly.tform, 
               metadata = BL_metadata, 
               var = 0.1)
# var 0.1, 6, 18, 12
# var 0.3, 6, 18, 12
# var 0.5 7, 18, 12
# var 0.7 8, 12, 12

all(rownames(BL.gonly.pca.obj$loadings) %in% rownames(gene_table))
rownames(BL.gonly.pca.obj$loadings) <- 
  gene_table[rownames(BL.gonly.pca.obj$loadings), 'display']

############################## BL BIPLOTS GENES ################################

## Biplot with projects (only genes)

biplot(BL.gonly.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "clinical_variant",
       colkey = c("Endemic BL" = wes_palette("Zissou1")[2], 
                  "Sporadic BL" = wes_palette("Zissou1")[4]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("plots/05b-BL_genes_biplot_pc1_pc2_clinvariant.pdf", height = 6, width = 6)

# Biplot with EBV status + clinical variant (only genes)

biplot(BL.gonly.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
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

ggsave("plots/05b-BL_genes_pc1_pc2_ebvclinvar.pdf", height = 6, width = 6)

# Gender

biplot(BL.gonly.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "gender",
       ellipse = FALSE,
       shape = "ebv_status",
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("plots/05b-BL_genes_biplot_pc1_pc2_gender.pdf", height = 6, width = 6)

############################## BL LOADINGS GENES ###############################

rangeRetain <- 0.01
PCAtools::plotloadings(BL.gonly.pca.obj,
                       title=paste0("BL Genes loadings"),
                       rangeRetain = rangeRetain,
                       caption = paste0('Top ', rangeRetain * 100, '% variables'),
                       subtitle = 'PC1-PC5',
                       shapeSizeRange = c(3,3),    
                       labSize = 3.0,    
                       shape = 24,
                       col = c('white', 'pink'),
                       drawConnectors = TRUE
)

eigencorplot(BL.gonly.pca.obj,
             metavars = c('ebv_status','clinical_variant','subgroup','gender',
                          'tissue_source_site','tumor_biopsy','MYC_SV_partner',
                          'subtype', 'Total_N_SSM', 'anatomic_site_classification'))


############################ BL DESEQ HERVs & GENES ############################

### DESeq2 (HERVs and Genes)

BL.g.countdat <- BL.filt.comb
cat(sprintf('%d variables\n', nrow(BL.g.countdat)))

stopifnot(all(colnames(BL.g.countdat) == rownames(BL_metadata)))

BL.g.dds <- DESeq2::DESeqDataSetFromMatrix(countData = BL.g.countdat,
                                         colData = BL_metadata,
                                         design = ~ clinical_variant +
                                           ebv_status)

BL.g.dds <- DESeq2::DESeq(BL.g.dds, parallel=T)
BL.g.tform <- DESeq2::varianceStabilizingTransformation(BL.g.dds, blind=FALSE)

## PCA
BL.g.pca.obj <-
  pca_standard(tform = BL.g.tform, 
               metadata = BL_metadata, 
               var = 0.3)
# var = 0.1, 5, 8, 14
# var = 0.3, 4, 19, 15
# var = 0.5, 5, 19, 14
# var - 0.7, 6, 18, 13
# 5 PCs for Elbow method
# 18 PCs for Horn method
# 4 PCs needed to explain 50 percent of variation

all(rownames(BL.g.pca.obj$loadings) %in% rownames(gene_table))
rownames(BL.g.pca.obj$loadings) <- 
  gene_table[rownames(BL.g.pca.obj$loadings), 'display']

########################## BL BIPLOTS HERVs + GENES ############################

## Biplot with projects (HERVS + GENES)

biplot(BL.g.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "clinical_variant",
       colkey = c("Endemic BL" = wes_palette("Zissou1")[2], 
                  "Sporadic BL" = wes_palette("Zissou1")[4]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("plots/05b-BL_hervsgenes_biplot_pc1_pc2_clinvariant.pdf", height = 6, width = 8)

# Biplot with COO call (only HERVs)

biplot(BL.g.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
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

ggsave("plots/05b-BL_hervsgenes_pc1_pc2_ebvclinvar.pdf", height = 6, width = 8)

biplot(BL.g.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "gender",
       ellipse = FALSE,
       shape = "ebv_status",
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("plots/05b-BL_genes_biplot_pc1_pc2_gender.pdf", height = 6, width = 6)


########################## BL LOADINGS HERVs & GENES ###########################

rangeRetain <- 0.01
PCAtools::plotloadings(BL.g.pca.obj,
                       title=paste0("BL HERV & Gene loadings"),
                       rangeRetain = rangeRetain,
                       caption = paste0('Top ', rangeRetain * 100, '% variables'),
                       subtitle = 'PC1-PC5',
                       shapeSizeRange = c(3,3),    
                       labSize = 3.0,    
                       shape = 24,
                       col = c('white', 'pink'),
                       drawConnectors = TRUE
)

ggsave("plots/05b-BL_hervsgenes_loadings_plot.pdf", height = 10, width = 10)

######################## BL HERVS ONLY WITHOUT Y CHROM #########################

### DESeq2 (HERVs Only)

hervs.to.keep <- intersect(rownames(BL.filt.herv), 
                           retro.annot$locus[retro.annot$chrom != "chrY"])

BL.filt.herv.no.y <- BL.filt.herv[hervs.to.keep,] 
BL.countdat <- BL.filt.herv.no.y
cat(sprintf('%d variables\n', nrow(BL.countdat)))

stopifnot(all(colnames(BL.countdat) == rownames(BL_metadata)))

BL.dds.no.y <- DESeq2::DESeqDataSetFromMatrix(countData = BL.countdat,
                                         colData = BL_metadata,
                                         design = ~ clinical_variant +
                                           ebv_status)

BL.dds.no.y <- DESeq2::DESeq(BL.dds.no.y, parallel=T)
BL.tform.no.y <- DESeq2::varianceStabilizingTransformation(BL.dds.no.y, blind=FALSE)

## PCA
BL.herv.no.y.pca.obj <-
  pca_standard(tform = BL.tform.no.y, 
               metadata = BL_metadata, 
               var = 0.5)

# 0.1, 4, 20, 21
# 0.2, 4, 20, 21
# 0.3, 4, 20, 20
# 0.5, 4, 19, 19
# 0.7, 5, 17, 18

######################### BL HERVS ONLY NO Y PCA PLOTS #########################

biplot(BL.herv.no.y.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       ellipse = FALSE,
       colby = "clinical_variant",
       colkey = c("Endemic BL" = wes_palette("Zissou1")[2], 
                  "Sporadic BL" = wes_palette("Zissou1")[4]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("plots/05b-BL_hervs_no_y_biplot_pc1_pc2_clinvariant.pdf", height = 6, width = 6)

biplot(BL.herv.no.y.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "ebv_status",
       ellipse = FALSE,
       shape = "clinical_variant", 
       shapekey = c("Endemic BL" = 15, "Sporadic BL" = 8),
       colkey = c("EBV-positive" = wes_palette("Darjeeling1")[2], 
                  "EBV-negative" = wes_palette("Darjeeling1")[4]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("plots/05b-BL_hervs_no_y_biplot_pc1_pc2_ebvclinvar.pdf", height = 6, width = 6)

biplot(BL.herv.no.y.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "gender",
       ellipse = FALSE,
       shape = "ebv_status",
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

ggsave("plots/05b-BL_hervs_no_y_biplot_pc1_pc2_gender.pdf", height = 6, width = 6)

###################### BL SEPARATED BY MALE AND FEMALE #########################

### DESeq2 (HERVs Only)

BL_metadata_y <- BL_metadata[BL_metadata$gender == "Male",]
BL_metadata_x <- BL_metadata[BL_metadata$gender == "Female",]

BL.filt.herv.y <- BL.filt.herv[rownames(BL_metadata_y)]
BL.filt.herv.x <- BL.filt.herv[rownames(BL_metadata_x)]

stopifnot(all(colnames(BL.filt.herv.y) == rownames(BL_metadata_y)))
stopifnot(all(colnames(BL.filt.herv.x) == rownames(BL_metadata_x)))

BL.dds.y <- DESeq2::DESeqDataSetFromMatrix(countData = BL.filt.herv.y,
                                         colData = BL_metadata_y,
                                         design = ~ clinical_variant +
                                           ebv_status)

BL.dds.y <- DESeq2::DESeq(BL.dds.y, parallel=T)
BL.tform.y <- DESeq2::varianceStabilizingTransformation(BL.dds.y, blind=FALSE)

## PCA
BL.herv.y.pca.obj <-
  pca_standard(tform = BL.tform.y, 
               metadata = BL_metadata_y, 
               var = 0.1)

# 0.1, 5, 13, 16
# 0.2, 5, 13, 16,
# 0.3, 5, 13, 15
# 0.5, 6, 13, 5
# 0.7, 6, 13, 14

BL.dds.x <- DESeq2::DESeqDataSetFromMatrix(countData = BL.filt.herv.x,
                                           colData = BL_metadata_x,
                                           design = ~ clinical_variant +
                                             ebv_status)

BL.dds.x <- DESeq2::DESeq(BL.dds.x, parallel=T)
BL.tform.x <- DESeq2::varianceStabilizingTransformation(BL.dds.x, blind=FALSE)

## PCA
BL.herv.x.pca.obj <-
  pca_standard(tform = BL.tform.x, 
               metadata = BL_metadata_x, 
               var = 0.1)

# 0.1, 2, 10, 12
# 0.2, 2, 12, 12
# 0.3, 2, 11, 11
# 0.5, 41, 10, 11

########################### BL HERVS ONLY Y/X BIPLOTS ##########################

biplot(BL.herv.y.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "ebv_status",
       ellipse = FALSE,
       shape = "clinical_variant", 
       shapekey = c("Endemic BL" = 15, "Sporadic BL" = 8),
       colkey = c("EBV-positive" = wes_palette("Darjeeling1")[2], 
                  "EBV-negative" = wes_palette("Darjeeling1")[4]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

eigencorplot(BL.herv.y.pca.obj,
             metavars = c('ebv_status','clinical_variant','subgroup',
                          'tissue_source_site','tumor_biopsy','MYC_SV_partner',
                          'subtype', 'Total_N_SSM', 'anatomic_site_classification'))

biplot(BL.herv.x.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "ebv_status",
       ellipse = FALSE,
       shape = "clinical_variant", 
       shapekey = c("Endemic BL" = 15, "Sporadic BL" = 8),
       colkey = c("EBV-positive" = wes_palette("Darjeeling1")[2], 
                  "EBV-negative" = wes_palette("Darjeeling1")[4]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)


eigencorplot(BL.herv.x.pca.obj,
             metavars = c('ebv_status','clinical_variant','subgroup','gender',
                          'tissue_source_site','tumor_biopsy','MYC_SV_partner',
                          'subtype', 'Total_N_SSM', 'anatomic_site_classification'))

############################# BL HERVS ONLY CHR 19 ############################# 

### DESeq2 (HERVs Only)

hervs.to.keep <- intersect(rownames(BL.filt.herv), 
                           retro.annot$locus[retro.annot$chrom == "chr19"])

BL.filt.herv.19 <- BL.filt.herv[hervs.to.keep,] 
BL.countdat <- BL.filt.herv.19
cat(sprintf('%d variables\n', nrow(BL.countdat)))

stopifnot(all(colnames(BL.countdat) == rownames(BL_metadata)))

BL.dds.19 <- DESeq2::DESeqDataSetFromMatrix(countData = BL.countdat,
                                            colData = BL_metadata,
                                            design = ~ clinical_variant +
                                              ebv_status)

BL.dds.19 <- DESeq2::DESeq(BL.dds.19, parallel=T)
BL.tform.19 <- DESeq2::varianceStabilizingTransformation(BL.dds.19, blind=FALSE)

## PCA
BL.herv.19.pca.obj <-
  pca_standard(tform = BL.tform.19, 
               metadata = BL_metadata, 
               var = 0.1)

# 0.01, 10, 11, 9
# 0.1, 10, 11, 9
# 0.3, 10, 9, 8


############################ BL HERVS ONLY 19 BIPLOTS ##########################

biplot(BL.herv.19.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 3, 
       encircle = FALSE,
       sizeLoadingsNames = 4,
       lengthLoadingsArrowsFactor = 1.5,
       drawConnectors = TRUE,
       colby = "ebv_status",
       ellipse = FALSE,
       shape = "clinical_variant", 
       shapekey = c("Endemic BL" = 15, "Sporadic BL" = 8),
       colkey = c("EBV-positive" = wes_palette("Darjeeling1")[2], 
                  "EBV-negative" = wes_palette("Darjeeling1")[4]),
       legendPosition = "right")  +
  theme_cowplot() +
  theme(aspect.ratio = 1)

eigencorplot(BL.herv.19.pca.obj,
             metavars = c('ebv_status','clinical_variant','subgroup','gender',
                          'tissue_source_site','tumor_biopsy','MYC_SV_partner',
                          'subtype', 'Total_N_SSM', 'anatomic_site_classification'))

################################ FL DESEQ HERVs ################################

### DESeq2 (HERVs Only)

FL.countdat <- FL.filt.herv
cat(sprintf('%d variables\n', nrow(FL.countdat)))

stopifnot(all(colnames(FL.countdat) == rownames(FL_metadata)))

FL.dds <- DESeq2::DESeqDataSetFromMatrix(countData = FL.countdat,
                                         colData = FL_metadata,
                                         design = ~ who_diagnosis)

FL.dds <- DESeq2::DESeq(FL.dds, parallel=T)
FL.tform <- DESeq2::varianceStabilizingTransformation(FL.dds, blind=FALSE)

## PCA
FL.herv.pca.obj <-
  pca_standard(tform = FL.tform, 
               metadata = FL_metadata, 
               var = 0.1)
# 2 PCs for Elbow method
# 4 PCs for Horn method
# 4 PCs needed to explain 50 percent of variation


############################## FL BIPLOTS HERVs ################################

## Biplot with projects (only HERVs)

biplot(FL.herv.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 4, 
       encircle = FALSE,
       sizeLoadingsNames =4,
       lengthLoadingsArrowsFactor = 3,
       drawConnectors = TRUE,
       colby = "site_progression",
       colkey = c(wes_palette("GrandBudapest1")[1:5]),
       shape = "who_diagnosis", 
       shapekey = c("FOLLICULAR GRADE 1" = 15, 
                    "FOLLICULAR GRADE 2" = 8,
                    "FOLLICULAR GRADE 3A" = 2),
       legendPosition = "right")  +
  theme_cowplot()

ggsave("plots/05b-FL_hervs_biplot_pc1_pc2_who.pdf", height = 6, width = 8)

eigencorplot(FL.herv.pca.obj,
             metavars = c('who_diagnosis','gender','stage','performance_status',
                          'LDH_ratio','extranodal_sites','IPI_score',
                          'tumor_size', 'primary_treatment'))

################################ FL DESEQ GENES ################################ 

### DESeq2 (HERVs and genes)

FL.g.countdat <- FL.filt.tx
cat(sprintf('%d variables\n', nrow(FL.g.countdat)))

stopifnot(all(colnames(FL.g.countdat) == rownames(FL_metadata)))

FL.g.dds <- DESeq2::DESeqDataSetFromMatrix(countData = FL.g.countdat,
                                         colData = FL_metadata,
                                         design = ~ who_diagnosis)

FL.g.dds <- DESeq2::DESeq(FL.g.dds, parallel=T)
FL.g.tform <- DESeq2::varianceStabilizingTransformation(FL.g.dds, blind=FALSE)

## PCA
FL.g.pca.obj <-
  pca_standard(tform = FL.g.tform, 
               metadata = FL_metadata, 
               var = 0.1)
# 6 PCs for Elbow method
# 5 PCs for Horn method
# 4 PCs needed to explain 50 percent of variation

all(rownames(FL.g.pca.obj$loadings) %in% rownames(gene_table))
rownames(FL.g.pca.obj$loadings) <- 
  gene_table[rownames(FL.g.pca.obj$loadings), 'display']

########################## FL BIPLOTS HERVs & GENES ############################

## Biplot with projects (only HERVs)

biplot(FL.g.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 3/4),
       pointSize = 4, 
       encircle = FALSE,
       sizeLoadingsNames =4,
       lengthLoadingsArrowsFactor = 3,
       drawConnectors = TRUE,
       colby = "stage",
       colkey = c(wes_palette("GrandBudapest1")[1:5]),
       shape = "who_diagnosis", 
       shapekey = c("FOLLICULAR GRADE 1" = 15, 
                    "FOLLICULAR GRADE 2" = 8,
                    "FOLLICULAR GRADE 3A" = 2),
       legendPosition = "right")  +
  theme_cowplot()

ggsave("plots/05b-FL_hervsgenes_biplot_pc1_pc2_who.pdf", height = 6, width = 8)


########################## ALL TOGETHER DESEQ HERVs ############################

### DESeq2 (HERVs Only)


all_metadata <- all_metadata %>% replace(is.na(.), "Missing")

hervs.to.keep <- intersect(rownames(all.counts.filt.herv), 
                           retro.annot$locus[retro.annot$chrom != "chrY"])

all.counts.filt.herv.y <- all.counts.filt.herv[hervs.to.keep,] 


all.countdat <- all.counts.filt.herv.y
cat(sprintf('%d variables\n', nrow(all.countdat)))

stopifnot(all(colnames(all.countdat) == rownames(all_metadata)))

all.dds <- DESeq2::DESeqDataSetFromMatrix(countData = all.countdat,
                                            colData = all_metadata,
                                            design = ~ cancer_type + 0)

all.dds <- DESeq2::DESeq(all.dds, parallel=T)
all.tform <- DESeq2::varianceStabilizingTransformation(all.dds, blind=FALSE)

## PCA
all.herv.pca.obj <-
  pca_standard(tform = all.tform, 
               metadata = all_metadata, 
               var = 0.7)
# No Y: 0.1, 9, 16, 1, With Y:
# No Y: 0.2, 8, 16, 1, With Y:
# No Y: 0.3, 8, 16, 1, With Y:
# No Y: 0.5, 8, 15, 1, With Y: 9, 16, 1
# No Y: 0.7, 8, 11, 1, With Y: 9, 12, 1

############################## ALL LYMPHOMA BIPLOTS HERVs ################################

## Biplot with projects (only HERVs)

biplot(all.herv.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 0.9),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames =4,
       lengthLoadingsArrowsFactor = 3,
       drawConnectors = TRUE,
       colby = "subtype",
       colkey = c("Endemic BL EBV-positive" = wes_palette("Darjeeling1")[2], 
                  "Sporadic BL EBV-positive" = "#006053",
                  "Sporadic BL EBV-negative" = wes_palette("Darjeeling1")[4],
                  "Endemic BL EBV-negative" = "#ae5c00",
                  "GCB" = "royalblue", 
                  "ABC" = "red3", 
                  "Unclass" = "lightblue", 
                  "Missing" = "grey",
                  "FOLLICULAR GRADE 1" = "#F1BB7B",
                  "FOLLICULAR GRADE 2" = "#FD6467",
                  "FOLLICULAR GRADE 3A" = "#5B1A18"),
       shape = "cancer_type", 
       shapekey = c("DLBCL" = 15, 
                    "BL" = 8,
                    "FL" = 2),
       legendPosition = "bottom")  +
  theme_cowplot() + 
  theme(aspect.ratio = 1)
  

ggsave("plots/05b-all_lymphoma_hervs_biplot_pc1_pc2_cancertype.pdf", height = 9, width = 9)

######################### ALL LYMPHOMA LOADINGS HERVs ##########################

rangeRetain <- 0.01
PCAtools::plotloadings(all.herv.pca.obj,
                       title=paste0("All lymphoma HERV loadings"),
                       rangeRetain = rangeRetain,
                       caption = paste0('Top ', rangeRetain * 100, '% variables'),
                       subtitle = 'PC1-PC5',
                       shapeSizeRange = c(3,3),    
                       labSize = 3.0,    
                       shape = 24,
                       col = c('white', 'pink'),
                       drawConnectors = TRUE
)

ggsave("plots/05b-all_lymphoma_hervs_loadings_plot.pdf", height = 10, width = 10)

######################## ALL TOGETHER DESEQ GENES ONLY #########################

### DESeq2 (HERVs and Genes)

all.g.countdat <- all.counts.filt.tx
cat(sprintf('%d variables\n', nrow(all.g.countdat)))

stopifnot(all(colnames(all.g.countdat) == rownames(all_metadata)))

all.g.dds <- DESeq2::DESeqDataSetFromMatrix(countData = all.g.countdat,
                                          colData = all_metadata,
                                          design = ~ cancer_type + 0)

all.g.dds <- DESeq2::DESeq(all.g.dds, parallel=T)
all.g.tform <- DESeq2::varianceStabilizingTransformation(all.g.dds, blind=FALSE)

## PCA
all.g.pca.obj <-
  pca_standard(tform = all.g.tform, 
               metadata = all_metadata, 
               var = 0.5)
# 16 PCs for Elbow method
# 43 PCs for Horn method
# 7 PCs needed to explain 50 percent of variation

all(rownames(all.g.pca.obj$loadings) %in% rownames(gene_table))
rownames(all.g.pca.obj$loadings) <- 
  gene_table[rownames(all.g.pca.obj$loadings), 'display']

####################### ALL LYMPHOMA BIPLOTS GENES ONLY ########################

## Biplot with projects (HERVs and genes)

biplot(all.g.pca.obj, 
       lab = NULL,
       showLoadings = FALSE,
       boxedLoadingsNames = TRUE,
       fillBoxedLoadings = alpha("white", 0.9),
       pointSize = 2, 
       encircle = FALSE,
       sizeLoadingsNames =4,
       lengthLoadingsArrowsFactor = 3,
       drawConnectors = TRUE,
       colby = "subtype",
       colkey = c("Endemic BL EBV-positive" = wes_palette("Darjeeling1")[2], 
                  "Sporadic BL EBV-positive" = "#006053",
                  "Sporadic BL EBV-negative" = wes_palette("Darjeeling1")[4],
                  "Endemic BL EBV-negative" = "#ae5c00",
                  "GCB" = "royalblue", 
                  "ABC" = "red3", 
                  "Unclass" = "lightblue", 
                  "Missing" = "grey",
                  "FOLLICULAR GRADE 1" = "#F1BB7B",
                  "FOLLICULAR GRADE 2" = "#FD6467",
                  "FOLLICULAR GRADE 3A" = "#5B1A18"),
       shape = "cancer_type", 
       shapekey = c("DLBCL" = 15, 
                    "BL" = 8,
                    "FL" = 2),
       legendPosition = "bottom")  +
  theme_cowplot() + 
  theme(aspect.ratio = 1)

ggsave("plots/05b-all_lymphoma_genes_biplot_pc1_pc2_cancertype.pdf", height = 9, width = 9)



################################# SAVE FILES ###################################

save(DLBCL.dds, DLBCL.tform, DLBCL.herv.pca.obj, 
     DLBCL.g.dds, DLBCL.g.tform, DLBCL.g.pca.obj, 
     DLBCL_metadata,
     file="r_outputs/05b-DLBCL_pca_dds.Rdata")

save(BL.dds, BL.tform, BL.herv.pca.obj, 
     BL.g.dds, BL.g.tform, BL.g.pca.obj, 
     BL_metadata,
     file="r_outputs/05b-BL_pca_dds.Rdata")

save(FL.dds, FL.tform, FL.herv.pca.obj, 
     FL.g.dds, FL.g.tform, FL.g.pca.obj, 
     FL_metadata,
     file="r_outputs/05b-FL_pca_dds.Rdata")

save(all.dds, all.tform, all.herv.pca.obj, 
     all.g.dds, all.g.tform, all.g.pca.obj, 
     all_metadata,
     file="r_outputs/05b-all_lymphoma_pca_dds.Rdata")

load("r_outputs/05b-all_lymphoma_pca_dds.Rdata")
load("r_outputs/05b-BL_pca_dds.Rdata")
