################################################################################
################################################################################
################################################################################
################################################################################
################################ LYMPHOMA PCAs #################################

## Plan:
## DLBCL individual DESEQ2 and PCA
## FL individual DESEQ2 and PCA
## BL individual DESEQ2 and PCA
## All together DESEQ2 and PCA

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


################################### LOAD DATA ##################################

load("r_outputs/02-DLBCL_filt_counts.Rdata")
load("r_outputs/02-BL_filt_counts.Rdata")
load("r_outputs/02-FL_filt_counts.Rdata")
load("r_outputs/01-metadata.Rdata")
load("r_outputs/01-refs.Rdata")

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
  
  return(pca.obj)
}

############################## DLBCL DESEQ HERVs ###############################

### DESeq2 (HERVs Only)

DLBCL.countdat <- DLBCL.counts.mfilt.herv
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
# 8 PCs for Elbow method
# 34 PCs for Horn method
# 29 PCs needed to explain 50 percent of variation

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

DLBCL.g.countDat <- DLBCL.counts.mfilt.comb
cat(sprintf('%d variables\n', nrow(DLBCL.g.countDat)))

stopifnot(all(colnames(DLBCL.g.countDat) == rownames(DLBCL_metadata)))

DLBCL.g.dds <- DESeq2::DESeqDataSetFromMatrix(countData = DLBCL.g.countDat,
                                      colData = DLBCL_metadata,
                                      design = ~ project)

DLBCL.g.dds <- DESeq2::DESeq(DLBCL.g.dds, parallel=T)
DLBCL.g.tform <- DESeq2::varianceStabilizingTransformation(DLBCL.g.dds, 
                                                           blind=FALSE)

## scree plot
DLBCL.g.pca.obj <-
  pca_standard(tform = DLBCL.g.tform, 
               metadata = DLBCL_metadata, 
               var = 0.9)

# 5 PCs for Elbow method
# 50 PCs for Horn method
# 20 PCs needed to explain 50 percent of variation

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

BL.countdat <- BL.counts.mfilt.herv
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
# 5 PCs for Elbow method
# 11 PCs for Horn method
# 12 PCs needed to explain 50 percent of variation


############################## BL BIPLOTS HERVs ################################

## Biplot with projects (only HERVs)

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
       colby = "clinical_variant",
       colkey = c("Endemic BL" = wes_palette("Zissou1")[2], 
                  "Sporadic BL" = wes_palette("Zissou1")[4]),
       legendPosition = "right")  +
  theme_cowplot()

ggsave("plots/05b-BL_hervs_biplot_pc1_pc2_clinvariant.pdf", height = 6, width = 8)

# Biplot with COO call (only HERVs)

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
  theme_cowplot()

ggsave("plots/05b-BL_hervs_biplot_pc1_pc2_ebvclinvar.pdf", height = 6, width = 8)

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

############################ BL DESEQ HERVs & GENES ############################

### DESeq2 (HERVs and Genes)

BL.g.countdat <- BL.counts.mfilt.comb
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
               var = 0.1)
# 8 PCs for Elbow method
# 11 PCs for Horn method
# 9 PCs needed to explain 50 percent of variation

all(rownames(BL.g.pca.obj$loadings) %in% rownames(gene_table))
rownames(BL.g.pca.obj$loadings) <- 
  gene_table[rownames(BL.g.pca.obj$loadings), 'display']

############################## BL BIPLOTS HERVs ################################

## Biplot with projects (only HERVs)

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
       colby = "clinical_variant",
       colkey = c("Endemic BL" = wes_palette("Zissou1")[2], 
                  "Sporadic BL" = wes_palette("Zissou1")[4]),
       legendPosition = "right")  +
  theme_cowplot()

ggsave("plots/05b-BL_hervsgenes_biplot_pc1_pc2_clinvariant.pdf", height = 6, width = 8)

# Biplot with COO call (only HERVs)

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
  theme_cowplot()

ggsave("plots/05b-BL_hervsgenes_pc1_pc2_ebvclinvar.pdf", height = 6, width = 8)

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

################################ FL DESEQ HERVs ################################

### DESeq2 (HERVs Only)

FL.countdat <- FL.counts.mfilt.herv
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
       showLoadings = TRUE,
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

ggsave("plots/05b-FL_hervs_biplot_pc1_pc2_who.pdf", height = 6, width = 8)


############################ FL DESEQ HERVs & GENEs ############################

### DESeq2 (HERVs and genes)

FL.g.countdat <- FL.counts.mfilt.comb
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
# 4 PCs for Elbow method
# 5 PCs for Horn method
# 4 PCs needed to explain 50 percent of variation

all(rownames(FL.g.pca.obj$loadings) %in% rownames(gene_table))
rownames(FL.g.pca.obj$loadings) <- 
  gene_table[rownames(FL.g.pca.obj$loadings), 'display']

########################## FL BIPLOTS HERVs & GENES ############################

## Biplot with projects (only HERVs)

biplot(FL.g.pca.obj, 
       lab = NULL,
       showLoadings = TRUE,
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
