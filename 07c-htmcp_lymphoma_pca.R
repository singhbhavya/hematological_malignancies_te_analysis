################################################################################
################################################################################
################################################################################
################################################################################
################################ LYMPHOMA PCAs #################################

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

load("r_outputs/07-htmcp_all_lymphoma_filt_counts.Rdata")
load("r_outputs/07-htmcp_DLBCL_filt_counts.Rdata")
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
  
  print(screeplot)
  
  return(pca.obj)
}

########################## ALL TOGETHER DESEQ HERVs ############################

### DESeq2 (HERVs Only)

all_metadata <- all_metadata_hiv %>% replace(is.na(.), "Missing")

hervs.to.keep <- intersect(rownames(all.hiv.filt.herv), 
                           retro.annot$locus[retro.annot$chrom != "chrY"])
all.counts.filt.herv.y <- all.hiv.filt.herv[hervs.to.keep,] 

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
                  "FOLLICULAR GRADE 3A" = "#5B1A18",
                  "HIV-positive" = "#fa7de7"),
       shape = "cancer_type", 
       shapekey = c("DLBCL" = 15, 
                    "BL" = 8,
                    "FL" = 2),
       legendPosition = "bottom")  +
  theme_cowplot() + 
  theme(aspect.ratio = 1)

ggsave("plots/07c-htmcp_all_lymphoma_hervs_biplot_pc1_pc2_cancertype.pdf", height = 9, width = 9)
