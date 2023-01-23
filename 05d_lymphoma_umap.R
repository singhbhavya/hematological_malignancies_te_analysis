################################################################################
################################################################################
################################################################################
################################################################################
################################ LYMPHOMA UMAP #################################

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
library(umap)
library(magrittr)

################################### LOAD DATA ##################################

load("r_outputs/05b-all_lymphoma_pca_dds.Rdata")
load("r_outputs/01-metadata.Rdata")

all_metadata$sample <- row.names(all_metadata)
all_metadata$subtype2 <- all_metadata$subtype
all_metadata$subtype2[all_metadata$cancer_type == "DLBCL"] <- 
  DLBCL_metadata$LymphGen_call[match(all_metadata$sample, DLBCL_metadata$case)]
all_metadata$subtype2[all_metadata$cancer_type == "BL"] <- 
  all_metadata$subtype[all_metadata$cancer_type == "BL"]
all_metadata$subtype2[all_metadata$cancer_type == "FL"] <-
  all_metadata$subtype[all_metadata$cancer_type == "FL"]
all_metadata$subtype2[is.na(all_metadata$subtype2)] <- "Unknown"

norm.counts.all <- assay(all.tform) %>%
  t() 

# Perform UMAP on the normalized data
custom <- umap::umap.defaults


custom$n_neighbors <- 9
custom$min_dist <- 0.01

umap.all <-
  umap::umap(norm.counts.all, custom)

umap.plot.all.df <- 
  data.frame(umap.all$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("sample") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(all_metadata, by = "sample")

ggplot(
  umap.plot.all.df,
  aes(
    x = X1,
    y = X2,
    color = subtype,
    shape = cancer_type)) +
  geom_point() +
  theme_cowplot() 

