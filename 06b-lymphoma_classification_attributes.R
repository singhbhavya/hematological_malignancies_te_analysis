################################################################################
################################################################################
################################################################################
###############################################################################
####################### LYMPHOMA CLASSIFICATION ATTRIBUTES ######################

## Plan:
## Analysis of what the lymphomas have been classified as. 

#################################### SETUP #####################################

library(knitr)
library(tidyverse)
library(matrixStats)
library(data.table)
library(DESeq2)
library(ggsci)
library(edgeR)
library(ashr)
library(cowplot)
library(readxl)
library(survival)
library(survminer)
library(GenomicDataCommons)
library(TCGAbiolinks)
library(cowplot)
library(wesanderson)

################################### LOAD DATA ##################################

load("r_outputs/02-BL_filt_counts.Rdata")
load("r_outputs/06-BL_clusters.Rdata")
load("r_outputs/02-FL_filt_counts.Rdata")
load("r_outputs/06-FL_clusters.Rdata")

################################### BL CLASSES #################################

BL_metadata$sample <- row.names(BL_metadata)
bl_clusters <- merge(bl_clusters, BL_metadata, by="sample")

bl_cluster_metadata <- bl_clusters[,c("sample", "cluster_name", 
                                      "clinical_variant", "ebv_status")]
bl_cluster_metadata$clinical_ebv <- paste(bl_cluster_metadata$clinical_variant,
                                          bl_cluster_metadata$ebv_status,
                                          sep=" ")
bl_cluster_metadata <- bl_cluster_metadata[,c("sample", "cluster_name", 
                                              "clinical_ebv")]

# Plot classifications, filled with clinical variant and EBV status
ggplot(bl_cluster_metadata, aes(x=cluster_name, fill = clinical_ebv)) +
  geom_bar(stat="count") +
  scale_fill_manual(values=c(wes_palette("Royal1"), wes_palette("Royal2"))) +
  coord_flip() +
  theme_cowplot()

# Plot clinical variant and EBV status, filled with classifications
pdf("plots/06b-bl_classifications_bar.pdf", height=4, width=10)
ggplot(bl_cluster_metadata, aes(x=clinical_ebv, fill = cluster_name)) +
  geom_bar(stat="count") +
  scale_fill_manual(values=c(wes_palette("Royal1"), wes_palette("Royal2"))) +
  coord_flip() +
  theme_cowplot() +
  labs(fill = "scCOO Cluster") +
  ylab("Number of BL Samples") +
  xlab("BL Clinical Variant and EBV Status")
dev.off()


################################### FL CLASSES #################################


FL_metadata$sample <- row.names(FL_metadata)
fl_clusters <- merge(fl_clusters, FL_metadata, by="sample")

fl_clusters$days_to_death <- fl_clusters$days_to_death[
  match(metadata$case, NCI_DLBCL_clinical_metadata$submitter_id)]

pdf("plots/06b-bf_classifications_bar.pdf", height=4, width=10)
# Plot classifications, filled with clinical variant and EBV status
ggplot(fl_clusters, aes(x=who_diagnosis, fill = cluster_name)) +
  geom_bar(stat="count") +
  scale_fill_manual(values=c(wes_palette("Royal1")[1:2], 
                             "lightblue",
                             wes_palette("Royal1")[4],
                             wes_palette("Rushmore1")[4],
                             wes_palette("Royal2")[1],
                             wes_palette("Rushmore1"))) +
  coord_flip() +
  theme_cowplot() +
  labs(fill = "scCOO Cluster") +
  ylab("Number of FL Samples") +
  xlab("FL WHO Diagnosis")

dev.off()


