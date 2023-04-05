################################################################################
################################################################################
################################################################################
###############################################################################
####################### LYMPHOMA CLASSIFICATION SURVIVAL #######################

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

load("r_outputs/05o-DLBCL_pca_ccp_clusters_metadata.Rdata")
load("r_outputs/01-refs.Rdata")

##################################### SETUP ####################################


# Add clusters to metadata
DLBCL_metadata$clust.retro.k2 <- clust.df$clust.retro.k2
DLBCL_metadata$clust.retro.k3 <- clust.df$clust.retro.k3
DLBCL_metadata$clust.retro.k4 <- clust.df$clust.retro.k4
DLBCL_metadata$clust.retro.k5 <- clust.df$clust.retro.k5
DLBCL_metadata$clust.retro.k7 <- clust.df$clust.retro.k7
DLBCL_metadata$clust.retro.k9 <- clust.df$clust.retro.k9

metadata <- DLBCL_metadata

metadata$status <- ifelse(metadata$vital_status == "Alive", 1, 2)

# Download TCGA metadata
NCI_DLBCL_clinical_metadata <- 
  GDCquery_clinic("NCICCR-DLBCL", type = "clinical", save.csv = FALSE)

metadata$days_to_death <- NCI_DLBCL_clinical_metadata$days_to_death[
  match(metadata$case, NCI_DLBCL_clinical_metadata$submitter_id)]

metadata$time <- ifelse(!is.na(metadata$days_to_death), 
                        metadata$days_to_death, 
                        metadata$days_to_last_follow_up)

################################### SURVIVAL ###################################

fit <- survfit(Surv(time, status) ~ clust.retro.k7, data = metadata)
print(fit)
res.sum <- surv_summary(fit)
print(res.sum)
surv_diff <- survival::survdiff(Surv(time, status) ~ clust.retro.k7, data = metadata)
print(surv_diff)

pdf("plots/06c-dlbcl_cluster.retro_survival.pdf", height = 10, width=10)

herv_clusters <-
  ggsurvplot(fit,
             risk.table =  "abs_pct",
             pval = TRUE,
             break.time.by = 500,
             conf.int = FALSE,
             risk.table.col = "strata",
             risk.table.y.text.col = T,
             risk.table.y.text = FALSE,
             linetype = "solid",
             surv.median.line = "hv",
             ncensor.plot = FALSE,
             ggtheme = theme_bw(), # Change ggplot2 theme
             palette =  ggsci::pal_npg(palette = c("nrc"))(7),
             legend.labs = 
               c("C1", "C2", "C3", "C4", "C5", "C6", "C7"),
             font.x = c(20),
             font.y = c(20),
             font.tickslab = c(20),
             size = 2,
             tables.theme = clean_theme()) 

herv_clusters$plot <- herv_clusters$plot + 
  theme(legend.text = element_text(size = 15, color = "black", face = "bold"))

herv_clusters
dev.off()

fit <- survfit(Surv(time, status) ~ COO_class, data = metadata)
print(fit)
res.sum <- surv_summary(fit)
print(res.sum)
surv_diff <- survdiff(Surv(time, status) ~ COO_class, data = metadata)
print(surv_diff)


pdf("plots/06c-survival_nci_coo.pdf", height = 7, width=15)

coo <-
  ggsurvplot(fit,
             risk.table =  "abs_pct",
             pval = TRUE,
             break.time.by = 500,
             conf.int = FALSE,
             risk.table.col = "strata",
             risk.table.y.text.col = T,
             risk.table.y.text = FALSE,
             linetype = "strata",
             surv.median.line = "hv",
             ncensor.plot = FALSE,
             ggtheme = theme_bw(), # Change ggplot2 theme
             palette = c("red3", "royalblue", "lightblue"),
             legend.labs = 
               c("ABC", "GCB", "Unclassified"),
             font.x = c(20),
             font.y = c(20),
             font.tickslab = c(20),
             size = 2,
             tables.theme = clean_theme()) 

coo$plot <- coo$plot + 
  theme(legend.text = element_text(size = 30, color = "black", face = "bold"))

coo

dev.off()
