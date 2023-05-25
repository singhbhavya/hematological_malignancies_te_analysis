################################################################################
################################################################################
################################################################################
################################################################################
################################# HTMCP DATA QC ################################


#################################### SETUP #####################################

library(tidyverse)
library(readxl)
library(GenomicDataCommons)
library(dplyr)
library(rtracklayer)
library(data.table)
library(edgeR)

################################## LOAD DATA ###################################

load("r_outputs/07-htmcp_all_lymphoma_counts.Rdata")
load("r_outputs/07-htmcp_dlbcl_counts.Rdata")

################################ SANITY CHECK ##################################

## Count the number of genes, HERVs, and L1s 
## (HERVs and L1s should add up to the number of retro genes AND annotation)
## Really, this should all be the same for combined, TCGA, and NCI. 

cat(sprintf("%d genes in all datasets\n", nrow(all.counts.hiv.tx)))
cat(sprintf("%d retro genes in all datasets\n", nrow(all.counts.hiv.rtx)))
cat(sprintf("%d HERVs in all datasets\n", nrow(all.counts.hiv.herv)))
cat(sprintf("%d L1s in all datasets\n", nrow(all.counts.hiv.l1)))
cat(sprintf("%d combined genes in all datasets\n", nrow(all.counts.hiv.comb)))

cat(sprintf("%d genes in all datasets\n", nrow(DLBCL.hiv.counts.tx)))
cat(sprintf("%d retro genes in all datasets\n", nrow(DLBCL.hiv.counts.rtx)))
cat(sprintf("%d HERVs in all datasets\n", nrow(DLBCL.counts.hiv.herv)))
cat(sprintf("%d L1s in all datasets\n", nrow(DLBCL.counts.hiv.l1)))
cat(sprintf("%d combined genes in all datasets\n", nrow(DLBCL.hiv.counts.comb)))

############################## FILTER ALL DATA #################################

# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(all.counts.hiv.rtx) * 0.005)

all.hiv.filt.tx <- all.counts.hiv.tx[rowSums(all.counts.hiv.tx > cutoff.count) > cutoff.samp,]
all.hiv.filt.rtx <- all.counts.hiv.rtx[rowSums(all.counts.hiv.rtx > cutoff.count) > cutoff.samp, ]
all.hiv.filt.herv <- all.counts.hiv.herv[rowSums(all.counts.hiv.herv > cutoff.count) > cutoff.samp, ]
all.hiv.filt.l1 <- all.counts.hiv.l1[rowSums(all.counts.hiv.l1 > cutoff.count) >cutoff.samp, ]
all.hiv.filt.comb <- all.counts.hiv.comb[rowSums(all.counts.hiv.comb > cutoff.count) > cutoff.samp,]

################################ FILTER DLBCL ##################################

# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- 3

DLBCL.hiv.filt.tx <- DLBCL.hiv.counts.tx[rowSums(DLBCL.hiv.counts.tx > cutoff.count) > cutoff.samp, ]
DLBCL.hiv.filt.rtx <- DLBCL.hiv.counts.rtx[rowSums(DLBCL.hiv.counts.rtx > cutoff.count) > cutoff.samp, ]
DLBCL.hiv.filt.herv <- DLBCL.counts.hiv.herv[rowSums(DLBCL.counts.hiv.herv > cutoff.count) > cutoff.samp, ]
DLBCL.hiv.filt.l1 <- DLBCL.counts.hiv.l1[rowSums(DLBCL.counts.hiv.l1 > cutoff.count) >cutoff.samp, ]
DLBCL.hiv.filt.comb <- DLBCL.hiv.counts.comb[rowSums(DLBCL.hiv.counts.comb > cutoff.count) >cutoff.samp, ]

################################ SANITY CHECK ##################################

## Count the number of genes, HERVs, and L1s 
## (HERVs and L1s should add up to the number of retro genes AND annotation)
## Really, this should all be the same for combined, TCGA, and NCI. 

cat(sprintf("%d genes in all datasets\n", nrow(all.hiv.filt.tx)))
cat(sprintf("%d retro genes in all datasets\n", nrow(all.hiv.filt.rtx)))
cat(sprintf("%d HERV in all datasets\n", nrow(all.hiv.filt.herv)))
cat(sprintf("%d L1 in all datasets\n\n", nrow(all.hiv.filt.l1)))
cat(sprintf("%d combined genes in all datasets\n\n", nrow(all.hiv.filt.comb)))

cat(sprintf("%d genes in all datasets\n", nrow(DLBCL.hiv.filt.tx)))
cat(sprintf("%d retro genes in all datasets\n", nrow(DLBCL.hiv.filt.rtx)))
cat(sprintf("%d HERV in all datasets\n", nrow(DLBCL.hiv.filt.herv)))
cat(sprintf("%d L1 in all datasets\n\n", nrow(DLBCL.hiv.filt.l1)))
cat(sprintf("%d combined genes in all datasets\n\n", nrow(DLBCL.hiv.filt.comb)))

############################# COUNT READS PER TYPE ############################# 

te_percent <- function(herv.df, rtx.df, comb.df, metadata, metadata.col) {
  
  herv.reads <- as.data.frame(colSums(herv.df))
  te.reads <- as.data.frame(colSums(rtx.df))
  all.reads <- as.data.frame(colSums(comb.df))
  
  colnames(herv.reads) <- c("reads")
  colnames(te.reads) <- c("reads")
  colnames(all.reads) <- c("reads")
  
  herv.reads$sample  <- rownames(herv.reads)
  te.reads$sample <- rownames(te.reads)
  all.reads$sample <- rownames(all.reads)
  
  herv.reads$type <- metadata[[metadata.col]]
  te.reads$type <- metadata[[metadata.col]]
  all.reads$type <- metadata[[metadata.col]]
  
  stopifnot(all(all.reads$sample == te.reads$sample))
  
  te.reads$proportion <- te.reads$reads/all.reads$reads*100
  herv.reads$proportion <- herv.reads$reads/all.reads$reads*100
  
  output <- list(herv.reads = herv.reads,
                 te.reads = te.reads)
  
  return(output)
}


####################### COUNT READS PER LYMPHOMA SUB-TYPE ###################### 

all_metadata_hiv$subtype <- replace_na(all_metadata$subtype, replace = "Unclass")
lymphoma.te.percent <-
  te_percent(all.hiv.filt.herv, all.hiv.filt.rtx, all.hiv.filt.comb,
             all_metadata_hiv, "subtype")

lymphoma.te.percent$herv.reads
lymphoma.te.percent$te.reads

lymphoma.te.percent$herv.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Cell Type") +
  ylab("% of TE Fragments") +
  ylim(0, 2.5) +
  scale_fill_manual(values = c("Endemic BL EBV-positive" = wes_palette("Darjeeling1")[2], 
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
                               "HIV-positive" = "#fa7de7")) + 
  theme(aspect.ratio = 1)


################################## SAVE FILES ##################################

save(all.hiv.filt.comb, all.hiv.filt.tx, all.hiv.filt.rtx, 
     all.hiv.filt.herv, all.hiv.filt.l1, 
     all_metadata_hiv,
     file="r_outputs/07-htmcp_all_lymphoma_filt_counts.Rdata")


save(DLBCL.hiv.filt.comb, DLBCL.hiv.filt.tx, DLBCL.hiv.filt.rtx, 
     DLBCL.hiv.filt.herv, DLBCL.hiv.filt.l1,
     DLBCL_metadata_hiv,
     file="r_outputs/07-htmcp_DLBCL_filt_counts.Rdata")
