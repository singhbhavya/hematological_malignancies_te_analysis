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

cat(sprintf("%d genes in all datasets\n", nrow(all.counts.hiv.rtx)))
cat(sprintf("%d retro genes in all datasets\n", nrow(all.counts.hiv.herv)))
cat(sprintf("%d combined genes in all datasets\n", nrow(all.counts.hiv.l1)))

cat(sprintf("%d genes in DLBCL\n", nrow(DLBCL.hiv.counts.rtx)))
cat(sprintf("%d retro genes in DLBCL\n", nrow(DLBCL.counts.hiv.herv)))
cat(sprintf("%d combined genes in DLBCL\n", nrow(DLBCL.counts.hiv.l1)))

############################## FILTER ALL DATA #################################

# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(all.counts.hiv.rtx) * 0.005)

all.hiv.filt.rtx <- all.counts.hiv.rtx[rowSums(all.counts.hiv.rtx > cutoff.count) > cutoff.samp, ]
all.hiv.filt.herv <- all.counts.hiv.herv[rowSums(all.counts.hiv.herv > cutoff.count) > cutoff.samp, ]
all.hiv.filt.l1 <- all.counts.hiv.l1[rowSums(all.counts.hiv.l1 > cutoff.count) >cutoff.samp, ]


################################ FILTER DLBCL ##################################

# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- 3

DLBCL.hiv.filt.rtx <- DLBCL.hiv.counts.rtx[rowSums(DLBCL.hiv.counts.rtx > cutoff.count) > cutoff.samp, ]
DLBCL.hiv.filt.herv <- DLBCL.counts.hiv.herv[rowSums(DLBCL.counts.hiv.herv > cutoff.count) > cutoff.samp, ]
DLBCL.hiv.filt.l1 <- DLBCL.counts.hiv.l1[rowSums(DLBCL.counts.hiv.l1 > cutoff.count) >cutoff.samp, ]


################################ SANITY CHECK ##################################

## Count the number of genes, HERVs, and L1s 
## (HERVs and L1s should add up to the number of retro genes AND annotation)
## Really, this should all be the same for combined, TCGA, and NCI. 

cat(sprintf("%d retro genes in all datasets\n", nrow(all.hiv.filt.rtx)))
cat(sprintf("%d HERV in all datasets\n", nrow(all.hiv.filt.herv)))
cat(sprintf("%d L1 in all datasets\n\n", nrow(all.hiv.filt.l1)))

cat(sprintf("%d retro genes in all datasets\n", nrow(DLBCL.hiv.filt.rtx)))
cat(sprintf("%d HERV in all datasets\n", nrow(DLBCL.hiv.filt.herv)))
cat(sprintf("%d L1 in all datasets\n\n", nrow(DLBCL.hiv.filt.l1)))

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

# need to do once we have combined counts

################################## SAVE FILES ##################################

save(all.hiv.filt.rtx, all.hiv.filt.herv, all.hiv.filt.l1, 
     all_metadata_hiv,
     file="r_outputs/07-htmcp_all_lymphoma_filt_counts.Rdata")


save(DLBCL.hiv.filt.rtx, DLBCL.hiv.filt.herv, DLBCL.hiv.filt.l1,
     DLBCL_metadata_hiv,
     file="r_outputs/07-htmcp_DLBCL_filt_counts.Rdata")
