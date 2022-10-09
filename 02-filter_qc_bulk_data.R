################################################################################
################################################################################
################################################################################
################################################################################
#################################  BULK DATA QC ################################

## Plan:
## - Count the number of genes, LINEs, HERVs.
## - Filter the datasets separately.
## - Filter the datasets together.

#################################### SETUP #####################################

library(tidyverse)
library(readxl)
library(GenomicDataCommons)
library(dplyr)
library(rtracklayer)
library(data.table)

################################## LOAD DATA ###################################

load("r_outputs/01-all_lymphoma_counts.Rdata")
load("r_outputs/01-DLBCL_counts.Rdata")
load("r_outputs/01-BL_counts.Rdata")
load("r_outputs/01-FL_counts.Rdata")

################################ SANITY CHECK ##################################

## Count the number of genes, HERVs, and L1s 
## (HERVs and L1s should add up to the number of retro genes AND annotation)
## Really, this should all be the same for combined, TCGA, and NCI. 

cat(sprintf("%d genes in all datasets\n", nrow(all.counts.tx)))
cat(sprintf("%d retro genes in all datasets\n", nrow(all.counts.rtx)))
cat(sprintf("%d combined genes in all datasets\n", nrow(all.counts.comb)))
cat(sprintf("%d HERV in all datasets\n", nrow(all.counts.herv)))
cat(sprintf("%d L1 in all datasets\n\n", nrow(all.counts.l1)))

cat(sprintf("%d genes in DLBCL\n", nrow(DLBCL.counts.tx)))
cat(sprintf("%d retro genes in DLBCL\n", nrow(DLBCL.counts.rtx)))
cat(sprintf("%d combined genes in DLBCL\n", nrow(DLBCL.counts.comb)))
cat(sprintf("%d HERV in DLBCL\n", nrow(DLBCL.counts.herv)))
cat(sprintf("%d L1 in DLBCL\n\n", nrow(DLBCL.counts.herv)))

cat(sprintf("%d genes in BL\n", nrow(BL.counts.tx)))
cat(sprintf("%d retro genes in BL\n", nrow(BL.counts.rtx)))
cat(sprintf("%d combined genes in BL\n", nrow(BL.counts.comb)))
cat(sprintf("%d HERV in BL\n", nrow(BL.counts.herv)))
cat(sprintf("%d L1 in BL\n\n", nrow(BL.counts.l1)))

cat(sprintf("%d genes in FL\n", nrow(FL.counts.tx)))
cat(sprintf("%d retro genes in FL\n", nrow(FL.counts.rtx)))
cat(sprintf("%d combined genes FL\n", nrow(FL.counts.comb)))
cat(sprintf("%d HERV in FL\n", nrow(FL.counts.herv)))
cat(sprintf("%d L1 in FL\n\n", nrow(FL.counts.l1)))

############################## FILTER ALL DATA #################################

# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(all.counts.comb) * 0.01)

all.counts.mfilt.tx <- all.counts.tx[rowSums(all.counts.tx) > cutoff.count, ]
all.counts.mfilt.rtx <- all.counts.rtx[rowSums(all.counts.rtx) > cutoff.count, ]
all.counts.mfilt.comb <- rbind(all.counts.mfilt.tx, all.counts.mfilt.rtx)
all.counts.mfilt.herv <- all.counts.herv[rowSums(all.counts.herv) > cutoff.count, ]
all.counts.mfilt.l1 <- all.counts.l1[rowSums(all.counts.l1) > cutoff.count, ]

all.counts.filt.tx <- all.counts.tx[rowSums(all.counts.tx > cutoff.count) > cutoff.samp, ]
all.counts.filt.rtx <- all.counts.rtx[rowSums(all.counts.rtx > cutoff.count) > cutoff.samp, ]
all.counts.filt.comb <- rbind(all.counts.filt.tx, all.counts.filt.rtx)
all.counts.filt.herv <- all.counts.herv[rowSums(all.counts.herv > cutoff.count) > cutoff.samp, ]
all.counts.filt.l1 <- all.counts.l1[rowSums(all.counts.l1 > cutoff.count) >cutoff.samp, ]

################################ FILTER DLBCL ##################################

# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(DLBCL.counts.comb) * 0.01)

DLBCL.counts.mfilt.tx <- DLBCL.counts.tx[rowSums(DLBCL.counts.tx) > cutoff.count, ]
DLBCL.counts.mfilt.rtx <- DLBCL.counts.rtx[rowSums(DLBCL.counts.rtx) > cutoff.count, ]
DLBCL.counts.mfilt.comb <- rbind(DLBCL.counts.mfilt.tx, DLBCL.counts.mfilt.rtx)
DLBCL.counts.mfilt.herv <- DLBCL.counts.herv[rowSums(DLBCL.counts.herv) > cutoff.count, ]
DLBCL.counts.mfilt.l1 <- DLBCL.counts.l1[rowSums(DLBCL.counts.l1) > cutoff.count, ]

DLBCL.filt.tx <- DLBCL.counts.tx[rowSums(DLBCL.counts.tx > cutoff.count) > cutoff.samp, ]
DLBCL.filt.rtx <- DLBCL.counts.rtx[rowSums(DLBCL.counts.rtx > cutoff.count) > cutoff.samp, ]
DLBCL.filt.comb <- rbind(DLBCL.filt.tx, DLBCL.filt.rtx)
DLBCL.filt.herv <- DLBCL.counts.herv[rowSums(DLBCL.counts.herv > cutoff.count) > cutoff.samp, ]
DLBCL.filt.l1 <- DLBCL.counts.l1[rowSums(DLBCL.counts.l1 > cutoff.count) >cutoff.samp, ]

################################## FILTER BL ###################################

# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(BL.counts.comb) * 0.05)

BL.counts.mfilt.tx <- BL.counts.tx[rowSums(BL.counts.tx) > cutoff.count, ]
BL.counts.mfilt.rtx <- BL.counts.rtx[rowSums(BL.counts.rtx) > cutoff.count, ]
BL.counts.mfilt.comb <- rbind(BL.counts.mfilt.tx, BL.counts.mfilt.rtx)
BL.counts.mfilt.herv <- BL.counts.herv[rowSums(BL.counts.herv) > cutoff.count, ]
BL.counts.mfilt.l1 <- BL.counts.l1[rowSums(BL.counts.l1) > cutoff.count, ]

BL.filt.tx <- BL.counts.tx[rowSums(BL.counts.tx > cutoff.count) > cutoff.samp, ]
BL.filt.rtx <- BL.counts.rtx[rowSums(BL.counts.rtx > cutoff.count) > cutoff.samp, ]
BL.filt.comb <- rbind(BL.filt.tx, BL.filt.rtx)
BL.filt.herv <- BL.counts.herv[rowSums(BL.counts.herv > cutoff.count) > cutoff.samp, ]
BL.filt.l1 <- BL.counts.l1[rowSums(BL.counts.l1 > cutoff.count) >cutoff.samp, ]

################################## FILTER FL ###################################

# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(FL.counts.comb) * 0.1)

FL.counts.mfilt.tx <- FL.counts.tx[rowSums(FL.counts.tx) > cutoff.count, ]
FL.counts.mfilt.rtx <- FL.counts.rtx[rowSums(FL.counts.rtx) > cutoff.count, ]
FL.counts.mfilt.comb <- rbind(FL.counts.mfilt.tx, FL.counts.mfilt.rtx)
FL.counts.mfilt.herv <- FL.counts.herv[rowSums(FL.counts.herv) > cutoff.count, ]
FL.counts.mfilt.l1 <- FL.counts.l1[rowSums(FL.counts.l1) > cutoff.count, ]

FL.filt.tx <- FL.counts.tx[rowSums(FL.counts.tx > cutoff.count) > cutoff.samp, ]
FL.filt.rtx <- FL.counts.rtx[rowSums(FL.counts.rtx > cutoff.count) > cutoff.samp, ]
FL.filt.comb <- rbind(FL.filt.tx, FL.filt.rtx)
FL.filt.herv <- FL.counts.herv[rowSums(FL.counts.herv > cutoff.count) > cutoff.samp, ]
FL.filt.l1 <- FL.counts.l1[rowSums(FL.counts.l1 > cutoff.count) >cutoff.samp, ]

################################## SAVE FILES ##################################

save(all.counts.filt.comb, all.counts.filt.tx, all.counts.filt.rtx, 
     all.counts.filt.herv, all.counts.filt.l1, all_metadata, 
     all.counts.mfilt.tx, all.counts.mfilt.rtx, all.counts.mfilt.comb,
     all.counts.mfilt.herv, all.counts.mfilt.l1,
     file="r_outputs/02-all_lymphoma_filt_counts.Rdata")


save(DLBCL.filt.comb, DLBCL.filt.tx, DLBCL.filt.rtx, 
     DLBCL.filt.herv, DLBCL.filt.l1, DLBCL_metadata, 
     DLBCL.counts.mfilt.tx, DLBCL.counts.mfilt.rtx,
     DLBCL.counts.mfilt.comb, DLBCL.counts.mfilt.herv,
     DLBCL.counts.mfilt.l1,
     file="r_outputs/02-DLBCL_filt_counts.Rdata")

save(BL.filt.comb, BL.filt.tx, BL.filt.rtx, 
     BL.filt.herv, BL.filt.l1, BL_metadata, 
     BL.counts.mfilt.tx, BL.counts.mfilt.rtx,
     BL.counts.mfilt.comb, BL.counts.mfilt.herv,
     BL.counts.mfilt.l1,
     file="r_outputs/02-BL_filt_counts.Rdata")

save(FL.filt.comb, FL.filt.tx, FL.filt.rtx, 
     FL.filt.herv, FL.filt.l1, FL_metadata, 
     FL.counts.mfilt.tx, FL.counts.mfilt.rtx,
     FL.counts.mfilt.comb, FL.counts.mfilt.herv,
     FL.counts.mfilt.l1,
     file="r_outputs/02-FL_filt_counts.Rdata")

