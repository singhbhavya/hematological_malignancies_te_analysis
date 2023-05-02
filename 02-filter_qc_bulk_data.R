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
library(edgeR)

################################## LOAD DATA ###################################

load("r_outputs/01-all_lymphoma_counts.Rdata")
load("r_outputs/01-DLBCL_counts.Rdata")
load("r_outputs/01-BL_counts.Rdata")
load("r_outputs/01-FL_counts.Rdata")
load("r_outputs/01-GCB_Bulk_counts.Rdata")
load("r_outputs/01-GCB_Agirre.Rdata")
load("r_outputs/01-refs.Rdata")

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

cat(sprintf("%d genes in GCB Bulk\n", nrow(GCB_Bulk.counts.tx)))
cat(sprintf("%d retro genes in GCB Bulk\n", nrow(GCB_Bulk.counts.rtx)))
cat(sprintf("%d combined genes GCB Bulk\n", nrow(GCB_Bulk.counts.comb)))
cat(sprintf("%d HERV in GCB Bulk\n", nrow(GCB_Bulk.counts.herv)))
cat(sprintf("%d L1 in GCB Bulk\n", nrow(GCB_Bulk.counts.l1)))

cat(sprintf("%d genes in GCB Agirre\n", nrow(GCB_Agirre.counts.tx)))
cat(sprintf("%d retro genes in GCB Agirre\n", nrow(GCB_Agirre.counts.rtx)))
cat(sprintf("%d combined genes GCB Agirre\n", nrow(GCB_Agirre.counts.comb)))
cat(sprintf("%d HERV in GCB Agirre\n", nrow(GCB_Agirre.counts.herv)))
cat(sprintf("%d L1 in GCB Agirre\n", nrow(GCB_Agirre.counts.l1)))

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

############################### FILTER GCB Bulk ################################

# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(GCB_Bulk.counts.comb) * 0.1)

GCB_Bulk.counts.mfilt.tx <- GCB_Bulk.counts.tx[rowSums(GCB_Bulk.counts.tx) > cutoff.count, ]
GCB_Bulk.counts.mfilt.rtx <- GCB_Bulk.counts.rtx[rowSums(GCB_Bulk.counts.rtx) > cutoff.count, ]
GCB_Bulk.counts.mfilt.comb <- rbind(GCB_Bulk.counts.mfilt.tx, GCB_Bulk.counts.mfilt.rtx)
GCB_Bulk.counts.mfilt.herv <- GCB_Bulk.counts.herv[rowSums(GCB_Bulk.counts.herv) > cutoff.count, ]
GCB_Bulk.counts.mfilt.l1 <- GCB_Bulk.counts.l1[rowSums(GCB_Bulk.counts.l1) > cutoff.count, ]

GCB_Bulk.filt.tx <- GCB_Bulk.counts.tx[rowSums(GCB_Bulk.counts.tx > cutoff.count) > cutoff.samp, ]
GCB_Bulk.filt.rtx <- GCB_Bulk.counts.rtx[rowSums(GCB_Bulk.counts.rtx > cutoff.count) > cutoff.samp, ]
GCB_Bulk.filt.comb <- rbind(GCB_Bulk.filt.tx, GCB_Bulk.filt.rtx)
GCB_Bulk.filt.herv <- GCB_Bulk.counts.herv[rowSums(GCB_Bulk.counts.herv > cutoff.count) > cutoff.samp, ]
GCB_Bulk.filt.l1 <- GCB_Bulk.counts.l1[rowSums(GCB_Bulk.counts.l1 > cutoff.count) >cutoff.samp, ]

############################## FILTER GCB Agirre ###############################

# Minimum count threshold
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold
cutoff.samp <- floor(ncol(GCB_Agirre.counts.comb) * 0.1)

GCB_Agirre.counts.mfilt.tx <- GCB_Agirre.counts.tx[rowSums(GCB_Agirre.counts.tx) > cutoff.count, ]
GCB_Agirre.counts.mfilt.rtx <- GCB_Agirre.counts.rtx[rowSums(GCB_Agirre.counts.rtx) > cutoff.count, ]
GCB_Agirre.counts.mfilt.comb <- rbind(GCB_Agirre.counts.mfilt.tx, GCB_Agirre.counts.mfilt.rtx)
GCB_Agirre.counts.mfilt.herv <- GCB_Agirre.counts.herv[rowSums(GCB_Agirre.counts.herv) > cutoff.count, ]
GCB_Agirre.counts.mfilt.l1 <- GCB_Agirre.counts.l1[rowSums(GCB_Agirre.counts.l1) > cutoff.count, ]

GCB_Agirre.filt.tx <- GCB_Agirre.counts.tx[rowSums(GCB_Agirre.counts.tx > cutoff.count) > cutoff.samp, ]
GCB_Agirre.filt.rtx <- GCB_Agirre.counts.rtx[rowSums(GCB_Agirre.counts.rtx > cutoff.count) > cutoff.samp, ]
GCB_Agirre.filt.comb <- rbind(GCB_Agirre.filt.tx, GCB_Agirre.filt.rtx)
GCB_Agirre.filt.herv <- GCB_Agirre.counts.herv[rowSums(GCB_Agirre.counts.herv > cutoff.count) > cutoff.samp, ]
GCB_Agirre.filt.l1 <- GCB_Agirre.counts.l1[rowSums(GCB_Agirre.counts.l1 > cutoff.count) >cutoff.samp, ]

################################ SANITY CHECK ##################################

## Count the number of genes, HERVs, and L1s 
## (HERVs and L1s should add up to the number of retro genes AND annotation)
## Really, this should all be the same for combined, TCGA, and NCI. 

cat(sprintf("%d genes in all datasets\n", nrow(all.counts.filt.tx)))
cat(sprintf("%d retro genes in all datasets\n", nrow(all.counts.filt.rtx)))
cat(sprintf("%d combined genes in all datasets\n", nrow(all.counts.filt.comb)))
cat(sprintf("%d HERV in all datasets\n", nrow(all.counts.filt.herv)))
cat(sprintf("%d L1 in all datasets\n\n", nrow(all.counts.filt.l1)))

cat(sprintf("%d genes in DLBCL\n", nrow(DLBCL.filt.tx)))
cat(sprintf("%d retro genes in DLBCL\n", nrow(DLBCL.filt.rtx)))
cat(sprintf("%d combined genes in DLBCL\n", nrow(DLBCL.filt.comb)))
cat(sprintf("%d HERV in DLBCL\n", nrow(DLBCL.filt.herv)))
cat(sprintf("%d L1 in DLBCL\n\n", nrow(DLBCL.filt.l1)))

cat(sprintf("%d genes in BL\n", nrow(BL.filt.tx)))
cat(sprintf("%d retro genes in BL\n", nrow(BL.filt.rtx)))
cat(sprintf("%d combined genes in BL\n", nrow(BL.filt.comb)))
cat(sprintf("%d HERV in BL\n", nrow(BL.filt.herv)))
cat(sprintf("%d L1 in BL\n\n", nrow(BL.filt.l1)))

cat(sprintf("%d genes in FL\n", nrow(FL.filt.tx)))
cat(sprintf("%d retro genes in FL\n", nrow(FL.filt.rtx)))
cat(sprintf("%d combined genes FL\n", nrow(FL.filt.comb)))
cat(sprintf("%d HERV in FL\n", nrow(FL.filt.herv)))
cat(sprintf("%d L1 in FL\n\n", nrow(FL.filt.l1)))

cat(sprintf("%d genes in GCB Bulk\n", nrow(GCB_Bulk.filt.tx)))
cat(sprintf("%d retro genes in GCB Bulk\n", nrow(GCB_Bulk.filt.rtx)))
cat(sprintf("%d combined genes GCB Bulk\n", nrow(GCB_Bulk.filt.comb)))
cat(sprintf("%d HERV in GCB Bulk\n", nrow(GCB_Bulk.filt.herv)))
cat(sprintf("%d L1 in GCB Bulk\n\n", nrow(GCB_Bulk.filt.l1)))

cat(sprintf("%d genes in GCB Agirre\n", nrow(GCB_Agirre.filt.tx)))
cat(sprintf("%d retro genes in GCB Agirre\n", nrow(GCB_Agirre.filt.rtx)))
cat(sprintf("%d combined genes Agirre\n", nrow(GCB_Agirre.filt.comb)))
cat(sprintf("%d HERV in GCB Agirre\n", nrow(GCB_Agirre.filt.herv)))
cat(sprintf("%d L1 in GCB Agirre\n\n", nrow(GCB_Agirre.filt.l1)))

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

####################### COUNT READS PER AGIRRE CELL TYPE ####################### 

agirre.te.percent <-
  te_percent(GCB_Agirre.filt.herv, GCB_Agirre.filt.rtx, GCB_Agirre.filt.comb,
           agirre_metadata, "source_name")

agirre.te.percent$herv.reads %>%
  group_by(type) %>%
  summarise_at(vars(proportion), list(name = mean))

agirre.te.percent$te.reads %>%
  group_by(type) %>%
  summarise_at(vars(proportion), list(name = mean))

pdf("plots/02-agirre_te_percent.pdf", height=4, width=4)
agirre.te.percent$te.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Cell Type") +
  ylab("% of TE Fragments") +
  ylim(0, 1) +
  scale_fill_manual(values = c("Naive" = pal_jco("default", alpha = 0.8)(7)[1],
                               "Bone Marrow plasma cell" = pal_jco("default", alpha = 0.8)(7)[7],
                               "Memory" = pal_jco("default", alpha = 0.8)(7)[3],
                               "Centroblast" = pal_jco("default", alpha = 0.8)(7)[4],
                               "Centrocyte" = pal_jco("default", alpha = 0.8)(7)[5],
                               "Tonsilar plasma cell" = pal_jco("default", alpha = 0.8)(7)[6])) + 
  scale_x_discrete(labels=c("Naive" = "NB", 
                            "Memory" = "MB",
                            "Centroblast" = "DZ",
                            "Centrocyte" = "LZ",
                            "Tonsilar plasma cell" = "PB",
                            "Bone Marrow plasma cell" = "BMPC" )) +
  theme(aspect.ratio = 1)

dev.off()

pdf("plots/02-agirre_herv_percent.pdf", height=4, width=4)
agirre.te.percent$herv.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Cell Type") +
  ylab("% of HERV Fragments") +
  ylim(0, 1) +
  scale_fill_manual(values = c("Naive" = pal_jco("default", alpha = 0.8)(7)[1],
                               "Bone Marrow plasma cell" = pal_jco("default", alpha = 0.8)(7)[7],
                               "Memory" = pal_jco("default", alpha = 0.8)(7)[3],
                               "Centroblast" = pal_jco("default", alpha = 0.8)(7)[4],
                               "Centrocyte" = pal_jco("default", alpha = 0.8)(7)[5],
                               "Tonsilar plasma cell" = pal_jco("default", alpha = 0.8)(7)[6])) + 
  scale_x_discrete(labels=c("Naive" = "NB", 
                            "Memory" = "MB",
                            "Centroblast" = "DZ",
                            "Centrocyte" = "LZ",
                            "Tonsilar plasma cell" = "PB",
                            "Bone Marrow plasma cell" = "BMPC"))  +
  theme(aspect.ratio = 1)

dev.off()

####################### COUNT READS PER HOLMES CELL TYPE ####################### 

holmes.te.percent <-
  te_percent(GCB_Bulk.filt.herv, GCB_Bulk.filt.rtx, GCB_Bulk.filt.comb,
             bulk_metadata, "source_name")

holmes.te.percent$herv.reads %>%
  group_by(type) %>%
  summarise_at(vars(proportion), list(name = mean))

holmes.te.percent$te.reads %>%
  group_by(type) %>%
  summarise_at(vars(proportion), list(name = mean))

pdf("plots/02-holmes_te_percent.pdf", height=4, width=4)
holmes.te.percent$te.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Cell Type") +
  ylab("% of TE Fragments") +
  ylim(0, 1) +
  scale_fill_manual(values = c("Naïve B cells" = pal_jco("default", alpha = 0.7)(5)[1],
                               "Germinal Center B cells" = pal_jco("default", alpha = 0.7)(5)[2],
                               "Memory B cells" = pal_jco("default", alpha = 0.7)(5)[3],
                               "Dark Zone Germinal Center B cells" = pal_jco("default", alpha = 0.7)(5)[4],
                               "Light Zone Germinal Center B cells" = pal_jco("default", alpha = 0.7)(5)[5])) + 
  scale_x_discrete(labels=c("Naïve B cells" = "NB", 
                            "Memory B cells" = "MB",
                            "Dark Zone Germinal Center B cells" = "DZ",
                            "Light Zone Germinal Center B cells" = "LZ",
                            "Germinal Center B cells" = "GCB")) +
  theme(aspect.ratio = 1)

dev.off()

pdf("plots/02-holmes_herv_percent.pdf", height=4, width=4)
holmes.te.percent$herv.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Cell Type") +
  ylab("% of HERV Fragments") +
  ylim(0, 1) +
  scale_fill_manual(values = c("Naïve B cells" = pal_jco("default", alpha = 0.7)(5)[1],
                               "Germinal Center B cells" = pal_jco("default", alpha = 0.7)(5)[2],
                               "Memory B cells" = pal_jco("default", alpha = 0.7)(5)[3],
                               "Dark Zone Germinal Center B cells" = pal_jco("default", alpha = 0.7)(5)[4],
                               "Light Zone Germinal Center B cells" = pal_jco("default", alpha = 0.7)(5)[5])) + 
  scale_x_discrete(labels=c("Naïve B cells" = "NB", 
                            "Memory B cells" = "MB",
                            "Dark Zone Germinal Center B cells" = "DZ",
                            "Light Zone Germinal Center B cells" = "LZ",
                            "Germinal Center B cells" = "GCB")) +
  theme(aspect.ratio = 1)

dev.off()

####################### COUNT READS PER LYMPHOMA SUB-TYPE ###################### 

all_metadata$subtype <- replace_na(all_metadata$subtype, replace = "Unclass")
lymphoma.te.percent <-
  te_percent(all.counts.filt.herv, all.counts.filt.rtx, all.counts.filt.comb,
             all_metadata, "subtype")

lymphoma.te.percent$herv.reads
lymphoma.te.percent$te.reads

pdf("plots/02-lymphoma_subtype_te_percent.pdf", height=4, width=4)
lymphoma.te.percent$te.reads %>%
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
                               "FOLLICULAR GRADE 3A" = "#5B1A18")) + 
  theme(aspect.ratio = 1)
dev.off()

pdf("plots/02-lymphoma_subtype_herv_percent.pdf", height=4, width=4)
lymphoma.te.percent$herv.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Cell Type") +
  ylab("% of HERV Fragments") +
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
                               "FOLLICULAR GRADE 3A" = "#5B1A18")) + 
  theme(aspect.ratio = 1)
dev.off()

####################### COUNT READS PER BURKITT LYMPHOMA ####################### 

BL.te.percent <-
  te_percent(BL.filt.herv, BL.filt.rtx, BL.filt.comb,
             BL_metadata, "ebv_status")

BL.te.percent$herv.reads %>%
  group_by(type) %>%
  summarise_at(vars(proportion), list(name = mean))

BL.te.percent$te.reads %>%
  group_by(type) %>%
  summarise_at(vars(proportion), list(name = mean))

pdf("plots/02-BL_ebv_te_percent.pdf", height=4, width=4)
BL.te.percent$te.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("EBV status") +
  ylab("% of TE Fragments") +
  ylim(0, 3) +
  scale_fill_manual(values = c("EBV-negative" = wes_palette("Zissou1")[1], 
                               "EBV-positive" = wes_palette("Zissou1")[4])) + 
  theme(aspect.ratio = 1)
dev.off()

pdf("plots/02-BL_ebv_herv_percent.pdf", height=4, width=4)
BL.te.percent$herv.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("EBV status") +
  ylab("% of HERV Fragments") +
  ylim(0, 3) +
  scale_fill_manual(values = c("EBV-negative" = wes_palette("Zissou1")[1], 
                               "EBV-positive" = wes_palette("Zissou1")[4])) + 
  theme(aspect.ratio = 1)
dev.off()
####################### COUNT READS PER LYMPHOMA SUB-TYPE ###################### 

all_metadata$subtype <- replace_na(all_metadata$subtype, replace = "Unclass")
lymphoma.te.percent <-
  te_percent(all.counts.filt.herv, all.counts.filt.rtx, all.counts.filt.comb,
             all_metadata, "cancer_type")

lymphoma.te.percent$herv.reads
lymphoma.te.percent$te.reads

pdf("plots/02-lymphoma_type_te_percent.pdf", height=4, width=4)
lymphoma.te.percent$te.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Cell Type") +
  ylab("% of TE Fragments") +
  ylim(0, 2.5) +
  scale_fill_manual(values = c("DLBCL" = "0073C2B2", 
                               "BL" = "#EFC000B2",
                               "FL" = "#868686B2")) + 
  theme(aspect.ratio = 1)
dev.off()

pdf("plots/02-lymphoma_type_herv_percent.pdf", height=4, width=4)
lymphoma.te.percent$herv.reads %>%
  ggplot(aes(x=type, y=proportion, fill=type))  +
  geom_boxplot(notch = TRUE) +
  theme_pubr() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=45, hjust=1)) +
  xlab("Cell Type") +
  ylab("% of HERV Fragments") +
  ylim(0, 2.5) +
  scale_fill_manual(values = c("DLBCL" = "0073C2B2", 
                               "BL" = "#EFC000B2",
                               "FL" = "#868686B2")) + 
  theme(aspect.ratio = 1)
dev.off()

################################## WRITE CSVs ##################################

write.csv(DLBCL.filt.comb, "r_outputs/dlbcl.filt.comb.csv", quote=FALSE)
write.csv(BL.filt.comb, "r_outputs/BL.filt.comb.csv", quote=FALSE)
write.csv(FL.filt.comb, "r_outputs/FL.filt.comb.csv", quote=FALSE)
write.csv(GCB_Bulk.filt.comb, "r_outputs/BHM.filt.comb.csv", quote=FALSE)
write.csv(GCB_Agirre.filt.comb, "r_outputs/BAG.filt.comb.csv", quote=FALSE)

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

save(GCB_Bulk.filt.comb, GCB_Bulk.filt.tx, GCB_Bulk.filt.rtx, 
     GCB_Bulk.filt.herv, GCB_Bulk.filt.l1, bulk_metadata, 
     GCB_Bulk.counts.mfilt.comb, GCB_Bulk.counts.mfilt.tx,
     GCB_Bulk.counts.mfilt.rtx, GCB_Bulk.counts.mfilt.herv, 
     GCB_Bulk.counts.mfilt.l1,
     file="r_outputs/02-GCB_Bulk_filt_counts.Rdata")

save(GCB_Agirre.filt.comb, GCB_Agirre.filt.tx, GCB_Agirre.filt.rtx, 
     GCB_Agirre.filt.herv, GCB_Agirre.filt.l1, agirre_metadata, 
     GCB_Agirre.counts.mfilt.comb, GCB_Agirre.counts.mfilt.tx,
     GCB_Agirre.counts.mfilt.rtx, GCB_Agirre.counts.mfilt.herv, 
     GCB_Agirre.counts.mfilt.l1,
     file="r_outputs/02-GCB_Agirre_filt_counts.Rdata")


