################################################################################
################################################################################
################################################################################
################################################################################
######################### ALL LYMPHOMA CLUSTER FEATURES ########################

## Plan:
##

#################################### SETUP #####################################

library(tidyverse)
library(DESeq2)
library(matrixStats)
library(glmnet)
library(Boruta)
library(UpSetR)
library(c060)
library(randomForest)
library(rpart)
library(rpart.plot)

################################### LOAD DATA ##################################

load("r_outputs/02-all_lymphoma_filt_counts.Rdata")
load("r_outputs/01-refs.Rdata")

#################################### DESEQ2 ####################################

all_metadata$subtype <- replace_na(all_metadata$subtype, replace = "Unclass")

hervs.to.keep <- intersect(rownames(all.counts.filt.herv), 
                           retro.annot$locus[retro.annot$chrom != "chrY"])

all.counts.filt.herv.y <- all.counts.filt.herv[hervs.to.keep,] 

countDat <- all.counts.filt.herv.y

cat(sprintf('%d variables\n', nrow(countDat)))

stopifnot(all(colnames(countDat) == rownames(all_metadata)))

dds <- DESeq2::DESeqDataSetFromMatrix(countDat, all_metadata, ~1 )
dds <- DESeq2::DESeq(dds, parallel=T)
tform <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)

##################################### LRT ######################################

p.cutoff <- 1e-3
lfc.cutoff <- 1.5
dds_lrt <- DESeq2::DESeqDataSetFromMatrix(countDat, all_metadata, ~ cancer_type)
dds_lrt <- DESeq2::DESeq(dds_lrt, test="LRT", reduced = ~ 1, parallel=T)
res_lrt <- DESeq2::results(dds_lrt, alpha=p.cutoff)
sig_lrt <- subset(res_lrt, padj < p.cutoff & abs(log2FoldChange) > lfc.cutoff)

selected_vars <- list()
selected_vars$lrt <- rownames(sig_lrt)

############################### FILTER VARIABLES ###############################

ntop <- nrow(tform) # No filtering
vars <- rowVars(assay(tform))
mat <- assay(tform)[order(-vars)[1:ntop], ]

cat(sprintf('%d variables\n', nrow(mat)))

#################################### BORUTA ####################################

all_metadata$cancer_type <- as.factor(all_metadata$cancer_type)
set.seed(12345)
bor.orig <- Boruta(x=t(mat), y=all_metadata$cancer_type, doTrace=2, ntree=1000, maxRuns=1000)
print(bor.orig)
bor.model <- TentativeRoughFix(bor.orig)
print(bor.model)

selected_vars$boruta <- names(bor.model$finalDecision)[bor.model$finalDecision == 'Confirmed']

#################################### LASSO #####################################

weakness <- 1    # For each subsample the features are reweighted by a random weight uniformly sampled in [weakness,1]
size <- 0.632      # proportion of samples drawn in every subsample
steps <- 200       # number of subsamples

error <- 0.01  # the desired type I error level w.r.t. to the chosen type I error rate.
pi_thr <- 0.6      # the threshold used for the stability selection

ncores <- 1        # force run on 1 core, prevents reseeding RNG on each core

set.seed(123)

res.stabpath.c <-
  c060::stabpath(
    all_metadata$cancer_type,
    t(mat),
    mc.cores = ncores,
    weakness=weakness,
    family = "multinomial",
    type.multinomial = "grouped"
  )

plot(res.stabpath.c)

res.stabsel.c <-
  c060::stabsel(
    res.stabpath.c,
    error = 0.1,
    pi_thr = 0.6,
    type = "pfer"
  )

length(res.stabsel.c$stable)

selected_vars$lasso <- names(res.stabsel.c$stable)

pdf("plots/05m-all_lymphoma_features_upset.pdf", height=5, width=7.5)
UpSetR::upset(
  fromList(selected_vars), 
  sets=c('lrt', 'boruta', 'lasso'), 
  order.by = "freq",
  text.scale = c(2, 2, 2, 2, 2, 2))
dev.off()

################################# FEATURE GRAPH ################################

mat.sel <- mat[selected_vars$lasso, ]
rpart.fit <- rpart(y ~ . , data=data.frame(y=all_metadata$cancer_type, 
                                           t(mat.sel)),
                   method = 'class')

pdf("plots/05m-all_lymphoma_features_rpart.pdf", height=7, width=9)
rpart.plot(rpart.fit, type=5, extra=1, digits=3)
dev.off()

############################### SELECTED FEATURES ##############################

# > selected_vars$lasso
# [1] "ERVL_1p34.2"  "ERVLB4_2p16.3" "MER4B_10q21.3" "ERVL_Xq21.1b"  "ERVLE_14q23.2"

plot.counts <- function(df, gene) {
  
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="subtype", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=subtype, y=count, fill=subtype))  +
    geom_boxplot(notch = TRUE) +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("subtype") +
    ylab("Counts") +
    ggtitle(gene) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1) +
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
    scale_y_log10(labels = label_comma()) 
}


p1 <- plot.counts(dds_lrt, "ERVL_1p34.2") + stat_compare_means(comparisons = list(c("ABC", "Endemic BL EBV-positive")))
p2 <- plot.counts(dds_lrt, "ERVLB4_2p16.3")
p3 <- plot.counts(dds_lrt, "MER4B_10q21.3")
p4 <- plot.counts(dds_lrt, "ERVL_Xq21.1b")
p5 <- plot.counts(dds_lrt, "ERVLE_14q23.2")

pdf("plots/05m-all_lymphoma_selected_features.pdf", height=4, width=18)
plot_grid(p1, p2, p3, p4, p5,
          nrow = 2, 
          ncol = 3,
          labels = "AUTO")
dev.off()

pdf("plots/05m-all_lymphoma_HARLEQUIN_1q32.1.pdf", height=5, width=7)
plot.counts(dds_lrt, "HARLEQUIN_1q32.1") +  
  stat_compare_means(comparisons = list(c("ABC", "GCB"),
                                        c("ABC", "Unclass"),
                                        c("ABC", "Endemic BL EBV-positive"),
                                        c("ABC", "Endemic BL EBV-negative"),
                                        c("ABC", "Sporadic BL EBV-negative")),
                     method = "t.test", 
                     label = "p.signif", 
                     hide.ns = TRUE) 
  # Weird one
dev.off()

pdf("plots/05m-all_lymphoma_LTR46_Xq11.1.pdf", height=4, width=6)
plot.counts(dds_lrt, "LTR46_Xq11.1") # Upregulated in FL
dev.off()

pdf("plots/05m-all_lymphoma_upreg_DZ.pdf", height=4, width=18)
plot_grid(plot.counts(dds_lrt, "HERVH_3q22.1e"), # Upregulated in DZ GCB, 
          plot.counts(dds_lrt, "HML5_1q22"), # Upregulated in DZ GCB 
          plot.counts(dds_lrt, "HERVFH21_1p36.31"),  # Upregulated in DZ Agirre
          nrow = 1, 
          ncol = 3)
dev.off()

pdf("plots/05m-all_lymphoma_upreg_PB.pdf", height=8, width=12)
plot_grid(plot.counts(dds_lrt, "HERVH_15q26.3b"), # Upregulated in BMPC Agirre
          plot.counts(dds_lrt, "MER61_3q24a"), # Upregulated in PB Agirre
          plot.counts(dds_lrt, "PABLB_2q31.1"), # Upregulated in PB Agirre
          plot.counts(dds_lrt, "HERVP71A_15q24.2"),# Upregulated in PB Agirre
          nrow = 2, 
          ncol = 2)
dev.off()

pdf("plots/05m-all_lymphoma_BL_lasso_features.pdf", height=8, width=18)
plot_grid(plot.counts(dds_lrt, "ERVLE_2p25.3c") +
            stat_compare_means(comparisons = list(c("ABC", "Endemic BL EBV-positive"))),
          plot.counts(dds_lrt, "MER61_4p16.3") +
            stat_compare_means(comparisons = list(c("ABC", "Endemic BL EBV-positive"))), 
          plot.counts(dds_lrt, "ERV316A3_2q21.2b") +
            stat_compare_means(comparisons = list(c("ABC", "Endemic BL EBV-positive"))),
          plot.counts(dds_lrt, "ERVLE_5p13.2c") +
            stat_compare_means(comparisons = list(c("ABC", "Endemic BL EBV-positive"))),
          nrow = 2, 
          ncol = 2) 
dev.off()

#################### SELECTED FEATURES BROAD LYMPHOMA TYPE ##################### 

# > selected_vars$lasso
# [1] "ERVL_1p34.2"  "ERVLB4_2p16.3" "MER4B_10q21.3" "ERVL_Xq21.1b"  "ERVLE_14q23.2"

plot.counts <- function(df, gene) {
  
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="cancer_type", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=cancer_type, y=count, fill=cancer_type))  +
    geom_boxplot(notch = TRUE) +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("subtype") +
    ylab("Counts") +
    ggtitle(gene) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1) +
    scale_fill_manual(values = c("DLBCL" = "0073C2B2", 
                                 "BL" = "#EFC000B2",
                                 "FL" = "#868686B2")) + 
    scale_y_log10(labels = label_comma())   +
    stat_compare_means(comparisons = list(c("DLBCL", "BL"),
                                          c("BL", "FL"),
                                          c("DLBCL", "FL")),
                       method = "t.test", 
                       label = "p.signif", 
                       hide.ns = TRUE) +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          title=element_text(size=14, face="bold")
}


p1 <- plot.counts(dds_lrt, "ERVL_1p34.2")
p2 <- plot.counts(dds_lrt, "ERVLB4_2p16.3")
p3 <- plot.counts(dds_lrt, "MER4B_10q21.3")
p4 <- plot.counts(dds_lrt, "ERVL_Xq21.1b")
p5 <- plot.counts(dds_lrt, "ERVLE_14q23.2")

pdf("plots/05m-all_broad_lymphoma_selected_features.pdf", height=8, width=12)
plot_grid(p1, p2, p3, p4, p5,
          nrow = 2, 
          ncol = 3,
          labels = "AUTO")
dev.off()

pdf("plots/05m-all_broad_lymphoma_HARLEQUIN_1q32.1.pdf", height=4, width=4)
plot.counts(dds_lrt, "HARLEQUIN_1q32.1") # Weird one
dev.off()

pdf("plots/05m-all_broad_lymphoma_LTR46_Xq11.1.pdf", height=4, width=4)
plot.counts(dds_lrt, "LTR46_Xq11.1") # Upregulated in FL
dev.off()

pdf("plots/05m-all_broad_lymphoma_upreg_DZ.pdf", height=8, width=8)
plot_grid(plot.counts(dds_lrt, "MER61_3q13.11"),
          plot.counts(dds_lrt, "HML5_1q22"), 
          plot.counts(dds_lrt, "HERV3_14q32.33"), 
          plot.counts(dds_lrt, "HARLEQUIN_19p12b"),
          nrow = 2, 
          ncol = 2)
dev.off()

pdf("plots/05m-all_broad_lymphoma_upreg_PB.pdf", height=8, width=8)
plot_grid(plot.counts(dds_lrt, "HERVH_15q26.3b"), # Upregulated in BMPC Agirre
          plot.counts(dds_lrt, "MER61_3q24a"), # Upregulated in PB Agirre
          plot.counts(dds_lrt, "PABLB_2q31.1"), # Upregulated in PB Agirre
          plot.counts(dds_lrt, "HERVP71A_15q24.2"),# Upregulated in PB Agirre
          nrow = 2, 
          ncol = 2)
dev.off()

pdf("plots/05m-all_broad_lymphoma_BL_lasso_features.pdf", height=8, width=8)
plot_grid(plot.counts(dds_lrt, "ERVLE_2p25.3c"),
          plot.counts(dds_lrt, "MER61_4p16.3"),
          plot.counts(dds_lrt, "ERV316A3_2q21.2b"),
          plot.counts(dds_lrt, "ERVLE_5p13.2c"),
          nrow = 2, 
          ncol = 2) 
dev.off()

pdf("plots/05m-all_broad_lymphoma_upreg_FL.pdf", height=4, width=8)
plot_grid(plot.counts(dds_lrt, "HERVW_6q22.1"),
          plot.counts(dds_lrt, "HARLEQUIN_4q13.2a")
          )
dev.off()

#################################### SAVE DATA #################################

save(selected_vars, mat.sel, sig_lrt, bor.model, res.stabsel.c, 
     file="r_outputs/05m-all_lymphoma_selected_features.Rdata")

load("r_outputs/05m-all_lymphoma_selected_features.Rdata")


