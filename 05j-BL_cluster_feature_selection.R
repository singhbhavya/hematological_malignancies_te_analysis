################################################################################
################################################################################
################################################################################
################################################################################
####################### BURKITT LYMPHOMA CLUSTER FEATURES ######################

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

load("r_outputs/02-BL_filt_counts.Rdata")
load("r_outputs/05i-BL_pca_ccp_clusters_metadata.Rdata")


################################# METADATA SETUP ###############################

BL_metadata$clust.retro.k3 <- clust.df$clust.retro.k3

#################################### DESEQ2 ####################################

countDat <- BL.filt.herv
cat(sprintf('%d variables\n', nrow(counts)))

stopifnot(all(colnames(countDat) == rownames(BL_metadata)))

dds <- DESeq2::DESeqDataSetFromMatrix(countDat, BL_metadata, ~1 )
dds <- DESeq2::DESeq(dds, parallel=T)
tform <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)

##################################### LRT ######################################

p.cutoff <- 1e-3
lfc.cutoff <- 1.5
dds_lrt <- DESeq2::DESeqDataSetFromMatrix(countDat, BL_metadata, ~ clust.retro.k3)
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

BL_metadata$clust.retro.k3 <- as.factor(BL_metadata$clust.retro.k3)
set.seed(12345)
bor.orig <- Boruta(x=t(mat), y=BL_metadata$clust.retro.k3, doTrace=2, ntree=1000, maxRuns=1000)
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
    BL_metadata$clust.retro.k3,
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

pdf("plots/05j-BL_features_upset.pdf", height=6, width=8)
upset(fromList(selected_vars), 
      sets=c('lrt', 'boruta', 'lasso'), 
      order.by = "freq",
      text.scale = c(1.5, 1.5, 1.3, 1.3, 1.3, 1.3))
dev.off()

################################# FEATURE GRAPH ################################

mat.sel <- mat[selected_vars$lasso, ]
rpart.fit <- rpart(y ~ . , data=data.frame(y=BL_metadata$clust.retro.k3, 
                                           t(mat.sel)),
                   method = 'class')

pdf("plots/05j-BL_features_rpart.pdf", height=7, width=9)
rpart.plot(rpart.fit, type=5, extra=1, digits=3)
dev.off()

############################### SELECTED FEATURES ##############################

# > selected_vars$lasso
# [1] "MER101_16p12.2a" "HERVH48_8q24.13" "HERV3_14q32.33"  "HERVL40_4q31.3b"

plot.counts <- function(df, gene) {
  
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="clust.retro.k3", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=clust.retro.k3, y=count, fill=clust.retro.k3))  +
    geom_boxplot() +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("clust.retro.k3") +
    ylab("Counts") +
    scale_fill_manual(values = c("C1" = wes_palette("Chevalier1")[1], 
                                 "C2" = wes_palette("Chevalier1")[2],
                                 "C3" = wes_palette("Chevalier1")[3])) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    ggtitle(gene) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)
}


p1 <- plot.counts(dds_lrt, "HERVH48_8q24.13")
p2 <- plot.counts(dds_lrt, "MER101_16p12.2a")
p3 <- plot.counts(dds_lrt, "HERV3_14q32.33")
p4 <- plot.counts(dds_lrt, "HERVL40_4q31.3b")

plot_grid(p1, p2, p3, p4,
          nrow = 2, 
          ncol = 2,
          labels = "AUTO")

