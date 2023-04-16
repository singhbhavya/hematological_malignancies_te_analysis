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
load("r_outputs/01-refs.Rdata")

################################################################################
################################################################################
#################################### CLUST 2 ###################################
################################################################################
################################################################################

################################# METADATA SETUP ###############################

BL_metadata$clust.retro.k2 <- clust.df$clust.retro.k2

#################################### DESEQ2 ####################################

hervs.to.keep <- intersect(rownames(BL.filt.herv), 
                           retro.annot$locus[retro.annot$chrom != "chrY"])

BL.filt.herv.y <- BL.filt.herv[hervs.to.keep,] 

countDat <- BL.filt.herv.y
cat(sprintf('%d variables\n', nrow(countDat)))

stopifnot(all(colnames(countDat) == rownames(BL_metadata)))

dds <- DESeq2::DESeqDataSetFromMatrix(countDat, BL_metadata, ~1 )
dds <- DESeq2::DESeq(dds, parallel=T)
tform <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)

##################################### LRT ######################################

p.cutoff <- 1e-3
lfc.cutoff <- 1.5
dds_lrt <- DESeq2::DESeqDataSetFromMatrix(countDat, BL_metadata, ~ clust.retro.k2)
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

BL_metadata$clust.retro.k2 <- as.factor(BL_metadata$clust.retro.k2)
set.seed(12345)
bor.orig <- Boruta(x=t(mat), y=BL_metadata$clust.retro.k2, doTrace=2, ntree=1000, maxRuns=1000)
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
    BL_metadata$clust.retro.k2,
    t(mat),
    mc.cores = ncores,
    weakness=weakness,
    family = "binomial"
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

mat.sel <- mat[Reduce(intersect, 
                      list(selected_vars$lrt, 
                           selected_vars$boruta, 
                           selected_vars$lasso)), ]

rpart.fit <- rpart(y ~ . , data=data.frame(y=BL_metadata$clust.retro.k2, 
                                           t(mat.sel)),
                   method = 'class')

pdf("plots/05j-BL_features_rpart.pdf", height=7, width=9)
rpart.plot(rpart.fit, type=5, extra=1, digits=3)
dev.off()

############################### SELECTED FEATURES ##############################

# > Reduce(intersect, list(selected_vars$lrt, selected_vars$boruta, selected_vars$lasso))
# [1] "HML6_Yp11.2" "HERV3_Yp11.2a" "HERVL_Yp11.2a" "MER4_Yq11.221b""HARLEQUIN_Yq11.221c"
# [6] "MER4_Yq11.221e"     

# > selected_vars$lasso
# [1] "ERVLE_2p25.3c"    "MER61_4p16.3"     "ERV316A3_2q21.2b" "ERVLE_5p13.2c"  

plot.counts <- function(df, gene) {
  
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="clust.retro.k2", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=clust.retro.k2, y=count, fill=clust.retro.k2))  +
    geom_boxplot() +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("clust.retro.k2") +
    ylab("Counts") +
    scale_fill_manual(values = c("C1" = wes_palette("Chevalier1")[1], 
                                 "C2" = wes_palette("Chevalier1")[2])) + 
    ggtitle(gene) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)  + 
    scale_y_log10(labels = label_comma()) 
}


p1 <- plot.counts(dds, "ERVLE_2p25.3c")
p2 <- plot.counts(dds, "MER61_4p16.3")
p3 <- plot.counts(dds, "ERV316A3_2q21.2b")
p4 <- plot.counts(dds, "ERVLE_5p13.2c")

pdf("plots/05j-BL_lasso_selected_features.pdf", height=6, width=9)
plot_grid(p1, p2, p3, p4, 
          nrow = 2, 
          ncol = 2)
dev.off()

################################# CLUSTER SIZES ################################

reads.left <- as.data.frame(colSums(BL.counts.mfilt.comb))
colnames(reads.left) <- c("reads")
reads.left$sample <- rownames(reads.left)

reads.left$clust <- BL_metadata$clust.retro.k2

reads.left %>%
  group_by(clust) %>%
  summarise_at(vars(reads), list(name = mean))

# # A tibble: 3 × 2
# clust      name
#<fct>     <dbl>
# C1    46688579.
# C2    43442928


