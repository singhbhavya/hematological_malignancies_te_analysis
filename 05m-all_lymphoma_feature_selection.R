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


#################################### DESEQ2 ####################################

countDat <- all.counts.filt.herv
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

pdf("plots/05m-all_lymphoma_features_upset.pdf", height=6, width=8)
upset(fromList(selected_vars), 
      sets=c('lrt', 'boruta', 'lasso'), 
      order.by = "freq",
      text.scale = c(1.5, 1.5, 1.3, 1.3, 1.3, 1.3))
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
                           intgroup="cancer_type", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=cancer_type, y=count, fill=cancer_type))  +
    geom_boxplot() +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("cancer_type") +
    ylab("Counts") +
    ggtitle(gene) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)
}


p1 <- plot.counts(dds_lrt, "ERVL_1p34.2")
p2 <- plot.counts(dds_lrt, "ERVLB4_2p16.3")
p3 <- plot.counts(dds_lrt, "MER4B_10q21.3")
p4 <- plot.counts(dds_lrt, "ERVL_Xq21.1b")
p5 <- plot.counts(dds_lrt, "ERVLE_14q23.2")

pdf("plots/05m-all_lymphoma_selected_features.pdf", height=6, width=9)
plot_grid(p1, p2, p3, p4, p5,
          nrow = 2, 
          ncol = 3,
          labels = "AUTO")
dev.off()


#################################### SAVE DATA #################################

save(selected_vars, mat.sel, sig_lrt, bor.model, res.stabsel.c, 
     file="r_outputs/05m-all_lymphoma_selected_features.Rdata")

load("r_outputs/05m-all_lymphoma_selected_features.Rdata")


