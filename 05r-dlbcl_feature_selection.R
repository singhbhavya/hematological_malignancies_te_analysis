################################################################################
################################################################################
################################################################################
################################################################################
############################ DLBCL CLUSTER FEATURES ############################ 

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

load("r_outputs/02-DLBCL_filt_counts.Rdata")
load("r_outputs/05o-DLBCL_pca_ccp_clusters_metadata.Rdata")
load("r_outputs/01-refs.Rdata")

################################################################################
################################################################################
#################################### CLUST 7 ###################################
################################################################################
################################################################################

################################# METADATA SETUP ###############################

DLBCL_metadata$clust.retro.k7 <- clust.df$clust.retro.k7


#################################### DESEQ2 ####################################

countDat <- DLBCL.filt.herv
cat(sprintf('%d variables\n', nrow(countDat)))

stopifnot(all(colnames(countDat) == rownames(DLBCL_metadata)))

dds <- DESeq2::DESeqDataSetFromMatrix(countDat, DLBCL_metadata, ~1 )
dds <- DESeq2::DESeq(dds, parallel=T)
tform <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)

##################################### LRT ######################################

p.cutoff <- 1e-3
lfc.cutoff <- 1.5
dds_lrt <- DESeq2::DESeqDataSetFromMatrix(countDat,DLBCL_metadata, ~ clust.retro.k7)
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

DLBCL_metadata$clust.retro.k7 <- as.factor(DLBCL_metadata$clust.retro.k7)
set.seed(12345)
bor.orig <- Boruta(x=t(mat), y=DLBCL_metadata$clust.retro.k7, doTrace=2, ntree=1000, maxRuns=1000)
print(bor.orig)
bor.model <- TentativeRoughFix(bor.orig)
print(bor.model)

selected_vars$boruta <- names(bor.model$finalDecision)[bor.model$finalDecision == 'Confirmed']

#################################### LASSO #####################################

weakness <- 1    # For each subsample the features are reweighted by a random weight uniformly sampled in [weakness,1]
size <- 1      # proportion of samples drawn in every subsample
steps <- 100       # number of subsamples

error <- 0.01  # the desired type I error level w.r.t. to the chosen type I error rate.
pi_thr <- 0.6      # the threshold used for the stability selection

ncores <- 1        # force run on 1 core, prevents reseeding RNG on each core

set.seed(12345)

res.stabpath.c <-
  c060::stabpath(
    DLBCL_metadata$clust.retro.k7,
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

pdf("plots/05r-DLBCL_features_upset.pdf", height=5, width=7)
UpSetR::upset(fromList(selected_vars), 
              sets=c('lrt', 'boruta', 'lasso'), 
              order.by = "freq",
              text.scale = c(2, 2, 2, 2, 2, 2))
dev.off()

################################# FEATURE GRAPH ################################

mat.sel <- mat[Reduce(intersect, list(selected_vars$lrt, selected_vars$boruta, selected_vars$lasso)), ]
mat.sel <- mat[selected_vars$lasso, ]
# "HERVW_2q32.3"   "HML2_7p22.1"    "HERVH_7q11.23a" "HERVH_16p13.2e"

rpart.fit <- rpart(y ~ . , data=data.frame(y=DLBCL_metadata$clust.retro.k7, 
                                           t(mat.sel)),
                   method = 'class')

pdf("plots/05r-DLBCL_features_rpart.pdf", height=7, width=9)
rpart.plot(rpart.fit, type=5, extra=1, digits=3)
dev.off()

################################# PLOT COUNTS ################################## 


plot.counts <- function(df, gene, title) {
  
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="clust.retro.k7", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=clust.retro.k7, y=count, fill=clust.retro.k7))  +
    geom_boxplot(notch = TRUE) +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("HERV Cluster") +
    ylab("Counts") +
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1) +
    ylim(0,10000) +
    scale_fill_manual(values = c("HC1" = "#E64B35FF",
                                 "HC2" = "#4DBBD5FF",
                                 "HC3" = "#00A087FF",
                                 "HC4" = "#3C5488FF",
                                 "HC5" = "#F39B7FFF",
                                 "HC6" = "#8491B4FF",
                                 "HC7" = "#91D1C2FF")) + 
    scale_y_log10(labels = label_comma()) 
}


p1 <- plot.counts(dds_lrt, "HML2_7p22.1", "HML2_7p22.1")+
  stat_compare_means(comparisons = list(c("HC4", "HC6"),
                                        c("HC4", "HC5"),
                                        c("HC4", "HC7")),
                     method = "t.test", 
                     label = "p.signif")

p2 <- plot.counts(dds_lrt, "HERVH_16p13.2e", "HERVH_16p13.2e") +
  stat_compare_means(comparisons = list(c("HC7", "HC6"),
                                        c("HC7", "HC5"),
                                        c("HC7", "HC4"),
                                        c("HC7", "HC3"),
                                        c("HC7", "HC2"),
                                        c("HC7", "HC1")),
                     method = "t.test", 
                     label = "p.signif")

p3 <- plot.counts(dds_lrt, "HERVW_2q32.3", "HERVW_2q32.3") +
  stat_compare_means(comparisons = list(c("HC1", "HC2"),
                                        c("HC2", "HC3"),
                                        c("HC2", "HC7")),
                     method = "t.test", 
                     label = "p.signif")

p4 <- plot.counts(dds_lrt, "HERVH_7q11.23a", "HERVH_7q11.23a") +
  stat_compare_means(comparisons = list(c("HC3", "HC5"),
                                        c("HC5", "HC7"),
                                        c("HC2", "HC5")),
                     method = "t.test", 
                     label = "p.signif")

pdf("plots/05r-DLBCL_lasso_selected_features.pdf", height=6, width=9)
plot_grid(p1, 
          p2, 
          p3, 
          p4, 
          nrow = 2, 
          ncol = 2)
dev.off()

plot.counts(dds_lrt, "HARLEQUIN_1q32.1", "HARLEQUIN_1q32.1")+
  stat_compare_means(comparisons = list(c("C1", "C2"),
                                        c("C2", "C3"),
                                        c("C2", "C4"),
                                        c("C2", "C5"),
                                        c("C2", "C7"),
                                        c("C2", "C6")),
                     method = "t.test", 
                     label = "p.signif")

#################################### SAVE DATA #################################

save(selected_vars, mat.sel, dds_lrt, sig_lrt, bor.model, res.stabsel.c, 
     file="r_outputs/05r-dlbcl_selected_features.Rdata")

load("r_outputs/05r-dlbcl_selected_features.Rdata")

