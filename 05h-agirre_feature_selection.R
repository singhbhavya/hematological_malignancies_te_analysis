################################################################################
################################################################################
################################################################################
################################################################################
####################### BULK AGIRRE / GCB COMMON HERVS #########################


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

load("r_outputs/02-GCB_Agirre_filt_counts.Rdata")

################################# METADATA SETUP ###############################

agirre_metadata$Cell <- agirre_metadata$source_name
agirre_metadata$Cell[agirre_metadata$Cell == "Naive"] <- "Naive B"
agirre_metadata$Cell[agirre_metadata$Cell == "Centroblast"] <- "Dark Zone Germinal Center B"
agirre_metadata$Cell[agirre_metadata$Cell == "Centrocyte"] <- "Light Zone Germinal Center B"
agirre_metadata$Cell[agirre_metadata$Cell == "Memory"] <- "Memory B"
agirre_metadata$Cell[agirre_metadata$Cell == "Tonsilar plasma cell"] <- "Plasma cells"
agirre_metadata$Cell[agirre_metadata$Cell == "Bone Marrow plasma cell"] <- "Plasma cells"

################################### DATA SETUP #################################

samples <- agirre_metadata[agirre_metadata$Cell != "Bone Marrow plasma cell",]
counts <- GCB_Agirre.filt.herv[,rownames(samples)]

#################################### DESEQ2 ####################################

countDat <- GCB_Agirre.filt.herv
cat(sprintf('%d variables\n', nrow(counts)))

stopifnot(all(colnames(countDat) == rownames(agirre_metadata)))

dds <- DESeq2::DESeqDataSetFromMatrix(countDat, agirre_metadata, ~1 )
dds <- DESeq2::DESeq(dds, parallel=T)
tform <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)

##################################### LRT ######################################

p.cutoff <- 1e-3
lfc.cutoff <- 1.5
dds_lrt <- DESeq2::DESeqDataSetFromMatrix(countDat, agirre_metadata, ~ Cell)
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

agirre_metadata$Cell <- as.factor(agirre_metadata$Cell)
set.seed(12345)
bor.orig <- Boruta(x=t(mat), y=agirre_metadata$Cell, doTrace=2, ntree=1000, maxRuns=1000)
print(bor.orig)
bor.model <- TentativeRoughFix(bor.orig)
print(bor.model)

selected_vars$boruta <- names(bor.model$finalDecision)[bor.model$finalDecision == 'Confirmed']

#################################### LASSO #####################################

weakness <- 0.6    # For each subsample the features are reweighted by a random weight uniformly sampled in [weakness,1]
size <- 1.0        # proportion of samples drawn in every subsample
steps <- 100       # number of subsamples

error <- 0.01  # the desired type I error level w.r.t. to the chosen type I error rate.
pi_thr <- 0.6      # the threshold used for the stability selection

ncores <- 1        # force run on 1 core, prevents reseeding RNG on each core

set.seed(12345)

res.stabpath.c <-
  c060::stabpath(
    agirre_metadata$Cell,
    t(mat),
    size = size,
    steps = steps,
    weakness = weakness,
    mc.cores = ncores,
    family = "multinomial",
    type.multinomial = "grouped"
    )

res.stabsel.c <-
  c060::stabsel(
    res.stabpath.c,
    error = error,
    pi_thr = pi_thr,
    type = "pfer"
    )

selected_vars$lasso.clust <- names(res.stabsel.c$stable)
plot(res.stabpath.c)

##################################### UPSET ####################################

pdf("plots/05h-agirre_features_upset.pdf", height=5.5, width=7.5)
upset(fromList(selected_vars), 
      sets=c('lrt', 'boruta', 'lasso.clust'), 
      order.by = "freq",
      text.scale = c(2, 2, 2, 2, 2, 2))
dev.off()

################################# FEATURE GRAPH ################################

mat.sel <- mat[selected_vars$lasso.clust, ]
rpart.fit <- rpart(y ~ . , data=data.frame(y=agirre_metadata$Cell, 
                                           t(mat.sel)),
                   method = 'class')

pdf("plots/05h-agirre_features_rpart.pdf", height=7, width=9)
rpart.plot(rpart.fit, type=5, extra=1, digits=3)
dev.off()

############################### SELECTED FEATURES ##############################

plot.counts <- function(df, gene) {
  
  as.data.frame(plotCounts(df, 
                           gene=gene, 
                           intgroup="Cell", 
                           returnData = TRUE)) %>%
    ggplot(aes(x=Cell, y=count, fill=Cell))  +
    geom_boxplot(notch = TRUE) +
    theme_pubr() +
    theme(legend.position="none",
          axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("Cell") +
    ylab("Counts") +
    scale_fill_manual(values = c("Naive B" = pal_jco("default", alpha = 0.8)(7)[1],
                                 "Memory B" = pal_jco("default", alpha = 0.8)(7)[3],
                                 "Dark Zone Germinal Center B" = pal_jco("default", alpha = 0.8)(7)[4],
                                 "Light Zone Germinal Center B" = pal_jco("default", alpha = 0.8)(7)[5],
                                 "Plasma cells" = pal_jco("default", alpha = 0.8)(7)[6])) + 
    scale_x_discrete(labels=c("Naive B" = "NB", 
                              "Memory B" = "MB",
                              "Dark Zone Germinal Center B" = "DZ",
                              "Light Zone Germinal Center B" = "LZ",
                              "Plasma cells" = "PB")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    ggtitle(gene) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1) + 
    scale_y_log10(labels = label_comma())
}


p1 <- plot.counts(dds_lrt, "ERVLB4_14q23.3") + 
  stat_compare_means(comparisons = list(c("Memory B", "Plasma cells")),
                     method = "t.test", 
                     label = "p.signif", 
                     hide.ns = TRUE)
p2 <- plot.counts(dds_lrt, "HERVL_2p12a") +
  stat_compare_means(comparisons = list(c("Memory B", "Naive B")),
                     method = "t.test", 
                     label = "p.signif", 
                     hide.ns = TRUE)
p3 <- plot.counts(dds_lrt, "HERVP71A_8q24.13") +
  stat_compare_means(comparisons = list(c("Memory B", "Plasma cells")),
                   method = "t.test", 
                   label = "p.signif", 
                   hide.ns = TRUE)
p4 <- plot.counts(dds_lrt, "MER61_19p12c") +
  stat_compare_means(comparisons = list(c("Dark Zone Germinal Center B", 
                                          "Memory B"),
                                        c("Dark Zone Germinal Center B",
                                          "Naive B")),
                     method = "t.test", 
                     label = "p.signif", 
                     hide.ns = TRUE)
p5 <- plot.counts(dds_lrt, "HARLEQUIN_19p12b") +
  stat_compare_means(comparisons = list(c("Dark Zone Germinal Center B", 
                                          "Light Zone Germinal Center B")),
                     method = "t.test", 
                     label = "p.signif", 
                     hide.ns = TRUE)
p6 <-  plot.counts(dds_lrt, "HERVFRD_2p12a") +
  stat_compare_means(comparisons = list(c("Memory B", "Naive B")),
                     method = "t.test", 
                     label = "p.signif", 
                     hide.ns = TRUE)
p7 <-  plot.counts(dds_lrt, "PABLB_7q11.21") +
  stat_compare_means(comparisons = list(c("Dark Zone Germinal Center B", 
                                          "Light Zone Germinal Center B")),
                     method = "t.test", 
                     label = "p.signif", 
                     hide.ns = TRUE)
p8 <-  plot.counts(dds_lrt, "HERVL_1q23.3a") +
  stat_compare_means(comparisons = list(c("Dark Zone Germinal Center B", 
                                          "Plasma cells")),
                     method = "t.test", 
                     label = "p.signif", 
                     hide.ns = TRUE)
p9 <-  plot.counts(dds_lrt, "HERVP71A_15q24.2") +
  stat_compare_means(comparisons = list(c("Light Zone Germinal Center B", 
                                          "Plasma cells"),
                                        c("Memory B",
                                          "Plasma cells"),
                                        c("Plasma cells",
                                          "Naive B")),
                     method = "t.test", 
                     label = "p.signif", 
                     hide.ns = TRUE)
p10 <-  plot.counts(dds_lrt, "HUERSP2_6p22.3") +
  stat_compare_means(comparisons = list(c("Light Zone Germinal Center B", 
                                          "Naive B"),
                                        c("Memory B",
                                          "Naive B"),
                                        c("Plasma cells",
                                          "Naive B")),
                     method = "t.test", 
                     label = "p.signif", 
                     hide.ns = TRUE)
p11 <- plot.counts(dds_lrt, "ERVLE_6p25.1b") +
  stat_compare_means(comparisons = list(c("Dark Zone Germinal Center B", 
                                          "Light Zone Germinal Center B"),
                                        c("Dark Zone Germinal Center B",
                                          "Memory B")),
                     method = "t.test", 
                     label = "p.signif", 
                     hide.ns = TRUE)

pdf("plots/05h-agirre_selected_features_lasso.pdf", height=10, width=12)
plot_grid(p1, p2, p3, p4, p5, p6, p7,
          p8, p9, p10, p11, 
          nrow = 3, 
          ncol = 4,
          labels = NULL)
dev.off()

#################################### SAVE DATA #################################

save(selected_vars, mat.sel, sig_lrt, bor.model, res.stabsel.c, 
     file="r_outputs/05h-agirre_selected_features.Rdata")

load("r_outputs/05h-agirre_selected_features.Rdata")




