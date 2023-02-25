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

upset(fromList(selected_vars), sets=c('lrt', 'boruta', 'lasso.clust'), order.by = "freq")

################################# FEATURE GRAPH ################################

mat.sel <- mat[selected_vars$lasso.clust, ]
rpart.fit <- rpart(y ~ . , data=data.frame(y=agirre_metadata$Cell, 
                                           t(mat.sel)),
                   method = 'class')

rpart.plot(rpart.fit, type=5, extra=1, digits=3)
