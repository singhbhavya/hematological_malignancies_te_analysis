################################################################################
################################################################################
################################################################################
################################################################################
################################# LYMPHOMA DE ##################################

## Plan:
## DLBCL top 20 and top 100 HERVs/genes heatmaps + features of interest
## BL top 20 and top 100 HERVs/genes heatmaps + features of interest
## FL top 20 and top 100 HERVs/genes heatmaps + features of interest
## All top 20 and top 100 HERVs/genes heatmaps + features of interest

#################################### SETUP #####################################

library(knitr)
library(tidyverse)
library(matrixStats)
library(data.table)
library(PCAtools)
library(DESeq2)
library(ggplot2)
library(ggsci)
library(edgeR)
library(ashr)
library(cowplot)
library(wesanderson)
library(pheatmap)

################################### LOAD DATA ##################################

load("r_outputs/01-metadata.Rdata")
load("r_outputs/01-refs.Rdata")
load("r_outputs/05b-DLBCL_pca_dds.Rdata")
load("r_outputs/05b-BL_pca_dds.Rdata")
load("r_outputs/05b-FL_pca_dds.Rdata")
load("r_outputs/05b-all_lymphoma_pca_dds.Rdata")

################################ SET THRESHOLDS ################################

fc=2 # fold change
l2fc=log2(fc) # log 2 fold change
pval=0.05 # p value threshold

########################## EXPLORE DESEQ DLBCL HERVs ###########################

dlbcl_res <- DESeq2::results(DLBCL.dds, name="COO_class_GCB_vs_ABC")
summary(dlbcl_res)
# Order by p value
dlbcl.res.ordered <- dlbcl_res[order(dlbcl_res$padj),]
# Filter by p value
dlbcl.res05 <- DESeq2::results(DLBCL.dds, alpha=pval)
summary(dlbcl.res05)
# How many HERVs left that are significant?
sum(dlbcl.res05$padj < 0.05, na.rm=TRUE)
# How many HERVs left after log2foldchange + padj?
sum(dlbcl.res05$padj < 0.05, na.rm =TRUE, lfcThreshold=2)

# Filter genes
dlbcl.res05lf2 <- DESeq2::results(DLBCL.dds, 
                                  independentFiltering = TRUE) %>% 
  as.data.frame(.) %>% tibble::rownames_to_column() %>% 
  as_tibble(., rowname = NA) %>%
  dplyr::select(everything()) %>%
  arrange(padj) %>%
  dplyr::filter(log2FoldChange >= l2fc | log2FoldChange <= -l2fc)

dlbcl.res05lf2[1:20,]

dlbcl.res05lf2 <-
  dlbcl.res05lf2 %>% 
  remove_rownames %>% tibble::column_to_rownames(var="rowname")

######################## DLBCL SPECIFIC HERV FEATURES ##########################

# Looking at some differentially expressed HERVs with lowest p value
pdf("plots/05c-dlbcl_HARLEQUIN_17q25.3a.pdf", height=4, width=4)
plotCounts(DLBCL.dds, gene=which.min(dlbcl_res$padj), intgroup="COO_class")
dev.off()
pdf("plots/05c-dlbcl_MER4_17q21.2d.pdf", height=4, width=4)
plotCounts(DLBCL.dds, gene="MER4_17q21.2d", intgroup="COO_class")
dev.off()
pdf("plots/05c-dlbcl_HERVE_19q12.pdf", height=4, width=4)
plotCounts(DLBCL.dds, gene="HERVE_19q12", intgroup="COO_class")
dev.off()
pdf("plots/05c-dlbcl_HERVH_17q25.3.pdf", height=4, width=4)
plotCounts(DLBCL.dds, gene="HERVH_17q25.3", intgroup="COO_class")
dev.off()
pdf("plots/05c-dlbcl_HARLEQUIN_1q32.1.pdf", height=4, width=4)
plotCounts(DLBCL.dds, gene="HARLEQUIN_1q32.1", intgroup="COO_class")
dev.off()
pdf("plots/05c-dlbcl_ERV316A3_9q21.13a.pdf", height=4, width=4)
plotCounts(DLBCL.dds, gene="ERV316A3_9q21.13a", intgroup="COO_class")
dev.off()
pdf("plots/05c-dlbcl_MER4_17q21.2c.pdf", height=4, width=4)
plotCounts(DLBCL.dds, gene="MER4_17q21.2c", intgroup="COO_class")
dev.off()
pdf("plots/05c-dlbcl_HARLEQUIN_4p15.2.pdf", height=4, width=4)
plotCounts(DLBCL.dds, gene="HARLEQUIN_4p15.2", intgroup="COO_class")
dev.off()
pdf("plots/05c-dlbcl_HML6_19q13.43b.pdf", height=4, width=4)
plotCounts(DLBCL.dds, gene="HML6_19q13.43b", intgroup="COO_class")
dev.off()
pdf("plots/05c-dlbcl_ERV316A3_8q13.3a.pdf", height=4, width=4)
plotCounts(DLBCL.dds, gene="ERV316A3_8q13.3a", intgroup="COO_class")
dev.off()

############################ DLBCL HERV HEATMAPS ###############################

select <- rownames(dlbcl.res05lf2[1:20,])
df <- as.data.frame(colData(DLBCL.dds)[,c("COO_class","project")])

cols <- rgb_gsea(palette = c("default"), n = 12, alpha = 0.7, reverse = FALSE)

## Top 20 HERVs
pdf("plots/05c-dlbcl_deseq_coo_top20_hervs_filtered.pdf", height=6, width=15)
pheatmap(assay(DLBCL.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         color = cols,
         breaks=seq(-3,3,length.out=14),
         scale="row",
         cluster_cols=TRUE, annotation_col=df)

dev.off()

## Top 100 HERVs
select <- rownames(dlbcl.res05lf2[1:100,])
pdf("plots/05c-dlbcl_deseq_coo_top100_hervs_filtered.pdf", height=16, width=18)
pheatmap(assay(DLBCL.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         scale="row",
         color = cols,
         breaks=seq(-3,3,length.out=14),
         cluster_cols=TRUE, annotation_col=df)

dev.off()

###################### EXPLORE DESEQ DLBCL HERVs & GENES #######################

dlbcl_res_all <- DESeq2::results(DLBCL.g.dds, name="COO_class_GCB_vs_ABC")
rownames(dlbcl_res_all) <- gene_table[rownames(dlbcl_res_all), 'display']
summary(dlbcl_res_all)
# Order by p value
dlbcl.res.ordered.all <- dlbcl_res_all[order(dlbcl_res_all$padj),]
dlbcl.res.ordered.all[1:20,]
# Filter by p value
dlbcl.res05.all <- DESeq2::results(DLBCL.g.dds, alpha=pval)
summary(dlbcl.res05.all)
# How many HERVs left that are significant?
sum(dlbcl.res05.all$padj < 0.05, na.rm=TRUE)
# How many HERVs left after log2foldchange + padj?
sum(dlbcl.res05.all$padj < 0.05, na.rm =TRUE, lfcThreshold=2)

# Filter genes
dlbcl.res05lf2.all <- DESeq2::results(DLBCL.g.dds, 
                                      independentFiltering = TRUE) %>% 
  as.data.frame(.) %>% tibble::rownames_to_column() %>% 
  as_tibble(., rowname = NA) %>%
  dplyr::select(everything()) %>%
  arrange(padj) %>%
  dplyr::filter(log2FoldChange >= l2fc | log2FoldChange <= -l2fc)

dlbcl.res05lf2.all <-
  dlbcl.res05lf2.all %>% 
  remove_rownames %>% tibble::column_to_rownames(var="rowname")


######################## DLBCL HERV & GENE HEATMAPS ############################

## Top 20 genes & HERVs
select <- rownames(dlbcl.res05lf2.all[1:20,])
df <- as.data.frame(colData(DLBCL.g.dds)[,c("COO_class","project")])

cols <- rgb_gsea(palette = c("default"), n = 12, alpha = 0.7, reverse = FALSE)

pdf("plots/05c-dlbcl_deseq_coo_top20_hervsgenes_filtered.pdf", height=6, width=15)
pheatmap(assay(DLBCL.g.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         scale="row",
         color = cols,
         breaks=seq(-3,3,length.out=14),
         labels_row = gene_table[rownames(dlbcl.res05lf2.all)[1:20], 'display'],
         cluster_cols=TRUE, annotation_col=df)

dev.off()

## Top 100 genes & HERVs

select <- rownames(dlbcl.res05lf2.all[1:100,])
pdf("plots/05c-dlbcl_deseq_coo_top100_hervsgenes_filtered.pdf", height=16, width=20)
pheatmap(assay(DLBCL.g.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         breaks=seq(-3,3,length.out=14),
         scale="row",
         color = cols,
         labels_row = gene_table[rownames(dlbcl.res05lf2.all)[1:100], 'display'],
         cluster_cols=TRUE, annotation_col=df)
dev.off()

########################### EXPLORE DESEQ BL HERVs #############################

bl_res <- DESeq2::results(BL.dds, 
                          name="ebv_status_EBV.positive_vs_EBV.negative")
summary(bl_res)
# Order by p value
bl.res.ordered <- bl_res[order(bl_res$padj),]
# Filter by p value
bl.res05 <- DESeq2::results(BL.dds, alpha=pval)
summary(bl.res05)
# How many HERVs left that are significant?
sum(bl.res05$padj < 0.05, na.rm=TRUE)
# How many HERVs left after log2foldchange + padj?
sum(bl.res05$padj < 0.05, na.rm =TRUE, lfcThreshold=2)

# Filter genes
bl.res05lf2 <- DESeq2::results(BL.dds, 
                               independentFiltering = TRUE) %>% 
  as.data.frame(.) %>% tibble::rownames_to_column() %>% 
  as_tibble(., rowname = NA) %>%
  dplyr::select(everything()) %>%
  arrange(padj) %>%
  dplyr::filter(log2FoldChange >= l2fc | log2FoldChange <= -l2fc)

bl.res05lf2[1:20,]

bl.res05lf2 <-
  bl.res05lf2 %>% 
  remove_rownames %>% tibble::column_to_rownames(var="rowname")

########################## BL SPECIFIC HERV FEATURES ###########################

# Looking at some differentially expressed HERVs with lowest p value
pdf("plots/05c-BL_HERVH_14q22.1a.pdf", height=4, width=4)
HERVH_14q22.1a <- plotCounts(BL.dds, gene="HERVH_14q22.1a", 
                             intgroup=c("ebv_status", "clinical_variant"),
                             returnData=TRUE)

HERVH_14q22.1a$clinical_variant_status <- paste(HERVH_14q22.1a$ebv_status,
                                                HERVH_14q22.1a$clinical_variant,
                                                sep=" ")

ggplot(HERVH_14q22.1a, aes(x = clinical_variant_status, 
                           y = count))  + 
  geom_bar(stat="identity") +
  ggtitle("HERVH_14q22.1a") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75))

dev.off()
 
pdf("plots/05c-BL_HERVH_11q13.4a.pdf", height=4, width=4)
HERVH_11q13.4a <-
  plotCounts(BL.dds, gene="HERVH_11q13.4a", 
             intgroup=c("ebv_status", "clinical_variant"),
             returnData = TRUE)

HERVH_11q13.4a$clinical_variant_status <- paste(HERVH_11q13.4a$ebv_status,
                                                HERVH_11q13.4a$clinical_variant,
                                                sep=" ")

ggplot(HERVH_11q13.4a, aes(x = clinical_variant_status, 
                           y = count))  + 
  geom_bar(stat="identity") +
  ggtitle("HERVH_11q13.4a") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75))
dev.off()

pdf("plots/05c-BL_HERVH_3q26.33a.pdf", height=4, width=4)

HERVH_3q26.33a <-
  plotCounts(BL.dds, gene="HERVH_3q26.33a", 
           intgroup=c("ebv_status", "clinical_variant"),
           returnData = TRUE)

HERVH_3q26.33a$clinical_variant_status <- paste(HERVH_3q26.33a$ebv_status,
                                                HERVH_3q26.33a$clinical_variant,
                                                sep=" ")

ggplot(HERVH_3q26.33a, aes(x = clinical_variant_status, 
                           y = count))  + 
  geom_bar(stat="identity") +
  ggtitle("HERVH_3q26.33a") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75))
dev.off()



pdf("plots/05c-BL_ERV316A3_3q13.31c.pdf", height=4, width=4)
ERV316A3_3q13.31c <-
  plotCounts(BL.dds, gene="ERV316A3_3q13.31c", intgroup=c("ebv_status", "clinical_variant"),
           returnData = TRUE)

ERV316A3_3q13.31c$clinical_variant_status <- paste(ERV316A3_3q13.31c$ebv_status,
                                                   ERV316A3_3q13.31c$clinical_variant,
                                                sep=" ")

ggplot(ERV316A3_3q13.31c, aes(x = clinical_variant_status, 
                           y = count))  + 
  geom_bar(stat="identity") +
  ggtitle("ERV316A3_3q13.31c") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75))


dev.off()

pdf("plots/05c-BL_HARLEQUIN_1q32.1.pdf", height=4, width=4)
HARLEQUIN_1q32.1 <-
  plotCounts(BL.dds, gene="HARLEQUIN_1q32.1", intgroup=c("ebv_status", "clinical_variant"),
             returnData = TRUE)

HARLEQUIN_1q32.1$clinical_variant_status <- paste(HARLEQUIN_1q32.1$ebv_status,
                                                  HARLEQUIN_1q32.1$clinical_variant,
                                                   sep=" ")

ggplot(HARLEQUIN_1q32.1, aes(x = clinical_variant_status, 
                              y = count))  + 
  geom_bar(stat="identity") +
  ggtitle("HARLEQUIN_1q32.1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75))


dev.off()


############################## BL HERV HEATMAPS ################################

select <- rownames(bl.res05lf2[1:20,])
df <- as.data.frame(colData(BL.dds)[,c("clinical_variant","ebv_status", "gender")])

cols <- rgb_gsea(palette = c("default"), n = 12, alpha = 0.7, reverse = FALSE)

## Top 20 HERVs
pdf("plots/05c-BL_deseq_ebvclinvar_top20_hervs_filtered.pdf", height=4, width=10)
pheatmap(assay(BL.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         color = cols,
         breaks=seq(-3,3,length.out=14),
         scale="row",
         cluster_cols=TRUE, annotation_col=df)

dev.off()

## Top 100 HERVs
select <- rownames(bl.res05lf2[1:100,])
pdf("plots/05c-BL_deseq_ebvclinvar_top100_hervs_filtered.pdf", height=16, width=10)

pheatmap(assay(BL.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         scale="row",
         color = cols,
         breaks=seq(-3,3,length.out=14),
         cluster_cols=TRUE, annotation_col=df)
dev.off()

######################## EXPLORE DESEQ BL HERVs & GENES ########################

bl_res_all <- DESeq2::results(BL.g.dds, 
                              name="ebv_status_EBV.positive_vs_EBV.negative")
rownames(bl_res_all) <- gene_table[rownames(bl_res_all), 'display']
summary(bl_res_all)
# Order by p value
bl.res.ordered.all <- bl_res_all[order(bl_res_all$padj),]
bl.res.ordered.all[1:20,]
# Filter by p value
bl.res05.all <- DESeq2::results(BL.g.dds, alpha=pval)
summary(bl.res05.all)
# How many HERVs left that are significant?
sum(bl.res05.all$padj < 0.05, na.rm=TRUE)
# How many HERVs left after log2foldchange + padj?
sum(bl.res05.all$padj < 0.05, na.rm =TRUE, lfcThreshold=2)

# Filter genes
bl.res05lf2.all <- DESeq2::results(BL.g.dds, 
                                   independentFiltering = TRUE) %>% 
  as.data.frame(.) %>% tibble::rownames_to_column() %>% 
  as_tibble(., rowname = NA) %>%
  dplyr::select(everything()) %>%
  arrange(padj) %>%
  dplyr::filter(log2FoldChange >= l2fc | log2FoldChange <= -l2fc)

bl.res05lf2.all <-
  bl.res05lf2.all %>% 
  remove_rownames %>% tibble::column_to_rownames(var="rowname")

########################## BL HERV & GENE HEATMAPS #############################

## Top 20 genes & HERVs
select <- rownames(bl.res05lf2.all[1:20,])
df <- as.data.frame(colData(BL.g.dds)[,c("clinical_variant","ebv_status", "gender")])

cols <- rgb_gsea(palette = c("default"), n = 12, alpha = 0.7, reverse = FALSE)

pdf("plots/05c-BL_deseq_ebvclinvar_top20_hervsgenes_filtered.pdf", height=4, width=10)
pheatmap(assay(BL.g.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         scale="row",
         color = cols,
         breaks=seq(-3,3,length.out=14),
         labels_row = gene_table[rownames(bl.res05lf2.all)[1:20], 'display'],
         cluster_cols=TRUE, annotation_col=df)

dev.off()

## Top 100 genes & HERVs

select <- rownames(bl.res05lf2.all[1:100,])
pdf("plots/05c-BL_deseq_ebvclinvar_top100_hervsgenes_filtered.pdf", height=16, width=10)
pheatmap(assay(BL.g.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         breaks=seq(-3,3,length.out=14),
         scale="row",
         color = cols,
         labels_row = gene_table[rownames(bl.res05lf2.all)[1:100], 'display'],
         cluster_cols=TRUE, annotation_col=df)
dev.off()

########################### EXPLORE DESEQ FL HERVs #############################

fl_res <- DESeq2::results(FL.dds, 
                          name="who_diagnosis_FOLLICULAR.GRADE.2_vs_FOLLICULAR.GRADE.1")
summary(fl_res)
# Order by p value
fl.res.ordered <- fl_res[order(fl_res$padj),]
# Filter by p value
fl.res05 <- DESeq2::results(FL.dds, alpha=pval)
summary(fl.res05)
# How many HERVs left that are significant?
sum(fl.res05$padj < 0.05, na.rm=TRUE)
# How many HERVs left after log2foldchange + padj?
sum(fl.res05$padj < 0.05, na.rm =TRUE, lfcThreshold=2)

# Filter genes
fl.res05lf2 <- DESeq2::results(FL.dds, 
                               independentFiltering = TRUE) %>% 
  as.data.frame(.) %>% tibble::rownames_to_column() %>% 
  as_tibble(., rowname = NA) %>%
  dplyr::select(everything()) %>%
  arrange(padj) %>%
  dplyr::filter(log2FoldChange >= l2fc | log2FoldChange <= -l2fc)

fl.res05lf2[1:20,]

fl.res05lf2 <-
  fl.res05lf2 %>% 
  remove_rownames %>% tibble::column_to_rownames(var="rowname")

########################## FL SPECIFIC HERV FEATURES ###########################

# Looking at some differentially expressed HERVs with lowest p value
pdf("plots/05c-fl_MER4_Yq11.221b.pdf", height=4, width=4)
plotCounts(FL.dds, gene="MER4_Yq11.221b", intgroup="who_diagnosis")
dev.off()

pdf("plots/05c-fl_HERVIP10FH_Xq24.pdf", height=4, width=4)
plotCounts(FL.dds, gene="HERVIP10FH_Xq24", intgroup="who_diagnosis")
dev.off()

############################## FL HERV HEATMAPS ################################

select <- rownames(fl.res05lf2[1:20,])
df <- as.data.frame(colData(FL.dds)[,c("who_diagnosis", "stage")])

cols <- rgb_gsea(palette = c("default"), n = 12, alpha = 0.7, reverse = FALSE)

## Top 20 HERVs
pdf("plots/05c-FL_deseq_diagnosis_top20_hervs_filtered.pdf", height=4, width=6)
pheatmap(assay(FL.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         color = cols,
         breaks=seq(-3,3,length.out=14),
         scale="row",
         cluster_cols=TRUE, annotation_col=df)

dev.off()

## Top 100 genes & HERVs

select <- rownames(fl.res05lf2[1:100,])
pdf("plots/05c-FL_deseq_diagnosis_top100_hervs_filtered.pdf", height=16, width=6)
pheatmap(assay(FL.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         breaks=seq(-3,3,length.out=14),
         scale="row",
         color = cols,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

######################## EXPLORE DESEQ FL HERVs & GENES ########################

fl_res_all <- DESeq2::results(FL.g.dds, 
                              name="who_diagnosis_FOLLICULAR.GRADE.2_vs_FOLLICULAR.GRADE.1")
rownames(fl_res_all) <- gene_table[rownames(fl_res_all), 'display']
summary(fl_res_all)
# Order by p value
fl.res.ordered.all <- fl_res_all[order(fl_res_all$padj),]
fl.res.ordered.all[1:20,]
# Filter by p value
fl.res05.all <- DESeq2::results(FL.g.dds, alpha=pval)
summary(fl.res05.all)
# How many HERVs left that are significant?
sum(fl.res05.all$padj < 0.05, na.rm=TRUE)
# How many HERVs left after log2foldchange + padj?
sum(fl.res05.all$padj < 0.05, na.rm =TRUE, lfcThreshold=2)

# Filter genes
fl.res05lf2.all <- DESeq2::results(FL.g.dds, 
                                   independentFiltering = TRUE) %>% 
  as.data.frame(.) %>% tibble::rownames_to_column() %>% 
  as_tibble(., rowname = NA) %>%
  dplyr::select(everything()) %>%
  arrange(padj) %>%
  dplyr::filter(log2FoldChange >= l2fc | log2FoldChange <= -l2fc)

fl.res05lf2.all <-
  fl.res05lf2.all %>% 
  remove_rownames %>% tibble::column_to_rownames(var="rowname")


########################## FL HERV & GENE HEATMAPS #############################

## Top 20 genes & HERVs
select <- rownames(fl.res05lf2.all[1:20,])
df <- as.data.frame(colData(FL.g.dds)[,c("who_diagnosis", "stage")])

cols <- rgb_gsea(palette = c("default"), n = 12, alpha = 0.7, reverse = FALSE)

pdf("plots/05c-FL_deseq_diagnosis_top20_hervsgenes_filtered.pdf", height=4, width=6)
pheatmap(assay(FL.g.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         scale="row",
         color = cols,
         breaks=seq(-3,3,length.out=14),
         labels_row = gene_table[rownames(fl.res05lf2.all)[1:20], 'display'],
         cluster_cols=TRUE, annotation_col=df)

dev.off()

## Top 100 genes & HERVs

select <- rownames(fl.res05lf2.all[1:100,])
pdf("plots/05c-FL_deseq_diagnosis_top100_hervsgenes_filtered.pdf", height=16, width=6)
pheatmap(assay(FL.g.tform)[select,], cluster_rows=TRUE,
         show_rownames=TRUE,
         show_colnames = FALSE,
         breaks=seq(-3,3,length.out=14),
         scale="row",
         color = cols,
         labels_row = gene_table[rownames(fl.res05lf2.all)[1:100], 'display'],
         cluster_cols=TRUE, annotation_col=df)
dev.off()
