################################################################################
################################################################################
################################################################################
################################################################################
####################### BULK AGIRRE / GCB COMMON HERVS #########################


#################################### SETUP #####################################

library(knitr)
library(tidyverse)
library(matrixStats)
library(data.table)
library(PCAtools)
library(DESeq2)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(edgeR)
library(ashr)
library(cowplot)
library(wesanderson)
library(UpSetR)
library(EnhancedVolcano)
library(ggvenn)

################################### LOAD DATA ##################################

load("r_outputs/05e-gcb_vars.Rdata")
load("r_outputs/05f-agirre_vars.Rdata")
load("r_outputs/05e_gcb_vars_no_gcb.Rdata")
load("r_outputs/05f-agirre_vars_no_pb.Rdata")
load("r_outputs/01-metadata.Rdata")
load("r_outputs/01-refs.Rdata")
load("r_outputs/05f-agirre_vars_naive.Rdata")
load("r_outputs/05_gcb_naive_vars.Rdata")

remove(BL_metadata, DLBCL_metadata, FL_metadata, all_metadata)

################################# COMMON HERVS #################################

common_up_hervs <- list(
  "DZ" = intersect(sig_gcb$DZ$display[sig_gcb$DZ$class == "LTR" &
                                        sig_gcb$DZ$log2FoldChange > 0], 
                   sig_agirre$DZ$display[sig_agirre$DZ$class == "LTR" &
                                           sig_agirre$DZ$log2FoldChange > 0]),
  "LZ" = intersect(sig_gcb$LZ$display[sig_gcb$LZ$class == "LTR" &
                                        sig_gcb$LZ$log2FoldChange > 0], 
                   sig_agirre$LZ$display[sig_agirre$LZ$class == "LTR" &
                                           sig_agirre$LZ$log2FoldChange > 0]),
  "MB" = intersect(sig_gcb$MB$display[sig_gcb$MB$class == "LTR" &
                                        sig_gcb$MB$log2FoldChange > 0], 
                   sig_agirre$MB$display[sig_agirre$MB$class == "LTR" &
                                           sig_agirre$MB$log2FoldChange > 0]),
  "NB" = intersect(sig_gcb$NB$display[sig_gcb$NB$class == "LTR" &
                                        sig_gcb$NB$log2FoldChange > 0], 
                   sig_agirre$NB$display[sig_agirre$NB$class == "LTR" &
                                           sig_agirre$NB$log2FoldChange > 0]),
  "DZvLZ" = intersect(sig_gcb$DZvLZ$display[sig_gcb$DZvLZ$class == "LTR" &
                                              sig_gcb$DZvLZ$log2FoldChange > 0], 
                      sig_agirre$DZvLZ$display[sig_agirre$DZvLZ$class == "LTR" &
                                                 sig_agirre$DZvLZ$log2FoldChange > 0]),
  "MBvsNB" = intersect(sig_gcb$MBvsNB$display[sig_gcb$MBvsNB$class == "LTR" &
                                                sig_gcb$MBvsNB$log2FoldChange > 0], 
                       sig_agirre$MBvsNB$display[sig_agirre$MBvsNB$class == "LTR" &
                                                   sig_agirre$MBvsNB$log2FoldChange > 0])
)


common_down_hervs <- list(
  "DZ" = intersect(sig_gcb$DZ$display[sig_gcb$DZ$class == "LTR" &
                                        sig_gcb$DZ$log2FoldChange < 0], 
                   sig_agirre$DZ$display[sig_agirre$DZ$class == "LTR" &
                                           sig_agirre$DZ$log2FoldChange < 0]),
  "LZ" = intersect(sig_gcb$LZ$display[sig_gcb$LZ$class == "LTR" &
                                        sig_gcb$LZ$log2FoldChange < 0], 
                   sig_agirre$LZ$display[sig_agirre$LZ$class == "LTR" &
                                           sig_agirre$LZ$log2FoldChange < 0]),
  "MB" = intersect(sig_gcb$MB$display[sig_gcb$MB$class == "LTR" &
                                        sig_gcb$MB$log2FoldChange < 0], 
                   sig_agirre$MB$display[sig_agirre$MB$class == "LTR" &
                                           sig_agirre$MB$log2FoldChange < 0]),
  "NB" = intersect(sig_gcb$NB$display[sig_gcb$NB$class == "LTR" &
                                        sig_gcb$NB$log2FoldChange < 0], 
                   sig_agirre$NB$display[sig_agirre$NB$class == "LTR" &
                                           sig_agirre$NB$log2FoldChange < 0]),
  "DZvLZ" = intersect(sig_gcb$DZvLZ$display[sig_gcb$DZvLZ$class == "LTR" &
                                              sig_gcb$DZvLZ$log2FoldChange < 0], 
                      sig_agirre$DZvLZ$display[sig_agirre$DZvLZ$class == "LTR" &
                                                 sig_agirre$DZvLZ$log2FoldChange < 0]),
  "MBvsNB" = intersect(sig_gcb$MBvsNB$display[sig_gcb$MBvsNB$class == "LTR" &
                                                sig_gcb$MBvsNB$log2FoldChange < 0], 
                       sig_agirre$MBvsNB$display[sig_agirre$MBvsNB$class == "LTR" &
                                                   sig_agirre$MBvsNB$log2FoldChange < 0])
)

sink(file = "r_outputs/05g-common_up_hervs_gcb_agirre.txt")
for (n in names(common_up_hervs)) {
  cat("\n----#Common upregulated HERVs in ", n, "---#\n")
  print(common_up_hervs[[n]])
}
sink(file=NULL)

sink(file = "r_outputs/05g-common_down_hervs_gcb_agirre.txt")
for (n in names(common_down_hervs)) {
  cat("\n----#Common downregulated HERVs in ", n, "---#\n")
  print(common_down_hervs[[n]])
}
sink(file=NULL)

########################## COMMON HERVS NO GCB NO PB ###########################

common_up_hervs_no_gcb_no_pb <- list(
  "DZ" = intersect(sig_no_gcb$DZ$display[sig_no_gcb$DZ$class == "LTR" &
                                           sig_no_gcb$DZ$log2FoldChange > 0], 
                   sig_agirre_no_pb$DZ$display[sig_agirre_no_pb$DZ$class == "LTR" &
                                                 sig_agirre_no_pb$DZ$log2FoldChange > 0]),
  "LZ" = intersect(sig_no_gcb$LZ$display[sig_no_gcb$LZ$class == "LTR" &
                                           sig_no_gcb$LZ$log2FoldChange > 0], 
                   sig_agirre_no_pb$LZ$display[sig_agirre_no_pb$LZ$class == "LTR" &
                                                 sig_agirre_no_pb$LZ$log2FoldChange > 0]),
  "MB" = intersect(sig_no_gcb$MB$display[sig_no_gcb$MB$class == "LTR" &
                                           sig_no_gcb$MB$log2FoldChange > 0], 
                   sig_agirre_no_pb$MB$display[sig_agirre_no_pb$MB$class == "LTR" &
                                           sig_agirre_no_pb$MB$log2FoldChange > 0]),
  "NB" = intersect(sig_no_gcb$NB$display[sig_no_gcb$NB$class == "LTR" &
                                           sig_no_gcb$NB$log2FoldChange > 0], 
                   sig_agirre_no_pb$NB$display[sig_agirre_no_pb$NB$class == "LTR" &
                                                 sig_agirre_no_pb$NB$log2FoldChange > 0]),
  "DZvLZ" = intersect(sig_no_gcb$DZvLZ$display[sig_no_gcb$DZvLZ$class == "LTR" &
                                                 sig_no_gcb$DZvLZ$log2FoldChange > 0], 
                      sig_agirre_no_pb$DZvLZ$display[sig_agirre_no_pb$DZvLZ$class == "LTR" &
                                                       sig_agirre_no_pb$DZvLZ$log2FoldChange > 0]),
  "MBvsNB" = intersect(sig_no_gcb$MBvsNB$display[sig_no_gcb$MBvsNB$class == "LTR" &
                                                   sig_no_gcb$MBvsNB$log2FoldChange > 0], 
                       sig_agirre_no_pb$MBvsNB$display[sig_agirre_no_pb$MBvsNB$class == "LTR" &
                                                   sig_agirre_no_pb$MBvsNB$log2FoldChange > 0])
)

common_down_hervs_no_gcb_no_pb <- list(
  "DZ" = intersect(sig_no_gcb$DZ$display[sig_no_gcb$DZ$class == "LTR" &
                                           sig_no_gcb$DZ$log2FoldChange < 0], 
                   sig_agirre_no_pb$DZ$display[sig_agirre_no_pb$DZ$class == "LTR" &
                                                 sig_agirre_no_pb$DZ$log2FoldChange < 0]),
  "LZ" = intersect(sig_no_gcb$LZ$display[sig_no_gcb$LZ$class == "LTR" &
                                           sig_no_gcb$LZ$log2FoldChange < 0], 
                   sig_agirre_no_pb$LZ$display[sig_agirre_no_pb$LZ$class == "LTR" &
                                                 sig_agirre_no_pb$LZ$log2FoldChange < 0]),
  "MB" = intersect(sig_no_gcb$MB$display[sig_no_gcb$MB$class == "LTR" &
                                           sig_no_gcb$MB$log2FoldChange < 0], 
                   sig_agirre_no_pb$MB$display[sig_agirre_no_pb$MB$class == "LTR" &
                                                 sig_agirre_no_pb$MB$log2FoldChange < 0]),
  "NB" = intersect(sig_no_gcb$NB$display[sig_no_gcb$NB$class == "LTR" &
                                           sig_no_gcb$NB$log2FoldChange < 0], 
                   sig_agirre_no_pb$NB$display[sig_agirre_no_pb$NB$class == "LTR" &
                                                 sig_agirre_no_pb$NB$log2FoldChange < 0]),
  "DZvLZ" = intersect(sig_no_gcb$DZvLZ$display[sig_no_gcb$DZvLZ$class == "LTR" &
                                                 sig_no_gcb$DZvLZ$log2FoldChange < 0], 
                      sig_agirre_no_pb$DZvLZ$display[sig_agirre_no_pb$DZvLZ$class == "LTR" &
                                                       sig_agirre_no_pb$DZvLZ$log2FoldChange < 0]),
  "MBvsNB" = intersect(sig_no_gcb$MBvsNB$display[sig_no_gcb$MBvsNB$class == "LTR" &
                                                   sig_no_gcb$MBvsNB$log2FoldChange < 0], 
                       sig_agirre_no_pb$MBvsNB$display[sig_agirre_no_pb$MBvsNB$class == "LTR" &
                                                         sig_agirre_no_pb$MBvsNB$log2FoldChange < 0])
)

sink(file = "r_outputs/05g-common_up_hervs_gcb_agirre_nogcbpb.txt")
for (n in names(common_up_hervs_no_gcb_no_pb)) {
  cat("\n----#Common upregulated HERVs in ", n, "---#\n")
  print(common_up_hervs_no_gcb_no_pb[[n]])
}
sink(file=NULL)

sink(file = "r_outputs/05g-common_down_hervs_gcb_agirre_nogcbpb.txt")
for (n in names(common_down_hervs_no_gcb_no_pb)) {
  cat("\n----#Common downregulated HERVs in ", n, "---#\n")
  print(common_down_hervs_no_gcb_no_pb[[n]])
}
sink(file=NULL)

######################## COMMON HERVS NAIVE AS BASELINE ########################

common_up_hervs_naive <- list(
  "DZ" = intersect(sig_naive$DZ$display[sig_naive$DZ$class == "LTR" &
                                        sig_naive$DZ$log2FoldChange > 0], 
                   sig_agirre_naive$DZ$display[sig_agirre_naive$DZ$class == "LTR" &
                                           sig_agirre_naive$DZ$log2FoldChange > 0]),
  "LZ" = intersect(sig_naive$LZ$display[sig_naive$LZ$class == "LTR" &
                                        sig_naive$LZ$log2FoldChange > 0], 
                   sig_agirre_naive$LZ$display[sig_agirre_naive$LZ$class == "LTR" &
                                           sig_agirre_naive$LZ$log2FoldChange > 0]),
  "MB" = intersect(sig_naive$MB$display[sig_naive$MB$class == "LTR" &
                                        sig_naive$MB$log2FoldChange > 0], 
                   sig_agirre_naive$MB$display[sig_agirre_naive$MB$class == "LTR" &
                                           sig_agirre_naive$MB$log2FoldChange > 0]),
  "NB" = intersect(sig_naive$NB$display[sig_naive$NB$class == "LTR" &
                                        sig_naive$NB$log2FoldChange > 0], 
                   sig_agirre_naive$NB$display[sig_agirre_naive$NB$class == "LTR" &
                                           sig_agirre_naive$NB$log2FoldChange > 0])
  )


common_down_hervs_naive <- list(
  "DZ" = intersect(sig_naive$DZ$display[sig_naive$DZ$class == "LTR" &
                                        sig_naive$DZ$log2FoldChange < 0], 
                   sig_agirre_naive$DZ$display[sig_agirre_naive$DZ$class == "LTR" &
                                           sig_agirre_naive$DZ$log2FoldChange < 0]),
  "LZ" = intersect(sig_naive$LZ$display[sig_naive$LZ$class == "LTR" &
                                        sig_naive$LZ$log2FoldChange < 0], 
                   sig_agirre_naive$LZ$display[sig_agirre_naive$LZ$class == "LTR" &
                                           sig_agirre_naive$LZ$log2FoldChange < 0]),
  "MB" = intersect(sig_naive$MB$display[sig_naive$MB$class == "LTR" &
                                        sig_naive$MB$log2FoldChange < 0], 
                   sig_agirre_naive$MB$display[sig_agirre_naive$MB$class == "LTR" &
                                           sig_agirre_naive$MB$log2FoldChange < 0]),
  "NB" = intersect(sig_naive$NB$display[sig_naive$NB$class == "LTR" &
                                        sig_naive$NB$log2FoldChange < 0], 
                   sig_agirre_naive$NB$display[sig_agirre_naive$NB$class == "LTR" &
                                           sig_agirre_naive$NB$log2FoldChange < 0])
  )


##################################### VENN ##################################### 

################################# DARK ZONE UP ################################# 

makevenn <- function(clust) {
  
  a <- ggvenn(list(
    `GCB` = sig_gcb[[clust]]$display[sig_gcb[[clust]]$class == "LTR" &
                                      sig_gcb[[clust]]$log2FoldChange > 0],
    `Agirre` = sig_agirre[[clust]]$display[sig_agirre[[clust]]$class == "LTR" &
                                            sig_agirre[[clust]]$log2FoldChange > 0]), 
    set_name_size = 3.5, 
    text_size = 3.5,
    c("GCB", "Agirre"),
    fill_color = c("#f8766d", "#00bfc4"))
  
  a_down <- ggvenn(list(
    `GCB` = sig_gcb[[clust]]$display[sig_gcb[[clust]]$class == "LTR" &
                                       sig_gcb[[clust]]$log2FoldChange < 0 ],
    `Agirre` = sig_agirre[[clust]]$display[sig_agirre[[clust]]$class == "LTR" &
                                             sig_agirre[[clust]]$log2FoldChange < 0 ]), 
    set_name_size = 3.5, 
    text_size = 3.5,
    c("GCB", "Agirre"),
    fill_color = c("#f8766d", "#00bfc4"))
  
  a_diff <- ggvenn(list(
    `GCB` = sig_gcb[[clust]]$display[sig_gcb[[clust]]$class == "LTR"],
    `Agirre` = sig_agirre[[clust]]$display[sig_agirre[[clust]]$class == "LTR"]), 
    set_name_size = 3.5, 
    text_size = 3.5,
    c("GCB", "Agirre"),
    fill_color = c("#f8766d", "#00bfc4"))
  
  
  a_gene <- ggvenn(list(
    `GCB` = sig_gcb[[clust]]$display[sig_gcb[[clust]]$class == "protein_coding" &
                                    sig_gcb[[clust]]$log2FoldChange > 0],
    `Agirre` = sig_agirre[[clust]]$display[sig_agirre[[clust]]$class == "protein_coding" &
                                          sig_agirre[[clust]]$log2FoldChange > 0]), 
    c("GCB", "Agirre"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  a_gene_down <- ggvenn(list(
    `GCB` = sig_gcb[[clust]]$display[sig_gcb[[clust]]$class == "protein_coding" &
                                       sig_gcb[[clust]]$log2FoldChange < 0 ],
    `Agirre` = sig_agirre[[clust]]$display[sig_agirre[[clust]]$class == "protein_coding" &
                                             sig_agirre[[clust]]$log2FoldChange < 0 ]), 
    c("GCB", "Agirre"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  a_gene_diff <- ggvenn(list(
    `GCB` = sig_gcb[[clust]]$display[sig_gcb[[clust]]$class == "protein_coding"],
    `Agirre` = sig_agirre[[clust]]$display[sig_agirre[[clust]]$class == "protein_coding"]), 
    c("GCB", "Agirre"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  b <- ggvenn(list(
    `GCB (no GCB)` = sig_no_gcb[[clust]]$display[sig_no_gcb[[clust]]$class == "LTR" &
                                       sig_no_gcb[[clust]]$log2FoldChange > 0],
    `Agirre (no PB)` = sig_agirre_no_pb[[clust]]$display[sig_agirre_no_pb[[clust]]$class == "LTR" &
                                                sig_agirre_no_pb[[clust]]$log2FoldChange > 0]), 
    c("GCB (no GCB)", "Agirre (no PB)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  b_down <- ggvenn(list(
    `GCB (no GCB)` = sig_no_gcb[[clust]]$display[sig_no_gcb[[clust]]$class == "LTR" &
                                                   sig_no_gcb[[clust]]$log2FoldChange < 0 ],
    `Agirre (no PB)` = sig_agirre_no_pb[[clust]]$display[sig_agirre_no_pb[[clust]]$class == "LTR" &
                                                           sig_agirre_no_pb[[clust]]$log2FoldChange < 0 ]), 
    c("GCB (no GCB)", "Agirre (no PB)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  b_diff <- ggvenn(list(
    `GCB (no GCB)` = sig_no_gcb[[clust]]$display[sig_no_gcb[[clust]]$class == "LTR"],
    `Agirre (no PB)` = sig_agirre_no_pb[[clust]]$display[sig_agirre_no_pb[[clust]]$class == "LTR"]), 
    c("GCB (no GCB)", "Agirre (no PB)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  b_gene <- ggvenn(list(
    `GCB (no GCB)` = sig_no_gcb[[clust]]$display[sig_no_gcb[[clust]]$class == "protein_coding" &
                                                sig_no_gcb[[clust]]$log2FoldChange > 0],
    `Agirre (no PB)` = sig_agirre_no_pb[[clust]]$display[sig_agirre_no_pb[[clust]]$class == "protein_coding" &
                                                        sig_agirre_no_pb[[clust]]$log2FoldChange > 0]), 
    c("GCB (no GCB)", "Agirre (no PB)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  b_gene_down <- ggvenn(list(
    `GCB (no GCB)` = sig_no_gcb[[clust]]$display[sig_no_gcb[[clust]]$class == "protein_coding" &
                                                   sig_no_gcb[[clust]]$log2FoldChange < 0],
    `Agirre (no PB)` = sig_agirre_no_pb[[clust]]$display[sig_agirre_no_pb[[clust]]$class == "protein_coding" &
                                                           sig_agirre_no_pb[[clust]]$log2FoldChange < 0]), 
    c("GCB (no GCB)", "Agirre (no PB)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  b_gene_diff <- ggvenn(list(
    `GCB (no GCB)` = sig_no_gcb[[clust]]$display[sig_no_gcb[[clust]]$class == "protein_coding"],
    `Agirre (no PB)` = sig_agirre_no_pb[[clust]]$display[sig_agirre_no_pb[[clust]]$class == "protein_coding"]), 
    c("GCB (no GCB)", "Agirre (no PB)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  c <- ggvenn(list(
    `GCB (v naive)` = sig_naive[[clust]]$display[sig_naive[[clust]]$class == "LTR" &
                                                sig_naive[[clust]]$log2FoldChange > 0],
    `Agirre (v naive)` = sig_agirre_naive[[clust]]$display[sig_agirre_naive[[clust]]$class == "LTR" &
                                                          sig_agirre_naive[[clust]]$log2FoldChange > 0]), 
    c("GCB (v naive)", "Agirre (v naive)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  c_down <- ggvenn(list(
    `GCB (v naive)` = sig_naive[[clust]]$display[sig_naive[[clust]]$class == "LTR" &
                                                   sig_naive[[clust]]$log2FoldChange < 0],
    `Agirre (v naive)` = sig_agirre_naive[[clust]]$display[sig_agirre_naive[[clust]]$class == "LTR" &
                                                             sig_agirre_naive[[clust]]$log2FoldChange < 0]), 
    c("GCB (v naive)", "Agirre (v naive)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  c_diff <- ggvenn(list(
    `GCB (v naive)` = sig_naive[[clust]]$display[sig_naive[[clust]]$class == "LTR"],
    `Agirre (v naive)` = sig_agirre_naive[[clust]]$display[sig_agirre_naive[[clust]]$class == "LTR"]), 
    c("GCB (v naive)", "Agirre (v naive)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  c_gene <- ggvenn(list(
    `GCB (v naive)` = sig_naive[[clust]]$display[sig_naive[[clust]]$class == "protein_coding" &
                                                sig_naive[[clust]]$log2FoldChange > 0],
    `Agirre (v naive)` = sig_agirre_naive[[clust]]$display[sig_agirre_naive[[clust]]$class == "protein_coding" &
                                                          sig_agirre_naive[[clust]]$log2FoldChange > 0]), 
    c("GCB (v naive)", "Agirre (v naive)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  c_gene_down <- ggvenn(list(
    `GCB (v naive)` = sig_naive[[clust]]$display[sig_naive[[clust]]$class == "protein_coding" &
                                                   sig_naive[[clust]]$log2FoldChange < 0],
    `Agirre (v naive)` = sig_agirre_naive[[clust]]$display[sig_agirre_naive[[clust]]$class == "protein_coding" &
                                                             sig_agirre_naive[[clust]]$log2FoldChange < 0]), 
    c("GCB (v naive)", "Agirre (v naive)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  c_gene_diff <- ggvenn(list(
    `GCB (v naive)` = sig_naive[[clust]]$display[sig_naive[[clust]]$class == "protein_coding"],
    `Agirre (v naive)` = sig_agirre_naive[[clust]]$display[sig_agirre_naive[[clust]]$class == "protein_coding"]), 
    c("GCB (v naive)", "Agirre (v naive)"),
    set_name_size = 3.5, 
    text_size = 3.5,
    fill_color = c("#f8766d", "#00bfc4"))
  
  clust.fig.up <- plot_grid(a, a_gene, b, b_gene, c, c_gene, 
            labels = "AUTO",
            ncol = 2)
  
  clust.fig.down <- plot_grid(a_down, a_gene_down, b_down, b_gene_down, c_down, c_gene_down, 
                              labels = "AUTO",
                              ncol = 2)
  
  clust.fig.diff <- plot_grid(a_diff, a_gene_diff, b_diff, b_gene_diff, c_diff, c_gene_diff, 
                              labels = "AUTO",
                              ncol = 2)
  
  pdf(paste0("plots/05g_gcb_v_agirre_", clust, "_up",".pdf"), height=8, width=10)
  print(clust.fig.up)
  dev.off()
  
  pdf(paste0("plots/05g_gcb_v_agirre_", clust, "_down",".pdf"), height=8, width=10)
  print(clust.fig.down)
  dev.off()
  
  pdf(paste0("plots/05g_gcb_v_agirre_", clust, "_diff",".pdf"), height=8, width=10)
  print(clust.fig.diff)
  dev.off()

}

makevenn("DZ")
makevenn("LZ")
makevenn("MB")
makevenn("NB")