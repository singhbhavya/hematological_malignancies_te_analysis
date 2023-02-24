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

################################### LOAD DATA ##################################

load("r_outputs/05e-gcb_vars.Rdata")
load("r_outputs/05f-agirre_vars.Rdata")
load("r_outputs/01-metadata.Rdata")
load("r_outputs/01-refs.Rdata")

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


