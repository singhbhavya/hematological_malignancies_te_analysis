################################################################################
################################################################################
################################################################################
################################################################################
########################### DLBCL CLASSIFIER ANALYSIS ##########################

## The purpose of this analysis is to see how the same NCICCR and TCGA samples 
## are classified by different classification schemes (COO, scCOO, etc), and to
## visualize data movement from COO to LymphGen to EcoTyper, etc.

#################################### SETUP #####################################

library(knitr)
library(tidyverse)
library(dplyr)
library(matrixStats)
library(data.table)
library(PCAtools)
library(DESeq2)
library(ggplot2)
library(ggsci)
library(edgeR)
library(ashr)
library(cowplot)
library(ggalluvial)
library(wesanderson)

################################## LOAD DATA ###################################

load("r_outputs/02-DLBCL_filt_counts.Rdata")
load("r_outputs/01-metadata.Rdata")

################################ ALLUVIAL PLOT #################################

DLBCL_metadata$age_at_diagnosis <- round(DLBCL_metadata$age_at_diagnosis / 365)

DLBCL_metadata %>% dplyr::count(project, gender)

#project gender   n
#1 NCICCR-DLBCL female 195
#2 NCICCR-DLBCL   male 286
#3    TCGA-DLBC female  26
#4    TCGA-DLBC   male  22

alluvial <- DLBCL_metadata %>% 
  dplyr::count(project, COO_class, DblHit_call, Chapuay_call, EcoTyper_call, 
               LymphGen_call, scCOO_group_call, EcoTyper_call)

alluvial <-
  alluvial %>% 
  mutate(LymphGen_call = strsplit(as.character(LymphGen_call), "/")) %>% 
  unnest(LymphGen_call) 

alluvial$scCOO_group_call <- 
  gsub('Group', 'G', alluvial$scCOO_group_call)

colors <- wes_palette("Zissou1")
ggplot(as.data.frame(alluvial),
       aes(y = n, 
           axis1 = COO_class, 
           axis2 = DblHit_call, 
           axis3 = scCOO_group_call, 
           axis4 = Chapuay_call, 
           axis5 = EcoTyper_call, 
           axis6 =LymphGen_call))  +
  geom_alluvium(aes(fill = COO_class),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/4, 
               reverse = FALSE, 
               fill = "lightgrey", 
               alpha = .7) +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)),
            reverse = FALSE) + coord_flip() + 
  scale_fill_manual(values = c("ABC" = colors[5], 
                               "GCB" = colors[3], 
                               "Unclass" = colors[1])) +
  scale_x_continuous(breaks = 1:6, labels = 
                       c("COO Classification\n (Alizadeh et al, 2000)", 
                         "DBL Hit\n(Ennishi et al, 2019)", 
                         "scCOO Group\n(Holmes et al, 2020)", 
                         "(Chapuy et al, 2018)", 
                         "EcoTyper\n(Luca et al, 2021)", 
                         "LymphGen\n(Steen et al, 2021)")) +
  theme_half_open() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=16))

ggsave("plots/05a-dlbcl_alluvial.pdf", height=5, width=10)
