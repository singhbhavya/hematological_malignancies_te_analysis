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
load("r_outputs/05o-DLBCL_pca_ccp_clusters_metadata.Rdata")


################################ ALLUVIAL PLOT #################################


# Add clusters to metadata
DLBCL_metadata$clust.retro.k2 <- clust.df$clust.retro.k2
DLBCL_metadata$clust.retro.k3 <- clust.df$clust.retro.k3
DLBCL_metadata$clust.retro.k4 <- clust.df$clust.retro.k4
DLBCL_metadata$clust.retro.k5 <- clust.df$clust.retro.k5
DLBCL_metadata$clust.retro.k7 <- clust.df$clust.retro.k7
DLBCL_metadata$clust.retro.k9 <- clust.df$clust.retro.k9

DLBCL_metadata$age_at_diagnosis <- round(DLBCL_metadata$age_at_diagnosis / 365)

DLBCL_metadata %>% dplyr::count(project, gender)

#project gender   n
#1 NCICCR-DLBCL female 195
#2 NCICCR-DLBCL   male 286
#3    TCGA-DLBC female  26
#4    TCGA-DLBC   male  22

alluvial <- DLBCL_metadata %>% 
  dplyr::count(project, COO_class, DblHit_call, Chapuay_call, EcoTyper_call, 
               LymphGen_call, scCOO_group_call, EcoTyper_call, clust.retro.k7)

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
           axis6 =LymphGen_call,
           axis7 = clust.retro.k7))  +
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
  scale_fill_manual(values = c("ABC" = "red3", 
                               "GCB" = "royalblue", 
                               "Unclass" = "lightblue")) +
  scale_x_continuous(breaks = 1:7, labels = 
                       c("COO Classification\n (Alizadeh et al, 2000)", 
                         "DBL Hit\n(Ennishi et al, 2019)", 
                         "scCOO Group\n(Holmes et al, 2020)", 
                         "(Chapuy et al, 2018)", 
                         "EcoTyper\n(Luca et al, 2021)", 
                         "LymphGen\n(Steen et al, 2021)",
                         "clust.retro.k7")) +
  theme_half_open() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=16))

ggsave("plots/05a-dlbcl_alluvial_k7.pdf", height=5, width=10)

################################ ALLUVIAL PLOT 2 #################################

alluvial <- DLBCL_metadata %>% 
  dplyr::count(COO_class, EcoTyper_call, LymphGen_call, EcoTyper_call, clust.retro.k7)

alluvial <-
  alluvial %>% 
  mutate(LymphGen_call = strsplit(as.character(LymphGen_call), "/")) %>% 
  unnest(LymphGen_call) 

ggplot(as.data.frame(alluvial),
       aes(y = n, 
           axis1 = COO_class, 
           axis2 = clust.retro.k7))  +
  guides(fill = FALSE) +
  geom_alluvium(aes(fill = COO_class),
                width = 0, knot.pos = 0, reverse = FALSE) +
  geom_stratum(width = 1/4, 
               reverse = FALSE, 
               fill = "lightgrey", 
               alpha = .7) +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)),
            reverse = FALSE) + coord_flip() + 
  scale_fill_manual(values = c("ABC" = "red3", 
                               "GCB" = "royalblue", 
                               "Unclass" = "lightblue")) +
  scale_x_continuous(breaks = 1:2, labels = 
                       c("COO Classification\n (Alizadeh et al, 2000)", 
                         "clust.retro.k7")) +
  theme_half_open() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=16))

ggsave("plots/05a-dlbcl_alluvial_k7_coo_only.pdf", height=3, width=10)


ggplot(as.data.frame(alluvial),
       aes(y = n, 
           axis1 = LymphGen_call, 
           axis2 = clust.retro.k7))  +
  guides(fill = FALSE) +
  geom_alluvium(aes(fill = LymphGen_call),
                width = 0, knot.pos = 0, reverse = FALSE) +
  geom_stratum(width = 1/4, 
               reverse = FALSE, 
               fill = "lightgrey", 
               alpha = .7) +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)),
            reverse = FALSE) + coord_flip() + 
  scale_fill_manual(values = c("A53" = "#e78bf0",
                               "BN2" = "#e28743",
                               "EZB" = "#3d3aaf",
                               "MCD" = "#4d8347",
                               "N1" = "#e84646",
                               "Other" = "#e9f190",
                               "ST2" = "#38baaa")) +
  scale_x_continuous(breaks = 1:2, labels = 
                       c("LymphGen\n(Steen et al, 2021)", 
                         "clust.retro.k7")) +
  theme_half_open() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=16))

ggsave("plots/05a-dlbcl_alluvial_k7_lymphgen_only.pdf", height=3, width=10)

ggplot(as.data.frame(alluvial),
       aes(y = n, 
           axis1 = EcoTyper_call, 
           axis2 = clust.retro.k7))  +
  guides(fill = FALSE) +
  geom_alluvium(aes(fill = EcoTyper_call),
                width = 0, knot.pos = 0, reverse = FALSE) +
  geom_stratum(width = 1/4, 
               reverse = FALSE, 
               fill = "lightgrey", 
               alpha = .7) +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)),
            reverse = FALSE) + coord_flip() + 
  scale_fill_manual(values = c("S1" = "#bcc779",
                               "S2" = "#60aaec",
                               "S3" = "#Ecaf60",
                               "S4" = "#B789EE",
                               "S5" = "#E0AF71",
                               "Unassigned" = "grey")) +
  scale_x_continuous(breaks = 1:2, labels = 
                       c("EcoTyper)", 
                         "clust.retro.k7")) +
  theme_half_open() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=16))

ggsave("plots/05a-dlbcl_alluvial_k7_ecotyper_only.pdf", height=3, width=10)

################################## BAR PLOTS ###################################


ggplot(DLBCL_metadata %>% 
         dplyr::count(clust.retro.k7,
                      COO_class, sort = TRUE), 
       aes(fill=COO_class, y=n, x=clust.retro.k7)) + 
  geom_bar(position="stack", stat="identity", alpha=0.7) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.line=element_blank()) + 
  scale_fill_manual(values = c("ABC" = "red3", 
                               "GCB" = "royalblue", 
                               "Unclass" = "lightblue")) +
  xlab("DLBCL cluster") +
  ylab("Number of samples") + 
  guides(fill=guide_legend(title="COO class"))

ggsave("plots/05a-dlbcl_barplot_k7_coo_only.pdf", height=5, width=5)

ggplot(DLBCL_metadata %>% 
         dplyr::count(clust.retro.k7,
                      LymphGen_call, sort = TRUE), 
       aes(fill=LymphGen_call, y=n, x=clust.retro.k7)) + 
  geom_bar(position="stack", stat="identity", alpha=0.7) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.line=element_blank()) + 
  scale_fill_manual(values = c("A53" = "#e78bf0",
                               "BN2" = "#e28743",
                               "BN2/A53" = "#AF6C3A",
                               "BN2/EZB" = "#c79875",
                               "BN2/MCD" = "#7a4c29",
                               "BN2/ST2" = "#58361d",
                               "EZB" = "#3d3aaf",
                               "EZB/A53" = "#7775c7",
                               "EZB/MCD" = "#8b89cf",
                               "EZB/N1/ST2/A53" = "#b1b0df",
                               "EZB/ST2/A53" = "#d8d8ef",
                               "MCD" = "#4d8347",
                               "MCD/A53" = "#82a87e",
                               "MCD/ST2" = "#b8cdb5",
                               "N1" = "#e84646",
                               "N1/ST2" = "#f19090",
                               "Other" = "#e9f190",
                               "ST2" = "#38baaa",
                               "NA" = "grey")) +
  xlab("DLBCL cluster") +
  ylab("Number of samples") + 
  guides(fill=guide_legend(title="LymphGen class")) 

ggsave("plots/05a-dlbcl_barplot_k7_lymphgen_only.pdf", height=5, width=5)

ggplot(DLBCL_metadata %>% 
         dplyr::count(clust.retro.k7,
                      EcoTyper_call, sort = TRUE), 
       aes(fill=EcoTyper_call, y=n, x=clust.retro.k7)) + 
  geom_bar(position="stack", stat="identity", alpha=0.7) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.line=element_blank()) + 
  scale_fill_manual(values = c("S1" = "#bcc779",
                               "S2" = "#60aaec",
                               "S3" = "#Ecaf60",
                               "S4" = "#B789EE",
                               "S5" = "lightpink",
                               "Unassigned" = "grey",
                               "Missing" = "grey")) +
  xlab("DLBCL cluster") +
  ylab("Number of samples") + 
  guides(fill=guide_legend(title="EcoTyper class"))

ggsave("plots/05a-dlbcl_barplot_k7_ecotyper_only.pdf", height=5.5, width=5.5)
