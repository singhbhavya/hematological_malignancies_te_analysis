################################################################################
################################################################################
################################################################################
################################################################################
################################## GCB FGSEA ###################################

## Plan:
## DESEq2 on cell types together
## FGEAS on cell types comparing between datasets

#################################### SETUP #####################################

library(knitr)
library(tidyverse)
library(matrixStats)
library(data.table)
library(dplyr)
library(PCAtools)
library(DESeq2)
library(ggplot2)
library(ggsci)
library(edgeR)
library(ashr)
library(cowplot)
library(wesanderson)
library(fgsea)

################################### LOAD DATA ##################################

load("r_outputs/01-GCB_Bulk_counts.Rdata")
load("r_outputs/01-GCB_Agirre.Rdata")
load("r_outputs/01-refs.Rdata")

################################### PATHWAYS ###################################

pathways.hallmark <- gmtPathways("gsea/h.all.v2023.1.Hs.symbols.gmt")
pathways.immune <- gmtPathways("gsea/c7.immunesigdb.v2023.1.Hs.symbols.gmt")

################################# METADATA SETUP ###############################

agirre_metadata$Cell_type <- agirre_metadata$source_name
agirre_metadata$Cell_type[agirre_metadata$Cell_type == "Naive"] <- "NB"
agirre_metadata$Cell_type[agirre_metadata$Cell_type == "Centroblast"] <- "DZ"
agirre_metadata$Cell_type[agirre_metadata$Cell_type == "Centrocyte"] <- "LZ"
agirre_metadata$Cell_type[agirre_metadata$Cell_type == "Memory"] <- "MB"
agirre_metadata$Cell_type[agirre_metadata$Cell_type == "Tonsilar plasma cell"] <- "PB"
agirre_metadata$Cell_type[agirre_metadata$Cell_type == "Bone Marrow plasma cell"] <- "BMPC"

bulk_metadata$Cell_type[bulk_metadata$Cell_type == "Naive B"] <- "NB"
bulk_metadata$Cell_type[bulk_metadata$Cell_type == "Germinal Center B"] <- "GCB"
bulk_metadata$Cell_type[bulk_metadata$Cell_type == "Memory B"] <- "MB"
bulk_metadata$Cell_type[bulk_metadata$Cell_type == "Dark Zone Germinal Center B"] <- "DZ"
bulk_metadata$Cell_type[bulk_metadata$Cell_type == "Light Zone Germinal Center B"] <- "LZ"


agirre_metadata <- subset(agirre_metadata, select = c("cancer_type", "Cell_type"))
colnames(agirre_metadata) <- c("Dataset", "Cell_type")

bulk_metadata <- subset(bulk_metadata, select = c("cancer_type", "Cell_type"))
colnames(bulk_metadata) <- c("Dataset", "Cell_type")

############################# FUNCTION SCREE PLOT ##############################

pca_standard <- function(tform, metadata, var) {
  
  removeVar <- var
  pca.obj <- PCAtools::pca(assay(tform), 
                           metadata=metadata, 
                           removeVar=removeVar)
  
  cat(sprintf('Removed %d pct low variance variables, %d retained\n', 
              removeVar*100, length(pca.obj$xvars)))
  
  varline <- 50
  varline.x <- min(which(cumsum(pca.obj$variance) >= varline))
  
  horn <- PCAtools::parallelPCA(assay(tform), removeVar = removeVar)
  elbow <- PCAtools::findElbowPoint(pca.obj$variance)
  
  screeplot <-PCAtools::screeplot(pca.obj,
                                  axisLabSize = 6,
                                  components = getComponents(pca.obj, 1:30),
                                  title=paste("Retrotranscriptome SCREE",
                                              metadata$cancer_type[1],
                                              sep=" "),
                                  hline=varline, vline=c(varline.x, horn$n, elbow)
  ) +
    geom_label(aes(x=varline.x+1, y=50, 
                   label = paste0(varline, '% var'), vjust = -1)) +
    geom_label(aes(x = horn$n + 1, y = 50,
                   label = 'Horn\'s', vjust = -1)) +
    geom_label(aes(x = elbow + 1, y = 50,
                   label = 'Elbow method', vjust = -1))
  
  
  cat(sprintf('%d PCs for Elbow method\n', elbow))
  cat(sprintf('%d PCs for Horn method\n', horn$n))
  cat(sprintf('%d PCs needed to explain %d percent of variation\n', 
              varline.x, varline))
  
  return(pca.obj)
}

################################ SET THRESHOLDS ################################

lfc.cutoff <- 1.5
pval=0.001 # p value threshold

################################### DARK ZONE ##################################

dz.metadata  <- rbind(agirre_metadata[agirre_metadata$Cell_type == "DZ",], 
                      bulk_metadata[bulk_metadata$Cell_type == "DZ",])

dz.counts <- cbind(GCB_Agirre.counts.comb[rownames(agirre_metadata[agirre_metadata$Cell_type == "DZ",])],
                   GCB_Bulk.counts.comb[rownames(bulk_metadata[bulk_metadata$Cell_type == "DZ",])])

stopifnot(all(colnames(dz.counts) == rownames(dz.metadata)))

# DESEQ2 

dz.dds <- DESeq2::DESeqDataSetFromMatrix(countData = dz.counts,
                                          colData = dz.metadata,
                                          design = ~ Dataset + 0)

dz.dds <- DESeq2::DESeq(dz.dds, parallel=T)
dz.tform <- DESeq2::varianceStabilizingTransformation(dz.dds, blind=FALSE)

# DESEQ2 results

dz.res <- DESeq2::results(dz.dds, contrast=c("Dataset", "GCB_Agirre", "GCB_Bulk"), 
                          alpha=pval)

dz.res$SYMBOL <- gene_table[rownames(dz.res),]$display

dz.res <- as.data.frame(dz.res) %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))

dz.ranks <- deframe(dz.res)

dz.fgsea <- fgsea(pathways=pathways.immune, stats=dz.ranks, eps=0)

dz.fgseaResTidy <- dz.fgsea %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
dz.fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

dz.pathways.subset <- dz.fgseaResTidy[dz.fgseaResTidy$pathway %like% "DARK", ]

ggplot(dz.pathways.subset, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

################################### LIGHT ZONE #################################

lz.metadata  <- rbind(agirre_metadata[agirre_metadata$Cell_type == "LZ",], 
                      bulk_metadata[bulk_metadata$Cell_type == "LZ",])

lz.counts <- cbind(GCB_Agirre.counts.comb[rownames(agirre_metadata[agirre_metadata$Cell_type == "LZ",])],
                   GCB_Bulk.counts.comb[rownames(bulk_metadata[bulk_metadata$Cell_type == "LZ",])])

stopifnot(all(colnames(lz.counts) == rownames(lz.metadata)))

# DESEQ2 

lz.dds <- DESeq2::DESeqDataSetFromMatrix(countData = lz.counts,
                                         colData = lz.metadata,
                                         design = ~ Dataset + 0)

lz.dds <- DESeq2::DESeq(lz.dds, parallel=T)
lz.tform <- DESeq2::varianceStabilizingTransformation(lz.dds, blind=FALSE)

# DESEQ2 results

lz.res <- DESeq2::results(lz.dds, contrast=c("Dataset", "GCB_Agirre", "GCB_Bulk"), 
                          alpha=pval)

lz.res$SYMBOL <- gene_table[rownames(lz.res),]$display

lz.res <- as.data.frame(lz.res) %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))

lz.ranks <- deframe(lz.res)

lz.fgsea <- fgsea(pathways=pathways.immune, stats=dz.ranks, eps=0)

lz.fgseaResTidy <- lz.fgsea %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
lz.fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

lz.pathways.subset <- lz.fgseaResTidy[lz.fgseaResTidy$pathway %like% "LIGHT", ]

ggplot(lz.pathways.subset, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

#################################### NAIVE B ###################################


nb.metadata  <- rbind(agirre_metadata[agirre_metadata$Cell_type == "NB",], 
                      bulk_metadata[bulk_metadata$Cell_type == "NB",])

nb.counts <- cbind(GCB_Agirre.counts.comb[rownames(agirre_metadata[agirre_metadata$Cell_type == "NB",])],
                   GCB_Bulk.counts.comb[rownames(bulk_metadata[bulk_metadata$Cell_type == "NB",])])

stopifnot(all(colnames(nb.counts) == rownames(nb.metadata)))

# DESEQ2 

nb.dds <- DESeq2::DESeqDataSetFromMatrix(countData = nb.counts,
                                         colData = nb.metadata,
                                         design = ~ Dataset + 0)

nb.dds <- DESeq2::DESeq(nb.dds, parallel=T)
nb.tform <- DESeq2::varianceStabilizingTransformation(nb.dds, blind=FALSE)

# DESEQ2 results

nb.res <- DESeq2::results(nb.dds, contrast=c("Dataset", "GCB_Agirre", "GCB_Bulk"), 
                          alpha=pval)

nb.res$SYMBOL <- gene_table[rownames(nb.res),]$display

nb.res <- as.data.frame(nb.res) %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))

nb.ranks <- deframe(nb.res)

nb.fgsea <- fgsea(pathways=pathways.immune, stats=dz.ranks, eps=0)

nb.fgseaResTidy <- nb.fgsea %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
nb.fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

nb.pathways.subset <- nb.fgseaResTidy[nb.fgseaResTidy$pathway %like% "NAIVE", ]

ggplot(nb.pathways.subset, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
