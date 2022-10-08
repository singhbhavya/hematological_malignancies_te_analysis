################################################################################
################################################################################
################################################################################
################################################################################
#############################  LOAD SINGLE CELL DATA ###########################

## Plan
## - Load data and create Seurat object
## - Seurat / Stellarscope QC
## - Seurat normalization
## - Merge data
## - Integration
## - Cell type annotation
## - Azimuth 

#################################### SETUP #####################################

library(tidyr)
library(scopetools)
library(scater)
library(Seurat)
library(edgeR)
library(Seurat)
library(SeuratObject)
library(rtracklayer)

################################### METADATA ###################################

# Import data
GCB_metadata <- read.csv("metadata/GCB/GCB_3_SRA.txt", sep = ",")
GCB_metadata$sample <- GCB_metadata$BioSample
GCB_metadata <- GCB_metadata %>% remove_rownames %>% column_to_rownames(var="BioSample")

# Remove the dark zone and light zone sorted cells (for now)
GCB_metadata <- GCB_metadata[!(row.names(GCB_metadata) %in% c("SAMN13191511",
                                                              "SAMN13191510",
                                                              "SAMN13191508",
                                                              "SAMN13191507")),]

################################# LOAD SEURAT ##################################

# Function to load Seurat objects for all samples
load_all_seurat <- function(i) {
  sample_name <- GCB_metadata$sample[GCB_metadata$sample == i]
  
  seurat_object <- 
    scopetools::load_stellarscope_seurat(stellarscope_dir = 
                                           paste("results/GCB/telescope_pseudobulk/", i, "_F/", sep = ""),
                                         starsolo_dir = 
                                           paste("results/GCB/starsolo_alignment/", i, "/", i, "_GDC38.Solo.out/Gene/filtered/", sep = ""),
                                         exp_tag = paste(i, "_F_pseudobulk", sep = ""),
                                         use.symbols = TRUE,
                                         TE_metadata = retro.hg38.v1)
  
  assign(paste(sample_name, "seurat", sep="."), seurat_object, envir=.GlobalEnv)
}

# Apply function to all samples 
invisible(lapply(GCB_metadata$sample, load_all_seurat))

# Add seurat obj names to sample sheet
GCB_metadata$seurat_obj_og <- paste(GCB_metadata$sample, ".seurat", sep="")

############################### STELLARSCOPE QC ################################

SAMN14979762.seurat.qc <- stellarscope_cell_qc(SAMN14979762.seurat)
SAMN14979745.seurat.qc <- stellarscope_cell_qc(SAMN14979745.seurat)
SAMN13191513.seurat.qc <- stellarscope_cell_qc(SAMN13191513.seurat)
SAMN13191512.seurat.qc <- stellarscope_cell_qc(SAMN13191512.seurat)

GCB_metadata$seurat_obj_qc <- paste(GCB_metadata$sample, ".seurat.qc", sep="")

  ############################### NORMALIZE DATA ###############################

SAMN14979762.seurat.norm <- Seurat::NormalizeData(SAMN14979762.seurat.qc)
SAMN14979745.seurat.norm <- Seurat::NormalizeData(SAMN14979745.seurat.qc)
SAMN13191513.seurat.norm <- Seurat::NormalizeData(SAMN13191513.seurat.qc)
SAMN13191512.seurat.norm <- Seurat::NormalizeData(SAMN13191512.seurat.qc)

################################## MERGE DATA ##################################

# Merge data that was first normalized
GCB.norm.merged <- merge(SAMN14979762.seurat.norm, 
                         y = c(SAMN14979745.seurat.norm, SAMN13191513.seurat.norm, 
                               SAMN13191512.seurat.norm), 
                         add.cell.ids = c("SAMN14979762", "SAMN14979745",
                                          "SAMN13191513", "SAMN13191512"), 
                         project = "GCB",
                         merge.data = TRUE)

# Merge un-normalized data, and then normalize 
GCB.merged.norm <- merge(SAMN14979762.seurat.qc, 
                         y = c(SAMN14979745.seurat.qc, SAMN13191513.seurat.qc, 
                               SAMN13191512.seurat.qc), 
                         add.cell.ids = c("SAMN14979762", "SAMN14979745",
                                          "SAMN13191513", "SAMN13191512"), 
                         project = "GCB")

GCB.merged.norm <- Seurat::NormalizeData(GCB.merged.norm)


################################## SAVE DATA ###################################

save(SAMN14979762.seurat.norm, SAMN14979745.seurat.norm, 
     SAMN13191513.seurat.norm, SAMN13191512.seurat.norm, GCB_metadata, 
     file="r_outputs/03-gcb_seurat_objects_individual.Rdata")

save(GCB.norm.merged, GCB.merged.norm, 
     file="r_outputs/03-gcb_seurat_objects_merged.Rdata")