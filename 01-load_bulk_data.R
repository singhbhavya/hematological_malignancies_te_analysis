################################################################################
################################################################################
################################################################################
################################################################################
##########################  BULK DATA / METADATA SETUP #########################

# Plan:
# - Download metadata from TCGA biolinks (for DLBCL, Burkitt, and Follicular)
# - Build metadata and sample files for all. 
# - Set up count matrices (all, rtx, and hervs only)
# - Set up annotations

#################################### SETUP #####################################

library(tidyverse)
library(readxl)
library(GenomicDataCommons)
library(TCGAbiolinks)
library(dplyr)
library(rtracklayer)
library(data.table)
library(scopetools)

############################# LOAD TE ANNOTATIONS ##############################

## load annotation
retro.hg38.v1 <- 
  readr::read_tsv(
    "https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/retro.hg38.v1/genes.tsv.gz", 
    na=c('.'))
retro.hg38.v1 <- retro.hg38.v1 %>%
  tidyr::separate(locus, c("family"), sep='_', remove=F, extra='drop') %>%
  dplyr::mutate(
    te_class = factor(ifelse(is.na(l1base_id), 'LTR', 'LINE'), levels=c('LTR','LINE')),
  )

# Annotation directory for scopetools
ddir <- system.file("extdata", package="scopetools")

# Remove the confounding LINE element (L1FLnI_Xq21.1db) that has a poly A tail
# in the middle of it:

retro.hg38.v1<-
  retro.hg38.v1[!(retro.hg38.v1$locus=="L1FLnI_Xq21.1db"),]

retro.annot <- retro.hg38.v1
row.names(retro.annot) <- retro.annot$locus
row.names(retro.annot) <- gsub("_", "-", row.names(retro.annot))

############################ LOAD GENE ANNOTATIONS #############################

gtf <- rtracklayer::import("refs/gencode.v38.annotation.gtf")
gtf_df=as.data.frame(gtf)
gtf_df <-
  gtf_df[,
         c("gene_id", "seqnames", "start", "end", "strand", "width", "gene_name",
           "gene_type")]


colnames(gtf_df) <- c("gene_id", "chrom", "start", "end", "strand", "length",
                      "gene_name", "gene_type")

gene_table <- 
  gtf_df[!duplicated(gtf_df[,c(1,7)]), ] %>% 
  dplyr::select('gene_id', 'gene_name', 'gene_type')

gene_table <- 
  rbind(gene_table, data.frame(gene_id=retro.annot$locus, 
                               gene_name=retro.annot$locus, 
                               gene_type=retro.annot$te_class))
rownames(gene_table) <- gene_table$gene_id

gene_table$display <- gene_table$gene_name
gene_table[duplicated(gene_table$gene_name), 'display'] <- 
  paste(gene_table[duplicated(gene_table$gene_name), 'display'], 
        gene_table[duplicated(gene_table$gene_name), 'gene_id'], sep='|')


############################## CLINICAL METADATA ###############################

# Download DLBCL NCICCR metadata
NCICCR_DLBCL_clinical_metadata <- 
  GDCquery_clinic("NCICCR-DLBCL", type = "clinical", save.csv = FALSE)

# Download DLBCL TCGA metadata
TCGA_DLBCL_clinical_metadata <- 
  GDCquery_clinic("TCGA-DLBC", type = "clinical", save.csv = FALSE)

# Import CCGI Burkitt Lymphoma metadata
BL_clinical_metadata <- 
  read.csv("metadata/BL/BL_samples_metadata.tsv", sep="\t")

# Import Follicular Lymphoma metadata
FL_clinical_metadata <- 
  read_excel("metadata/FL/CGCI_NHL_ClinicalDataSet_20110710.xlsx")

FL_SRA_run_table <- 
  read.csv("metadata/FL/FL_SRA.txt")

########################## DLBCL SUB-CLASSIFICATIONS ###########################

# Import paper-specific metadata for DLBCL classifications 

LymphGen <- 
  read_excel("metadata/DLBCL/Full_Wright_et_al_2021_LymphGen.xlsx", skip = 1)

EcoTyper <- 
  read_excel("metadata/DLBCL/Full_Steen_et_al_2021_EcoTyper.xlsx", sheet = "S2F", skip=1)

Holmes_scCOO <- 
  read_excel("metadata/DLBCL/Full_Holmes_et_al_2020.xlsx", skip=1)

########################## CREATE DLBCL METADATA FILES #########################

# Put TCGA and NCICCR metadata together
# Add the following metadata columns: ann_arbor_clinical_stage, age_at_diagnosis,
# gender, vital_status, days_to_last_follow_up, prior_malignancy, prior_treatment, 
# international_prognostic_index, last_known_disease_status, 
# progression_or_recurrence,year_of_diagnosis, year_of_death, tissue_or_organ_of_origin

DLBCL_clinical_metadata <- 
  plyr::rbind.fill(TCGA_DLBCL_clinical_metadata, NCICCR_DLBCL_clinical_metadata)
DLBCL_clinical_metadata <-
  DLBCL_clinical_metadata[,
                                             c("submitter_id",
                                               "gender", "vital_status",
                                               "ann_arbor_clinical_stage",
                                               "age_at_diagnosis",
                                               "international_prognostic_index",
                                               "days_to_last_follow_up",
                                               "tissue_or_organ_of_origin")]

# Filter metadata from all other studies so that we only have the NCI and/or 
# TCGA datasets together

EcoTyper <- EcoTyper[EcoTyper$Cohort == "Schmitz et al.",]
Holmes_scCOO <- Holmes_scCOO[Holmes_scCOO$Dataset == "NCI-DLBCL",]
LymphGen <- LymphGen[LymphGen$Cohort == "NCI",]
LymphGen <- LymphGen[!grepl("CTSP",LymphGen$Cohort),]

# Subset our original sample metadata sheet to only contain essentials
# We'll pull out updated metadata from the other files

DLBCL_metadata <- read.table("metadata/DLBCL/DLBCL_samples.tsv", sep = "\t", header = TRUE)
DLBCL_metadata <- DLBCL_metadata[, c("sample", "project", "case")]

# Add essential metadata from other files

DLBCL_metadata <- merge(DLBCL_metadata, 
                                 DLBCL_clinical_metadata, 
                                 by.x = "case", by.y = "submitter_id")

DLBCL_metadata$IPI_score <- LymphGen$`IPI\r\nScore`[match(DLBCL_metadata$case, 
                                                                   LymphGen$`Donor name`)]

DLBCL_metadata$RCHOP_like_chemo <- LymphGen$`R-CHOp-like\r\nChemo`[match(DLBCL_metadata$case, 
                                                                          LymphGen$`Donor name`)]

DLBCL_metadata$COO_class <- LymphGen$`COO\r\nClass`[match(DLBCL_metadata$case, 
                                                          LymphGen$`Donor name`)]

DLBCL_metadata$LymphGen_call <- LymphGen$`LymphGen\r\ncall`[match(DLBCL_metadata$case, 
                                                                  LymphGen$`Donor name`)]

DLBCL_metadata$DblHit_call <- LymphGen$`Dbl.Hit\r\nCall`[match(DLBCL_metadata$case, 
                                                               LymphGen$`Donor name`)]

DLBCL_metadata$EcoTyper_call <- 
  EcoTyper$`B cell state`[match(DLBCL_metadata$case, EcoTyper$`Sample ID`)]

DLBCL_metadata$Schmitz_call <- 
  Holmes_scCOO$`Genetic Subtype (Schmitz et al., NEJM 2018)`[match(DLBCL_metadata$case, 
                                                                   Holmes_scCOO$`Sample ID`)]

DLBCL_metadata$Chapuay_call <- 
  Holmes_scCOO$`Genetic Subtype (Chapuy et al., Nat. Medicine 2018)`[match(DLBCL_metadata$case, 
                                                                           Holmes_scCOO$`Sample ID`)]

DLBCL_metadata$scCOO_class_call <- 
  Holmes_scCOO$`sc-COO Class`[match(DLBCL_metadata$case, 
                                    Holmes_scCOO$`Sample ID`)]

DLBCL_metadata$scCOO_group_call <- 
  Holmes_scCOO$`sc-COO Group`[match(DLBCL_metadata$case,
                                    Holmes_scCOO$`Sample ID`)]

remove(EcoTyper, Holmes_scCOO, LymphGen, NCICCR_DLBCL_clinical_metadata,
       TCGA_DLBCL_clinical_metadata, DLBCL_clinical_metadata)

############################ CREATE FL METADATA FILES ##########################

# Import Follicular Lymphoma metadata
FL_clinical_metadata <- 
  FL_clinical_metadata[FL_clinical_metadata$`WHO diagnosis` %like% "FOLLICULAR",]

FL_clinical_metadata$sample_id <- FL_SRA_run_table$BioSample[match(FL_clinical_metadata$Patient_ID,
                                                                   FL_SRA_run_table$submitted_subject_id)]

FL_metadata <- FL_clinical_metadata

remove(FL_clinical_metadata, FL_SRA_run_table)

############################ CREATE BL METADATA FILES ##########################

# Import Follicular Lymphoma metadata

#BL_metadata <- BL_clinical_metadata[BL_clinical_metadata$pilot %like% "True",]

BL_metadata <- BL_clinical_metadata

BL_metadata <- BL_metadata[, c("case", "project_id", "submitter_id", 
                                        "sample_type", "sample", "tissue_type",
                                        "tumor_descriptor", "cohort", "clinical_variant",
                                        "ebv_status", "ebv_genome_type", "sex",
                                        "age_at_diagnosis", "anatomic_site_classification",
                                        "tissue_source_site")]

remove(BL_clinical_metadata)

######################## RENAME COLUMNS & MERGE METADATA #######################

colnames(BL_metadata) <- c("case", "project_id", "submitter_id", 
                           "sample_type", "sample", "tissue_type",
                           "tumor_descriptor", "cohort", "clinical_variant",
                           "ebv_status", "ebv_genome_type", "gender",
                           "age_at_diagnosis", "anatomic_site_classification",
                           "tissue_source_site")

colnames(FL_metadata) <- c("patient_id", "who_diagnosis", 
                           "days_to_birth_from_date_of_diagnosis",
                           "gender", "stage", "stage_group", "performance_status",
                           "LDH_ratio", "extranodal_sites", "tumor_size",
                           "IPI_score", "primary_treatment", 
                           "days_to_primary_treatment_start_from_date_of_diagnosis",
                           "response", "days_to_response_assessment_from_date_of_diagnosis",
                           "days_to_progression_from_date_of_diagnosis",
                           "site_progression", "secondary_treatment", "HSCT",
                           "days_to_HSCT_from_date_of_diagnosis", 
                           "days_to_transformed_lymphoma_diagnosed_from_date_of_diagnosis",
                           "days_to_last_follow_up_from_date_of_diagnosis", 
                           "status_lasrt_follow_up", "cause_of_death",
                           "cause_of_death_ICD10", "sample")

DLBCL_metadata <- DLBCL_metadata %>% remove_rownames %>% column_to_rownames(var="sample")
BL_metadata <- BL_metadata %>% remove_rownames %>% column_to_rownames(var="sample")
FL_metadata <- FL_metadata %>% remove_rownames %>% column_to_rownames(var="sample")
FL_metadata <- FL_metadata[!(row.names(FL_metadata) %in% "SAMN05182469"),]

# Add cancer type & COO

DLBCL_metadata$cancer_type <- "DLBCL"
DLBCL_metadata$subtype <- DLBCL_metadata$COO_class
FL_metadata$cancer_type <- "FL"
FL_metadata$subtype <- FL_metadata$who_diagnosis
BL_metadata$cancer_type <- "BL"
BL_metadata$subtype <- BL_metadata$clinical_variant

all_metadata <- rbind(
  DLBCL_metadata[, c("cancer_type", "subtype")],
  BL_metadata[, c("cancer_type", "subtype")],
  FL_metadata[, c("cancer_type", "subtype")]
)


################################ LOAD TELESCOPE ################################

# Old telescope reports

load_all_lymphoma_old <- function(df) {
  
  sample_names <- rownames(df)
  
  t_files <- file.path(paste("results/", df$cancer_type[1], "/telescope",
                             sep = ""), 
                       paste0(sample_names, '/', sample_names, 
                              '_telescope.report.tsv'))
  names(t_files) <- df$bulk_RNAseq
  
  counts.rtx <- load_telescope_reports(t_files, 
                                       all_locs=retro.hg38.v1$locus, 
                                       count_column = "count")
  
  assign(paste(df$cancer_type[1], "counts", "rtx", sep="."), 
         counts.rtx, 
         envir=.GlobalEnv)
  
}

load_all_lymphoma_old(DLBCL_metadata)

# New telescope reports

load_all_lymphoma_new <- function(df) {
  
  sample_names <- rownames(df)
  
  t_files <- file.path(paste("results/", df$cancer_type[1], "/telescope",
                             sep = ""), 
                       paste0(sample_names, '/', sample_names, 
                              '-telescope_report.tsv'))
  names(t_files) <- df$bulk_RNAseq
  
  counts.rtx <- load_telescope_reports(t_files, 
                                       all_locs=retro.hg38.v1$locus, 
                                       count_column = "final_count")
  
  assign(paste(df$cancer_type[1], "counts", "rtx", sep="."), 
         counts.rtx, 
         envir=.GlobalEnv)
  
  }

load_all_lymphoma_new(BL_metadata)
load_all_lymphoma_new(FL_metadata)

################################## LOAD STAR ###################################

DLBCL_files <- Sys.glob(file.path("results/DLBCL/star_alignment/*", '*.ReadsPerGene.out.tab'))
BL_files <- Sys.glob(file.path("results/BL/star_alignment/*", '*.ReadsPerGene.out.tab'))
FL_files <- Sys.glob(file.path("results/FL/star_alignment/*", '*.ReadsPerGene.out.tab'))

DLBCL.counts.tx <- load_star_counts(DLBCL_files)
BL.counts.tx <- load_star_counts(BL_files)
FL.counts.tx <- load_star_counts(FL_files)

################################# SANITY CHECK #################################

stopifnot(all(rownames(DLBCL.counts.rtx) == retro.hg38.v1$locus))
stopifnot(all(rownames(BL.counts.rtx) == retro.hg38.v1$locus))
stopifnot(all(rownames(FL.counts.rtx) == retro.hg38.v1$locus))

########################### ORDER SAMPLES / METADATA ###########################

# reorder counts.tx by metadata rowname
reorder_idx_counts.tx <- match(rownames(DLBCL_metadata), colnames(DLBCL.counts.tx))
DLBCL.counts.tx <- DLBCL.counts.tx[,reorder_idx_counts.tx]

reorder_idx_counts.tx <- match(rownames(BL_metadata), colnames(BL.counts.tx))
BL.counts.tx <- BL.counts.tx[,reorder_idx_counts.tx]

reorder_idx_counts.tx <- match(rownames(FL_metadata), colnames(FL.counts.tx))
FL.counts.tx <- FL.counts.tx[,reorder_idx_counts.tx]

# reorder counts.rtx by metadata rowname
reorder_idx_counts.rtx <- match(rownames(DLBCL_metadata), colnames(DLBCL.counts.rtx))
DLBCL.counts.rtx <- DLBCL.counts.rtx[,reorder_idx_counts.rtx]

reorder_idx_counts.rtx <- match(rownames(BL_metadata), colnames(BL.counts.rtx))
BL.counts.rtx <- BL.counts.rtx[,reorder_idx_counts.rtx]

reorder_idx_counts.rtx <- match(rownames(FL_metadata), colnames(FL.counts.rtx))
reorder_idx_counts.tx <- match(rownames(FL_metadata), colnames(FL.counts.tx))
FL.counts.rtx <- FL.counts.rtx[,reorder_idx_counts.rtx]
FL.counts.tx <- FL.counts.tx[,reorder_idx_counts.tx]

# sanity check
stopifnot(all(names(DLBCL.counts.tx) == names(DLBCL.counts.rtx)))
stopifnot(all(names(BL.counts.tx) == names(BL.counts.rtx)))
stopifnot(all(names(FL.counts.tx) == names(FL.counts.rtx)))

################################ COMBINE SAMPLES ###############################

# combine .tx and .rtx counts for all lymphoma samples

all.counts.rtx <- cbind(DLBCL.counts.rtx, BL.counts.rtx, FL.counts.rtx)
all.counts.tx <- cbind(DLBCL.counts.tx, BL.counts.tx, FL.counts.tx)

stopifnot(all(names(all.counts.tx) == names(all.counts.rtx)))

# combine .tx and .rtx counts in the same matrices

DLBCL.counts.comb <- rbind(DLBCL.counts.tx, DLBCL.counts.rtx)
BL.counts.comb <- rbind(BL.counts.tx, BL.counts.rtx)
FL.counts.comb <- rbind(FL.counts.tx, FL.counts.rtx)
all.counts.comb <- rbind(all.counts.tx, all.counts.rtx)

############################# SUBSET HERVs and L1s #############################

retro.hg38.v1 <- retro.hg38.v1 %>% remove_rownames %>% column_to_rownames(var="locus")

DLBCL.counts.herv <- DLBCL.counts.rtx[retro.hg38.v1$te_class == 'LTR',]
BL.counts.herv <- BL.counts.rtx[retro.hg38.v1$te_class == 'LTR',]
FL.counts.herv <- FL.counts.rtx[retro.hg38.v1$te_class == 'LTR',]
all.counts.herv <- all.counts.rtx[retro.hg38.v1$te_class == 'LTR',]

DLBCL.counts.l1 <- DLBCL.counts.rtx[retro.hg38.v1$te_class == 'LINE',]
BL.counts.l1 <- BL.counts.rtx[retro.hg38.v1$te_class == 'LINE',]
FL.counts.l1 <- FL.counts.rtx[retro.hg38.v1$te_class == 'LINE',]
all.counts.l1 <- all.counts.rtx[retro.hg38.v1$te_class == 'LINE',]

################################## SAVE FILES ##################################


save(all.counts.comb, all.counts.tx, all.counts.rtx, all.counts.herv, 
     all.counts.l1, all_metadata, file="r_outputs/01-all_lymphoma_counts.Rdata")

save(DLBCL.counts.comb, DLBCL.counts.tx, DLBCL.counts.rtx, DLBCL.counts.herv, 
     DLBCL.counts.l1, DLBCL_metadata, file="r_outputs/01-DLBCL_counts.Rdata")

save(BL.counts.comb, BL.counts.tx, BL.counts.rtx, BL.counts.herv, 
     BL.counts.l1, BL_metadata, file="r_outputs/01-BL_counts.Rdata")

save(FL.counts.comb, FL.counts.tx, FL.counts.rtx, FL.counts.herv, 
     FL.counts.l1, FL_metadata, file="r_outputs/01-FL_counts.Rdata")

save(all_metadata, DLBCL_metadata, BL_metadata, FL_metadata, 
     file="r_outputs/01-metadata.Rdata")

save(retro.hg38.v1, retro.annot, gene_table, 
     file="r_outputs/01-refs.Rdata")
