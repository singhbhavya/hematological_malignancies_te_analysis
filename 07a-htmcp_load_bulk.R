################################################################################
################################################################################
################################################################################
################################################################################
##########################  HTMCP DATA / METADATA SETUP ########################


#################################### SETUP #####################################

library(tidyverse)
library(readxl)
library(GenomicDataCommons)
library(TCGAbiolinks)
library(dplyr)
library(rtracklayer)
library(data.table)
library(scopetools)
library(tximport)

################################### LOAD DATA ##################################

load("r_outputs/01-refs.Rdata")
load("r_outputs/01-all_lymphoma_counts.Rdata")
load("r_outputs/01-DLBCL_counts.Rdata")

############################ LOAD GENE ANNOTATIONS #############################

ttg <- read.table('refs/ttg.tsv', 
                  sep = '\t',
                  header = T,
                  stringsAsFactors = F
)

gsym <- read.table('refs/gsym.tsv', 
                   sep = '\t',
                   header = T,
                   stringsAsFactors = F
)


HTMCP_metadata <- read.csv('metadata/HTMCP/samples_hiv_dlbcl.csv', stringsAsFactors = F) %>%
  dplyr::select(sample_id=BioSample, primary_site=body_site) %>%
  dplyr::mutate(
    primary_site=factor(recode(primary_site, `lymph node`='LN', `blood`='BL'),
                        levels=c('LN','BL')),
    HIV=factor("positive", levels=c('negative', 'positive'))
  )

row.names(HTMCP_metadata) <- HTMCP_metadata$sample_id

################################## LOAD HTMCP ##################################

files <- file.path('results/HTMCP/kallisto_alignment', rownames(HTMCP_metadata), 
                   paste0(rownames(HTMCP_metadata), ".", 'kallisto.abundance.tsv'))
names(files) <- rownames(HTMCP_metadata)
txi <- tximport(files, type = "kallisto", tx2gene=ttg)

################################ LOAD TELESCOPE  ###############################
#--- Load Telescope counts
files <- file.path('results/HTMCP/telescope', rownames(HTMCP_metadata),  
                   paste0(rownames(HTMCP_metadata), "-", 'telescope_report.tsv'))
names(files) <- rownames(HTMCP_metadata)

cts <- lapply(1:length(files), 
              function(i) {
                tmp <- read.table(files[i],
                                  sep='\t', header=T, stringsAsFactors=F)
                ret <- data.frame(
                  transcript=retro.annot$locus, 
                  stringsAsFactors = F
                ) %>%
                  left_join(tmp, by='transcript') %>%
                  dplyr::select(gene_id=transcript, count=final_count)
                
                ret[is.na(ret)] <- 0
                stopifnot(all(ret$gene_id == retro.annot$locus))
                ret$gene_id <- NULL
                names(ret) <- c(names(files)[i])
                ret
              }) %>% bind_cols

row.names(cts) <- retro.annot$locus
rxi <- list(counts=cts)

################################ COMBINE SAMPLES ###############################

# all counts
stopifnot(all(rownames(cts) == rownames(all.counts.rtx)))
all.counts.hiv.rtx <- cbind(cts, all.counts.rtx)

# dlbcl counts
stopifnot(all(rownames(cts) == rownames(DLBCL.counts.rtx)))
DLBCL.hiv.counts.rtx <- cbind(cts, DLBCL.counts.rtx)

############################# SUBSET HERVs and L1s #############################

# all counts
all.counts.hiv.herv <- all.counts.hiv.rtx[retro.hg38.v1$te_class == 'LTR',]
all.counts.hiv.l1 <- all.counts.hiv.rtx[retro.hg38.v1$te_class == 'LINE',]

# dlbcl counts
DLBCL.counts.hiv.herv <- DLBCL.hiv.counts.rtx[retro.hg38.v1$te_class == 'LTR',]
DLBCL.counts.hiv.l1 <- DLBCL.hiv.counts.rtx[retro.hg38.v1$te_class == 'LINE',]

################################ COMBINE METADATA ##############################

# all counts
HTMCP_metadata_formatted <- HTMCP_metadata
HTMCP_metadata_formatted$cancer_type <- "DLBCL"
HTMCP_metadata_formatted$subtype <- "HIV-positive"


all_metadata_hiv <- 
  rbind(
  subset(HTMCP_metadata_formatted, select = c("cancer_type", "subtype")), 
  all_metadata
)


# dlbcl
DLBCL_metadata_hiv <- dplyr::bind_rows(
  subset(HTMCP_metadata_formatted, select = c("cancer_type", "subtype")),
  DLBCL_metadata)
  
##########################@###### SANITY CHECK ##########@###################### 

stopifnot(all(names(DLBCL.hiv.counts.rtx) == rownames(DLBCL_metadata_hiv)))
stopifnot(all(names(all.counts.hiv.rtx) == rownames(all_metadata_hiv)))

################################## SAVE FILES ##################################

save(all.counts.hiv.rtx, all.counts.hiv.herv, all.counts.hiv.l1, all_metadata_hiv, 
     file="r_outputs/07-htmcp_all_lymphoma_counts.Rdata")

save(DLBCL.hiv.counts.rtx, DLBCL.counts.hiv.herv, DLBCL.counts.hiv.l1, DLBCL_metadata_hiv, 
     file="r_outputs/07-htmcp_dlbcl_counts.Rdata")
