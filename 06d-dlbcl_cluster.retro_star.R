################################################################################
################################################################################
################################################################################
################################################################################
############################# DLBCL CLUSTER STATS ############################## 

#################################### SETUP #####################################

library(tidyverse)
library(fossil)

################################### LOAD DATA ##################################

load("r_outputs/05o-DLBCL_pca_ccp_clusters_metadata.Rdata")

################################################################################
################################################################################
#################################### CLUST 7 ###################################
################################################################################
################################################################################
################################################################################

# Add clusters to metadata
DLBCL_metadata$clust.retro.k2 <- clust.df$clust.retro.k2
DLBCL_metadata$clust.retro.k3 <- clust.df$clust.retro.k3
DLBCL_metadata$clust.retro.k4 <- clust.df$clust.retro.k4
DLBCL_metadata$clust.retro.k5 <- clust.df$clust.retro.k5
DLBCL_metadata$clust.retro.k6 <- clust.df$clust.retro.k6
DLBCL_metadata$clust.retro.k7 <- clust.df$clust.retro.k7
DLBCL_metadata$clust.retro.k9 <- clust.df$clust.retro.k9

herv.vname <- names(clust.df)[grep('^clust', names(clust.df))]
alt.vname <- names(DLBCL_metadata)[grep('^clust', names(DLBCL_metadata))]
clin.vname <- c('vital_status', 'international_prognostic_index', 'gender',
                'age_at_diagnosis', 'tissue_or_organ_of_origin', 
                'RCHOP_like_chemo', "COO_class", "LymphGen_call",
                "DblHit_call", "EcoTyper_call", "Chapuay_call", 
                "scCOO_group_call")

catdat <- data.frame(
  clust.herv = clust.df$clust.retro.k7,
  DLBCL_metadata[,alt.vname],
  DLBCL_metadata[,clin.vname]
)

# Rand index
rand.df <- lapply(herv.vname, function(c1) {
  sapply(alt.vname, function(c2) {
    cass1 <- as.numeric(clust.df[,c1])
    cass2 <- as.numeric(DLBCL_metadata[,c2])
    fossil::adj.rand.index(cass1, cass2)
  })
}) %>% bind_rows %>% data.frame
rownames(rand.df) <- herv.vname


rand.df %>%
  rownames_to_column("rclust") %>%
  tidyr::separate(rclust, c("x1","x2","k"), sep="\\.") %>%
  select(-c(x1,x2)) %>%
  tidyr::pivot_longer(2:7) %>%
  ggplot(aes(x=k, y=value, group=name, color=name)) + geom_line() + 
  labs(title="Cluster assignment similarity", ylab="Rand index") + theme_minimal()


format(rand.df["clust.retro.k7",], digits=3)

ft <- lapply(alt.vname, function(acn) {
  fisher.test(table(catdat$clust.herv, catdat[,acn]), workspace=1e9, simulate.p.value=TRUE)
})
ft.p <- sapply(ft, function(fit) fit$p.value)
names(ft.p) <- alt.vname
format(data.frame(t(ft.p)), digits=3, scientific=T)

# gender

table(catdat$clust.herv, DLBCL_metadata$gender)
fisher.test(table(catdat$clust.herv, DLBCL_metadata$gender), workspace=1e7)

table(catdat$clust.herv, DLBCL_metadata$vital_status)
fisher.test(table(catdat$clust.herv, DLBCL_metadata$vital_status), workspace=1e9, simulate.p.value=TRUE)

table(catdat$clust.herv, DLBCL_metadata$ann_arbor_clinical_stage)
fisher.test(table(catdat$clust.herv, DLBCL_metadata$ann_arbor_clinical_stage), workspace=1e9, simulate.p.value=TRUE)

table(catdat$clust.herv, DLBCL_metadata$COO_class)
fisher.test(table(catdat$clust.herv, DLBCL_metadata$COO_class), workspace=1e9, simulate.p.value=TRUE)

table(catdat$clust.herv, DLBCL_metadata$LymphGen_call)
fisher.test(table(catdat$clust.herv, DLBCL_metadata$LymphGen_call), workspace=1e9, simulate.p.value=TRUE)

table(catdat$clust.herv, DLBCL_metadata$DblHit_call)
fisher.test(table(catdat$clust.herv, DLBCL_metadata$DblHit_call), workspace=1e9, simulate.p.value=TRUE)

table(catdat$clust.herv, DLBCL_metadata$Chapuay_call)
fisher.test(table(catdat$clust.herv, DLBCL_metadata$Chapuay_call), workspace=1e9, simulate.p.value=TRUE)

table(catdat$clust.herv, DLBCL_metadata$scCOO_group_call)
fisher.test(table(catdat$clust.herv, DLBCL_metadata$scCOO_group_call), workspace=1e9, simulate.p.value=TRUE)

