################################################################################
################################################################################
################################################################################
################################################################################
############################### BL CLUSTER STATS ############################### 

#################################### SETUP #####################################

library(tidyverse)
library(fossil)

################################### LOAD DATA ##################################

load("r_outputs/05i-BL_pca_ccp_clusters_metadata.Rdata")

################################################################################
################################################################################
#################################### CLUST 7 ###################################
################################################################################
################################################################################
################################################################################

herv.vname <- names(clust.df)[grep('^clust', names(clust.df))]
alt.vname <- names(BL_metadata)[grep('^clust', names(BL_metadata))]
clin.vname <- c('clinical_variant', 'ebv_status', 'ebv_genome_type',
                'gender', 'age_at_diagnosis', 
                'anatomic_site_classification', "tissue_source_site", "tumor_biopsy",
                "MYC_SV", "MYC_SV_partner", "subgroup", 
                "subtype")

catdat <- data.frame(
  clust.herv = clust.df$clust.retro.k2,
  BL_metadata[,alt.vname],
  BL_metadata[,clin.vname]
)

# Rand index
rand.df <- lapply(herv.vname, function(c1) {
  sapply(alt.vname, function(c2) {
    cass1 <- as.numeric(clust.df[,c1])
    cass2 <- as.numeric(BL_metadata[,c2])
    fossil::adj.rand.index(cass1, cass2)
  })
}) %>% bind_rows %>% data.frame
rownames(rand.df) <- herv.vname


rand.df %>%
  rownames_to_column("rclust") %>%
  tidyr::separate(rclust, c("x1","x2","k"), sep="\\.") %>%
  dplyr::select(-c(x1,x2)) %>%
  tidyr::pivot_longer(2:4) %>%
  ggplot(aes(x=k, y=value, group=name, color=name)) + geom_line() + 
  labs(title="Cluster assignment similarity", ylab="Rand index") + theme_minimal()

format(rand.df["clust.retro.k2",], digits=3)

ft <- lapply(alt.vname, function(acn) {
  fisher.test(table(catdat$clust.herv, catdat[,acn]), workspace=1e9, simulate.p.value=TRUE)
})
ft.p <- sapply(ft, function(fit) fit$p.value)
names(ft.p) <- alt.vname
format(data.frame(t(ft.p)), digits=3, scientific=T)

table(catdat$clust.herv, BL_metadata$gender)
fisher.test(table(catdat$clust.herv, BL_metadata$gender), workspace=1e7)

table(catdat$clust.herv, BL_metadata$clinical_variant)
fisher.test(table(catdat$clust.herv, BL_metadata$clinical_variant), workspace=1e7)

table(catdat$clust.herv, BL_metadata$ebv_status)
fisher.test(table(catdat$clust.herv, BL_metadata$ebv_status), workspace=1e7)

table(catdat$clust.herv, BL_metadata$MYC_SV_partner)
fisher.test(table(catdat$clust.herv, BL_metadata$MYC_SV_partner), workspace=1e7)

table(catdat$clust.herv, BL_metadata$ebv_genome_type)
fisher.test(table(catdat$clust.herv, BL_metadata$ebv_genome_type), workspace=1e7)


