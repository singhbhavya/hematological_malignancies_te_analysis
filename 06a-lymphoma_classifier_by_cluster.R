################################################################################
################################################################################
################################################################################
################################################################################
########################### LYMPHOMA CLASSIFICATION ############################

## Plan:
## Prelim lymphoma classification, based on previous work from Holmes et al, 2020.

#################################### SETUP #####################################

library(knitr)
library(tidyverse)
library(matrixStats)
library(data.table)
library(DESeq2)
library(ggsci)
library(edgeR)
library(ashr)
library(cowplot)
library(readxl)
library(boot)
library(boot.pval)

################################### LOAD DATA ##################################

load("r_outputs/01-refs.Rdata")
load("r_outputs/02-BL_filt_counts.Rdata")
load("r_outputs/05b-BL_pca_dds.Rdata")
load("r_outputs/04-gcb.norm.merged.azimuth.Rdata")
load("r_outputs/02-FL_filt_counts.Rdata")
load("r_outputs/05b-FL_pca_dds.Rdata")
load("r_outputs/02-DLBCL_filt_counts.Rdata")
load("r_outputs/05b-DLBCL_pca_dds.Rdata")

################################# BL Z SCORES ##################################

# Change the vst to z-scores
BL.vsd <- assay(BL.g.tform)
BL.z.tform <- t(scale(t(BL.vsd)))

# Change all gene IDs to gene names
rownames(BL.z.tform) <- gene_table[rownames(BL.z.tform), 'display']

# Sanity check
all(rownames(BL.z.tform) %in% gene_table$display)

################################# FL Z SCORES ##################################

# Change the vst to z-scores
FL.vsd <- assay(FL.g.tform)
FL.z.tform <- t(scale(t(FL.vsd)))

# Change all gene IDs to gene names
rownames(FL.z.tform) <- gene_table[rownames(FL.z.tform), 'display']

# Sanity check
all(rownames(FL.z.tform) %in% gene_table$display)

############################### DLBCL Z SCORES #################################

# Change the vst to z-scores
DLBCL.vsd <- assay(DLBCL.g.tform)
DLBCL.z.tform <- t(scale(t(DLBCL.vsd)))

# Change all gene IDs to gene names
rownames(DLBCL.z.tform) <- gene_table[rownames(DLBCL.z.tform), 'display']

# Sanity check
all(rownames(DLBCL.z.tform) %in% gene_table$display)


############################### CLASSIFICATION #################################

classify_lymphomas <- 
  
  # Function takes as input a file of Z scores, and a file of pre-defined clusters,
  # which should have a list of upregulated and downregulated genes per cluster.
  
  function(z_scores_file, cluster_file, count_file) {
    
    # Create an empty dataframe where values will be stored later 
    classifications <<- data.frame(matrix(ncol = 6, nrow = 0))
    cols_classifications <<- c("sample", "cluster", "score", "p_val", "ci_l",
                                     "ci_u")
    colnames(classifications) <<- cols_classifications
    
    # Convert z_file to a dataframe
    z_file <- as.data.frame(z_scores_file)
    
    # For loop to iterate over each sample (samp ID derived from z_scores_file)
    for (samps in 1:length(colnames(z_scores_file))) {
      
      # Assign sample ID name
      sample_id <- colnames(z_file)[samps]
      
      # Print update
      message(paste0("Currently working on ", sample_id))
      
      # Subset z_scores_file to one sample 
      sample_data <- z_file[c(sample_id)]
      
      # Add z-scores for the sample to the cluster_file, match by gene
      cluster_file$samp_z <-  sample_data[[`sample_id`]][
        match(cluster_file$gene, rownames(sample_data))]
      
      # For each z-score, compare to the average z-score, and assign a sample score 
      # of either -1 or 1 depending on whether the two z-scores are going in the 
      # same direction
      cluster_file$samp_score <- 
        ifelse(
          cluster_file$avg_log2FC > 0 & cluster_file$samp_z > 0, 1, 
          ifelse(
            cluster_file$avg_log2FC > 0 & cluster_file$samp_z < 0, -1,
            ifelse(
              cluster_file$avg_log2FC < 0 & cluster_file$samp_z < 0, 1,
              ifelse(
                cluster_file$avg_log2FC < 0 & cluster_file$samp_z > 0, -1, 
                ifelse(
                  cluster_file$avg_log2FC == 0 & cluster_file$samp_z == 0, 1,
                  ifelse(
                    cluster_file$avg_log2FC & cluster_file$samp_z > 1 | cluster_file$samp_z < 1, -1, 
                    "Error")))))) 
      
      # For loop to iterate over each cluster in the cluster dataframe. Basically,
      # for each cluster, we want to 1) calculate the normalzed sum (which I think
      # is the mean?), 2) bootstrap the sample score values
      
      
      
      for (clust in 1:length(unique(cluster_file$cluster))) { 
        
        # Get cluster name
        clust_name <- unique(cluster_file$cluster)[clust]
        
        # Print message
        message(paste0("Currently working on ", clust_name))
        
        # Subset the cluster_file to have just the clusters and sample scores
        cluster_score_df <- cluster_file[,c("cluster", "samp_score")]
        
        # Subset the cluster_file for the cluster that the loop is currently on.
        cluster_score_df <- subset(cluster_file, cluster==clust_name)
        
        # Convert to dataframe
        cluster_score_df <- as.data.frame(cluster_score_df[,c("gene", "samp_score")])
        
        # Create function to get the normalized sum (mean?) of the sample scores
        # per cluster. Remove NAs (genes that were not expressed in this sample)
        my.mean = function(x, indices) {return( mean(x[indices], na.rm=TRUE) ) }
        
        # Bootstrap the scores for each sample 10,000 times. Calculate mean.
        score_boot <- boot(cluster_score_df$samp_score, my.mean, 10000)
        
        # Confidence interval for normalized sum.
        score_ci <- boot.ci(score_boot)
        
        # p-value for normalized sum. 
        score_pval <-boot.pval(score_boot)
        
        # Add row to the classifications dataframe with sample name, 
        # clusters, scores for each cluster, p value, and confidence interval.
        classifications[nrow(classifications) + 1,] <<-
          c(sample_id, clust_name, as.numeric(score_boot$t0), as.numeric(score_pval), 
            as.numeric(score_ci$normal[2]), as.numeric(score_ci$normal[3]))
        
      }
      
    }
    
    # Convert to numeric
    classifications[, 3:6] <<- sapply(classifications[, 3:6], as.numeric)
    
    # Subset cluster assignments that are higher than an absolute score of 0.3
    classifications_filtered <<- classifications[which(classifications[,3]>0.3),]
    
    # pick scores for each sample with highest p values
    final_classes <<- unique(setDT(classifications_filtered)[order(p_val)], 
                                   by = "sample")
    
    
  }

################################ CLASSIFY BL  ##################################

# Run classifier on BL
classify_lymphomas(z_scores_file = BL.z.tform,
                   cluster_file = l2.markers,
                   count_file = BL.counts.mfilt.comb)

# Assign dfs
bl_clusters <- final_classes
bl_classifications <- classifications
bl_classifications_filtered <- classifications_filtered

# Assign cluster names

bl_clusters$cluster_name <- 
  unique(l2.markers$cluster)[as.numeric(bl_clusters$cluster)]

# Save files

save(bl_clusters, bl_classifications, bl_classifications_filtered,
     file="r_outputs/06-BL_clusters.Rdata")

################################ CLASSIFY FL  ##################################

# Run classifier on FL
classify_lymphomas(z_scores_file = FL.z.tform,
                   cluster_file = l2.markers,
                   count_file = FL.counts.mfilt.comb)

# Assign dfs
fl_clusters <- final_classes
fl_classifications <- classifications
fl_classifications_filtered <- classifications_filtered

# Assign cluster names

fl_clusters$cluster_name <- 
  unique(l2.markers$cluster)[as.numeric(fl_clusters$cluster)]

# Save files

save(fl_clusters, fl_classifications, fl_classifications_filtered,
     file="r_outputs/06-FL_clusters.Rdata")

############################## CLASSIFY DLBCL  #################################

# Run classifier on FL
classify_lymphomas(z_scores_file = DLBCL.z.tform,
                   cluster_file = l2.markers,
                   count_file = DLBCL.counts.mfilt.comb)

# Assign dfs
dlbcl_clusters <- final_classes
dlbcl_classifications <- classifications
dlbcl_classifications_filtered <- classifications_filtered

# Assign cluster names

dlbcl_clusters$cluster_name <- 
  unique(l2.markers$cluster)[as.numeric(dlbcl_clusters$cluster)]

# Save files

save(dlbcl_clusters, dlbcl_classifications, dlbcl_classifications_filtered,
     file="r_outputs/06-DLBCL_clusters.Rdata")
