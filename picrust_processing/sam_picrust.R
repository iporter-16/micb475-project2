#### Nov 1, 2023 - SKA ####
#!/usr/bin/env Rscript

# load packages
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ggh4x)

### load abundance data ###
abundance_file <- "picrust_processing/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) %>% as.data.frame()
#abundance_data$pathway = rownames(abundance_data)
#abundance_data = abundance_data[,-1]
### load metadata ###
metadata <- read_delim("colombia/metadata_categorized_CL.txt",
                       delim = "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)

# samples = col_names()
# results_file_input <- ggpicrust2(file = abundance_file,metadata = metadata,group = "smoker", pathway = "MetaCyc", daa_method = "LinDA", ko_to_kegg = FALSE, order = "pathway_class", p_values_bar = TRUE, x_lab = "pathway_name")
# kegg_abundance <- ko2kegg_abundance("picrust_processing/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv") 

### Perform pathway differential abundance analysis (DAA) using LinDA method ###
daa_results_df <- pathway_daa(abundance = abundance_data %>% column_to_rownames("pathway"), 
                              metadata = metadata, 
                              group = "smoker", 
                              daa_method = "LinDA", 
                              select = NULL, reference = NULL) 
### Generate pathway heatmap ###
feature_with_p_0.05 <- daa_results_df %>% filter(p_values < 0.05)
heatmap_smoker <- pathway_heatmap(abundance = abundance_data %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "smoker")
# heatmap_smoker

### Generate pathway PCA plot ###
pca_plot_smoker <- pathway_pca(abundance = abundance_data %>% column_to_rownames("pathway"), metadata = metadata, group = "smoker")+ggtitle("PCA plot comparing smoking in all participants")
# pca_plot_smoker
pca_plot_LDL <- pathway_pca(abundance = abundance_data %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category")+ggtitle("PCA plot comparing LDL in all participants")
# pca_plot_LDL

# setwd("/Users/saman/Desktop/micb475-project2/pcoa_plots")
# ggsave("pca_smoker_all.png", pca_plot_smoker, scale=1)
# ggsave("pca_LDL_all.png", pca_plot_LDL, scale=1)

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("DESeq2", "edgeR") #"ALDEx2"
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data %>% column_to_rownames("pathway"), metadata = metadata, group = "smoker", daa_method = method)
})
