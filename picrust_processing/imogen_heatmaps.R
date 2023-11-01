#From picrust2 analysis, produces heatmaps of pathways signififantly (p<0.05 up or downregulated)
#in smoking and nonsmoking adults, displaying their fiber consumption (high or low)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

#read abundance data
abundance_file <- "picrust_processing/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) %>% as.data.frame()
#read metadata
metadata <- read_delim("colombia/metadata_categorized_CL.txt",delim = "\t",escape_double = FALSE,trim_ws = TRUE)

#isolate smoking and nonsmoking from abundance data and metadata
metadata_smoking <- filter(metadata, metadata$smoker=="Yes")
metadata_nonsmoking <- filter(metadata, metadata$smoker=="No")
abundance_data_smoking <- select(abundance_data, metadata_smoking$`#SampleID`)
abundance_data_smoking$pathway = rownames(abundance_data_smoking)
rownames(abundance_data_smoking) = NULL

abundance_data_nonsmoking <- select(abundance_data, metadata_nonsmoking$`#SampleID`)
abundance_data_nonsmoking$pathway = rownames(abundance_data_nonsmoking)
rownames(abundance_data_nonsmoking) = NULL

# Perform pathway differential abundance analysis (DAA) using LinDA
# This requires 
daa_results_df <- pathway_daa(abundance = abundance_data[,-1], metadata = metadata, 
                              group = "smoker", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_smoking <- pathway_daa(abundance = abundance_data_smoking[,-59], metadata = metadata_smoking, 
                                      group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_nonsmoking <- pathway_daa(abundance = abundance_data_nonsmoking[,-384], metadata = metadata_nonsmoking, 
                                         group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL)
# Generate pathway heatmap
feature_with_p_0.05_smok <- daa_results_df_smoking %>% filter(p_values < 0.05)
heat_smok <- pathway_heatmap(abundance = abundance_data_smoking %>% filter(pathway %in% feature_with_p_0.05_smok$feature) %>% 
                  column_to_rownames("pathway"), metadata = metadata, group = "fibre_category")

feature_with_p_0.05_nonsmok <- daa_results_df_nonsmoking %>% filter(p_values < 0.05)
heat_nonsmok <- pathway_heatmap(abundance = abundance_data_nonsmoking %>% filter(pathway %in% feature_with_p_0.05_nonsmok$feature) %>% 
                               column_to_rownames("pathway"), metadata = metadata, group = "fibre_category")
