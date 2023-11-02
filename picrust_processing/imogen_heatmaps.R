#From picrust2 analysis, produces heatmaps of pathways signififantly (p<0.05 up or downregulated)
#in smoking and nonsmoking adults, displaying their fiber consumption (high or low)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ggh4x)

#read abundance data
abundance_file <- "picrust_processing/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) %>% as.data.frame()
#read metadata
metadata <- read_delim("colombia/metadata_categorized_CL.txt",delim = "\t",escape_double = FALSE,trim_ws = TRUE)

#isolate smoking and nonsmoking from abundance data and metadata
metadata_smoking <- filter(metadata, metadata$smoker=="Yes")
metadata_nonsmoking <- filter(metadata, metadata$smoker=="No")

abundance_data_smoking <- select(abundance_data, metadata_smoking$`#SampleID`)
abundance_data_smoking$pathway = abundance_data$pathway
abundance_data_smoking <- abundance_data_smoking %>% select(pathway, everything())

abundance_data_nonsmoking <- select(abundance_data, metadata_nonsmoking$`#SampleID`)
abundance_data_nonsmoking$pathway = abundance_data$pathway
abundance_data_nonsmoking <- abundance_data_nonsmoking %>% select(pathway, everything())

#Run PCA plots
# pathway_pca(abundance = abundance_data %>% column_to_rownames("pathway"), metadata = metadata, group = "smoker")
# pathway_pca(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata_nonsmoking, group = "LDL_category")
# pathway_pca(abundance = abundance_data_smoking %>% column_to_rownames("pathway"), metadata = metadata_smoking, group = "LDL_category")

# Perform pathway differential abundance analysis (DAA) using LinDA
daa_results_df <- pathway_daa(abundance = abundance_data %>% column_to_rownames("pathway"), metadata = metadata, 
                              group = "smoker", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_smoking <- pathway_daa(abundance = abundance_data_smoking%>% column_to_rownames("pathway"), metadata = metadata_smoking, 
                                      group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_nonsmoking <- pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata_nonsmoking, 
                                         group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL)
# Generate pathway heatmap
feature_with_p_0.05_smok <- daa_results_df_smoking %>% filter(p_values < 0.05)
heat_smok <- pathway_heatmap(abundance = abundance_data_smoking %>% filter(pathway %in% feature_with_p_0.05_smok$feature) %>% 
                  column_to_rownames("pathway"), metadata = metadata, group = "fibre_category") + 
                  ggtitle("Significant differential expression in smokers")

feature_with_p_0.05_nonsmok <- daa_results_df_nonsmoking %>% filter(p_values < 0.05)
heat_nonsmok <- pathway_heatmap(abundance = abundance_data_nonsmoking %>% filter(pathway %in% feature_with_p_0.05_nonsmok$feature) %>% 
                  column_to_rownames("pathway"), metadata = metadata, group = "fibre_category") + 
                  ggtitle("Significant differential expression in nonsmokers")

### LDL ###
# Perform pathway differential abundance analysis (DAA) using LinDA
daa_results_df <- pathway_daa(abundance = abundance_data %>% column_to_rownames("pathway"), metadata = metadata, 
                              group = "smoker", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_smoking_LDL <- pathway_daa(abundance = abundance_data_smoking%>% column_to_rownames("pathway"), metadata = metadata_smoking, 
                                      group = "LDL_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_nonsmoking_LDL <- pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata_nonsmoking, 
                                         group = "LDL_category", daa_method = "LinDA", select = NULL, reference = NULL)
# Generate pathway heatmap for LDL groups
feature_with_p_0.05_smok_LDL <- daa_results_df_smoking_LDL %>% filter(p_values < 0.05)
heat_smok_LDL <- pathway_heatmap(abundance = abundance_data_smoking %>% filter(pathway %in% feature_with_p_0.05_smok$feature) %>% 
                  column_to_rownames("pathway"), metadata = metadata, group = "LDL_category") + 
                  ggtitle("Significant differential expression in smokers with high/low LDL")

feature_with_p_0.05_nonsmok_LDL <- daa_results_df_nonsmoking_LDL %>% filter(p_values < 0.05)
heat_nonsmok_LDL <- pathway_heatmap(abundance = abundance_data_nonsmoking %>% filter(pathway %in% feature_with_p_0.05_nonsmok$feature) %>% 
                   column_to_rownames("pathway"), metadata = metadata, group = "LDL_category") + 
                  ggtitle("Significant differential expression in nonsmokers with high/low LDL")

ggsave("/heatmaps/heatmap_smok_fibre.png",heat_smok)
