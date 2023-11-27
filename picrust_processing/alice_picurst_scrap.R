#S

library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library("ggh4x")

abundance_file <- "picrust_processing/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
abundance_data  =as.data.frame(abundance_data)

#rownames(abundance_data) = abundance_data$pathway
#abundance_data = abundance_data[,-1]

metadata <- read_delim(
  "colombia/metadata_categorized_CL.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

#select smoking and nonsmoking from abundance data and metadata
metadata_nonsmoking <- filter(metadata, metadata$smoker=="No")
abundance_data_nonsmoking <- select(abundance_data, metadata_nonsmoking$`#SampleID`)
abundance_data_nonsmoking$pathway <- abundance_data$pathway
#abundance_data_nonsmoking <- tibble::rownames_to_column(abundance_data_nonsmoking, var = "pathway")
#abundance_data_nonsmoking$pathway = rownames(abundance_data_nonsmoking)
#rownames(abundance_data_nonsmoking) = NULL


# Perform pathway DAA using LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
#daa_results_df <- pathway_daa(abundance = abundance_data_smoking, metadata = metadata, group = "LDL_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_nonsmoking <- pathway_daa(abundance = abundance_data_nonsmoking %>% 
                                           column_to_rownames("pathway"), metadata = metadata_nonsmoking, group = "LDL_category", daa_method = "LinDA", select = NULL, reference = NULL)

feature_with_p_0.05 <- daa_results_df_nonsmoking %>% filter(p_values < 0.05)
pathway_heatmap(abundance = abundance_data_nonsmoking %>% filter(pathway %in% feature_with_p_0.05$feature) %>%
                  column_to_rownames("pathway"), metadata = metadata, group = "LDL_category") + ggtitle("Significant Differential expression in non-smokers with high/low LDL")

# Generate pathway PCA plot 
pathway_pca(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata_nonsmoking, group = "LDL_category")



