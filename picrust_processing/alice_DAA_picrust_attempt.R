#### Nov 06, 2023 - AW ####

#Attempt at differential expression analysis (LinDa, ALDEx2, DESeq2, edgeR)

#install & load packages
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library("ggh4x")
library(ALDEx2)
library(DESeq2)
library(dplyr)

#read abundance and metadata
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


#Filter smoking and nonsmoking from metadata
metadata_smoking <- filter(metadata, metadata$smoker=="Yes")
metadata_nonsmoking <- filter(metadata, metadata$smoker=="No")

#Select smoking and nonsmoking from abundance data
abundance_data_smoking <- select(abundance_data, metadata_smoking$`#SampleID`)
abundance_data_smoking$pathway = abundance_data$pathway
abundance_data_smoking <- abundance_data_smoking %>% select(pathway, everything())

abundance_data_nonsmoking <- select(abundance_data, metadata_nonsmoking$`#SampleID`)
abundance_data_nonsmoking$pathway = abundance_data$pathway
abundance_data_nonsmoking <- abundance_data_nonsmoking %>% select(pathway, everything())

# Filter data to remove 0 SD
row_sds_nonsmok <- abundance_data_nonsmoking %>% column_to_rownames("pathway") %>% apply(1, sd)
filtered_LDL_nonsmok <- abundance_data_nonsmoking[row_sds_nonsmok != 0, ]
row.names(filtered_LDL_nonsmok) <- NULL
row_sds_smok <- abundance_data_smoking %>% column_to_rownames("pathway") %>% apply(1, sd)
filtered_LDL_smok <- abundance_data_smoking[row_sds_smok != 0, ]
row.names(filtered_LDL_smok) <- NULL

### Perform pathway differential abundance analysis (DAA) ###
# LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
daa_results_df <- pathway_daa(abundance = abundance_data %>% column_to_rownames("pathway"), 
                              metadata = metadata, 
                              group = "smoker", 
                              daa_method = "LinDA", 
                              select = NULL, reference = NULL)

daa_results_df_smoking_LDL <- pathway_daa(abundance = abundance_data_smoking%>% column_to_rownames("pathway"), metadata = metadata_smoking, group = "LDL_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_nonsmoking_LDL <- pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata_nonsmoking, group = "LDL_category", daa_method = "LinDA", select = NULL, reference = NULL)

daa_results_df_smoking_fibre <- pathway_daa(abundance = abundance_data_smoking%>% column_to_rownames("pathway"), metadata = metadata_smoking, group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_nonsmoking_fibre <- pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata_nonsmoking, group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL)

# Run pathway DAA for other methods (ALDEx2, DESeq2, edgeR)
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR")

daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data %>% column_to_rownames("pathway"), metadata = metadata, group = "smoker", daa_method = method)
})

#LDL
daa_results_list_smoking_LDL <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data_smoking %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category", daa_method = method)
})

daa_results_list_nonsmoking_LDL <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category", daa_method = method)
})

#Fibre
daa_results_list_smoking_fibre <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data_smoking %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category", daa_method = method)
})

daa_results_list_nonsmoking_fibre <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category", daa_method = method)
})


