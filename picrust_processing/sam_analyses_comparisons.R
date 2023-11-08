#### Nov 7, 2023 - SKA ####
#!/usr/bin/env Rscript
pkgs <- c("ALDEx2", "edgeR", "DESeq2")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}
### Load packages
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ggh4x)
library(ALDEx2)
library(edgeR)
library(DESeq2)

### Read abundance data
abundance_file <- "picrust_processing/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) %>% as.data.frame()
### Read metadata
metadata <- read_delim("colombia/metadata_categorized_CL.txt",delim = "\t",escape_double = FALSE,trim_ws = TRUE)

### Isolate smoking and nonsmoking from abundance data and metadata
metadata_smoking <- filter(metadata, metadata$smoker=="Yes")
metadata_nonsmoking <- filter(metadata, metadata$smoker=="No")

abundance_data_smoking <- select(abundance_data, metadata_smoking$`#SampleID`)
abundance_data_smoking$pathway = abundance_data$pathway
abundance_data_smoking <- abundance_data_smoking %>% select(pathway, everything())

abundance_data_nonsmoking <- select(abundance_data, metadata_nonsmoking$`#SampleID`)
abundance_data_nonsmoking$pathway = abundance_data$pathway
abundance_data_nonsmoking <- abundance_data_nonsmoking %>% select(pathway, everything())

### Filter data to remove 0 standard deviation
row_sds_nonsmok <- abundance_data_nonsmoking %>% column_to_rownames("pathway") %>% apply(1, sd)
filtered_LDL_nonsmok <- abundance_data_nonsmoking[row_sds_nonsmok != 0, ]
row.names(filtered_LDL_nonsmok) <- NULL
row_sds_smok <- abundance_data_smoking %>% column_to_rownames("pathway") %>% apply(1, sd)
filtered_LDL_smok <- abundance_data_smoking[row_sds_smok != 0, ]
row.names(filtered_LDL_smok) <- NULL

### Perform pathway differential abundance analysis (DAA) ###
## LinDA method
daa_results_df <- pathway_daa(abundance = abundance_data %>% column_to_rownames("pathway"), 
                              metadata = metadata, 
                              group = "smoker", 
                              daa_method = "LinDA", 
                              select = NULL, reference = NULL)
daa_results_df_smoking_fibre <- pathway_daa(abundance = abundance_data_smoking%>% column_to_rownames("pathway"), metadata = metadata_smoking, group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_nonsmoking_fibre <- pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata_nonsmoking, group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL)
daa_results_df_smoking_LDL <- pathway_daa(abundance = abundance_data_smoking%>% column_to_rownames("pathway"), metadata = metadata_smoking, group = "LDL_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_nonsmoking_LDL <- pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata_nonsmoking, group = "LDL_category", daa_method = "LinDA", select = NULL, reference = NULL)
## Multiple methods: ALDEx2, DESeq2, edgeR
methods <- c("ALDEx2", "DESeq2", "edgeR") #"ALDEx2"
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data %>% column_to_rownames("pathway"), metadata = metadata, group = "smoker", daa_method = method)
}) # smoker only
daa_results_list_smoking_fibre <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data_smoking %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category", daa_method = method)
}) # smoking-fibre
daa_results_list_nonsmoking_fibre <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category", daa_method = method)
}) # nonsmoking-fibre
daa_results_list_smoking_LDL <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data_smoking %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category", daa_method = method)
}) # smoking-LDL
daa_results_list_nonsmoking_LDL <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category", daa_method = method)
}) # nonsmoking-LDL

# print(daa_results_list[[1]])

### Generate pathway heatmap using multiple methods###
## smokers only
sig_feature_list <- lapply(daa_results_list, function(results) {
  results %>% filter(p_values < 0.05)
})
heatmap_list_smoker <- lapply(sig_feature_list, function(sig) {
  pathway_heatmap(abundance = abundance_data %>% filter(pathway %in% sig$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "smoker")
})
heatmap_ALDEx2_smoker <- heatmap_list_smoker[[1]]
heatmap_DESeq2_smoker <- heatmap_list_smoker[[2]]
heatmap_edgeR_smoker <- heatmap_list_smoker[[3]]
# heatmap_ALDEx2_smoker
# heatmap_DESeq2_smoker
# heatmap_edgeR_smoker

## smoking-fibre
sig_feature_list_smoking_fibre <- lapply(daa_results_list_smoking_fibre, function(results) {
  results %>% filter(p_values < 0.05)
})
heatmap_list_smoking_fibre <- lapply(sig_feature_list_smoking_fibre, function(sig) {
  pathway_heatmap(abundance = abundance_data %>% filter(pathway %in% sig$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category")
})
heatmap_ALDEx2_smoking_fibre <- heatmap_list_smoking_fibre[[1]]
heatmap_DESeq2_smoking_fibre <- heatmap_list_smoking_fibre[[2]]
heatmap_edgeR_smoking_fibre <- heatmap_list_smoking_fibre[[3]]
heatmap_ALDEx2_smoking_fibre
heatmap_DESeq2_smoking_fibre
heatmap_edgeR_smoking_fibre

## nonsmoking-fibre
sig_feature_list_nonsmoking_fibre <- lapply(daa_results_list_nonsmoking_fibre, function(results) {
  results %>% filter(p_values < 0.05)
})
heatmap_list_nonsmoking_fibre <- lapply(sig_feature_list_nonsmoking_fibre, function(sig) {
  pathway_heatmap(abundance = abundance_data %>% filter(pathway %in% sig$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category")
})
heatmap_ALDEx2_nonsmoking_fibre <- heatmap_list_nonsmoking_fibre[[1]]
heatmap_DESeq2_nonsmoking_fibre <- heatmap_list_nonsmoking_fibre[[2]]
heatmap_edgeR_nonsmoking_fibre <- heatmap_list_nonsmoking_fibre[[3]]
heatmap_ALDEx2_nonsmoking_fibre
heatmap_DESeq2_nonsmoking_fibre
heatmap_edgeR_nonsmoking_fibre

## smoking-LDL
sig_feature_list_smoking_LDL <- lapply(daa_results_list_smoking_LDL, function(results) {
  results %>% filter(p_values < 0.05)
})
heatmap_list_smoking_LDL <- lapply(sig_feature_list_smoking_LDL, function(sig) {
  pathway_heatmap(abundance = abundance_data %>% filter(pathway %in% sig$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category")
})
heatmap_ALDEx2_smoking_LDL <- heatmap_list_smoking_LDL[[1]]
heatmap_DESeq2_smoking_LDL <- heatmap_list_smoking_LDL[[2]]
heatmap_edgeR_smoking_LDL <- heatmap_list_smoking_LDL[[3]]
heatmap_ALDEx2_smoking_LDL
heatmap_DESeq2_smoking_LDL
heatmap_edgeR_smoking_LDL

## nonsmoking-LDL
sig_feature_list_nonsmoking_LDL <- lapply(daa_results_list_nonsmoking_LDL, function(results) {
  results %>% filter(p_values < 0.05)
})
heatmap_list_nonsmoking_LDL <- lapply(sig_feature_list_nonsmoking_LDL, function(sig) {
  pathway_heatmap(abundance = abundance_data %>% filter(pathway %in% sig$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category")
})
heatmap_ALDEx2_nonsmoking_LDL <- heatmap_list_nonsmoking_LDL[[1]]
heatmap_DESeq2_nonsmoking_LDL <- heatmap_list_nonsmoking_LDL[[2]]
heatmap_edgeR_nonsmoking_LDL <- heatmap_list_nonsmoking_LDL[[3]]
heatmap_ALDEx2_nonsmoking_LDL
heatmap_DESeq2_nonsmoking_LDL
heatmap_edgeR_nonsmoking_LDL