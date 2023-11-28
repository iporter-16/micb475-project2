#### Nov 7, 2023 - SKA ####
#!/usr/bin/env Rscript

### Load packages
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ggh4x)
# library(ALDEx2)
# library(edgeR)
library(DESeq2)

### Read abundance data
abundance_file <- "picrust_processing/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) %>% as.data.frame()
### Read metadata
metadata <- read_delim("colombia/metadata_categorized_CL.txt",delim = "\t",escape_double = FALSE,trim_ws = TRUE)

### Isolate smoking and nonsmoking from abundance data and metadata
metadata_smoking <- filter(metadata, metadata$smoker=="Yes")
metadata_nonsmoking <- filter(metadata, metadata$smoker=="No")

abundance_data_smoking <- dplyr::select(abundance_data, metadata_smoking$`#SampleID`) # if MASS package is installed, it will mask select() from dplyr - using dplyr:: forces R to use select() directly from dplyr package
abundance_data_smoking$pathway = abundance_data$pathway
abundance_data_smoking <- abundance_data_smoking %>% dplyr::select(pathway, everything())

abundance_data_nonsmoking <- dplyr::select(abundance_data, metadata_nonsmoking$`#SampleID`)
abundance_data_nonsmoking$pathway = abundance_data$pathway
abundance_data_nonsmoking <- abundance_data_nonsmoking %>% dplyr::select(pathway, everything())

### Filter data to remove 0 standard deviation
row_sds_nonsmok <- abundance_data_nonsmoking %>% column_to_rownames("pathway") %>% apply(1, sd)
filtered_LDL_nonsmok <- abundance_data_nonsmoking[row_sds_nonsmok != 0, ]
row.names(filtered_LDL_nonsmok) <- NULL
row_sds_smok <- abundance_data_smoking %>% column_to_rownames("pathway") %>% apply(1, sd)
filtered_LDL_smok <- abundance_data_smoking[row_sds_smok != 0, ]
row.names(filtered_LDL_smok) <- NULL

### Perform pathway differential abundance analysis (DAA) ### # tbd: combine LinDA method with other methods
## Multiple methods: DESeq2, edgeR, LinDA, ALDEx2
methods <- c("DESeq2", "edgeR", "LinDA") #"ALDEx2", <- error when running ALDEx2, worked initially but not anymore
daa_results_list <- lapply(methods, function(method) {
  pathway_daa(abundance = abundance_data %>% column_to_rownames("pathway"), metadata = metadata, group = "smoker", daa_method = method)
}) # smoker only
# daa_results_list_smoking_fibre <- lapply(methods, function(method) {pathway_daa(abundance = abundance_data_smoking %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category", daa_method = method)
# }) # smoking-fibre
# daa_results_list_nonsmoking_fibre <- lapply(methods, function(method) {pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category", daa_method = method)
# }) # nonsmoking-fibre
daa_results_list_smoking_LDL <- lapply(methods, function(method) {pathway_daa(abundance = abundance_data_smoking %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category", daa_method = method)
}) # smoking-LDL
daa_results_list_nonsmoking_LDL <- lapply(methods, function(method) {pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category", daa_method = method)
}) # nonsmoking-LDL

### Significnat Features, p<0.05 ###
sig_feature_list <- lapply(daa_results_list, function(results) {
  results %>% filter(p_values < 0.05)
}) # smokers only
# sig_feature_list_smoking_fibre <- lapply(daa_results_list_smoking_fibre, function(results) {results %>% filter(p_values < 0.05)
# }) # smoking-fibre
# sig_feature_list_nonsmoking_fibre <- lapply(daa_results_list_nonsmoking_fibre, function(results) {results %>% filter(p_values < 0.05)
# }) # nonsmoking-fibre
sig_feature_list_smoking_LDL <- lapply(daa_results_list_smoking_LDL, function(results) {results %>% filter(p_values < 0.05)
}) # smoking-LDL
sig_feature_list_nonsmoking_LDL <- lapply(daa_results_list_nonsmoking_LDL, function(results) {results %>% filter(p_values < 0.05)
}) # nonsmoking-LDL

view(sig_feature_list_smoking_LDL[[3]])

### Generate pathway heatmap using multiple methods ###
heatmap_list_smoker <- lapply(sig_feature_list, function(sig) {
  pathway_heatmap(abundance = abundance_data %>% filter(pathway %in% sig$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "smoker")
}) # smokers only
# heatmap_list_smoking_fibre <- lapply(sig_feature_list_smoking_fibre, function(sig) {pathway_heatmap(abundance = abundance_data_smoking %>% filter(pathway %in% sig$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category")
# }) # smoking-fibre
# heatmap_list_nonsmoking_fibre <- lapply(sig_feature_list_nonsmoking_fibre, function(sig) {pathway_heatmap(abundance = abundance_data_nonsmoking %>% filter(pathway %in% sig$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category")
# }) # nonsmoking-fibre
heatmap_list_smoking_LDL <- lapply(sig_feature_list_smoking_LDL, function(sig) {pathway_heatmap(abundance = abundance_data_smoking %>% filter(pathway %in% sig$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category")
}) # smoking-LDL
heatmap_list_nonsmoking_LDL <- lapply(sig_feature_list_nonsmoking_LDL, function(sig) {pathway_heatmap(abundance = abundance_data_nonsmoking %>% filter(pathway %in% sig$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category")
}) # nonsmoking-LDL

### Assigning each heatmap ###
## smoker only
heatmap_DESeq2_smoker <- heatmap_list_smoker[[1]]
heatmap_edgeR_smoker <- heatmap_list_smoker[[2]]
heatmap_LinDA_smoker <- heatmap_list_smoker[[3]]
# heatmap_ALDEx2_smoker <- heatmap_list_smoker[[4]]
# ## smoking-fibre
# heatmap_DESeq2_smoking_fibre <- heatmap_list_smoking_fibre[[1]]
# heatmap_edgeR_smoking_fibre <- heatmap_list_smoking_fibre[[2]]
# heatmap_LinDA_smoking_fibre <- heatmap_list_smoking_fibre[[3]]
# # heatmap_ALDEx2_smoking_fibre <- heatmap_list_smoking_fibre[[4]]
# ## nonsmoking-fibre
# heatmap_DESeq2_nonsmoking_fibre <- heatmap_list_nonsmoking_fibre[[1]]
# heatmap_edgeR_nonsmoking_fibre <- heatmap_list_nonsmoking_fibre[[2]]
# heatmap_LinDA_nonsmoking_fibre <- heatmap_list_nonsmoking_fibre[[3]]
# # heatmap_ALDEx2_nonsmoking_fibre <- heatmap_list_nonsmoking_fibre[[4]]
## smoking-LDL
heatmap_DESeq2_smoking_LDL <- heatmap_list_smoking_LDL[[1]]
heatmap_edgeR_smoking_LDL <- heatmap_list_smoking_LDL[[2]]
heatmap_LinDA_smoking_LDL <- heatmap_list_smoking_LDL[[3]]
# heatmap_ALDEx2_smoking_LDL <- heatmap_list_smoking_LDL[[4]]
## nonsmoking-LDL
heatmap_DESeq2_nonsmoking_LDL <- heatmap_list_nonsmoking_LDL[[1]]
heatmap_edgeR_nonsmoking_LDL <- heatmap_list_nonsmoking_LDL[[2]]
heatmap_LinDA_nonsmoking_LDL <- heatmap_list_nonsmoking_LDL[[3]]
# heatmap_ALDEx2_nonsmoking_LDL <- heatmap_list_nonsmoking_LDL[[4]]

### Interesting Pathways ###
interesting_pathways_list <- lapply(sig_feature_list, function(sig) {
  abundance_data %>% filter(pathway %in% sig$feature)
}) # smoker only
# interesting_pathways_list_smoking_fibre <- lapply(sig_feature_list_smoking_fibre, function(sig) {abundance_data_smoking %>% filter(pathway %in% sig$feature)
# }) # smoking-fibre
# interesting_pathways_list_nonsmoking_fibre <- lapply(sig_feature_list_nonsmoking_fibre, function(sig) {abundance_data_nonsmoking %>% filter(pathway %in% sig$feature)
# }) # nonsmoking-fibre
interesting_pathways_list_smoking_LDL <- lapply(sig_feature_list_smoking_LDL, function(sig) {abundance_data_smoking %>% filter(pathway %in% sig$feature)
}) # smoking-LDL
interesting_pathways_list_nonsmoking_LDL <- lapply(sig_feature_list_nonsmoking_LDL, function(sig) {abundance_data_nonsmoking %>% filter(pathway %in% sig$feature)
}) # nonsmoking-LDL


## isolate pathway category
all_interesting_list <- lapply(interesting_pathways_list, function(interesting) {
  data.frame(feature=c(interesting$pathway))
}) # smoker only
# all_interesting_list_smoking_fibre <- lapply(interesting_pathways_list_smoking_fibre, function(interesting) {data.frame(feature=c(interesting$pathway))
# }) # smoking-fibre
# all_interesting_list_nonsmoking_fibre <- lapply(interesting_pathways_list_nonsmoking_fibre, function(interesting) {data.frame(feature=c(interesting$pathway))
# }) # nonsmoking-fibre
all_interesting_list_smoking_LDL <- lapply(interesting_pathways_list_smoking_LDL, function(interesting) {data.frame(feature=c(interesting$pathway))
}) # smoking-LDL
all_interesting_list_nonsmoking_LDL <- lapply(interesting_pathways_list_nonsmoking_LDL, function(interesting) {data.frame(feature=c(interesting$pathway))
}) # nonsmoking-LDL

## add pathway descriptions
all_interesting_list <- lapply(all_interesting_list, function(all_interesting) {
  pathway_annotation(pathway = "MetaCyc", daa_results_df = all_interesting, ko_to_kegg = FALSE)
}) # smoker only
# all_interesting_list_smoking_fibre <- lapply(all_interesting_list_smoking_fibre, function(all_interesting) {pathway_annotation(pathway = "MetaCyc", daa_results_df = all_interesting, ko_to_kegg = FALSE)
# }) # smoking-fibre
# all_interesting_list_nonsmoking_fibre <- lapply(all_interesting_list_nonsmoking_fibre, function(all_interesting) {pathway_annotation(pathway = "MetaCyc", daa_results_df = all_interesting, ko_to_kegg = FALSE)
# }) # nonsmoking-fibre
all_interesting_list_smoking_LDL <- lapply(all_interesting_list_smoking_LDL, function(all_interesting) {pathway_annotation(pathway = "MetaCyc", daa_results_df = all_interesting, ko_to_kegg = FALSE)
}) # smoking-LDL
all_interesting_list_nonsmoking_LDL <- lapply(all_interesting_list_nonsmoking_LDL, function(all_interesting) {pathway_annotation(pathway = "MetaCyc", daa_results_df = all_interesting, ko_to_kegg = FALSE)
}) # nonsmoking-LDL

## smoker only
all_interesting_DESeq2_smoker <- all_interesting_list[[1]]
all_interesting_edgeR_smoker <- all_interesting_list[[2]]
all_interesting_LinDA_smoker <- all_interesting_list[[3]]
# all_interesting_ALDEx2_smoker <- all_interesting_list[[4]]

# ## smoking-fibre
# all_interesting_DESeq2_smoking_fibre <- all_interesting_list_smoking_fibre[[1]]
# all_interesting_edgeR_smoking_fibre <- all_interesting_list_smoking_fibre[[2]]
# all_interesting_LinDA_smoking_fibre <- all_interesting_list_smoking_fibre[[3]]
# # all_interesting_ALDEx2_smoking_fibre <- all_interesting_list_smoking_fibre[[4]]
# 
# ## nonsmoking-fibre
# all_interesting_DESeq2_nonsmoking_fibre <- all_interesting_list_nonsmoking_fibre[[1]]
# all_interesting_edgeR_nonsmoking_fibre <- all_interesting_list_nonsmoking_fibre[[2]]
# all_interesting_LinDA_nonsmoking_fibre <- all_interesting_list_nonsmoking_fibre[[3]]
# # all_interesting_ALDEx2_nonsmoking_fibre <- all_interesting_list_nonsmoking_fibre[[4]]

## smoking-LDL
all_interesting_DESeq2_smoking_LDL <- all_interesting_list_smoking_LDL[[1]]
all_interesting_edgeR_smoking_LDL <- all_interesting_list_smoking_LDL[[2]]
all_interesting_LinDA_smoking_LDL <- all_interesting_list_smoking_LDL[[3]]
# all_interesting_ALDEx2_smoking_LDL <- all_interesting_list_smoking_LDL[[4]]

## nonsmoking-LDL
all_interesting_DESeq2_nonsmoking_LDL <- all_interesting_list_nonsmoking_LDL[[1]]
all_interesting_edgeR_nonsmoking_LDL <- all_interesting_list_nonsmoking_LDL[[2]]
all_interesting_LinDA_nonsmoking_LDL <- all_interesting_list_nonsmoking_LDL[[3]]
# all_interesting_ALDEx2_nonsmoking_LDL <- all_interesting_list_nonsmoking_LDL[[4]]