#!/usr/bin/env Rscript
# pulling log2foldchange information

#Installing all packages required
pkgs <- c("phyloseq", "ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "BiocManager", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

devtools::install_github('cafferychen777/ggpicrust2')

#Loading packages
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)

#loading a custom function to provide all DESeq2 results
source("picrust_processing/DESeq2_function.R")

#Importing the pathway PICrust2
abundance_file <- "picrust_processing/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
abundance_data  =as.data.frame(abundance_data)

#Read in your metadata file
metadata <- read_delim(
  "colombia/metadata_categorized_CL.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

#Removing all columns that include NA's for the pathways_daa() function
metadata = metadata[ , colSums(is.na(metadata)) == 0]
#Filter your metadata as needed to look at specific comparisons
metadata_smoking = metadata %>%filter(smoker == "Yes")
metadata_nonsmoking = metadata %>%filter(smoker == "No")

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names_smoking = metadata_smoking$`#SampleID`
sample_names_smoking = append(sample_names_smoking, "pathway")
abundance_data_smoking = abundance_data[, colnames(abundance_data) %in% sample_names_smoking] #This step is the actual filtering

sample_names_nonsmoking = metadata_nonsmoking$`#SampleID`
sample_names_nonsmoking = append(sample_names_nonsmoking, "pathway")
abundance_data_nonsmoking = abundance_data[, colnames(abundance_data) %in% sample_names_nonsmoking] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_smoking =  abundance_data_smoking[, colSums(abundance_data_smoking != 0) > 0]
abundance_data_nonsmoking =  abundance_data_nonsmoking[, colSums(abundance_data_nonsmoking != 0) > 0]

#Filtering out columns that represent a total abundance < 10000
#Filtering out rows (pathways) that have a total count < 100
abundance_data_smoking = abundance_data_smoking[,colSums(abundance_data_smoking[,-1]) > 10000]
abundance_data_smoking = abundance_data_smoking[rowSums(abundance_data_smoking[,-1]) > 100,]

abundance_data_nonsmoking = abundance_data_nonsmoking[,colSums(abundance_data_nonsmoking[,-1]) > 10000]
abundance_data_nonsmoking = abundance_data_nonsmoking[rowSums(abundance_data_nonsmoking[,-1]) > 100,]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_smoking) = NULL
rownames(abundance_data_nonsmoking) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples_smoking = rownames(t(abundance_data_smoking[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadata_smoking = metadata_smoking[metadata_smoking$`#SampleID` %in% abun_samples_smoking,] #making sure the filtered metadata only includes these samples

abun_samples_nonsmoking = rownames(t(abundance_data_nonsmoking[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadata_nonsmoking = metadata_nonsmoking[metadata_nonsmoking$`#SampleID` %in% abun_samples_nonsmoking,] #making sure the filtered metadata only includes these samples

# Workflow for MetaCyc Pathway and Enzyme classification (EC)

# Perform pathway DAA using DESeq2 method
metacyc_daa_results_df_smoking <- pathway_daa(abundance = abundance_data_smoking %>% column_to_rownames("pathway"), metadata = metadata_smoking, group = "LDL_category", daa_method = "DESeq2")
# Annotate MetaCyc pathway results without KO to KEGG conversion. This line provides a useful description for your pathway accessions. 
metacyc_daa_annotated_results_df_smoking <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df_smoking, ko_to_kegg = FALSE)

metacyc_daa_results_df_nonsmoking <- pathway_daa(abundance = abundance_data_nonsmoking %>% column_to_rownames("pathway"), metadata = metadata_nonsmoking, group = "LDL_category", daa_method = "DESeq2")
metacyc_daa_annotated_results_df_nonsmoking <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df_nonsmoking, ko_to_kegg = FALSE)

# # Generate pathway heatmap
# # Please change column_to_rownames() to the feature column if you are not using example dataset
# # Please change group to "your_group_column" if you are not using example dataset
# feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05)
# pathway_heatmap(abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = healthy_only_metadata, group = "Species")
# 
# # Generate pathway PCA plot
# # Please change column_to_rownames() to the feature column if you are not using example dataset
# # Please change group to "your_group_column" if you are not using example dataset
# pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = healthy_only_metadata, group = "Species")

#Generating the full DESeq2 results dataframe
#Here is where we can see HOW up or downregulated a certain pathway is.
#The reference group is selected alphabetically. In this case, Dogs are the reference and anything with a positive log2foldchange means it is higher in humans.

#This function takes two variables.
# 1-the filtered abundance table with pathways in the first column and no rownames
# 2-the filtered metadata
#This dataframe shows the same results as the "metacyc_daa_results_df" above, but it also includes the log2foldchange information that is most valuable to us.

#Here is an example
smoking_res = DEseq2_function(abundance_data_smoking, metadata_smoking)
View(smoking_res)








