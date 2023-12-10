#!/bin/R
#PURPOSE: From picrust2 analysis, produces PCA plots for smokers and nonsmokers
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(ggh4x)

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

### Pulling the interesting unexpressed stuff
# interesting_nonsmok <- abundance_data_nonsmoking[row_sds_nonsmok== 0, ]
# interesting_smok <- abundance_data_smoking[row_sds_smok== 0, ]
# all_interesting <- data.frame(feature=c(interesting_nonsmok$pathway,interesting_smok$pathway))
# all_interesting <- pathway_annotation(pathway = "MetaCyc", daa_results_df = all_interesting, ko_to_kegg = FALSE)

### Create plots
PCA_smoking_all <- pathway_pca(abundance = abundance_data %>% column_to_rownames("pathway"), 
                         metadata = metadata, group = "smoker")+
                        ggtitle("PCA plot comparing smoking in all participants")
PCA_LDL_all <- pathway_pca(abundance = abundance_data %>% column_to_rownames("pathway"), 
                           metadata = metadata, group = "LDL_category")+
                        ggtitle("PCA plot comparing LDL in all participants")
PCA_nonsmok_LDL <- pathway_pca(abundance = filtered_LDL_nonsmok %>% column_to_rownames("pathway"), 
                               metadata = metadata_nonsmoking, group = "LDL_category")+
                                ggtitle("            Non-smokers")
PCA_smok_LDL <- pathway_pca(abundance = filtered_LDL_smok %>% column_to_rownames("pathway"), 
                            metadata = metadata_smoking, group = "LDL_category")+
                                ggtitle("            Smokers")

# additional plots for smoker category only and LDL category only among all participants. SKA - Nov 7, 2023 #
pca_plot_smoker <- pathway_pca(abundance = abundance_data %>% column_to_rownames("pathway"), metadata = metadata, group = "smoker")+
  ggtitle("PCA plot comparing smoking in all participants")
pca_plot_smoker+ 
# pca_plot_smoker
pca_plot_LDL <- pathway_pca(abundance = abundance_data %>% column_to_rownames("pathway"), metadata = metadata, group = "LDL_category")+
  ggtitle("PCA plot comparing LDL in all participants")
# pca_plot_LDL

### Save plots. Only change wd to save plots! Keep commented ut
# setwd("/Users/porte16049/Desktop/MICB475/micb475-project2/pcoa_plots")
# ggsave("pca_all_LDL.png",PCA_total,scale=1)
ggsave("pca_smokers_LDL.png",PCA_smok_LDL,scale=1)
ggsave("pca_nonsmokers_LDL.png",PCA_nonsmok_LDL,scale=1)

# setwd("/Users/saman/Desktop/micb475-project2/pcoa_plots")
ggsave("pca_smoker_all.png", pca_plot_smoker, scale=1)
ggsave("pca_LDL_all.png", pca_plot_LDL, scale=1)