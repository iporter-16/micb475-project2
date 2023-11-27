#### October 27, 2023, AW #####
#Non-smokers and LDL PICRUst Tests 

#load libraries 
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library()

#load in pathway file & convert pathway 
abundance_file <- "picrust_processing/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
abundance_data  =as.data.frame(abundance_data)
rownames(abundance_data) = abundance_data$pathway
abundance_data = abundance_data[,-1]

#read metadata
metadata <- read_delim(
  "colombia/metadata_categorized_CL.txt",
  delim = "\t",escape_double = FALSE,
  trim_ws = TRUE
  )

samples = colnames()
results_file_input <- ggpicrust2(file = abundance_data,
                                 metadata = metadata,
                                 group = "smoker", 
                                 pathway = "MetaCyc",
                                 daa_method = "LinDA",
                                 ko_to_kegg = FALSE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

#select smoking and nonsmoking from abundance data and metadata
metadata_smoking <- filter(metadata, metadata$smoker=="Yes")
metadata_nonsmoking <- filter(metadata, metadata$smoker=="No")
abundance_data_smoking <- select(abundance_data, metadata_smoking$`#SampleID`)
abundance_data_nonsmoking <- select(abundance_data, metadata_nonsmoking$`#SampleID`)
# abundance_data_smoking <- cbind(abundance_data$pathway,select(abundance_data, metadata_smoking$`#SampleID`))
# abundance_data_nonsmoking <- cbind(abundance_data$pathway,select(abundance_data, metadata_nonsmoking$`#SampleID`))
# names(abundance_data_smoking)[names(abundance_data_smoking) == 'abundance_data$pathway'] <- 'pathway'
# names(abundance_data_nonsmoking)[names(abundance_data_nonsmoking) == 'abundance_data$pathway'] <- 'pathway'

# Perform pathway differential abundance analysis (DAA) using LinDA
daa_results_df <- pathway_daa(abundance = abundance_data, metadata = metadata, 
                              group = "LDL_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_smoking <- pathway_daa(abundance = abundance_data_smoking, metadata = metadata_smoking, 
                                      group = "LDL_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_nonsmoking <- pathway_daa(abundance = abundance_data_nonsmoking, metadata = metadata_nonsmoking, 
                                         group = "LDL_category", daa_method = "LinDA", select = NULL, reference = NULL)

