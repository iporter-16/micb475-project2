#Simplified walkthrough for questions 

library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

#read abundance data
abundance_file <- "picrust_processing/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
abundance_data  =as.data.frame(abundance_data)
#change pathways from a column to rownames, required for pathway_daa() function
rownames(abundance_data) = abundance_data$pathway
abundance_data = abundance_data[,-1]

#read metadata
metadata <- read_delim("colombia/metadata_categorized_CL.txt",delim = "\t",escape_double = FALSE,trim_ws = TRUE)

#isolate smoking and nonsmoking from abundance data and metadata
metadata_smoking <- filter(metadata, metadata$smoker=="Yes")
metadata_nonsmoking <- filter(metadata, metadata$smoker=="No")
abundance_data_smoking <- select(abundance_data, metadata_smoking$`#SampleID`)
abundance_data_nonsmoking <- select(abundance_data, metadata_nonsmoking$`#SampleID`)
# abundance_data_smoking <- cbind(abundance_data$pathway,select(abundance_data, metadata_smoking$`#SampleID`))
# abundance_data_nonsmoking <- cbind(abundance_data$pathway,select(abundance_data, metadata_nonsmoking$`#SampleID`))
# names(abundance_data_smoking)[names(abundance_data_smoking) == 'abundance_data$pathway'] <- 'pathway'
# names(abundance_data_nonsmoking)[names(abundance_data_nonsmoking) == 'abundance_data$pathway'] <- 'pathway'

# Perform pathway differential abundance analysis (DAA) using LinDA
# This requires 
daa_results_df <- pathway_daa(abundance = abundance_data, metadata = metadata, 
                              group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_smoking <- pathway_daa(abundance = abundance_data_smoking, metadata = metadata_smoking, 
                                      group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_nonsmoking <- pathway_daa(abundance = abundance_data_nonsmoking, metadata = metadata_nonsmoking, 
                                         group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL)
# Generate pathway heatmap
### ISSUE: there is no column called pathway, because we removed it at line 15. 
### Removing it appears to subsequently make pathway_daa() fail, but then it's needed here?
feature_with_p_0.1 <- daa_results_df_nonsmoking %>% filter(p_adjust < 0.1)
pathway_heatmap(abundance = abundance_data_smoking %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category")
