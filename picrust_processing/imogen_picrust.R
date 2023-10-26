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
#change pathways from a column to rownames:
rownames(abundance_data) = abundance_data$pathway
abundance_data = abundance_data[,-1]

#read metadata
metadata <- read_delim("colombia/metadata_categorized_CL.txt",
                    delim = "\t",escape_double = FALSE,trim_ws = TRUE)

#isolate smoking and nonsmoking from abundance data and metadata
metadata_smoking <- filter(metadata, metadata$smoker=="Yes")
metadata_nonsmoking <- filter(metadata, metadata$smoker=="No")
abundance_data_smoking <- select(abundance_data, metadata_smoking$`#SampleID`)
abundance_data_nonsmoking <- select(abundance_data, metadata_nonsmoking$`#SampleID`)
# abundance_data_smoking <- cbind(abundance_data$pathway,select(abundance_data, metadata_smoking$`#SampleID`))
# abundance_data_nonsmoking <- cbind(abundance_data$pathway,select(abundance_data, metadata_nonsmoking$`#SampleID`))
# names(abundance_data_smoking)[names(abundance_data_smoking) == 'abundance_data$pathway'] <- 'pathway'
# names(abundance_data_nonsmoking)[names(abundance_data_nonsmoking) == 'abundance_data$pathway'] <- 'pathway'

# write.table(metadata_smoking,file="picrust_processing/imogen_picrust_tests/metadata_smoking.tsv",row.names=FALSE,sep="\t")
# write.table(metadata_nonsmoking,file="picrust_processing/imogen_picrust_tests/metadata_nonsmoking.tsv",row.names=FALSE,sep="\t")
# write.table(abundance_data_smoking,file="picrust_processing/imogen_picrust_tests/abundance_smoking.tsv",row.names=FALSE,sep="\t")
# write.table(abundance_data_nonsmoking,file="picrust_processing/imogen_picrust_tests/abundance_nonmoking.tsv",row.names=FALSE,sep="\t")

# Perform pathway differential abundance analysis (DAA) using LinDA
# Please change group to "your_group_column" if you are not using example dataset
rownames(abundance_data) = abundance_data$pathway
abundance_data = abundance_data[,-1]
daa_results_df <- pathway_daa(abundance = abundance_data, metadata = metadata, 
                              group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_smoking <- pathway_daa(abundance = abundance_data_smoking, metadata = metadata_smoking, 
                                      group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL) 
daa_results_df_nonsmoking <- pathway_daa(abundance = abundance_data_nonsmoking, metadata = metadata_nonsmoking, 
                                      group = "fibre_category", daa_method = "LinDA", select = NULL, reference = NULL)
# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_with_p_0.05 <- daa_results_df_nonsmoking %>% filter(p_adjust < 0.1)
pathway_heatmap(abundance = abundance_data_smoking %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = metadata, group = "fibre_category")

# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pathway_pca(abundance = abundance_data_smoking %>% column_to_rownames("pathway"), metadata = metadata_smoking, group = "fibre_category")

