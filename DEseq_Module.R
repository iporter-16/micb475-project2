##### Nov 14, 2023 - AW ####

#DESeq module - relative abundance of LDL genus 

#load libraries 
library(tidyverse)
library(phyloseq)
library(DESeq2)

#load phyloseq data 
load("phyloseq_object_final.RData") #do phyloseq final 

#add one to reads 
phyloseq_object_final_plus1 <- transform_sample_counts(phyloseq_object_final, function(x) x+1)
#sep into smoker categories like in picrust

#isolate smoking and nonsmoking from abundance data and metadata
phyloseq_plus1_smoking <- filter(phyloseq_object_final_plus1, ~smoker=="Yes")
phyloseq_plus1_nonsmoking <- filter(metadata, metadata$smoker=="No")

abundance_data_smoking <- select(abundance_data, metadata_smoking$`#SampleID`)
abundance_data_smoking$pathway = abundance_data$pathway
abundance_data_smoking <- abundance_data_smoking %>% select(pathway, everything())

abundance_data_nonsmoking <- select(abundance_data, metadata_nonsmoking$`#SampleID`)
abundance_data_nonsmoking$pathway = abundance_data$pathway
abundance_data_nonsmoking <- abundance_data_nonsmoking %>% select(pathway, everything())

#make DESeq object
phyloseq_object_final_deseq <- phyloseq_to_deseq2(phyloseq_object_final_plus1, ~LDL_category) #what category? 
phyloseq_final <- DESeq(phyloseq_object_final_deseq)

# Make sure that the Healthy group is your reference group
res <- results(phyloseq_final, tidy=TRUE)
View(res)

## Volcano plot: effect size VS significance
ggplot(res) + #show number genes increasing/decreasing abundance compared to no group
  geom_point(aes(x=log2FoldChange, y=-log10(padj))) #mising values due to NAs 

## Make variable to color by whether it is significant + large change
vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>% #new column in results table 
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot

ggsave(filename="vol_plot.png",vol_plot)

# Create bar plot
# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get a vector of ASV names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
phyloseq_filt <- prune_taxa(sigASVs_vec,phyloseq_object_final)
# Add taxonomy onto DESeq results table
merged_results <- tax_table(phyloseq_filt) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

# Make DESeq plot
first_DESeq <- (ggplot(merged_results) +
                  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
                  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
                  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)))
first_DESeq
# Make sure that you have a line that saves the bar plot as a png and this file is present within your project folder

ggsave("phyloseq_DeSeq.png", first_DESeq)