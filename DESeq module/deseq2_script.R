##### Nov 14, 2023 - AW ####
##### Nov 27, 2023 - AW - updated to add labels####

#DESeq module - relative abundance of LDL genus 

#load libraries 
library(tidyverse)
library(phyloseq)
library(DESeq2)

#load phyloseq data 
load("phyloseq_object_final.RData") #do phyloseq final 


#filter to smoking or non smoking 
phylo_smoking = subset_samples(phyloseq_object_final, smoker == "Yes")


################################################################ CHRIS ADDED
taxa_info = as.data.frame(tax_table(phylo_smoking))
taxa_info$ASV = rownames(taxa_info)

################################################################
#add one to reads 
phyloseq_object_final_plus1 <- transform_sample_counts(phylo_smoking, function(x) x+1)
phyloseq_object_final_deseq <- phyloseq_to_deseq2(phyloseq_object_final_plus1, ~LDL_category) #what category? 
phyloseq_final <- DESeq(phyloseq_object_final_deseq)

# Make sure that the Healthy group is your reference group
res <- results(phyloseq_final, tidy=TRUE)


#################################################################################### CHRIS ADDED
#Change first column name from "row" to "ASV" beacuse thats what I did on line 18
colnames(res)[1] = "ASV"
View(res)

#Merge taxa information with "res" 
res_with_taxa = inner_join(taxa_info,res, by = "ASV" )

#take significant hits
res_sig = res_with_taxa %>%
  filter( padj<0.01 & abs(log2FoldChange)>2)

#Order by log2FC
res_sig <- res_sig[order(res_sig$log2FoldChange),]
#Make the plot
ggplot(data = res_sig, aes(y = reorder(Genus, -(as.numeric(log2FoldChange))), x = log2FoldChange, fill = pvalue))+
  geom_col()

#Function to combined Pvalues
combine_pvalues = function(p){
  return(1-pchisq(-2*sum(log(p),na.rm=T),2*sum(!is.na(p))))
}

#Merging the data for rows that have the same genus. Need to take the average for log2FC and combined the Pvalues
res_genus_combined = res_sig %>%
  group_by(Genus) %>%
  summarize(log2FoldChange_avg = mean(log2FoldChange), pvalues =  combine_pvalues(pvalue))
#Remove NAs
res_genus_combined = na.omit(res_genus_combined)

#Remove g__  
res_genus_combined$Genus <- gsub("g__", "", res_genus_combined$Genus)


res_genus_combined <- res_genus_combined[order(res_genus_combined$log2FoldChange_avg),]
sighits = ggplot(data = res_genus_combined, aes(x= log2FoldChange_avg,y=reorder(Genus, -(as.numeric(log2FoldChange_avg))), fill = pvalues))+
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA) + 
  labs(x = "Average Log2 Fold Change", y = "Genus")

ggsave("smoking_LDL_phyloseq_DeSeq.png", sighits)

##########################################################################################


## Volcano plot: effect size VS significance
ggplot(res) + #show number genes increasing/decreasing abundance compared to no group
  geom_point(aes(x=log2FoldChange, y=-log10(padj))) #mising values due to NAs 

## Make variable to color by whether it is significant + large change
vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>% #new column in results table 
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) + 
  labs(x = "Log2 Fold Change", y = "-Log10 (padj)")
vol_plot

ggsave(filename="vol_plot.png",vol_plot)



##### Dec 01, 2023 - AW - add volcano point labels using EnhancedVolcano package ####


if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

# BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)

#Annotate the Ensembl gene IDs to gene symbols

#IP changes: included all data, reduced log2foldchange to 1 in accordance with bar plots, 
#fixed labels so they were the same throughout allowing genera of interest to be labelled
res_with_taxa_noNA <- res_with_taxa %>% filter(!is.na(res_with_taxa$padj))
res_with_taxa_noNA$Genus <- gsub("g__", "", res_with_taxa_noNA$Genus)
labels <- res_genus_combined$Genus

volcano_plot <- EnhancedVolcano(
  res_with_taxa_noNA,
  lab = res_with_taxa_noNA$Genus, selectLab = labels,
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-5, 5), # Adjust xlim based on your data distribution
  ylim = c(-1, 15),
  title = 'EnhancedVolcano Plot',
  pCutoff = 0.05, # Adjust p-value cutoff
  FCcutoff = 2,   # Adjust log2 fold change cutoff
  pointSize = 3   # Adjust point size
)
print(volcano_plot)

#Issue is that many ASVs are within that genus. As a solution, only fix the genus names (remove g__) 
#of the significant ASVs so that they are the only ones labelled. IP 03/12
res_with_taxa_noNA <- res_with_taxa %>% filter(!is.na(res_with_taxa$padj))
res_with_taxa_noNA <-  dplyr::mutate(res_with_taxa_noNA, Genus= if_else(ASV %in% res_sig$ASV, gsub("g__","",Genus),Genus))
res_sig$Genus <- gsub("g__", "", res_sig$Genus)

EnhancedVolcano(
  res_with_taxa_noNA,
  lab = res_with_taxa_noNA$Genus, selectLab = labels,
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-5, 5), # Adjust xlim based on your data distribution
  ylim = c(-1, 15),
  title = paste('EnhancedVolcano Plot - ',labels[l]),
  pCutoff = 0.05, # Adjust p-value cutoff
  FCcutoff = 2,   # Adjust log2 fold change cutoff
  pointSize = 2,   # Adjust point size
  labCol = 'black',
  labFace = 'bold',
  # boxedLabels = TRUE,
  drawConnectors = TRUE,
  widthConnectors = 1.0,
  colConnectors = 'black',
  max.overlaps=50,
  labSize = 4)




















