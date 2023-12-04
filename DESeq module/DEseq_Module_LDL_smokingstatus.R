##### Nov 14, 2023 - AW ####
#### Dec 03, 2023 - AW - make low LDL the reference group #### 
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


res_genus_combined <- res_genus_combined[order(res_genus_combined$log2FoldChange_avg),]
sighits = ggplot(data = res_genus_combined, aes(x= log2FoldChange_avg,y=reorder(Genus, -(as.numeric(log2FoldChange_avg))), fill = pvalues))+
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)

ggsave("smoking_LDL_phyloseq_DeSeq.png", sighits)

##########################################################################################



## Volcano plot: effect size VS significance
ggplot(res) + #show number genes increasing/decreasing abundance compared to no group
  geom_point(aes(x=log2FoldChange, y=-log10(padj))) #mising values due to NAs 

## Make variable to color by whether it is significant + large change
vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>% #new column in results table 
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot

ggsave(filename="vol_plot_smokingLDL.png",vol_plot)

###########################################################################

##### Nov 16, 2023 - AW ####
#repeat for nonsmoking
phylo_nonsmoking = subset_samples(phyloseq_object_final, smoker == "No")
taxa_info_nonsmoking = as.data.frame(tax_table(phylo_nonsmoking))
taxa_info_nonsmoking$ASV = rownames(taxa_info)

phyloseq_nonsmoke_plus1 <- transform_sample_counts(phylo_nonsmoking, function(x) x+1)
#phyloseq_nonsmoke_deseq <- phyloseq_to_deseq2(phyloseq_nonsmoke_plus1, ~LDL_category) #what category? 
phyloseq_nonsmoke_deseq <- phyloseq_to_deseq2(
  phyloseq_nonsmoke_plus1, 
  ~ relevel(LDL_category, "low"))
phyloseq_nonsmoke_final <- DESeq(phyloseq_nonsmoke_deseq)

res_nonsmoke <- results(phyloseq_nonsmoke_final, tidy=TRUE)

colnames(res_nonsmoke)[1] = "ASV"
View(res_nonsmoke)

res_nonsmoke_taxa = inner_join(taxa_info,res_nonsmoke, by = "ASV" )

res_nonsmoke_sig = res_nonsmoke_taxa %>%
  filter( padj<0.01 & abs(log2FoldChange)>2)

upregulated_count <- sum(res_nonsmoke_sig$log2FoldChange > 0)
downregulated_count <- sum(res_nonsmoke_sig$log2FoldChange < 0)

cat("Number of upregulated nonsmoking ASVs:", upregulated_count, "\n")
cat("Number of downregulated nonsmoking ASVs:", downregulated_count, "\n")

res_nonsmoke_sig <- res_nonsmoke_sig[order(res_nonsmoke_sig$log2FoldChange),]

ggplot(data = res_nonsmoke_sig, aes(y = reorder(Genus, -(as.numeric(log2FoldChange))), x = log2FoldChange, fill = pvalue))+
  geom_col()

combine_pvalues = function(p){
  return(1-pchisq(-2*sum(log(p),na.rm=T),2*sum(!is.na(p))))
}

res_nonsmoke_genus_combined = res_nonsmoke_sig %>%
  group_by(Genus) %>%
  summarize(log2FoldChange_avg = mean(log2FoldChange), pvalues =  combine_pvalues(pvalue))

res_nonsmoke_genus_combined = na.omit(res_nonsmoke_genus_combined)

res_nonsmoke_genus_combined <- res_nonsmoke_genus_combined[order(res_nonsmoke_genus_combined$log2FoldChange_avg),]
sighits = ggplot(data = res_nonsmoke_genus_combined, aes(x= log2FoldChange_avg,y=reorder(Genus, -(as.numeric(log2FoldChange_avg))), fill = pvalues))+
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA) + 
  labs(x = "Average Log2 Fold Change", y = "Genus")

ggsave("nonsmoking_LDL_phyloseq_DeSeq.png", sighits)


# Volcano plot
ggplot(res_nonsmoke) + #show number genes increasing/decreasing abundance compared to no group
  geom_point(aes(x=log2FoldChange, y=-log10(padj))) #mising values due to NAs 

## Make variable to color by whether it is significant + large change
vol_plot <- res_nonsmoke %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>% #new column in results table 
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) + 
  labs(x = "Log2 Fold Change", y = "-Log10 (padj)")
vol_plot

ggsave(filename="vol_plot_nonsmokingLDL.png",vol_plot)
