##2023 Nov 15 TW##

library(DESeq2)
library(phyloseq)
library(tidyverse)

load("phyloseq_object_final.RData")

sample_metadata <- sample_data(phyloseq_object_final)

filtered_smoking <- subset_samples(sample_metadata, `smoker` %in% c('yes', 'no'))

mouse_plus1 <- transform_sample_counts(phyloseq_object_final, function(x) x+1)
mouse_deseq <- phyloseq_to_deseq2(mouse_plus1, ~`LDL_category`)
mouse_fecal <- DESeq(mouse_deseq)

# Make sure that the Healthy group is your reference group
res <- results(mouse_fecal, tidy=TRUE)
View(res)

##Volcano plot
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave(filename="vol_plot.png",vol_plot)



# trial
ggplot(res) + #show number genes increasing/decreasing abundance compared to no group
  geom_point(aes(x=log2FoldChange, y=-log10(padj))) #mising values due to NAs 

## Make variable to color by whether it is significant + large change
vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>% #new column in results table 
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) + 
  labs(x = "Log2 Fold Change", y = "-Log10 (padj)")
vol_plot


# Create bar plot
# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get a vector of ASV names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
mouse_filt <- prune_taxa(sigASVs_vec,phyloseq_object_final)
# Add taxonomy onto DESeq results table
merged_results <- tax_table(mouse_filt) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

# Make DESeq plot
ggplot(merged_results) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Make sure that you have a line that saves the bar plot as a png and this file is present within your project folder
ggsave("phyloseq_DeSeq.png", first_DESeq)


### trial
#load libraries 
library(tidyverse)
library(phyloseq)
library(DESeq2)

#load phyloseq data 
load("phyloseq_object_final.RData") #do phyloseq final 

phylo_smoking = subset_samples(phyloseq_object_final, smoker == "No")


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


## non-smoking LDL
#repeat for nonsmoking
phylo_nonsmoking = subset_samples(phyloseq_object_final, smoker == "No")
taxa_info_nonsmoking = as.data.frame(tax_table(phylo_nonsmoking))
taxa_info_nonsmoking$ASV = rownames(taxa_info_nonsmoking)

phyloseq_nonsmoke_plus1 <- transform_sample_counts(phylo_nonsmoking, function(x) x+1)
phyloseq_nonsmoke_deseq <- phyloseq_to_deseq2(phyloseq_nonsmoke_plus1, ~LDL_category) #what category? 
phyloseq_nonsmoke_final <- DESeq(phyloseq_nonsmoke_deseq)

res_nonsmoke <- results(phyloseq_nonsmoke_final, tidy=TRUE)

colnames(res_nonsmoke)[1] = "ASV"
View(res_nonsmoke)

res_nonsmoke_taxa = inner_join(taxa_info_nonsmoking,res_nonsmoke, by = "ASV" )

res_nonsmoke_sig = res_nonsmoke_taxa %>%
  filter( padj<0.01 & abs(log2FoldChange)>2)

res_nonsmoke_sig <- res_nonsmoke_sig[order(res_nonsmoke_sig$log2FoldChange),]

res_nonsmoke_sig$Genus = c("RF39_ASV1","RF39_ASV2")

ggplot(data = res_nonsmoke_sig, aes(y = reorder(Genus, -(as.numeric(log2FoldChange))), x = log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)+
  geom_col()
  
  
  res_genus_combined <- res_genus_combined[order(res_genus_combined$log2FoldChange_avg),]
  sighits = ggplot(data = res_genus_combined, aes(x= log2FoldChange_avg,y=reorder(Genus, -(as.numeric(log2FoldChange_avg))), fill = pvalues))+
    geom_bar(stat = "identity") +
    theme_bw()+
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA)

combine_pvalues = function(p){
  return(1-pchisq(-2*sum(log(p),na.rm=T),2*sum(!is.na(p))))
}

res_nonsmoke_genus_combined = res_nonsmoke_sig %>%
  group_by(Species) %>%
  summarize(log2FoldChange_avg = mean(log2FoldChange), pvalues =  combine_pvalues(pvalue))

res_nonsmoke_genus_combined = na.omit(res_nonsmoke_genus_combined)

res_nonsmoke_genus_combined <- res_nonsmoke_genus_combined[order(res_nonsmoke_genus_combined$log2FoldChange_avg),]
sighits = ggplot(data = res_nonsmoke_genus_combined, aes(x= log2FoldChange_avg,y=reorder(Species, -(as.numeric(log2FoldChange_avg))), fill = pvalues))+
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA)

ggsave("nonsmoking_LDL_phyloseq_DeSeq.png", sighits)


# Volcano plot
ggplot(res_nonsmoke) + #show number genes increasing/decreasing abundance compared to no group
  geom_point(aes(x=log2FoldChange, y=-log10(padj))) #mising values due to NAs 

## Make variable to color by whether it is significant + large change
vol_plot <- res_nonsmoke %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>% #new column in results table 
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
vol_plot

ggsave(filename="vol_plot_nonsmokingLDL.png",vol_plot)

