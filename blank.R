#blank file for testing
#making modifications

#### Oct 19, 2023 - AW ####

#load libraries 
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#load data 
metafp <- "colombia/metadata_categorized_CL.txt"
meta <- read_delim(metafp, delim="\t")
meta

otufp <- "colombia/colombia_export/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "colombia/colombia_export/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "colombia/colombia_export/rooted_tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

### Format OTU table ###

otu_mat <- as.matrix(otu[,-1]) 
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU) 

### Format updated metadata ###
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'#SampleID'
SAMP <- sample_data(samp_df)
class(SAMP)

### Formatting taxonomy ###
tax_mat <- tax %>% select(-Confidence)%>% 
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)

### Create phyloseq object ###
colombia <- phyloseq(OTU, SAMP, TAX, phylotree)

### Looking at phyloseq object ####
otu_table(colombia)
sample_data(colombia)
tax_table(colombia)
phy_tree(colombia)