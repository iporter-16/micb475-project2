#### Nov 9, 2023 - SKA ####
#!/usr/bin/env Rscript

### Beta Diversity ###

# Load packages
library(vegan)
library(phyloseq)
library(tidyverse)

# Load phyloseq object
load("phyloseq_object_final_highcount_highcount.RData")
samp_dat_wdiv <- data.frame(sample_data(phyloseq_object_final_highcount), estimate_richness(phyloseq_object_final_highcount))

### PERMANOVA (Permutational ANOVA) ####
# Use phyloseq to calculate distance matrix
dm_weighted_unifrac <- UniFrac(phyloseq_object_final_highcount, weighted=TRUE)#weighted unifrac
dm_unweighted_unifrac <- UniFrac(phyloseq_object_final_highcount_highcount, weighted=FALSE)#unweighted unifrac
dm_bray <- vegdist(t(otu_table(phyloseq_object_final_highcount)), method="bray")#bray-curtis
dm_jaccard <- vegdist(t(otu_table(phyloseq_object_final_highcount)), method="jaccard")#jaccard
# permanova: smoker-fibre_category
adonis2(dm_weighted_unifrac ~ `smoker`*fibre_category, data=samp_dat_wdiv)#weighted unifrac
adonis2(dm_unweighted_unifrac ~ `smoker`*fibre_category, data=samp_dat_wdiv)#unweighted unifrac
adonis2(dm_bray ~ `smoker`*fibre_category, data=samp_dat_wdiv)#bray-curtis
adonis2(dm_jaccard ~ `smoker`*fibre_category, data=samp_dat_wdiv)#jaccard
# permanova: smoker-LDL_category
adonis2(dm_weighted_unifrac ~ `smoker`*LDL_category, data=samp_dat_wdiv)#weighted unifrac
adonis2(dm_unweighted_unifrac ~ `smoker`*LDL_category, data=samp_dat_wdiv)#unweighted unifrac
adonis2(dm_bray ~ `smoker`*LDL_category, data=samp_dat_wdiv)#bray-curtis
adonis2(dm_jaccard ~ `smoker`*LDL_category, data=samp_dat_wdiv)#jaccard

### PCoA Plot - Final Figure ###
# Beta Diversity figure panel will include Weighted Unifrac PCoA Plot of the smoker and LDL categories, 
# divided into 2 plots, one with smokers only and the other with nonsmokers only
## Filter smokers only and nonsmokers only ##
view(sample_data(phyloseq_object_final_highcount_highcount))
phyloseq_object_smoking <- subset_samples(phyloseq_object_final_highcount_highcount, smoker == "Yes")
phyloseq_object_nonsmoking <- subset_samples(phyloseq_object_final_highcount_highcount, smoker == "No")
view(sample_data(phyloseq_object_smoking))
## Generate PCoA plots using weighted unifrac ##
ord.weighted_unifrac_smoking <- ordinate(phyloseq_object_smoking, method="PCoA", distance="unifrac", weighted=TRUE)
ord.weighted_unifrac_nonsmoking <- ordinate(phyloseq_object_nonsmoking, method="PCoA", distance="unifrac", weighted=TRUE)
gg_pcoa_wunifrac_smoking_LDL <- plot_ordination(phyloseq_object_smoking, ord.weighted_unifrac_smoking, color = "LDL_category") +
  labs(col="LDL category") + stat_ellipse(type = "norm")#smoking-LDL
gg_pcoa_wunifrac_nonsmoking_LDL <- plot_ordination(phyloseq_object_nonsmoking, ord.weighted_unifrac_nonsmoking, color = "LDL_category") +
  labs(col="LDL category") + stat_ellipse(type = "norm")#nonsmoking-LDL

gg_pcoa_wunifrac_smoking_LDL
gg_pcoa_wunifrac_nonsmoking_LDL
### Save plots ###
# setwd("/Users/saman/Desktop/micb475-project2/pcoa_plots/beta_diversity")
ggsave("plot_pcoa_wunifrac_smoking_LDL.png", gg_pcoa_wunifrac_smoking_LDL, height=2, width=5)
ggsave("plot_pcoa_wunifrac_nonsmoking_LDL.png", gg_pcoa_wunifrac_nonsmoking_LDL, height=2, width=5)
# setwd("/Users/saman/Desktop/micb475-project2")

#### Oct 26, 2023 - SKA ####
#Making a new column in the metadata to combined fibre and smoking
sample_data(phyloseq_object_final_highcount)[[40]] = paste(sample_data(phyloseq_object_final_highcount)[[32]],sample_data(phyloseq_object_final_highcount)[[38]])
#Making a new column in the metadata to combined LDL and smoking
sample_data(phyloseq_object_final_highcount)[[41]] = paste(sample_data(phyloseq_object_final_highcount)[[32]],sample_data(phyloseq_object_final_highcount)[[37]])

### PCoA plot ###
# weighted unifrac
ord.weighted_unifrac <- ordinate(phyloseq_object_final_highcount, method="PCoA", distance="unifrac", weighted=TRUE)
gg_pcoa_wunifrac_fibre <- plot_ordination(phyloseq_object_final_highcount, ord.weighted_unifrac, color = "V40") +
  labs(col="smoker, fibre category") + stat_ellipse(type = "norm")#smoker-fibre
gg_pcoa_wunifrac_LDL <- plot_ordination(phyloseq_object_final_highcount, ord.weighted_unifrac, color = "V41") +
  labs(col="smoker, LDL category") + stat_ellipse(type = "norm")#smoker-LDL
# unweighted unifrac
ord.unweighted_unifrac <- ordinate(phyloseq_object_final_highcount, method="PCoA", distance="unifrac", weighted=FALSE)
gg_pcoa_unwunifrac_fibre <- plot_ordination(phyloseq_object_final_highcount, ord.unweighted_unifrac, color = "V40") +
  labs(col="smoker, fibre category") + stat_ellipse(type = "norm")#smoker-fibre
gg_pcoa_unwunifrac_LDL <- plot_ordination(phyloseq_object_final_highcount, ord.unweighted_unifrac, color = "V41") +
  labs(col="smoker, LDL category") + stat_ellipse(type = "norm")#smoker-LDL
# bray-curtis
ord.bray <- ordinate(phyloseq_object_final_highcount, method="PCoA", distance="bray")
gg_pcoa_bray_fibre <- plot_ordination(phyloseq_object_final_highcount, ord.bray, color = "V40") +
  labs(col="smoker, fibre category") + stat_ellipse(type = "norm")#smoker-fibre
gg_pcoa_bray_LDL <- plot_ordination(phyloseq_object_final_highcount, ord.bray, color = "V41") +
  labs(col="smoker, LDL category") + stat_ellipse(type = "norm")#smoker-LDL

### Save plots ###
# setwd("/Users/saman/Desktop/micb475-project2/pcoa_plots/beta_diversity")
# ggsave("plot_pcoa_wunifrac_fibre.png", gg_pcoa_wunifrac_fibre, height=4, width=5)
# ggsave("plot_pcoa_wunifrac_LDL.png", gg_pcoa_wunifrac_LDL, height=4, width=5)
# ggsave("plot_pcoa_unwunifrac_fibre.png", gg_pcoa_unwunifrac_fibre, height=4, width=5)
# ggsave("plot_pcoa_unwunifrac_LDL.png", gg_pcoa_unwunifrac_LDL, height=4, width=5)
# ggsave("plot_pcoa_bray_fibre.png", gg_pcoa_bray_fibre, height=4, width=5)
# ggsave("plot_pcoa_bray_LDL.png", gg_pcoa_bray_LDL, height=4, width=5)
# setwd("/Users/saman/Desktop/micb475-project2")

#### Oct 24, 2023 - SKA ####
# # PCoA Plot
# # bray-curtis
# bc_dm <- distance(phyloseq_object_final_highcount, method="bray")
# pcoa_bc <- ordinate(phyloseq_object_final_highcount, method="PCoA", distance=bc_dm)
# 
# gg_pcoa_bc_fibre <- plot_ordination(phyloseq_object_final_highcount, pcoa_bc, color = "fibre_category", shape="smoker") +
#   labs(pch="smoking status", col = "fibre level") + stat_ellipse()
# gg_pcoa_bc_fibre
# gg_pcoa_bc_LDL <- plot_ordination(phyloseq_object_final_highcount, pcoa_bc, color = "LDL_category", shape="smoker") +
#   labs(pch="smoking status", col = "LDL level") + stat_ellipse()
# gg_pcoa_bc_LDL
# 
# # weighted unifrac
# wu_dm <- distance(phyloseq_object_final_highcount, method="unifrac", weighted=TRUE)
# pcoa_wu <- ordinate(phyloseq_object_final_highcount, method="PCoA", wu_dm)
# gg_pcoa_wu_fibre <- plot_ordination(phyloseq_object_final_highcount, pcoa_wu, color = "fibre_category", shape="smoker") +
#   labs(pch="smoking status", col = "fibre level") + stat_ellipse()
# gg_pcoa_wu_fibre
# gg_pcoa_wu_LDL <- plot_ordination(phyloseq_object_final_highcount, pcoa_wu, color = "LDL_category", shape="smoker") +
#   labs(pch="smoking status", col = "fibre level") + stat_ellipse()
# gg_pcoa_wu_LDL
# 
# # unweighted unifrac
# unwu_dm <- distance(phyloseq_object_final_highcount, method="unifrac", weighted=FALSE)
# pcoa_unwu <- ordinate(phyloseq_object_final_highcount, method="PCoA", unwu_dm)
# gg_pcoa_unwu_fibre <- plot_ordination(phyloseq_object_final_highcount, pcoa_unwu, color = "fibre_category", shape="smoker") +
#   labs(pch="smoking status", col = "fibre level") + stat_ellipse()
# gg_pcoa_unwu_fibre
# gg_pcoa_unwu_LDL <- plot_ordination(phyloseq_object_final_highcount, pcoa_unwu, color = "LDL_category", shape="smoker") +
#   labs(pch="smoking status", col = "fibre level") + stat_ellipse()
# gg_pcoa_unwu_LDL
# 
# # Save plots
# ggsave("plot_pcoa_fibre.png"
#        , gg_pcoa_fibre
#        , height=4, width=5)
# 
# ggsave("plot_pcoa_LDL.png"
#        , gg_pcoa_LDL
#        , height=4, width=5)