#### Oct 24, 2023 - SKA ####
#!/usr/bin/env Rscript

### Beta Diversity ###

# Load packages
library(vegan)
library(phyloseq)
library(tidyverse)

# Load phyloseq object
load("phyloseq_object_final.RData")

# PCoA Plot
# bray-curtis
bc_dm <- distance(phyloseq_object_final, method="bray")
pcoa_bc <- ordinate(phyloseq_object_final, method="PCoA", distance=bc_dm)
plot_ordination(phyloseq_object_final, pcoa_bc, color = "fibre_category", shape="smoker")
plot_ordination(phyloseq_object_final, pcoa_bc, color = "LDL_category", shape="smoker")

gg_pcoa_fibre <- plot_ordination(phyloseq_object_final, pcoa_bc, color = "fibre_category", shape="smoker") +
  labs(pch="smoking status", col = "fibre level")
gg_pcoa_fibre
gg_pcoa_LDL <- plot_ordination(phyloseq_object_final, pcoa_bc, color = "LDL_category", shape="smoker") +
  labs(pch="smoking status", col = "LDL level")
gg_pcoa_LDL

# weighted unifrac
wu_dm <- distance(phyloseq_object_final, method="unifrac", weighted=TRUE)
pcoa_wu <- ordinate(phyloseq_object_final, method="PCoA", wu_dm)
plot_ordination(phyloseq_object_final, pcoa_wu, color = "fibre_category", shape="smoker")
plot_ordination(phyloseq_object_final, pcoa_wu, color = "LDL_category", shape="smoker")

# unweighted unifrac
unwu_dm <- distance(phyloseq_object_final, method="unifrac", weighted=FALSE)
pcoa_unwu <- ordinate(phyloseq_object_final, method="PCoA", unwu_dm)
plot_ordination(phyloseq_object_final, pcoa_unwu, color = "fibre_category", shape="smoker")
plot_ordination(phyloseq_object_final, pcoa_unwu, color = "LDL_category", shape="smoker")

# Save plots
ggsave("plot_pcoa_fibre.png"
       , gg_pcoa_fibre
       , height=4, width=5)

ggsave("plot_pcoa_LDL.png"
       , gg_pcoa_LDL
       , height=4, width=5)