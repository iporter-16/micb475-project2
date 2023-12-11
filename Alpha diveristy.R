library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggplot2)

#### Load in RData ####
load("phyloseq_object_final_highcount.RData")

plot_richness(phyloseq_object_final_highcount) 

plot_richness(phyloseq_object_final_highcount, measures = c("Shannon","Chao1")) 

gg_richness <- plot_richness(phyloseq_object_final_highcount, x = "smoker", measures = c("Shannon","Chao1")) +
  xlab("Smoker") +
  geom_boxplot()
gg_richness

estimate_richness(phyloseq_object_final_highcount)

phylo_dist <- pd(t(otu_table(phyloseq_object_final_highcount)), phy_tree(phyloseq_object_final_highcount),
                 include.root=F)

sample_data(phyloseq_object_final_highcount)$PD <- phylo_dist$PD

plot.pd <- ggplot(sample_data(phyloseq_object_final_highcount), aes(smoker, PD)) + 
  geom_boxplot() +
  xlab("Smoker") +
  ylab("Phylogenetic Diversity")

plot.pd

sample_data(phyloseq_object_final_highcount)[["legend2"]] = paste(sample_data(phyloseq_object_final_highcount)[[32]],sample_data(phyloseq_object_final_highcount)[[38]])
#sample colors
sample_colors <- c("lightblue", "maroon", "olivedrab","orange")



# plot (fiber)
phyloseq_object_final_highcount %>%                                                              #phyloseq object
  plot_richness(
    x = "legend2",                                                    #compare diversity of datatype
    measures = c("Observed", "Shannon", "Chao1", "ACE", "Simpson", "InvSimpson", "Fisher")) +                              #choose diversity measures
  geom_boxplot(aes(fill = 'legend2'), show.legend = FALSE)+             #make violin plot, set fill aes to sampletype
  theme_linedraw()+                                                     #change theme to classic
  xlab(NULL)+                                                           #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                     #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = -90),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                        #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+           #adjust headings
  scale_fill_manual(values = sample_colors)+                            #set fill colors
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5))

meta = sample_data(phyloseq_object_final_highcount)

alphadiv <- estimate_richness(phyloseq_object_final_highcount, measures = c("Observed", "Shannon", "Chao1", "ACE", "Simpson", "InvSimpson", "Fisher")) %>%
  rownames_to_column(var = "legend2") 
colnames(alphadiv)[1] = "uniqueID"
meta_alpha = inner_join(meta, alphadiv, by = "uniqueID")

meta_alpha = do.call(cbind.data.frame, meta_alpha)


res.aov20 <- aov(Shannon ~ legend2, data = meta_alpha)
summary(res.aov20)
TukeyHSD(res.aov20, which = "legend2")

res.aov21 <- aov(Observed ~ legend2, data = meta_alpha)
summary(res.aov21)
TukeyHSD(res.aov21, which = "legend2")

res.aov22 <- aov(Chao1 ~ legend2, data = meta_alpha)
summary(res.aov22)
TukeyHSD(res.aov22, which = "legend2")

res.aov23 <- aov(ACE ~ legend2, data = meta_alpha)
summary(res.aov23)
TukeyHSD(res.aov23, which = "legend2")

res.aov24 <- aov(Simpson ~ legend2, data = meta_alpha)
summary(res.aov24)
TukeyHSD(res.aov24, which = "legend2")

res.aov25 <- aov(Simpson ~ legend2, data = meta_alpha)
summary(res.aov25)
TukeyHSD(res.aov25, which = "legend2")

#head(alphadiv)

#kruskal.test(Observed ~ V41, data = alphadiv)
#kruskal.test(Shannon ~ V41, data = alphadiv)

##### Tiffany
sample_data(phyloseq_object_final_highcount)[["legend"]] = paste(sample_data(phyloseq_object_final_highcount)[[32]],sample_data(phyloseq_object_final_highcount)[[37]])

# plot (LDL)
phyloseq_object_final_highcount %>%                                                              #phyloseq object
  plot_richness(
    x = "legend",                                                    #compare diversity of datatype
    measures = c("Observed", "Shannon", "Chao1", "ACE", "Simpson", "InvSimpson", "Fisher")) +                              #choose diversity measures
  geom_boxplot(aes(fill = legend), show.legend = FALSE)+             #make violin plot, set fill aes to sampletype
  theme_linedraw()+                                                     #change theme to classic
  xlab(NULL)+                                                           #no label on x-axis
  theme(axis.text.y.left = element_text(size = 20),                     #adjust y-axis text
        axis.text.x = element_text(size = 20, hjust = 0.5, vjust = 1, angle = -90),           #adjust x-axis label position
        axis.title.y = element_text(size = 20))+                        #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 25))+           #adjust headings
  scale_fill_manual(values = sample_colors)+                            #set fill colors
  theme(plot.title=element_text(size = 25, face = "bold", hjust = 0.5))

meta = sample_data(phyloseq_object_final_highcount)

alphadiv <- estimate_richness(phyloseq_object_final_highcount, measures = c("Observed", "Shannon", "Chao1", "ACE", "Simpson", "InvSimpson", "Fisher")) %>%
  rownames_to_column(var = "legend") 
colnames(alphadiv)[1] = "uniqueID"
meta_alpha = inner_join(meta, alphadiv, by = "uniqueID")

meta_alpha = do.call(cbind.data.frame, meta_alpha)


res.aov2 <- aov(Shannon ~ legend, data = meta_alpha)
summary(res.aov2)

TukeyHSD(res.aov2, which = "legend")

res.aov3 <- aov(Observed ~ legend, data = meta_alpha)
summary(res.aov3)
TukeyHSD(res.aov3, which = "legend")

res.aov4 <- aov(Chao1 ~ legend, data = meta_alpha)
summary(res.aov4)
TukeyHSD(res.aov4, which = "legend")

res.aov5 <- aov(ACE ~ legend, data = meta_alpha)
summary(res.aov5)
TukeyHSD(res.aov5, which = "legend")

res.aov6 <- aov(Simpson ~ legend, data = meta_alpha)
summary(res.aov6)
TukeyHSD(res.aov6, which = "legend")

res.aov7 <- aov(InvSimpson ~ legend, data = meta_alpha)
summary(res.aov7)
TukeyHSD(res.aov7, which = "legend")

res.aov8 <- aov(Fisher ~ legend, data = meta_alpha)
summary(res.aov8)
TukeyHSD(res.aov8, which = "legend")

# For the manuscript
length(meta$legend2[meta$legend2 == 'Yes high'])
length(meta$legend2[meta$legend2 == 'Yes low'])
length(meta$smoker[meta$smoker == 'Yes'])
length(meta$smoker[meta$smoker == 'No'])
length(meta$legend2[meta$legend2 == 'No high'])
length(meta$legend2[meta$legend2 == 'No low'])

length(meta$legend[meta$legend == 'Yes high'])
length(meta$legend[meta$legend == 'Yes low'])



### haven't used
alphadiv2 <- estimate_richness(phyloseq_object_final_highcount, measures = c("Observed", "Shannon")) %>%
  rownames_to_column(var = "V42") %>%
  left_join(as.data.frame(sample_data(phyloseq_object_final_highcount)), by = "V42")

head(alphadiv)

kruskal.test(Observed ~ V42, data = alphadiv2)
kruskal.test(Shannon ~ V42, data = alphadiv2)
t.test(dat$response_vector ~ dat$predictor_vector)
t.test(response_vector ~ predictor_vector, data=dat)

result <- t.test(data$alpha_diversity ~ data$group)


alpha_diversity_data <- phyloseq_object_final_highcount

# Perform Wilcoxon rank-sum test for Shannon diversity between categories in V41
shannon_test <- by(alpha_diversity_data@otu_table, alpha_diversity_data@sam_data$LDL_category, function(x) wilcox.test(x, alpha_diversity_data@sam_data$LDL_cateogry))

# Extract the p-values
p_values <- sapply(shannon_test, function(test) test$p.value)

