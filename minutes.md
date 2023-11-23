# Meeting minutes
Organised by date, descending.
**Shared server:** root@10.19.139.116 ; Biome1376

**Proposal** https://docs.google.com/document/d/1oZlL-Ai-Lp_nA-4W0G4WqwVJkbWYVR4hhuk1EpJy-nc/edit

[Oct 2](#oct-2) ; [Oct 5](#oct-5) ; [Oct 12](#oct-12) ; [Oct 19](#oct-19) ; [Oct 26](#oct-26) ; [Nov 2](#nov-2) ; [Nov 9](#nov-9)

## Nov 23

- Fig1 Alpha diversity: present observed and shannon in 2-panel figure (LDL). Relabel to be specific about groups.
- Fig2 Beta diversity: present smoking + nonsmoking in 2-panel figure (LDL)
- Fig3 Picrust barplots and heatmaps in 4-panel figure
- Fig4 Volcano and DESeq bar plot for smokers
- Supp1: fibre diversity plots,
- Supp2: PCA plots from picrust analysis
- Supp3: volcano plot (nonsmokers sig results)

Can't make huge claims/suggestions, only specifics.

Takeaway: smokers are more susceptible and volatile to differences in LDL, shown both in picrust and deseq plots.

Mention that there are two points in the nonsmokers volcano plots – possible NA, but go back and BLAST the sequence to link it to something.

**Presentation:** stick to picrust barplots and significance figures. Present a clear takeaway (smokers need to watch their LDL!)


## Nov 9
| Agenda item                                                                        	| Conclusion 	|
|-----------------------------------------	|---------	|
| Alpha diversity: no significant data is shown between the groups                        |  |
| Distribute tasks for manuscript writing?                        |  |
| Analyses Comparisons: how do we get the top 10 pathways? |  |  

### Next steps:

Tiffany + Alice –> **DESeq analysis** to figure out if there are a few specific genus / species that are affected in abundance, resulting in changes in metabolic information during picrust. Run on taxonomic-level table (phyloseq object w/OTU table + taxonomy file, prior to picrust analyses), then lit review for metabolic functions to see if it aligns with metabolic results

- Test smokers high/low LDL, check which bacteria are upregulated.
- Hope to see less difference between high/low LDL in nonsmokers.

Sean + Sam -> Perform using **Chris's DESeq function (will be emailed)** on smoker/nonsmoker high/low LDL groups.

Start manuscript drafting!

- Imogen –> Introduction
- Tiffany + Sam –> Figure 1 (diversity)
- Sam –> Methods
- Sean + Alice –> Discussion later...


## Nov 2

### To-do:

- **Imogen:** Export all plots onto Github _(done!)_
- Combine plots 
- **Imogen & ?:** [Lit review](https://docs.google.com/document/d/1gAc4GGjX9uzj-aR70NqgcJvRnQcLFoWcT7jLYDcZJVs/edit) on results: pathways up/down regulated in LDL samples
- **Imogen:**  Lit review on results: pathways needed to be excluded in PCA plots _(done! dead end -i)_
- Compare top 10 results between LinDA, deseq, aldex2, edgeR
- **Tiffany:** [Alpha diversity signifiance comparison](https://docs.google.com/presentation/d/1pIplkOJ19TjLork48qhRfuhgrWncWUlC/edit?usp=sharing&ouid=103207874254987114452&rtpof=true&sd=true)

## Oct 26

| Agenda item                                                                        	| Conclusion 	|
|-----------------------------------------	|---------	|
| Picrust tutorial                        |See file  |
| Looking over alpha metrics – LDL (Tiffany)    |LDL appears to have no effect on a-diversity so has no big impact on compositional changes, but don't disregard LDL as there may still be smaller-level functional changes occuring (DESeq2 will tell us if specific species are up/downregulated)|
| Looking over alpha metrics – fiber (Tiffany) | Seems that without smoking fiber intake doesn't affect alpha diversity, but with smoking there is an impact: smoking + high fiber appears to correspond to lower alpha diversity (specifically no difference in evenness, but change in richness)* | 
| Looking over ß diversity | Axes are very small % meaning they show a very small amount of diversity difference between the two. Thus there is not a huge overarching impact of fiber or LDL on ß diveristy|

*For the paper, use Observed, Chao1, Shannon.

### Action items

Picrust analyses: run using pathways table found in **picrust2-2.5.2/picrust2_out_pipeline/pathways_out/**

- Sam –> Pure smoking
- Imogen –> Fiber, nonsmoking
- Sean –> Fiber, smoking
- Alice –> LDL, nonsmoking
- Tiffany –> LDL, smoking


## Oct 19

### Agenda

| Agenda item                                                                        	| Conclusion 	|
|------------------------------------------------------------------------------------	|---------	|
| Draft proposal feedback                                                         	|x |
| Assign data analysis roles based on approach table, below                                                         	|x |

### Action items

- Fix Chris's suggestions in proposal, email to Chris by Friday noon (ALL)
- Run picrust on dataset (Imogen)
- Set up R project on github (Imogen)
- First milestone: make phyloseq object (Alice)
- Rarefy data after phyloseq (Sean) using 22700 -> phyloseq_object_final
- Alpha diversity analyses after phyloseq_final (Tiffany)
- Beta diversity analyses after phyloseq_final (Sam)

## Oct 12

### Agenda

| Agenda item                                                                        	| Conclusion 	|
|------------------------------------------------------------------------------------	|---------	|
| Confirm project objectives                                                         	|[Done](#project-aims) |
| Assign proposal roles                                                              	|[Done](#proposal-allocations)|
| Discuss categorisation cutoffs for fibre/cholesterol / other variables of interest 	| Fibre cutoffs: above/below 21g(f), 30g(m). LDL cutoffs: <100mg(mf). |
| Determine ^ associated data manipulation (i.e. categorisation) 	|         	|
| Assign metadata filtering / wrangling role(s)                                      	|Tiffany|

### Project Aims

1. Identify changes in microbiome composition (a/ß diversity) in smokers vs. nonsmokers, and by diet categories

   1.1. bin data based on fiber and cholesterol levels (H/L)

   *Figure outputs: All a-diversity metrics (Shannons, Faith's PD) in boxplots by non/smoke and diet categories*
   
2. Analyse changes in metabolic function between smokers and non-smokers, using picrust.

   2.1. Install picrust, process data to remove all features in table.qza file that are mitochondria/chloroplasts and also <5 (increase efficiency) through QIIME

   *Figure outputs: Pathway tables, enzyme tables. Produce table of top 10 upregulated pathways, and top 10 downregulated. Produce heatmap & volcano  to visualise.*

3. Within the smoking population, how how does fibre intake affect metabolic gut function?

   Split into four groups (nonsmoker high-fibre, nonsmoker low-fibre, smoker high-fibre, smoker low-fibre), generate plot showing which pathways are up/downregulated. Control is low-fibre group.

   *Figure outputs: same as A2*

4. Within the smoking population, how how does LDL intake affect metabolic gut function?

   Same as A3.


### Proposal allocations

- Intro & Background –> Imogen
- Objectives, Aims & Rationale –> Alice
- Approach & Timeline –> Sean
- Overview flowchart –> Tiffany/Sam
- QIIME processing (cutoff @248) –> Tiffany/Sam

### Action items

- Submit proposal to Chris by Oct 17th for feedback
- Metadata filtering asap –> Tiffany
- Look at qza repseqs table dmux files, rarefaction curve on **Oct 17th** –> Sam/Tiffany
- If possible, a-diversity analysis by **Oct 17t** for smoking/nonsmoking populations.

**Timeline** 
- Picrust analysis by Oct 22nd
- Coding finished by early Nov
- Start manuscript by mid Nov

  
## Oct 5 

### Agenda
Our options: look for confounding variables (above), or use picrust for metabolic analysis (fatty acid metabolism)

Things to decide upon: timeline/process, assigning roles and draft deadlines. Using truncation of 225 as this is what the cardiometabolic paper used, and we are investigating their confounding variable(?)

### Meeting notes
*notetaker: Imogen*

*Designing question* – confounding variable is a good first aim, but what is the question itself? Need to have background in proposal as to why you are studying this field, drawing from literature.

- Eg. if city, perform lit review into consumption and environment. How is the microbiome of Columbians undergoing diet westernisation affected by geographical location?

*How to find confounding variable* – if PCA plot shows the separation by those two groups!

*Using picrust* – introduces dataset with metabolic function, to draw connections between microbiome and pathways. Instead of calling a read based on species, calls based on metabolism. Can perform analysis grouping based on function / pathway.

- Would use the same reads, manifest, sequence etc.
- Produces feature tables, one for enzymes and one for pathways.
- Compare the changes in pathways and function with the changes in microbiome highlighted in the student paper on cardiometabolic status.

**Question: How does the microbial composition reflect metabolic function between smokers and non-smokers? How significantly does diet affect function?**

**Aims:**

1. Changes in microbiome composition (a/ß diversity) in smokers vs. nonsmokers, and by diet categories (see below).
2. Changes in metabolic function between smokers and non-smokers, using picrust. 
3. Within the smoking population, how how does **fibre intake** affect metabolic gut function? Split into four groups (nonsmoker high-fibre, nonsmoker low-fibre, smoker high-fibre, smoker low-fibre), generate plot showing which pathways are up/downregulated. *Control is low-fibre group?*
4. Within the smoking population, how does **cholesterol intake** affect metabolic function? As above.


**Workflow suggestion**

- Use QIIME processing to generate feature tables and phylogenetic tree –> look at bacterial abundance and compare to picrust metabolic data. Eg. Get subset of upregulated pathways, the bacteria associated with these pathways, and whether any of those bacteria are upregulated.
- Within the smoker category, create dietary pattern categories. Which pattern results in the least dysbiosis? 
- Investigate the change in short-chain fatty acid metabolism and the link to firmicute population changes.
- Link above to the gut-brain axis, is there a connection to neurotransmitters?
- Investigate the impact of geographic location on smoking/nonsmoking metabolic statuses

*Note: picrust is making huge assumptions as it does not directly correlate bacteria to function*

### Action items

- Chris will send picrust information –> download software, learn how to run
- Decide who is doing what for the proposal.



## Oct 2
Overall project goal: perform weighted unifrac analysis on the Columbian dataset, looking at metadata categories previously [shown](https://ojs.library.ubc.ca/index.php/UJEMI/article/view/198186/192791?fbclid=IwAR0iTZopMvDnj4u4ff_Y713ByjeSGnvi86pGAkuLxliEXvQDzXDXm4_k-OA) to be significant in ß diversity analysis. 

Hope to identify the confounding variable demonstrated in a prior student [paper](https://ojs.library.ubc.ca/index.php/UJEMI/article/view/198169?fbclid=IwAR2YNekj_Vk7d4nClvOmYG3rJulVkaPYJ1Y8h8GdXCUovs8sDr9Kn98OMaA). Will begin by looking at city, fibre intake, sex, medication. Perform unweighted unifrac analysis as the students did for the original paper, and assess clustering based on these variables (previously noted to have a significant impact on ß diversity metrics) to see if one is responsible for the pattern observed.

Question: what is the confounding variable that causes the population clustering seen above?

**Action items:** start QIIME processing (Sam, with group helping!) + begin project proposal

