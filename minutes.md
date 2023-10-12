# Meeting minutes
Organised by date, descending.
**Shared server:** root@10.19.139.116 ; Biome1376

**Proposal** https://docs.google.com/document/d/1oZlL-Ai-Lp_nA-4W0G4WqwVJkbWYVR4hhuk1EpJy-nc/edit

## Oct 12

### Agenda

| Agenda item                                                                        	| Outcome 	|
|------------------------------------------------------------------------------------	|---------	|
| Confirm project objectives                                                         	|         	|
| Assign proposal roles                                                              	|         	|
| Discuss categorisation cutoffs for fibre/cholesterol / other variables of interest 	|         	|
| Determine ^ associated data manipulation (i.e. categorisation) 	|         	|
| Assign metadata filtering / wrangling role(s)                                      	|         	|

### Meeting notes

### Action items

- Submit proposal to Chris by 17th for feedback
- Have QIIME done by next meeting

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

