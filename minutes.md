# Meeting minutes
Organised with date, topic, notetaker and action items.

**Shared server:** root@10.19.139.116 ; Biome1376

## Oct 2
Overall project goal: perform weighted unifrac analysis on the Columbian dataset, looking at metadata categories previously [shown](https://ojs.library.ubc.ca/index.php/UJEMI/article/view/198186/192791?fbclid=IwAR0iTZopMvDnj4u4ff_Y713ByjeSGnvi86pGAkuLxliEXvQDzXDXm4_k-OA) to be significant in ß diversity analysis. 

Hope to identify the confounding variable demonstrated in a prior student [paper](https://ojs.library.ubc.ca/index.php/UJEMI/article/view/198169?fbclid=IwAR2YNekj_Vk7d4nClvOmYG3rJulVkaPYJ1Y8h8GdXCUovs8sDr9Kn98OMaA). Will begin by looking at city, fibre intake, sex, medication. Perform unweighted unifrac analysis as the students did for the original paper, and assess clustering based on these variables (previously noted to have a significant impact on ß diversity metrics) to see if one is responsible for the pattern observed.

Question: what is the confounding variable that causes the population clustering seen above?

**Action items:** start QIIME processing (Sam, with group helping!) + begin project proposal

## Oct 5 (Agenda)
Our options: look for confounding variables (above), or use piecrust for metabolic analysis (fatty acid metabolism)

Things to decide upon: timeline/process, assigning roles and draft deadlines. Using truncation of 225 as this is what the cardiometabolic paper used, and we are investigating their confounding variable(?)

## Oct 5 (Meeting notes)
*Designing question* – confounding variable is a good first aim, but what is the question itself? Need to have background in proposal as to why you are studying this field, drawing from literature.

- Eg. if city, perform lit review into consumption and environment. How is the microbiome of Columbians undergoing diet westernisation affected by geographical location?

*How to find confounding variable* – if PCA plot shows the separation by those two groups!

*Using piecrust* – introduces dataset with metabolic function, to draw connections between microbiome and pathways. Instead of calling a read based on species, calls based on metabolism. Can perform analysis grouping based on function / pathway.

- Would use the same reads, manifest, sequence etc.
- Produces feature tables, one for enzymes and one for pathways.

