# Microbial assemblages are structured by host taxonomy and host morphology in a diverse seaweed community 

This page describes the statistical analysis performed in the manuscript entitled "Microbial assemblages are structured by host taxonomy and host morphology in a diverse seaweed community" by Lemay et al. 

### Summary

This page focuses on the analysis of the factor shaping microbiota composition in a diverse seaweed community. Composition dissimilarities are described with beta-diversity metrics and effect of enironmental or host factors are tested using a PerMANOVA approach using the vegan R package with the adonis2 function. When multiple factors are included, we used a type III sum of squares (param "by" set to "margin" in adonis2) in order to test the effect of each factor while controlling for the others. Here, we conduct the main analysis using the Bray-Curtis indice (on rarefied data) on raw OTUs but tested the sensity of our results to various microbial phylogenetic scale using the [BDTT approach](https://github.com/FloMazel/BDTT).

R Script are provided ("Scripts" folder) along with figures ("Outputs" folder)

Data will be provided upon acceptance of the manuscript

### Overall flow of the analysis of microbiota composition

We first produce the OTU table by rarefaction (1000 reads, step #1), then compute beta-diversities (step #2) and then relate these beta-diversities to environmental and host factors (step #3-6). Step #3 test if seaweeds have a microbiome that is distinct from the environment (rock and water). They do, so we explore in more detail the seaweed microbiome. In step #4, we test if the seaweed microbiome is related to host species identities and if it is also related to the tidal zone (low, mid and high) for those species that occur in multiple zones. Both factor are significants, but host species is more important (for 16S data). To vizualize the data, we group samples by host species X tidal zone combinations and produce a dendrogram. To further test the effect of species characteristics (taxonomy and morphology), we group samples per species (mixing different heights) to test for host taxonomic effects ("phylosymbiosis", step #5) and host morphology effects (step #6).


### Summary of the different steps

#### 1. Rarefaction of the OTU tables

Corresponding R script is "Rarefaction_Step.R" (in the beta-diversity computation folder)


#### 2. Beta-diversity computation 

Corresponding R script is "Beta_diversity_Computations.R" (in the beta-diversity computation folder)


#### 3. Effect of substrate on microbiota composition (rock vs. seawater vs. seaweed substrates)

Corresponding R script is "Substrate_Zone_effects.R". Because the sampling of the different substrate is higly imbalanced, and is is known that PerMANOVA could be sensible to this, we also performed PerMANOVA tests on a balanced sampling (results are similars to the full sampling). Corresponding plots and tables are in the "Output" folder


#### 4. Effect of tidal height vs. host species identity

Corresponding R script is "Host_Species_vs_TidalHeight_effects.R". Corresponding plots and tables are in the "Output" folder.

#### 5. Effect of host taxonomy 

Corresponding R script is "Host_taxonomy_effects.R". Corresponding plots and tables are in the "Output" folder.


#### 6. Effect of host morphology

Corresponding R script is "Host_taxonomy_effects.R" Corresponding plots and tables are in the "Output" folder.
