# ------------------------------#
# ------------------------------#
#    Effect of Host Morphology  #
# ------------------------------#
# ------------------------------#

rm(list=ls())

# Question: Does hosts with similar morphologies have similar microbiome, independently of host taxonomy/phylogeny, i.e. convergences? 

# ----------------------
# 1. Loading data and functions 
# -----------------------------

#Loading R packages
library(vegan)
library(dplyr)
library(ggplot2)
library(ape)
library(BiodiversityR)
library(RVAideMemoire)
library(data.tree)
source("Documents/GitHub/Sea_weed_Microbiome/Scripts/Useful_functions.R")

# Load metadata and OTU_table for different host taxonomic ranks
Bray_B_HostTaxo=readRDS("Documents/GitHub/Sea_weed_Microbiome/Outputs/Intermediate_Rfiles/OTUs_Bacteria_table_byHost_Taxonomic_rank_BDTT.RDS")
Bray_Euk_HostTaxo=readRDS("Documents/GitHub/Sea_weed_Microbiome/Outputs/Intermediate_Rfiles/OTUs_Euk_table_byHost_Taxonomic_rank_BDTT.RDS")
Metadata_Host=readRDS("Documents/GitHub/Sea_weed_Microbiome/Outputs/Intermediate_Rfiles/Metadata_byHost_Taxonomic_rank_BDTT.RDS")

# Define the unit of the analysis 
genera=rownames(Bray_B_HostTaxo[[1]][["Genus"]])
genera=genera[-match("Mastocarpus",genera)] #removing Mastocarpus
Metadata_Genus=Metadata_Host[['0']][["Genus"]][genera,]

#Add the morphological data
Metadata=read.table("~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Metadata/Westbeach_Metadata_16s_REVISED_TAXONOMY.txt",header=T,stringsAsFactors = F)
Metadata=Metadata[!Metadata$Genus%in%c("Rock","rock","water","Mastocarpus","Phylospadix"),]
Metadata$Genus=as.character(Metadata$Genus);Metadata$morphology_revised2018=as.character(Metadata$morphology_revised2018)
SummaryMorpho=data.frame(Metadata %>% group_by(Genus) %>% summarise(Morphology=unique(morphology_revised2018))) # retrieve morphology for genus
Metadata_Genus[as.character(SummaryMorpho$Genus),'Morphology']=factor(as.character(SummaryMorpho$Morphology))
Metadata_Genus$Genus=row.names(Metadata_Genus)

# 2. Plot the taxonomic distribution of morphologies
# -----------------------------------------------

TaxonomyH=Metadata_Genus # creating taxonomic tree
TaxonomyH$pathString <- paste(TaxonomyH$Phylum, 
                              TaxonomyH$Class, 
                              TaxonomyH$Order, 
                              TaxonomyH$Family, 
                              TaxonomyH$Genus, 
                              sep = "/")
TaxoTree=as.phylo(as.Node(TaxonomyH))

colors=c('darkgreen',"brown","red","darkblue","orange") #Create a color palette
names(colors)=as.character(unique(TaxonomyH$Morphology))

pdf("Documents/GitHub/Sea_weed_Microbiome/Outputs/Host_taxonomy_tree_withMorphology.pdf",width = 10,height = 10)
plot.phylo(TaxoTree,label.offset=.1,show.node.label = T,type = "fan",open.angle = 45,cex=.5,tip.color = colors[as.character(TaxonomyH[TaxoTree$tip.label,'Morphology'])])
legend(x=70,y=-10,legend=names(colors),fill=colors,cex = .5)
dev.off()

# ----------------------
# 3. Statistics
# ----------------------

# 3.1. Simple PERMANOVA with all samples/all microbial phylo.scales
# -----------------------------------------------------------------

slices=names(Bray_B_HostTaxo)
Permanova_Euk_Morpho=Permanova_Bacteria_Morpho=Permanova_Bacteria_Morpho_family=Permanova_Bacteria_Morpho_Order=Permanova_Euk_Morpho_family=Permanova_Euk_Morpho_Order=list()
nperms=1000
for (s in slices)
{
  # Sub sample metadata and OTU table to this unit
  Bray_B_Genus=Bray_B_HostTaxo[[s]][["Genus"]][genera,genera]
  Bray_Euk_Genus=Bray_Euk_HostTaxo[[s]][["Genus"]][genera,genera]
  
  # Main tatistical analysis
  Permanova_Bacteria_Morpho_family[[s]]=adonis2(Bray_B_Genus~Family+Morphology,data=Metadata_Genus,by="margin",permutations = nperms)
  Permanova_Bacteria_Morpho_Order[[s]]=adonis2(Bray_B_Genus~Order+Morphology,data=Metadata_Genus,by="margin",permutations = nperms)
  Permanova_Bacteria_Morpho[[s]]=adonis2(Bray_B_Genus~Morphology,data=Metadata_Genus,permutations = nperms)
  
  Permanova_Euk_Morpho_family[[s]]=adonis2(Bray_Euk_Genus~Family+Morphology,data=Metadata_Genus,by="margin",permutations = nperms)
  Permanova_Euk_Morpho_Order[[s]]=adonis2(Bray_Euk_Genus~Order+Morphology,data=Metadata_Genus,by="margin",permutations = nperms)
  Permanova_Euk_Morpho[[s]]=adonis2(Bray_Euk_Genus~Morphology,data=Metadata_Genus,permutations = nperms)
  
  print(s)
}

# Create dataframe with results for plotting in ggplot2
testGroup=c("Family","Order","No_cofactor")
F_values_summary=expand.grid(slices=slices,Controlling_for=testGroup,Data=c("16S","18S"))
F_values_summary$F_value=F_values_summary$P_value=NA

for (s in slices)
{
  
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="Family")&(F_values_summary$Data=="16S"),"F_value"]=Permanova_Bacteria_Morpho_family[[s]]$F[2]
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="Order")&(F_values_summary$Data=="16S"),"F_value"]=Permanova_Bacteria_Morpho_Order[[s]]$F[2]
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="No_cofactor")&(F_values_summary$Data=="16S"),"F_value"]=Permanova_Bacteria_Morpho[[s]]$F[1]
  
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="Family")&(F_values_summary$Data=="16S"),"P_value"]=Permanova_Bacteria_Morpho_family[[s]]$`Pr(>F)`[2]
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="Order")&(F_values_summary$Data=="16S"),"P_value"]=Permanova_Bacteria_Morpho_Order[[s]]$`Pr(>F)`[2]
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="No_cofactor")&(F_values_summary$Data=="16S"),"P_value"]=Permanova_Bacteria_Morpho[[s]]$`Pr(>F)`[1]
  
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="Family")&(F_values_summary$Data=="18S"),"F_value"]=Permanova_Euk_Morpho_family[[s]]$F[2]
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="Order")&(F_values_summary$Data=="18S"),"F_value"]=Permanova_Euk_Morpho_Order[[s]]$F[2]
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="No_cofactor")&(F_values_summary$Data=="18S"),"F_value"]=Permanova_Euk_Morpho[[s]]$F[1]
  
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="Family")&(F_values_summary$Data=="18S"),"P_value"]=Permanova_Euk_Morpho_family[[s]]$`Pr(>F)`[2]
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="Order")&(F_values_summary$Data=="18S"),"P_value"]=Permanova_Euk_Morpho_Order[[s]]$`Pr(>F)`[2]
  F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Controlling_for=="No_cofactor")&(F_values_summary$Data=="18S"),"P_value"]=Permanova_Euk_Morpho[[s]]$`Pr(>F)`[1]
  
}

#   3.2. Pairwise PERMANOVA
# ---------------------------

s='0'
Bray_B_Genus=Bray_B_HostTaxo[[s]][["Genus"]][genera,genera]
Bray_Euk_Genus=Bray_Euk_HostTaxo[[s]][["Genus"]][genera,genera]
pairwisePERMANOVA=pairwise.perm.manova(as.dist(Bray_B_Genus),Metadata_Genus$Morphology)
pairwisePERMANOVA_F=pairwise.perm.manovaF(as.dist(Bray_B_Genus),Metadata_Genus$Morphology)

# 3.3. Simple PERMANOVA without crust seaweeds
# --------------------------------------------

s='0'
Non_Crust_Genera=row.names(Metadata_Genus)[!as.character(Metadata_Genus$Morphology)=="crust"]
Metadata_Genus_S1=Metadata_Genus[Non_Crust_Genera,]
Bray_B_Genus_S1=Bray_B_HostTaxo[[s]][["Genus"]][Non_Crust_Genera,Non_Crust_Genera]
Bray_Euk_Genus_S1=Bray_Euk_HostTaxo[[s]][["Genus"]][Non_Crust_Genera,Non_Crust_Genera]

Permanova_Bacteria_Morpho_Order_S1=adonis2(Bray_B_Genus_S1~Order+Morphology,data=Metadata_Genus_S1,by="margin",permutations = 10000)
Permanova_Euk_Morpho_Order_S1=adonis2(Bray_Euk_Genus_S1~Order+Morphology,data=Metadata_Genus_S1,by="margin",permutations = 10000)

Permanova_Bacteria_Morpho_S1=adonis2(Bray_B_Genus_S1~Morphology,data=Metadata_Genus_S1)
Permanova_Euk_Morpho_S1=adonis2(Bray_Euk_Genus_S1~Morphology,data=Metadata_Genus_S1)



# ----------------------
# 4. Plotting and tables
# ----------------------

# 4.1. Export main results at slice "0" as a table
# -----------------------------------------------
s='0'
toexport=rbind(format_aov_table(Permanova_Bacteria_Morpho_Order[[s]],"16S"),format_aov_table(Permanova_Euk_Morpho_Order[[s]],"18S"))
write.csv(toexport,"Documents/GitHub/Sea_weed_Microbiome/Outputs/Tables/Morphology/Host_Morphology-controlling_for_Host_Order_PerMANOVA.csv",na="")
toexport=rbind(format_aov_table(Permanova_Bacteria_Morpho_family[[s]],"16S"),format_aov_table(Permanova_Euk_Morpho_family[[s]],"18S"))
write.csv(toexport,"Documents/GitHub/Sea_weed_Microbiome/Outputs/Tables/Morphology/Host_Morphology-controlling_for_Host_Family_PerMANOVA.csv",na="")

# 4.2. Export BDTT results as a figure 
# ------------------------------------

BDTT_taxo=ggplot(data = F_values_summary, aes(x = slices, y = F_value,colour=Controlling_for,shape=P_value<.05,group=Controlling_for)) +geom_point(cex=5)+facet_wrap(~Data)

pdf("Documents/GitHub/Sea_weed_Microbiome/Outputs/F_value_plots/Host_Morphology_Fvalue_BrayC.pdf",width = 12,height = 7)
BDTT_taxo+stat_summary(geom="line")+scale_shape_manual(values=c(1, 16 ))+ggtitle("Effect of host morphology")#+ylim=c(0.5,1.8)
dev.off()

# 4.3. Export sensitivity analsyis
# --------------------------------

toexport=rbind(format_aov_table(Permanova_Bacteria_Morpho_Order_S1,"16S"),format_aov_table(Permanova_Euk_Morpho_Order_S1,"18S"))
write.csv(toexport,"Documents/GitHub/Sea_weed_Microbiome/Outputs/Tables/Morphology/Host_Morphology_withoutCrust-controlling_for_Host_Order_PerMANOVA.csv",na="")

toexport=rbind(format_aov_table(Permanova_Bacteria_Morpho_S1,"16S"),format_aov_table(Permanova_Euk_Morpho_S1,"18S"))
write.csv(toexport,"Documents/GitHub/Sea_weed_Microbiome/Outputs/Tables/Morphology/Host_Morphology_withoutCrust-No_Phylo_control.csv",na="")


# 4.4. Export pairwise Permanova as tables 
# ----------------------------------------

toexport=round(as.data.frame(pairwisePERMANOVA$p.value),3)
write.csv(toexport,"Documents/GitHub/Sea_weed_Microbiome/Outputs/Tables/Morphology/Host_Morphology_PairwisePermanova_Tests-No_Phylo_control_Pvalues.csv",na="")

toexport=round(as.data.frame(pairwisePERMANOVA_F$F),3)
write.csv(toexport,"Documents/GitHub/Sea_weed_Microbiome/Outputs/Tables/Morphology/Host_Morphology_PairwisePermanova_Tests-No_Phylo_control_Fvalues.csv",na="")


