##-----------------------------##
## Beta-diversity computations ##
##-----------------------------##

rm(list=ls())

source("~/Documents/GitHub/BDTT/BDTT_functions.R")
library(GUniFrac)
library(betapart)
library(vegan)
library(biomformat)
library(ape)
library(phytools)
library(abind)
library(Matrix)


EukTree=read.tree("~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Phylogenies/RAxML_labelledTree.EPA_placement_rooted.figtree.tre")
EukTree

BactTree=read.tree("~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Phylogenies/Westbeach_16s_MEDS_tree.tre")
BactTree=reroot(tree=BactTree, node.number=sample(x=c(BactTree$edge),size=1))
BactTree

# Get Betas
#----------

slices=c(seq(from=0,to=.7,by=0.05)) #BDTT


#--------------#
#  16 S data   #
#--------------#

Bacteria_OTU_t_raref=read.table(file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Rarefied_Full_OTUs_Tables/WB_16s_OTUs_NoSeagrass_raref1000.txt")

commonBact=(intersect(row.names(Bacteria_OTU_t_raref),BactTree$tip.label))
BactTree=drop.tip(BactTree,BactTree$tip.label[!BactTree$tip.label%in%commonBact])
Hnodes=getHnodes(BactTree) #get Nodes height 
unifracs <- aperm(GUniFrac(t(Bacteria_OTU_t_raref), BactTree, alpha=c(0, 0.5, 1))$unifracs,c(3,1,2)) #Unifrac
BDTT_results=BDTT(similarity_slices=slices,tree=BactTree,sampleOTUs=Bacteria_OTU_t_raref,onlyBeta=T,metric=c("jac","bc"))
  
saveRDS(BDTT_results,file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Beta_diversity/BDTT_Rarefs1000_16S.RDS")
saveRDS(unifracs,file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Beta_diversity/unifracs_Rarefs1000_16S.RDS")
  

#--------------#
#   18S data   #
#--------------#

Euk_OTU_t_raref=read.table(file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Rarefied_Full_OTUs_Tables/WB_18_OTUs_NoSeagrass_raref1000.txt")

commonEuk=(intersect(row.names(Euk_OTU_t_raref),EukTree$tip.label))
EukTree=drop.tip(EukTree,EukTree$tip.label[!EukTree$tip.label%in%commonEuk])
Hnodes=getHnodes(EukTree) #get Nodes height 
unifracs <- aperm(GUniFrac(t(Euk_OTU_t_raref), EukTree, alpha=c(0, 0.5, 1))$unifracs,c(3,1,2)) #Unifrac
BDTT_results=BDTT(similarity_slices=slices,tree=EukTree,sampleOTUs=Euk_OTU_t_raref,onlyBeta=T,metric=c("jac","bc"))

saveRDS(BDTT_results,file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Beta_diversity/BDTT_Rarefs1000_18S.RDS")
saveRDS(unifracs,file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Beta_diversity/unifracs_Rarefs1000_18S.RDS")


