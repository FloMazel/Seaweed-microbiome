####### -----------------------------------------####### 
#######  Effect of substrate and tidal height    ####### 
####### -----------------------------------------#######

rm(list=ls())

# 1. Loading data and functions 
# -----------------------------

#Loading R packages and some functions
library(vegan)
library(dplyr)
library(ggplot2)
library(ape)
source("~/Documents/GitHub/Sea_weed_Microbiome/Scripts/Useful_functions.R")

#Define my ggplot theme
MyTheme=theme_bw()+theme(axis.text=element_text(size=9),axis.text.x =element_text(size=9),axis.title=element_text(size=14),legend.text=element_text(size=11),axis.ticks = element_line(size = 1.5))

#Loading Metadata
Metadata=read.table("~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Metadata/Westbeach_Metadata_16s.txt",header=T,row.names = "SampleID");Metadata$SampleID=rownames(Metadata)

#Loading Beta-diversity matrices
Bacterial_BDTT=readRDS(file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Beta_diversity/BDTT_Rarefs1000_16S.RDS")
Euk_BDTT=readRDS(file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Beta_diversity/BDTT_Rarefs1000_18S.RDS")

#Selectiong samples that are shared between 16S and 18S data only
shared_sample=intersect(rownames(Metadata),intersect(dimnames(Bacterial_BDTT)[[3]],dimnames(Euk_BDTT)[[3]]))
Metadata=Metadata[shared_sample,]


# 2.1 Statistical tests: raw OTUS 
# ------------------------------

#Selecting the tips of the ASV phylogenies (i.e. the ASV themselves) as microbial units and bray-curtis ('bc') as a beta-diversity metric
Bray_Bacterial=Bacterial_BDTT['0','bc',shared_sample,shared_sample]
Bray_Euk=Euk_BDTT['0','bc',shared_sample,shared_sample]

# 16S -- all samples 
Permanova_16S=adonis2(Bray_Bacterial~Tidal_height+Substrate,data=Metadata,by="margin")
betadisper_16S_substrate=anova(betadisper(as.dist(Bray_Bacterial),Metadata$Substrate)) # higly significant 

#Pairwise PERMANOVA
Pairwise_Permanova_16S=My_pairwise.perm.manova(resp=as.dist(Bray_Bacterial),fact=Metadata$Substrate)


# 18S -- all samples 
Permanova_18S=adonis2(Bray_Euk~Tidal_height+Substrate,data=Metadata,by="margin")
betadisper_16S_substrate=anova(betadisper(as.dist(Bray_Euk),Metadata$Substrate)) # higly significant 
#Pairwise PERMANOVA
Pairwise_Permanova_18S=My_pairwise.perm.manova(resp=as.dist(Bray_Euk),fact=Metadata$Substrate)

# Create a balanced design to test to remove bias from dispersion effect in adonis2
table(Metadata$Substrate) # 25 samples for water
Metadata_Balanced=Metadata %>% group_by(Substrate) %>% sample_n(size = 25) %>% data.frame()
Bal_sample=Metadata_Balanced$SampleID

Permanova_16S_bal=adonis2(Bray_Bacterial[Bal_sample,Bal_sample]~Tidal_height+Substrate,data=Metadata_Balanced,by="margin")
Permanova_18S_bal=adonis2(Bray_Euk[Bal_sample,Bal_sample]~Tidal_height+Substrate,data=Metadata_Balanced,by="margin")

# Export the anova tables 
All_tables=rbind(format_aov_table(Permanova_16S,"16S_allSamples"),format_aov_table(Permanova_16S_bal,"16S_balanced"),format_aov_table(Permanova_18S,"18S_allSamples"),format_aov_table(Permanova_18S_bal,"18S_balanced"))
write.csv(All_tables,"Documents/GitHub/Sea_weed_Microbiome/Outputs/Tables/Substrate_TidalHeight_PerMANOVAs.csv")
write.csv(Pairwise_Permanova_16S,"/Users/fmazel/Documents/GitHub/Sea_weed_Microbiome/Outputs/Tables/Pairwise_Substrate_TidalHeight_PerMANOVAs_16S.csv")
write.csv(Pairwise_Permanova_18S,"/Users/fmazel/Documents/GitHub/Sea_weed_Microbiome/Outputs/Tables/Pairwise_Substrate_TidalHeight_PerMANOVAs_18S.csv")

# 2.2 Plotting  raw OTUS 
# -----------------------------

Pcoa_B=pcoa(Bray_Bacterial)
Pcoa_Euk=pcoa(Bray_Euk)

PC_scores=data.frame(Axe1=Pcoa_B$vectors[,1],Axe2=Pcoa_B$vectors[,2],Morphology=Metadata$morphology_revised2018,Gene="16S")
PC_scores=rbind(PC_scores,data_frame(Axe1=Pcoa_Euk$vectors[,1],Axe2=Pcoa_Euk$vectors[,2],Morphology=Metadata$morphology_revised2018,Gene="18S"))

SubstratePlot=ggplot(aes(y=Axe1,x=Axe2,colour=Morphology),data=PC_scores)+geom_point()+MyTheme+ theme(legend.text=element_text(size=10))+facet_wrap(~Gene)

pdf("Documents/GitHub/Sea_weed_Microbiome/Outputs/PcoA_plots/Substrate_type_PcoA_BrayC.pdf",width = 12,height = 5)
SubstratePlot
dev.off()

# 3.1. Statistical tests: BDTT 
# ----------------------------
# Run the Permanova tests on different microbial phylogenetic resolutions

slices=dimnames(Bacterial_BDTT)[[1]]

Permanova_substrate_BDTT_B=list()
Permanova_substrate_BDTT_Euk=list()

for (i in slices)
{
  print(i)
  Bray_Bacterial=Bacterial_BDTT[i,'bc',shared_sample,shared_sample]
  Bray_Euk=Euk_BDTT[i,'bc',shared_sample,shared_sample]
  Permanova_substrate_BDTT_B[[i]]=adonis2(Bray_Bacterial~Tidal_height*Substrate,data=Metadata)
  Permanova_substrate_BDTT_Euk[[i]]=adonis2(Bray_Euk~Tidal_height*Substrate,data=Metadata) 
}



# 3.2. Plotting BDTT
# -------------------

#format results to plot in ggplot

F_values_summary=expand.grid(predictor=c("Tidal_height","Substrate","Substrate:Tidal_height"),slices=slices,Data=c("16S","18S"))
F_values_summary$F_value=F_values_summary$P_value=NA

Permanova_substrate_BDTT_B_Fvalues=Permanova_substrate_BDTT_Euk_Fvalues=list()
for (i in slices)
{
  F_values_summary[(F_values_summary$slices==i)&(F_values_summary$Data=="16S"),"F_value"]=Permanova_substrate_BDTT_B[[i]]$F[1:3]
  F_values_summary[F_values_summary$slices==i&F_values_summary$Data=="16S","P_value"]=Permanova_substrate_BDTT_B[[i]]$`Pr(>F)`[1:3]
  F_values_summary[F_values_summary$slices==i&F_values_summary$Data=="18S","F_value"]=Permanova_substrate_BDTT_Euk[[i]]$F[1:3]
  F_values_summary[F_values_summary$slices==i&F_values_summary$Data=="18S","P_value"]=Permanova_substrate_BDTT_Euk[[i]]$`Pr(>F)`[1:3]
}


x=10
SubstratePlot=ggplot(aes(y=F_value,x=slices,colour=predictor,shape=P_value<.05,group=predictor),data=F_values_summary)+geom_point(cex=5)+MyTheme+ theme(legend.text=element_text(size=x))+facet_wrap(~Data)+scale_shape_manual(values=c(16, 1 ))
SubstratePlot=SubstratePlot+stat_summary(geom="line")+xlab("Phylogenetic scale")+ggtitle("Effect of substrate and tidal height on microbiota composition")
SubstratePlot

pdf("~/Documents/GitHub/Sea_weed_Microbiome/Outputs/F_value_plots/Substrate_TidalHeight_BrayC_BDTT.pdf",width = 12,height = 5)
SubstratePlot
dev.off()



