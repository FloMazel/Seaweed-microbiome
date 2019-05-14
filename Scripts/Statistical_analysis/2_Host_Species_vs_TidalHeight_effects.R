# --------------------------------------------#
#    Effect of Host Species vs Tidal height   #
# --------------------------------------------#

rm(list=ls())

# 1. Loading data and functions 
# -----------------------------

#Define my ggplot theme
MyTheme=theme_bw()+theme(axis.text=element_text(size=9),axis.text.x =element_text(size=9),axis.title=element_text(size=14),legend.text=element_text(size=11),axis.ticks = element_line(size = 1.5))

#Loading R packages and some functions
library(vegan)
library(dplyr)
library(ggplot2)
library(ape)
source("~/Documents/GitHub/Sea_weed_Microbiome/Scripts/Useful_functions.R")

#Loading Metadata
Metadata=read.table("~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Metadata/Westbeach_Metadata_16s.txt",header=T)
rownames(Metadata)=Metadata$SampleID

#Loading Beta-diversity matrices
Bacterial_BDTT=readRDS(file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Beta_diversity/BDTT_Rarefs1000_16S.RDS")
Euk_BDTT=readRDS(file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Beta_diversity/BDTT_Rarefs1000_18S.RDS")
dim(Bacterial_BDTT)
dim(Euk_BDTT)

# Selectiong samples that are shared between 16S and 18S data only
shared_sample=intersect(rownames(Metadata),intersect(dimnames(Bacterial_BDTT)[[3]],dimnames(Euk_BDTT)[[3]]))

# And without water/rock samples
Host_samples=intersect(shared_sample,Metadata$SampleID[!Metadata$Species%in%c("Water","Rock")])
Metadata_Host_raw=Metadata[Host_samples,]

# Select those Species that are only found in all tidal heigth to contrast Species host identity Vs. tidal height effects
Species_Oc=apply(table(Metadata$Species,Metadata$Tidal_height)>0,1,sum)
Generalists_sp=names(Species_Oc)[Species_Oc==3]
Generalists_sp=Generalists_sp[!Generalists_sp=="Rock"] #removing rocks
Generalists_sp_sample=as.character(Metadata_Host_raw$SampleID[as.character(Metadata_Host_raw$Species)%in%Generalists_sp])


# 2. Statistical tests: raw OTUS 
# ------------------------------

#Selecting the tips of the ASV phylogenies (i.e. the ASV themselves) as microbial units and bray-curtis ('bc') as a beta-diversity metric
Bray_Bacterial=Bacterial_BDTT['0','bc',Host_samples,Host_samples]
Bray_Euk=Euk_BDTT['0','bc',Host_samples,Host_samples]

Permanova_Bact=adonis2(Bray_Bacterial[Generalists_sp_sample,Generalists_sp_sample]~Species+Tidal_height,data=Metadata_Host_raw[Generalists_sp_sample,],by="margin")
Permanova_Euk=adonis2(Bray_Euk[Generalists_sp_sample,Generalists_sp_sample]~Species+Tidal_height,data=Metadata_Host_raw[Generalists_sp_sample,],by="margin")

All_tables=rbind(format_aov_table(Permanova_Bact,"16S"),format_aov_table(Permanova_Euk,"18S"))
write.csv(All_tables,"Documents/GitHub/Sea_weed_Microbiome/Outputs/Tables/HostSpecies_TidalHeight_PerMANOVAs.csv")

# 3.1. Statistical tests: BDTT
# ------------------------------

slices=dimnames(Bacterial_BDTT)[[1]]
Permanova_BDTT_B=Permanova_BDTT_Euk=list()

for (i in slices)
{
  print(i)
  Bray_Bacterial=Bacterial_BDTT[i,'bc',Generalists_sp_sample,Generalists_sp_sample]
  Bray_Euk=Euk_BDTT[i,'bc',Generalists_sp_sample,Generalists_sp_sample]
  Permanova_BDTT_B[[i]]=adonis2(Bray_Bacterial~Species+Tidal_height,data=Metadata_Host_raw[Generalists_sp_sample,],by="margin")
  Permanova_BDTT_Euk[[i]]=adonis2(Bray_Euk~Species+Tidal_height,data=Metadata_Host_raw[Generalists_sp_sample,],by="margin")
}

# 3.2. Plotting BDTT
# -------------------

#format results to plot in ggplot

F_values_summary=expand.grid(predictor=c("Species","Tidal_height"),slices=slices,Data=c("16S","18S"))
F_values_summary$F_value=F_values_summary$P_value=NA

for (i in slices)
{
  F_values_summary[(F_values_summary$slices==i)&(F_values_summary$Data=="16S"),"F_value"]=Permanova_BDTT_B[[i]]$F[1:2]
  F_values_summary[F_values_summary$slices==i&F_values_summary$Data=="16S","P_value"]=Permanova_BDTT_B[[i]]$`Pr(>F)`[1:2]
  F_values_summary[F_values_summary$slices==i&F_values_summary$Data=="18S","F_value"]=Permanova_BDTT_Euk[[i]]$F[1:2]
  F_values_summary[F_values_summary$slices==i&F_values_summary$Data=="18S","P_value"]=Permanova_BDTT_Euk[[i]]$`Pr(>F)`[1:2]
}

x=10
TidalPlot=ggplot(aes(y=F_value,x=slices,colour=predictor,shape=P_value<.05,group=predictor),data=F_values_summary)+geom_point(cex=5)+MyTheme+ theme(legend.text=element_text(size=x))+facet_wrap(~Data)+scale_shape_manual(values=c(1, 16 ))
TidalPlot=TidalPlot+stat_summary(geom="line")+xlab("Phylogenetic scale")+ggtitle("Effect of substrate and tidal height on microbiota composition")
TidalPlot

pdf("~/Documents/GitHub/Sea_weed_Microbiome/Outputs/F_value_plots/F_value_Tidal_Vs_SPecies_BDTT.pdf",width = 12,height = 5)
TidalPlot
dev.off()
