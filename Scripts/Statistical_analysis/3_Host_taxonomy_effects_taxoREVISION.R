# ----------------------------#
# ----------------------------#
#    Effect of Host Taxonomy  #
# ----------------------------#
# ----------------------------#

rm(list=ls())

# -----------------------------
# 1. Loading data and functions 
# -----------------------------

#Loading R packages
library(vegan)
library(dplyr)
library(ggplot2)
library(ape)
source("~/Documents/GitHub/Sea_weed_Microbiome/Scripts/Useful_functions.R")


#Loading REVISED TAXONOMY Metadata
Metadata=read.table("~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Metadata/Westbeach_Metadata_16s_REVISED_TAXONOMY.txt",header=T,stringsAsFactors = F)
#Metadata=read.table("~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Metadata/Westbeach_Metadata_16s.txt",header=T,stringsAsFactors = F)

rownames(Metadata)=Metadata$SampleID

#Loading Beta-diversity matrices
Bacterial_BDTT=readRDS(file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Beta_diversity/BDTT_Rarefs1000_16S.RDS")
Euk_BDTT=readRDS(file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Beta_diversity/BDTT_Rarefs1000_18S.RDS")

# Selectiong samples that are shared between 16S and 18S data only
shared_sample=intersect(rownames(Metadata),intersect(dimnames(Bacterial_BDTT)[[3]],dimnames(Euk_BDTT)[[3]]))


# And without water/rock samples
Host_samples=intersect(shared_sample,Metadata$SampleID[!Metadata$Species%in%c("Water","Rock")])
length(Host_samples)

#Selecting the tips of the ASV phylogenies (i.e. the ASV themselves) as microbial units and bray-curtis ('bc') as a beta-diversity metric
Metadata_Host_raw=Metadata[Host_samples,]

# Defining units used in the analysis
Taxonomy=c("Species","Genus","Family","Order","Class","Phylum")
data_units=c("Individuals","Species","Genus","Family","Order","Class")

#Theme for plots 
MyTheme=theme_bw()+theme(axis.text=element_text(size=9),axis.text.x =element_text(size=9),axis.title=element_text(size=14),legend.text=element_text(size=11),axis.ticks = element_line(size = 1.5)) #Define my ggplot theme


#Define the microbial phylogenetic slices 
slices=dimnames(Bacterial_BDTT)[[1]]
#slices="0"

# -------------
# 2. Analysis 
# -------------

# 2.1. Main analysis looped over slices
# -------------------------------------

# 2.1.1. Loop 

Metadata_Host=Bray_B_HostTaxo=Bray_Euk_HostTaxo=Taxo_groups=Host_Permanova_B=Host_Permanova_Euk=Pcoa_B=Pcoa_Euk=TaxoPlot=list()

for (s in slices)
{
  Bray_Bacterial=Bacterial_BDTT[s,'bc',Host_samples,Host_samples]
  Bray_Euk=Euk_BDTT[s,'bc',Host_samples,Host_samples]
  
  # Create list of objects
  Metadata_Host[[s]]=Bray_B_HostTaxo[[s]]=Bray_Euk_HostTaxo[[s]]=Taxo_groups[[s]]=Host_Permanova_B[[s]]=Host_Permanova_Euk[[s]]=Pcoa_B[[s]]=Pcoa_Euk[[s]]=TaxoPlot[[s]]=list()
  
  #Initialize metadata 
  Metadata_Host[[s]][["Individuals"]]=Metadata_Host_raw[,c("Position_m",Taxonomy)]
  Bray_B_HostTaxo[[s]][["Individuals"]]=Bray_Bacterial
  Bray_Euk_HostTaxo[[s]][["Individuals"]]=Bray_Euk
  
  for (n in 1:5) #loop over taxonomic ranks 
  {
    
    #define taxonomic units
    k=data_units[n-1] # sub-units (k) to combine  in units (i)
    i=data_units[n] # units (i) 
    j=Taxonomy[n] # group to test (j)
    
    # Create corresponding metadata
    Metadata_Host[[s]][[j]]=cbind(data.frame(Metadata_Host[[s]][[i]] %>% group_by_(as.name(j)) %>% summarise(total=n())),data.frame(Metadata_Host[[s]][[i]][,-1] %>% group_by_(as.name(j)) %>% summarise_all(funs(unique))))
    row.names(Metadata_Host[[s]][[j]])=Metadata_Host[[s]][[j]][,j]; Metadata_Host[[s]][[j]]=Metadata_Host[[s]][[j]][,-c(1,3)]
    
    # Average Beta_diversity for units different from individuals
    if (!i=="Individuals") {
      Bray_B_HostTaxo[[s]][[i]]=as.matrix(as.dist(meandist(as.dist(Bray_B_HostTaxo[[s]][[k]]),group=as.character(Metadata_Host[[s]][[k]][rownames(Bray_B_HostTaxo[[s]][[k]]),i]))))
      Bray_Euk_HostTaxo[[s]][[i]]=as.matrix(as.dist(meandist(as.dist(Bray_Euk_HostTaxo[[s]][[k]]),group=as.character(Metadata_Host[[s]][[k]][rownames(Bray_B_HostTaxo[[s]][[k]]),i]))))
    } 
    
    # Select groups and units of the analysis (i.e. groups with more than 1 unit, e.g. genera with more than one Species)
    Taxo_groups[[s]][[j]]=rownames(Metadata_Host[[s]][[j]])[Metadata_Host[[s]][[j]]$total>1]
    samplesUnit=rownames(Metadata_Host[[s]][[i]])[Metadata_Host[[s]][[i]][,j]%in%Taxo_groups[[s]][[j]]]
    
    # Permanova test all samples
    Host_taxonomy=Metadata_Host[[s]][[i]][samplesUnit,j]
    
    Bray_curtis=Bray_B_HostTaxo[[s]][[i]][samplesUnit,samplesUnit]
    Host_Permanova_B[[s]][[j]]=adonis2(Bray_curtis~Host_taxonomy)
    
    Bray_curtis=Bray_Euk_HostTaxo[[s]][[i]][samplesUnit,samplesUnit]
    Host_Permanova_Euk[[s]][[j]]=adonis2(Bray_curtis~Host_taxonomy)
    
    # Permanova test low Zone samples
    
    
    #Corresponding PcoA plot
    Pcoa_B[[s]][[j]]=pcoa(Bray_B_HostTaxo[[s]][[i]][samplesUnit,samplesUnit])
    Pcoa_Euk[[s]][[j]]=pcoa(Bray_Euk_HostTaxo[[s]][[i]][samplesUnit,samplesUnit])
    
    PC_scores=data.frame(Taxonomic_ranks=j,Axe1=Pcoa_B[[s]][[j]]$vectors[,1],Axe2=Pcoa_B[[s]][[j]]$vectors[,2],Taxonomic_groups=Metadata_Host[[s]][[i]][samplesUnit,j],Gene="16S")
    PC_scores=rbind(PC_scores,data_frame(Taxonomic_ranks=j,Axe1=Pcoa_Euk[[s]][[j]]$vectors[,1],Axe2=Pcoa_Euk[[s]][[j]]$vectors[,2],Taxonomic_groups=Metadata_Host[[s]][[i]][samplesUnit,j],Gene="18S"))
    
    x=10
    if (i=="Individuals"){x=2}
    TaxoPlot[[s]][[j]]=ggplot(aes(y=Axe1,x=Axe2,colour=Taxonomic_groups),data=PC_scores)+geom_point()+MyTheme+ theme(legend.text=element_text(size=x))+facet_wrap(~Gene)+ggtitle(j)
  } 
  print(s)
}

# SAVE OTUs and metadata tables (for the morphological analysis)
saveRDS(Bray_B_HostTaxo,"Documents/GitHub/Sea_weed_Microbiome/Outputs/Intermediate_Rfiles/OTUs_Bacteria_table_byHost_Taxonomic_rank_BDTT.RDS")
saveRDS(Bray_Euk_HostTaxo,"Documents/GitHub/Sea_weed_Microbiome/Outputs/Intermediate_Rfiles/OTUs_Euk_table_byHost_Taxonomic_rank_BDTT.RDS")
saveRDS(Metadata_Host,"Documents/GitHub/Sea_weed_Microbiome/Outputs/Intermediate_Rfiles/Metadata_byHost_Taxonomic_rank_BDTT.RDS")


Bray_B_HostTaxo=readRDS("Documents/GitHub/Sea_weed_Microbiome/Outputs/Intermediate_Rfiles/OTUs_Bacteria_table_byHost_Taxonomic_rank_BDTT.RDS")
Bray_Euk_HostTaxo=readRDS("Documents/GitHub/Sea_weed_Microbiome/Outputs/Intermediate_Rfiles/OTUs_Euk_table_byHost_Taxonomic_rank_BDTT.RDS")
Metadata_Host=readRDS("~/Documents/GitHub/Sea_weed_Microbiome/Outputs/Intermediate_Rfiles/Metadata_byHost_Taxonomic_rank_BDTT.RDS")


# 2.1.2. GGplot formatting
testGroup=Taxonomy[-6]
F_values_summary=expand.grid(slices=slices,Taxonomy=testGroup,Data=c("16S","18S"))
F_values_summary$F_value=F_values_summary$P_value=NA

for (s in slices)
{
  for (i in testGroup)
  {
    F_values_summary[(F_values_summary$slices==s)&(F_values_summary$Taxonomy==i)&(F_values_summary$Data=="16S"),"F_value"]=Host_Permanova_B[[s]][[i]]$F[1]
    F_values_summary[(F_values_summary$slices==s)&F_values_summary$Taxonomy==i&F_values_summary$Data=="16S","P_value"]=Host_Permanova_B[[s]][[i]]$`Pr(>F)`[1]
    F_values_summary[(F_values_summary$slices==s)&F_values_summary$Taxonomy==i&F_values_summary$Data=="18S","F_value"]=Host_Permanova_Euk[[s]][[i]]$F[1]
    F_values_summary[(F_values_summary$slices==s)&F_values_summary$Taxonomy==i&F_values_summary$Data=="18S","P_value"]=Host_Permanova_Euk[[s]][[i]]$`Pr(>F)`[1]
  }
}

# 2.2 Sensitivity analysis on one unique tidal zone (low)
# --------------------------------------------------------

# 2.2.1 Loop
Low_samples=as.character(Metadata_Host_raw$SampleID[Metadata_Host_raw$Tidal_height=="low"])
Low_Metadata_Host=Low_Bray_B_HostTaxo=Low_Bray_Euk_HostTaxo=Low_Taxo_groups=Low_Host_Permanova_B=Low_Host_Permanova_Euk=Low_Pcoa_B=Low_Pcoa_Euk=Low_TaxoPlot=list()

for (s in c('0')) #only for classical OTUs
{
  Low_Bray_Bacterial=Bacterial_BDTT[s,'bc',Low_samples,Low_samples]
  Low_Bray_Euk=Euk_BDTT[s,'bc',Low_samples,Low_samples]
  
  # Create list of objects
  Low_Metadata_Host[[s]]=Low_Bray_B_HostTaxo[[s]]=Low_Bray_Euk_HostTaxo[[s]]=Low_Taxo_groups[[s]]=Low_Host_Permanova_B[[s]]=Low_Host_Permanova_Euk[[s]]=Low_Pcoa_B[[s]]=Low_Pcoa_Euk[[s]]=Low_TaxoPlot[[s]]=list()
  
  #Initialize metadata 
  Low_Metadata_Host[[s]][["Individuals"]]=Metadata_Host_raw[Low_samples,c("Position_m",Taxonomy)]
  Low_Bray_B_HostTaxo[[s]][["Individuals"]]=  Low_Bray_Bacterial
  Low_Bray_Euk_HostTaxo[[s]][["Individuals"]]=  Low_Bray_Euk
  
  for (n in 1:5) #loop over taxonomic ranks 
  {
    
    #define taxonomic units
    k=data_units[n-1] # sub-units (k) to combine  in units (i)
    i=data_units[n] # units (i) 
    j=Taxonomy[n] # group to test (j)
    
    # Create corresponding metadata
    Low_Metadata_Host[[s]][[j]]=cbind(data.frame(Low_Metadata_Host[[s]][[i]] %>% group_by_(as.name(j)) %>% summarise(total=n())),data.frame(Low_Metadata_Host[[s]][[i]][,-1] %>% group_by_(as.name(j)) %>% summarise_all(funs(unique))))
    row.names(Low_Metadata_Host[[s]][[j]])=Low_Metadata_Host[[s]][[j]][,j]; Low_Metadata_Host[[s]][[j]]=Low_Metadata_Host[[s]][[j]][,-c(1,3)]
    
    # Average Beta_diversity for units different from individuals
    if (!i=="Individuals") {
      Low_Bray_B_HostTaxo[[s]][[i]]=as.matrix(as.dist(meandist(as.dist(Low_Bray_B_HostTaxo[[s]][[k]]),group=as.character(Low_Metadata_Host[[s]][[k]][rownames(Low_Bray_B_HostTaxo[[s]][[k]]),i]))))
      Low_Bray_Euk_HostTaxo[[s]][[i]]=as.matrix(as.dist(meandist(as.dist(Low_Bray_Euk_HostTaxo[[s]][[k]]),group=as.character(Low_Metadata_Host[[s]][[k]][rownames(Low_Bray_B_HostTaxo[[s]][[k]]),i]))))
    } 
    
    # Select groups and units of the analysis (i.e. groups with more than 1 unit, e.g. genera with more than one Species)
    Low_Taxo_groups[[s]][[j]]=rownames(Low_Metadata_Host[[s]][[j]])[Low_Metadata_Host[[s]][[j]]$total>1]
    samplesUnit=rownames(Low_Metadata_Host[[s]][[i]])[Low_Metadata_Host[[s]][[i]][,j]%in%Low_Taxo_groups[[s]][[j]]]
    
    # Permanova test all samples
    Low_Host_taxonomy=Low_Metadata_Host[[s]][[i]][samplesUnit,j]
    
    Low_Bray_curtis=Low_Bray_B_HostTaxo[[s]][[i]][samplesUnit,samplesUnit]
    Low_Host_Permanova_B[[s]][[j]]=adonis2(Low_Bray_curtis~Low_Host_taxonomy)
    
    Low_Bray_curtis=Low_Bray_Euk_HostTaxo[[s]][[i]][samplesUnit,samplesUnit]
    Low_Host_Permanova_Euk[[s]][[j]]=adonis2(Low_Bray_curtis~Low_Host_taxonomy)
    
    #Corresponding PcoA plot
    Low_Pcoa_B[[s]][[j]]=pcoa(Low_Bray_B_HostTaxo[[s]][[i]][samplesUnit,samplesUnit])
    Low_Pcoa_Euk[[s]][[j]]=pcoa(Low_Bray_Euk_HostTaxo[[s]][[i]][samplesUnit,samplesUnit])
    
    PC_scores=data.frame(Taxonomic_ranks=j,Axe1=Low_Pcoa_B[[s]][[j]]$vectors[,1],Axe2=Low_Pcoa_B[[s]][[j]]$vectors[,2],Taxonomic_groups=Low_Metadata_Host[[s]][[i]][samplesUnit,j],Gene="16S")
    PC_scores=rbind(PC_scores,data_frame(Taxonomic_ranks=j,Axe1=Low_Pcoa_Euk[[s]][[j]]$vectors[,1],Axe2=Low_Pcoa_Euk[[s]][[j]]$vectors[,2],Taxonomic_groups=Low_Metadata_Host[[s]][[i]][samplesUnit,j],Gene="18S"))
    
    x=10
    if (i=="Individuals"){x=2}
    Low_TaxoPlot[[s]][[j]]=ggplot(aes(y=Axe1,x=Axe2,colour=Taxonomic_groups),data=PC_scores)+geom_point()+MyTheme+ theme(legend.text=element_text(size=x))+facet_wrap(~Gene)+ggtitle(j)
  } 
  print(s)
}

# 2.2.2. GGplot formatting
testGroup=Taxonomy[-6]
F_values_summaryLow=expand.grid(slices=slices,Taxonomy=testGroup,Data=c("16S","18S"))
F_values_summaryLow$F_value=F_values_summaryLow$P_value=NA

for (s in c('0'))
{
  for (i in testGroup)
  {
    F_values_summaryLow[(F_values_summaryLow$slices==s)&(F_values_summaryLow$Taxonomy==i)&(F_values_summaryLow$Data=="16S"),"F_value"]=Low_Host_Permanova_B[[s]][[i]]$F[1]
    F_values_summaryLow[(F_values_summaryLow$slices==s)&F_values_summaryLow$Taxonomy==i&F_values_summaryLow$Data=="16S","P_value"]=Low_Host_Permanova_B[[s]][[i]]$`Pr(>F)`[1]
    F_values_summaryLow[(F_values_summaryLow$slices==s)&F_values_summaryLow$Taxonomy==i&F_values_summaryLow$Data=="18S","F_value"]=Low_Host_Permanova_Euk[[s]][[i]]$F[1]
    F_values_summaryLow[(F_values_summaryLow$slices==s)&F_values_summaryLow$Taxonomy==i&F_values_summaryLow$Data=="18S","P_value"]=Low_Host_Permanova_Euk[[s]][[i]]$`Pr(>F)`[1]
  }
}


# 2.3. Sensitivity analysis for balanced data on raw OTUs

Sp=names(table(Metadata_Host_raw$Species))[table(Metadata_Host_raw$Species)>1]
Metadata_Host_Balanced=Metadata_Host_raw[Metadata_Host_raw$Species%in%Sp,] %>% group_by(Species) %>% sample_n(size = 2) %>% data.frame()
Bal_sample=as.character(Metadata_Host_Balanced$SampleID)
Permanova_Bacteria_BalSamp=adonis2(Bacterial_BDTT['0','bc',Bal_sample,Bal_sample]~Species,data=Metadata_Host_Balanced)
Permanova_Euk_BalSamp=adonis2(Euk_BDTT['0','bc',Bal_sample,Bal_sample]~Species,data=Metadata_Host_Balanced)





# ----------------------------
# 3. Plot and save the objects 
# ----------------------------

# 3.1. Plot and table for slice "0": tips of the phylogenies
# -------------------------------------------------------------

# Plot within vs between taxonomic ranks with raw beta-values
library(reshape2)
library(grid)
vp <- viewport(width = 0.4, height = 0.4, x = 0.75, y = 0.7)
s="0"
n=1
LL=c(expression("Within \n\n    Between"),expression("Within \n\n    Between"),expression("Within \n\n    Between"),expression("Within \n\n    Between"),expression("Within \n\n    Between"))
MyTheme1=theme_bw()+theme(axis.text.x=element_text(size=9,angle=45, vjust=.5),axis.title=element_text(size=14),legend.text=element_text(size=11),axis.ticks = element_line(size = 1.5)) #Define my ggplot theme

#Bacteria 
Compa=list()
for (n in 1:5)
{
  i=data_units[n] # units (i) 
  j=Taxonomy[n] # group to test (j)

  t=rownames(Metadata_Host[[s]][[j]])[Metadata_Host[[s]][[j]]$total>1]
  samplesUnit=rownames(Metadata_Host[[s]][[i]])[Metadata_Host[[s]][[i]][,j]%in%t]
  Host_taxonomy=Metadata_Host[[s]][[i]][samplesUnit,j]
  
  Bray_curtis=as.data.frame(melt(Bray_B_HostTaxo[[s]][[i]][samplesUnit,samplesUnit]))
  Bray_curtis=Bray_curtis[!Bray_curtis$Var1==Bray_curtis$Var2,] #remove self comparisons
  Bray_curtis$F1=Metadata_Host[[s]][[i]][as.character(Bray_curtis$Var1),j];Bray_curtis$F2=Metadata_Host[[s]][[i]][as.character(Bray_curtis$Var2),j]
  Bray_curtis$Type=Bray_curtis$F2==Bray_curtis$F1
  Bray_curtis$Type[Bray_curtis$Type==T]="Within"
  Bray_curtis$Type[Bray_curtis$Type==F]="Between"
  Bray_curtis$Rank=j
  
  Compa[[j]]=Bray_curtis[,c("Type","value","Rank")]
}

SummaryBeta=do.call(rbind,Compa)
SummaryBeta$Rank <- factor(SummaryBeta$Rank, levels=unique(SummaryBeta$Rank))
SummaryBeta$Type=factor(SummaryBeta$Type,levels=c("Within","Between"))

Bac_RawBeta=ggplot(SummaryBeta, aes(x=Rank, y=1-value, fill=Type)) + 
  geom_boxplot(show.legend = FALSE) +
  scale_x_discrete(name = NULL, labels = LL) + 
  scale_y_continuous(name = "Microbiome similarity (1 - Bray-Curtis)") +
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  ggtitle("(A) Bacteria")+
  MyTheme1
Bac_RawBeta

#Add f ratio
F_ratio_B=ggplot(aes(y=F_value,x=Taxonomy,shape=P_value<.05,group=Data),data=F_values_summary[(F_values_summary$slices==s)&(as.character(F_values_summary$Data)=="16S"),])+geom_point(size=4,show.legend = FALSE)+scale_shape_manual(values=c(1, 16 ))+MyTheme1+geom_line()+xlab("")+ylab("Pseudo F")

pdf("Documents/GitHub/Sea_weed_Microbiome/Outputs/F_value_plots/Host_Taxonomy_Fvalue_Beta_Bact_BrayC.pdf",width = 6,height = 5)
print(Bac_RawBeta)
print(F_ratio_B, vp = vp)
dev.off()

#Eukariotes
Compa=list()
for (n in 1:5)
{
  i=data_units[n] # units (i) 
  j=Taxonomy[n] # group to test (j)
  
  t=rownames(Metadata_Host[[s]][[j]])[Metadata_Host[[s]][[j]]$total>1]
  samplesUnit=rownames(Metadata_Host[[s]][[i]])[Metadata_Host[[s]][[i]][,j]%in%t]
  Host_taxonomy=Metadata_Host[[s]][[i]][samplesUnit,j]
  
  Bray_curtis=as.data.frame(melt(Bray_Euk_HostTaxo[[s]][[i]][samplesUnit,samplesUnit]))
  Bray_curtis=Bray_curtis[!Bray_curtis$Var1==Bray_curtis$Var2,] #remove self comparisons
  Bray_curtis$F1=Metadata_Host[[s]][[i]][as.character(Bray_curtis$Var1),j];Bray_curtis$F2=Metadata_Host[[s]][[i]][as.character(Bray_curtis$Var2),j]
  Bray_curtis$Type=Bray_curtis$F2==Bray_curtis$F1
  Bray_curtis$Type[Bray_curtis$Type==T]="Within"
  Bray_curtis$Type[Bray_curtis$Type==F]="Between"
  Bray_curtis$Rank=j
  
  Compa[[j]]=Bray_curtis[,c("Type","value","Rank")]
}

SummaryBeta=do.call(rbind,Compa)
SummaryBeta$Rank <- factor(SummaryBeta$Rank, levels=unique(SummaryBeta$Rank))
SummaryBeta$Type=factor(SummaryBeta$Type,levels=c("Within","Between"))

Euk_RawBeta=ggplot(SummaryBeta, aes(x=Rank, y=1-value, fill=Type)) + 
  geom_boxplot(show.legend = FALSE) +
  scale_x_discrete(name = NULL, labels = LL) + 
  scale_y_continuous(name = "Microbiome similarity (1 - Bray-Curtis)") +
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  ggtitle("(B) Eukaryotes")+
  MyTheme1
Euk_RawBeta

#Add f ratio
F_ratio_Euk=ggplot(aes(y=F_value,x=Taxonomy,shape=P_value<.05,group=Data),data=F_values_summary[(F_values_summary$slices==s)&(as.character(F_values_summary$Data)=="18S"),])+geom_point(size=4,show.legend = FALSE)+scale_shape_manual(values=c(1, 16 ))+MyTheme1+geom_line()+xlab("")+ylab("Pseudo F")

pdf("Documents/GitHub/Sea_weed_Microbiome/Outputs/F_value_plots/Host_Taxonomy_Fvalue_Beta_Euk_BrayC.pdf",width = 6,height = 5)
print(Euk_RawBeta)
print(F_ratio_Euk, vp = vp)
dev.off()

# PcoA Plots 
pdf("Documents/GitHub/Sea_weed_Microbiome/Outputs/PcoA_plots/Host_Taxonomy_PcoA_BrayC.pdf",width = 20,height = 10)
multiplot(TaxoPlot[[s]][[1]],TaxoPlot[[s]][[2]],TaxoPlot[[s]][[3]],TaxoPlot[[s]][[4]],TaxoPlot[[s]][[5]],cols=2)
dev.off()

# F_value along taxomomic scale
pdf("Documents/GitHub/Sea_weed_Microbiome/Outputs/F_value_plots/Host_Taxonomy_Fvalue_BrayC.pdf",width = 8,height = 5)
ggplot(aes(y=F_value,x=Taxonomy,color=Data,shape=P_value<.05,group=Data),data=F_values_summary[(F_values_summary$slices==s),])+geom_point()+scale_shape_manual(values=c(1, 16 ))+MyTheme+stat_summary(geom="line")+ggtitle("Effect of host taxonomy on microbiota composition")
dev.off()

# F_value along taxomomic scale (low zone)
pdf("Documents/GitHub/Sea_weed_Microbiome/Outputs/F_value_plots/Host_Taxonomy_Fvalue_BrayC_lowZone.pdf",width = 8,height = 5)
ggplot(aes(y=F_value,x=Taxonomy,color=Data,shape=P_value<.05,group=Data),data=F_values_summaryLow[(F_values_summaryLow$slices==s),])+geom_point()+scale_shape_manual(values=c(1, 16 ))+MyTheme+stat_summary(geom="line")+ggtitle("Effect of host taxonomy on microbiota composition")
dev.off()

# Permanova tables 
s='0'
B=Euk=list()
for (i in names(Host_Permanova_B[[s]]))
{
  B[[i]]=format_aov_table(Host_Permanova_B[[s]][[i]],"16S")
  Euk[[i]]=format_aov_table(Host_Permanova_Euk[[s]][[i]],"18S")
}

B[["Species_BalancedSamping"]]=format_aov_table(Permanova_Bacteria_BalSamp,"16S")
Euk[["Species_BalancedSamping"]]=format_aov_table(Permanova_Euk_BalSamp,"18S")

write.csv(rbind(do.call(rbind,B),do.call(rbind,Euk)),"Documents/GitHub/Sea_weed_Microbiome/Outputs/Tables/Host_Taxonomy_PerMANOVA.csv",na="")


# 3.3. BDTT Plot 

F_values_summarySIGNIF=F_values_summary[F_values_summary$P_value<.05,]
BDTT_taxo=ggplot(data = F_values_summarySIGNIF, aes(x = Taxonomy, y = slices,colour=F_value)) +geom_point(cex=15)+facet_wrap(~Data)

pdf("Documents/GitHub/Sea_weed_Microbiome/Outputs/F_value_plots/Host_Taxonomy_Fvalue_BrayC_BDTT.pdf",width = 12,height = 5)
BDTT_taxo+xlab("Host taxonomy (explanatory variable)")+ylab("Microbial phylogenetic scale")+scale_color_gradient(low="grey", high="red")
dev.off()



