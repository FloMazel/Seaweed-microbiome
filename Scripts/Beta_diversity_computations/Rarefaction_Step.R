#Rarefaction

rm(list=ls())
library(biomformat)
library(vegan)

setwd("~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Full_OTUs_Tables/")

#--------------#
#  16 S data   #
#--------------#

Bacteria_OTU_t=t(as.matrix(biom_data(read_biom("WB_16s_OTUs_NoSeagrass.biom"))))

hist(apply(Bacteria_OTU_t,1,sum))
summary(apply(Bacteria_OTU_t,1,sum))

Bacteria_OTU_t_raref=t(rrarefy(Bacteria_OTU_t,sample = 1000))
colnames(Bacteria_OTU_t_raref)[apply(Bacteria_OTU_t_raref,2,sum)<1000] #"Odonthalia.471"
Bacteria_OTU_t_raref=Bacteria_OTU_t_raref[,apply(Bacteria_OTU_t_raref,2,sum)==1000] #removing sampels with ≤ 999 reads

write.table(Bacteria_OTU_t_raref,file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Rarefied_Full_OTUs_Tables/WB_16s_OTUs_NoSeagrass_raref1000.txt")

#--------------#
#   18S data   #
#--------------#

Euk_OTU_t=t(as.matrix(biom_data(read_biom("WB_18s_OTUs_NoSeagrass.biom"))))

hist(apply(Euk_OTU_t,1,sum))
summary(apply(Euk_OTU_t,1,sum))

Euk_OTU_t_raref=t(rrarefy(Euk_OTU_t,sample = 1000))
write.table(Euk_OTU_t_raref,file="~/Desktop/Recherche/En_cours/Analyses_en_cours/WestBeach/Data/V3/Rarefied_Full_OTUs_Tables/WB_18_OTUs_NoSeagrass_raref1000.txt")
