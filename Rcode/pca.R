#PCoA


setwd("/Users/jerry/Dropbox/CSBQ/shapiro")
 
#packages
library(pals)
library(ggplot2)
library(gridExtra)
library(vegan)
library(ape)
library(lubridate)

system("ls -1 results/org_results/*RefSeq_annot_organism.tsv >all.tsv")
all.tsv = read.table("all.tsv", stringsAsFactors = F)
all_spp = NULL
all_spp_shortnames = NULL
all_spp_five = NULL
#first loop is to subset the major species (over 1%)
for(i in 1:nrow(all.tsv))
{
  refseq.all = read.table(all.tsv[i,1],sep = "\t",stringsAsFactors = F)
  all_spp = c(all_spp,refseq.all[1:100,3]) #keep only top 10
  all_spp_five = c(all_spp_five,refseq.all[refseq.all[,1]>1,3]) #keep only the ones that are at more than 1% in at least one sample
}
all_spp_five_m = unique(sort(all_spp_five))
all_spp_m = data.frame(unique(sort(all_spp)),stringsAsFactors = F)

for(i in 1:nrow(all.tsv))
{
  refseq.all = read.table(all.tsv[i,1],sep = "\t",stringsAsFactors = F)
  for(j in 1:nrow(all_spp_m))
  {
    temp = refseq.all[refseq.all[,3] == all_spp_m[j,1],1]
    if(length(temp) == 1) all_spp_m[j,i+1] = temp
    if(length(temp) == 0) all_spp_m[j,i+1] = 0
  }
  #get a proper shortname for the graph
  replicate = substring(strsplit(all.tsv[i,1],split = "WatPhotz_")[[1]][2],4,4)
  location = substring(strsplit(all.tsv[i,1],split = "Champ")[[1]][2],1,3)
  date = strsplit(all.tsv[i,1],split = "-")[[1]][2]
  colnames(all_spp_m)[i+1] = paste(location,date,replicate,sep = "_")
  colnames(all_spp_m)[1] = "species" 
  all_spp_shortnames = c(all_spp_shortnames,rep(paste(date,location,replicate,sep = "_"),nrow(all_spp_m)))
}

#standardization using hellinger transform
#asv.filt.abundants.norm.hel <-decostand(asv.filt.abundants.norm, "hel")

#PERMANOVA
#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
#permanova = adonis(formula=asv.filt.abundants.norm.hel~fertilization*species,strata=(design.keep$bloc/design.keep$replicate), data=design.keep, permutations=9999, method="bray")
#permanova$aov.tab$comparison = "root_fungi"
#all_spp_m = all_spp_m[,-c(6:8)]
all_spp_m.pcoa = data.frame(all_spp_m[,-1])
rownames(all_spp_m.pcoa) = all_spp_m[,1]

#transpose
all_spp_m.pcoa=t(all_spp_m.pcoa)

#standardize
all_spp_m.pcoa.hel <-decostand(all_spp_m.pcoa, "hel")

#dissimilarity
all_spp_m.pcoa.hel.bray <-vegdist(all_spp_m.pcoa.hel, method="bray")

#Calculating PCoA
all_spp_m.pcoa.hel.bray.pcoa<-pcoa(dist(all_spp_m.pcoa.hel.bray))

#How many axes represent more variability (17)
bs = all_spp_m.pcoa.hel.bray.pcoa$values$Broken_stick
length(bs[bs>mean(bs)])

#PVE of first 2 axes (4.7% & 3.8%)
axis.1.2 = round((all_spp_m.pcoa.hel.bray.pcoa$values$Broken_stick/sum(all_spp_m.pcoa.hel.bray.pcoa$values$Broken_stick))[1:2],4)*100

#Ploting the PCoAs - with fertilization as empty circles
#crops are "darkred","darkblue","darkorange
pcoa.plot = data.frame(all_spp_m.pcoa.hel.bray.pcoa$vectors[,1:2])
pcoa.plot[regexpr("St1",rownames(pcoa.plot))>0,3] = "St1" #circle
pcoa.plot[regexpr("PRM",rownames(pcoa.plot))>0,3] = "PRM" #triangle
pcoa.plot[regexpr("St2",rownames(pcoa.plot))>0,3] = "St2" #square
date = (unlist(strsplit(rownames(pcoa.plot),split = "_"))[seq(2,195,by =3)])
pcoa.plot[,4] = ymd(date)
colnames(pcoa.plot) = c("axis1","axis2","Sampling","date")


p1=ggplot(pcoa.plot,aes(axis1,axis2)) +
  labs(title = "Lake Champlain - PcoA (temp. sps, present at >1%)",shape = "Sampling Site", colour = "Sampling Date") +
  geom_point(aes(colour=as.factor(date),shape = factor(Sampling)),stroke=3) +
  scale_shape(solid = F) +
  ylab(paste("Axis 2 (PVE:",axis.1.2[2],"%)",sep = "")) + xlab(paste("Axis 1 (PVE:",axis.1.2[1],"%)",sep = "")) +
  scale_color_brewer(palette= "Spectral")

#

dev.new()
p1
dev.print(device=pdf, "figures/Champlain_pcoa.pdf", onefile=FALSE)
dev.off()


