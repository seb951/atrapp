#!/usr/local/bin/Rscript

setwd("/Users/jerry/Dropbox/CSBQ/shapiro")

#packages
library(ggplot2)
library(gridExtra)


system("ls -1 results/org_results/*RefSeq_annot_organism.tsv >all.tsv")
all.tsv = read.table("all.tsv", stringsAsFactors = F)
all_spp = NULL
all_spp_shortnames = NULL
all_spp_five = NULL
#first loop is to subset the major species (over 5%)
for(i in 1:nrow(all.tsv))
{
refseq.all = read.table(all.tsv[i,1],sep = "\t",stringsAsFactors = F)
all_spp = c(all_spp,refseq.all[1:10,3]) #keep only top 10
all_spp_five = c(all_spp_five,refseq.all[refseq.all[,1]>5,3]) #keep only the ones that are at more than 5%
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
  all_spp_shortnames = c(all_spp_shortnames,rep(paste(date,location,replicate,sep = "_"),nrow(all_spp_m)))
}
all_spp_m_ggplot = data.frame(unlist(all_spp_m[,2:66]))
temp = all_spp_m_ggplot
all_spp_m_ggplot[,2] = rep(all_spp_m[,1],65)
all_spp_m_ggplot[,3] = all_spp_shortnames
colnames(all_spp_m_ggplot) = c("fraction","species","sample")

#top12
all_spp_m_ggplot_top12 = NULL

for(i in 1:12)
{
  all_spp_m_ggplot_top12 = rbind(all_spp_m_ggplot_top12,all_spp_m_ggplot[all_spp_m_ggplot[,2] == all_spp_five_m[i],])
}
#these are percentage values
all_spp_m_ggplot_top12[,1] = all_spp_m_ggplot_top12[,1]/100

#split into 3 locations.
all_spp_m_ggplot_top12_st1 = all_spp_m_ggplot_top12[regexpr("St1",all_spp_m_ggplot_top12[,3])>0,]
all_spp_m_ggplot_top12_st2 = all_spp_m_ggplot_top12[regexpr("St2",all_spp_m_ggplot_top12[,3])>0,]
all_spp_m_ggplot_top12_PRM = all_spp_m_ggplot_top12[regexpr("PRM",all_spp_m_ggplot_top12[,3])>0,]


###plot - St1
p1=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampSt1 samples)",fill = "Taxonomy") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_bar(aes(y = fraction, x = sample, fill = species),
           data = all_spp_m_ggplot_top12_st1,stat="identity") + ylab("fraction of annotated species") + scale_fill_brewer(palette = "Set3") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#scale_fill_brewer(palette = "Blues") +

###plot - St2
p2=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampSt2 samples)",fill = "Taxonomy") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_bar(aes(y = fraction, x = sample, fill = species),
           data = all_spp_m_ggplot_top12_st2,stat="identity") + ylab("fraction of annotated species")  + scale_fill_brewer(palette = "Set3") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

###plot - PRM
p3=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampPRM samples)",fill = "Taxonomy") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_bar(aes(y = fraction, x = sample, fill = species),
           data = all_spp_m_ggplot_top12_PRM,stat="identity") + ylab("fraction of annotated species") + scale_fill_brewer(palette = "Set3") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#+ scale_fill_brewer(palette = "Blues")

#all three (dimensions in inches)
dev.new(width=14, height=15,noRStudioGD = TRUE)
#pdf('figures/Champlain_barplot.pdf',width=14, height=15)
grid.arrange(p1,p2,p3, ncol = 1)
dev.print(device=pdf, "figures/Champlain_species_barplot.pdf", onefile=FALSE)
dev.off()

