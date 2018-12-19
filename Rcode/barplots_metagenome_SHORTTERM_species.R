#!/usr/local/bin/Rscript

setwd("/Users/jerry/Documents/CSBQ/shapiro")

#packages
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(dplyr)
library(lubridate)

system("ls -1 results/Short_time_course_Lac_Champlain/org_results/*RefSeq_annot_organism.tsv >all.tsv")
all.tsv = read.table("all.tsv", stringsAsFactors = F)
all_spp = NULL
dates = NULL
locations = NULL
all_spp_five = NULL
indexes = NULL
#first loop is to subset the major species (over 5%)
for(i in 1:nrow(all.tsv))
{
refseq.all = read.table(all.tsv[i,1],sep = "\t",stringsAsFactors = F)
#all_spp = c(all_spp,refseq.all[1:25,3]) #keep only top 10
#all_spp_five = c(all_spp_five,refseq.all[refseq.all[,1]>2,3]) #keep only the ones that are at more than 5%
#simplify by genus
#all_spp_five = c(all_spp_five,refseq.all[refseq.all[,2]>5,3]) #keep only the ones that are at more than 5%

#simplify by genus
refseq.all$genus = unlist(lapply(strsplit(refseq.all[,3],split = " "), `[[`, 1))
refseq.all = refseq.all[,c(1,2,4)] %>% group_by(genus) %>% summarise_all(sum)
refseq.all = refseq.all[order(refseq.all$V1,decreasing =T),]
all_spp = c(all_spp,refseq.all[1:25,1]) #keep only top 10
all_spp_five = c(all_spp_five,refseq.all[refseq.all[,2]>5,1]) #keep only the ones that are at more than 5%
}
#simplify by genus
all_spp_five_m = unique(sort(unlist(all_spp_five)))
all_spp_m = data.frame(unique(sort(unlist(all_spp))),stringsAsFactors = F)

#all_spp_five_m = unique(sort(all_spp_five))
#all_spp_m = data.frame(unique(sort(all_spp)),stringsAsFactors = F)

for(i in 1:nrow(all.tsv))
{
  refseq.all = read.table(all.tsv[i,1],sep = "\t",stringsAsFactors = F)
  
  #simplify by genus
  refseq.all$genus = unlist(lapply(strsplit(refseq.all[,3],split = " "), `[[`, 1))
  refseq.all = refseq.all[,c(1,2,4)] %>% group_by(genus) %>% summarise_all(sum)
  
  for(j in 1:nrow(all_spp_m))
  {
    #simplify by genus
    temp = refseq.all[refseq.all[,1] == all_spp_m[j,1],2]
    #temp = refseq.all[refseq.all[,3] == all_spp_m[j,1],1]
    if(length(temp) == 1) all_spp_m[j,i+1] = temp
    if(length(temp) == 0) all_spp_m[j,i+1] = 0
  }
  #get a proper shortname for the graph
  shortname = strsplit(all.tsv[i,1],split = "_merged")[[1]][1]
  shortname_date = strsplit(all.tsv[i,1],split = "_")[[1]][7] 
  date = strsplit(shortname_date,split = ".",fixed = T)[[1]][2]
  location = substring(shortname,nchar(shortname)-2,nchar(shortname))
  index = strsplit(shortname,split=".",fixed =T)[[1]][4]
  
  indexes = c(indexes,rep(index,nrow(all_spp_m)))
  dates = c(dates,rep(date,nrow(all_spp_m)))
  locations = c(locations,rep("st1",nrow(all_spp_m))) #shortterm are all st1
}

#create your ggplot object
all_spp_m_ggplot = data.frame(unlist(all_spp_m[,2:24]))
temp = all_spp_m_ggplot
all_spp_m_ggplot[,2] = rep(all_spp_m[,1],23)
all_spp_m_ggplot[,3] = dmy(dates)
all_spp_m_ggplot[,4] = locations
all_spp_m_ggplot[,5] = factor(all_spp_m_ggplot[,4], levels=c('st1'))
all_spp_m_ggplot[,6] = unlist(lapply(strsplit(all_spp_m_ggplot[,2],split = " "), `[[`, 1))
all_spp_m_ggplot[,7] = indexes 
all_spp_m_ggplot[,8] = paste(as.character(all_spp_m_ggplot[,3]),gsub("Index","",all_spp_m_ggplot[,7]),sep = "") 
colnames(all_spp_m_ggplot) = c("fraction","species","samples","locations","locations_f","genus","indexes","unique_ID")

#top percentage according to "all_spp_five_m" vector
all_spp_m_ggplot_top12 = NULL
for(i in 1:length(all_spp_five_m))
{
  all_spp_m_ggplot_top12 = rbind(all_spp_m_ggplot_top12,all_spp_m_ggplot[all_spp_m_ggplot[,2] == all_spp_five_m[i],])
}
#these are percentage values
all_spp_m_ggplot_top12[,1] = all_spp_m_ggplot_top12[,1]/100

#plot
x = colorRampPalette(brewer.pal(12,"Paired"))
p1=ggplot() + labs(title = "Lake Champlain - annotated genera - Short term",fill = "Taxonomy") +
  theme_bw() + 
  #  scale_x_date(date_breaks = "weeks") +
  theme(legend.text = element_text(face="italic"),plot.title = element_text(hjust = 0.5, size=14, face="bold")) + scale_fill_manual(values = x(length(all_spp_five_m))) +
  geom_bar(aes(y = fraction, x = unique_ID, fill = genus),
           data = all_spp_m_ggplot_top12,stat="identity") + ylab("fraction of annotated genera")  + theme(axis.text.x = element_text(angle = 90, hjust = 1))


#PDF (dimensions in inches)
dev.new(width=10, height=7,noRStudioGD = TRUE)
p1 + facet_grid(rows=vars(locations_f))
dev.print(device=pdf,"figures/Metagenomic_Champlain_short_term_genera_barplot.pdf", onefile=FALSE)
dev.off()


png("figures/Metagenomic_Champlain_short_genera_barplot.png",width=10, res =400,height=7,units= 'in')
p1 + facet_grid(rows=vars(locations_f))
dev.off()


#save
write.table(all_spp_m_ggplot_top12,"results/all_spp_m_ggplot_top12_short_genome",row.names = F, col.names = T, quote = F)
