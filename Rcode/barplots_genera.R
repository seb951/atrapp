#!/usr/local/bin/Rscript

setwd("/Users/jerry/Documents/CSBQ/shapiro")

#packages
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(dplyr)
library(lubridate)

system("ls -1 results/org_results/*RefSeq_annot_organism.tsv >all.tsv")
all.tsv = read.table("all.tsv", stringsAsFactors = F)
all_spp = NULL
samples = NULL
locations = NULL
all_spp_five = NULL
#first loop is to subset the major species (over 5%)
for(i in 1:nrow(all.tsv))
{
refseq.all = read.table(all.tsv[i,1],sep = "\t",stringsAsFactors = F)
#simplify by genus
refseq.all$genus = unlist(lapply(strsplit(refseq.all[,3],split = " "), `[[`, 1))
refseq.all = refseq.all[,c(1,2,4)] %>% group_by(genus) %>% summarise_all(sum)
refseq.all = refseq.all[order(refseq.all$V1,decreasing =T),]
all_spp = c(all_spp,refseq.all[1:25,1]) #keep only top 10
all_spp_five = c(all_spp_five,refseq.all[refseq.all[,2]>5,1]) #keep only the ones that are at more than 5%

#all_spp = c(all_spp,refseq.all[1:10,3]) #keep only top 10
#all_spp_five = c(all_spp_five,refseq.all[refseq.all[,1]>5,3]) #keep only the ones that are at more than 5%
}
#all_spp_five_m = unique(sort(all_spp_five))
#all_spp_m = data.frame(unique(sort(all_spp)),stringsAsFactors = F)
#simplify by genus
all_spp_five_m = unique(sort(unlist(all_spp_five)))
all_spp_m = data.frame(unique(sort(unlist(all_spp))),stringsAsFactors = F)


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
  replicate = substring(strsplit(all.tsv[i,1],split = "WatPhotz_")[[1]][2],4,4)
  location = substring(strsplit(all.tsv[i,1],split = "Champ")[[1]][2],1,3)
  date = strsplit(all.tsv[i,1],split = "-")[[1]][2]
  samples = c(samples,rep(paste(date,replicate,sep = "_"),nrow(all_spp_m)))
  locations = c(locations,rep(location,nrow(all_spp_m)))
}

all_spp_m_ggplot = data.frame(unlist(all_spp_m[,2:66]))
temp = all_spp_m_ggplot
all_spp_m_ggplot[,2] = rep(all_spp_m[,1],65)
all_spp_m_ggplot[,3] = samples
all_spp_m_ggplot[,4] = locations
all_spp_m_ggplot[,5] = factor(all_spp_m_ggplot[,4], levels=c('St1','St2','PRM'))
colnames(all_spp_m_ggplot) = c("fraction","species","samples","locations","locations_f")

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
p1=ggplot() + labs(title = "Lake Champlain - annotated genera",fill = "Taxonomy") +
  theme_bw() + 
  #scale_x_date(date_breaks = "months" , date_labels = "%b") +
  theme(legend.text = element_text(face="italic"),plot.title = element_text(hjust = 0.5, size=14, face="bold")) + scale_fill_manual(values = x(length(all_spp_five_m))) +
  geom_bar(aes(y = fraction, x = samples, fill = species),
           data = all_spp_m_ggplot_top12,stat="identity") + ylab("fraction of annotated genera")  + theme(axis.text.x = element_text(angle = 90, hjust = 1))


#PDF (dimensions in inches)
dev.new(width=10, height=7,noRStudioGD = TRUE)
p1 + facet_grid(rows=vars(locations_f))
dev.print(device=pdf,"figures/Champlain_genera_barplot.pdf", onefile=FALSE)
dev.off()

#PNG
png("figures/Champlain_genera_barplot.png",width=10, res =400,height=7,units= 'in')
p1 + facet_grid(rows=vars(locations_f))
dev.off()


###
###Plot the Fold change over time...
###
#add dates to the dataframe to calculate mean per sampling date replicate per location
all_spp_m_ggplot_top12$dates = unlist(strsplit(all_spp_m_ggplot_top12$samples,"_"))[seq(1,(nrow(all_spp_m_ggplot_top12)*2),by = 2)]
all_spp_m_ggplot_top12$dates = ymd(all_spp_m_ggplot_top12$dates)
#Keep only fractions
fraction = as.data.frame(all_spp_m_ggplot_top12$fraction)

#dplyr fractions data per replicate, per location and per species
all_spp_m_ggplot_top12.replicates_merged = fraction %>% group_by(all_spp_m_ggplot_top12$locations,all_spp_m_ggplot_top12$dates,all_spp_m_ggplot_top12$species) %>% summarise_all(mean)

#add columns for FoldChanges.
all_spp_m_ggplot_top12.replicates_merged$fold_change_mean = 0
all_spp_m_ggplot_top12.replicates_merged$fold_change = 0

#add the locations as factors as well + give column names
all_spp_m_ggplot_top12.replicates_merged$locations_f = 0
colnames(all_spp_m_ggplot_top12.replicates_merged) = c("locations","dates","species","fraction","fold_change_mean","fold_change","locations_f")
all_spp_m_ggplot_top12.replicates_merged = as.data.frame(all_spp_m_ggplot_top12.replicates_merged)
all_spp_m_ggplot_top12.replicates_merged$locations_f = factor(all_spp_m_ggplot_top12.replicates_merged[,1], levels=c('St1','St2','PRM'))

####Fold change progression from Time 0 (logFC = 1)
for(i in 1: nrow(all_spp_m_ggplot_top12.replicates_merged))
{
  #per species
  temp = all_spp_m_ggplot_top12.replicates_merged[all_spp_m_ggplot_top12.replicates_merged$species == all_spp_m_ggplot_top12.replicates_merged[i,3],]
  
  #per location
  temp = as.data.frame(temp)
  temp2 = temp[temp$locations == all_spp_m_ggplot_top12.replicates_merged[i,1],]
  
  #value divided by mean value
  all_spp_m_ggplot_top12.replicates_merged$fold_change_mean[i] = log(all_spp_m_ggplot_top12.replicates_merged[i,4] / mean(temp2[,4], na.rm =T),10)
  #value divided by Time 0
  all_spp_m_ggplot_top12.replicates_merged$fold_change[i] = log(all_spp_m_ggplot_top12.replicates_merged[i,4] / temp2[temp2[,2]=="2016-06-01",4],10)
    }

#add cyano column 
all_spp_m_ggplot_top12.replicates_merged$cyano = "0"
all_spp_m_ggplot_top12.replicates_merged$cyano[all_spp_m_ggplot_top12.replicates_merged[,3] == "Anabaena"] = "cyano"
all_spp_m_ggplot_top12.replicates_merged$cyano[all_spp_m_ggplot_top12.replicates_merged[,3] == "Dolichospermum"] = "cyano"
all_spp_m_ggplot_top12.replicates_merged$cyano[all_spp_m_ggplot_top12.replicates_merged[,3] == "Microcystis"] = "cyano"
all_spp_m_ggplot_top12.replicates_merged$cyano[all_spp_m_ggplot_top12.replicates_merged[,3] == "Synechococcus"] = "cyano"
all_spp_m_ggplot_top12.replicates_merged$cyano[all_spp_m_ggplot_top12.replicates_merged[,3] == "Nostocales"] = "cyano"

#ggplot object
x = colorRampPalette(brewer.pal(12,"Paired"))
p1=ggplot() + labs(title = "Lake Champlain - annotated genera",fill = "Taxonomy") +
  theme_bw() + 
  theme(legend.text = element_text(face="italic"),plot.title = element_text(hjust = 0.5, size=14, face="bold")) + scale_colour_manual(aesthetics = "colour",values = x(length(all_spp_five_m))) +
  geom_line(aes(y = fold_change, x = dates, linetype = cyano, colour = species,group=species),size = 2,data = all_spp_m_ggplot_top12.replicates_merged,stat="identity") +
  ylab("log10 Fold change since Time 0")  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(colour="Genera (dashed = cyanobacteria)") +
  scale_linetype_discrete(name="cyano",guide=F) +
  scale_x_date(date_breaks = "months" , date_labels = "%b")

#PDF (dimensions in inches)
dev.new(width=10, height=7,noRStudioGD = TRUE)
p1 + facet_grid(rows=vars(locations_f))
dev.print(device=pdf,"figures/Champlain_genera_foldchange.pdf", onefile=FALSE)
dev.off()

##PNG
png("figures/Champlain_genera_foldchange.png",width=10, res =400,height=7,units= 'in')
p1 + facet_grid(rows=vars(locations_f))
dev.off()



#save
write.table(all_spp_m_ggplot_top12.replicates_merged,"results/all_spp_m_ggplot_top12_metatranscriptome_genera",row.names = F, col.names = T, quote = F)


