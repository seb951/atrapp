#metagenome vs metatranscriptome

setwd("/Users/jerry/Documents/CSBQ/shapiro")

library(lubridate)
library(ggplot2)
library(RColorBrewer)

metat_genera = read.table("results/all_spp_m_ggplot_top12_metatranscriptome_genera",stringsAsFactors = F,header = T)
short_genome = read.table("results/all_spp_m_ggplot_top12_short_genome",stringsAsFactors = F, header = T)


metat_genera$dates = ymd(metat_genera$dates)
metat_genera_st1 = metat_genera[metat_genera$locations_f=="St1",]
colnames(metat_genera_st1)[7] = "locations"
metat_genera_st1[,7] = "Littoral (transcriptome-2016)"
short_genome$dates = ymd(short_genome$samples)

y2015 = interval(ymd("2015-01-01"), ymd("2016-01-01"))
y2016 = interval(ymd("2016-01-01"), ymd("2017-01-01"))

#two years...
short_genome$locations[short_genome$dates %within% y2015] = "Littoral (genome-2015)"
short_genome$locations[short_genome$dates %within% y2016] = "Littoral (genome-2016)"
short_genome$dates = gsub("2015","2016",short_genome$dates)
short_genome$dates = ymd(short_genome$dates)

metat_metag = rbind(short_genome[,c(1,2,9,4)],metat_genera_st1[,c(4,3,2,7)])
metat_metag$dates = ymd(metat_metag$dates)
metat_metag$dates_c = as.character(metat_metag$dates)

#ggplot object
x = colorRampPalette(brewer.pal(12,"Paired"))
p1=ggplot() + labs(title = "Lake Champlain - metatranscriptomes & metagenomes",fill = "Taxonomy (genera)") +
  theme_bw() + 
  #  scale_x_date(date_breaks = "weeks") +
  facet_grid(rows=vars(locations),scales="free") +
  theme(legend.text = element_text(face="italic"),axis.title = element_text(size = 20),plot.title = element_text(hjust = 0.5, size=20, face="bold")) + scale_fill_manual(values = x(length(unique(sort(metat_metag[,2]))))) +
  geom_bar(aes(y = fraction, x = dates_c, fill = species),
           data = metat_metag,stat="identity") + xlab("Sampling dates") + ylab("Fraction of annotated reads")  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


#ggplot
dev.new(width=10, height=7,noRStudioGD = TRUE)
p1
dev.print(device=pdf,"figures/metagenome_vs_metatranscriptome_pretty.pdf", onefile=FALSE)
dev.off()


#####
system("ls -1 results/cyano/func_results/*annot_function.tsv >all.tsv")
all.tsv = read.table("all.tsv", stringsAsFactors = F)
all_spp = NULL
samples = NULL
locations = NULL
all_spp_five = NULL
#first loop is to subset the major species (over 1%)
for(i in 1:nrow(all.tsv))
{
  refseq.all = read.table(all.tsv[i,1],sep = "\t",stringsAsFactors = F)
  
  #unknowns
  refseq.all[refseq.all[,3] == "Not",3] = c("unknown unknown unknown unknown unknown unknown")
  
  #simplify name
  refseq.all$genes1 = unlist(lapply(strsplit(refseq.all[,3],split = " "), `[[`, 1))
  refseq.all$genes2 = unlist(lapply(strsplit(refseq.all[,3],split = " "), `[[`, 2))
  refseq.all$genes3 = unlist(lapply(strsplit(refseq.all[,3],split = " "), `[[`, 3))
  refseq.all$genes = paste(refseq.all$genes1,refseq.all$genes2,refseq.all$genes3,sep = " ")
  refseq.all = refseq.all[,c(1,2,7)] %>% group_by(genes) %>% summarise_all(sum)
  refseq.all = refseq.all[order(refseq.all$V1,decreasing =T),]
  
  all_spp = c(all_spp,refseq.all[,1]) #keep only top 10
  all_spp_five = c(all_spp_five,refseq.all[refseq.all[,3]>150,1]) #keep only the ones that are at more than a count of 1....
}
#simplify by genus
all_spp_five_m = unique(sort(unlist(all_spp_five)))
all_spp_m = data.frame(unique(sort(unlist(all_spp))),stringsAsFactors = F)



for(i in 1:nrow(all.tsv))
{
  refseq.all = read.table(all.tsv[i,1],sep = "\t",stringsAsFactors = F)
  
  #unknowns
  refseq.all[refseq.all[,3] == "Not",3] = c("unknown unknown unknown unknown")
  
  #simplify by genus
  refseq.all$genes1 = unlist(lapply(strsplit(refseq.all[,3],split = " "), `[[`, 1))
  refseq.all$genes2 = unlist(lapply(strsplit(refseq.all[,3],split = " "), `[[`, 2))
  refseq.all$genes3 = unlist(lapply(strsplit(refseq.all[,3],split = " "), `[[`, 3))
  refseq.all$genes = paste(refseq.all$genes1,refseq.all$genes2,refseq.all$genes3,sep = " ")
  refseq.all = refseq.all[,c(1,2,7)] %>% group_by(genes) %>% summarise_all(sum)
  
  for(j in 1:nrow(all_spp_m))
  {
    temp = refseq.all[refseq.all[,1] == all_spp_m[j,1],3]
    if(length(temp$V2) == 1) all_spp_m[j,i+1] = temp
    if(length(temp$V2) == 0) all_spp_m[j,i+1] = 0
  }
  #get a proper shortname for the graph
  replicate = substring(strsplit(all.tsv[i,1],split = "WatPhotz_")[[1]][2],4,4)
  location = substring(strsplit(all.tsv[i,1],split = "Champ")[[1]][2],1,3)
  date = strsplit(all.tsv[i,1],split = "-")[[1]][2]
  samples = c(samples,rep(paste(date,replicate,sep = "_"),nrow(all_spp_m)))
  locations = c(locations,rep(location,nrow(all_spp_m)))
}

all_spp_m_ggplot = data.frame(unlist(all_spp_m[,2:66]))
all_spp_m_ggplot[,2] = rep(all_spp_m[,1],65)
all_spp_m_ggplot[,3] = samples
all_spp_m_ggplot[,4] = locations
all_spp_m_ggplot[,5] = factor(all_spp_m_ggplot[,4], levels=c('St1','St2','PRM'))
colnames(all_spp_m_ggplot) = c("fraction","species","samples","locations","locations_f")
all_spp_m_ggplot[,5] = factor(all_spp_m_ggplot$locations, levels=c('St1','St2','PRM'))
all_spp_m_ggplot[,6] = unlist(lapply(strsplit(all_spp_m_ggplot$samples,split = "_"), `[[`, 1))
all_spp_m_ggplot[,6] = ymd(all_spp_m_ggplot[,6])
colnames(all_spp_m_ggplot)[6] = "dates"


###incorporate the transcriptome datastructure
metag_false = data.frame(fraction=0,species="alpha alpha alpha",samples=unique(sort(metat_metag$dates)),locations="transcriptome",locations_f="transcriptome",dates=unique(sort(metat_metag$dates)))
colnames(metag_false) = colnames(all_spp_m_ggplot)

all_spp_m_ggplot= rbind(all_spp_m_ggplot,metag_false)

all_spp_m_ggplot[,7] = as.character(all_spp_m_ggplot$dates)
colnames(all_spp_m_ggplot)[7] = "dates_c"

#top percentage according to "all_spp_five_m" vector
all_spp_five_m = c(all_spp_five_m,"alpha alpha alpha")
all_spp_m_ggplot_top12 = NULL
for(i in 1:length(all_spp_five_m))
{
  all_spp_m_ggplot_top12 = rbind(all_spp_m_ggplot_top12,all_spp_m_ggplot[all_spp_m_ggplot[,2] == all_spp_five_m[i],])
}

#just get a better name
temp1 = unlist(lapply(strsplit(all_spp_m_ggplot_top12[,2],split = " "), `[[`, 1))
temp2 = unlist(lapply(strsplit(all_spp_m_ggplot_top12[,2],split = " "), `[[`, 2))
temp3 = unlist(lapply(strsplit(all_spp_m_ggplot_top12[,2],split = " "), `[[`, 3))
all_spp_m_ggplot_top12$names = paste(temp3, " (",temp1," ",temp2,")",sep = "")
all_spp_m_ggplot_top12$names = temp3

#plot
x = colorRampPalette(brewer.pal(9,"Set1"))
p_cyano=ggplot() + labs(title = "Lake Champlain - cyanotoxic genes",fill = "Gene names (species)") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) + scale_fill_manual(values = x(length(all_spp_five_m))) +
  geom_bar(aes(y = fraction, x = dates_c, fill = names),
           data = all_spp_m_ggplot_top12,stat="identity") + ylab("Counts of annotated cyanotoxic genes")  + theme(axis.text.x = element_text(angle = 90, hjust = 1))


#PDF (dimensions in inches)
dev.new(width=10, height=5.3333,noRStudioGD = TRUE)
p_cyano
dev.print(device=pdf,"figures/Champlain_cyano_barplot_pretty.pdf", onefile=FALSE)
dev.off()
