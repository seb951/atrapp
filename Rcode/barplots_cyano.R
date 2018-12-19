
setwd("/Users/jerry/Documents/CSBQ/shapiro")

#packages
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(dplyr)

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




#top percentage according to "all_spp_five_m" vector
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

#plot
x = colorRampPalette(brewer.pal(12,"Paired"))
p1=ggplot() + labs(title = "Lake Champlain - cyanotoxic genes",fill = "Gene names (species)") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) + scale_fill_manual(values = x(length(all_spp_five_m))) +
  geom_bar(aes(y = fraction, x = samples, fill = names),
           data = all_spp_m_ggplot_top12,stat="identity") + ylab("Counts of annotated cyanotoxic genes")  + theme(axis.text.x = element_text(angle = 90, hjust = 1))


#PDF (dimensions in inches)
dev.new(width=10, height=8,noRStudioGD = TRUE)
p1 + facet_grid(rows=vars(locations_f))
dev.print(device=pdf,"figures/Champlain_cyano_barplot.pdf", onefile=FALSE)
dev.off()

#PNG
png("figures/Champlain_cyano_barplot.png",width=10, res =400,height=8,units= 'in')
p1 + facet_grid(rows=vars(locations_f))
dev.off()


  
  
