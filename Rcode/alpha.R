#!/usr/local/bin/Rscript

setwd("/Users/jerry/Documents/CSBQ/shapiro")

#packages
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(vegan)
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
  all_spp = c(all_spp,refseq.all[1:10,3]) #keep only top 10
  all_spp_five = c(all_spp_five,refseq.all[refseq.all[,1]>3,3]) #keep only the ones that are at more than 5%
}
all_spp_five_m = unique(sort(all_spp_five))
all_spp_m = data.frame(unique(sort(all_spp)),stringsAsFactors = F)

alpha_div = data.frame(replicate = rep(0,65),location = 0, date = 0,samples = 0,alpha = 0)
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
  alpha_div[i,1] = substring(strsplit(all.tsv[i,1],split = "WatPhotz_")[[1]][2],4,4) #replicate
  alpha_div[i,2] = substring(strsplit(all.tsv[i,1],split = "Champ")[[1]][2],1,3) #location
  alpha_div[i,3] = strsplit(all.tsv[i,1],split = "-")[[1]][2] #date
  alpha_div[i,4] = paste(date,replicate,sep = "_") #samples
#  locations = c(locations,rep(location,nrow(all_spp_m)))
}
alpha_div[,2] = factor(alpha_div[,2],levels=c('St1','St2','PRM')) 
alpha_div[,3] = ymd(alpha_div[,3])
#alpha diversity
alpha_div[,5] = vegan::diversity(t(all_spp_m[,2:66]), index= "invsimpson")


#plot
p1 = ggplot() + labs(title = "Lake Champlain - alpha diversity") +
  geom_point(aes(x = date, y = alpha),na.rm=T,data = alpha_div) + 
  theme_bw() + 
  ylab("alpha diversity (Inverse Simpson Index)") +
  xlab("Sampling date") +
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  scale_x_date(date_breaks = "months" , date_labels = "%b") +
  facet_grid(rows=vars(location),scales="free")

#PDF (dimensions in inches)
dev.new(width=10, height=7,noRStudioGD = TRUE)
p1
dev.print(device=pdf,"figures/alpha.pdf", onefile=FALSE)
dev.off()

#PNG
png("figures/alpha.png",width=10, res =400,height=7,units= 'in')
p1
dev.off()


###
###Plot the Fold change over time...
###
#add dates to the dataframe to calculate mean per sampling date replicate per location
all_spp_m_ggplot_top12$dates = unlist(strsplit(all_spp_m_ggplot_top12$samples,"_"))[seq(1,(nrow(all_spp_m_ggplot_top12)*2),by = 2)]

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
  all_spp_m_ggplot_top12.replicates_merged$fold_change[i] = log(all_spp_m_ggplot_top12.replicates_merged[i,4] / temp2[temp2[,2]=="20160601",4],10)
}

#add cyano column 
all_spp_m_ggplot_top12.replicates_merged$cyano = "0"
all_spp_m_ggplot_top12.replicates_merged$cyano[all_spp_m_ggplot_top12.replicates_merged[,3] == "Anabaena sp."] = "cyano"
all_spp_m_ggplot_top12.replicates_merged$cyano[all_spp_m_ggplot_top12.replicates_merged[,3] == "Dolichospermum circinale"] = "cyano"
all_spp_m_ggplot_top12.replicates_merged$cyano[all_spp_m_ggplot_top12.replicates_merged[,3] == "Microcystis aeruginosa"] = "cyano"
all_spp_m_ggplot_top12.replicates_merged$cyano[all_spp_m_ggplot_top12.replicates_merged[,3] == "Synechococcus sp."] = "cyano"


#ggplot object
x = colorRampPalette(brewer.pal(12,"Paired"))
p1=ggplot() + labs(title = "Lake Champlain - annotated species",fill = "Taxonomy") +
  theme_bw() + 
  theme(legend.text = element_text(face="italic"),plot.title = element_text(hjust = 0.5, size=14, face="bold")) + scale_colour_manual(aesthetics = "colour",values = x(length(all_spp_five_m))) +
  geom_line(aes(y = fold_change, x = dates, linetype = cyano, colour = species,group=species),size = 2,data = all_spp_m_ggplot_top12.replicates_merged,stat="identity") +
  ylab("log10 Fold change since Time 0")  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(colour="Species (dashed = cyanobacteria)") +
  scale_linetype_discrete(name="cyano",guide=F)

#dev.new(width=10, height=7,noRStudioGD = TRUE)
p1 + facet_grid(rows=vars(locations_f),scales="free")

#PDF (dimensions in inches)
dev.new(width=10, height=7,noRStudioGD = TRUE)
p1 + facet_grid(rows=vars(locations_f))
dev.print(device=pdf,"figures/Champlain_species_foldchange.pdf", onefile=FALSE)
dev.off()

#PNG
png("figures/Champlain_species_foldchange.png",width=10, res =400,height=7,units= 'in')
p1 + facet_grid(rows=vars(locations_f))
dev.off()