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
metat_genera_st1[,7] = "St1 (transcriptome-2016)"
short_genome$dates = ymd(short_genome$samples)

y2015 = interval(ymd("2015-01-01"), ymd("2016-01-01"))
y2016 = interval(ymd("2016-01-01"), ymd("2017-01-01"))

#two years...
short_genome$locations[short_genome$dates %within% y2015] = "St1 (genome-2015)"
short_genome$locations[short_genome$dates %within% y2016] = "St1 (genome-2016)"
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
  facet_grid(rows=vars(locations)) +
  #scales="free"
  theme(legend.text = element_text(face="italic"),axis.title = element_text(size = 20),plot.title = element_text(hjust = 0.5, size=20, face="bold")) + scale_fill_manual(values = x(length(unique(sort(metat_metag[,2]))))) +
  geom_bar(aes(y = fraction, x = dates, fill = species),
           data = metat_metag,stat="identity") + xlab("Sampling dates") + ylab("Fraction of annotated reads")  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


#ggplot
dev.new(width=10, height=7,noRStudioGD = TRUE)
p1
dev.print(device=pdf,"figures/metagenome_vs_metatranscriptome.pdf", onefile=FALSE)
dev.off()

png("figures/metagenome_vs_metatranscriptome.png",width=10, res =400,height=7,units= 'in')
p1
dev.off()


####fold changes for cyano metag and metat.
#remove 2015
metat_metag = metat_metag[metat_metag$locations != "St1 (genome-2015)",]

#Keep only fractions
fraction = as.data.frame(metat_metag$fraction)

#add columns for FoldChanges.
metat_metag$fold_change_mean = 0
metat_metag$fold_change = 0

####Fold change progression from Time 0 (logFC = 1)
for(i in 1: nrow(metat_metag))
{
 #   if(metat_metag$locations[i] == "St1 (genome-2015)") loc_specific_time_0 = "2016-07-16"
    if(metat_metag$locations[i] == "St1 (genome-2016)") loc_specific_time_0 = "2016-07-21"
    if(metat_metag$locations[i] == "St1 (transcriptome-2016)") loc_specific_time_0 = "2016-06-01"
  
  #per species
  temp = metat_metag[metat_metag$species == metat_metag$species[i],]
  
  #per location
  temp = as.data.frame(temp)
  temp2 = temp[temp$locations == metat_metag$locations[i],]
  
  #value divided by mean value
  if(nrow(temp2)>0) metat_metag$fold_change_mean[i] = log(metat_metag$fraction[i] / mean(temp2$fraction, na.rm =T),10)
  #value divided by Time 0
  if((nrow(temp2)>0) & (length(temp2$fraction[temp2$dates==loc_specific_time_0]) >0)) metat_metag$fold_change[i] = log(metat_metag$fraction[i] / temp2$fraction[temp2$dates==loc_specific_time_0],10)
}

#Keep only cyano
metat_metag$cyano = "0"
metat_metag$cyano[metat_metag[,2] == "Anabaena"] = "cyano"
metat_metag$cyano[metat_metag[,2] == "Dolichospermum"] = "cyano"
metat_metag$cyano[metat_metag[,2] == "Microcystis"] = "cyano"
metat_metag$cyano[metat_metag[,2] == "Synechococcus"] = "cyano"
metat_metag$cyano[metat_metag[,2] == "Nostocales"] = "cyano"
#metat_metag = metat_metag[metat_metag$cyano != "0",]

#ggplot object
x = colorRampPalette(brewer.pal(12,"Paired"))
p1=ggplot() + labs(title = "Lake Champlain - annotated cyanobacteria",fill = "Taxonomy") +
  theme_bw() + 
  theme(legend.text = element_text(face="italic"),plot.title = element_text(hjust = 0.5, size=14, face="bold")) + scale_colour_manual(aesthetics = "colour",values = x(length(unique(sort(metat_metag$species))))) +
  geom_line(aes(y = fold_change, x = dates, linetype = cyano, colour = species,group=species),size = 2,data = metat_metag,stat="identity") +
  ylab("log10 Fold change since Time 0")  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(colour="cyanobacteria") +
  #scale_linetype_discrete(name="cyano",guide=F) +
  scale_x_date(date_breaks = "weeks")

#PDF (dimensions in inches)
dev.new(width=10, height=7,noRStudioGD = TRUE)
p1 + facet_grid(rows=vars(locations))
dev.print(device=pdf,"figures/Champlain_genera_foldchange_metag_metat.pdf", onefile=FALSE)
dev.off()




####Mantel test
#To compare both the metag and metat data matrices.
library(ape)
library(vegan)

#mantel.test()

#PERMANOVA
#see how well metag-2016 explains the metat-2015... 
#must keep only the samples with metag and metat data...
metat_metag_2016 = metat_metag[metat_metag[,4]!="St1 (genome-2015)",]

metat_metag_2016_permanova = metat_metag_2016[metat_metag_2016[,3]=="2016-07-21",]
metat_metag_2016_permanova = rbind(metat_metag_2016_permanova,metat_metag_2016[metat_metag_2016[,3]=="2016-08-02",])
metat_metag_2016_permanova = rbind(metat_metag_2016_permanova,metat_metag_2016[metat_metag_2016[,3]=="2016-08-08",])
metat_metag_2016_permanova = rbind(metat_metag_2016_permanova,metat_metag_2016[metat_metag_2016[,3]=="2016-08-09",])
metat_metag_2016_permanova = rbind(metat_metag_2016_permanova,metat_metag_2016[metat_metag_2016[,3]=="2016-09-01",])
metat_metag_2016_permanova = rbind(metat_metag_2016_permanova,metat_metag_2016[metat_metag_2016[,3]=="2016-09-19",])
metat_metag_2016_permanova = rbind(metat_metag_2016_permanova,metat_metag_2016[metat_metag_2016[,3]=="2016-10-06",])

metat_metag_2016_permanova[metat_metag_2016_permanova[,3]=="2016-08-08",3] = "2016-08-09"

metag_2016_permanova = metat_metag_2016_permanova[metat_metag_2016_permanova[,4]=="St1 (genome-2016)",]
metat_2016_permanova = metat_metag_2016_permanova[metat_metag_2016_permanova[,4]=="St1 (transcriptome-2016)",]


#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
asv.filt.abundants.norm.hel <-decostand(asv.filt.abundants.norm, "hel")
permanova=adonis(formula=asv.filt.abundants.norm.hel~fertilization*species,strata=(design.keep$bloc/design.keep$replicate), data=design.keep, permutations=9999, method="bray")



