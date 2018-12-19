#metagenome vs metatranscriptome


library(lubridate)

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

#ggplot object
x = colorRampPalette(brewer.pal(12,"Paired"))
p1=ggplot() + labs(title = "Lake Champlain - annotated genera",fill = "Taxonomy") +
  theme_bw() + 
  scale_x_date(date_breaks = "weeks") +
  facet_grid(rows=vars(locations),scales="free") +
  theme(legend.text = element_text(face="italic"),plot.title = element_text(hjust = 0.5, size=14, face="bold")) + scale_fill_manual(values = x(length(unique(sort(metat_metag[,2]))))) +
  geom_bar(aes(y = fraction, x = dates, fill = species),
           data = metat_metag,stat="identity") + ylab("fraction of annotated genera")  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


#ggplot
dev.new(width=10, height=7,noRStudioGD = TRUE)
p1
dev.print(device=pdf,"figures/metagenome_vs_metatranscriptome.pdf", onefile=FALSE)
dev.off()

png("figures/metagenome_vs_metatranscriptome.png",width=10, res =400,height=7,units= 'in')
p1
dev.off()
