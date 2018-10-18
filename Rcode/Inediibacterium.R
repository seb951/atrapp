
###Inediibacterium genome hits ----

setwd("/Users/jerry/Dropbox/CSBQ/shapiro")

#blast out
blastn.out = read.table("results/Inediibacterium_ChampSt1-20160915-WatPhotz_RNAa.fasta100k.blastn.out",stringsAsFactors = F)

hist(blastn.out[,9],breaks = 1000)

dev.new(noRStudioGD = TRUE)
#par(mfrow = c(2,2))

#region1
region.1.1t = blastn.out[blastn.out[,9]>1637200,9]
region.1.2t = blastn.out[blastn.out[,9]>1637200,10]
region.1.1 = region.1.1t[region.1.1t<1638000]
region.1.1.o = region.1.1[order(region.1.2)]/1000
region.1.2 = region.1.2t[region.1.1t<1638000]
region.1.2.o = region.1.2[order(region.1.2)]/1000

plot(region.1.2.o,1:length(region.1.1),col = "darkblue",main = "I. massiliense",xlab = "genomic position",yaxt ="n",ylab="",pch =20,cex =0.8)
points(region.1.1.o,1:length(region.1.1),col = "darkred",pch =20,cex =0.8)
segments(region.1.1.o,1:length(region.1.1),region.1.2.o,1:length(region.1.1),col="darkgrey")

dev.print(device=pdf, "figures/Imassiliense_genome_mapping.pdf", onefile=FALSE)
dev.off()
