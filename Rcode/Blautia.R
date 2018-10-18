
###Blautia refseq hits -----
setwd("/home/renaut/scratch/blast_database/Blautia")
#check paths...
system("grep 'Blautia' ../refseq/RefSeq_bac.fa | grep 'massiliensis' >Blautia_refseq")

#blautia names
blautia = read.table("Blautia_refseq",stringsAsFactors = F,sep = '\t',na.strings="noneNAs",quote = "")
blautia_names = rep(0,nrow(blautia))

#simplify names
for(i in 1:nrow(blautia))
{
blautia_names[i] = strsplit(blautia[i,1]," ")[[1]][1]
blautia_names[i] = sub(">","",blautia_names[i],fixed = T)
}

#Diamond annotation against refseq
annot = read.table("/home/renaut/scratch/ATRAPP_Champlain_2016_Metat/samsa2_run_St1/step_4_annotation/HI.4752.003.NEBNext_Index_13.000576_ChampSt1-20160915-WatPhotz_RNAa_merged.RefSeq_annotated", stringsAsFactors = F)
annot = cbind(annot,0)

#which ones are Blautia?
for(i in 1:nrow(annot))
{
temp = blautia_names %in% annot[i,2]
if(length(temp[temp==T])==1) annot[i,13] = blautia_names[temp==T]
if(i%%50000==0) print(paste("Done match ",i, " of",nrow(annot)," ;The time is: ",Sys.time(),sep = ""))
}

#A single hit...
rle(sort(annot[,13]))

###Blautia genome hits ----

setwd("/Users/jerry/Dropbox/CSBQ/shapiro")

#blast out
blastn.out = read.table("results/ChampSt1-20160915-WatPhotz_RNAa.fasta100k.blastn.out",stringsAsFactors = F)

hist(blastn.out[,9],breaks = 1000)

dev.new(noRStudioGD = TRUE)
par(mfrow = c(2,2))

#region1
region.1.1t = blastn.out[blastn.out[,9]>1500000,9]
region.1.2t = blastn.out[blastn.out[,9]>1500000,10]
region.1.1 = region.1.1t[region.1.1t<1510000]
region.1.1.o = region.1.1[order(region.1.1)]/1000000
region.1.2 = region.1.2t[region.1.1t<1510000]
region.1.2.o = region.1.2[order(region.1.1)]/1000000

plot(region.1.1.o,1:length(region.1.1),col = "darkred",main = "B. massiliensis Region 1",xlab = "genomic position (Mb)",yaxt ="n",ylab="",pch =20,cex =0.8)
points(region.1.2.o,1:length(region.1.1),col = "darkblue",pch =20,cex =0.8)
segments(region.1.1.o,1:length(region.1.1),region.1.2.o,1:length(region.1.1),col="darkgrey")

#region2
region.1.1t = blastn.out[blastn.out[,9]>2300000,9]
region.1.2t = blastn.out[blastn.out[,9]>2300000,10]
region.1.1 = region.1.1t[region.1.1t<2365000]
region.1.1.o = region.1.1[order(region.1.1)]/1000000
region.1.2 = region.1.2t[region.1.1t<2365000]
region.1.2.o = region.1.2[order(region.1.1)]/1000000

plot(region.1.1.o,1:length(region.1.1),col = "darkred",main = "B. massiliensis Region 2",xlab = "genomic position (Mb)",yaxt ="n",ylab="",pch =20,cex =0.8)
points(region.1.2.o,1:length(region.1.1),col = "darkblue",pch =20,cex =0.8)
segments(region.1.1.o,1:length(region.1.1),region.1.2.o,1:length(region.1.1),col="darkgrey")

#region2 & 3
region.1.1t = blastn.out[blastn.out[,9]>2380000,9]
region.1.2t = blastn.out[blastn.out[,9]>2380000,10]
region.1.1 = region.1.1t[region.1.1t<2400000]
region.1.1.o = region.1.1[order(region.1.1)]/1000000
region.1.2 = region.1.2t[region.1.1t<2400000]
region.1.2.o = region.1.2[order(region.1.1)]/1000000

plot(region.1.1.o,1:length(region.1.1),col = "darkred",main = "B. massiliensis Region 3",xlab = "genomic position (Mb)",yaxt ="n",ylab="",pch =20,cex =0.8)
points(region.1.2.o,1:length(region.1.1),col = "darkblue",pch =20,cex =0.8)
segments(region.1.1.o,1:length(region.1.1),region.1.2.o,1:length(region.1.1),col="darkgrey")

#region4 2632560
region.1.1t = blastn.out[blastn.out[,9]>2600000,9]
region.1.2t = blastn.out[blastn.out[,9]>2600000,10]
region.1.1 = region.1.1t[region.1.1t<2640000]
region.1.1.o = region.1.1[order(region.1.1)]/1000000
region.1.2 = region.1.2t[region.1.1t<2640000]
region.1.2.o = region.1.2[order(region.1.1)]/1000000

plot(region.1.1.o,1:length(region.1.1),col = "darkred",main = "B. massiliensis Region 4",xlab = "genomic position (Mb)",yaxt ="n",ylab="",pch =20,cex =0.8)
points(region.1.2.o,1:length(region.1.1),col = "darkblue",pch =20,cex =0.8)
segments(region.1.1.o,1:length(region.1.1),region.1.2.o,1:length(region.1.1),col="darkgrey")

dev.print(device=pdf, "figures/Bmassiliensis_genome_mapping.pdf", onefile=FALSE)
dev.off()
