

setwd("/home/renaut/scratch/blast_database/Blautia")

grep 'Blautia' ../refseq/RefSeq_bac.fa | grep 'massiliensis' >Blautia_refseq


#blautia names
blautia = read.table("Blautia_refseq",stringsAsFactors = F,sep = '\t',na.strings="noneNAs",quote = "")
blautia_names = rep(0,nrow(blautia))

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
