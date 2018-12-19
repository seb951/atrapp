

#get the names in proper format for SAMSA2 DIAMOND annotation
fasta = read.table("~/Desktop/New_Toxin_Database_OneLine_1.faa",stringsAsFactors = F)
annot = read.table("~/Desktop/New_Toxin_Database_OneLine_1.fxn",sep="\t",stringsAsFactors = F)

fasta_names = fasta[seq(1,nrow(fasta),by = 2),1]
functions = rep("n",nrow(annot))
species  = rep("n",nrow(annot))
for(i in 1:nrow(annot))
  {
  temp = strsplit(annot[i,3],split = " ")[[1]]
  functions[i] = paste(temp[-c(length(temp))],collapse = " ")
  species[i] = paste(temp[c(1,2)],collapse = " ")
  
  fasta_names[i] = paste(fasta_names[i]," ",functions[i]," [",species[i],"]",collapse = " ",sep = "")
  }

fasta[seq(1,nrow(fasta),by = 2),1] = fasta_names

write.table(fasta,"~/Desktop/New_Toxin_Database_OneLine_1.faa2",quote = F, row.names =F, col.names = F)


