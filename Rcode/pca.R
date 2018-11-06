#PCoA


setwd("/Users/jerry/Dropbox/CSBQ/shapiro")
 
#packages
#library(pals)
library(ggplot2)
#library(gridExtra)
#library(vegan)
#library(ape)
#library(lubridate)
library(phyloseq)

system("ls -1 results/org_results/*RefSeq_annot_organism.tsv >all.tsv")
all.tsv = read.table("all.tsv", stringsAsFactors = F)
all_spp = NULL
dates =  NULL
replicates = NULL
locations = NULL
full = NULL
#first loop is to subset the major species (over 1%)
for(i in 1:nrow(all.tsv))
{
  refseq.all = read.table(all.tsv[i,1],sep = "\t",stringsAsFactors = F)
  all_spp = c(all_spp,refseq.all[refseq.all[,1]>1,3]) #keep only top 1%
}
all_spp_m = data.frame(unique(sort(all_spp)),stringsAsFactors = F)

for(i in 1:nrow(all.tsv))
{
  refseq.all = read.table(all.tsv[i,1],sep = "\t",stringsAsFactors = F)
  for(j in 1:nrow(all_spp_m))
  {
    temp = refseq.all[refseq.all[,3] == all_spp_m[j,1],1]
    if(length(temp) == 1) all_spp_m[j,i+1] = temp
    if(length(temp) == 0) all_spp_m[j,i+1] = 0
  }
  replicate = substring(strsplit(all.tsv[i,1],split = "WatPhotz_")[[1]][2],4,4)
  location = substring(strsplit(all.tsv[i,1],split = "Champ")[[1]][2],1,3)
  date = strsplit(all.tsv[i,1],split = "-")[[1]][2]
  dates = c(dates,date)
  replicates = c(replicates,replicate)
  locations = c(locations,location)
  full = c(full,paste(location,date,replicate,sep = "_"))
  }

###get rid of samples st1 from sept15th, they are weid
all_spp_m = all_spp_m[,-c(6:8)]
dates = dates[-c(6:8)]
replicates = replicates[-c(6:8)]
locations = locations[-c(6:8)]
full = full[-c(6:8)]

#####phyloseq
####otu table
colnames(all_spp_m)[-1] = full
otumat = otu_table(all_spp_m[,-1],taxa_are_rows=T)
rownames(otumat) = all_spp_m[,1]

###FAKE tax table
taxmat = matrix(sample(letters,nrow(all_spp_m)*7, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat = tax_table(taxmat)

###design
design = sample_data(data.frame(dates,replicates,locations,full,row.names=colnames(all_spp_m)[-1],stringsAsFactors = F))
  

###phyloseq object
physeq = phyloseq(otumat, taxmat,design)


###ordination (NMDS is the same as pcoa?)

physeq.ord <- ordinate(physeq, "PCoA", "bray")
p1 = plot_ordination(physeq, physeq.ord, type="samples", color="dates", title="taxa",shape="locations") +
  theme_bw() +
  geom_point(stroke = 2) + 
  labs(title = "Lake Champlain - PcoA \n(36 species present in > 1% of annotated sequences)") +
  scale_shape(solid = T) +
  scale_color_brewer(palette= "PuOr")

#PDF
dev.new(width=6, height=5,noRStudioGD = TRUE)
p1
dev.print(device=pdf, "figures/Champlain_pcoa.pdf", onefile=FALSE)
dev.off()

#PNG
png("figures/Champlain_pcoa.png",width=6, res =400,height=5,units= 'in')
p1
dev.off()




