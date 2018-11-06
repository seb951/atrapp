setwd("/Users/jerry/Documents/CSBQ/shapiro")

##############################-----
#########Fixing the missing data
##############################

system("ls -1 results/Subsystems_results/*tated.hierarchy.reduced >all.tsv")
all.tsv = read.table("all.tsv", stringsAsFactors = F)
#first loop is to subset the major species (over 1%)
for(i in 1:nrow(all.tsv))
{
refseq.all.temp = read.table(all.tsv[i,1],stringsAsFactors = F,sep = "#",quote = "")
#
#this loop is to fix the NO HIERARCHY and "" problems...
refseq.all = data.frame(rep(0,nrow(refseq.all.temp)))
refseq.all$percent = 0
refseq.all$count = 0
refseq.all$hier1 = 0
refseq.all$hier2 = 0
refseq.all$hier3 = 0
refseq.all$hier4 = 0

  for(n in 1:nrow(refseq.all.temp))
  {
    temp = strsplit(refseq.all.temp[n,],"\t")
    refseq.all$percent[n] = temp[[1]][1]
    refseq.all$count[n] = temp[[1]][2]
    refseq.all$hier1[n] = temp[[1]][3]
    refseq.all$hier2[n] = temp[[1]][4]
    refseq.all$hier3[n] = temp[[1]][5]
    refseq.all$hier4[n] = temp[[1]][6]
    refseq.all[n,is.na(refseq.all[n,])] = "NO HIERARCHY" #no hiercharchy (nothing at ll)
    refseq.all[n,nchar(refseq.all[n,])==0] = refseq.all[n,7] #some files have nothing quoted as ""
  }

write.table(refseq.all[,2:7],paste(all.tsv[i,1],"_modif",sep = ""),quote = T, row.names = F, col.names =T)

if(i%%5==0) print(paste("file nb: ",i," of:",nrow(all.tsv),"; The time is: ",Sys.time(),sep = ""))
}


#######################################
###hierachical level info -----
###get all the data
#######################################
all.tsv = read.table("all.tsv", stringsAsFactors = F)
#packages
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
level=3 #choose levels 1,2,3,4

samples = NULL
locations = NULL
all_spp = NULL
refseq.all.list = as.list(1)
for(i in 1:nrow(all.tsv))
{
  refseq.all = read.table(paste(all.tsv[i,1],"_modif",sep = ""),header = T,stringsAsFactors = F)
  refseq.all.list[[i]] = refseq.all
  #get a proper shortname for the graph
  replicate = substring(strsplit(all.tsv[i,1],split = "WatPhotz_")[[1]][2],4,4)
  location = substring(strsplit(all.tsv[i,1],split = "Champ")[[1]][2],1,3)
  date = strsplit(all.tsv[i,1],split = "-")[[1]][2]
  samples = c(samples,paste(date,replicate,sep = "_"))
  locations = c(locations,location)
  all_spp = c(all_spp,refseq.all[,level+2])
}

#unique.
all_spp_m = unique(sort(all_spp))

#level 4 uniques
refseq.all.uniques = data.frame(all_spp_m,stringsAsFactors = F)

for(j in 1:65)
  {
    for(i in 1:length(all_spp_m))
    {
      temp = refseq.all.list[[j]][refseq.all.list[[j]][,level+2] == all_spp_m[i],]
  #summing the percentages per categories for all categories...
      refseq.all.uniques[i,j+1] = sum(temp[,1])
    }
  }
colnames(refseq.all.uniques)[-1] = samples


#top X percent of functions
if(level == 1) X=0.5
if(level == 2|3|4) X=1
refseq.all.order = refseq.all.uniques[order(rowMeans(refseq.all.uniques[,-1]),decreasing=T),]
refseq.all.top1 = refseq.all.uniques[rowMeans(refseq.all.uniques[,-1])>X,]
dim(refseq.all.top1)

#ggplot object
all_spp_m_ggplot = data.frame(unlist(refseq.all.top1[,2:66]))
temp = all_spp_m_ggplot
all_spp_m_ggplot[,2] = rep(refseq.all.top1[,1],65)
all_spp_m_ggplot[,3] = as.vector(t(matrix(rep(samples,nrow(refseq.all.top1)),nrow=length(samples),ncol=nrow(refseq.all.top1))))
all_spp_m_ggplot[,4] = as.vector(t(matrix(rep(locations,nrow(refseq.all.top1)),nrow=length(locations),ncol=nrow(refseq.all.top1))))
all_spp_m_ggplot[,5] = factor(all_spp_m_ggplot[,4], levels=c('St1','St2','PRM'))
colnames(all_spp_m_ggplot) = c("fraction","species","samples","locations","locations_f")

#these are percentage values
all_spp_m_ggplot[,1] = all_spp_m_ggplot[,1]/100

#plot
x = colorRampPalette(brewer.pal(12,"Paired"))
p1=ggplot() + labs(title = "Lake Champlain - annotated functions",fill = "Taxonomy") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) + scale_fill_manual(values = x(nrow(refseq.all.top1))) +
  geom_bar(aes(y = fraction, x = samples, fill = species),
           data = all_spp_m_ggplot,stat="identity") + ylab("fraction of annotated species")  + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#PDF (dimensions in inches)
dev.new(width=10, height=7,noRStudioGD = TRUE)
p1 + facet_grid(rows=vars(locations_f))
dev.print(device=pdf,paste("figures/Champlain_functions_barplot_level",level,".pdf",sep =""), onefile=FALSE)
dev.off()

#PNG
png(paste("figures/Champlain_functions_barplot_level",level,".png",sep =""),width=10, res =400,height=7,units= 'in')
p1 + facet_grid(rows=vars(locations_f))
dev.off()



