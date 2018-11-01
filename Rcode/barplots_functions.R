setwd("/Users/jerry/Dropbox/CSBQ/shapiro")

#packages
library(ggplot2)
library(gridExtra)

##############################-----
#########Fixing the missing data
##############################

system("ls -1 results/Subsystems_results/*tated.hierarchy.reduced >all.tsv")
all.tsv = read.table("all.tsv", stringsAsFactors = F)
all_spp = NULL
all_spp_shortnames = NULL
all_spp_one = NULL
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

all_spp = c(all_spp,refseq.all$hier1[1:30]) #keep only top 10
all_spp_one = c(all_spp_one,refseq.all$hier1[as.numeric(refseq.all$percent)>1]) #keep only the ones that are at more than 1% in at least one of the samples
}


#######################################
###hierachical level info -----
###get all the data
#######################################
level = 3 #choose levels 1,2,3,4
all_spp = NULL
all_spp_shortnames = NULL
all_spp_one = NULL
refseq.all.list = as.list(1)
for(i in 1:nrow(all.tsv))
{
  refseq.all = read.table(paste(all.tsv[i,1],"_modif",sep = ""),header = T,stringsAsFactors = F)
  refseq.all.list[[i]] = refseq.all
  #get a proper shortname for the graph
  replicate = substring(strsplit(all.tsv[i,1],split = "WatPhotz_")[[1]][2],4,4)
  location = substring(strsplit(all.tsv[i,1],split = "Champ")[[1]][2],1,3)
  date = strsplit(all.tsv[i,1],split = "-")[[1]][2]
  all_spp_shortnames = c(all_spp_shortnames,paste(date,location,replicate,sep = "_"))
  all_spp_one = c(all_spp_one,refseq.all[,level+2])
}

#unique.
all_spp_one_m = unique(sort(all_spp_one))

#level 4 uniques
refseq.all.uniques = data.frame(all_spp_one_m,stringsAsFactors = F)

for(j in 1:65)
  {
    for(i in 1:length(all_spp_one_m))
    {
      temp = refseq.all.list[[j]][refseq.all.list[[j]][,level+2] == all_spp_one_m[i],]
      refseq.all.uniques[i,j+1] = sum(temp[,1])
    }
  }
colnames(refseq.all.uniques)[-1] = all_spp_shortnames


#top X percent of functions
refseq.all.order = refseq.all.uniques[order(rowMeans(refseq.all.uniques[,-1]),decreasing=T),]
refseq.all.top1 = refseq.all.uniques[rowMeans(refseq.all.uniques[,-1])>1,]
dim(refseq.all.top1)

#ggplot object
all_spp_m_ggplot = data.frame(unlist(refseq.all.top1[,2:66]))
temp = all_spp_m_ggplot
all_spp_m_ggplot[,2] = rep(refseq.all.top1[,1],65)
all_spp_m_ggplot[,3] = as.vector(t(matrix(rep(all_spp_shortnames,nrow(refseq.all.top1)),nrow=length(all_spp_shortnames),ncol=nrow(refseq.all.top1))))
colnames(all_spp_m_ggplot) = c("fraction","species","sample")

#these are percentage values
all_spp_m_ggplot[,1] = all_spp_m_ggplot[,1]/100

#split into 3 locations.
all_spp_m_ggplot_st1 = all_spp_m_ggplot[regexpr("St1",all_spp_m_ggplot[,3])>0,]
all_spp_m_ggplot_st2 = all_spp_m_ggplot[regexpr("St2",all_spp_m_ggplot[,3])>0,]
all_spp_m_ggplot_PRM = all_spp_m_ggplot[regexpr("PRM",all_spp_m_ggplot[,3])>0,]

###plot - St1
x = colorRampPalette(brewer.pal(12,"Paired"))
p1=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampSt1 samples)",fill = "Taxonomy (>1%)") +
  theme_bw() + 
  theme(legend.position="none",plot.margin = margin(0.2, 8, 0.2, 0.2, "cm"),
  plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_bar(aes(y = fraction, x = sample, fill = species),
           data = all_spp_m_ggplot_st1,stat="identity") + ylab("fraction of annotated functions") + scale_fill_manual(values = x(nrow(refseq.all.top1))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#scale_fill_brewer(palette = "Blues") +

###plot - St2
p2=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampSt2 samples)",fill = "Taxonomy (>1%)") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_bar(aes(y = fraction, x = sample, fill = species),
           data = all_spp_m_ggplot_st2,stat="identity") + ylab("fraction of annotated functions") + scale_fill_manual(values = x(nrow(refseq.all.top1))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

###plot - PRM
p3=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampPRM samples)",fill = "Taxonomy (>1%)") +
  theme_bw() +
  theme(legend.position="none",plot.margin = margin(0.2, 8, 0.2, 0.2, "cm"),
        plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  scale_color_discrete(guide = guide_legend(override.aes = list(color = "white"))) +
  geom_bar(aes(y = fraction, x = sample, fill = species),
           data = all_spp_m_ggplot_PRM,stat="identity") + ylab("fraction of annotated functions") + scale_fill_manual(values = x(nrow(refseq.all.top1))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#+ scale_fill_brewer(palette = "Blues")

#all three (dimensions in inches)
dev.new(width=14, height=15,noRStudioGD = TRUE)
#pdf('figures/Champlain_functions_barplot.pdf',width=14, height=15)
grid.arrange(p1,p2,p3, ncol = 1)
dev.print(device=pdf, paste("figures/Champlain_functions_barplot_level",level,".pdf",sep =""), onefile=FALSE)
dev.off()

###Sandbox ----


p1 = plot
p2 = plot + theme(legend.position="none",plot.margin = margin(0, 2, 0, 0, "cm"))


if(1==2)
{
  
  #these must be above the 2% in at least one comparison.
  all_spp_one_m = unique(sort(all_spp_one))
  
  #these must be present in the topTEN
  all_spp_m = data.frame(unique(sort(all_spp)),stringsAsFactors = F)
  
  ###get the data for the topTEN
  for(i in 1:nrow(all.tsv))
  {
    refseq.all = read.table(paste(all.tsv[i,1],"_modif",sep = ""),header = T,stringsAsFactors = F)
    for(j in 1:nrow(all_spp_m))
    {
      temp = refseq.all[refseq.all[,3] == all_spp_m[j,1],1]
      if(length(temp) == 1) all_spp_m[j,i+1] = temp
      if(length(temp) == 0) all_spp_m[j,i+1] = 0
    }
    #get a proper shortname for the graph
    replicate = substring(strsplit(all.tsv[i,1],split = "WatPhotz_")[[1]][2],4,4)
    location = substring(strsplit(all.tsv[i,1],split = "Champ")[[1]][2],1,3)
    date = strsplit(all.tsv[i,1],split = "-")[[1]][2]
    all_spp_shortnames = c(all_spp_shortnames,rep(paste(date,location,replicate,sep = "_"),nrow(all_spp_m)))
  }
  all_spp_m_ggplot = data.frame(unlist(all_spp_m[,2:66]))
  temp = all_spp_m_ggplot
  all_spp_m_ggplot[,2] = rep(all_spp_m[,1],65)
  all_spp_m_ggplot[,3] = all_spp_shortnames
  colnames(all_spp_m_ggplot) = c("fraction","species","sample")
  
  #GET ONLY THE DATA FOR THE TOP % 
  all_spp_m_ggplot_top12 = NULL
  
  for(i in 1:length(all_spp_one_m))
  {
    all_spp_m_ggplot_top12 = rbind(all_spp_m_ggplot_top12,all_spp_m_ggplot[all_spp_m_ggplot[,2] == all_spp_one_m[i],])
  }
  #these are percentage values
  all_spp_m_ggplot_top12[,1] = all_spp_m_ggplot_top12[,1]/100
  
  #split into 3 locations.
  all_spp_m_ggplot_top12_st1 = all_spp_m_ggplot_top12[regexpr("St1",all_spp_m_ggplot_top12[,3])>0,]
  all_spp_m_ggplot_top12_st2 = all_spp_m_ggplot_top12[regexpr("St2",all_spp_m_ggplot_top12[,3])>0,]
  all_spp_m_ggplot_top12_PRM = all_spp_m_ggplot_top12[regexpr("PRM",all_spp_m_ggplot_top12[,3])>0,]
  
  
  ###plot - St1
  x = colorRampPalette(brewer.pal(12,"Paired"))
  p1=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampSt1 samples)",fill = "Taxonomy") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
    geom_bar(aes(y = fraction, x = sample, fill = species),
             data = all_spp_m_ggplot_top12_st1,stat="identity") + ylab("fraction of annotated functions") + scale_fill_manual(values = x(length(all_spp_one_m))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  #scale_fill_brewer(palette = "Blues") +
  
  ###plot - St2
  p2=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampSt2 samples)",fill = "Taxonomy") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
    geom_bar(aes(y = fraction, x = sample, fill = species),
             data = all_spp_m_ggplot_top12_st2,stat="identity") + ylab("fraction of annotated functions") + scale_fill_manual(values = x(length(all_spp_one_m))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ###plot - PRM
  p3=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampPRM samples)",fill = "Taxonomy") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
    geom_bar(aes(y = fraction, x = sample, fill = species),
             data = all_spp_m_ggplot_top12_PRM,stat="identity") + ylab("fraction of annotated functions") + scale_fill_manual(values = x(length(all_spp_one_m))) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  #+ scale_fill_brewer(palette = "Blues")
  
  #all three (dimensions in inches)
  dev.new(width=14, height=15,noRStudioGD = TRUE)
  #pdf('figures/Champlain_functions_barplot.pdf',width=14, height=15)
  grid.arrange(p1,p2,p3, ncol = 1)
  dev.print(device=pdf, "figures/Champlain_functions_barplot_level1.pdf", onefile=FALSE)
  dev.off()
}
