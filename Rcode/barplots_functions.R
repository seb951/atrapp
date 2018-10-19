


setwd("/Users/jerry/Dropbox/CSBQ/shapiro")

#packages
library(ggplot2)
library(gridExtra)


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
#this loop is to fix the NO HIERARCHY FUCKING problem...
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
  }

write.table(refseq.all[,2:7],paste(all.tsv[i,1],"_modif",sep = ""),quote = T, row.names = F, col.names =T)

print(paste("file nb: ",i," of:",nrow(all.tsv),"; The time is: ",Sys.time(),sep = ""))

all_spp = c(all_spp,refseq.all$hier1[1:10]) #keep only top 10
all_spp_one = c(all_spp_one,refseq.all$hier1[as.numeric(refseq.all$percent)>2]) #keep only the ones that are at more than 2% IN AT LEAST ONE OF THE SAMPLES....
}

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

#GET ONLY THE DATA FOR THE TOP 2% 
all_spp_m_ggplot_top12 = NULL

for(i in 1:11)
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
p1=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampSt1 samples)",fill = "Taxonomy") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_bar(aes(y = fraction, x = sample, fill = species),
           data = all_spp_m_ggplot_top12_st1,stat="identity") + ylab("fraction of annotated functions") + scale_fill_brewer(palette = "Set3") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#scale_fill_brewer(palette = "Blues") +

###plot - St2
p2=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampSt2 samples)",fill = "Taxonomy") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_bar(aes(y = fraction, x = sample, fill = species),
           data = all_spp_m_ggplot_top12_st2,stat="identity") + ylab("fraction of annotated functions") + scale_fill_brewer(palette = "Set3") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

###plot - PRM
p3=ggplot() + labs(title = "Lake Champlain - all annotated species (ChampPRM samples)",fill = "Taxonomy") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_bar(aes(y = fraction, x = sample, fill = species),
           data = all_spp_m_ggplot_top12_PRM,stat="identity") + ylab("fraction of annotated functions") + scale_fill_brewer(palette = "Set3") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#+ scale_fill_brewer(palette = "Blues")

#all three (dimensions in inches)
dev.new(width=14, height=15,noRStudioGD = TRUE)
#pdf('figures/Champlain_functions_barplot.pdf',width=14, height=15)
grid.arrange(p1,p2,p3, ncol = 1)
dev.print(device=pdf, "figures/Champlain_functions_barplot.pdf", onefile=FALSE)
dev.off()






###Sandbox ---
if(1==2)
{
  
  refseq.all.temp = read.table(all.tsv[i,1],stringsAsFactors = F,quote = "",sep = "\t",comment.char = "$",na.strings="none")
  
#this loop is to fix the NO HIERARCHY FUCKING problem...
no_hier = refseq.all.temp[regexpr("NO HIERARCHY",refseq.all.temp[,1])>1, 1]
no_hier = paste(no_hier,"\tNO HIERARCHY\tNO HIERARCHY\tNO HIERARCHY",sep = "")
no_misc = refseq.all.temp[regexpr("Miscellaneous",refseq.all.temp[,1])>1, 1]
no_misc = paste(no_misc,"\tMiscellaneous\tMiscellaneous\tMiscellaneous",sep = "")
refseq.all.temp[regexpr("NO HIERARCHY",refseq.all.temp[,1])>1, 1] = no_hier
refseq.all.temp[regexpr("Miscellaneous",refseq.all.temp[,1])>1, 1] = no_misc
write.table(refseq.all.temp,"results/Subsystems_results/temp",quote = F, row.names = F, col.names =F)

#this is to solve the function COMMENT CHARACTER problem:
command = paste("awk '{gsub(/&#963;-/,\"\")}; 1' results/Subsystems_results/temp  >","results/Subsystems_results/temp1",sep = "")
system(command)

#double FUCKING spaces
command = paste("awk '{gsub(/\t\t/,\"\t\")}; 1' results/Subsystems_results/temp1  >","results/Subsystems_results/temp2",sep = "")
system(command)

}