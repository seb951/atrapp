

setwd("/Users/jerry/Documents/CSBQ/shapiro")

library(lubridate)
library(ggplot2)
#load env. data
env = read.table("results/GRDI_ECO_formated_Metadata_CNRC_2016_for_stats_Jesse_v2.Final.txt",
                 header = F, stringsAsFactors = F,sep = "\t",comment.char = "",skip = 1)

header  = read.table("results/GRDI_ECO_formated_Metadata_CNRC_2016_for_stats_Jesse_v2.Final.txt",
                     header = F, stringsAsFactors = F,sep = "\t",comment.char = "",nrows = 1)


env = env[,1:91]

#this is because some decimals are seperated by commas....
for(i in 1:91)
{
  env[,i] = gsub(",",".",env[,i])
}

colnames(env) = header[1,1:91]
env[,4] = mdy(paste(env[,4]," 2016"))

##remove this funky data point...
#env = env[env[,4] != "2016-09-15",]


#merge the replicated (get the mean)
env_merged = data.frame(matrix(NA,ncol = 91, nrow = length(unique(env[,4]))))

env_merged[,4] = unique(env[,4])
colnames(env_merged) = colnames(env)
for(i in 1:nrow(env_merged))
    {
        temp = env[env[,4] == env_merged[i,4],]  
      for(j in c(1,2,3,5:91))
        {
        temp2 = as.numeric(temp[,j]) #change to numerics
        temp3 = temp2[!is.na(temp2)] #keep only numerics
        if(length(temp3)>0) env_merged[i,j] = mean(temp3) # mean of replicates if possible...
        }
}

colnames(env_merged) = colnames(env) 

env_ggplot = data.frame(values = as.numeric(unlist(env_merged[,c(6,21,27,30,33)])),date = env_merged[,4],var1 = c(rep("N",14),rep("P",14),rep("Chl",14),rep("intrac_tox",14),rep("extrac_tox",14)))
env_ggplot$var2 <- factor(env_ggplot$var1, labels = c("Chlorophyl A\n(ug/L)", "Nitrogen\n(mg N/L)", "Phosphorus\n(ug P/L)","Toxins (intracellular)\n(ug/L)","Toxins (extracellular)\n(ug/L)"))


graph = "pdf"
#graph = "png"

if(graph == "pdf") { dev.new(width=6, height=6,noRStudioGD = TRUE) }
if(graph == "png") { png("figures/environmental_data.png",width=6, res =400,height=4,units= 'in') }
  
  p = ggplot() + labs(title = "Lake Champlain") +
  geom_point(aes(x = date, y = values),na.rm=T,data = env_ggplot) + 
    theme_bw() + 
    ylab("Concentration") +
    xlab("Sampling date") +
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"),strip.text.y = element_text(size = 6)) +
    scale_x_date(date_breaks = "months" , date_labels = "%b") +
    facet_grid(rows=vars(var2),scales="free")
 
  p
  
#PDF
if(graph == "pdf")  {  dev.print(device=pdf, "figures/environmental_data.pdf", onefile=FALSE)  }

dev.off()

###chlorophyll only.
chl = env_ggplot[env_ggplot$var1=="Chl",]
chl = chl[!is.na(chl[,1]),] 
#chl[,1] = log(chl[,1],10)

dev.new(width=6, height=6,noRStudioGD = TRUE)
#graph parameters
p = ggplot() + theme(title = element_text(size = 26), axis.text.x = element_text(size = 16),axis.text.y = element_text(size = 16),axis.title =element_text(size = 26)) +
  geom_point(aes(x = date, y = values),na.rm=T,data = chl,col = "#808000ff",size = 5) +
  geom_line(aes(x = date, y = values),na.rm=T,data = chl,col = "#808000ff",size = 2) +
  labs(x = "Sampling date", y = "Chlorophyll A (ug/L)") +
  scale_y_continuous(trans='log10') +
  scale_x_date(date_breaks = "2 week" , date_labels = "%b %d")

  p
  

dev.print(device=pdf, "figures/chloro.pdf", onefile=FALSE) 

###cyano vs. env ----
#(recall that next section is contigent on running the barplots_cyano.R script...)

#In this section, I will first split the cyano annotation data into the 3 locations
#Then, merge the replicate dates
#Then, add the N,P,Cl,Toxin env. data
#With this in mind, I'll look for correlations...

##split the cyano annotation data into the 3 locations
all_spp_cyano = t(all_spp_m[,2:66]) #transpose matrix to get the same format as env. data.
all_location = list() #empty list
for(sample in c("St1","St2","PRM"))
  {
  single_location = all_spp_cyano[regexpr(sample,all.tsv[,1])>0,] #samples
  single_location_dates = dates[regexpr(sample,all.tsv[,1])>0] #dates
  single_location_all.tsv = all.tsv[regexpr(sample,all.tsv[,1])>0,] #names
  unique_dates = unique(single_location_dates) #unique_dates

  #merge the replicates (get the mean)
  single_location_merged = data.frame(matrix(NA,ncol = 36, nrow = length(unique_dates)))

    for(i in 1:nrow(single_location_merged))
      {
        temp = single_location[single_location_dates == unique_dates[i],]  
        single_location_merged[i,] = colMeans(temp) #mean
    }
  
  #add 5 columns for env. data.
  single_location_merged = cbind(single_location_merged,0,0,0,0,0)
  
  colnames(single_location_merged) = c(all_spp_m[,1],colnames(env_merged[,c(6,21,27,30,33)]))
  rownames(single_location_merged) = unique_dates
  
  ####put the resulting dataframes into a list for easy handling.
  all_location[[sample]] = single_location_merged
  
  #add env. data in last 5 columns.
  for(i in 1:nrow(all_location[[sample]]))
    {
        temp = env_merged[rownames(all_location[[sample]])[i] == env_merged[,4],c(6,21,27,30,33)]
        if(nrow(temp)==1) all_location[[sample]][i,37:41] = temp
    }
  }





###Look for correlations with sum of all cyanotoxins
#bof
#There just isn't much, in part because there are so few data points AND data is skewed by the 20160915 sample...

png("figures/cyano_env_correl.png",width=10, res =400,height=10,units= 'in')
par(mfrow=c(4,4))
st2 = all_location[["St2"]]
prm = all_location[["PRM"]]
st1 = all_location[["St1"]]
plot(rowSums(st2[,1:36]),st2[,37])
plot(rowSums(st2[,1:36]),st2[,38])
plot(rowSums(st2[,1:36]),st2[,39])
plot(rowSums(st2[,1:36]),st2[,40])
plot(rowSums(st2[,1:36]),st2[,41])

plot(rowSums(st1[,1:36]),st1[,37])
plot(rowSums(st1[,1:36]),st1[,38])
plot(rowSums(st1[,1:36]),st1[,39])
plot(rowSums(st1[,1:36]),st1[,40])
plot(rowSums(st1[,1:36]),st1[,41])

plot(rowSums(prm[,1:36]),prm[,37])
plot(rowSums(prm[,1:36]),prm[,38])
plot(rowSums(prm[,1:36]),prm[,39])
plot(rowSums(prm[,1:36]),prm[,40])
plot(rowSums(prm[,1:36]),prm[,41])

dev.off()
