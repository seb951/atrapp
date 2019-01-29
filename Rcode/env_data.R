

setwd("/Users/jerry/Documents/CSBQ/shapiro")

library(lubridate)
library(ggplot2)
#load env. data
env_old = read.table("results/GRDI_ECO_formated_Metadata_CNRC_2016_for_stats_Jesse_v1seb.txt",
                 header = F, stringsAsFactors = F,sep = "\t",comment.char = "",skip = 1)

env = read.table("results/GRDI_ECO_formated_Metadata_CNRC_2016_for_stats_Jesse_v2.Final.txt",
                 header = F, stringsAsFactors = F,sep = "\t",comment.char = "",skip = 1)


header_old  = read.table("results/GRDI_ECO_formated_Metadata_CNRC_2016_for_stats_Jesse_v1seb.txt",
                   header = F, stringsAsFactors = F,sep = "\t",comment.char = "",nrows = 1)

header  = read.table("results/GRDI_ECO_formated_Metadata_CNRC_2016_for_stats_Jesse_v2.Final.txt",
                     header = F, stringsAsFactors = F,sep = "\t",comment.char = "",nrows = 1)


env = env[,1:91]
colnames(env) = header[1,1:91]
env[,4] = mdy(paste(env[,4]," 2016"))

env_ggplot = data.frame(values = as.numeric(unlist(env[,c(6,21,27)])),date = env[,4],var1 = c(rep("N",44),rep("P",44),rep("Chl",44)))
env_ggplot$var2 <- factor(env_ggplot$var1, labels = c("Chlorophyl A\n(ug/L)", "Nitrogen\n(mg N/L)", "Phosphorus\n(ug P/L)"))



#graph = "pdf"
graph = "png"

if(graph == "pdf") { dev.new(width=6, height=4,noRStudioGD = TRUE) }
if(graph == "png") { png("figures/environmental_data.png",width=6, res =400,height=4,units= 'in') }
  
  p = ggplot() + labs(title = "Lake Champlain") +
  geom_point(aes(x = date, y = values),na.rm=T,data = env_ggplot) + 
    theme_bw() + 
    ylab("Concentration") +
    xlab("Sampling date") +
    theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
    scale_x_date(date_breaks = "months" , date_labels = "%b") +
    facet_grid(rows=vars(var2),scales="free")
 
  p
  
#PDF
if(graph == "pdf")  {  dev.print(device=pdf, "figures/environmental_data.pdf", onefile=FALSE)  }

dev.off()




