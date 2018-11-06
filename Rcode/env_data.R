

setwd("/Users/jerry/Documents/CSBQ/shapiro")

library(lubridate)

#load env. data
env = read.table("results/GRDI_ECO_formated_Metadata_CNRC_2016_for_stats_Jesse_v1seb.txt",
                 header = F, stringsAsFactors = F,sep = "\t",comment.char = "",skip = 1)

header  = read.table("results/GRDI_ECO_formated_Metadata_CNRC_2016_for_stats_Jesse_v1seb.txt",
                   header = F, stringsAsFactors = F,sep = "\t",comment.char = "",nrows = 1)

env = env[,1:91]
colnames(env) = header[1,1:91]
env[,4] = mdy(paste(env[,4]," 2016"))

for(graph in c("pdf", "png"))
{
if(graph == "pdf") { dev.new(width=7, height=7,noRStudioGD = TRUE) }
if(graph == "png") { png("figures/environmental_data.png",width=7, res =400,height=7,units= 'in') }
par(mfrow = c(3,1))
#N (mg N/L)
plot(env[env[,3]=="station1",c(4,6)],type = "p",col = "darkblue",pch =20,cex =2,
     ylim = c(0,max(env[env[,3]=="station1",6])*1.1))
#points(env[env[,3]=="station1",c(4,6)],type = "l",col = "black",pch =20)
x = env[env[,3]=="station1",c(4)]
axis(3,at=x,labels=x,cex.axis=0.8,las=2,tick=F,line=-7)
#text(env[40,4]-15,env[40,6],env[40,4],xpd = T,col="darkblue")

#P (µg P/L)
plot(env[env[,3]=="station1",c(4,21)],type = "p",col = "darkblue",pch =20,cex =2,
     ylim = c(0,max(env[env[,3]=="station1",21])*1.1))
axis(3,at=x,labels=x,cex.axis=0.8,las=2,tick=F,line=-7)
#points(env[env[,3]=="station1",c(4,21)],type = "l",col = "black",pch =20)
#text(env[40,4]-15,env[40,21],env[40,4],xpd = T,col="darkblue")

#Chlorophyll (filtration; ug/L)
plot(env[env[,3]=="station1",c(4,27)],type = "p",col = "darkblue",pch =20,cex =2,
     ylim = c(0,600))
axis(3,at=x,labels=x,cex.axis=0.8,las=2,tick=F,line=-7)
#points(env[env[,3]=="station1",c(4,27)],type = "l",col = "black",pch =20)
#text(env[40,4]-15,as.numeric(env[40,27]),env[40,4],xpd = T,col="darkblue")


#PDF
if(graph == "pdf")  {  dev.print(device=pdf, "figures/environmental_data.pdf", onefile=FALSE)  }

dev.off()
}




