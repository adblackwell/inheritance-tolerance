#Models for publication. 
library(doParallel)
library(parallel)
library(foreach)

#Load functions first
source("1 Model Functions 3.0.R")

#color palette to use for figures.
library(wesanderson)
cols<-c("black",wes_palette("Darjeeling1")[c(3,2,1,4,5)])

## ---- modelflu ----
#set for flu. Medium cost immunity, same cost tolerance, high chance clearing, immunity
paramsflu<-makeparams(months=120,h2current=1,h2past=0.1,PI=20,infectedT0=10,Pmr=0.10,Pmt=0.05,Pct=0,Pcr=0.8,Pim=0.90)

#Visualize the parasite parameters. A is Resist, B is Tolerate
paramscheck(params2,xmax1=10,xmax2=10)
paramscheck(paramsflu,xmax1=2,xmax2=2)


set.seed(2367)
popsflu<-runmodel(paramsflu)


pdf("figures/flufig.pdf",width=3.43, height=1.7,pointsize=6)
n<-12
pops<-popsflu$pops[1:n]
infections<-popsflu$infections
params<-popsflu$params

layout(matrix(c(1,2, 0,0),ncol=2,byrow=TRUE),widths=c(4,4),heights=c(1,0.2))
bottommar<-1

par(mar=c(bottommar,4,1,1))
infected<-Reduce(rbind,lapply(pops,function(x) c(sum(x$infected==0 & x$immune==0),sum(x$immune==1),sum(x$infected>0))))
popsize<-infected[,1]+infected[,2]+infected[,3]
plot(1,1,ylim=c(0,max(popsize)),xlim=c(1,n),type="n",xaxt="n",ylab=NA)
lines(1:nrow(infected),popsize,lty=2,lwd=2,col=cols[1])
lines(1:nrow(infected),infected[,1],lty=1,lwd=2,col=cols[2])
lines(1:nrow(infected),infected[,2],lty=1,lwd=2,col=cols[3])
lines(1:nrow(infected),infected[,3],lty=1,lwd=2,col=cols[4])
ax<-axTicks(1,axp=c(0,n/12,10))
axis(1,at=ax*12,labels=ax)
mtext("Individuals",side=2,line=2.5)
mtext("Years",side=1,line=2.5)

legend("right",c("All","Susceptible","Immune","Infected"),col=cols[1:4],lty=c(2,1,1,1),lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0),inset=c(0.05,0),title=expression(bold("Status")))
mtext("A",side=3,line=-1,at=-1.7,cex=2)

par(mar=c(bottommar,4,1,1))
strategy<-Reduce(rbind,lapply(pops,function(x) c(sum(x$strategy==1 & x$immune==0),sum(x$strategy==1 & x$immune==1),sum(x$strategy==2))))
plot(1,1,ylim=c(0,1),xlim=c(1,n),type="n",ylab=NA,xaxt="n")
lines(1:nrow(strategy),strategy[,1]/popsize,lty=1,lwd=2,col=cols[1])
lines(1:nrow(strategy),strategy[,2]/popsize,lty=1,lwd=2,col=cols[3])
lines(1:nrow(strategy),strategy[,3]/popsize,lty=1,lwd=2,col=cols[5])
axis(1,at=ax*12,labels=ax)
mtext("Proportion",side=2,line=2.5)
mtext("Years",side=1,line=2.5)

legend("right",c("Fight","Immune","Tolerate"),col=cols[c(1,3,5)],lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0),inset=c(0.05,0),title=expression(bold("Strategy")))
mtext("B",side=3,line=-1,at=-1.7,cex=2)

dev.off()
