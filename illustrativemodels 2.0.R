library(wesanderson)
cols<-c("black",wes_palette("Darjeeling1")[c(3,2,1,4,5)])

## ---- modelflu ----
#set for flu. Medium cost immunity, same cost tolerance, high chance clearing, immunity
paramsflu<-list(n0=1000, 
                months=120,
                h2inf=1,
                h2past=0.1,
                K=5000,
                reprodrate=0.02,
                R0=20,
                infectedT0=10,
                startingstrategy=1,
                Pmr=0.10,
                Pmt=0.05,
                Pcr=0.8,
                Pct=0,
                Pim=0.90,
                PIdecline=0.98,
                dieexp=0.0015,
                priorMon=12,
                infectionweight=2,
                treat=NA,
                treatefficacy=0,
                inheritstrategy=FALSE,
                noflip=FALSE)
#seed is just to produce a consistent graph each time. Randomizing the seed will produce similar, but not 100% identical results.


set.seed(2367)
popsflu<-runmodel(paramsflu)
#tiff("flufig.tif",width=2000,height=1000,res=400,compression="lzw",pointsize=8)
#pdf("figures/flufig.pdf",width=3.43, height=1.7,pointsize=6)

library(wesanderson)
cols<-c("black",wes_palette("Darjeeling1")[c(3,2,1,4,5)])
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

#par(mar=c(bottommar,0,0.2,1))
#plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
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

#par(mar=c(bottommar,0,0.2,1))
#plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
legend("right",c("Fight","Immune","Tolerate"),col=cols[c(1,3,5)],lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0),inset=c(0.05,0),title=expression(bold("Strategy")))
mtext("B",side=3,line=-1,at=-1.7,cex=2)

dev.off()

## ---- model123 ----

#set for malaria-like. High cost of immunity, low chance of clearing, median cost tolerance, low prob immunity
params<-list(  n0=1000,
               months=1200,
               h2inf=1,
               h2past=0,
               K=10000,
               reprodrate=0.02,
               R0=2,
               infectedT0=1,
               startingstrategy=1,
               Pmr=0.02,
               Pmt=0.002,
               Pcr=0.20,
               Pct=0.05,
               Pim=0.02,
               flip=1,
               PIdecline=0.98,
               inheritstrategy=FALSE)
#paramscheck(params,2,2) 
set.seed(65625)
pops1<-runmodel(params)
#strategyreport3(pops1,300,cols=cols)

#same as 1 with high immunity
params2<-list(  n0=1000,
               months=1200,
               h2inf=1,
               h2past=0,
               K=10000,
               reprodrate=0.02,
               R0=2,
               infectedT0=1,
               startingstrategy=1,
               Pmr=0.02,
               Pmt=0.002,
               Pcr=0.20,
               Pct=0.05,
               Pim=0.5,
               flip=1,
               PIdecline=0.98,
               treat=NA,
               treatefficacy=0,
               inheritstrategy=FALSE)
#paramscheck(params,2,2) 
set.seed(98201)
pops2<-runmodel(params2)
#strategyreport3(pops2,1200,cols=cols)


#same as 1 with low chance of tolerance
params3<-list(  n0=1000,
               months=1200,
               h2inf=1,
               h2past=0, #think about this!!
               K=10000,
               reprodrate=0.02,
               R0=2,
               infectedT0=1,
               startingstrategy=1,
               Pmr=0.02,
               Pmt=0.002,
               Pcr=0.20,
               Pct=0.05,
               Pim=0.02,
               flip=0.437,
               PIdecline=0.98,
               treat=NA,
               treatefficacy=0,
               inheritstrategy=FALSE)
#paramscheck(params,2,2) 
set.seed(98201)
pops3<-runmodel(params3)
#strategyreport3(pops3,1200,cols=cols)
#strategyreport3compare(list(pops2,pops1,pops3),1200,cols=cols)

#same as 1 with low chance of tolerance, inheritance
params4<-list(  n0=1000,
                months=1200,
                h2inf=1,
                h2past=0,
                K=10000,
                reprodrate=0.02,
                R0=2,
                infectedT0=1,
                startingstrategy=1,
                Pmr=0.02,
                Pmt=0.002,
                Pcr=0.20,
                Pct=0.05,
                Pim=0.02,
                flip=0.437,
                PIdecline=0.98,
                treat=NA,
                treatefficacy=0,
                inheritstrategy=TRUE)
#paramscheck(params,2,2) 
set.seed(98201)
pops4<-runmodel(params4)
strategyreportcompare(list(pops2,pops1,pops3,pops4),1200,cols=cols)

## ---- modelinheritcompare ----
#0 just own experience
paramsi0<-list(  n0=1000,
                 months=1200,
                 h2inf=0,
                 h2past=0,
                 K=5000,
                 reprodrate=0.02,
                 R0=2,
                 infectedT0=1,
                 startingstrategy=1,
                 Pmr=0.02,
                 Pmt=0.002,
                 Pcr=0.20,
                 Pct=0.05,
                 Pim=0.02,
                 flip=0.437,
                 PIdecline=0.98,
                 treat=NA,
                 treatefficacy=0,
                 inheritstrategy=FALSE)
#paramscheck(params,2,2) 
set.seed(65625)
popsi0<-runmodel(paramsi0)
strategyreport3(popsi0,600,cols=cols)

#1st just inf during pregnancy
paramsi1<-list(  n0=1000,
               months=1200,
               h2inf=1,
               h2past=0,
               K=5000,
               reprodrate=0.02,
               R0=2,
               infectedT0=1,
               startingstrategy=1,
               Pmr=0.02,
               Pmt=0.002,
               Pcr=0.20,
               Pct=0.05,
               Pim=0.02,
               flip=0.437,
               PIdecline=0.98,
               treat=NA,
               treatefficacy=0,
               inheritstrategy=FALSE)
#paramscheck(params,2,2) 
set.seed(65625)
popsi1<-runmodel(paramsi1)
strategyreport3(popsi1,600,cols=cols)

#2nd inf during pregnancy and mom's longer history
paramsi2<-list(  n0=1000,
                 months=1200,
                 h2inf=1,
                 h2past=0.25,
                 K=5000,
                 reprodrate=0.02,
                 R0=2,
                 infectedT0=1,
                 startingstrategy=1,
                 Pmr=0.02,
                 Pmt=0.002,
                 Pcr=0.20,
                 Pct=0.05,
                 Pim=0.02,
                 flip=0.437,
                 PIdecline=0.98,
                 treat=NA,
                 treatefficacy=0,
                 inheritstrategy=FALSE)
#paramscheck(params,2,2) 
set.seed(65625)
popsi2<-runmodel(paramsi2)
strategyreport3(popsi2,300,cols=cols)


#3 iinherit strategy
paramsi3<-list(  n0=1000,
                 months=1200,
                 h2inf=0,
                 h2past=0,
                 K=5000,
                 reprodrate=0.02,
                 R0=2,
                 infectedT0=1,
                 startingstrategy=1,
                 Pmr=0.02,
                 Pmt=0.002,
                 Pcr=0.20,
                 Pct=0.05,
                 Pim=0.02,
                 flip=0.437,
                 PIdecline=0.98,
                 treat=NA,
                 treatefficacy=0,
                 inheritstrategy=TRUE)
#paramscheck(params,2,2) 
set.seed(65625)
popsi3<-runmodel(paramsi3)
strategyreport3(popsi3,300,cols=cols)

#4 all of the above
paramsi4<-list(  n0=1000,
                 months=1200,
                 h2inf=1,
                 h2past=0.25,
                 K=5000,
                 reprodrate=0.02,
                 R0=2,
                 infectedT0=1,
                 startingstrategy=1,
                 Pmr=0.02,
                 Pmt=0.002,
                 Pcr=0.20,
                 Pct=0.05,
                 Pim=0.02,
                 flip=0.437,
                 PIdecline=0.98,
                 treat=NA,
                 treatefficacy=0,
                 inheritstrategy=TRUE)
#paramscheck(params,2,2) 
set.seed(65625)
popsi4<-runmodel(paramsi4)
strategyreport3(popsi4,300,cols=cols)

strategyreportcompare(list(popsi0,popsi1,popsi2,popsi3,popsi4),1200,cols=cols)



## ---- modelinheritcomparehighprob ----
#0 just own experience
paramsi0<-list(  n0=1000,
                 months=1200,
                 h2inf=0,
                 h2past=0,
                 K=5000,
                 reprodrate=0.02,
                 R0=2,
                 infectedT0=1,
                 startingstrategy=1,
                 Pmr=0.02,
                 Pmt=0.002,
                 Pcr=0.20,
                 Pct=0.05,
                 Pim=0.02,
                 flip=1,
                 PIdecline=0.98,
                 treat=NA,
                 treatefficacy=0,
                 inheritstrategy=FALSE)
#paramscheck(params,2,2) 
set.seed(65625)
popsi0<-runmodel(paramsi0)
#strategyreport3(popsi0,600,cols=cols)

#1st just inf during pregnancy
paramsi1<-list(  n0=1000,
                 months=1200,
                 h2inf=1,
                 h2past=0,
                 K=5000,
                 reprodrate=0.02,
                 R0=2,
                 infectedT0=1,
                 startingstrategy=1,
                 Pmr=0.02,
                 Pmt=0.002,
                 Pcr=0.20,
                 Pct=0.05,
                 Pim=0.02,
                 flip=1,
                 PIdecline=0.98,
                 treat=NA,
                 treatefficacy=0,
                 inheritstrategy=FALSE)
#paramscheck(params,2,2) 
set.seed(65625)
popsi1<-runmodel(paramsi1)
#strategyreport3(popsi1,600,cols=cols)

#2nd inf during pregnancy and mom's longer history
paramsi2<-list(  n0=1000,
                 months=1200,
                 h2inf=1,
                 h2past=0.25,
                 K=5000,
                 reprodrate=0.02,
                 R0=2,
                 infectedT0=1,
                 startingstrategy=1,
                 Pmr=0.02,
                 Pmt=0.002,
                 Pcr=0.20,
                 Pct=0.05,
                 Pim=0.02,
                 flip=1,
                 PIdecline=0.98,
                 treat=NA,
                 treatefficacy=0,
                 inheritstrategy=FALSE)
#paramscheck(params,2,2) 
set.seed(65625)
popsi2<-runmodel(paramsi2)
#strategyreport3(popsi2,300,cols=cols)


#3 iinherit strategy
paramsi3<-list(  n0=1000,
                 months=1200,
                 h2inf=0,
                 h2past=0,
                 K=5000,
                 reprodrate=0.02,
                 R0=2,
                 infectedT0=1,
                 startingstrategy=1,
                 Pmr=0.02,
                 Pmt=0.002,
                 Pcr=0.20,
                 Pct=0.05,
                 Pim=0.02,
                 flip=1,
                 PIdecline=0.98,
                 treat=NA,
                 treatefficacy=0,
                 inheritstrategy=TRUE)
#paramscheck(params,2,2) 
set.seed(65625)
popsi3<-runmodel(paramsi3)
#strategyreport3(popsi3,300,cols=cols)

#4 all of the above
paramsi4<-list(  n0=1000,
                 months=1200,
                 h2inf=1,
                 h2past=0.25,
                 K=5000,
                 reprodrate=0.02,
                 R0=2,
                 infectedT0=1,
                 startingstrategy=1,
                 Pmr=0.02,
                 Pmt=0.002,
                 Pcr=0.20,
                 Pct=0.05,
                 Pim=0.02,
                 flip=1,
                 PIdecline=0.98,
                 treat=NA,
                 treatefficacy=0,
                 inheritstrategy=TRUE)
#paramscheck(params,2,2) 
set.seed(65625)
popsi4<-runmodel(paramsi4)
#strategyreport3(popsi4,300,cols=cols)

strategyreportcompare(list(popsi0,popsi1,popsi2,popsi3,popsi4),1200,cols=cols)


## ---- modelmicrobiota ----
params4<-list( n0=1000,
               months=1200,
               h2epi=0.5,
               K=5000,
               reprodrate=0.01,
               R0=1,
               infectedT0=20,
               startingstrategy=1,
               Pmr=0.005,
               Pmt=0,
               Pcr=0.3,
               Pct=0.01,
               Pim=0.5,
               flip=0.1,
               PIdecline=0.98,
               treat=NA,
               treatefficacy=0)
paramscheck(params,10,10)
set.seed(98201)
pops5<-runmodel(params4)
strategyreport3(pops5,1200,cols=cols,statlim=c(0,1))


## ---- modelhelminth ----
#set for helminth. Medium cost immunity, very low chance clearing, low chance immunity
params<-list(  n0=1000,
               months=1200,
               h2epi=0.25,
               K=5000,
               reprodrate=0.02,
               R0=0.25,
               infectedT0=1,
               startingstrategy=1,
               Pmr=0.005,
               Pmt=0.001,
               Pcr=0.04,
               Pct=0.01,
               Pim=0.01,
               flip=0.5,
               PIdecline=0.01,
               treat=NA,
               treatefficacy=0)
paramscheck(params,10,10) 
pops2<-runmodel(params)
strategyreport3(pops2,1200,cols=cols)

