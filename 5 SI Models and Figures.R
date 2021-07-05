#SI Models
library(doParallel)
library(parallel)
library(foreach)

#Load functions first
source("1 Model Functions 3.0.R")

#color palette to use for figures.
library(wesanderson)
cols<-c("black",wes_palette("Darjeeling1")[c(3,2,1,4,5)])



#helminth like example
#very low transmssion

params1H0<-makeparams(h2current=0,h2past=0,Pim=0,infectedT0=0,PI=0.075,Pmr=0,Pmt=0,Pcr=0.2,Pct=0.1,K=5000,dieexp=0.004,reprodrate=0.010) #Base with no infection
params2H0<-makeparams(h2current=0,h2past=0,Pim=0,noflip=TRUE,PI=0.075,Pmr=0.0005,Pmt=0.0001,Pcr=0.04,Pct=0.02,K=5000,dieexp=0.004,reprodrate=0.010) #No Tolerance
params3H0<-makeparams(h2current=0,h2past=0,Pim=0,PI=0.075,Pmr=0.0005,Pmt=0.0001,Pcr=0.04,Pct=0.02,K=5000,dieexp=0.004,reprodrate=0.010) #Tolerance but no inheritance
params4H0<-makeparams(h2current=1,h2past=0.1,Pim=0,PI=0.075,Pmr=0.0005,Pmt=0.0001,Pcr=0.04,Pct=0.02,K=5000,dieexp=0.004,reprodrate=0.010) #Tolerance and inheritance
allparamsH0<-list(params1H0,params2H0,params3H0,params4H0)

#This is set up to run these parameter sets in parallel, with output to the viewer. Hit refresh in the viewer to update model progress.
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

htmlFile <- file.path(tempdir(),"index.html")
zz <- file(htmlFile, open = "wt")
cat("Running ",length(allparamsH)," Models<br>Hit Refresh To See Progress<br>",file=zz)
close(zz)
# (code to write some content to the file)
viewer <- getOption("viewer")
viewer(htmlFile)
allpopsH0<-foreach(p=1:length(allparamsH0)) %dopar% {
  zz <- file(htmlFile, open = "at")
  sink(file = zz, append = TRUE, type ="message", split = FALSE)
  set.seed(236) #Seed is arbitrary. Change to see other versions of the same.
  pops1<-runmodel(allparamsH0[[p]],model=p)
  sink()
  pops1
}

stopCluster(cl)
#Optional save and load of models, so they don't need to be rerun each time
save(allpopsH0,file="allpopsHelminth0",compress=TRUE)
#load(file="allpopsHelminth0")


#Note, baseline survival and reproduction have been changed to keep the population size more managable in the simulation.

params1H<-makeparams(h2current=0,h2past=0,Pim=0,infectedT0=0,PI=0.1,Pmr=0,Pmt=0,Pcr=0.2,Pct=0.1,K=5000,dieexp=0.004,reprodrate=0.010) #Base with no infection
params2H<-makeparams(h2current=0,h2past=0,Pim=0,noflip=TRUE,PI=0.1,Pmr=0.0005,Pmt=0.0001,Pcr=0.04,Pct=0.02,K=5000,dieexp=0.004,reprodrate=0.010) #No Tolerance
params3H<-makeparams(h2current=0,h2past=0,Pim=0,PI=0.1,Pmr=0.0005,Pmt=0.0001,Pcr=0.04,Pct=0.02,K=5000,dieexp=0.004,reprodrate=0.010) #Tolerance but no inheritance
params4H<-makeparams(h2current=1,h2past=0.1,Pim=0,PI=0.1,Pmr=0.0005,Pmt=0.0001,Pcr=0.04,Pct=0.02,K=5000,dieexp=0.004,reprodrate=0.010) #Tolerance and inheritance
allparamsH<-list(params1H,params2H,params3H,params4H)

#Visualize the parasite parameters. A is Resist, B is Tolerate
paramscheck(params2H,xmax1=10,xmax2=10)

#This is set up to run these parameter sets in parallel, with output to the viewer. Hit refresh in the viewer to update model progress.
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

htmlFile <- file.path(tempdir(),"index.html")
zz <- file(htmlFile, open = "wt")
cat("Running ",length(allparamsH)," Models<br>Hit Refresh To See Progress<br>",file=zz)
close(zz)
# (code to write some content to the file)
viewer <- getOption("viewer")
viewer(htmlFile)
allpopsH<-foreach(p=1:length(allparamsH)) %dopar% {
  zz <- file(htmlFile, open = "at")
  sink(file = zz, append = TRUE, type ="message", split = FALSE)
  set.seed(236) #Seed is arbitrary. Change to see other versions of the same.
  pops1<-runmodel(allparamsH[[p]],model=p)
  sink()
  pops1
}

stopCluster(cl)
#Optional save and load of models, so they don't need to be rerun each time
save(allpopsH,file="allpopsHelminth",compress=TRUE)
#load(file="allpopsHelminth")


#helminth like example 2, higher virulence
#Note, baseline survival and reproduction have been changed to keep the population size more managable in the simulation.

params1H2<-makeparams(h2current=0,h2past=0,Pim=0,infectedT0=0,PI=0.2,Pmr=0,Pmt=0,Pcr=0.2,Pct=0.1,K=5000,dieexp=0.004,reprodrate=0.010) #Base with no infection
params2H2<-makeparams(h2current=0,h2past=0,Pim=0,noflip=TRUE,PI=0.2,Pmr=0.0005,Pmt=0.0001,Pcr=0.04,Pct=0.02,K=5000,dieexp=0.004,reprodrate=0.010) #No Tolerance
params3H2<-makeparams(h2current=0,h2past=0,Pim=0,PI=0.2,Pmr=0.0005,Pmt=0.0001,Pcr=0.04,Pct=0.02,K=5000,dieexp=0.004,reprodrate=0.010) #Tolerance but no inheritance
params4H2<-makeparams(h2current=1,h2past=0.1,Pim=0,PI=0.2,Pmr=0.0005,Pmt=0.0001,Pcr=0.04,Pct=0.02,K=5000,dieexp=0.004,reprodrate=0.010) #Tolerance and inheritance
allparamsH2<-list(params1H2,params2H2,params3H2,params4H2)

#Visualize the parasite parameters. A is Resist, B is Tolerate
paramscheck(params2H2,xmax1=10,xmax2=10)

#This is set up to run these parameter sets in parallel, with output to the viewer. Hit refresh in the viewer to update model progress.
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

htmlFile <- file.path(tempdir(),"index.html")
zz <- file(htmlFile, open = "wt")
cat("Running ",length(allparamsH)," Models<br>Hit Refresh To See Progress<br>",file=zz)
close(zz)
# (code to write some content to the file)
viewer <- getOption("viewer")
viewer(htmlFile)
allpopsH2<-foreach(p=1:length(allparamsH2)) %dopar% {
  zz <- file(htmlFile, open = "at")
  sink(file = zz, append = TRUE, type ="message", split = FALSE)
  set.seed(236) #Seed is arbitrary. Change to see other versions of the same.
  pops1<-runmodel(allparamsH2[[p]],model=p)
  sink()
  pops1
}

stopCluster(cl)
#Optional save and load of models, so they don't need to be rerun each time
save(allpopsH2,file="allpopsHelminth2",compress=TRUE)

pdf("figures/helminthfig.pdf",width=7, height=4,pointsize=12)
strategyreport(c(allpopsH0[4],allpopsH[4],allpopsH2[4]),n=50*12,cols=cols,popmax=4500)
dev.off()

pdf("figures/helmsurvivalfig.pdf",width=7, height=5.2,pointsize=10)
layout(matrix(c(1,2,0,3,4,0,5,6,0),ncol=3,byrow=FALSE),heights=c(1,1,0.13))
par(mar=c(1.75,4,0.75,1))
survplot(allpopsH0[1:4],n1=20*12+1,n2=120*12,cols=c("black",cols[3:5]),labs=c("No Infection","No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("A",side=3,line=-1,at=-25,cex=1.5)
mortalityplot(allpopsH0[1:4],n1=20*12+1,n2=120*12,cols=c("black",cols[3:5]),labs=c("No Infection","No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("B",side=3,line=-1,at=-25,cex=1.5)


survplot(allpopsH[1:4],n1=20*12+1,n2=120*12,cols=c("black",cols[3:5]),labs=c("No Infection","No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("C",side=3,line=-1,at=-25,cex=1.5)
mortalityplot(allpopsH[1:4],n1=20*12+1,n2=120*12,cols=c("black",cols[3:5]),labs=c("No Infection","No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("D",side=3,line=-1,at=-25,cex=1.5)

survplot(allpopsH2[1:4],n1=20*12+1,n2=120*12,cols=c("black",cols[3:5]),labs=c("No Infection","No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("E",side=3,line=-1,at=-25,cex=1.5)
mortalityplot(allpopsH2[1:4],n1=20*12+1,n2=120*12,cols=c("black",cols[3:5]),labs=c("No Infection","No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("F",side=3,line=-1,at=-25,cex=1.5)

dev.off()


#SI Figure showing these three pathogen sets
paramsflu<-makeparams(months=120,h2current=1,h2past=0.1,PI=20,infectedT0=10,Pmr=0.10,Pmt=0.05,Pct=0,Pcr=0.8,Pim=0.90)
params2<-makeparams(h2current=0,h2past=0,Pim=0.05,noflip=TRUE) #No Tolerance

pdf("figures/parasiteprofiles.pdf",width=7, height=4,pointsize=12)
layout(matrix(c(1,2,0,3,4,7,5,6,0),nrow=3,byrow=FALSE),widths=c(4,4,4),heights=c(1,1,0.1))
paramscheck(paramsflu,xmax1=4,xmax2=4,nolayout=TRUE,legend="topright")
paramscheck(params2,xmax1=10,xmax2=10,nolayout=TRUE,legend=NA,letters=c("C","D"))
paramscheck(params2H,xmax1=10,xmax2=10,nolayout=TRUE,legend=NA,letters=c("E","F"))
par(mar=c(0,5,0,0.5))
frame()
text(0.5,0.7,"Years",cex=2,adj=c(0.5,0.5))
dev.off()

