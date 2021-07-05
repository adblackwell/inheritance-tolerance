#Models for publication. 
library(doParallel)
library(parallel)
library(foreach)

#Load functions first
source("1 Model Functions 3.0.R")

#color palette to use for figures.
library(wesanderson)
cols<-c("black",wes_palette("Darjeeling1")[c(3,2,1,4,5)])

#primary comparison

params1<-makeparams(h2current=0,h2past=0,Pim=0.9,infectedT0=0,Pmr=0,Pmt=0,K=5000,dieexp=0.004,reprodrate=0.010) #Base with no infection
params2<-makeparams(h2current=0,h2past=0,Pim=0.9,noflip=TRUE,K=5000,dieexp=0.004,reprodrate=0.010) #No Tolerance
params3<-makeparams(h2current=0,h2past=0,Pim=0.9,K=5000,dieexp=0.004,reprodrate=0.010) #Tolerance but no inheritance
params4<-makeparams(h2current=1,h2past=0.1,Pim=0.9,K=5000,dieexp=0.004,reprodrate=0.010) #Tolerance and inheritance
allparams<-list(params1,params2,params3,params4)

#Visualize the parasite parameters. A is Resist, B is Tolerate
paramscheck(params2)


#This is set up to run these parameter sets in parallel, with output to the viewer. Hit refresh in the viewer to update model progress.
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

htmlFile <- file.path(tempdir(),"index.html")
zz <- file(htmlFile, open = "wt")
cat("Running ",length(allparams)," Models<br>Hit Refresh To See Progress<br>",file=zz)
close(zz)
# (code to write some content to the file)
viewer <- getOption("viewer")
viewer(htmlFile)
allpops<-foreach(p=1:length(allparams)) %dopar% {
  zz <- file(htmlFile, open = "at")
  sink(file = zz, append = TRUE, type ="message", split = FALSE)
  set.seed(65629) #Seed is arbitrary. Change to see other versions of the same.
  pops1<-runmodel(allparams[[p]],model=p)
  sink()
  pops1
}
stopCluster(cl)
#Optional save and load of models, so they don't need to be rerun each time
save(allpops,file="allpopsImmunity",compress=TRUE)
load(file="allpopsImmunity")

#Plot figures
pdf("figures/immunefigSup.pdf",width=7, height=4,pointsize=12)
strategyreport(allpops[2:4],n=100*12,cols=cols,popmax=4500)
dev.off()

pdf("figures/survivalfigSup.pdf",width=3.43, height=5.2,pointsize=12)
layout(matrix(c(1,2,3),ncol=1),heights=c(1,1,0.13))
par(mar=c(1.75,4,0.75,1))
survplot(allpops[1:4],n1=20*12+1,n2=120*12,cols=c("black",cols[3:5]),labs=c("No Infection","No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("A",side=3,line=-1,at=-20,cex=1.5)
mortalityplot(allpops[1:4],n1=20*12+1,n2=120*12,cols=c("black",cols[3:5]),labs=c("No Infection","No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("B",side=3,line=-1,at=-20,cex=1.5)
dev.off()

pdf("figures/prevalenceSup2.pdf",width=3.43, height=6.2,pointsize=12)
layout(matrix(c(1,2,3,4,5),ncol=1),heights=c(1,1,1,1,0.13))
par(mar=c(1.75,4,0.75,1))
prevplot(allpops[2:4],n1=0*12+1,n2=10*12,cols=cols[3:5],labs=c("No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("A",side=3,line=-1,at=-20,cex=1.5)
prevplot(allpops[2:4],n1=10*12+1,n2=20*12,cols=cols[3:5],labs=c("No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("B",side=3,line=-1,at=-20,cex=1.5)
prevplot(allpops[2:4],n1=40*12+1,n2=50*12,cols=cols[3:5],labs=c("No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("C",side=3,line=-1,at=-20,cex=1.5)
prevplot(allpops[2:4],n1=90*12+1,n2=100*12,cols=cols[3:5],labs=c("No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
mtext("D",side=3,line=-1,at=-20,cex=1.5)
axis(1,at=seq(5,95,10))
dev.off()
