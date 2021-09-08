#Malaria-like example (Figure 6 & 7) 
library(doParallel)
library(parallel)
library(foreach)

#Load functions first
source("1 Model Functions.R")

#color palette to use for figures.
library(wesanderson)
cols<-c("black",wes_palette("Darjeeling1")[c(3,2,1,4,5)])

#primary comparison

paramsflu<-makeparams(months=120,h2current=1,h2past=0.2,PI=20,infectedT0=10,Pmr=0.10,Pmt=0.05,Pct=0,Pcr=0.8,Pim=0.90) #flu-like
params4H2<-makeparams(h2current=1,h2past=0.2,Pim=0,PI=0.2,Pmr=0.0005,Pmt=0.0001,Pcr=0.04,Pct=0.02,K=5000,dieexp=0.004,reprodrate=0.010) #Helminth-like
params4<-makeparams(h2current=1,h2past=0.2,Pim=0.05) #Malaria-like
params4I<-makeparams(h2current=1,h2past=0.2,Pim=0.2,K=5000,dieexp=0.004,reprodrate=0.010) #Malaria-like high immunity
allparams<-list(paramsflu,params4H2,params4,params4I)

#Visualize the parasite parameters. A is Resist, B is Tolerate
paramscheck(params2)


#This is set up to run these parameter sets in parallel, with output to the viewer. Hit refresh in the viewer to update model progress.
Sys.sleep(0.1) #keeps makecluster from crashing
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

htmlFile <- file.path(tempdir(),"index.html")
zz <- file(htmlFile, open = "wt")
cat("Running ",length(allparams)," Models<br>Hit Refresh To See Progress<br>",file=zz)
close(zz)
# (code to write some content to the file)
viewer <- getOption("viewer")
viewer(htmlFile)
allpopsM<-foreach(p=1:length(allparams)) %dopar% {
  zz <- file(htmlFile, open = "at")
  sink(file = zz, append = TRUE, type ="message", split = FALSE)
  set.seed(65629) #Seed is arbitrary. Change to see other versions of the same.
  pops1<-runmodel(allparams[[p]],model=p)
  sink()
  pops1
}
stopCluster(cl)
#Optional save and load of models, so they don't need to be rerun each time
save(allpopsM,file="savedruns/allpopsMaternal",compress=TRUE)
#load(file="savedruns/allpopsMaternal")
load(file="savedruns/allpops")
load(file="savedruns/allpopsHelminth2")

#Plot figures
pdf("../inheritance-immunity paper/figures/sup/strategyfigMaternalExperience.pdf",width=7, height=4,pointsize=12)
strategyreport(c(allpops[4],allpopsM[3],allpopsH2[4],allpopsM[2]),n=50*12,cols=cols,popmax=c(2000,2000,4500,4500))
dev.off()

pdf("../inheritance-immunity paper/figures/sup/maternalsurvivalfigSup.pdf",width=3.43, height=5.2,pointsize=12)
layout(matrix(c(1,2,3),ncol=1),heights=c(1,1,0.13))
par(mar=c(1.75,4,0.75,1))
survplot(c(allpops[2:4],allpopsM[3]),n1=20*12+1,n2=120*12,cols=c("black",cols[3:5]),labs=c("No Tolerance","Tolerance","Maternal Effect","Maternal Experience"))
mtext("A",side=3,line=-1,at=-20,cex=1.5)
survplot(c(allpopsH2[2:4],allpopsM[2]),n1=20*12+1,n2=120*12,cols=c("black",cols[3:5]),labs=c("No Tolerance","Tolerance","Maternal Effect","Maternal Experience"))

mtext("B",side=3,line=-1,at=-20,cex=1.5)
axis(1,at=seq(0,100,20))
mtext("Age (years)",side=1,line=2,cex=0.5)
dev.off()
