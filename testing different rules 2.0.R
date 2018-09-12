makeparams<-function(n0=1000,
                     months=1200,
                     h2inf=0,
                     h2past=0,
                     K=5000,
                     reprodrate=0.02,
                     R0=1,
                     infectedT0=1,
                     startingstrategy=1,
                     Pmr=0.02,
                     Pmt=0.002,
                     Pcr=0.10,
                     Pct=0.05,
                     Pim=0.4,
                     flip=1,
                     PIdecline=0.98,
                     treat=NA,
                     treatefficacy=0,
                     inheritstrategy=FALSE,
                     dieexp=0.0015,
                     noflip=FALSE){
  params<-list(n0=n0, months=months, h2inf=h2inf, h2past=h2past, K=K, reprodrate=reprodrate, R0=R0, infectedT0=infectedT0, startingstrategy=startingstrategy, Pmr=Pmr, Pmt=Pmt, Pcr=Pcr, Pct=Pct, Pim=Pim, flip=flip, PIdecline=PIdecline, treat=treat, treatefficacy=treatefficacy, inheritstrategy=inheritstrategy,dieexp=dieexp,noflip=noflip)
}

#primary comparison?
params1<-makeparams(h2inf=0,h2past=0,Pim=0.05,flip=0,noflip=TRUE)
params2<-makeparams(h2inf=0,h2past=0,Pim=0.05,flip=1)
params3<-makeparams(h2inf=1,h2past=0.1,Pim=0.05,flip=1)
params4<-makeparams(h2inf=0,h2past=0,Pim=0.8,flip=0,noflip=TRUE)
params5<-makeparams(h2inf=1,h2past=0.1,Pim=0.8,flip=1)
allparams<-list(params1,params2,params3,params4,params5)

library(doParallel)
library(parallel)
library(foreach)
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

allpops<-foreach(p=1:length(allparams)) %dopar% {
  set.seed(65628)
  pops1<-runmodel(allparams[[p]])
}
#stopCluster(cl)
save(allpops,file="allpops4.2",compress=TRUE)
#exit()
load(file="allpops4.2")
library(wesanderson)
cols<-c("black",wes_palette("Darjeeling")[c(3,2,1,4,5)])
strategyreportcompare(allpops,1200,cols=cols)

params4<-makeparams(h2inf=0,h2past=0,Pim=0.95,flip=0,noflip=TRUE)
params5<-makeparams(h2inf=1,h2past=0.1,Pim=0.95,flip=1)
allparams<-list(params4,params5)

allpops2<-foreach(p=1:length(allparams)) %dopar% {
  set.seed(65628)
  pops1<-runmodel(allparams[[p]])
}
strategyreportcompare(allpops2,1200,cols=cols)

#primary comparison?
params1<-makeparams(h2inf=1,h2past=0.1,Pim=0.1,flip=1.5)
params2<-makeparams(h2inf=1,h2past=0.1,Pim=0.4,flip=1.5)
params3<-makeparams(h2inf=1,h2past=0.1,Pim=0.8,flip=1.5)
params4<-makeparams(h2inf=1,h2past=0,Pim=0.1,flip=1.5)
params5<-makeparams(h2inf=1,h2past=0,Pim=0.4,flip=1.5)
params6<-makeparams(h2inf=1,h2past=0,Pim=0.8,flip=1.5)

allparams<-list(params1,params2,params3,params4,params5,params6)
#allparams<-list(params1,params2,params3)
# 
# remoter::client("128.111.165.102")
# c2s(runmodel)
# c2s(allparams)
# setwd("D:/Dropbox/My Files/Anthropology/UCSB Studies/Inheritance of Immunity/Immunity Inheritance Model")
library(doParallel)
library(parallel)
library(foreach)
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

allpops<-foreach(p=1:length(allparams)) %dopar% {
  set.seed(65625)
  pops1<-runmodel(allparams[[p]])
}
#stopCluster(cl)
save(allpops,file="allpops4.2",compress=TRUE)
#exit()
load(file="allpops4.2")
library(wesanderson)
cols<-c("black",wes_palette("Darjeeling")[c(3,2,1,4,5)])
strategyreportcompare(allpops,1200,cols=cols)

#secondary comparison?
params1<-makeparams(h2inf=0,h2past=0,Pim=0.2,flip=0.01)
params1<-makeparams(h2inf=0,h2past=0,Pim=0.2)
params2<-makeparams(h2inf=1,h2past=0.1,Pim=0.2)
#params3<-makeparams(h2inf=0,flip=1.5,Pim=0.2,Pmr=0.06)
#params4<-makeparams(h2inf=1,flip=1.5,Pim=0.2,Pmr=0.06)

allparams<-list(params1,params2)
#allparams<-list(params1,params2,params3)

library(doParallel)
library(parallel)
library(foreach)
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

allpops<-foreach(p=1:length(allparams)) %dopar% {
  set.seed(65625)
  pops1<-runmodel(allparams[[p]])
}
stopCluster(cl)

library(wesanderson)
cols<-c("black",wes_palette("Darjeeling")[c(3,2,1,4,5)])
strategyreportcompare(allpops,1200,cols=cols)

#the above is interesting. Shows that population growth is double with the inheritance, mostly by reducing infant mortality I think (can get mortality graphs to see). Infants are fast tracked to immunity or tolerance, and not the moddle ground. 

save(allpops,file="allpopsinheritcompare")

