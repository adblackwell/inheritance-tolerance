makeparams<-function(n0=1000, #starting population size
                     months=1500, #how long to run the model
                     h2inf=0, #percent of mom's current infection value inherited
                     h2past=0, #percent of mom's past infection history inherited
                     K=3000, #maximum population size
                     reprodrate=0.015, #Monthy risk of reproduction for reproductive age individual
                     R0=0.5, #Poisson distribution lambda (mean) number of indivduals an infection infects each month
                     infectedT0=10, #number of individuals to infect at time zero
                     startingstrategy=1, #1 = resist, 2=tolerate
                     Pmr=0.025, #Monthly probablity of mortality when resisting
                     Pmt=0.0001, #Monthly probablity of mortality when tolerating
                     Pcr=0.20, #Monthly probablity of clearance when resisting
                     Pct=0.01, #Monthly probablity of clearance when tolerating
                     Pim=0.2, #Probability of gaining immunity after clearance
                     PIdecline=0.98, #Proportion of immunological memory retained each month. when not infected
                     PNIdecline=1,
                     treat=NA, #cut from final, allows treatment
                     treatefficacy=0, #cut from final, allows treatment
                     inheritstrategy=FALSE, #whether offspring directly inherit mom's strategy (rather than just her history)
                     flipthreshold=0.75, #proportion of posterior past this threshold determines prob of flipping
                     dieexp=0.002, #exponent for non-infection mortality curve
                     priorMon=0, #infants are born with an inate prior equal to this many months not infected
                     infectionweight=1, #how much infection months count relative to unifected in updating prior
                     noflip=FALSE, #When TRUE, there is no startegy switching
                     trackinfections=TRUE){ #When true, infections are recorded, but slows down models
  params<-as.list(environment())
}

library(doParallel)
library(parallel)
library(foreach)

#primary comparison?
params1<-makeparams(h2inf=0,h2past=0,Pim=0.05,noflip=TRUE)
params2<-makeparams(h2inf=0,h2past=0,Pim=0.05)
params3<-makeparams(h2inf=1,h2past=0.1,Pim=0.05)
params4<-makeparams(h2inf=0,h2past=0,Pim=0.05,infectedT0=0,Pmr=0,Pmt=0)
allparams<-list(params1,params2,params3,params4)
# 

params5<-makeparams(h2inf=0,h2past=0,Pim=0.05,noflip=TRUE,Pmr=0.015)
params6<-makeparams(h2inf=0,h2past=0,Pim=0.05,Pmr=0.015)
params7<-makeparams(h2inf=1,h2past=0.1,Pim=0.05,Pmr=0.015)
# #params4<-makeparams(h2inf=0,h2past=0,Pim=0.8,noflip=TRUE)
# #params5<-makeparams(h2inf=1,h2past=0,Pim=0.8)
# 
allparams<-list(params1,params2,params3,params4,params5,params6,params7)

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

htmlFile <- file.path(tempdir(),"index.html")
zz <- file(htmlFile, open = "wt")
cat("Running ",length(allparams)," Models<br>",file=zz)
close(zz)
# (code to write some content to the file)
viewer <- getOption("viewer")
viewer(htmlFile)
allpops<-foreach(p=1:length(allparams)) %dopar% {
  zz <- file(htmlFile, open = "at")
  sink(file = zz, append = TRUE, type ="message", split = FALSE)
  set.seed(65629)
  pops1<-runmodel(allparams[[p]],model=p)
  sink()
  pops1
}
stopCluster(cl)
save(allpops,file="allpops4.3",compress=TRUE)
#exit()
load(file="allpops4.3")
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

