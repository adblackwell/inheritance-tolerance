#optimal strategy should depend on how much time spent in tolerance vs. how much time spent fighting/free of infection
#remoter::client("128.111.165.102",password="Credoelvumipsum")
#should it just depend on the ratio between the two and the ratio of time infected? One might also be more stochastic...

#simple model
params<-c(
  costfight=0.04,
  costtolerate=0.004,
  probkill=0.10,
  probclear=0.01,
  probimmunity=0.05,
  reinfectionprob=0.135
)
#reinfectionprob will be R0*infections/popsize. Actually it'll be 1-(1-R0/popsize)^infected
#so for R0=2.5 and 80% infected 1-(1-2.5/1000)^800 = 0.865
#cumulative sum of mortality risk for each

library(doParallel)
library(parallel)
library(foreach)
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

#at each point, should be cumulative prob dying - cumulative prob of clearing? 
survcheck<-function(costfight,probkill,probimmunity,reinfectionprob,n){
  #lengths<-NA
  lengths<-foreach(j=1:n, .combine=c) %dopar% {
    infected<-1
    immune<-0
    previnfect<-0
    for(i in 1:840){
      if(infected==1) {
        infected <- 1-rbinom(1,1,probkill)
        if(infected==0) {
          previnfect<-previnfect+1
          immune <- rbinom(1,1,probimmunity/previnfect)
        }
      }
      if(infected==0 & immune==0) {
        infected <- rbinom(1,1,reinfectionprob)
      }
      if(infected==1) {
        if(rbinom(1,1,costfight)==1){
          lengthi<-i
          break()
        }
      }
      if(i==840) lengthi<-i
    }
    lengthi
    #lengths<-c(lengths,lengthi)
  }
return(lengths)
}

s1<-survcheck(0.04,0.10,0.05,0.865,1000)
s2<-survcheck(0.004,0.01,0,0.865,1000)

s1<-survcheck(0.04,0.20,0.02,0,1000)
s2<-survcheck(0.01,0.01,0,0,1000)


#go back and modify my model. If mom was able to be immune, should never go for tolerance strategy...
immunity<-seq(0,1,0.05)
reinf<-seq(0,1,0.025)
check<-expand.grid(immunity,reinf)
n<-2000
#mm<-function(nums) c(mean(nums),median(nums))

tol<-t(sapply(reinf,function(r) mean(survcheck(0.002,0.05,0,r,n))))

compare1<-apply(check,1,function(r) mean(survcheck(0.004,0.10,r[1],r[2],n)))

compare2<-apply(check,1,function(r) mean(survcheck(0.004,0.50,r[1],r[2],n)))

compare3<-apply(check,1,function(r) mean(survcheck(0.02,0.10,r[1],r[2],n)))

compare4<-apply(check,1,function(r) mean(survcheck(0.02,0.50,r[1],r[2],n)))

save(tol,compare1,compare2,compare3,compare4,file="optimalstrategies2")


plot(reinf,compare1[1:(length),1],type="l",ylim=c(0,840))
lines(reinf,compare1[,2],col="red")

# 
# costtol<-0.004
# costfight<-costtol*seq(1,10,1)
# probkilltol<-0.05
# probkillfight<-probkilltol*seq(1,10,1)
# immunity<-0
# #1-(1-0.207)^3  0.207 median 3 months, 0.11 means median 6 months, 0.056 means median 12 months, 0.0115 5 years,0.0058 10 years)
#reinf<-c(0,0.0058,0.0115,0.056,0.11,0.207)
# check<-expand.grid(reinf,costfight,probkillfight)
# n<-2000
# 
# tolerance<-sapply(reinf,function(r) mean(survcheck(costtol,probkilltol,immunity,r,n)))
# 
# compare1<-t(apply(check,1,function(r) c(mean(survcheck(0.004,0.10,r[1],r[2],n)), mean(survcheck(0.002,0.05,0,r[2],n)))))
# 
# compare1<-t(apply(check,1,function(r) c(mean(survcheck(0.004,0.10,r[1],r[2],n)), mean(survcheck(0.002,0.05,0,r[2],n)))))
# 
# compare2<-t(apply(check,1,function(r) c(mean(survcheck(0.004,0.50,r[1],r[2],n)), mean(survcheck(0.002,0.05,0,r[2],n)))))
# 
# compare3<-t(apply(check,1,function(r) c(mean(survcheck(0.02,0.10,r[1],r[2],n)), mean(survcheck(0.002,0.05,0,r[2],n)))))
# 
# compare4<-t(apply(check,1,function(r) c(mean(survcheck(0.02,0.50,r[1],r[2],n)), mean(survcheck(0.002,0.05,0,r[2],n)))))
# 
# save(compare1,compare2,compare3,compare4,file="optimalstrategies2")
# 

