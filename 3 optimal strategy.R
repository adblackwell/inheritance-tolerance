#optimal strategy should depend on how much time spent in tolerance vs. how much time spent fighting/free of infection
#should it just depend on the ratio between the two and the ratio of time infected? One might also be more stochastic...

# #simple model
# params<-c(
#   costfight=0.04,
#   costtolerate=0.004,
#   probkill=0.10,
#   probclear=0.01,
#   probimmunity=0.05,
#   reinfectionprob=0.135
# )
#reinfectionprob will be R0*infections/popsize. Actually it'll be 1-(1-R0/popsize)^infected
#so for R0=2.5 and 80% infected 1-(1-2.5/1000)^800 = 0.865
#cumulative sum of mortality risk for each

library(doParallel)
library(parallel)
library(foreach)

#at each point, should be cumulative prob dying - cumulative prob of clearing? 
survcheck<-function(costfight,probkill,probimmunity,reinfectionprob,n,cores){
  #lengths<-NA
  lengths<-foreach(j=1:cores, .combine=c, .packages="foreach") %dopar% {
    foreach(h=1:ceiling((n/cores)), .combine=c) %do% {
      #want to know survival if infected so start each infected
      infected<-1
      immune<-0
      previnfect<-0
      for(i in 1:840){
        if(infected==1) {
          #check if infection clears
          infected <- 1-rbinom(1,1,probkill)
          if(infected==0) {
            #check if gains immunity after clearance
            previnfect<-previnfect+1
            immune <- rbinom(1,1,probimmunity/previnfect)
          }
        }
        if(infected==0 & immune==0) {
          #check if gains infection
          infected <- rbinom(1,1,reinfectionprob)
        }
        if(infected==1) {
          #check if dies from infection
          if(rbinom(1,1,costfight)==1){
            lengthi<-i
            break()
          }
        }
        #if makes it to 70, dies at 70
        if(i==840) lengthi<-i
      }
      lengthi
      #lengths<-c(lengths,lengthi)
    }
  }
  return(lengths)
  }


s1<-survcheck(0.04,0.10,0.05,0.865,1000)
s2<-survcheck(0.004,0.01,0,0.865,1000)

s1<-survcheck(0.04,0.20,0.02,0,1000)
s2<-survcheck(0.01,0.01,0,0,1000)

# 
# #go back and modify my model. If mom was able to be immune, should never go for tolerance strategy...
# immunity<-seq(0,1,0.05)
# reinf<-seq(0,1,0.025)
# check<-expand.grid(immunity,reinf)
# n<-2000
# #mm<-function(nums) c(mean(nums),median(nums))
# cores<-(detectCores()-1)
# cl <- makeCluster(cores)
# registerDoParallel(cl)
# tol<-t(sapply(reinf,function(r) mean(survcheck(0.002,0.05,0,r,n,cores))))
# 
# compare1<-apply(check,1,function(r) mean(survcheck(0.004,0.10,r[1],r[2],n,cores)))
# 
# compare2<-apply(check,1,function(r) mean(survcheck(0.004,0.50,r[1],r[2],n,cores)))
# 
# compare3<-apply(check,1,function(r) mean(survcheck(0.02,0.10,r[1],r[2],n,cores)))
# 
# compare4<-apply(check,1,function(r) mean(survcheck(0.02,0.50,r[1],r[2],n,cores)))
# stopCluster(cl)
# save(tol,compare1,compare2,compare3,compare4,file="optimalstrategies2")
# 
# 
# plot(reinf,compare1[1:(length),1],type="l",ylim=c(0,840))
# lines(reinf,compare1[,2],col="red")

#different version
#go back and modify my model. If mom was able to be immune, should never go for tolerance strategy...
#immunity<-seq(0,1,0.05)
reinf<-c(0,0.01,0.05,0.25)
costs2<-0.05*1:10
clear<-0.01*seq(1,51,length=10)
check2<-expand.grid(costs2,clear)
n<-2000
#mm<-function(nums) c(mean(nums),median(nums))
#costfight,probkill,probimmunity,reinfectionprob,n
cores<-(detectCores()-1)
cl <- makeCluster(cores)
registerDoParallel(cl)

#tol<-t(sapply(reinf,function(r) mean(survcheck(0.002,0.05,0,r,n,cores))))

compare1<-apply(check2,1,function(r) mean(survcheck(r[1],r[2],0,reinf[1],n,cores)))

compare2<-apply(check2,1,function(r) mean(survcheck(r[1],r[2],0,reinf[2],n,cores)))

compare3<-apply(check2,1,function(r) mean(survcheck(r[1],r[2],0,reinf[3],n,cores)))

compare4<-apply(check2,1,function(r) mean(survcheck(r[1],r[2],0,reinf[4],n,cores)))
stopCluster(cl)
save(compare1,compare2,compare3,compare4,file="optimalstrategies3")


load("optimalstrategies3")

library(wesanderson)
colp<-colorRampPalette(c("white",wes_palette("Darjeeling1"))[c(3,1,5)])
#colp<-colorRampPalette(wes_palette("Darjeeling")[c(3,5)])

ratio<-matrix(log(compare1[1]/compare1),nrow=length(costs1),byrow=FALSE)
ratio2<-matrix(log(compare2[1]/compare2),nrow=length(costs1),byrow=FALSE)
ratio3<-matrix(log(compare3[1]/compare3),nrow=length(costs2),byrow=FALSE)
ratio4<-matrix(log(compare4[1]/compare4),nrow=length(costs2),byrow=FALSE)

v<-round(max(abs(quantile(c(ratio,ratio2,ratio3,ratio4),c(0.05,0.95)))),1)
ratio[ratio>v]<-v
ratio[ratio<(v*-1)]<-(v*-1)
ratio2[ratio2>v]<-v
ratio2[ratio2<(v*-1)]<-(v*-1)
ratio3[ratio3>v]<-v
ratio3[ratio3<(v*-1)]<-(v*-1)
ratio4[ratio4>v]<-v
ratio4[ratio4<(v*-1)]<-(v*-1)

ratios<-list(ratio,ratio2,ratio3,ratio4)

pdf("figures/tolerateorresist2.pdf",width=3.43, height=2.75,pointsize=8)
layout(matrix(c(0,0,0,6,7, 12,12,10,1,2,  12,12,10,8,9, 0,11,10,3,4, 0,0,0,5,5),ncol=5,byrow=TRUE),widths=c(0.5,0.5,0.15,1,1),heights=c(0.15,1.5,0.15,1.5,0.15))
mars<-c(2,1.8,0.2,1)
par(mar=mars)
ncols<-100
for(i in 1:length(ratios)){
  plot(0,0,type="n",xlab=NA,ylab=NA,xlim=c(1,10),ylim=c(1,10),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
  .filled.contour(1:10,1:10,t(ratios[[i]]),levels=seq(-v,v,length=(ncols+1)),col=colp(ncols))
  text(1+0.02*10,10-0.02*10,LETTERS[i+1],cex=2,adj=c(0,1),font=2)
  box()
  axis(1,at=c(0.10,0.20,0.30,0.40,0.50)*10/0.51,labels=FALSE)
  axis(2,at=c(2,4,6,8,10),labels=FALSE)
  if(i %in% c(1,3)) {
    axis(2,cex.axis=1,at=c(2,4,6,8,10),labels=c("2x","4x","6x","8x","10x"))
  }
  if(i %in% c(3,4)){
    axis(1,cex.axis=1,at=c(0.10,0.20,0.30,0.40,0.50)*10/0.51,labels=c("10x","20x","30x","40x","50x"))
  }
}
par(mar=c(0,mars[2],0,mars[4]))
frame()
text(0.5,0.5,"Relative Clearance from Resistance",cex=1)
frame()
text(0.5,0.5,"Zero reinfection",cex=1,font=2)
frame()
text(0.5,0.5,"1% Reinfection",cex=1,font=2)
frame()
text(0.5,0.5,"5% Reinfection",cex=1,font=2)
frame()
text(0.5,0.5,"25% Reinfection",cex=1,font=2)

par(mar=c(mars[1],0,mars[3],0))
frame()
text(0.5,0.5,"Relative Mortality Cost of Resistance",srt=90,cex=1)


#11
par(mar=c(mars[1]+0.5,3,mars[3]+0.5,1))
plot(0,0,ylim=c(-v,v),xlim=c(0,1),type="n",xlab=NA,ylab=NA,main=NA,yaxt="n",yaxs="i",xaxs="i",xaxt="n")
ints<-seq(-v,v,length=ncols+1)

for(i in 1:ncols){
  polygon(c(0,0,1,1),c(ints[i],ints[i+1],ints[i+1],ints[i]),col=colp(ncols)[i],border=NA)
}
axis(2,at=c(log(0.2),0,log(5)),labels=c("5x Resistance\nSurvival","Neutral","5x Tolerance\nSurvival"),las=2)
box()

#12
par(mar=c(6,3,0.2,1))
frame()
box()
lines(c(0,1),c(0.8,0),col=colp(ncols)[1],lwd=2)
lines(c(0,1),c(0.6,0.2),col=colp(ncols)[ncols],lwd=2)
text(0,1,LETTERS[1],cex=2,adj=c(0,1),font=2)
mtext("Total Time Infected",1,1.2,cex=0.75)
mtext("(prevalence, transmission,\nclearance, immunity)",1,3.2,cex=0.5)
mtext("Host Survival",2,1.5,cex=0.75)

arrows(0, -0.12, 1, -0.12, length=0.06, xpd = TRUE)
arrows(-0.12, 0, -0.12, 1, length=0.06, xpd = TRUE)
legend("topright",c("Resisting","Tolerating"),col=colp(ncols)[c(1,ncols)],lty=1,lwd=2,bty="o",box.col=NA,bg=NA,inset=c(0.07,0.07),title=expression(bold("Strategy")),cex=0.75)
dev.off()









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

