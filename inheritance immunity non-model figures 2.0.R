
#<<trade-offfig>>=
pdf("figures/trade-off fig.pdf",width=3.43, height=3.2,pointsize=10)
library(wesanderson)
cols<-c("black",wes_palette("Darjeeling")[c(3,2,1,4,5)])
par(mar=c(0,0,0,0))
plot(0,0,type="n",ylim=c(0,0.9),xlim=c(0,1),bty="n",xaxt="n",yaxt="n")
res<-1000
xs<-seq(0,1,length=res)
frac<-0.9
mid<-0.6
left<-0.15
top<-(xs-0.5)^2
top<-top/max(top)*(frac-mid)+mid
bottom<-(1-top)

middle<-(xs-0.5)^2
middle<-middle/max(middle)/2+0.5
middle[(res/2):res]<- 1-middle[(res/2):res]
middle<-(middle)*(top[res*left]-(1-top[res*left]))+(1-top[res*left])
xs2<-(xs-0.5)*(1-left*2)+0.5

px<-c(xs[1:(res-res*left)],rev(xs2)[2:(res-1)],rev(xs[1:res*left]))
py<-c(bottom[1:(res-res*left)],rev(middle)[2:(res-1)],rev(top[1:res*left]))
polygon(px,py,col=cols[6],border="black")
polygon(1-px,1-py,col=cols[3],border="black")


a1<-0.0
a2<-0.05
arrows(1, a1, 0, a1, length=0.1, xpd = TRUE,col=cols[6],lwd=2)

sw<-strwidth("Parasite Replication",font=2)
sh<-strheight("Parasite Replication",font=2)
polygon(c(0.5-sw/2,0.5-sw/2,0.5+sw/2,0.5+sw/2),c(a1-sh/2,a1+sh/2,a1+sh/2,a1-sh/2),border=NA,col="white")
text(0.5,a1,"Parasite Replication",col=cols[6],font=2)

arrows(0, a2, 1, a2, length=0.1, xpd = TRUE,col=cols[3],lwd=2)
sw<-strwidth("Immune Response",font=2)
sh<-strheight("Immune Response",font=2)
polygon(c(0.5-sw/2,0.5-sw/2,0.5+sw/2,0.5+sw/2),c(a2-sh/2,a2+sh/2,a2+sh/2,a2-sh/2),border=NA,col="white")
text(0.5,a2,"Immune Response",col=cols[3],font=2)

text(left-0.02,0.5,"Parasite\n pathology",font=2)
text(1-left+0.02,0.5,"Immuno-\n pathology",font=2)
tol1<-0.37
tol2<-0.63
y1<-top[which.min(abs(xs-(1-left)))]+0.05
y4<-top[which.min(abs(xs-(tol2)))]+0.05
y2<-bottom[which.min(abs(xs-(tol2)))]-0.05
y3<-bottom[which.min(abs(xs-(tol1)))]-0.05
lines(c(1-left,1-left),c(y1+0.07,y1))
lines(c(tol2,tol2),c(y1+0.07,y4))
lines(c(tol2,1-left),c(y1+0.07,y1+0.07))

text((tol2+(1-left))/2,y1+0.1,"Resistance",font=2)

lines(c(tol2,tol2),c(y2-0.1,y2))
lines(c(tol1,tol1),c(y2-0.1,y3))
lines(c(tol1,tol2),c(y2-0.1,y2-0.1))
text(0.5,y2-0.13,"Tolerance",font=2)

dev.off()
#@



#<<modelfig2, echo=FALSE, fig.pos="h", fig.height=3.5, fig.width=3.5, fig.cap="Host survival for each strategy depends on the time an individual is infected. Expected time infected is determined by prevalence and transmission rate, but also by effectiveness of clearance and probability of immunity.", dev.args=list(pointsize=10)>>=
#@

#<<modelfigtol, echo=FALSE, fig.pos="h", fig.height=6, fig.width=6, fig.cap="", dev.args=list(pointsize=10)>>=
  
load("optimalstrategies2")

library(wesanderson)
colp<-colorRampPalette(c("white",wes_palette("Darjeeling"))[c(3,1,5)])
#colp<-colorRampPalette(wes_palette("Darjeeling")[c(3,5)])
immunity<-seq(0,1,0.05)
reinf<-seq(0,1,0.025)
check<-expand.grid(immunity,reinf)
n<-2000
ratio<-matrix(log(rep(tol,each=length(immunity))/compare1),nrow=length(immunity),byrow=FALSE)
ratio2<-matrix(log(rep(tol,each=length(immunity))/compare2),nrow=length(immunity),byrow=FALSE)
ratio3<-matrix(log(rep(tol,each=length(immunity))/compare3),nrow=length(immunity),byrow=FALSE)
ratio4<-matrix(log(rep(tol,each=length(immunity))/compare4),nrow=length(immunity),byrow=FALSE)

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

pdf("figures/tolerateorresist.pdf",width=3.43, height=3,pointsize=8)
layout(matrix(c(0,0,6,7,0, 12,8,1,2,9, 11,8,3,4,10, 0,0,5,5,0),ncol=5,byrow=TRUE),widths=c(1,0.15,1,1,0.15),heights=c(0.15,1.5,1.5,0.15))
mars<-c(2,1.8,0.2,0.2)
par(mar=mars)
ncols<-100
for(i in 1:length(ratios)){
  plot(0,0,type="n",xlab=NA,ylab=NA,xlim=range(reinf),ylim=range(immunity),xaxs="i",yaxs="i",xaxt="n",yaxt="n")
  .filled.contour(reinf,immunity,t(ratios[[i]]),levels=seq(-v,v,length=(ncols+1)),col=colp(ncols))
  text(0+0.02*max(reinf),max(immunity)-0.02*max(immunity),LETTERS[i+1],cex=2,adj=c(0,1),font=2)
  box()
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  if(i %in% c(1,3)) {
    axis(2,cex.axis=1)
  }
  if(i %in% c(3,4)){
    axis(1,cex.axis=1)
  }
}
par(mar=c(0,mars[2],0,mars[4]))
frame()
text(0.5,0.5,"Monthly Probability of Repeat Exposure",cex=1)
frame()
text(0.5,0.5,"Clearance x2",cex=1.2,font=2)
frame()
text(0.5,0.5,"Clearance x10",cex=1.2,font=2)

par(mar=c(mars[1],0,mars[3],0))
frame()
text(0.5,0.5,"Chance of Immunity",srt=90,cex=1)
frame()
text(0.5,0.5,"Cost x2",cex=1.2,font=2,srt= -90)
frame()
text(0.5,0.5,"Cost x10",cex=1.2,font=2,srt= -90)

#11
par(mar=c(mars[1]+1,10,mars[3]+1,1))
plot(0,0,ylim=c(-v,v),xlim=c(0,1),type="n",xlab=NA,ylab=NA,main=NA,yaxt="n",yaxs="i",xaxs="i",xaxt="n")
ints<-seq(-v,v,length=ncols+1)

for(i in 1:ncols){
  polygon(c(0,0,1,1),c(ints[i],ints[i+1],ints[i+1],ints[i]),col=colp(ncols)[i],border=NA)
}
axis(2,at=c(log(0.5),0,log(2)),labels=c("2x Resistance\nSurvival","Neutral","2x Tolerance\nSurvival"),las=2)
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
#@


#<<fig1, echo=FALSE, fig.pos="h", fig.height=4, fig.width=7, fig.cap="Host survival and clearance of two pathogens (A/C and B/D) when resisted (A and B) or tolerated (C and D). The first pathogen (A/C) has a high mortality from immunopathology when fought (A), but low pathology when tolerated (C). Resisting does not affect clearance (A and C). The second pathogen has a high mortality under either condition (B and D), but can be cleared quickly by fighting (B).", dev.args=list(pointsize=10),cache=TRUE>>=
  ## ---- fig1 ----
xmax1<-4
xmax2<-4
costfight<-0.04
costtolerate<-0.004
probkill<-0.02
probclear=0.02

costfight2<-0.10
costtolerate2<-0.10
probkill2<-0.4
probclear2<-0.02

#graph survival with fight and tolerate
months<-seq(0,70*12,0.1)
survivalfight<-sapply(months,function(x) (1-costfight)^x)
mediansurvivalfight<-months[which.min(abs(0.50-survivalfight))]
survivaltolerate<-sapply(months,function(x) (1-costtolerate)^x)
mediansurvivaltolerate<-months[which.min(abs(0.50-survivaltolerate))]

#graph clearance with and without fight
clearfight<-1-sapply(months,function(x) (1-probkill)^x)
medianclearfight<-months[which.min(abs(0.50-clearfight))]
cleartolerate<-1-sapply(months,function(x) (1-probclear)^x)
mediancleartolerate<-months[which.min(abs(0.50-cleartolerate))]

#graph survival with fight and tolerate 2
survivalfight2<-sapply(months,function(x) (1-costfight2)^x)
mediansurvivalfight2<-months[which.min(abs(0.50-survivalfight2))]
survivaltolerate2<-sapply(months,function(x) (1-costtolerate2)^x)
mediansurvivaltolerate2<-months[which.min(abs(0.50-survivaltolerate2))]

#graph clearance with and without fight 2
clearfight2<-1-sapply(months,function(x) (1-probkill2)^x)
medianclearfight2<-months[which.min(abs(0.50-clearfight2))]
cleartolerate2<-1-sapply(months,function(x) (1-probclear2)^x)
mediancleartolerate2<-months[which.min(abs(0.50-cleartolerate2))]


layout(matrix(c(0,0,0,0,7,1,3,4,8,2,5,6,0,0,0,0),ncol=4,byrow=TRUE),widths=c(0.3,4,4,1.2),heights=c(0.2,1,1,0.1))
par(mar=c(3,5,1,1))
plot(1,1,ylim=c(0,1),xlim=c(0,xmax1),type="n",ylab=NA,xlab=NA)
lines(months/12,survivalfight,col=cols[1],lwd=2)
lines(rep(mediansurvivalfight,2)/12,c(0,0.5),col=cols[1],lty=2,lwd=2)

lines(months/12,clearfight,col=cols[4],lwd=2)
lines(rep(medianclearfight,2)/12,c(0,0.5),col=cols[4],lty=2,lwd=2)
mtext(expression(bold("A")),side=3,line=-1,at=-1,cex=1.5)
mtext(expression(bold("Pathogen 1")),side=3,line=0.5)
mtext(expression(bold("Fight")),side=2,line=5,cex=1.5)
mtext("Proportion",side=2,line=3)

par(mar=c(3,5,1,1))
plot(1,1,ylim=c(0,1),xlim=c(0,xmax2),type="n",ylab=NA,xlab="Years")
lines(months/12,survivaltolerate,col=cols[3],lwd=2)
lines(rep(mediansurvivaltolerate,2)/12,c(0,0.5),col=cols[3],lty=2,lwd=2)
lines(months/12,cleartolerate,col=cols[2],lwd=2)
lines(rep(mediancleartolerate,2)/12,c(0,0.5),col=cols[2],lty=2,lwd=2)
mtext(expression(bold("C")),side=3,line=-1,at=-1,cex=1.5)
mtext("Years",side=1,line=3)
mtext("Proportion",side=2,line=3)
mtext(expression(bold("Tolerate")),side=2,line=5,cex=1.5)

#pathogen 2
par(mar=c(3,5,1,1))
plot(1,1,ylim=c(0,1),xlim=c(0,xmax1),type="n",ylab=NA,xlab=NA)
lines(months/12,survivalfight2,col=cols[1],lwd=2)
lines(rep(mediansurvivalfight2,2)/12,c(0,0.5),col=cols[1],lty=2,lwd=2)

lines(months/12,clearfight2,col=cols[4],lwd=2)
lines(rep(medianclearfight2,2)/12,c(0,0.5),col=cols[4],lty=2,lwd=2)
mtext(expression(bold("B")),side=3,line=-1,at=-1,cex=1.5)
mtext(expression(bold("Pathogen 2")),side=3,line=0.5)

par(mar=c(3,0,1,0.2))
plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
legend(0,1,legend=c("Survival","Clearance"),col=cols[c(1,4)],lty=c(1,1),bty="n",lwd=2)

par(mar=c(3,5,1,1))
plot(1,1,ylim=c(0,1),xlim=c(0,xmax2),type="n",ylab=NA,xlab="Years")
lines(months/12,survivaltolerate2,col=cols[3],lwd=2)
lines(rep(mediansurvivaltolerate2,2)/12,c(0,0.5),col=cols[3],lty=2,lwd=2)
lines(months/12,cleartolerate2,col=cols[2],lwd=2)
lines(rep(mediancleartolerate2,2)/12,c(0,0.5),col=cols[2],lty=2,lwd=2)
mtext(expression(bold("D")),side=3,line=-1,at=-1,cex=1.5)
mtext("Years",side=1,line=3)

par(mar=c(3,0,1,0.2))
plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
legend(0,1,legend=c("Survival","Clearance"),col=cols[c(3,2)],lty=c(1,1),bty="n",lwd=2)
#@

