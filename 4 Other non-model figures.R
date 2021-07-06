
#<<trade-offfig>>= (Fig 1)
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


#Figure S1 showing profiles for three pathogen sets
paramsflu<-makeparams(months=120,h2current=1,h2past=0.1,PI=20,infectedT0=10,Pmr=0.10,Pmt=0.05,Pct=0,Pcr=0.8,Pim=0.90)
params2<-makeparams(h2current=0,h2past=0,Pim=0.05,noflip=TRUE)
params2H<-makeparams(h2current=0,h2past=0,Pim=0,noflip=TRUE,PI=0.1,Pmr=0.0005,Pmt=0.0001,Pcr=0.04,Pct=0.02,K=5000,dieexp=0.004,reprodrate=0.010)

pdf("figures/parasiteprofiles.pdf",width=7, height=3.5,pointsize=12)
layout(matrix(c(8,9,0, 1,2,0, 3,4,7, 5,6,0, 10,11,0),nrow=3,byrow=FALSE),widths=c(0.75,4,4,4,2.5),heights=c(1,1,0.1))
paramscheck(paramsflu,xmax1=4,xmax2=4,nolayout=TRUE,legend=NA)
paramscheck(params2,xmax1=10,xmax2=10,nolayout=TRUE,legend=NA,letters=c("C","D"))
paramscheck(params2H,xmax1=10,xmax2=10,nolayout=TRUE,legend=NA,letters=c("E","F"))
par(mar=c(0,5,0,0.5))
frame()
text(0.5,0.7,"Years",cex=2,adj=c(0.5,0.5))
par(mar=c(3,0,0.5,0))
frame()
text(0.5,0.5,expression(bold("Resist")),cex=2,adj=c(0.5,0.5),srt=90)
frame()
text(0.5,0.5,expression(bold("Tolerate")),cex=2,adj=c(0.5,0.5),srt=90)


par(mar=c(3,0,1,0.2))
plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
legend(0,1,legend=c("Survival","Clearance"),col=cols[c(1,4)],lty=c(1,1),bty="n",lwd=2)

par(mar=c(3,0,1,0.2))
plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
legend(0,1,legend=c("Survival","Clearance"),col=cols[c(3,2)],lty=c(1,1),bty="n",lwd=2)

dev.off()


#Figure 4 Model description
library(wesanderson)
pdf("figures/modeldiagram.pdf",width=3.43, height=2.2,pointsize=7)
cols<-c("lightgrey",wes_palette("Darjeeling1")[c(3,2,1,4,5)])

par(mar=c(0,0,0,0))
plot(0,0,type="n",xlim=c(0.1,1),ylim=c(0.15,0.78),xaxt="n",yaxt="n",bty="n")
#three columns
boxwidth<-0.2
boxheight<-0.12
boxx<-c(0.2,0.6,0.6,0.3,0.9,0.9)
boxy<-c(0.4,0.4,0.65,0.65,0.275,0.525)
boxtext<-c("Susceptible\n(Resist)","Infected\nResisting","Infected\nTolerating","Susceptible\n(Tolerate)","Immune","Deceased")
boxcols<-c(3,2,2,3,5,1)
left<-boxx-boxwidth/2
bottom<-boxy-boxheight/2
right<-boxx+boxwidth/2
top<-boxy+boxheight/2



#last set of xs and ys will be used as arrow
straight<-matrix(c(right[1],boxy[1],left[2],boxy[2],
                   right[2],boxy[2],left[5],boxy[5],
                   right[2],boxy[2],left[6],boxy[6],
                   right[3],boxy[3],left[6],boxy[6],
                   boxx[2],top[2],boxx[3],bottom[3],
                   right[4],boxy[4],left[3],boxy[3]
),ncol=4,byrow=TRUE)
straightcols<-c(3,2,1,1,2,3)
straightlwds<-c(3,3,3,1,3,3)
library(shape)
Arrows(straight[,1],straight[,2],straight[,3],straight[,4],arr.adj = 1,arr.type="triangle",lwd=straightlwds,col=cols[straightcols],arr.length = 0.2,arr.lwd=)

off1<-0.05

alllines<-list(matrix(c(boxx[2],bottom[2], 
                        boxx[2],bottom[2]-off1,
                        boxx[1],bottom[2]-off1,
                        boxx[1],bottom[1]
),ncol=2,byrow=TRUE),
matrix(c(boxx[3],top[3], 
         boxx[3],top[3]+off1,
         boxx[4],top[3]+off1,
         boxx[4],top[4]
),ncol=2,byrow=TRUE),
matrix(c(boxx[5],bottom[5], 
         boxx[5],bottom[5]-off1,
         boxx[1]-boxwidth/4,bottom[5]-off1,
         boxx[1]-boxwidth/4,bottom[1]
),ncol=2,byrow=TRUE),
matrix(c(boxx[4],bottom[4],
         boxx[4],top[1]+off1,
         boxx[1]-boxwidth/4,top[1]+off1,
         boxx[1]-boxwidth/4,top[1]
),ncol=2,byrow=TRUE)
)

arrs<-lapply(alllines,function(x) c(x[nrow(x)-1,],x[nrow(x),]))
ltys<-c(1,1,2,2)
lcols<-c(2,2,5,3)
lwds<-c(3,1,1,1)
for (i in 1:length(alllines)){               
  lines(alllines[[i]],lty=ltys[i],col=cols[lcols[i]],lwd=lwds[i])
  Arrows(arrs[[i]][1],arrs[[i]][2],arrs[[i]][3],arrs[[i]][4],arr.adj = 1,arr.type="triangle",col=cols[lcols[i]],lwd=lwds[i],arr.length = 0.2,arr.lwd=2,segment=FALSE)  
  
}

rect(left, bottom, right, top,col=cols[boxcols])
text(boxx,boxy,boxtext)
labels<-c(expression('P'['I']),expression('P'['IM']),expression('P'['MR']),expression('P'['MT']),expression('P'['R->T']),expression('P'['I']))
labels2<-c(expression('P'['CR']),expression('P'['CT']),expression('E'['w']~'='~'0'),expression('E'['w']~'='~'0'))

placement<-t(apply(straight,1,function(x) c(mean(x[c(1,3)]),mean(x[c(2,4)]))))
placement2<-Reduce(rbind,lapply(alllines,function(x) c(mean(x[2:3,1]),mean(x[2:3,2]))))


text(placement[,1],placement[,2],labels)
text(placement2[c(1,2),1],placement2[c(1,2),2],labels2[c(1,2)])
text(placement2[c(3,4),1],placement2[c(3,4),2],labels2[c(3,4)],pos=3)

text(0.9,0.75,expression('P'['MR'] ~ '>' ~ 'P'['MT']))
text(0.9,0.7,expression('P'['CR'] ~ '>' ~ 'P'['CT']))
text(0.9,0.65,expression('P'['R->T']~'~'~ 'f(E'['w']~')'))
dev.off()