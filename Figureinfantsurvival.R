#model figure survival
load("allpopsinheritcompare")
tiff("survplot.tif",width=2000,height=800,res=400,compression="lzw",pointsize=8)
survplot(allpops[1:3],cols=cols)
dev.off()

dieexp=0.0015
pdf("figures/survivalfig.pdf",width=3.43, height=3,pointsize=10)
library(wesanderson)
cols<-c("black",wes_palette("Darjeeling1")[c(3,2,1,4,5)])
par(mfrow=c(1,1),mar=c(5,4,1,1))
survplot2(allpops[c(4,1,2,3)],n1=20*12+1,n2=120*12,cols=c("black",cols[3:5]),labs=c("No Infection","No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))
dev.off()

survplot2(allpops[c(1,2,4)],n1=1,n2=50*12,cols=cols[3:5],labs=c("No Tolerance (A-B)","No Maternal Effect (C-D)","Maternal Effect (E-F)"))

survplot<-function(allpops,n=length(allpops[[1]]$pops),ymax=4000,cols=c("black","blue","green","red","orange","purple")){
  require(survival)
  par(mfrow=c(1,3),mar=c(5,4,1,1))
  bottommar<-1
  ii<-1
  sss<-list()
  for(pops1 in allpops){
    pops<-pops1$pops[1:n]
    infections<-pops1$infections
    params<-pops1$params
    infected<-Reduce(rbind,lapply(pops,function(x) c(sum(x$infected==0 & x$immune==0),sum(x$immune==1),sum(x$infected>0))))
    popsize<-infected[,1]+infected[,2]+infected[,3]
    plot(1,1,ylim=c(0,ymax),xlim=c(1,n),type="n",ylab="Individuals (x 1000)",xaxt="n",yaxt="n",xlab=NA)
    lines(1:nrow(infected),popsize,lty=2,lwd=2,col=cols[1])
    lines(1:nrow(infected),infected[,1],lty=1,lwd=2,col=cols[2])
    lines(1:nrow(infected),infected[,2],lty=1,lwd=2,col=cols[3])
    lines(1:nrow(infected),infected[,3],lty=1,lwd=2,col=cols[4])
    ax<-axTicks(1,axp=c(0,n/12,10))
    #ay<-axTicks(2,axp=c(0,max(popsize)/1000,5))
    axis(2,at=seq(0,10,0.5)*1000,labels=seq(0,10,0.5))
    mtext("Years",side=1,line=3)
    axis(1,at=ax*12,labels=ax)

    legend("topleft",c("All","Susceptible","Immune","Infected"),col=cols[1:4],lty=c(2,1,1,1),lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0.05),title=expression(bold("Status")))
    text(0,ymax,LETTERS[ii],cex=2,adj=c(0,1),font=2)

    
    vitalrate2<-Reduce(rbind,lapply(1:n,function(x) expand.grid(month=x,agedeath=pops[[x]]$age[pops[[x]]$die==1],die=1,mod=1)))
    vitalrate2<-rbind(vitalrate2,expand.grid(month=n,agedeath=pops[[n]]$age[pops[[n]]$die==0],die=0,mod=1))
    vitalrate2$monthborn<-round(vitalrate2$month-vitalrate2$agedeath*12)
    vitalrate2<-vitalrate2[vitalrate2$monthborn>1,]
    s<-survfit(Surv(agedeath,die)~1,data=vitalrate2)
    sss[[ii]]<-s
    ii<-ii+1
  }
    
  plot(1,1,ylim=c(0,1),xlim=c(1,90),type="n",ylab="Survival (%)",xlab=NA)
  ii<-5
  for(s in sss){
    lines(s,col=cols[ii])
    ii<-ii+1
  }
  mtext("Age",side=1,line=3)
  legend("topright",c("No Maternal Effect","Maternal Effect"),col=cols[5:6],lty=c(1,1),lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0.05),title=expression(bold("Maternal Effects")))
  text(0,1,"C",cex=2,adj=c(0,1),font=2)
}


survplot2<-function(allpops,n1=1,n2=length(allpops[[1]]$pops),cols=rainbow(length(allpops)),labs=seq(1:length(allpops))){
  require(survival)
  ii<-1
  sss<-list()
  for(pops1 in allpops){
    pops <- pops1$pops[n1:n2]
    n <- (n2-n1)+1
    vitalrate2<-Reduce(rbind,lapply(1:n,function(x) expand.grid(month=x,agedeath=pops[[x]]$age[pops[[x]]$die==1],die=1,mod=1)))
    vitalrate2<-rbind(vitalrate2,expand.grid(month=n,agedeath=pops[[n]]$age[pops[[n]]$die==0],die=0,mod=1))
    vitalrate2$monthborn<-round(vitalrate2$month-vitalrate2$agedeath*12)
    #vitalrate2<-vitalrate2[vitalrate2$monthborn>1,]
    s<-survfit(Surv(agedeath,die)~1,data=vitalrate2)
    sss[[ii]]<-s
    ii<-ii+1
  }
  
  plot(1,1,ylim=c(0,1),xlim=c(1,100),type="n",ylab="Survival (%)",xlab=NA)
  ii<-1
  for(s in sss){
    lines(s,conf.int=FALSE,col=cols[ii])
    ii<-ii+1
  }
  # #add line for no infection
  # if(noinfection){
  #   dierisk <- abs((1-((960-seq(0,1200,1)+240)/960)^0.002))
  #   lines(seq(0,1200,1)/12, cumprod(1-dierisk))
  #   labs<-c("No Infection",labs)
  #   cols<-c("black",cols)
  # }
  mtext("Age",side=1,line=3)
  legend("topright",legend=labs,col=cols,lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0.05),cex=0.75)
}

pops <- pops1$pops[n1:n2]
n <- (n2-n1)+1
p2<-lapply(pops,)
vitalrate2<-Reduce(rbind,lapply(1:n,function(x) expand.grid(month=x,agegroup=pops[[x]]$age[pops[[x]]$die==1],die=1,mod=1)))

table(floor(pops[[800]]$age/10)*10,ifelse(pops[[800]]$immune==1,3,pops[[800]]$strategy))
