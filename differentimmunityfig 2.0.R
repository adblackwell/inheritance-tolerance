load(file="allpops4.2")
library(wesanderson)
cols<-c("black",wes_palette("Darjeeling")[c(3,2,1,4,5)])
pdf("figures/immunefig.pdf",width=7, height=4,pointsize=12)
strategyreport4(allpops[1:3],n=50*12,cols=cols,popmax=2000)
dev.off()

strategyreport4(allpops2,n=100*12,cols=cols,popmax=4000)

tiff("zoomfig.tif",width=1000,height=1000,res=400,compression="lzw",pointsize=8)
srzoom(allpops[[4]],240,cols=cols)
dev.off()

tiff("singlefig.tif",width=900,height=1500,res=400,compression="lzw",pointsize=9)
strategyreportsingle(allpops[[1]],cols=cols)
dev.off()

strategyreport4<-function(allpops,n=length(allpops[[1]]$pops),cols=c("black","blue","green","red","orange","purple"),popmax=NA){
  
  layout(rbind(matrix(seq(1,length(allpops)*2,1),ncol=length(allpops),byrow=FALSE),rep(0,length(allpops))),widths=rep(4,length(allpops)),heights=c(1,1,0.3))
  bottommar<-1
  ii<-1
  for(pops1 in allpops){
    pops<-pops1$pops[1:n]
    infections<-pops1$infections
    params<-pops1$params
    
    par(mar=c(bottommar,4.5,1,1))
    infected<-Reduce(rbind,lapply(pops,function(x) c(sum(x$infected==0 & x$immune==0),sum(x$immune==1),sum(x$infected>0))))
    popsize<-infected[,1]+infected[,2]+infected[,3]
    if(is.na(popmax)) popmax<-max(popsize)
    plot(1,1,ylim=c(0,popmax),xlim=c(1,n),type="n",ylab="Individuals (x 1000)",xaxt="n",yaxt="n")
    lines(1:nrow(infected),popsize,lty=2,lwd=2,col=cols[1])
    lines(1:nrow(infected),infected[,1],lty=1,lwd=2,col=cols[2])
    lines(1:nrow(infected),infected[,2],lty=1,lwd=2,col=cols[3])
    lines(1:nrow(infected),infected[,3],lty=1,lwd=2,col=cols[4])
    ax<-axTicks(1,axp=c(0,n/12,10))
    #ay<-axTicks(2,axp=c(0,max(popsize)/1000,5))
    axis(2,at=seq(0,10,0.5)*1000,labels=seq(0,10,0.5))
    axis(1,at=ax*12,labels=NA)
    
    #par(mar=c(bottommar,0,0.2,1))
    #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    legend("topleft",c("All","Susceptible","Immune","Infected"),col=cols[1:4],lty=c(2,1,1,1),lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0.05),title=expression(bold("Status")))
    mtext(LETTERS[ii],side=3,line=-1,at=-200,cex=1.5)
    ii<-ii+1
    
    par(mar=c(bottommar,4.5,1,1))
    strategy<-Reduce(rbind,lapply(pops,function(x) c(sum(x$strategy==1 & x$immune==0),sum(x$immune==1),sum(x$strategy==2 & x$immune==0))))
    plot(1,1,ylim=c(0,1),xlim=c(1,n),type="n",ylab="Proportion",xaxt="n")
    lines(1:nrow(strategy),strategy[,1]/popsize,lty=1,lwd=2,col=cols[1])
    lines(1:nrow(strategy),strategy[,2]/popsize,lty=1,lwd=2,col=cols[3])
    lines(1:nrow(strategy),strategy[,3]/popsize,lty=1,lwd=2,col=cols[5])
    mtext("Years",side=1,line=3)
    axis(1,at=ax*12,labels=ax)

    #par(mar=c(bottommar,0,0.2,1))
    #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    legend("right",c("Resist","Immune","Tolerate"),col=cols[c(1,3,5)],lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0),title=expression(bold("Strategy")))
    mtext(LETTERS[ii],side=3,line=-1,at=-200,cex=1.5)
    ii<-ii+1
  }
}


srzoom<-function(allpops,n=length(allpops$pops),cols=c("black","blue","green","red","orange","purple"),statlim=c(0.2,1)){
    pops<-allpops$pops[1:n]
    par(mar=c(4,4,1,1))
    infected<-Reduce(rbind,lapply(pops,function(x) c(sum(x$infected==0 & x$immune==0),sum(x$immune==1),sum(x$infected>0))))
    popsize<-infected[,1]+infected[,2]+infected[,3]
    plot(1,1,ylim=c(0,max(popsize)),xlim=c(1,n),type="n",ylab="Individuals (x 1000)",xaxt="n",yaxt="n")
    lines(1:nrow(infected),popsize,lty=2,lwd=2,col=cols[1])
    lines(1:nrow(infected),infected[,1],lty=1,lwd=2,col=cols[2])
    lines(1:nrow(infected),infected[,2],lty=1,lwd=2,col=cols[3])
    lines(1:nrow(infected),infected[,3],lty=1,lwd=2,col=cols[4])
    ax<-axTicks(1,axp=c(0,n/12,10))
    #ay<-axTicks(2,axp=c(0,max(popsize)/1000,5))
    axis(2,at=seq(0,10,0.5)*1000,labels=seq(0,10,0.5))
    mtext("Years",side=1,line=3)
    axis(1,at=ax*12,labels=ax)
}


strategyreportsingle<-function(pops1,n=length(pops1$pops),cols=c("black","blue","green","red","orange","purple"),statlim=c(0.2,1)){
  
  layout(matrix(c(1,2,3,4,0,0),ncol=2,byrow=TRUE),widths=c(0.5,4),heights=c(1,1,0.3))
  bottommar<-1

    pops<-pops1$pops[1:n]
    infections<-pops1$infections
    params<-pops1$params
    
    par(mar=c(bottommar,1,0.2,0))
    
    frame()
    text(0,1,"A",cex=2,adj=c(0,1),font=2)
    
    par(mar=c(bottommar,4,1,1))
    infected<-Reduce(rbind,lapply(pops,function(x) c(sum(x$infected==0 & x$immune==0),sum(x$immune==1),sum(x$infected>0))))
    popsize<-infected[,1]+infected[,2]+infected[,3]
    plot(1,1,ylim=c(0,max(popsize)),xlim=c(1,n),type="n",ylab="Individuals (x 1000)",xaxt="n",yaxt="n")
    lines(1:nrow(infected),popsize,lty=2,lwd=2,col=cols[1])
    lines(1:nrow(infected),infected[,1],lty=1,lwd=2,col=cols[2])
    lines(1:nrow(infected),infected[,2],lty=1,lwd=2,col=cols[3])
    lines(1:nrow(infected),infected[,3],lty=1,lwd=2,col=cols[4])
    ax<-axTicks(1,axp=c(0,n/12,10))
    #ay<-axTicks(2,axp=c(0,max(popsize)/1000,5))
    axis(2,at=seq(0,10,0.5)*1000,labels=seq(0,10,0.5))
    axis(1,at=ax*12,labels=NA)
    
    #par(mar=c(bottommar,0,0.2,1))
    #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    legend("topleft",c("All","Susceptible","Immune","Infected"),col=cols[1:4],lty=c(2,1,1,1),lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0.05),title=expression(bold("Status")))
    
    par(mar=c(bottommar,1,0.2,0))
    
    frame()
    text(0,1,"B",cex=2,adj=c(0,1),font=2)
    
    par(mar=c(bottommar,4,1,1))
    strategy<-Reduce(rbind,lapply(pops,function(x) c(sum(x$strategy==1 & x$immune==0),sum(x$immune==1),sum(x$strategy==2 & x$immune==0))))
    plot(1,1,ylim=c(0,1),xlim=c(1,n),type="n",ylab="Proportion",xaxt="n")
    lines(1:nrow(strategy),strategy[,1]/popsize,lty=1,lwd=2,col=cols[1])
    lines(1:nrow(strategy),strategy[,2]/popsize,lty=1,lwd=2,col=cols[3])
    lines(1:nrow(strategy),strategy[,3]/popsize,lty=1,lwd=2,col=cols[5])
      mtext("Years",side=1,line=3)
      axis(1,at=ax*12,labels=ax)

    #par(mar=c(bottommar,0,0.2,1))
    #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    legend("right",c("Resist","Immune","Tolerate"),col=cols[c(1,3,5)],lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0),title=expression(bold("Strategy")))
}
