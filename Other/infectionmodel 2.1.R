##changelog
##2.2 Added Bayesian flip function
## ---- modelfunctions ----

runmodel<-function(params,ageburn=2000,model=NA){
  with(params,{    
    #starting population
    pop<-data.frame(age=abs(rnorm(n0,0,40)),strategy=startingstrategy,infected=0,immune=0,previousinfections=0,PIsum=0,PIweighted=0,PNIweighted=0,dieinfect=0,dieage=0,die=0,reprod=0,RS=0,R0=0,RR=0)
    #baseline immune history for uninfected population
    pop$PNIweighted<-sapply(pop$age,function(x) sum(PNIdecline^(seq(1,x*12+priorMon))))
    #burn in age structure
    for(i in 1:ageburn){
      pop$age<-pop$age+(1/12)
      pop$PNIweighted<-pop$PNIweighted*PNIdecline+1
      #count time since last reprod, know if eligible to reproduce again
      pop$reprod[pop$reprod>=1]<-pop$reprod[pop$reprod>=1]+1
      pop$reprod[pop$reprod>9]<-0
      #reproduce
      pop$reprod[pop$reprod==0 & pop$age>=18 & pop$age<=45] <- rbinom(length(pop$reprod[pop$reprod==0 & pop$age>=18 & pop$age<=45]),1,reprodrate)
      pop$RS[pop$reprod==1]<-pop$RS[pop$reprod==1]+1
      if(sum(pop$reprod==1)>0){
        babies<-data.frame(age=rep(0,sum(pop$reprod==1)),strategy=startingstrategy,infected=0,immune=0,previousinfections=0,PIsum=0,PIweighted=0,PNIweighted=priorMon,dieinfect=0,dieage=0,die=0,reprod=0,RS=0,R0=0,RR=0)
        #baby gets mom's strategy perfectly for this start, just in case I mix them later and don't want proportion to change with reproduction
        babies$strategy<-pop$strategy[pop$reprod==1]
        pop<-rbind(pop,babies)
      }
      dierisk <- (1-((100-pop$age)/100)^dieexp)
      dierisk[pop$age>100]<-1
      pop$dieage<-rbinom(length(pop$age),1,dierisk)
      #pop$dieage<-as.numeric(pop$age>70)
      pop$die<-as.numeric(pop$dieage==1)
      pop<-pop[pop$die==0,]
      #keep sample at starting size
      if(nrow(pop)>n0) pop<-pop[sample(nrow(pop),n0),]
    }
    #assign ID numbers
    pop$ID<-as.character(seq(1,nrow(pop),1))
    
    message(model, " Finished Age Structure Burn <br>")

    
    #infect starting number
    pop$infected[sample(nrow(pop),infectedT0)]<-1
    
    pops<-list(pop)
    
    #track infections by adding to this whenever an infection ends
    if(trackinfections) infections<-data.frame(month=1,length=NA,R0=NA,RR=NA,immunekilled=NA,strategy=NA,hostdeath=NA,dieinfect=NA)
    
    for(i in 2:months){
      #age
      pop$age<-pop$age+(1/12)
      #cumulative past history gets downweighted each round, by same decline as loss of tolerance
      pop$PNIweighted<-pop$PNIweighted*PNIdecline
      pop$PIweighted<-pop$PIweighted*PIdecline
      pop$PIweighted[pop$PIweighted<(PIdecline^2)]<-0
      
      #chance to lose tolerance/immunity (immunological memory) if not infected and PIweighted reaches zero
      pop$immune[pop$infected==0 & pop$PIweighted==0]<-0
      pop$strategy[pop$infected==0 & pop$PIweighted==0]<-1
      
      #check for strategy switching, can only switch if infected (do this after infection so infants can switch right after infected)
      # define grid
      p_grid <- seq( from=0 , to=1 , length.out=1000 )
      uniprior <- rep( 1 , 100 )
      #posterior is based on weighted months infected x a multiplier (since months infected should update prior more than unifected), and weighted months uninfected. Weighting is to simulate gradual decline in immunological memory
      #posterior is the liklihood that an infection is endemic and cannot or should not be resisted
      
      if(noflip) flipprob<-0 
      else flipprob <- apply(pop[pop$infected>0,c("infected","PIweighted","PNIweighted")],1,function(pp){
          infmonths<-round(pp[1]+pp[2])*infectionweight
          notinfmonths<-round(pp[3])
          likelihood <-dbinom( infmonths , size=(infmonths+notinfmonths), prob=p_grid )
          # compute product of likelihood and prior
          unstd.posterior <- likelihood * uniprior
          # standardize the posterior, so it sums to 1
          posterior <- unstd.posterior / sum(unstd.posterior)
          #sample from the posterior and see if it is > the set threshold for flipping
          #sample(p_grid,1,prob=posterior)
          post<-sample(p_grid,1000,prob=posterior,replace=TRUE)
          sum(post>=flipthreshold)/1000
      })
          
      tolerate<-rbinom(length(pop$infected[pop$infected>0]),1,flipprob)
      pop$strategy[pop$infected>0][tolerate==1]<-2
      
      #count time since last reprod, know if eligible to reproduce again
      pop$reprod[pop$reprod>=1]<-pop$reprod[pop$reprod>=1]+1
      pop$reprod[pop$reprod>9]<-0
      
      #fight infection
      
      #this allows me to add treatment into the model, maybe cut for final version
      # if(i %in% treat) {
      #   treatm<-rbinom(length(pop$infected[pop$infected>0]),1,treatefficacy)
      # } else treatm<-rep(0,length(pop$infected[pop$infected>0]))
      # clearR<-rbinom(length(pop$infected[pop$infected>0]),1,Pcr)*(pop$strategy[pop$infected>0]==1)
      # clearT<-as.numeric(rbinom(length(pop$infected[pop$infected>0]),1,Pct)==1 | treatm==1)
      clearR<-rbinom(length(pop$infected[pop$infected>0]),1,Pcr)*(pop$strategy[pop$infected>0]==1)
      clearT<-rbinom(length(pop$infected[pop$infected>0]),1,Pct)*(pop$strategy[pop$infected>0]==2)
      #spontanous clearance doesn't produce immunity and chance of immunity declines with each subsequent infection (assumes if you were going to develop immunity it would happen with the first or second infection...)
      immune<-rbinom(length(clearR),1,Pim/(pop$previousinfections+1))*clearR
      pop$immune[pop$infected>0]<-immune
      
      pop$previousinfections[pop$infected>0][(clearR==1 | clearT==1)]<- pop$previousinfections[pop$infected>0][(clearR==1 | clearT==1)] + 1
      pop$PIsum[pop$infected>0][(clearR==1 | clearT==1)]<- pop$PIsum[pop$infected>0][(clearR==1 | clearT==1)] + pop$infected[pop$infected>0][(clearR==1 | clearT==1)]
      pop$PIweighted[pop$infected>0][(clearR==1 | clearT==1)]<- pop$PIweighted[pop$infected>0][(clearR==1 | clearT==1)] + pop$infected[pop$infected>0][(clearR==1 | clearT==1)]
      
      if(trackinfections) infections<-rbind(infections,data.frame(month=rep(i,sum(clearR==1 | clearT==1)),length=pop$infected[pop$infected>0][(clearR==1 | clearT==1)],R0=pop$R0[pop$infected>0][(clearR==1 | clearT==1)],RR=pop$RR[pop$infected>0][(clearR==1 | clearT==1)],immunekilled=clearR[(clearR==1 | clearT==1)],strategy=pop$strategy[pop$infected>0][(clearR==1 | clearT==1)],hostdeath=rep(0,sum(clearR==1 | clearT==1)),dieinfect=rep(0,sum(clearR==1 | clearT==1))))
      pop$RR[pop$infected>0][(clearR==1 | clearT==1)]<- 0
      pop$R0[pop$infected>0][(clearR==1 | clearT==1)]<- 0
      pop$infected[pop$infected>0][(clearR==1 | clearT==1)]<- 0
      
      #if still infected, mark infection duration for longer
      pop$infected[pop$infected>0] <- pop$infected[pop$infected>0]+1
      #if not infected add to counter of uninfected time
      pop$PNIweighted[pop$infected==0]<-pop$PNIweighted[pop$infected==0]+1

      #see if infection kills you
      pop$dieinfect<-rbinom(nrow(pop),1,Pmr)*(pop$strategy==1)*(pop$infected>0)+rbinom(nrow(pop),1,Pmt)*(pop$strategy==2)*(pop$infected>0)
      dierisk <- abs((1-((80-pop$age+20)/80)^dieexp))
      dierisk[pop$age>100]<-1
      pop$dieage<-rbinom(length(pop$age),1,dierisk)
      #pop$dieage<-as.numeric(pop$age>70)
      pop$die<-as.numeric(pop$dieinfect==1 | pop$dieage==1)
      addinfect<-pop[pop$die==1 & pop$infected>0,]
      #record infections of dead as ending
      if(trackinfections) infections<-rbind(infections,data.frame(month=rep(i,nrow(addinfect)),length=addinfect$infected,R0=addinfect$R0,RR=addinfect$RR,immunekilled=rep(0,nrow(addinfect)),strategy=addinfect$strategy,hostdeath=rep(1,nrow(addinfect)),dieinfect=addinfect$dieinfect),make.row.names=FALSE)
      
      #infection spreads
      if(sum(pop$infected>0)>0){
        infects <- rpois(sum(pop$infected>0),R0)
        pop$R0[pop$infected>0]<-pop$R0[pop$infected>0]+infects
        RR1<-rep(1:length(infects),infects)
        #select those to infect and infect if not already infected or immune
        samp<-unlist(sapply(infects,function(x) sample(nrow(pop),x,replace=TRUE)))
        #mark immune as exposed, so don't lose immunity
        pop$previousinfections[samp][pop$immune[samp]==1]<-pop$previousinfections[samp][pop$immune[samp]==1]+1
        pop$PIsum[samp][pop$immune[samp]==1]<-pop$PIsum[samp][pop$immune[samp]==1]+2
        pop$PIweighted[samp][pop$immune[samp]==1]<-pop$PIweighted[samp][pop$immune[samp]==1]+2
        
        infectable<-apply(pop[samp,c("infected","immune")],1,function(x) x[1]==0 & x[2]==0)
        RR2<-table(RR1[infectable])
        pop$RR[pop$infected>0][as.numeric(names(RR2))]<-pop$RR[pop$infected>0][as.numeric(names(RR2))]+RR2
        pop$infected[samp[infectable]]<-1
      #add mutation to R0 later? if so should probably be cost related to R0
      }
      
      
      #who will reproduce?
      pop$reprod[pop$reprod==0 & pop$age>=18 & pop$age<=45] <- rbinom(length(pop$reprod[pop$reprod==0 & pop$age>=18 & pop$age<=45]),1,reprodrate)
      pop$RS[pop$reprod==1]<-pop$RS[pop$reprod==1]+1
      #baby has chance of getting mom's strategy if mom is infected
      if(sum(pop$reprod==1)>0){
        babies<-data.frame(ID=paste0(pop$ID[pop$reprod==1],".",pop$RS[pop$reprod==1]),age=rep(0,sum(pop$reprod==1)),strategy=startingstrategy,infected=0,immune=0,previousinfections=0,PIsum=0,PIweighted=0,PNIweighted=priorMon,dieinfect=0,dieage=0,die=0,reprod=0,RS=0,R0=0,RR=0)
        
        #baby gets immunity if mom is infected with Pim and history with weight h2inf
        babies$immune<-rbinom(length(babies$immune),1,Pim*h2inf)*(pop$infected[pop$reprod==1]>0)
        #babies get weighted current infection duration if mom is infected and also her past exposure if h2past>0, regardless of current infection
        babies$PIweighted <- pop$infected[pop$reprod==1]*h2inf+(pop$PIweighted[pop$reprod==1])*h2past
        #babies can start with mom's strategy
        if(inheritstrategy) babies$strategy<-pop$strategy[pop$reprod==1]
        #babies can start with tolerate if right balance of recent mom infection and mom history. Baby has age 9 months.
        # flipprob<-1/(1 + exp(-(((pop$infected[pop$reprod==1]+babies$PIweighted)*flip-infantsample)/infantsample)*log(0.99/0.01)))
        # tolerate<-rbinom(nrow(babies),1,flipprob)
        # babies$strategy[tolerate==1]<-2
        #babies get ID numbers
        
        pop<-rbind(pop,babies,make.row.names=FALSE)
      }
      

      
      pops[[i]]<-pop
      pop<-pop[pop$die==0,]
      #if exceeds K, sample just K individuals to keep
      if(nrow(pop)>K) pop<-pop[sample(nrow(pop),K),]
      if(nrow(pop)==0) break()
      #if(sum(pop$infected>0)==0 & i<5) break()
      if(i %in% round(seq(120,months,120))) message(paste0(model," Month ",i,"/", months,"<br>"))
    }
    if(trackinfections) infections<-infections[!is.na(infections$length),] else infections<-NULL
    return(list(pops=pops,infections=infections,params=params))
  })
}


strategyreport<-function(pops,n=length(pops),...){
  pops<-pops[1:n]
  par(mfrow=c(3,1),mar=c(3,4,1,1))
  infected<-Reduce(rbind,lapply(pops,function(x) c(sum(x$infected==0 & x$immune==0),sum(x$immune==1),sum(x$infected>0))))
  popsize<-infected[,1]+infected[,2]+infected[,3]
  plot(1,1,ylim=c(0,max(popsize)),xlim=c(1,n),type="n",ylab="Individuals",xaxt="n",xaxs="i",...)
  lines(1:nrow(infected),popsize)
  lines(1:nrow(infected),infected[,1],lty=1,col="blue")
  lines(1:nrow(infected),infected[,2],lty=3,col="green")
  lines(1:nrow(infected),infected[,3],lty=3,col="red")
  ax<-axTicks(1,axp=c(0,n/12,10))
  axis(1,at=ax*12,labels=ax)
  
  plot(1,1,ylim=c(0,1),xlim=c(1,n),type="n",ylab="Proportion",xaxt="n",xaxs="i",...)
  lines(1:nrow(infected),infected[,1]/popsize,lty=1,col="blue")
  lines(1:nrow(infected),infected[,2]/popsize,lty=3,col="green")
  lines(1:nrow(infected),infected[,3]/popsize,lty=3,col="red")
  axis(1,at=ax*12,labels=ax)
  
  strategy<-Reduce(rbind,lapply(pops,function(x) c(sum(x$strategy==1),sum(x$strategy==2))))
  plot(1,1,ylim=c(0,1),xlim=c(1,n),type="n",ylab="Proportion",xaxt="n",xaxs="i",...)
  lines(1:nrow(strategy),strategy[,1]/popsize,lty=1,col="orange")
  lines(1:nrow(strategy),strategy[,2]/popsize,lty=1,col="purple")
  axis(1,at=ax*12,labels=ax)
}


strategyreport2<-function(pops,n=length(pops),cols=c("black","blue","green","red","orange","purple"),...){
  pops<-pops[1:n]
  
  layout(matrix(c(1,2,3,4,0,0),ncol=2,byrow=TRUE),widths=c(4,0.75),heights=c(1,1,0.1))
  par(mar=c(3,5,0.2,0))
  
  infected<-Reduce(rbind,lapply(pops,function(x) c(sum(x$infected==0 & x$immune==0),sum(x$immune==1),sum(x$infected>0))))
  popsize<-infected[,1]+infected[,2]+infected[,3]
  plot(1,1,ylim=c(0,max(popsize)),xlim=c(0,n),type="n",ylab="Individuals",xaxt="n",...)
  lines(1:nrow(infected),popsize,lty=2,lwd=2,col=cols[1])
  lines(1:nrow(infected),infected[,1],lty=1,lwd=2,col=cols[2])
  lines(1:nrow(infected),infected[,2],lty=1,lwd=2,col=cols[3])
  lines(1:nrow(infected),infected[,3],lty=1,lwd=2,col=cols[4])
  ax<-axTicks(1,axp=c(0,n/12,10))
  axis(1,at=ax*12,labels=NA)
  text(-(n/7.5),max(popsize), "A",cex=3,xpd=TRUE)
  
  par(mar=c(3,0,0.2,1))
  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  legend(0,1,c("All","Naive","Immune","Infected"),col=cols[1:4],lty=c(2,1,1,1),lwd=2,bty="n",title=expression(bold("Status")))
  
  par(mar=c(3,5,0.2,0))
  strategy<-Reduce(rbind,lapply(pops,function(x) c(sum(x$strategy==1),sum(x$strategy==2))))
  plot(1,1,ylim=c(0,1),xlim=c(1,n),type="n",ylab="Proportion",xaxt="n",...)
  lines(1:nrow(strategy),strategy[,1]/popsize,lty=1,lwd=2,col=cols[5])
  lines(1:nrow(strategy),strategy[,2]/popsize,lty=1,lwd=2,col=cols[6])
  axis(1,at=ax*12,labels=ax)
  mtext("Years",side=1,line=3)
  text(-(n/7.5),1, "B",cex=3,xpd=TRUE)
  
  par(mar=c(3,0,0.2,1))
  plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  legend(0,1,c("Fight","Tolerate"),col=cols[5:6],lty=1,lwd=2,bty="n",title=expression(bold("Strategy")))
}


paramscheck<-function(params,xmax1,xmax2){
  with(params,{
    #graph survival with fight and tolerate
    months<-seq(0,70*12,0.1)
    survivalfight<-sapply(months,function(x) (1-Pmr)^x)
    mediansurvivalfight<-months[which.min(abs(0.50-survivalfight))]
    survivaltolerate<-sapply(months,function(x) (1-Pmt)^x)
    mediansurvivaltolerate<-months[which.min(abs(0.50-survivaltolerate))]
    
    #graph clearance with and without fight
    clearRight<-1-sapply(months,function(x) (1-Pcr)^x)
    medianclearRight<-months[which.min(abs(0.50-clearRight))]
    cleartolerate<-1-sapply(months,function(x) (1-Pct)^x)
    mediancleartolerate<-months[which.min(abs(0.50-cleartolerate))]
    
    #R0 with and without, should just be based on clearance times R0 param
    
    
    layout(matrix(c(1,2,3,4,0,0),ncol=2,byrow=TRUE),widths=c(4,0.75),heights=c(1,1,0.1))
    par(mar=c(3,5,0.2,0))
    plot(1,1,ylim=c(0,1),xlim=c(0,xmax1),type="n",ylab="Proportion",xlab=NA)
    lines(months/12,survivalfight,col=cols[1])
    lines(rep(mediansurvivalfight,2)/12,c(0,0.5),col=cols[1],lty=2)
    
    lines(months/12,clearRight,col=cols[4])
    lines(rep(medianclearRight,2)/12,c(0,0.5),col=cols[4],lty=2)
    text(-(xmax1/7.5),1, "A",cex=3,xpd=TRUE)
    
    par(mar=c(3,0,0.2,1))
    plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    legend(0,1,legend=c("Survival","Clearance"),col=cols[c(1,4)],lty=1,bty="n")
    
    par(mar=c(3,5,0.2,0))
    plot(1,1,ylim=c(0,1),xlim=c(0,xmax2),type="n",ylab="Proportion",xlab="Years")
    lines(months/12,survivaltolerate,col=cols[3])
    lines(rep(mediansurvivaltolerate,2)/12,c(0,0.5),col=cols[3],lty=2)
    lines(months/12,cleartolerate,col=cols[2])
    lines(rep(mediancleartolerate,2)/12,c(0,0.5),col=cols[2],lty=2)
    text(-(xmax2/7.5),1, "B",cex=3,xpd=TRUE)
    
    par(mar=c(3,0,0.2,1))
    plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    legend(0,1,legend=c("Survival","Clearance"),col=cols[c(3,2)],lty=1,bty="n")
  })
}


strategyreport3<-function(pops1,n=length(pops),cols=c("black","blue","green","red","orange","purple"),statlim=c(0.2,1)){
  
  pops<-pops1$pops[1:n]
  infections<-pops1$infections
  params<-pops1$params
  if(nrow(infections)>0){
    infections$monthinf<-infections$month-infections$length
    infstats<-aggregate(cbind(R0,length,dieinfect)~strategy,data=infections[infections$monthinf>=(n*statlim[1]) & infections$monthinf<=(n*statlim[2]),],mean)
    infstats2<-aggregate(cbind(R0,length)~strategy+dieinfect,data=infections[infections$monthinf>=(n*statlim[1]) & infections$monthinf<=(n*statlim[2]),],mean)
    lifeinfect<-Reduce(rbind,lapply(pops[(n*statlim[1]):(n*statlim[2])],function(x) cbind(x$previousinfections[x$die==1]+(x$infected>0)[x$die==1],x$strategy[x$die==1])))
    lifeinfect<-aggregate(lifeinfect[,1]~lifeinfect[,2],FUN=mean)
  }
  
  layout(matrix(c(7,10, 1,2, 8,11, 3,4, 9,12, 5,6, 0,0),ncol=2,byrow=TRUE),widths=c(3,4),heights=c(0.2,1,0.2,1,0.2,1,0.3))
  bottommar<-1
  par(mar=c(bottommar,1,0.2,1))

  
  frame()
  if(nrow(infections)>0){
    xpad<-0.1
    ypad<-0.3
    adj<-c(1,1)
    adj2<-c(0.5,1)
    ybase<-1
    checkwidth<-strwidth(expression(bold("Probability Immunity (%) Tolerate Tolerate")),cex=1)
    cex<-1/checkwidth
    width<-strwidth(expression(bold("Tolerate")),cex=cex)*(1+xpad)
    height<-strheight(expression(bold('R'[0])),cex=cex)*(1+ypad)
    col1<-strwidth(expression(bold("Probability Immunity (%)")),cex=cex)
    adjwidth<-strwidth(expression(bold("8.8")),cex=cex)
    text(col1+width-adjwidth,ybase,expression(bold("Fight")),adj=adj2,cex=cex)
    text(col1+width*2-adjwidth,ybase,expression(bold("Tolerate")),adj=adj2,cex=cex)
    text(col1,ybase-height,expression(bold('R'[0])),cex=cex,adj=adj)
    text(col1+width,ybase-height,format(round(infstats[1,2],1),nsmall=1),adj=adj,cex=cex)
    text(col1+width*2,ybase-height,format(round(infstats[2,2],1),nsmall=1),adj=adj,cex=cex)
    text(col1,ybase-height*2,expression(bold("Duration (Months)")),adj=adj,cex=cex)
    text(col1+width,ybase-height*2,format(round(infstats[1,3],1),nsmall=1),adj=adj,cex=cex)
    text(col1+width*2,ybase-height*2,format(round(infstats[2,3],1),nsmall=1),adj=adj,cex=cex)
    text(col1,ybase-height*3,expression(bold("Case Fatality Rate (%)")),adj=adj,cex=cex)
    text(col1+width,ybase-height*3,format(round(infstats[1,4]*100,1),nsmall=1),adj=adj,cex=cex)
    text(col1+width*2,ybase-height*3,format(round(infstats[2,4]*100,1),nsmall=1),adj=adj,cex=cex)
    text(col1,ybase-height*4,expression(bold("Lifetime Infections")),adj=adj,cex=cex)
    text(col1+width,ybase-height*4,format(round(lifeinfect[1,2],1),nsmall=1),adj=adj,cex=cex)
    text(col1+width*2,ybase-height*4,format(round(lifeinfect[2,2],1),nsmall=1),adj=adj,cex=cex)
    text(col1,ybase-height*5,expression(bold("Probability Immunity (%)")),adj=adj,cex=cex)
    text(col1+width,ybase-height*5,format(round(params$Pim*100,1),nsmall=1),adj=adj,cex=cex)
    text(col1+width*2,ybase-height*5,"--",adj=adj,cex=cex)
  }

  par(mar=c(bottommar,5,0.2,1))
  infected<-Reduce(rbind,lapply(pops,function(x) c(sum(x$infected==0 & x$immune==0),sum(x$immune==1),sum(x$infected>0))))
  popsize<-infected[,1]+infected[,2]+infected[,3]
  plot(1,1,ylim=c(0,max(popsize)),xlim=c(0,n),type="n",ylab="Individuals",xaxt="n")
  lines(1:nrow(infected),popsize,lty=2,lwd=2,col=cols[1])
  lines(1:nrow(infected),infected[,1],lty=1,lwd=2,col=cols[2])
  lines(1:nrow(infected),infected[,2],lty=1,lwd=2,col=cols[3])
  lines(1:nrow(infected),infected[,3],lty=1,lwd=2,col=cols[4])
  ax<-axTicks(1,axp=c(0,n/12,10))
  axis(1,at=ax*12,labels=NA)
  
  #par(mar=c(bottommar,0,0.2,1))
  #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  legend("right",c("All","Susceptible","Immune","Infected"),col=cols[1:4],lty=c(2,1,1,1),lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.7),inset=c(0.05,0),title=expression(bold("Status")))

  
  par(mar=c(bottommar,5,0.2,0))
  vitalrate<-Reduce(rbind,lapply(pops,function(x) c(sum(x$dieinfect==1),sum(x$dieage==1),sum(x$die==1),sum(x$reprod==1))/nrow(x) ))
  vitalrate<-cbind(1:nrow(vitalrate),vitalrate)
  colnames(vitalrate)<-c("month","dieinfect","dieage","die","reprod")
  row.names(vitalrate)<-vitalrate[,1]
  vitalrate<-data.frame(vitalrate)
  vitalrate$year<-round((vitalrate$month+6)/12)
  vitalratey<-aggregate(cbind(dieinfect,dieage,die,reprod)~year,data=vitalrate,sum)
  vitalratey[,2:5]<-vitalratey[,2:5]*1000
  plot(vitalratey$die,type="l",ylab="Mortality Rate",xlab=NA,xaxt="n")
  axis(1,at=ax,labels=NA)

  par(mar=c(bottommar,5,0.2,1))
  plot(1,1,ylim=c(0,1),xlim=c(0,n),type="n",ylab="Individuals",xaxt="n")
  lines(1:nrow(infected),infected[,1]/popsize,lty=1,lwd=2,col=cols[2])
  lines(1:nrow(infected),infected[,2]/popsize,lty=1,lwd=2,col=cols[3])
  lines(1:nrow(infected),infected[,3]/popsize,lty=1,lwd=2,col=cols[4])
  axis(1,at=ax*12,labels=NA)

  par(mar=c(bottommar,5,0.2,0))
  plot(vitalratey$reprod,type="l",ylab="Birth Rate",xlab=NA)
  mtext("Years",side=1,line=3)
  
  par(mar=c(bottommar,5,0.2,1))
  strategy<-Reduce(rbind,lapply(pops,function(x) c(sum(x$strategy==1 & x$immune==0),sum(x$strategy==1 & x$immune==1),sum(x$strategy==2))))
  plot(1,1,ylim=c(0,1),xlim=c(1,n),type="n",ylab="Proportion",xaxt="n")
  lines(1:nrow(strategy),strategy[,1]/popsize,lty=1,lwd=2,col=cols[1])
  lines(1:nrow(strategy),strategy[,2]/popsize,lty=1,lwd=2,col=cols[3])
  lines(1:nrow(strategy),strategy[,3]/popsize,lty=1,lwd=2,col=cols[5])
  axis(1,at=ax*12,labels=ax)
  mtext("Years",side=1,line=3)

  #par(mar=c(bottommar,0,0.2,1))
  #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  legend("right",c("Fight","Immune","Tolerate"),col=cols[c(1,3,5)],lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.7),inset=c(0.05,0),title=expression(bold("Strategy")))
  
  par(mar=c(0,0,0,0))
  for(i in 1:6){
    frame()
    text(0,0,LETTERS[i],cex=2,adj=c(0,0),font=2)
  }
  
}


strategyreportflu<-function(pops1,n=length(pops),cols=c("black","blue","green","red","orange","purple"),statlim=c(0.2,1)){
  
  pops<-pops1$pops[1:n]
  infections<-pops1$infections
  params<-pops1$params

  layout(matrix(c(3,4, 1,2, 0,0),ncol=2,byrow=TRUE),widths=c(4,4),heights=c(0.2,1,0.3))
  bottommar<-1
  
  par(mar=c(bottommar,5,0.2,1))
  infected<-Reduce(rbind,lapply(pops,function(x) c(sum(x$infected==0 & x$immune==0),sum(x$immune==1),sum(x$infected>0))))
  popsize<-infected[,1]+infected[,2]+infected[,3]
  plot(1,1,ylim=c(0,max(popsize)),xlim=c(1,n),type="n",ylab="Individuals",xaxt="n")
  lines(1:nrow(infected),popsize,lty=2,lwd=2,col=cols[1])
  lines(1:nrow(infected),infected[,1],lty=1,lwd=2,col=cols[2])
  lines(1:nrow(infected),infected[,2],lty=1,lwd=2,col=cols[3])
  lines(1:nrow(infected),infected[,3],lty=1,lwd=2,col=cols[4])
  ax<-axTicks(1,axp=c(0,n/12,10))
  axis(1,at=ax*12,labels=ax)
  mtext("Years",side=1,line=3)
  
  #par(mar=c(bottommar,0,0.2,1))
  #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  legend("right",c("All","Susceptible","Immune","Infected"),col=cols[1:4],lty=c(2,1,1,1),lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.7),inset=c(0.05,0),title=expression(bold("Status")))
  
  par(mar=c(bottommar,5,0.2,1))
  strategy<-Reduce(rbind,lapply(pops,function(x) c(sum(x$strategy==1 & x$immune==0),sum(x$strategy==1 & x$immune==1),sum(x$strategy==2))))
  plot(1,1,ylim=c(0,1),xlim=c(1,n),type="n",ylab="Proportion",xaxt="n")
  lines(1:nrow(strategy),strategy[,1]/popsize,lty=1,lwd=2,col=cols[1])
  lines(1:nrow(strategy),strategy[,2]/popsize,lty=1,lwd=2,col=cols[3])
  lines(1:nrow(strategy),strategy[,3]/popsize,lty=1,lwd=2,col=cols[5])
  axis(1,at=ax*12,labels=ax)
  mtext("Years",side=1,line=3)
  
  #par(mar=c(bottommar,0,0.2,1))
  #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
  legend("right",c("Fight","Immune","Tolerate"),col=cols[c(1,3,5)],lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.7),inset=c(0.05,0),title=expression(bold("Strategy")))
  
  par(mar=c(0,0,0,0))
  for(i in 1:2){
    frame()
    text(0,0,LETTERS[i],cex=2,adj=c(0,0),font=2)
  }
  
}


strategyreportcompare<-function(popsAll,n=length(pops),cols=c("black","blue","green","red","orange","purple"),statlim=c(0.2,1)){
  
  
  layout(matrix(c(seq(1,length(popsAll)*3,1), 0,0,0),ncol=3,byrow=TRUE),widths=c(4,4,4),heights=c(rep(1,length(popsAll)),0.3))
  bottommar<-1
  ii<-1
  for(pops1 in popsAll){
    pops<-pops1$pops[1:n]
    infections<-pops1$infections
    params<-pops1$params
    
    if(nrow(infections)>0){
      infections$monthinf<-infections$month-infections$length
      infstats<-aggregate(cbind(R0,length,dieinfect)~strategy,data=infections[infections$monthinf>=(n*statlim[1]) & infections$monthinf<=(n*statlim[2]),],mean)
      infstats2<-aggregate(cbind(R0,length)~strategy+dieinfect,data=infections[infections$monthinf>=(n*statlim[1]) & infections$monthinf<=(n*statlim[2]),],mean)
      lifeinfect<-Reduce(rbind,lapply(pops[(n*statlim[1]):(n*statlim[2])],function(x) cbind(x$previousinfections[x$die==1]+(x$infected>0)[x$die==1],x$strategy[x$die==1])))
      lifeinfect<-aggregate(lifeinfect[,1]~lifeinfect[,2],FUN=mean)
    }
    
    par(mar=c(bottommar,1,0.2,0))
    
    frame()
    text(0,1,LETTERS[ii],cex=2,adj=c(0,1),font=2)
    if(nrow(infections)>0){
      xpad<-0.1
      ypad<-0.3
      adj<-c(1,1)
      adj2<-c(0.5,1)
      ybase<-0.9
      checkwidth<-strwidth(expression(bold("Probability Immunity (%) Tolerate Tolerate")),cex=1)
      cex<-1/checkwidth
      width<-strwidth(expression(bold("Tolerate")),cex=cex)*(1+xpad)
      height<-strheight(expression(bold('R'[0])),cex=cex)*(1+ypad)
      col1<-strwidth(expression(bold("Probability Immunity (%)")),cex=cex)
      adjwidth<-strwidth(expression(bold("8.8")),cex=cex)
      text(col1+width-adjwidth,ybase,expression(bold("Fight")),adj=adj2,cex=cex)
      text(col1+width*2-adjwidth,ybase,expression(bold("Tolerate")),adj=adj2,cex=cex)
      text(col1,ybase-height,expression(bold('R'[0])),cex=cex,adj=adj)
      text(col1+width,ybase-height,format(round(infstats[1,2],1),nsmall=1),adj=adj,cex=cex)
      text(col1+width*2,ybase-height,format(round(infstats[2,2],1),nsmall=1),adj=adj,cex=cex)
      text(col1,ybase-height*2,expression(bold("Mean Duration (Months)")),adj=adj,cex=cex)
      text(col1+width,ybase-height*2,format(round(infstats[1,3],1),nsmall=1),adj=adj,cex=cex)
      text(col1+width*2,ybase-height*2,format(round(infstats[2,3],1),nsmall=1),adj=adj,cex=cex)
      text(col1,ybase-height*3,expression(bold("Case Fatality Rate (%)")),adj=adj,cex=cex)
      text(col1+width,ybase-height*3,format(round(infstats[1,4]*100,1),nsmall=1),adj=adj,cex=cex)
      text(col1+width*2,ybase-height*3,format(round(infstats[2,4]*100,1),nsmall=1),adj=adj,cex=cex)
      text(col1,ybase-height*4,expression(bold("Lifetime Infections")),adj=adj,cex=cex)
      text(col1+width,ybase-height*4,format(round(lifeinfect[1,2],1),nsmall=1),adj=adj,cex=cex)
      text(col1+width*2,ybase-height*4,format(round(lifeinfect[2,2],1),nsmall=1),adj=adj,cex=cex)
      text(col1,ybase-height*5,expression(bold("Probability Immunity (%)")),adj=adj,cex=cex)
      text(col1+width,ybase-height*5,format(round(params$Pim*100,1),nsmall=1),adj=adj,cex=cex)
      text(col1+width*2,ybase-height*5,"--",adj=adj,cex=cex)
      text(col1,ybase-height*6,expression(bold("Probability Switch (%)")),adj=adj,cex=cex)
      #flipprob<-1/(1 + exp(-((9*params$flip-9)/9)*log(0.999/0.001)))
      #text(col1+width,ybase-height*6,format(round(flipprob*100,1),nsmall=1),adj=adj,cex=cex)
      text(col1+width*2,ybase-height*6,"--",adj=adj,cex=cex)
      text(col1,ybase-height*7,expression(bold("Inherit Strategy")),adj=adj,cex=cex)
      text(col1+width,ybase-height*7,ifelse(params$inheritstrategy,"Yes","No"),adj=adj,cex=cex)
      text(col1+width*2,ybase-height*7,ifelse(params$inheritstrategy,"Yes","No"),adj=adj,cex=cex)
    }
  
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
    if(ii==length(popsAll)) {
      mtext("Years",side=1,line=3)
      axis(1,at=ax*12,labels=ax)
    } else axis(1,at=ax*12,labels=NA)
    
    #par(mar=c(bottommar,0,0.2,1))
    #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    legend("topleft",c("All","Susceptible","Immune","Infected"),col=cols[1:4],lty=c(2,1,1,1),lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0.05),title=expression(bold("Status")))
    
    par(mar=c(bottommar,4,1,1))
    strategy<-Reduce(rbind,lapply(pops,function(x) c(sum(x$strategy==1 & x$immune==0),sum(x$immune==1),sum(x$strategy==2 & x$immune==0))))
    plot(1,1,ylim=c(0,1),xlim=c(1,n),type="n",ylab="Proportion",xaxt="n")
    lines(1:nrow(strategy),strategy[,1]/popsize,lty=1,lwd=2,col=cols[1])
    lines(1:nrow(strategy),strategy[,2]/popsize,lty=1,lwd=2,col=cols[3])
    lines(1:nrow(strategy),strategy[,3]/popsize,lty=1,lwd=2,col=cols[5])
    if(ii==length(popsAll)) {
      mtext("Years",side=1,line=3)
      axis(1,at=ax*12,labels=ax)
    } else axis(1,at=ax*12,labels=NA)
    
    #par(mar=c(bottommar,0,0.2,1))
    #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    legend("right",c("Fight","Immune","Tolerate"),col=cols[c(1,3,5)],lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0),title=expression(bold("Strategy")))
    ii<-ii+1
  }
}


strategyreportcompare2<-function(popsAll,n=length(pops),cols=c("black","blue","green","red","orange","purple"),statlim=c(0.2,1)){
  
  layout(matrix(c(seq(1,length(popsAll)*3,1), 0,0,0),ncol=3,byrow=TRUE),widths=c(0.5,4,4),heights=c(rep(1,length(popsAll)),0.3))
  bottommar<-1
  ii<-1
  for(pops1 in popsAll){
    pops<-pops1$pops[1:n]
    infections<-pops1$infections
    params<-pops1$params
    
    par(mar=c(bottommar,1,0.2,0))
    
    frame()
    text(0,1,LETTERS[ii],cex=2,adj=c(0,1),font=2)
    
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
    if(ii==length(popsAll)) {
      mtext("Years",side=1,line=3)
      axis(1,at=ax*12,labels=ax)
    } else axis(1,at=ax*12,labels=NA)
    
    #par(mar=c(bottommar,0,0.2,1))
    #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    legend("topleft",c("All","Susceptible","Immune","Infected"),col=cols[1:4],lty=c(2,1,1,1),lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0.05),title=expression(bold("Status")))
    
    par(mar=c(bottommar,4,1,1))
    strategy<-Reduce(rbind,lapply(pops,function(x) c(sum(x$strategy==1 & x$immune==0),sum(x$immune==1),sum(x$strategy==2 & x$immune==0))))
    plot(1,1,ylim=c(0,1),xlim=c(1,n),type="n",ylab="Proportion",xaxt="n")
    lines(1:nrow(strategy),strategy[,1]/popsize,lty=1,lwd=2,col=cols[1])
    lines(1:nrow(strategy),strategy[,2]/popsize,lty=1,lwd=2,col=cols[3])
    lines(1:nrow(strategy),strategy[,3]/popsize,lty=1,lwd=2,col=cols[5])
    if(ii==length(popsAll)) {
      mtext("Years",side=1,line=3)
      axis(1,at=ax*12,labels=ax)
    } else axis(1,at=ax*12,labels=NA)
    
    #par(mar=c(bottommar,0,0.2,1))
    #plot(1,1,xlim=c(0,1),ylim=c(0,1),bty="n",type="n",xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    legend("right",c("Fight","Immune","Tolerate"),col=cols[c(1,3,5)],lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0),title=expression(bold("Strategy")))
    ii<-ii+1
  }
}