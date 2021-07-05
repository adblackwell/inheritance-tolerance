##changelog
##2.2 Added Bayesian flip function
##3.0 Clean up for publication
## ---- modelfunctions ----

#function to quickly repeat default parameters for multiple simulations
makeparams<-function(n0=1000, #starting population size
                     months=1500, #how long to run the model in months
                     h2current=0, #percent of mom's current infection value inherited
                     h2past=0, #percent of mom's past infection history inherited
                     K=3000, #maximum population size
                     reprodrate=0.015, #Monthly risk of reproduction for reproductive age individual, F in the paper
                     gestation=12, #Minimum time between reproductions, representing gestation and lactation
                     PI=0.5, #Poisson distribution lambda (mean) number of individuals an infection infects each month, referred to as P_I in the paper
                     infectedT0=10, #number of individuals to infect at time zero. Starting with just 1 infection often results in a failure for the infection to take hold in the population.
                     startingstrategy=1, #1 = resist, 2=tolerate
                     Pmr=0.025, #Monthly probability of mortality when resisting
                     Pmt=0.0001, #Monthly probability of mortality when tolerating
                     Pcr=0.20, #Monthly probability of clearance when resisting
                     Pct=0.01, #Monthly probability of clearance when tolerating
                     Pim=0.2, #Probability of gaining immunity after clearance
                     Ewdecline=0.98, #Proportion of immunological memory retained each month. when not infected, 'delta' in the paper
                     NIwdecline=1, #Proportion of uninfected time to retain each month. Effectively determines how far back 'memory' of not being infected goes.
                     inheritstrategy=FALSE, #whether offspring directly inherit mom's strategy (rather than just her history)
                     flipthreshold=0.75, #proportion of posterior past this threshold determines prob of flipping, Tflip in the paper
                     dieexp=0.002, #exponent for baseline non-infection mortality curve
                     priorMon=0, #infants are born with an innate prior equal to this many months not infected
                     infectionweight=1, #how much infection months count relative to uninfected in updating prior
                     noflip=FALSE, #When TRUE, there is no strategy switching
                     trackinfections=TRUE){ #When true, infections are recorded, but makes the final dataset larger
  params<-as.list(environment())
}

#Function to run model with a set of parameters
#Parameter 'params' is a set of parameters from makeparams
#Parameter 'model' is just for html output to indicate which model is running
#Parameter 'ageburn' determines the number of months used to establish the stable population structure before introducing the infection.
runmodel<-function(params,ageburn=2000,model=NA){
  with(params,{    
    #starting population
    pop<-data.frame(age=abs(rnorm(n0,0,40)),strategy=startingstrategy,infected=0,immune=0,previousinfections=0,PIsum=0,Eweighted=0,NIweighted=0,dieinfect=0,dieage=0,die=0,reprod=0,RS=0,PI=0,RR=0)
    #baseline immune history for uninfected population
    pop$NIweighted<-sapply(pop$age,function(x) sum(NIwdecline^(seq(1,x*12+priorMon))))
    #burn in age structure
    for(i in 1:ageburn){
      pop$age<-pop$age+(1/12)
      pop$NIweighted<-pop$NIweighted*NIwdecline+1
      #count time since last reprod, know if eligible to reproduce again
      pop$reprod[pop$reprod>=1]<-pop$reprod[pop$reprod>=1]+1
      pop$reprod[pop$reprod>gestation]<-0
      #reproduce
      pop$reprod[pop$reprod==0 & pop$age>=18 & pop$age<=45] <- rbinom(length(pop$reprod[pop$reprod==0 & pop$age>=18 & pop$age<=45]),1,reprodrate)
      pop$RS[pop$reprod==1]<-pop$RS[pop$reprod==1]+1
      if(sum(pop$reprod==1)>0){
        babies<-data.frame(age=rep(0,sum(pop$reprod==1)),strategy=startingstrategy,infected=0,immune=0,previousinfections=0,PIsum=0,Eweighted=0,NIweighted=priorMon,dieinfect=0,dieage=0,die=0,reprod=0,RS=0,PI=0,RR=0)
        #baby gets mom's strategy perfectly for this start, just in case I mix them later and don't want proportion to change with reproduction
        babies$strategy<-pop$strategy[pop$reprod==1]
        pop<-rbind(pop,babies)
      }
      dierisk <- (1-((100-pop$age)/100)^dieexp)
      dierisk[pop$age>100]<-1
      pop$dieage<-rbinom(length(pop$age),1,dierisk)
      pop$die<-as.numeric(pop$dieage==1)
      pop<-pop[pop$die==0,]
      #keep sample at starting size
      if(nrow(pop)>n0) pop<-pop[sample(nrow(pop),n0),]
    }
    #assign ID numbers
    pop$ID<-as.character(seq(1,nrow(pop),1))
    
    message("Model", model, " Finished Age Structure Burn <br>")

    #infect starting number
    pop$infected[sample(nrow(pop),infectedT0)]<-1
    
    pops<-list(pop)
    
    #track infections by adding to this whenever an infection ends
    if(trackinfections) infections<-data.frame(month=1,length=NA,PI=NA,RR=NA,immunekilled=NA,strategy=NA,hostdeath=NA,dieinfect=NA)
    
    for(i in 2:months){
      #age
      pop$age<-pop$age+(1/12)
      #cumulative past history gets down-weighted each round, by same decline as loss of tolerance
      pop$NIweighted<-pop$NIweighted*NIwdecline
      pop$Eweighted<-pop$Eweighted*Ewdecline
      pop$Eweighted[pop$Eweighted<(Ewdecline^2)]<-0
      
      #chance to lose tolerance/immunity (immunological memory) if not infected and Eweighted reaches zero
      pop$immune[pop$infected==0 & pop$Eweighted==0]<-0
      pop$strategy[pop$infected==0 & pop$Eweighted==0]<-1
      
      #check for strategy switching, can only switch if infected (do this after infection so infants can switch right after infected)
      # define grid
      p_grid <- seq( from=0 , to=1 , length.out=1000 )
      uniprior <- rep( 1 , 100 )
      #posterior is based on weighted months infected x a multiplier and weighted months uninfected. Weighting is to simulate gradual decline in immunological memory
      #posterior is the likelihood that an infection is endemic/chronic and cannot or should not be resisted
      
      if(noflip) flipprob<-0 
      else flipprob <- apply(pop[pop$infected>0,c("infected","Eweighted","NIweighted")],1,function(pp){
          infmonths<-round(pp[1]+pp[2])*infectionweight
          notinfmonths<-round(pp[3])
          likelihood <-dbinom( infmonths , size=(infmonths+notinfmonths), prob=p_grid )
          unstd.posterior <- likelihood * uniprior
          posterior <- unstd.posterior / sum(unstd.posterior)
          post<-sample(p_grid,1000,prob=posterior,replace=TRUE)
          sum(post>=flipthreshold)/1000
      })
          
      tolerate<-rbinom(length(pop$infected[pop$infected>0]),1,flipprob)
      pop$strategy[pop$infected>0][tolerate==1]<-2
      
      #count time since last reprod, know if eligible to reproduce again
      pop$reprod[pop$reprod>=1]<-pop$reprod[pop$reprod>=1]+1
      pop$reprod[pop$reprod>gestation]<-0
      
      #fight infection
      clearR<-rbinom(length(pop$infected[pop$infected>0]),1,Pcr)*(pop$strategy[pop$infected>0]==1)
      clearT<-rbinom(length(pop$infected[pop$infected>0]),1,Pct)*(pop$strategy[pop$infected>0]==2)
      #spontaneous clearance doesn't produce immunity and chance of immunity declines with each subsequent infection (assumes if you were going to develop immunity it would happen with the first or second infection...)
      immune<-rbinom(length(clearR),1,Pim/(pop$previousinfections+1))*clearR
      pop$immune[pop$infected>0]<-immune
      
      pop$previousinfections[pop$infected>0][(clearR==1 | clearT==1)]<- pop$previousinfections[pop$infected>0][(clearR==1 | clearT==1)] + 1
      pop$PIsum[pop$infected>0][(clearR==1 | clearT==1)]<- pop$PIsum[pop$infected>0][(clearR==1 | clearT==1)] + pop$infected[pop$infected>0][(clearR==1 | clearT==1)]
      pop$Eweighted[pop$infected>0][(clearR==1 | clearT==1)]<- pop$Eweighted[pop$infected>0][(clearR==1 | clearT==1)] + pop$infected[pop$infected>0][(clearR==1 | clearT==1)]
      
      if(trackinfections) infections<-rbind(infections,data.frame(month=rep(i,sum(clearR==1 | clearT==1)),length=pop$infected[pop$infected>0][(clearR==1 | clearT==1)],PI=pop$PI[pop$infected>0][(clearR==1 | clearT==1)],RR=pop$RR[pop$infected>0][(clearR==1 | clearT==1)],immunekilled=clearR[(clearR==1 | clearT==1)],strategy=pop$strategy[pop$infected>0][(clearR==1 | clearT==1)],hostdeath=rep(0,sum(clearR==1 | clearT==1)),dieinfect=rep(0,sum(clearR==1 | clearT==1))))
      pop$RR[pop$infected>0][(clearR==1 | clearT==1)]<- 0
      pop$PI[pop$infected>0][(clearR==1 | clearT==1)]<- 0
      pop$infected[pop$infected>0][(clearR==1 | clearT==1)]<- 0
      
      #if still infected, mark infection duration for longer
      pop$infected[pop$infected>0] <- pop$infected[pop$infected>0]+1
      #if not infected add to counter of uninfected time
      pop$NIweighted[pop$infected==0]<-pop$NIweighted[pop$infected==0]+1

      #see who infection kills
      pop$dieinfect<-rbinom(nrow(pop),1,Pmr)*(pop$strategy==1)*(pop$infected>0)+rbinom(nrow(pop),1,Pmt)*(pop$strategy==2)*(pop$infected>0)
      #baseline mortality curve for death due to other causes
      dierisk <- abs((1-((100-pop$age)/80)^dieexp))
      #max lifespan 100
      dierisk[pop$age>100]<-1
      pop$dieage<-rbinom(length(pop$age),1,dierisk)
      pop$die<-as.numeric(pop$dieinfect==1 | pop$dieage==1)
      addinfect<-pop[pop$die==1 & pop$infected>0,]
      #record infections of dead as ending
      if(trackinfections) infections<-rbind(infections,data.frame(month=rep(i,nrow(addinfect)),length=addinfect$infected,PI=addinfect$PI,RR=addinfect$RR,immunekilled=rep(0,nrow(addinfect)),strategy=addinfect$strategy,hostdeath=rep(1,nrow(addinfect)),dieinfect=addinfect$dieinfect),make.row.names=FALSE)
      
      #infection spreads
      if(sum(pop$infected>0)>0){
        infects <- rpois(sum(pop$infected>0),PI)
        pop$PI[pop$infected>0]<-pop$PI[pop$infected>0]+infects
        RR1<-rep(1:length(infects),infects)
        #select those to infect and infect if not already infected or immune
        samp<-unlist(sapply(infects,function(x) sample(nrow(pop),x,replace=TRUE)))
        #mark immune as exposed again, so don't lose immunity
        pop$previousinfections[samp][pop$immune[samp]==1]<-pop$previousinfections[samp][pop$immune[samp]==1]+1
        pop$PIsum[samp][pop$immune[samp]==1]<-pop$PIsum[samp][pop$immune[samp]==1]+2
        pop$Eweighted[samp][pop$immune[samp]==1]<-pop$Eweighted[samp][pop$immune[samp]==1]+2
        
        infectable<-apply(pop[samp,c("infected","immune")],1,function(x) x[1]==0 & x[2]==0)
        RR2<-table(RR1[infectable])
        pop$RR[pop$infected>0][as.numeric(names(RR2))]<-pop$RR[pop$infected>0][as.numeric(names(RR2))]+RR2
        pop$infected[samp[infectable]]<-1
      }
      
      
      #who will reproduce?
      pop$reprod[pop$reprod==0 & pop$age>=18 & pop$age<=45] <- rbinom(length(pop$reprod[pop$reprod==0 & pop$age>=18 & pop$age<=45]),1,reprodrate)
      pop$RS[pop$reprod==1]<-pop$RS[pop$reprod==1]+1
      #baby has chance of getting mom's strategy if mom is infected
      if(sum(pop$reprod==1)>0){
        babies<-data.frame(ID=paste0(pop$ID[pop$reprod==1],".",pop$RS[pop$reprod==1]),age=rep(0,sum(pop$reprod==1)),strategy=startingstrategy,infected=0,immune=0,previousinfections=0,PIsum=0,Eweighted=0,NIweighted=priorMon,dieinfect=0,dieage=0,die=0,reprod=0,RS=0,PI=0,RR=0)
        
        #baby gets immunity if mom is infected and resisting with Pim and history with weight h2current
        babies$immune<-rbinom(length(babies$immune),1,Pim*h2current)*(pop$infected[pop$reprod==1]>0 & pop$strategy[pop$reprod==1]==1)
        #babies get weighted current infection duration if mom is infected and also her past exposure if h2past>0, regardless of current infection
        babies$Eweighted <- pop$infected[pop$reprod==1]*h2current+(pop$Eweighted[pop$reprod==1])*h2past
        #babies can start with mom's strategy if inheritstrategy=TRUE
        if(inheritstrategy) babies$strategy<-pop$strategy[pop$reprod==1]
        #babies get ID numbers
        pop<-rbind(pop,babies,make.row.names=FALSE)
      }
      
      pops[[i]]<-pop
      pop<-pop[pop$die==0,]
      #if exceeds K, sample just K individuals to keep. Used to keep simulation from becoming unmanagable
      if(nrow(pop)>K) pop<-pop[sample(nrow(pop),K),]
      if(nrow(pop)==0) break()
      #if(sum(pop$infected>0)==0 & i<5) break()
      if(i %in% round(seq(120,months,120))) message(paste0("Model", model," Month ",i,"/", months,"<br>"))
    }
    if(trackinfections) infections<-infections[!is.na(infections$length),] else infections<-NULL
    return(list(pops=pops,infections=infections,params=params))
  })
}


#Function to graph infection characteristics for survival and clearance.
paramscheck<-function(params,xmax1=70,xmax2=70,nolayout=FALSE,legend="topright",letters=c("A","B")){
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
    
    #PI with and without, should just be based on clearance times PI param
    
    
    if(!nolayout) layout(matrix(c(1,2,0),ncol=1,byrow=TRUE),heights=c(1,1,0.1))
    par(mar=c(3,5,0.5,0.5))
    plot(1,1,ylim=c(0,1),xlim=c(0,xmax1),type="n",ylab="Proportion",xlab=NA)
    lines(months/12,survivalfight,col=cols[1])
    lines(rep(mediansurvivalfight,2)/12,c(0,0.5),col=cols[1],lty=2)
    
    lines(months/12,clearRight,col=cols[4])
    lines(rep(medianclearRight,2)/12,c(0,0.5),col=cols[4],lty=2)
    #text(-(xmax1/7.5),1, letters[1],cex=3,xpd=TRUE)
    mtext(letters[1],side=3,line=-1.5,at=-(xmax1/2.65),cex=1.5)
    
    if(!is.na(legend)){
      legend(legend,legend=c("Survival","Clearance"),col=cols[c(1,4)],lty=1,bty="n")
    }
    
    par(mar=c(3,5,0.5,0.5))
    plot(1,1,ylim=c(0,1),xlim=c(0,xmax2),type="n",ylab="Proportion",xlab="Years")
    lines(months/12,survivaltolerate,col=cols[3])
    lines(rep(mediansurvivaltolerate,2)/12,c(0,0.5),col=cols[3],lty=2)
    lines(months/12,cleartolerate,col=cols[2])
    lines(rep(mediancleartolerate,2)/12,c(0,0.5),col=cols[2],lty=2)
    #text(-(xmax2/7.5),1, letters[2],cex=3,xpd=TRUE)
    mtext(letters[2],side=3,line=-1.5,at=-(xmax2/2.65),cex=1.5)
    
    if(!is.na(legend)){
      legend(legend,legend=c("Survival","Clearance"),col=cols[c(3,2)],lty=1,bty="n")
    }
  })
}

#Function to make a plot comparing models
strategyreport<-function(allpops,n=length(allpops[[1]]$pops),cols=c("black","blue","green","red","orange","purple"),popmax=NA){
  
  layout(rbind(matrix(seq(1,length(allpops)*2,1),ncol=length(allpops),byrow=FALSE),rep(0,length(allpops))),widths=rep(4,length(allpops)),heights=c(1,1,0.3))
  bottommar<-1
  ii<-1
  for(pops1 in allpops){
    pops<-pops1$pops[1:n]
    #infections<-pops1$infections
    #params<-pops1$params
    
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

#function to plot survival curves from simulations
survplot<-function(allpops,n1=1,n2=length(allpops[[1]]$pops),cols=rainbow(length(allpops)),labs=seq(1:length(allpops))){
  require(survival)
  ii<-1
  sss<-list()
  for(pops1 in allpops){
    pops <- pops1$pops[n1:n2]
    n <- (n2-n1)+1
    vitalrate2<-Reduce(rbind,lapply(1:n,function(x) expand.grid(month=x,agedeath=pops[[x]]$age[pops[[x]]$die==1],die=1,mod=1)))
    vitalrate2<-rbind(vitalrate2,expand.grid(month=n,agedeath=pops[[n]]$age[pops[[n]]$die==0],die=0,mod=1))
    vitalrate2$monthborn<-round(vitalrate2$month-vitalrate2$agedeath*12)
    s<-survfit(Surv(agedeath,die)~1,data=vitalrate2)
    sss[[ii]]<-s
    ii<-ii+1
  }
  
  plot(1,1,ylim=c(0,1),xlim=c(1,100),type="n",ylab="Survival (%)",xlab=NA,xaxt="n")
  ii<-1
  for(s in sss){
    lines(s,conf.int=FALSE,col=cols[ii])
    ii<-ii+1
  }
  axis(1,at=seq(1,100,10),labels=NA)
  #mtext("Age",side=1,line=3)
  legend("topright",legend=labs,col=cols,lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0.05),cex=1)
}



#function to plot survival curves from simulations
mortalityplot<-function(allpops,n1=1,n2=length(allpops[[1]]$pops),cols=rainbow(length(allpops)),labs=seq(1:length(allpops))){
  require(survival)
  ii<-1
  sss<-list()
  for(pops1 in allpops){
    pops <- pops1$pops[n1:n2]
    n <- (n2-n1)+1
    vitalrate2<-Reduce(rbind,lapply(1:n,function(x) expand.grid(month=x,agedeath=pops[[x]]$age[pops[[x]]$die==1],die=1,mod=1)))
    vitalrate2<-rbind(vitalrate2,expand.grid(month=n,agedeath=pops[[n]]$age[pops[[n]]$die==0],die=0,mod=1))
    vitalrate2$monthborn<-round(vitalrate2$month-vitalrate2$agedeath*12)
    s<-survfit(Surv(agedeath,die)~1,data=vitalrate2)
    sss[[ii]]<-s
    ii<-ii+1
  }
  
  plot(1,1,ylim=c(0,20),xlim=c(1,100),type="n",ylab="Mortality (%)",xlab=NA)
  ii<-1
  for(s in sss){
    group<-floor(s$time)
    mrt<-aggregate(list(mrt=s$n.event/s$n.risk),by=list(Time=group),sum)
    #drate<-diff((1-s$surv)[seq(1,length(s$cumhaz),12)])
    #dtime<-s$time[seq(1,length(s$surv),12)][-1]
    #lines(drate[-length(drate)]~dtime[-length(dtime)],col=cols[ii])
    lines(I(mrt$mrt*100)~I(mrt$Time+c(diff(mrt$Time)/2,0.5)),col=cols[ii])
    ii<-ii+1
  }
  mtext("Age",side=1,line=3)
  legend("topleft",legend=labs,col=cols,lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0.05),cex=1)
}


#function to plot infection prevalence by age
prevplot<-function(allpops,n1=1,n2=length(allpops[[1]]$pops),cols=rainbow(length(allpops)),labs=seq(1:length(allpops))){
  ii<-1
  sss<-list()
  for(pops1 in allpops){
    pops <- pops1$pops[n1:n2]
    n <- (n2-n1)+1
    vitalrate2<-Reduce(rbind,lapply(1:n,function(x) data.frame(month=x,age=floor(pops[[x]]$age/5)*5,infected=as.numeric(pops[[x]]$infected>1))))
    vitalrate3<-table(vitalrate2[,2:3])
    if(ncol(vitalrate3)==1) vitalrate3<-cbind(vitalrate3,rep(0,nrow(vitalrate3)))
    sss[[ii]]<-vitalrate3[,2]/rowSums(vitalrate3)
    ii<-ii+1
  }
  #ymax<-max(unlist(lapply(sss,max)))*100
  ymax<-100
  plot(1,1,ylim=c(0,ymax),xlim=c(0,95),type="n",ylab="Prevalence (%)",xlab=NA,xaxt="n")
  ii<-1
  for(s in sss){
    s<-s[as.numeric(names(s))<95] #exclude 95 since very few make it to this age and so it jumps around
    lines(as.numeric(names(s))+2.5,s*100,col=cols[ii])
    ii<-ii+1
  }
  axis(1,at=seq(0,95,5),labels=NA)
  #mtext("Age",side=1,line=3)
  legend("topright",legend=labs,col=cols,lty=1,lwd=2,bty="o",box.col=NA,bg=rgb(1,1,1,0.5),inset=c(0.05,0.05),cex=1)
}
