x<-allpops[[3]]$pops[[400]]
infected<-x$infected
PIweighted<-x$PIweighted
PNIweighted<-x$PNIweighted
flipthreshold<-0.75
infectionweight<-1
PIdecline<-0.98
priorMon<-0

plot(aggregate(strategy~age,data=x,mean))


flipprob <- apply(x[,c("infected","PIweighted","PNIweighted")],1,function(pp){
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

PNIweightedbase<-sapply(x$age,function(x) sum(PIdecline^(seq(1,x*12+priorMon))))
