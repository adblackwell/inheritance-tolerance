
x<-allpops[[2]]$pops[[1]]
p_grid <- seq( from=0 , to=1 , length.out=1000 )
uniprior <- rep( 1 , 100 )
infectionweight<-2
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
  sum(post>=0.75)/1000
})

plot(flipprob~x$age,xlim=c(0,50))
