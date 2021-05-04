

## --we want the true probability resistance will work given past experience ##

#months without/total months



pri<-c(0,8)
infectionweight<-2
flipprob<-0.5
## --now for example, show how updates with infection or not
## -- this is the probability that an infection will be cleared by resistance
# define grid
p_grid <- seq( from=0 , to=1 , length.out=100 )
# define prior as equivalent to pri[1]/pri[2] months uninfected
prior <- rep( 1 , 100 )
likelihood <- dbinom( pri[1] , size=pri[2] , prob=p_grid )
# compute product of likelihood and prior
unstd.posterior <- likelihood * prior
# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

par(mfrow=c(3,3))
plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )

post<-sample(p_grid,1000,prob=posterior,replace=TRUE)
flips<-sum(post<0.4)/length(post)

for(i in 1:18){
  prior2<-posterior
  likelihood <- dbinom( infweight , size=infweight , prob=p_grid )
  # compute product of likelihood and prior
  unstd.posterior <- likelihood * prior2
  # standardize the posterior, so it sums to 1
  posterior <- unstd.posterior / sum(unstd.posterior)
  plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
  abline(v=0.5)
  post<-sample(p_grid,1000,prob=posterior,replace=TRUE)
  flips2<-sum(post<0.4)/length(post)
  flips<-c(flips,flips2)
}

flip<-sample(p_grid,1,prob=posterior)<flipprob

#compare with old method
times<-1:18
ages<-rep(0,18)
flipw<-1.5

flipprob<-1/(1 + exp(-(((times)*flipw-(ages*12+9))/(ages*12+9))*log(0.999/0.001)))

plot(times,flipprob)
abline(v=9)

lines(times,flips[1:18])


#check method used in model


likelihood <-dbinom( 48 , size=127, prob=p_grid )
unstd.posterior <- likelihood * prior
# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)
plot( p_grid , posterior , type="b" , xlab="probability of water" , ylab="posterior probability" )

p_grid <- seq( from=0 , to=1 , length.out=100 )
uniprior <- rep( 1 , 100 )
pop<-data.frame(infected=0:8,PIweighted=0,PNIweighted=12*0.98^(0:8))
infectionweight<-3
flipprob <- apply(pop[,c("infected","PIweighted","PNIweighted")],1,function(pp){
  infmonths<-round(pp[1]+pp[2])*infectionweight
  notinfmonths<-round(pp[3])
  likelihood <-dbinom( infmonths , size=(infmonths+notinfmonths), prob=p_grid )
  # compute product of likelihood and prior
  unstd.posterior <- likelihood * uniprior
  # standardize the posterior, so it sums to 1
  posterior <- unstd.posterior / sum(unstd.posterior)
  #sample from the posterior and see if it is > the set threshold for flipping
  post<-sample(p_grid,1000,prob=posterior,replace=TRUE)
  sum(post>=0.5)/1000
})

plot(1:9,flipprob)
abline(h=9)

par(mfrow=c(3,3))
for(i in 1:9) plot( p_grid , flipprob[i] , type="b" , xlab="probability of flipping" , ylab="posterior probability" )
