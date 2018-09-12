times<-1:18
ages<-rep(0,18)
flip<-0

flipprob<-1/(1 + exp(-(((times)*flip-(ages*12+9))/(ages*12+9))*log(0.999/0.001)))

plot(times,flipprob)
abline(v=9)
