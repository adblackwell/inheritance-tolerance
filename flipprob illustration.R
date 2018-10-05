times<-1:18
ages<-rep(0,18)
flip<-1.5

flipprob<-1/(1 + exp(-(((times)*flip-(ages*12+9))/(ages*12+9))*log(0.999/0.001)))

plot(times,flipprob)
abline(v=9)

lines(times,flips[1:18])

flipprob<-1/(1 + exp(-(((pop$infected[pop$infected>0]+pop$PIweighted[pop$infected>0])*flip-(pop$age[pop$infected>0]*12+9))/(pop$age[pop$infected>0]*12+9))*log(0.999/0.001)))

