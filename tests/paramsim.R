### Simulation to test parameter update functions to see if they recover the true values of the parameters
popsize = 200
beta.prior = c(.012, 1.1)
mu.prior = c(0.96, 1.96)
p.prior = c(0.022, 0.084)

param.res <- matrix(nrow = 1000, ncol = 3)

for(j in 1:1000){
    print(j)
    
    SIRres<-SIRsim(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 20, censusInterval=.25, sampprob = 0.25, returnX = TRUE, binomsamp = TRUE)
    
    if(dim(SIRres$results)[1] < 20){
        while(dim(SIRres$results)[1] < 20){
            SIRres<-SIRsim(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 20, censusInterval=.25, sampprob = 0.25, returnX = TRUE, binomsamp = TRUE)
            
        }
    }
    
    param.res[j,1:2] <- apply(replicate(2001, update_rates(SIRres$trajectory, beta.prior, mu.prior, popsize =200))[1:2,],1,median)
    param.res[j,3] <-  median(replicate(2001, update_prob(SIRres$results, p.prior)))    
     
}

