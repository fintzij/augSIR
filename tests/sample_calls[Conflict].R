#### Call to SIR simulation function
library(GillespieSSA)


initcounts <- rmultinom(1, 200, prob = c(0.95, 0.05, 0)); x0 = c(S = initcounts[1], I = initcounts[2], R = initcounts[3])
parms <- c(b = 0.0075, m = 0.5)
a = c("b*S*I", "m*I")
nu <- matrix(c(-1,+1,0, 0, -1, +1),nrow=3)
SSAres <- ssa(x0, a, nu, parms, tf = 25)

cl <- makeCluster(3)
registerDoParallel(cl)

mytrajecs <- foreach(t = 1:250000, .packages=c('augSIR', 'GillespieSSA')) %dopar%{
    # call the function
    SIRres<-SIRsim(popsize = 200, initdist = c(0.95, 0.05, 0), b = 0.0075, mu = 0.5, a=0, tmax = 30, censusInterval=0.25, sampprob = 0.25, binomsamp = TRUE, returnX = TRUE)
    
    # If the epidemic dies out, keep calling the function
    if(dim(SIRres$results)[1] < 50){
        while(dim(SIRres$results)[1] < 50){
            SIRres<-SIRsim(popsize = 200, initdist = c(0.95, 0.05, 0), b = 0.0075, mu=0.5, a=0, tmax = 30, censusInterval=0.25, sampprob = 0.25, binomsamp = TRUE, returnX = TRUE)
            
        }
    }
    
#     
# #     # get data 
#     dat <- SIRres$results
#     dat.m <- melt(dat,id.vars="time")
# #     
#     ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + theme_bw()
#     
    
    # matrix with event times, subject id and event codes. 
    # Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
#     X.cur <- SIRres$trajectory
    
    # Build the matrix with counts of number of infecteds and susceptibles at each event time
#     Xcount.cur <- build_countmat(X = X.cur, popsize = 200)
    
#     popsize - min(Xcount.cur[,3])
#     # Build matrix with counts at observation times
#     W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$Observed, augmented = 0))
#     W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
#     
#     # calculate log-likelihood for SIR alone, and log-likelihood for trajectory and data
#     pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = max(W.cur[,1])+1, b = 0.005, m = 0.5, a = 0, initdist = c(0.95, 0.05,0), popsize = 200)
#     
#     pop_prob.cur

    initcounts <- rmultinom(1, 200, prob = c(0.95, 0.05, 0))
    
    if(initcounts[2] == 0){
        while(initcounts[2] == 0){
            initcounts <- rmultinom(1, 200, prob = c(0.95, 0.05, 0))
            
        }
    }

    x0 = c(S = initcounts[1], I = initcounts[2], R = initcounts[3])
    list(SIRres$trajectory,ssa(x0, a, nu, parms, tf = 25))
    
}
stopCluster(cl)

save(parlist, file="parlist.Rdata")
