#### Comparison with the pomp package
library(pomp)




#### Comparison with the GillespieSSA package
library(GillespieSSA)

# initialize ssa settings
initcounts <- rmultinom(1, 200, prob = c(0.95, 0.05, 0)); x0 = c(S = initcounts[1], I = initcounts[2], R = initcounts[3])
parms <- c(b = 0.0075, m = 0.5)
a = c("b*S*I", "m*I")
nu <- matrix(c(-1,+1,0, 0, -1, +1),nrow=3)

# initialize doPar
cl <- makeCluster(2)
registerDoParallel(cl)

niter <- 250000; tmax = 30

# run simulation
mytrajecs <- foreach(t = 1:niter, .packages=c('augSIR', 'GillespieSSA')) %dopar%{
    # call the function
    SIRres<-SIRsim(popsize = 200, initdist = c(0.95, 0.05, 0), b = 0.0075, mu = 0.5, a=0, tmax = tmax, censusInterval=0.25, sampprob = 0.25, trim = FALSE, returnX = FALSE)
    
    # If the epidemic dies out, keep calling the function
    if(max(SIRres$Truth) < 5){
        while(max(SIRres$Truth) < 5){
            SIRres<-SIRsim(popsize = 200, initdist = c(0.95, 0.05, 0), b = 0.0075, mu=0.5, a=0, tmax = tmax, censusInterval=0.25, sampprob = 0.25, trim = FALSE, returnX = FALSE)
            
        }
    }
    
    SIRres <- SIRres[,c(1,3)]
    
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
    
#     Build the matrix with counts of number of infecteds and susceptibles at each event time
#     Xcount.cur <- build_countmat(X = X.cur, popsize = 200)
#     

    # Build matrix with counts at observation times
#     W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$Observed, augmented = 0))
#     W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
    
#     # calculate log-likelihood for SIR alone, and log-likelihood for trajectory and data
#     pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = 0.005, m = 0.5, a = 0, initdist = c(0.95, 0.05,0), popsize = 200)
#     
#     pop_prob.cur

    initcounts <- rmultinom(1, 200, prob = c(0.95, 0.05, 0))
    
    if(initcounts[2] == 0){
        while(initcounts[2] == 0){
            initcounts <- rmultinom(1, 200, prob = c(0.95, 0.05, 0))
            
        }
    }

    x0 = c(S = initcounts[1], I = initcounts[2], R = initcounts[3])
    
    gSSA_res <- ssa(x0, a, nu, parms, tf = 25)$data[,c(1,3,2)]
    
    # Build matrix with counts at observation times for the GillespieSSA matrix
    gSSA <- as.matrix(data.frame(time = seq(0,tmax,by=0.25), Truth= 0))
    
    for(s in 1:dim(gSSA)[1]){
        gSSA[s,2] <- gSSA_res[sum(gSSA_res[,1] <= gSSA[s,1]), 2]
    }

    list(SIRres, gSSA)
    
}
stopCluster(cl)

save(mytrajecs, file="mytrajecs")


# Plot results
obstimes <- seq(0,tmax,by=0.25)

trajecs_SIRsim <- matrix(0, nrow = length(obstimes), ncol = niter+1); trajecs_SIRsim[,1] <- obstimes
trajecs_SSA <- matrix(0, nrow = length(obstimes), ncol = niter+1); trajecs_SSA[,1] <- obstimes

# extract counts of numbers of infecteds

for(s in 1:niter){
    trajecs_SIRsim[,s+1] <- mytrajecs[[s]][[1]][,2]
    trajecs_SSA[,s+1] <- mytrajecs[[s]][[2]][,2]
}

# construct matrices to compare simulation methods

# compute the mean, 25th and 75th quantiles, and the standard deviation
comp.summary <- data.frame(Method = c(rep("SIRsim", length(obstimes)), rep("GillespieSSA", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Median_Infected = 0, 
                           Q25_Infected = 0,
                           Q75_Infected = 0,
                           mcsd = 0, 
                           mean = 0,
                           lowermcci = 0,
                           uppermcci = 0)

for(n in 1:length(obstimes)){

    comp.summary[n, 3:5] <- quantile(trajecs_SIRsim[n, 2:dim(trajecs_SIRsim)[2]], probs = c(0.5, 0.25, 0.75))
    comp.summary[n + length(obstimes), 3:5] <- quantile(trajecs_SSA[n, 2:dim(trajecs_SSA)[2]], probs = c(0.5, 0.25, 0.75))
    
    comp.summary[n, 6] <- sd(trajecs_SIRsim[n, 2:dim(trajecs_SIRsim)[2]])/sqrt(niter)
    comp.summary[n + length(obstimes), 6] <- sd(trajecs_SSA[n, 2:dim(trajecs_SSA)[2]])/sqrt(niter)
    
    comp.summary[n, 7] <- mean(trajecs_SIRsim[n, 2:dim(trajecs_SIRsim)[2]])
    comp.summary[n + length(obstimes), 7] <- mean(trajecs_SSA[n, 2:dim(trajecs_SSA)[2]])
    
    comp.summary[n, 8] <- comp.summary[n,7] - 2*comp.summary[n,6]
    comp.summary[n + length(obstimes), 8] <- comp.summary[n+ length(obstimes),7] - 1.96*comp.summary[n+ length(obstimes),6]
    
    comp.summary[n, 9] <- comp.summary[n,7] + 2*comp.summary[n,6]
    comp.summary[n + length(obstimes), 9] <- comp.summary[n+ length(obstimes),7] + 1.96*comp.summary[n+ length(obstimes),6]
    
}

comp.summary$lower_demean <- comp.summary$lowermcci - comp.summary$mean
comp.summary$upper_demean <- comp.summary$uppermcci - comp.summary$mean



ggplot(comp.summary, aes(x = Time, y = Median_Infected, colour = Method)) + geom_line() + geom_ribbon(data=comp.summary, aes(x=Time, ymin = Q25_Infected, ymax = Q75_Infected, fill=Method), alpha=0.25) +
    labs(y="Count", title = "Median, 25th and 75th Quantiles")

ggplot(comp.summary, aes(x = Time, y = mean, colour = Method)) + geom_line() + geom_ribbon(data=comp.summary, aes(x=Time, ymin = lowermcci, ymax = uppermcci, fill=Method), alpha=0.25) +
    labs(y="Count",  title = "Mean, 95% Monte Carlo CI")

ggplot(comp.summary, aes(x = Time, y=sd, colour = Method)) + geom_line() + labs(y="", title = "Standard deviation")

ggplot(comp.summary, aes(x=Time,fill=Method)) + geom_ribbon(aes(ymin=lower_demean, ymax=upper_demean), alpha=0.25)

ggplot(subset(comp.summary, Time = 2.5), aes(x = Time, fill = method)) + (data = subset(comp.summary, Time = 2.5), aes(x=Time, ymin = lowermcci, ymax = uppermcci, fill=Method), alpha=0.25)

write.csv(comp.summary, file = "comp.summary.csv")
