######################################################################
### eval.R contains commands to test functionality of the package. ###
######################################################################
library(profr)
library(proftools)
# augSIR simulation -------------------------------------------------------

# simulate data
SIRres<-SIRsim(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 20, censusInterval=0.25, sampprob = 0.25, returnX=TRUE)
if(dim(SIRres$results)[1] < 10){
    while(dim(SIRres$results)[1] < 10){
        SIRres<-SIRsim(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 200, censusInterval=.25, sampprob = 0.25, returnX = TRUE)
        
    }
}

# get data 
dat <- SIRres$results
dat.m <- melt(dat,id.vars="time")

ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + theme_bw()

sim.settings <- list(popsize = 200,
                     tmax = 20,
                     niter = 500,
                     amplify = 5,
                     initdist = c(0.995, 0.005, 0))

inits <- list(beta.init = 0.01 + runif(1,-0.005, 0.005),
              mu.init = 0.5 + runif(1, -0.05, 0.05),
              alpha.init = 0, 
              probs.init = 0.2 + runif(1,-0.1, 0.1))

priors <- list(beta.prior = c(.012, 1.1),
               mu.prior = c(0.96, 1.96),
               alpha.prior = NULL,
               p.prior = c(0.022, 0.084))

# run sampler

results <- augSIR(dat, sim.settings, priors, inits, returnX=TRUE)

# #profile
# Rprof("~/School/UW/Year 3 +/Dissertation/Code/augSIR/tests/profile/augSIRprofile.out")
# results <- augSIR(dat, sim.settings, priors, inits)
# Rprof()
# summaryRprof("~/School/UW/Year 3 +/Dissertation/Code/augSIR/tests/profile/augSIRprofile.out")
# # results.prof <- lineprof(augSIR(dat, sim.settings, priors, inits))
# plotProfileCallGraph(readProfileData("~/School/UW/Year 3 +/Dissertation/Code/augSIR/tests/profile/augSIRprofile.out"),score = "total")

# Plots -------------------------------------------------------------------

censusInterval <- 0.25; p <- 0.2
trajectories <- list(); observations <- list(); likelihoods <- list()

for(k in 1:(length(results[[4]]))){
    if ((k%%10)==0){
        traj <- results[[4]][[k]]
        Xobs <- data.frame(time = unique(traj[,1]), 
                           infected = c(sum(traj[traj[,1]==0,3]),sum(traj[traj[,1]==0,3]) + cumsum(traj[traj[,1]!=0,3])), 
                           simnum = k)
        trajectories[[k]] <- Xobs
        
        obs <- data.frame(time = seq(0,max(Xobs$time),by=0.25), 
                          truth = 0, 
                          sampled = 0, 
                          simnum = k)
        
        for(j in 1:dim(obs)[1]){
            obs$truth[j] <- sum(traj[traj[,1]<=obs$time[j],3])
        }
        
        obs$sampled <- rbinom(dim(obs)[1], obs$truth, prob = 0.2)
        
        observations[[k]] <- obs
    }
    
}

trajecs <- do.call(rbind,trajectories)
samples <- do.call(rbind,observations)
likelihoods <- do.call(rbind,likelihoods)

truetrajec <- data.frame(time = unique(SIRres$trajectory[,1]), 
                         infected = c(sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]), sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]) + cumsum(SIRres$trajectory[SIRres$trajectory[,1]!=0,3])),
                         simnum = 0)

trajecs.gg <- ggplot(data = trajecs, aes(x = time, y = infected, group = simnum)) + geom_line(alpha = 0.1) + geom_line(data = truetrajec, aes(x = time, y = infected),colour="red", size = 1) + 
    geom_boxplot(data = samples, aes(x=time, y=sampled, group= time),outlier.shape=NA) + geom_point(data=data.frame(dat,simnum=0), aes(x=time,y=Observed),size=4,colour="blue") + theme_bw()
