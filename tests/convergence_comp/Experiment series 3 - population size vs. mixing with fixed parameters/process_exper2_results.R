library(ggplot2)
library(dplyr)

setwd("C:/Users/Jonathan/Google Drive/UW/Year 3 +/Dissertation/Code/augSIR/tests/convergence_comp/Experiment series 2 - population size and binomial prob - assessing mixing/results files")


# Collect data for each R0 ------------------------------------------------

files <- list.files()

simsettings <- expand.grid(simnum = 1:3, popsize = c(50, 100, 150, 200, 300, 400, 500), samp_prob = c(0.05, 0.2, 0.5), R0 = c(2,5,10))

results <- cbind(simsettings,data.frame(loglik = NA, accepts = NA, time = NA, beta = NA, mu = NA, p = NA))

for(k in 1:nrow(results)) {
    filepattern <- paste("augSIR__",results$popsize[k],"_0.2_",results$samp_prob[k],"_",results$R0[k],"_",results$simnum[k],".Rdata", sep = "")
    objpattern <-  paste("augSIR__",results$popsize[k],"_0.2_",results$samp_prob[k],"_",results$R0[k],"_",results$simnum[k], sep = "")
    
    load(filepattern)
    
    dat <- eval(parse(text = objpattern))
    
    eval(parse(text = paste("rm(",objpattern,")")))
    
    results$time[k] <- ifelse(dat[[1]] < 8, dat[[1]] * 60, dat[[1]])
    results$loglik[k] <- mean(dat[[2]][,"loglik"])
    results$accepts[k] <- mean(dat[[2]][,"accepts"])/results$popsize[k]
    
    results$beta[k] <- mean(dat[[2]][,"Beta"])
    results$mu[k] <- mean(dat[[2]][,"Mu"])
    results$p[k] <- mean(dat[[2]][,"probs"])
    
}

# # R0 vs. accepts
# ggplot(results[results$R0 == 5,], aes(x = accepts, y = beta*popsize, colour = as.factor(popsize))) + geom_point() + labs(y = "R0")

pdf(file = "convcomp_popsize_results.pdf")
for(R0 in c(2,5,10)){
    dat <- results[results$R0 == R0, ]
    
    print(ggplot(dat, aes(x = popsize, y = accepts, colour = as.factor(samp_prob))) + geom_point(size = 3)+ theme_bw() + labs(title = paste("Popsize vs. proportion of trajectories accepted; R0 =", R0)))
    
    print(ggplot(dat, aes(x = popsize, y = loglik, colour = as.factor(samp_prob))) + geom_point(size = 3)+ theme_bw() + labs(title = paste("Popsize vs. proportion of trajectories accepted; R0 =", R0)))
    
    print(ggplot(dat, aes(x = popsize, y = time, colour = as.factor(samp_prob))) + geom_point(size = 3)+ theme_bw() + labs(title = paste("Popsize vs. proportion of trajectories accepted; R0 =", R0)))
    
}
dev.off()

# to plot a set of trajectories
popsize <- 50; censusInterval <- 0.05; R0 = 5; b = R0/popsize; m = 1; tmax = 10; initdist <- c(0.95, 0.05, 0); samp_prob <- 0.05; initialization_num <- 1


set.seed(183427)
SIRres <- SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval = censusInterval, sampprob = samp_prob, returnX = TRUE)
if(all(SIRres$results$Observed == 0)){
    while(all(SIRres$results$Observed == 0)){
        SIRres <- SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval = censusInterval, sampprob = samp_prob, returnX = TRUE)
    }
}

truetrajec <- data.frame(time = unique(SIRres$trajectory[,1]), infected = c(sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]), sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]) + cumsum(SIRres$trajectory[SIRres$trajectory[,1]!=0,3])), simnum = 0)

trajectories <- res[[3]]
trajectories2 <- list(); observations2 <- list()
for(k in 1:length(res[[3]])){
    if ((k%%1)==0){
        traj <- trajectories[[k]]
        Xobs <- data.frame(time = traj[,1], 
                           infected = traj[,2], 
                           simnum = k)
        trajectories2[[k]] <- Xobs
        
        obs <- data.frame(time = seq(0,max(Xobs$time),by=censusInterval), 
                          truth = 0, 
                          sampled = 0, 
                          simnum = k)
        
        for(j in 1:dim(obs)[1]){
            obs$truth[j] <- traj[traj[,1] <= obs[j,1],2][sum(traj[,1] <= obs[j,1])]
        }
        
        obs$sampled <- rbinom(dim(obs)[1], obs$truth, prob = samp_prob)
        
        observations2[[k]] <- obs
    }
    
}
trajecs <- data.frame(do.call(rbind,trajectories2))

print(ggplot(data = trajecs, aes(x = time, y = infected, group = simnum)) + geom_path(alpha = 0.3) + geom_path(data = truetrajec, aes(x = time, y = infected),colour="red", size = 1)  + geom_point(data=data.frame(SIRres$res,simnum=0), aes(x=time,y=Observed),size=4,colour="blue") + theme_bw() + scale_x_continuous(limits = c(0, tmax)))
