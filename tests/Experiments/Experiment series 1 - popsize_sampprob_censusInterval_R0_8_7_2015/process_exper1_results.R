library(ggplot2)
library(dplyr)
library(gridExtra)

setwd("C:/Users/Jonathan/Google Drive/UW/Year 3 +/Dissertation/Code/augSIR/tests/convergence_comp/Experiment series 4 - reparameterized model/Results data")

# Collect data for each R0 ------------------------------------------------

files <- list.files()

simsettings <- expand.grid(simnum = 1:3, popsize = c(50, 150, 300), censusInterval = c(0.05, 0.2), samp_prob = c(0.1, 0.4), R0 = c(4, 8))

results <- cbind(simsettings,data.frame(loglik = NA, traj_accepts = NA, R0_accepts = NA, Mu_accepts = NA, Samp_prob_accepts = NA, p_init_accepts = NA, time = NA))

datafiles <- list.files()

for(k in 1:nrow(results)) {
    filepattern <- paste(paste("augSIR_",results$popsize[k],results$censusInterval[k],results$samp_prob[k],results$R0[k],results$simnum[k], sep = "_"), ".Rdata", sep = "")
    objpattern <-  paste("augSIR_",results$popsize[k],results$censusInterval[k],results$samp_prob[k],results$R0[k],results$simnum[k], sep = "_")
    
    if(! filepattern %in% datafiles) next
    
    load(filepattern)
    
    dat <- eval(parse(text = objpattern))
    
    eval(parse(text = paste("rm(",objpattern,")")))
    
    results$time[k] <- ifelse(dat[[1]] < 8, dat[[1]] * 60, dat[[1]])
    results$loglik[k] <- mean(dat[[2]][,"loglik"])
    results$traj_accepts[k] <- mean(dat[[2]][,"traj_accepts"])/results$popsize[k]
    results$R0_accepts[k] <- mean(dat[[2]][,"R0_accepts"])/results$popsize[k]
    results$Mu_accepts[k] <- mean(dat[[2]][,"Mu_accepts"])/results$popsize[k]
    results$Samp_prob_accepts[k] <- mean(dat[[2]][,"Samp_prob_accepts"])/results$popsize[k]
    results$p_init_accepts[k] <- mean(dat[[2]][,"p_init_accepts"])/results$popsize[k]
    
    
    
    # to plot a set of trajectories
    popsize <- results$popsize[k]; censusInterval <- results$censusInterval[k]; R0 = results$R0[k]; b = results$R0[k]/results$popsize[k]; m = 1; tmax = 10; initdist <- c(0.95, 0.05, 0); samp_prob <- results$samp_prob[k]; initialization_num <- results$simnum[k]
    
    
    set.seed(183427)
    SIRres <- SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval = censusInterval, sampprob = samp_prob, returnX = TRUE)
    if(all(SIRres$results$Observed == 0)){
        while(all(SIRres$results$Observed == 0)){
            SIRres <- SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval = censusInterval, sampprob = samp_prob, returnX = TRUE)
        }
    }
    
    truetrajec <- data.frame(time = unique(SIRres$trajectory[,1]), infected = c(sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]), sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]) + cumsum(SIRres$trajectory[SIRres$trajectory[,1]!=0,3])), simnum = 0)
    
    trajectories <- dat[[3]]
    trajectories2 <- list(); observations2 <- list()
    for(k in 1:length(dat[[3]])){
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
    
    assign(paste(objpattern,"_plot",sep=""), ggplot(data = trajecs, aes(x = time, y = infected, group = simnum)) + geom_path(alpha = 0.3) + geom_path(data = truetrajec, aes(x = time, y = infected),colour="red", size = 1)  + geom_point(data=data.frame(SIRres$results,simnum=0), aes(x=time,y=Observed),size=1,colour="blue") + theme_bw() + scale_x_continuous(limits = c(0, tmax)) + labs(title = objpattern))
    
}

plotcombns <- expand.grid(popsize = c(50, 150, 300), censusInterval = c(0.05, 0.2), samp_prob = c(0.1, 0.4), R0 = c(4, 8))
for(k in 1:nrow(plotcombns)){
    pdf(paste("augSIR_",plotcombns$popsize[k],plotcombns$censusInterval[k], plotcombns$samp_prob[k], plotcombns$R0[k],"plot.pdf",sep = "_"))
    plotnames <- paste("augSIR_",plotcombns$popsize[k],plotcombns$censusInterval[k], plotcombns$samp_prob[k], plotcombns$R0[k],1:3,"plot",sep = "_")
    plotnames <- plotnames[sapply(plotnames, exists)]
    
    eval(parse(text = paste("grid.arrange(",paste(plotnames,collapse=","),")")))
    dev.off()
}

# # R0 vs. accepts
pdf(file = "convcomp_popsize_results.pdf")
for(R0 in c(4, 8)){
    dat <- results[results$R0 == R0, ]
    
    print(ggplot(dat, aes(x = popsize, y = accepts, colour = as.factor(samp_prob))) + geom_point(size = 3)+ theme_bw() + labs(title = paste("Popsize vs. proportion of trajectories accepted; R0 =", R0)))
    
    print(ggplot(dat, aes(x = popsize, y = loglik, colour = as.factor(samp_prob))) + geom_point(size = 3)+ theme_bw() + labs(title = paste("Popsize vs. proportion of trajectories accepted; R0 =", R0)))
    
    print(ggplot(dat, aes(x = popsize, y = time, colour = as.factor(samp_prob))) + geom_point(size = 3)+ theme_bw() + labs(title = paste("Popsize vs. proportion of trajectories accepted; R0 =", R0)))
    
}
dev.off()

