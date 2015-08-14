library(ggplot2)
library(dplyr)
library(gridExtra)

setwd("C:/Users/Jonathan/Google Drive/UW/Year 3 +/Dissertation/Code/augSIR/tests/Experiments/Experiment series 1 - popsize_sampprob_censusInterval_R0_8_7_2015/Data")

# Collect data for each R0 ------------------------------------------------

files <- list.files()

simsettings <- expand.grid(simnum = 1, popsize = c(50, 200, 500), censusInterval = c(0.05, 0.2, 0.4), samp_prob = c(0.1, 0.4), R0 = c(4, 8))

results <- cbind(simsettings,data.frame(loglik = NA, traj_accepts = NA, R0_accepts = NA, Mu_accepts = NA, Samp_prob_accepts = NA, p_init_accepts = NA, time = NA))

datafiles <- list.files()

for(s in 1:nrow(results)) {
    print(s)
    filepattern <- paste(paste("augSIR_",results$popsize[s],results$censusInterval[s],results$samp_prob[s],results$R0[s],results$simnum[s], sep = "_"), ".Rdata", sep = "")
    objpattern <-  paste("augSIR_",results$popsize[s],results$censusInterval[s],results$samp_prob[s],results$R0[s],results$simnum[s], sep = "_")
    
    if(! filepattern %in% datafiles) break
    
    load(filepattern)
    
    dat <- eval(parse(text = objpattern))
    
    eval(parse(text = paste("rm(",objpattern,")")))
    
    results$time[s] <- ifelse(dat[[1]] < 8, dat[[1]] * 60, dat[[1]])
    results$loglik[s] <- mean(dat[[2]][,"loglik"])
    results$traj_accepts[s] <- mean(dat[[2]][,"traj_accepts"])/results$popsize[s]
    results$R0_accepts[s] <- mean(dat[[2]][,"R0_accepts"])/results$popsize[s]
    results$Mu_accepts[s] <- mean(dat[[2]][,"Mu_accepts"])/results$popsize[s]
    results$Samp_prob_accepts[s] <- mean(dat[[2]][,"Samp_prob_accepts"])/results$popsize[s]
    results$p_init_accepts[s] <- mean(dat[[2]][,"p_init_accepts"])/results$popsize[s]
    
    
    # to plot a set of trajectories
    popsize <- results$popsize[s]; censusInterval <- results$censusInterval[s]; R0 = results$R0[s]; b = results$R0[s]/results$popsize[s]; m = 1; tmax = 10; initdist <- c(0.95, 0.05, 0); samp_prob <- results$samp_prob[s]; initialization_num <- results$simnum[s]
    
    
    set.seed(183427)
    SIRres <- SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval = censusInterval, sampprob = samp_prob, returnX = TRUE)
    if(all(SIRres$results$Observed == 0)){
        while(all(SIRres$results$Observed == 0)){
            SIRres <- SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval = censusInterval, sampprob = samp_prob, returnX = TRUE)
        }
    }
    
    truetrajec <- data.frame(time = unique(SIRres$trajectory[,1]), infected = c(sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]), sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]) + cumsum(SIRres$trajectory[SIRres$trajectory[,1]!=0,3])), simnum = 0)
    
    trajectories <- dat[[3]]
    trajectories2 <- list(); #observations2 <- list()
    for(j in 1:length(dat[[3]])){
        if ((j%%1)==0){
            traj <- trajectories[[j]]
            Xobs <- data.frame(time = traj[,1], 
                               infected = traj[,2], 
                               simnum = j)
            trajectories2[[j]] <- Xobs
            
#             obs <- data.frame(time = seq(0,max(Xobs$time),by=censusInterval), 
#                               truth = 0, 
#                               sampled = 0, 
#                               simnum = j)
#             
#             for(t in 1:dim(obs)[1]){
#                 obs$truth[t] <- traj[traj[,t] <= obs[j,t],2][sum(traj[,t] <= obs[j,t])]
#             }
#             
#             obs$sampled <- rbinom(dim(obs)[1], obs$truth, prob = samp_prob)
#             
#             observations2[[j]] <- obs
        }
        
    }
    trajecs <- data.frame(do.call(rbind,trajectories2))
    
    
    params <- data.frame(dat[[2]])
    loglik_plot <- ggplot(params, aes(x=seq(1,nrow(dat[[2]])), y=loglik)) + geom_line() + labs(title = "log likelihood", xlab = "Iteration", ylab = "log-likelihood")
    R0_plot <- ggplot(params, aes(R0)) + geom_density()+ labs(title = paste("R0; % accepted = ",results$R0_accepts[s],sep="")) + geom_vline(xintercept = results$R0[s], colour="red")
    mu_plot <- ggplot(params, aes(Mu)) +  geom_density() + labs(title = paste("Mu; % accepted = ",results$Mu_accepts[s],sep="")) + geom_vline(xintercept = 1, colour="red")
    probs_plot <- ggplot(params, aes(probs)) + geom_density() + labs( main = paste("Binomial probability; % accepted = ",results$Samp_prob_accepts[s],sep="")) + geom_vline(xintercept = results$samp_prob[s], colour="red")
    pinit_plot <- ggplot(params, aes(p_initinfec)) + geom_density() + labs( main = paste("Prob of infection at time 0; % accepted = ",results$Samp_prob_accepts[s],sep="")) + geom_vline(xintercept = 0.05, colour="red")

    
    assign(paste(objpattern,"_params",sep=""),grid.arrange(loglik_plot, R0_plot, mu_plot, probs_plot, pinit_plot, nrow=5))
    
    assign(paste(objpattern,"_plot",sep=""), grid.arrange((ggplot(data = trajecs, aes(x = time, y = infected, group = simnum)) + geom_path(alpha = 0.3) + geom_path(data = truetrajec, aes(x = time, y = infected),colour="red", size = 1)  + geom_point(data=data.frame(SIRres$results,simnum=0), aes(x=time,y=Observed),size=1,colour="blue") + theme_bw() + scale_x_continuous(limits = c(0, tmax)) + labs(title = objpattern)), eval(parse(text = paste(objpattern,"_params",sep=""))), ncol = 2))
}

plotcombns <- expand.grid(popsize = c(50, 200, 500), censusInterval = c(0.05, 0.2, 0.4), samp_prob = c(0.1, 0.4), R0 = c(4, 8))
for(k in 1:nrow(plotcombns)){
    pdf(paste("augSIR_",plotcombns$popsize[k],plotcombns$censusInterval[k], plotcombns$samp_prob[k], plotcombns$R0[k],"plot.pdf",sep = "_"))
    plotnames <- paste("augSIR_",plotcombns$popsize[k],plotcombns$censusInterval[k], plotcombns$samp_prob[k], plotcombns$R0[k],1,"plot",sep = "_")
    plotnames <- plotnames[sapply(plotnames, exists)]
    
    eval(parse(text = paste("grid.arrange(",paste(plotnames,collapse=","),")")))
    dev.off()
}