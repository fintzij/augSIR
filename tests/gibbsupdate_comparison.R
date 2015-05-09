library(mcmc)

# set up log unnormalized posterior functions
logupost <- function(Xcount, W, theta, priors, initdist, popsize, tmax){    
    
    b = theta[1]
    m = theta[2]
    p = theta[3]
    
    # unpack priors
    beta.prior <- priors[[1]]
    mu.prior <- priors[[2]]
    p.prior <- priors[[3]]
    
    # set up log likelihood
    
    notpossible <- (b <=0 )| (m <= 0) | (p<=0 | p >= 1)
    
    if(notpossible == TRUE){
        
        lup <- -Inf
        
    } else{
        
        indend <- dim(Xcount)[1]
        
        numinf <- Xcount[,2]
        numsusc <- Xcount[,3]
        
        events <- diff(Xcount[,2], lag = 1)
        
        infec.rates <- (b * numinf) * numsusc
        recov.rates <- m * numinf
        
        hazards <- infec.rates + recov.rates
        
        rates <- ifelse(events==1, infec.rates, recov.rates)
        
        logl <- dbinom(sum(W[,2]), sum(W[,3]), prob=p, log=TRUE) + 
            dmultinom(c(Xcount[1,3], Xcount[1,2], 0), prob = initdist, log=TRUE) + 
            sum(log(rates[1:(indend - 1)])) - sum(hazards[1:(indend - 1)]*diff(Xcount[,1], lag = 1)) - 
            hazards[indend]*max(0,tmax - Xcount[indend,1])
        
        prior_lprob <- dgamma(b, shape = beta.prior[1], rate = beta.prior[2], log = TRUE) + dgamma(m, shape = mu.prior[1], rate = mu.prior[2], log = TRUE) +
            dbeta(p, p.prior[1], p.prior[2], log = TRUE)
        
        # return log unnormalized posterior
        
        lup <- ifelse(notpossible == TRUE, -Inf, logl + prior_lprob + dnorm(theta, c(0.005, 0.5, 0.25), c(0.001, 0.1, 0.1),log=TRUE))
        
    }    
    
    return(lup)
}

# function to get mcmc statistics from a batch

get.stats <- function(batch){
    
    # posterior means
    batch_mean <- apply(batch, 2, mean)
    
    # posterior second moments
    batch_mean_ofsq <- apply(batch^2, 2, mean) # second moment
    
    # posterior variances
    batch_sigmasq <- batch_mean_ofsq - batch_mean^2
    
    # monte carlo se for posterior mean
    mu_mcse <- apply(batch, 2, sd)
    
    # monte carlo se for posterior variance
    u <- batch
    v <- batch^2
    ubar <- apply(u, 2, mean)
    vbar <- apply(v, 2, mean)
    deltau <- sweep(u, 2, ubar)
    deltav <- sweep(v, 2, vbar)
    foo <- sweep(deltau, 2, ubar, "*")
    sigmasq_mcse <- sqrt(apply((deltav - 2 * foo)^2, 2, mean))
    
    list(batch_mean = batch_mean, batch_sigmasq = batch_sigmasq, mu_mcse = mu_mcse, sigmasq_mcse = sigmasq_mcse)
}

# simulate data and run mcmc

prior_list <- list(priors_1 = list(beta.prior = c(1, 200), mu.prior <- c(1, 0.001), p.prior <- c(1,1)), # moderately informative
                    priors_2 = list(beta.prior = c(1, 2), mu.prior <- c(1, 0.001), p.prior <- c(1,1)), # slightly informative
                    priors_3 = list(beta.prior = c(1e-4, 2e-2), mu.prior <- c(1, 0.001), p.prior <- c(1,1)), # diffuse
                    priors_4 = list(beta.prior = c(6.25e-6, 1.25e-3), mu.prior <- c(1, 0.001), p.prior <- c(1,1)), # more diffuse
                    priors_1 = list(beta.prior = c(4e-7, 8e-5), mu.prior <- c(1, 0.001), p.prior <- c(1,1))) # highly diffuse

prior_labels <- list("beta ~ Gamma(1, 200); mu ~ Gamma(1, 200), p ~ Beta(1,1)",
                      "beta ~ Gamma(1, 2); mu ~ Gamma(1, 2), p ~ Beta(1,1)",
                      "beta ~ Gamma(1e-4, 2e-2); mu ~ Gamma(1e-4, 2e-2), p ~ Beta(1,1)",
                      "beta ~ Gamma(6.25e-6, 1.25e-3); mu ~ Gamma(6.25e-6, 1.25e-3), p ~ Beta(1,1)",
                      "beta ~ Gamma(4e-7, 8e-5); mu ~ Gamma(4e-7, 8e-5), p ~ Beta(1,1)")

for(k in 1:5){
    
    pdf(file = paste("mcmc_gibbs_comparison_Prior",k,".pdf",sep=""))
    
    par(mfrow = c(3,2))
    
    for(j in 1:5){
        print(paste("Priorset",k," Dataset",j))
        
        # simulate data
        SIRres<-SIRsim(popsize = 200, initdist = c(0.95, 0.05, 0), b = 0.005, mu=0.5, a=0, tmax = 25, censusInterval=0.25, sampprob = 0.25, returnX = TRUE)
        
        if(max(SIRres$results[,3]) < 5){
            while(max(SIRres$results[,3]) < 5){
                SIRres<-SIRsim(popsize = 200, initdist = c(0.95, 0.05, 0), b = 0.005, mu=0.5, a=0, tmax = 25, censusInterval=0.25, sampprob = 0.25, returnX = TRUE)
                
            }
        }
        
        # build count matrix
        Xcount <- build_countmat(SIRres$trajectory, popsize = 200)
        
        # build observation matrix
        W <- as.matrix(data.frame(time = SIRres$results[,1], sampled = SIRres$results[,2], augmented = 0))
        
        W <- updateW(W, Xcount)
        
        
        # run mcmc and gibbs samplers
        
        theta.init <- c(0.004, 0.6, 0.2); results <- list()
        
        priors <- prior_list[[k]]
        beta.prior <- priors[[1]]
        mu.prior <- priors[[2]]
        p.prior <- priors[[3]]
        
        mcmc.out <- metrop(obj = logupost, initial = theta.init, nbatch = 300000, scale = c(0.0316^2, 0.1^2, 0.1^2), Xcount = Xcount, W=W, priors = priors, initdist = c(0.95, 0.05, 0), popsize=200, tmax = 25)
        
        rateparams <- replicate(50000, update_rates(Xcount = Xcount, beta.prior = beta.prior, mu.prior = mu.prior, alpha.prior = NULL, popsize = 200)[1:2])
        sampparam <- replicate(50000, update_prob(W = W, p.prior = p.prior))
        
        results[[j]] <- list(mcmcstats = c(get.stats(mcmc.out$batch), mcmc.out$accept), gibbsstats = get.stats(cbind(t(rateparams),sampparam)))
        
        # extract results
        mcmc_stats <- results[[j]][[1]]
        gibbs_stats <- results[[j]][[2]]
        
        # generate plots for beta parameter
        ts.plot(mcmc.out$batch[,1], main = paste("Metropolis - Beta - Prior",k," - ", "Dataset",j,sep=""),xlab = prior_labels[[k]])
        abline(a = mcmc_stats[[1]][1], b = 0, col = "red", lwd = 2)
        abline(a = c(mcmc_stats[[1]][1] - 1.96 * mcmc_stats[[3]][1]), b = 0, col = "red", lwd = 2)
        abline(a = mcmc_stats[[1]][1] + 1.96 * mcmc_stats[[3]][1], b = 0, col = "red", lwd = 2)
        
        ts.plot(rateparams[1,], main = paste("Gibbs - Beta - Prior",k," - ", "Dataset",j,sep=""),xlab = prior_labels[[k]])
        abline(a = gibbs_stats[[1]][1], b = 0, col = "red", lwd = 2)
        abline(a = c(gibbs_stats[[1]][1] - 1.96 * gibbs_stats[[3]][1]), b = 0, col = "red", lwd = 2)
        abline(a = gibbs_stats[[1]][1] + 1.96 * gibbs_stats[[3]][1], b = 0, col = "red", lwd = 2)
        
        plot(density(mcmc.out$batch[,1]), main = paste("Metropolis - Beta - Prior",k," - ", "Dataset",j,sep=""),xlab = prior_labels[[k]]) 
        
        plot(density(rateparams[1,]), main = paste("Gibbs - Beta - Prior",k," - ", "Dataset",j,sep=""),xlab = prior_labels[[k]])
        
        
        # generate plots for mu paramete
        ts.plot(mcmc.out$batch[,2], main = paste("Metropolis - Mu - Prior",k," - ", "Dataset",j,sep=""),xlab = prior_labels[[k]])
        abline(a = mcmc_stats[[1]][2], b = 0, col = "red", lwd = 2)
        abline(a = c(mcmc_stats[[1]][2] - 1.96 * mcmc_stats[[3]][2]), b = 0, col = "red", lwd = 2)
        abline(a = mcmc_stats[[1]][2] + 1.96 * mcmc_stats[[3]][2], b = 0, col = "red", lwd = 2)
        
        ts.plot(rateparams[2,], main = paste("Gibbs - Mu - Prior",k," - ", "Dataset",j,sep=""),xlab = prior_labels[[k]])
        abline(a = gibbs_stats[[1]][2], b = 0, col = "red", lwd = 2)
        abline(a = c(gibbs_stats[[1]][2] - 1.96 * gibbs_stats[[3]][2]), b = 0, col = "red", lwd = 2)
        abline(a = gibbs_stats[[1]][2] + 1.96 * gibbs_stats[[3]][2], b = 0, col = "red", lwd = 2)
        
        plot(density(mcmc.out$batch[,2]), main = paste("Metropolis - Mu - Prior",k," - ", "Dataset",j,sep=""),xlab = prior_labels[[k]])
        
        plot(density(rateparams[2,]), main = paste("Gibbs - Mu - Prior",j,sep=""),xlab = prior_labels[[k]])
        
        
        # generate plots for binomial parameter
        ts.plot(mcmc.out$batch[,3], main = paste("Metropolis - Sampling prob. - Prior",k," - ", "Dataset",j,sep=""),xlab = prior_labels[[k]])
        abline(a = mcmc_stats[[1]][3], b = 0, col = "red", lwd = 2)
        abline(a = c(mcmc_stats[[1]][3] - 1.96 * mcmc_stats[[3]][3]), b = 0, col = "red", lwd = 2)
        abline(a = mcmc_stats[[1]][3] + 1.96 * mcmc_stats[[3]][3], b = 0, col = "red", lwd = 2)
        
        ts.plot(sampparam, main = paste("Gibbs - Sampling prob. - Prior",k," - ", "Dataset",j,sep=""),xlab = prior_labels[[k]])
        abline(a = gibbs_stats[[1]][3], b = 0, col = "red", lwd = 2)
        abline(a = c(gibbs_stats[[1]][3] - 1.96 * gibbs_stats[[3]][3]), b = 0, col = "red", lwd = 2)
        abline(a = gibbs_stats[[1]][3] + 1.96 * gibbs_stats[[3]][3], b = 0, col = "red", lwd = 2)
        
        plot(density(mcmc.out$batch[,3]), main = paste("Metropolis - Sampling prob. - Prior",k," - ", "Dataset",j,sep=""),xlab = prior_labels[[k]])
        
        plot(density(sampparam), main = paste("Gibbs - Sampling prob. - Prior",j,sep=""),xlab = prior_labels[[k]])        
        
    }
    dev.off()
    
}

