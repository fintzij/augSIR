library(batch)
library(ggplot2)
source("augSIR.R")
source("auxilliary_functions.R")
source("SIRsimulation.R")
source("rjmcmc_functions.R")
source("matrix_build_update.R")
source("metropolis_hastings_functions.R")
source("path_sampling_functions.R")

# args <- commandArgs(TRUE)
# print(args)
# 
# popsize <- eval(parse(text = args[1]))
# censusInterval <- eval(parse(text = args[2]))
# samp_prob <- eval(parse(text = args[3]))
# resample_prop <- eval(parse(text = args[4]))
# initialization_num <- eval(parse(text = args[5]))

parseCommandArgs()

tmax = 10
b <- R0/popsize
m <- 1
# samp_prob <- 0.5 
initdist <- c(0.95, 0.05, 0)
# insert.prob = 1/3; remove.prob = 1/3; shift.prob = 1/3
# shift.int <- 0.1
niter <- 10000
trajecs_every <- 1
gillespie_trajecs <- vector("list", length = niter)
gillespie_extended <- vector("list", length = niter)

obstimes <- seq(0, tmax, by = censusInterval)

set.seed(initialization_num + 9235)
SIRres <- SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval = censusInterval, sampprob = samp_prob, returnX = TRUE, trim = FALSE)
if(sum(SIRres$results$Observed) < 5){
    while(sum(SIRres$results$Observed) < 5){
        SIRres <- SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval = censusInterval, sampprob = samp_prob, returnX = TRUE, trim = FALSE)
    }
}

gillespie_trajecs[[1]] <- build_countmat(SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval = censusInterval, sampprob = samp_prob, returnX = TRUE, trim = FALSE)$trajectory, popsize)
gillespie_extended[[1]] <- build_countmat(SIRsim_extended(popsize = popsize, initdist = initdist, b = b, m = m, obstimes = obstimes, samp_prob = samp_prob, trim = FALSE)$trajectory, popsize)

# initdist.shift <- runif(1,-0.01,0.01)
initdist.shift <- 0

inits <- list(beta.init = b - runif(1,0,b/3),
              mu.init = m - runif(1,0,m/3),
              alpha.init = 0, 
              probs.init = samp_prob + runif(1,-0.05,0.05),
              initdist.init = initdist + c(initdist.shift,-initdist.shift,0))

# priors
priors <- list(beta.prior = c(5/popsize + runif(1,-0.1/popsize,0.1/popsize), 0.001),
               mu.prior = c(0.001, 0.001),
               alpha.prior = NULL,
               p.prior = c(1,1),
               init.prior = c(1,1))

tmax <- max(SIRres$results[,1])


# vectors for parameters
Beta <- vector(length = niter); Beta[1] <- b + rnorm(1, 0, 1e-5)
Mu <- vector(length = niter); Mu[1] <- m + rnorm(1, 0, 1e-5)
Alpha <- vector(length = niter); Alpha[1] <- 0
probs <- vector(length = niter); probs[1] <- samp_prob
# Beta <- vector(length=niter); Beta[1] <- inits$beta.init
# Mu <- vector(length = niter); Mu[1] <- inits$mu.init
# Alpha <- vector(length = niter); Alpha[1] <- inits$alpha.init 
# probs <- vector(length = niter); probs[1] <- inits$probs.init
p_initinfec <- vector(length = niter)
suff.stats <- matrix(nrow = niter, ncol = 4); colnames(suff.stats) <- c("num_infec", "num_recov", "beta_suff.stat", "mu_suff.stat")
R_0 <- vector(length = niter)

accepts <- vector(length = niter)


# vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
beta.prior <- priors$beta.prior
mu.prior <- priors$mu.prior
# alpha.prior <- c(6, 12000)
p.prior <- priors$p.prior
init.prior <- priors$init.prior


# log-likelihood vector
loglik <- vector(length=niter)

# list to store trajectories
trajectories <- vector(mode = "list", length = 1 + niter/trajecs_every)

# observation matrix
W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Observed, augmented = 0))

# set initial distribution
# init.prob <- rbeta(1, shape1 = init.prior[1] + W.cur[1,2]/probs[1], shape2 = init.prior[2] + popsize - W.cur[1,2]/probs[1])
# initdist <- c(1 - init.prob, init.prob, 0)
initdist <- initdist
p_initinfec[1] <- initdist[2]

# individual trajectories
X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], samp_prob = probs[1], initdist = initdist, censusInterval = censusInterval, tmax = tmax, popsize = popsize)

# count matrix
Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)

trajectories[[1]] <- Xcount.cur

# update observation matrix
W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)

# calculate likelihood for initial trajectory
loglik[1] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = Beta[1], m = Mu[1], a = Alpha[1], p = probs[1], initdist = initdist, popsize = popsize)

# calculate sufficient statistics for initialization
suff.stats[1,] <- c(numinfec = sum(diff(Xcount.cur[,2], lag = 1)>0), 
                    numrecov =  sum(diff(Xcount.cur[,2], lag = 1)<=0), 
                    beta_suffstat = sum(Xcount.cur[1:(nrow(Xcount.cur) - 1),2]*Xcount.cur[1:(nrow(Xcount.cur)-1),3]*diff(Xcount.cur[,1], lag = 1)),
                    mu_suffstat = sum(Xcount.cur[1:(nrow(Xcount.cur) - 1),2]*diff(Xcount.cur[,1], lag = 1)))
R_0[1] <- Beta[1] * popsize / Mu[1]

writeLines(c(""), paste(paste("augSIR_log",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"),".txt", sep=""))

start.time <- Sys.time()
# start sampler
for(k in 2:niter){

    if(k %% 100 == 0){
        sink(paste(paste("log",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"),".txt", sep=""), append=TRUE)  
        cat(paste("Starting iteration",k,"\n"))  
        sink()
    }
    
    subjects <- sample(unique(X.cur[,2]), floor(popsize*resample_prop), replace=TRUE)

    pathirm.cur <- build_irm(Xcount = Xcount.cur, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], pop = FALSE)
    patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
    
    pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], initdist = initdist, popsize = popsize)
    
    accepts.k <- 0
    
    # redraw individuals
    for(j in 1:length(subjects)){
        # get current path
        path.cur <- getpath(X.cur, subjects[j])
        
        # get .other matrices
        Xother <- X.cur[X.cur[,2]!=subjects[j],]
        Xcount.other <- build_countmat(Xother, popsize - 1)
        W.other <- get_W_other(W.cur = W.cur, path = path.cur)
        
        # draw new path
        path.new <- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = samp_prob, initdist = initdist, tmax = tmax)
        
        # generate .new objects
        X.new <- updateX(X = X.cur, path = path.new, j = subjects[j])
        Xcount.new <- update_Xcount(Xcount.other = Xcount.other, path = path.new)
        W.new <- updateW(W = W.other, path = path.new)
        
        # check if update is needed to path.irm
        if(max(Xcount.new[,2]) == pathirm.cur[4,4,dim(pathirm.cur)[3]]){
            
            new.numinf <- pathirm.cur[4,4,dim(pathirm.cur)[3]]+1
            
            pathirm.cur <- update_irm(irm = pathirm.cur, new.numinf = new.numinf, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize)
            patheigen.cur <- update_eigen(patheigen = patheigen.cur, pathirm = pathirm.cur)
            
        }
        
        # Metropolis-Hastings 
        # population trajectory likelihoods 
        pop_prob.new <- pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
        
        # path likelihoods
        path_prob.new <- path_prob(path = path.new, Xcount.other = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
        path_prob.cur <- path_prob(path = path.cur, Xcount.other = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
        
        # compute log acceptance ratio
        a.prob <- ifelse(!is.na(pop_prob.new),accept_prob(pop_prob.new = pop_prob.new, pop_prob.cur = pop_prob.cur, path_prob.cur = path_prob.cur, path_prob.new = path_prob.new), -Inf)
        
        # accept or reject
        if(a.prob >=0 || exp(a.prob) > runif(1)) {
            X.cur <- X.new
            Xcount.cur <- Xcount.new
            W.cur <- W.new
            pop_prob.cur <- pop_prob.new
            accepts.k <- accepts.k + 1
        } 
    }
    
    # save proportion accepted
    accepts[k-1] <- mean(accepts.k)
    
    # draw new parameters
    # probs[k] <- update_prob(W = W.cur, p.prior = p.prior)
    probs[k] <- probs[k-1]
    
    # new rate parameters 
    params.update <- update_rates2(Xcount = Xcount.cur, beta.prior = beta.prior, mu.prior = mu.prior, init.prior = init.prior, popsize = popsize)
    
    params.new <- params.update[[1]]
    suff.stats[k, ] <- params.update[[2]]
    
    # Beta[k] <- params.new[1]
    Beta[k] <- Beta[k-1]
    
    # Mu[k] <- params.new[2]
    Mu[k] <- Mu[k-1]
    
    # Alpha[k] <- params.new[3]
    Alpha[k] <- Alpha[k-1]
    
    # p_initinfec[k] <- params.new[5]
    p_initinfec[k] <- p_initinfec[k-1]
    # initdist <- params.new[4:6]

    loglik[k] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = Beta[k], m = Mu[k], a = Alpha[k], p = probs[k], initdist = initdist, popsize = popsize)  
    
    R_0[k] <- Beta[k] * popsize / Mu[k]
    
    if(k%%trajecs_every == 0) {
        trajectories[[k/trajecs_every]] <- Xcount.cur
    }
    
    W.cur[,2] <- rbinom(nrow(W.cur), W.cur[,3], probs[1])
    
    gillespie_trajecs[[k]] <- build_countmat(SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval = censusInterval, sampprob = samp_prob, returnX = TRUE, trim = FALSE)$trajectory, popsize)
    
    gillespie_extended[[k]] <- build_countmat(SIRsim_extended(popsize = popsize, initdist = initdist, b = b, m = m, obstimes = obstimes, samp_prob = samp_prob, trim = FALSE)$trajectory, popsize)
}
end.time <- Sys.time()

total.time <- difftime(end.time, start.time, units = "mins")

results <- list(trajectories, gillespie_trajecs, gillespie_extended)

assign(paste("augSIR_",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"),results)

save(list = paste("augSIR_",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"), file = paste(paste("augSIR_",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"),".Rdata", sep=""))


# plot stuff
pdf(file = paste(paste("convcomp_plots_augSIR",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"),".pdf", sep=""))

    times <- seq(0, tmax, by = 0.05)
    augSIR <- matrix(0, nrow = length(times), ncol = niter)
    gillespie <- matrix(0, nrow = length(times), ncol = niter)
    extended <- matrix(0, nrow = length(times), ncol = niter)
    
    for(s in 1:niter){
        
        augSIR_traj <- results[[1]][[s]]
        gillespie_traj <- results[[2]][[s]]
        gillespie_ext <- results[[3]][[s]]
        
        for(t in 1:length(times)){
            augSIR[t,s] <- augSIR_traj[sum(augSIR_traj[,"time"] <= times[t]), "numsick"]
            gillespie[t,s] <- gillespie_traj[sum(gillespie_traj[,"time"] <= times[t]), "numsick"]
            extended[t,s] <- gillespie_ext[sum(gillespie_ext[,"time"] <= times[t]), "numsick"]
        }
    }
    
    results_comp <- data.frame(method = rep(c("augSIR", "Gillespie Counts", "Gillespie Extended"), each = length(times)), time = rep(times, 3), count = 0)
    results_comp[results_comp$method == "augSIR","count"] <- rowMeans(augSIR)
    results_comp[results_comp$method == "Gillespie Counts","count"] <- rowMeans(gillespie)
    results_comp[results_comp$method == "Gillespie Extended","count"] <- rowMeans(extended)
    
    ggplot(results_comp, aes(x = time, y = count, colour = method)) + geom_line() + theme_bw() + labs(title = "augSIR vs. Gillespie")

dev.off()
