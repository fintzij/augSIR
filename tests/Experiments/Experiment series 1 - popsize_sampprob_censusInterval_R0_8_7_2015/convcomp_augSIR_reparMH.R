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
initdist_true <- initdist
# insert.prob = 1/3; remove.prob = 1/3; shift.prob = 1/3
# shift.int <- 0.1
niter <- 10000
trajecs_every <- 50

obstimes <- seq(0, tmax, by = censusInterval)

set.seed(183427)
SIRres <- SIRsim_extended(popsize = popsize, initdist = initdist, b = b, m = m, obstimes = obstimes, samp_prob = samp_prob)
if(sum(SIRres$results$Observed) < 5){
    while(sum(SIRres$results$Observed) < 5){
        SIRres <- SIRsim_extended(popsize = popsize, initdist = initdist, b = b, m = m, obstimes = obstimes, samp_prob = samp_prob)
    }
}

tmax <- max(SIRres$results[,1])

set.seed(initialization_num + 9235)

# initdist.shift <- runif(1,-0.01,0.01)
initdist.shift <- 0

inits <- list(beta.init = b - runif(1,-b/6,b/6),
              mu.init = m - runif(1,-m/6,m/6),
              alpha.init = 0, 
              probs.init = samp_prob + runif(1,-0.05,0.05),
              initdist.init = initdist + c(initdist.shift,-initdist.shift,0))

# functions to evaluate prior densities
d_priors <- list(R0_prior = function(r) dnorm(r, 0, 1.8, log = TRUE),
               mu_prior = function(u) dnorm(u, log(m + runif(1, -0.01, 0.01)), 3, log = TRUE),
               p_prior = function(b) dnorm(b, logit(samp_prob + runif(1, -0.01, 0.01)), 3, log = TRUE),
               p_init_prior = function(i) dnorm(i, logit(initdist_true[2] + runif(1, -0.005, 0.005)), 3, log = TRUE))

# lists of functions to move between model and estimation scales for parameters
toEst <- list(function(params) log(params[1]*popsize/params[2]),
              function(params) log(params[2]),
              function(params) logit(params[3]),
              function(params) logit(params[4]))

fromEst <- list(function(params_est) exp(params_est[2])/popsize * exp(params_est[1]),
                function(params_est) exp(params_est[2]),
                function(params_est) expit(params_est[3]),
                function(params_est) expit(params_est[4]))

# transition kernel for MH proposals
kernel <- function(params_est, which_par){
    
    inds <- rep(0, length(params_est)); inds[which_par] <- 1
    
    params_new <- params_est + rnorm(length(params_est), 0, c(0.2, 0.2, 0.2, 0.2)) * inds
    
    return(params_new)
    
}

# vectors for parameters
Beta <- vector(length=niter); Beta[1] <- inits$beta.init
Mu <- vector(length = niter); Mu[1] <- inits$mu.init
Alpha <- vector(length = niter); Alpha[1] <- inits$alpha.init 
probs <- vector(length = niter); probs[1] <- inits$probs.init
p_initinfec <- vector(length = niter)
R_0 <- vector(length = niter); R_0[1] <- Beta[1] * popsize / Mu[1]

# suff.stats <- matrix(nrow = niter, ncol = 4); colnames(suff.stats) <- c("num_infec", "num_recov", "beta_suff.stat", "mu_suff.stat")

# vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
# beta.prior <- priors$beta.prior
# mu.prior <- priors$mu.prior
# # alpha.prior <- c(6, 12000)
# p.prior <- priors$p.prior
# init.prior <- priors$init.prior


# log-likelihood vector
loglik <- vector(length=niter)

# list to store trajectories
trajectories <- vector(mode = "list", length = 1 + niter/trajecs_every)

# observation matrix
W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Observed, augmented = 0))

# set initial distribution
init.prob <- rbeta(1, shape1 = 1 + W.cur[1,2]/probs[1], shape2 = 1 + popsize - W.cur[1,2]/probs[1])
initdist <- c(1 - init.prob, init.prob, 0)
p_initinfec[1] <- initdist[2]

# set vector of current parameters on model and estimation scale
params.cur <- c(Beta[1], Mu[1], probs[1], p_initinfec[1])
params_est.cur <- sapply(toEst, function(f) f(params.cur)) 

# set vector of prior probabilities
prior_probs.cur <- evaluate_priors(params_est.cur, d_priors)


# initialize vector of estimation scale parameters
l_R0 <- vector(length = niter); l_R0[1] <- params_est.cur[1]
l_Mu <- vector(length = niter); l_Mu[1] <- params_est.cur[2]
l_probs <- vector(length = niter); l_probs[1] <- params_est.cur[3]
l_p_initinfec <- vector(length = niter); l_p_initinfec[1] <- params_est.cur[4]

# vectors to track acceptances
traj_accepts <- vector(length = niter)
param_accepts <- matrix(0, nrow = niter, ncol = length(params.cur))
colnames(param_accepts) <- c("R0_accepts", "Mu_accepts", "Samp_prob_accepts", "p_init_accepts")

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
# suff.stats[1,] <- c(numinfec = sum(diff(Xcount.cur[,2], lag = 1)>0), 
#                     numrecov =  sum(diff(Xcount.cur[,2], lag = 1)<=0), 
#                     beta_suffstat = sum(Xcount.cur[1:(nrow(Xcount.cur) - 1),2]*Xcount.cur[1:(nrow(Xcount.cur)-1),3]*diff(Xcount.cur[,1], lag = 1)),
#                     mu_suffstat = sum(Xcount.cur[1:(nrow(Xcount.cur) - 1),2]*diff(Xcount.cur[,1], lag = 1)))

writeLines(c(""), paste(paste("augSIR_log",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"),".txt", sep=""))

start.time <- Sys.time()
# start sampler
for(k in 2:niter){

    if(k %% 500 == 0){
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
        a.prob <- accept_prob(pop_prob.new = pop_prob.new, pop_prob.cur = pop_prob.cur, path_prob.cur = path_prob.cur, path_prob.new = path_prob.new)
        
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
    traj_accepts[k-1] <- mean(accepts.k)
    
    # get current log-likelihood
    log_lik.cur <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], p = probs[k-1], initdist = initdist, popsize = popsize)
    
    # draw new parameters
    for(p in 1:length(params.cur)){
        
        # propose new parameter
        param_prop <- propose_params(p, params.cur, kernel, toEst, fromEst)
        
        # update vector prior probabilities
        prior_probs.new <- update_prior_prob(p, prior_probs = prior_probs.cur, params_est = param_prop[[1]], d_priors = d_priors)
        
        # calculate new log_likelihood
        log_lik.new <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = param_prop[[1]][1], m = param_prop[[1]][2], a = 0, p = param_prop[[1]][3], initdist = c(1-param_prop[[1]][4], param_prop[[1]][4], 0), popsize = popsize)
        
        a.prob_MH <- accept_prob_param_MH(log_lik.new = log_lik.new, log_lik.cur = log_lik.cur, prior_probs.new = prior_probs.new, prior_probs.cur = prior_probs.cur)
        
        if(a.prob_MH >=0 || exp(a.prob_MH) > runif(1)) {
            params.cur <- param_prop[[1]] # update vector of parameters on model scale
            params_est.cur <- param_prop[[2]] # update parameters on estimation scale
            log_lik.cur <- log_lik.new # retain new log likelihood
            prior_probs.cur <- prior_probs.new # retain new vector of prior probs 
            param_accepts[k-1, p] <- 1
        } 
        
    }
    
    # record parameter values
    Beta[k] <- params.cur[1]; l_R0[k] <- params_est.cur[1] 
    Mu[k] <- params.cur[2]; l_Mu[k] <- params_est.cur[2]
    probs[k] <- params.cur[3]; l_probs[k] <- params_est.cur[3]
    p_initinfec[k] <- params.cur[4]; l_p_initinfec[k] <- params_est.cur[4]
    
    initdist <- c(1 - params.cur[3], params.cur[3], 0)

    loglik[k] <- log_lik.cur 
    
    if(k%%trajecs_every == 0) {
        trajectories[[1 + k/trajecs_every]] <- Xcount.cur
    }
    
    if(k == 200 && all(traj_accepts == 0)){
        sink(paste(paste("BAD_INITIALIZATION",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"),".txt", sep=""), append=TRUE)  
        cat("SAMPLER HAS NOT ACCEPTED ANY TRAJECTORIES!")  
        sink()
    }
}
end.time <- Sys.time()

total.time <- difftime(end.time, start.time, units = "mins")

results <- list(total.time, quantities = cbind(loglik, Beta, R0 = (Beta * popsize / Mu), Mu, probs, p_initinfec, traj_accepts,param_accepts), trajectories)

assign(paste("augSIR_",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"),results)

save(list = paste("augSIR_",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"), file = paste(paste("augSIR_",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"),".Rdata", sep=""))


# plot stuff
pdf(file = paste(paste("convcomp_plots_augSIR",popsize,censusInterval,samp_prob,R0,initialization_num,sep="_"),".pdf", sep=""))
par(mfrow = c(5,2))
ts.plot(results[[2]][,"Beta"]); plot(density(results[[2]][,"Beta"]), main = "Beta"); abline(v = b, col = "red")
ts.plot(results[[2]][,"Mu"]); plot(density(results[[2]][,"Mu"]), main = "Mu"); abline(v = m, col = "red")
ts.plot(results[[2]][,"R_0"]); plot(density(results[[2]][,"R_0"]), main = "R0"); abline(v = 0.05, col = "red")
ts.plot(results[[2]][,"probs"]); plot(density(results[[2]][,"probs"]), main = "Binomial probability"); abline(v = samp_prob, col = "red")
ts.plot(results[[2]][,"p_initinfec"]); plot(density(results[[2]][,"p_initinfec"]), main = "Prob. of infection at time 0"); abline(v = 0.05, col = "red")

par(mfrow = c(1,1))
pairs(results[[2]][seq(1,niter,by=5),1:5])
ts.plot(loglik, main = "log-likelihood", xlab = "iteration")
    
    truetrajec <- data.frame(time = unique(SIRres$trajectory[,1]), 
                             infected = c(sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]), sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]) + cumsum(SIRres$trajectory[SIRres$trajectory[,1]!=0,3])),
                             simnum = 0)
    
    trajectories <- results[[3]]
    trajectories2 <- list(); observations2 <- list()
    for(k in 1:length(results[[3]])){
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
    
    ggplot(data = trajecs, aes(x = time, y = infected, group = simnum)) + geom_path(alpha = 0.3) + geom_path(data = truetrajec, aes(x = time, y = infected),colour="red", size = 1)  + geom_point(data=data.frame(SIRres$results,simnum=0), aes(x=time,y=Observed),size=4,colour="blue") + theme_bw() + scale_x_continuous(limits = c(0, tmax))

dev.off()
