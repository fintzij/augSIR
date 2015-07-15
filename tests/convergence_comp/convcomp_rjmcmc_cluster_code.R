source("augSIR.R")
source("auxilliary_functions.R")
source("SIRsimulation.R")
source("rjmcmc_functions.R")
source("matrix_build_update.R")
source("metropolis_hastings_functions.R")
source("path_sampling_functions.R")

args <- commandArgs(TRUE)
print(args)

popsize <- as.numeric(args[1])
censusInterval <- as.numeric(args[2])
samp_prob <- as.numeric(args[3])
resample_prop <- as.numeric(args[4])
initialization_num <- as.numeric(args[5])

tmax = 10
b <- 5/popsize
m <- 1
samp_prob <- 0.5 
initdist <- c(0.95, 0.05, 0)
insert.prob = 1/3; remove.prob = 1/3; shift.prob = 1/3
shift.int <- 0.1
niter <- 25000

obstimes <- seq(0, tmax, by = censusInterval)

set.seed(183427)
SIRres<-SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval=censusInterval, sampprob = samp_prob, returnX = TRUE)

if(initialization_num == 1){
    set.seed(561916)
} else if(initialization_num == 2){
    set.seed(683629)
} else if(initialization_num == 3){
    set.seed(724318)
}

initdist.shift <- runif(1,-0.01,0.01)

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
Beta <- vector(length=niter); Beta[1] <- inits$beta.init
Mu <- vector(length = niter); Mu[1] <- inits$mu.init
Alpha <- vector(length = niter); Alpha[1] <- inits$alpha.init 
probs <- vector(length = niter); probs[1] <- inits$probs.init
p_initinfec <- vector(length = niter)
suff.stats <- matrix(nrow = niter, ncol = 4); colnames(suff.stats) <- c("num_infec", "num_recov", "beta_suff.stat", "mu_suff.stat")

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
trajectories <- vector(mode = "list", length = 100)

# observation matrix
W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Observed, augmented = 0))

# individual trajectories
X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = 0, p = probs[1], amplify = 2, tmax = tmax, popsize = popsize)

# count matrix
Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)

# update observation matrix
W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)

# initialize population likelihood
pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = Beta[1], m = Mu[1], a = Alpha[1], initdist = initdist, popsize = popsize)

# calculate likelihood for initial trajectory
loglik[1] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = Beta[1], m = Mu[1], a = Alpha[1], p = probs[1], initdist = initdist, popsize = 200)

# set initial distribution
initdist <- c(1-W.cur[1,2]/probs[1]/popsize, W.cur[1,2]/probs[1]/popsize, 0)
p_initinfec[1] <- initdist[2]

# calculate sufficient statistics for initialization
suff.stats[1,] <- c(numinfec = sum(diff(Xcount.cur[,2], lag = 1)>0), 
                    numrecov =  sum(diff(Xcount.cur[,2], lag = 1)<=0), 
                    beta_suffstat = sum(Xcount.cur[1:(nrow(Xcount.cur) - 1),2]*Xcount.cur[1:(nrow(Xcount.cur)-1),3]*diff(Xcount.cur[,1], lag = 1)),
                    mu_suffstat = sum(Xcount.cur[1:(nrow(Xcount.cur) - 1),2]*diff(Xcount.cur[,1], lag = 1)))

writeLines(c(""), paste(paste("log",popsize,censusInterval,samp_prob,resample_prop,initialization_num,sep="_"),".txt", sep=""))

start.time <- Sys.time()
# start sampler
for(k in 2:niter){
    
    if(k %% 2500 == 0){
        sink(paste(paste("log",popsize,censusInterval,samp_prob,resample_prop,initialization_num,sep="_"),".txt", sep=""), append=TRUE)  
        cat(paste("Starting iteration",k,"\n"))  
        sink()
    }
    
    subjects <- sample(unique(X.cur[,2]), floor(popsize*resample_prop), replace=TRUE)
    # subjects <- 1:popsize
    
    pathirm.cur <- build_irm(Xcount = Xcount.cur, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = FALSE)
    patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
    
    accepts.k <- 0
    
    # redraw individuals
    for(j in 1:length(subjects)){
        
        path.cur <- X.cur[X.cur[,2] == subjects[j], ]
        path.new <- rjmcmc_draw(path.cur = path.cur, Xcount.cur, j = subjects[j], initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, tmax = max(W.cur[,1]), b = b, m = m, p = samp_prob)
        
        # update bookkeeping objects
        X.new <- X.cur; X.new[X.new[,2] == subjects[j], ] <- path.new
        Xcount.new <- build_countmat(X = X.new, popsize = popsize)
        W.new <- updateW(W = W.cur, Xcount = Xcount.new)
        
        # calculate acceptance ratio
        rjmcmc.ratio <- rjmcmc_ratio(W.cur = W.cur, W.new = W.new, X.cur = X.cur, X.new = X.new, Xcount.cur = Xcount.cur, Xcount.new = Xcount.new, path.cur = path.cur, path.new = path.new, initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, b = b, m = m, samp_prob = samp_prob, tmax = tmax, popsize = popsize)
        
        # decide whether to accept. if accept, update current trajectories
        if(rjmcmc.ratio > log(runif(1))){
            
            X.cur <- X.new
            Xcount.cur <- Xcount.new
            W.cur <- W.new
            accepts.k <- accepts.k + 1
            
        }
        
    }
    
    # save proportion accepted
    accepts[k-1] <- mean(accepts.k)
    
    # draw new parameters
    probs[k] <- update_prob(W = W.cur, p.prior = p.prior)
    
    # new rate parameters 
    params.update <- update_rates2(Xcount = Xcount.cur, beta.prior = beta.prior, mu.prior = mu.prior, init.prior = init.prior, popsize = popsize)
    
    params.new <- params.update[[1]]
    suff.stats[k, ] <- params.update[[2]]
    
    Beta[k] <- params.new[1]
    
    Mu[k] <- params.new[2]
    
    Alpha[k] <- params.new[3]
    
    p_initinfec[k] <- params.new[5]
    initdist <- params.new[4:6]
    
    loglik[k] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = Beta[k], m = Mu[k], a = Alpha[k], p = probs[k], initdist = initdist, popsize = popsize)  
    
    if(k%%250 == 0){
        trajectories[[k/50]] <- Xcount.cur
    }
}
end.time <- Sys.time()

total.time <- difftime(end.time, start.time)

results <- list(total.time, quantities = cbind(loglik, Beta, Mu, probs, p_initinfec, suff.stats, accepts), trajectories)

assign(paste("rjmcmc",popsize,censusInterval,samp_prob,resample_prop,initialization_num,sep="_"),results)

save(list = paste("rjmcmc",popsize,censusInterval,samp_prob,resample_prop,initialization_num,sep="_"), file = paste(paste("rjmcmc",popsize,censusInterval,samp_prob,resample_prop,initialization_num,sep="_"),".Rdata", sep=""))
