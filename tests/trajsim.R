# evaluating how well the method reconstructs the true trajectory

# set simulation parameters
niter <- 50
popsize = 50; tmax = 20
b <- 0.05 + runif(1, -0.0001, 0.0001)
m <- 1 + runif(1, -0.0001, 0.0001)
samp_prob <- 0.5 
initdist <- c(0.98, 0.02, 0)
censusInterval <- 0.05
obstimes <- seq(0, tmax, by = censusInterval)

# generate dataset
SIRres<-SIRsim(popsize = popsize, initdist = initdist, b = b, mu = m, a=0, tmax = tmax, censusInterval=censusInterval, sampprob = samp_prob, returnX = TRUE)

# plot dataset
# get data 
# dat <- SIRres$results[SIRres$results[,1]%%0.1 == 0, ]
dat <- SIRres$results
dat.m <- melt(dat,id.vars="time")

ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_line() + theme_bw()

# Settings
# sim.settings <- list(popsize = popsize,
#                      tmax = max(SIRres$results[,1]),
#                      niter = niter,
#                      amplify = 10,
#                      initdist = initdist)
# 
inits <- list(beta.init = b,
              mu.init = m,
              alpha.init = 0, 
              probs.init = samp_prob)

# priors
priors <- list(beta.prior = c(0.0005, 0.001),
               mu.prior = c(0.001, 0.001),
               alpha.prior = NULL,
               p.prior = c(1,1))

tmax <- max(SIRres$results[,1])

# vectors for parameters
Beta <- vector(length=niter); Beta[1] <- inits$beta.init
Mu <- vector(length = niter); Mu[1] <- inits$mu.init
Alpha <- vector(length = niter); Alpha[1] <- inits$alpha.init 
probs <- vector(length = niter); probs[1] <- inits$probs.init
accepts <- vector(length = niter)


# vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
beta.prior <- priors$beta.prior
mu.prior <- priors$mu.prior
# alpha.prior <- c(6, 12000)
p.prior <- priors$p.prior


# log-likelihood vector
loglik <- vector(length=niter)

# list to store trajectories
trajectories <- vector(mode = "list", length = niter)

# observation matrix
W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Observed, augmented = 0))

# individual trajectories
X.cur <- SIRres$trajectory

# count matrix
Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)

# update observation matrix
W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)

# initialize population likelihood
pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = Beta[1], m = Mu[1], a = Alpha[1], initdist = initdist, popsize = popsize)

# calculate likelihood for initial trajectory
loglik[1] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = Beta[1], m = Mu[1], a = Alpha[1], p = probs[1], initdist = initdist, popsize = 200)


# start sampler
for(k in 2:niter){
    
    print(k)
    subjects <- sample(unique(X.cur[,2]), popsize, replace=TRUE)
    # subjects <- 1:popsize
    
    pathirm.cur <- build_irm(Xcount = Xcount.cur, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = FALSE)
    patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
    
    accepts.k <- 0
    
    # redraw individuals
    for(j in 1:length(subjects)){
        # get current path
        path.cur <- getpath(X.cur, subjects[j])
        
        # get .other matrices
        Xother <- X.cur[X.cur[,2]!=subjects[j],]
        Xcount.other <- build_countmat(Xother, popsize - 1)
        W.other <-get_W_other(W.cur = W.cur, path = path.cur)
        
        # draw new path
        path.new<- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = samp_prob, initdist = initdist, tmax = tmax)
        
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
        
        # Metropolis-Hastings (should always accept, but include check to make sure)
        # population trajectory likelihoods 
        # pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
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
    accepts[k-1] <- mean(accepts.k)
    
    # draw new parameters
        probs[k] <- update_prob(W = W.cur, p.prior = p.prior)
    # probs[k] <- probs[k-1]
    # new rate parameters 
        params.new <- update_rates(Xcount = Xcount.cur, beta.prior = beta.prior, mu.prior = mu.prior, popsize = popsize)
    # params.new <- c(Beta[k-1], Mu[k-1], Alpha[k-1])
    
    Beta[k] <- params.new[1]
    
    Mu[k] <- params.new[2]
    
    Alpha[k] <- params.new[3]
    
    loglik[k] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = Beta[k], m = Mu[k], a = Alpha[k], p = probs[k], initdist = initdist, popsize = popsize)  
    
    trajectories[[k]] <- Xcount.cur
    
}

results <- list(Beta = Beta, Mu = Mu, probs=probs, loglik = loglik, accepts = accepts, trajectories = trajectories) 



censusInterval <- censusInterval; p <- samp_prob
trajectories2 <- list(); observations2 <- list(); likelihoods <- list()

for(k in 2:(length(trajectories))){
    if ((k%%1)==0){
        traj <- results$trajectories[[k]]
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
        
        obs$sampled <- rbinom(dim(obs)[1], obs$truth, prob = p)
        
        observations2[[k]] <- obs
    }
    
}

trajecs <- do.call(rbind,trajectories2)
samples <- do.call(rbind,observations2)
likelihoods <- do.call(rbind,likelihoods)

truetrajec <- data.frame(time = unique(SIRres$trajectory[,1]), 
                         infected = c(sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]), sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]) + cumsum(SIRres$trajectory[SIRres$trajectory[,1]!=0,3])),
                         simnum = 0)


ggplot(data = trajecs, aes(x = time, y = infected, group = simnum)) + geom_path(alpha = 0.3) + geom_path(data = truetrajec, aes(x = time, y = infected),colour="red", size = 1) + geom_boxplot(data = samples, aes(x=time, y=sampled, group= time),outlier.shape=NA) + geom_point(data=data.frame(dat,simnum=0), aes(x=time,y=Observed),size=4,colour="blue") + theme_bw() + scale_x_continuous(limits = c(0, tmax))
