
# simulate data

SIRres<-SIRsim(popsize = 3, initdist = c(0.7, 0.3, 0), b = 0.5, mu=1, a=0, tmax = 3, censusInterval=0.05, sampprob = 0.25, returnX = TRUE, trim = FALSE)


# get data 
dat <- SIRres$results
dat.m <- melt(dat,id.vars="time")
# 
ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + ylim(c(0,3)) + xlim(0,3) + theme_bw() 

colnames(dat) <- c("time", "Observed", "Truth")


# Settings
sim.settings <- list(popsize = 3,
                     tmax = max(dat[,1]),
                     niter = 10000,
                     amplify = 10,
                     initdist = c(0.7, 0.3, 0))

inits <- list(beta.init = 0.5 + runif(1,-0.0005, 0.0005),
              mu.init = 1 + runif(1, -0.005, 0.005),
              alpha.init = 0, 
              probs.init = 0.5 + runif(1,-0.005, 0.005))

priors <- list(beta.prior = c(1e-4, 2e-2),
               mu.prior = c(1e-4, 1e-2),
               alpha.prior = NULL,
               p.prior = c(1,1))


popsize <- sim.settings$popsize # size of the population; 
tmax <- sim.settings$tmax # maximum time of observation
amplify <- sim.settings$amplify # amplification parameter for initialization
niter <- sim.settings$niter # number of iterations in the sampler
initdist <- sim.settings$initdist # initial distribution for individual infection status

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
trajectories <- list()

# observation matrix
W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$Observed, augmented = 0))

# matrix with event times, subject id and event codes. 
# Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
X.cur <- SIRres$trajectory
# X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], p=probs[1], amplify = amplify, tmax=tmax, popsize = popsize)

# update Xcount.cur and W.cur
Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)# update observation matrix
W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)

if(!checkpossible(X=X.cur, W=W.cur)) {
    while(!checkpossible(X=X.cur,W=W.cur)){
        X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], p=probs[1], amplify = amplify, tmax=tmax, popsize = popsize)
        Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
        W.cur <- updateW(W = W.cur, X = X.cur)
    }
}

pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = Beta[1], m = Mu[1], a = Alpha[1], initdist = initdist, popsize = popsize)

loglik[1] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = Beta[1], m = Mu[1], a = Alpha[1], p = probs[1], initdist = initdist, popsize = popsize)

keep.going <- TRUE

# M-H sampler
for(k in 2:niter){
    if(keep.going == FALSE) break
    # Update trajectories
    if(k%%100 == 0) print(k)
    subjects <- sample(unique(X.cur[,2]), popsize, replace=TRUE)
    
    pathirm.cur <- build_irm(Xcount = Xcount.cur, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = FALSE)
    patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
    
    accepts.k <- 0; keep.going <- TRUE
    
    for(j in 1:length(subjects)){
        
        Xother <- X.cur[X.cur[,2]!=subjects[j],]
        
        path.cur <- getpath(X.cur, subjects[j])
        
        if(!all(path.cur == Xcount.cur[,1])){
        
        Xcount.other <- get_Xcount_other(Xcount = Xcount.cur, path = path.cur)
        
        if(nrow(Xcount.other)!=1){
        
        if(check_subj_pop_mats(Xother, Xcount.other) == FALSE | Xcount.other[1,3] == popsize){
            
            keep.going <- FALSE
            print("STOP AND CHECK .other!")
            break
        }
        
        W.other <-get_W_other(W.cur = W.cur, path = path.cur)
        
        path.new<- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = probs[k-1], initdist = initdist, tmax = tmax)
        
        X.new <- updateX(X = X.cur, path = path.new, j = subjects[j])
        Xcount.new <- update_Xcount(Xcount.other = Xcount.other, path = path.new)
        
        if(check_subj_pop_mats(X.new, Xcount.new) == FALSE | any(Xcount.new[,c(2:3)] > popsize) | any(rowSums(Xcount.new[,c(2:3)]) > popsize)){
            
            keep.going <- FALSE
            print("STOP AND CHECK .new!")
            break
            
        }
        
        W.new <- updateW(W = W.other, path = path.new)
        
        if(max(Xcount.new[,2]) == pathirm.cur[4,4,dim(pathirm.cur)[3]]){
            
            new.numinf <- pathirm.cur[4,4,dim(pathirm.cur)[3]]+1
            
            pathirm.cur <- update_irm(irm = pathirm.cur, new.numinf = new.numinf, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize)
            patheigen.cur <- update_eigen(patheigen = patheigen.cur, pathirm = pathirm.cur)
            
        } 
        
        
        pop_prob.new <- pop_prob(Xcount = Xcount.new, tmax = tmax, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], initdist = initdist, popsize = popsize)
        
        path_prob.new <- path_prob(path = path.new, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
        path_prob.cur <- path_prob(path = path.cur, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
        
        a.prob <- accept_prob(pop_prob.new = pop_prob.new, pop_prob.cur = pop_prob.cur, path_prob.cur = path_prob.cur, path_prob.new = path_prob.new)
        
        if(min(a.prob, 0) > log(runif(1))) {
            X.cur <- X.new
            Xcount.cur <- Xcount.new
            W.cur <- W.new
            pop_prob.cur <- pop_prob.new
            accepts.k <- accepts.k + 1
        }
        }
        }
        
    }
    
    # save proportion accepted
    accepts[k-1] <- mean(accepts.k)
    
    # draw new parameters
#         probs[k] <- update_prob(W = W.cur, p.prior = p.prior)
    probs[k] <- probs[k-1]
    # new rate parameters 
#         params.new <- update_rates(Xcount = Xcount.cur, beta.prior = beta.prior, mu.prior = mu.prior, popsize = popsize)
    params.new <- c(Beta[k-1], Mu[k-1], Alpha[k-1])
    
    Beta[k] <- params.new[1]
    
    Mu[k] <- params.new[2]
    
    Alpha[k] <- params.new[3]
    
#     loglik[k] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = Beta[k], m = Mu[k], a = Alpha[k], p = probs[k], initdist = initdist, popsize = popsize)  
    
    trajectories[[k]] <- Xcount.cur
    
}

results2 <- list(Beta = Beta, Mu = Mu, probs=probs, loglik = loglik, accepts = accepts, trajectories = trajectories) 


# Results for the case with error, p=0.25  -----------------------------------------------------------------

censusInterval <- 0.05; p <- 0.5
trajectories2 <- list(); observations2 <- list(); likelihoods <- list()

# for(k in 1:(length(results2[[6]]))){
for(k in 2:(length(trajectories))){
    if ((k%%50)==0){
        traj <- results2$trajectories[[k]]
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
trajecs$infected <- trajecs$infected + rnorm(nrow(trajecs), sd = 0.03)
samples <- do.call(rbind,observations2)
likelihoods <- do.call(rbind,likelihoods)

truetrajec <- data.frame(build_countmat(SIRres$trajectory, popsize), simnum = 0)

sample_wide <- reshape(samples, timevar = "simnum", idvar = c("time"), direction = "wide", drop = "truth")
sample_wide <- data.frame(time = sample_wide$time[1:(nrow(sample_wide) - 2)], avg = apply(sample_wide[1:(nrow(sample_wide) - 2), 2:ncol(sample_wide)], 1, mean))

# truetrajec <- data.frame(time = gSSA_res[,1], 
#                          infected = gSSA_res[,2],
#                          simnum = 0)

ggplot(data = trajecs, aes(x = time, y = infected, group = simnum)) + geom_step(alpha = 0.1) + geom_step(data = truetrajec, aes(x = time, y = numsick),colour="red", size = 1) + geom_point(data=data.frame(dat,simnum=0), aes(x=time,y=Observed),size=4,colour="blue") + theme_bw()

#     geom_boxplot(data = samples, aes(x=time, y=sampled, group= time),outlier.shape=NA) + 


# Series of simulations ---------------------------------------------------

probs_seq <- seq(0.3,0.9,by=0.15)

pdf(file="smallsamplesime.pdf")

for(s in 1:length(probs_seq)){
    
    p <- probs_seq[s]
    
    for(t in 1:3){
        print(paste("p = ",probs_seq[s],", Iter = ",t, sep=""))
        
        SIRres<-SIRsim(popsize = 3, initdist = c(0.7, 0.3, 0), b = 0.5, mu=1, a=0, tmax = 3, censusInterval=0.05, sampprob = probs_seq[s], returnX = TRUE, trim = FALSE)
        if(sum(SIRres$trajectory[,3]!=0)<=2){
            while(sum(SIRres$trajectory[,3]!=0)<=2){
                SIRres<-SIRsim(popsize = 3, initdist = c(0.7, 0.3, 0), b = 0.5, mu=1, a=0, tmax = 3, censusInterval=0.05, sampprob = probs_seq[s], returnX = TRUE, trim = FALSE) 
            }
        }
        
        # Settings
        # set initial values
        inits <- list(beta.init = 0.5 + runif(1,-0.0005, 0.0005),
                      mu.init = 1 + runif(1, -0.005, 0.005),
                      alpha.init = 0, 
                      probs.init = 0.5 + runif(1,-0.005, 0.005))
        
        
        # vectors for parameters of prior distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
        beta.prior <- c(1e-4, 2e-2)
        mu.prior <- c(1e-4, 1e-2)
        p.prior <- c(1,1)
        
        popsize <- 3 # size of the population; 
        tmax <- 3 # maximum time of observation
        niter <- 50000 # number of iterations in the sampler
        initdist <- c(0.7, 0.3, 0) # initial distribution for individual infection status
        
        # vectors for parameters
        Beta <- vector(length=niter); Beta[1] <- inits$beta.init
        Mu <- vector(length = niter); Mu[1] <- inits$mu.init
        probs <- vector(length = niter); probs[1] <- inits$probs.init
        accepts <- vector(length = niter)
        
        # list to store trajectories
        trajectories <- list()
        
        # observation matrix
        W.cur <- SIRres$results; colnames(W.cur) <- c("time", "sampled", "augmented")
        
        # matrix with event times, subject id and event codes. 
        # Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
        X.cur <- SIRres$trajectory
        # X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], p=probs[1], amplify = amplify, tmax=tmax, popsize = popsize)
        
        # update Xcount.cur and W.cur
        Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)# update observation matrix
        W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
        
        trajectories[[1]] <- Xcount.cur
        
        pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = Beta[1], m = Mu[1], a = 0, initdist = initdist, popsize = popsize)
        
        
        # M-H sampler
        for(k in 2:niter){
            if(k%%1000 == 0) print(k)
            
            # Update trajectories
            subjects <- sample(unique(X.cur[,2]), popsize, replace=TRUE)
            
            pathirm.cur <- build_irm(Xcount = Xcount.cur, b = Beta[k-1], m = Mu[k-1], a = 0, popsize = popsize, pop = FALSE)
            patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
            
            
            for(j in 1:length(subjects)){
                
                Xother <- X.cur[X.cur[,2]!=subjects[j],]
                
                path.cur <- getpath(X.cur, subjects[j])
                
                Xcount.other <- build_countmat(Xother, (popsize - 1))
                
                if(nrow(Xcount.other)!=1){
                    W.other <-get_W_other(W.cur = W.cur, path = path.cur)
                    
                    path.new<- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = probs[k-1], initdist = initdist, tmax = tmax)
                    
                    X.new <- updateX(X = X.cur, path = path.new, j = subjects[j])
                    Xcount.new <- build_countmat(X = X.new, popsize = popsize)
#                     Xcount.new <- update_Xcount(Xcount.other = Xcount.other, path = path.new)
                    
                    W.new <- updateW(W = W.other, path = path.new)
                    
                    if(max(Xcount.new[,2]) == pathirm.cur[4,4,dim(pathirm.cur)[3]]){
                        
                        new.numinf <- pathirm.cur[4,4,dim(pathirm.cur)[3]]+1
                        
                        pathirm.cur <- update_irm(irm = pathirm.cur, new.numinf = new.numinf, b = Beta[k-1], m = Mu[k-1], a = 0, popsize = popsize)
                        patheigen.cur <- update_eigen(patheigen = patheigen.cur, pathirm = pathirm.cur)
                        
                    } 
                    
                    
                    pop_prob.new <- pop_prob(Xcount = Xcount.new, tmax = tmax, b = Beta[k-1], m = Mu[k-1], a = 0, initdist = initdist, popsize = popsize)
                    
                    path_prob.new <- path_prob(path = path.new, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
                    path_prob.cur <- path_prob(path = path.cur, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
                    
                    a.prob <- accept_prob(pop_prob.new = pop_prob.new, pop_prob.cur = pop_prob.cur, path_prob.cur = path_prob.cur, path_prob.new = path_prob.new)
                    
                    if(min(a.prob, 0) > log(runif(1))) {
                        X.cur <- X.new
                        Xcount.cur <- Xcount.new
                        W.cur <- W.new
                        pop_prob.cur <- pop_prob.new
                    }
                    
                }                
                
            }
            
            # sample new parameters 
            
#             probs[k] <- update_prob(W = W.cur, p.prior = p.prior)
            probs[k] <- probs[k-1]
            
#             params.new <- update_rates(Xcount = Xcount.cur, beta.prior = beta.prior, mu.prior = mu.prior, popsize = popsize)
            
#             Beta[k] <- params.new[1]
#             Mu[k] <- params.new[2]

            Beta[k] <- Beta[k-1]
            Mu[k] <- Mu[k-1]
            
                                
            trajectories[[k]] <- Xcount.cur
            
        }
        
        results2 <- list(Beta = Beta, Mu = Mu, probs=probs_seq[s], loglik = 0, accepts = 0, trajectories = trajectories) 
        
        
        # Plot results -----------------------------------------------------------------
        
        censusInterval <- 0.05; p <- probs_seq[s]
        trajectories2 <- list(); observations2 <- list(); likelihoods <- list()
        
        # for(k in 1:(length(results2[[6]]))){
        for(k in 2:(length(trajectories))){
            if ((k%%50)==0){
                traj <- results2$trajectories[[k]]
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
                
                obs$sampled <- rbinom(dim(obs)[1], obs$truth, prob = probs_seq[s])
                
                observations2[[k]] <- obs
            }
            
        }
        
        trajecs <- do.call(rbind,trajectories2)
        trajecs$infected <- trajecs$infected + rnorm(nrow(trajecs), sd = 0.05)
        samples <- do.call(rbind,observations2)
        
        truetrajec <- data.frame(build_countmat(SIRres$trajectory, popsize), simnum = 0)
        
        sample_wide <- reshape(samples, timevar = "simnum", idvar = c("time"), direction = "wide", drop = "truth")
        sample_wide <- data.frame(time = sample_wide$time[1:(nrow(sample_wide) - 2)], avg = apply(sample_wide[1:(nrow(sample_wide) - 2), 2:ncol(sample_wide)], 1, mean, na.rm=TRUE))
        
        
        print(ggplot(data = trajecs, aes(x = time, y = infected, group = simnum)) + geom_step(alpha = 0.1) + geom_step(data = truetrajec, aes(x = time, y = numsick),colour="red", size = 1) +geom_point(data=data.frame(SIRres$results,simnum=0), aes(x=time,y=Observed),size=4,colour="blue", alpha = 0.6) + labs(title = paste("b = 0.5, m = 1, N=3, p = ",probs_seq[s],", Iter = ",t, sep="")) + theme_bw())
        
    }
}

dev.off()
