
# augSIR - Main Wrapper ---------------------------------------------------

augSIR <- function(dat, sim.settings, priors, inits, returnX = FALSE) {
    
    # initialize simulation settings
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
    accepts <- vector(length = (niter - 1))
    
    # vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
    beta.prior <- priors$beta.prior 
    mu.prior <- priors$mu.prior
    # alpha.prior <- c(6, 12000)
    p.prior <- priors$p.prior
    
    # log-likelihood vector
    loglik <- vector(length=niter)
    
    # observation matrix
    W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$Observed, augmented = 0)); 
    
    # matrix with event times, subject id and event codes. 
    # Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range    
    X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], p=probs[1], amplify = amplify, tmax=20, popsize = popsize)
    Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
    
    # update observation matrix
    W.new <- updateW(W = W.other, path = path.new)
    
    if(returnX == TRUE) {
        trajectories <- list(length = niter)
        trajectories[[1]] <- X.cur
    }
    
    if(!checkpossible(X=X.cur, W=W.cur)) {
        while(!checkpossible(X=X.cur,W=W.cur)){
            X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], p=probs[1], amplify = amplify, tmax=20, popsize = popsize)
            Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
            W.cur <- updateW(W = W.cur, X = X.cur)
        }
    }
    
    popirm.cur <- build_irm(Xcount = Xcount.cur, b = Beta[1], m = Mu[1], a = Alpha[1], popsize = popsize, pop = TRUE)
    pop_prob.cur <- pop_prob(Xcount = Xcount.cur, irm = popirm.cur, initdist = initdist, popsize = popsize)
    
    # M-H sampler
    for(k in 2:niter){
        # Update trajectories
        print(k)
        subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
        
        pathirm.cur <- build_irm(Xcount = Xcount.cur, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = FALSE)
        patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
        
        accepts.k <- 0
        
        for(j in 1:length(subjects)){
            
            Xother <- X.cur[X.cur[,2]!=subjects[j],]
            
            path.cur <- getpath(X.cur, subjects[j])
            
            Xcount.other <- get_Xcount_other(Xcount = Xcount.cur, path = path.cur)
            W.other <-get_W_other(W.cur = W.cur, path = path.cur)
            
            path.new<- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = probs[k-1], initdist = initdist, tmax = tmax)
            
            X.new <- updateX(X = X.cur, path = path.new, j = subjects[j]); path.new <- getpath(X = X.new, j = subjects[j])
            Xcount.new <- update_Xcount(Xcount.other = Xcount.other, path = path.new)
            
            W.new <- updateW(W = W.other, path = path.new)
            
            if(max(Xcount.new[,2]) == pathirm.cur[4,4,dim(pathirm.cur)[3]]){
                
                new.numinf <- pathirm.cur[4,4,dim(pathirm.cur)[3]]+1
                
                pathirm.cur <- update_irm(irm = pathirm.cur, new.numinf = new.numinf, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize)
                patheigen.cur <- update_eigen(patheigen = patheigen.cur, pathirm = pathirm.cur)
                
            } 
            
            popirm.new <- build_irm(Xcount = Xcount.new, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = TRUE)
            pop_prob.new <- pop_prob(Xcount = Xcount.new, irm = popirm.new, initdist = initdist, popsize = popsize)
            
            path_prob.new <- path_prob(path = path.new, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
            path_prob.cur <- path_prob(path = path.cur, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
            
            obs_prob.new <- obs_prob(W.new, probs[k-1])
            obs_prob.cur <- obs_prob(W.cur, probs[k-1])
            
            a.prob <- accept_prob(pop_prob.new = pop_prob.new, pop_prob.cur = pop_prob.cur, path_prob.cur = path_prob.cur, path_prob.new = path_prob.new)
            
            if(min(a.prob, 0) > log(runif(1))) {
                X.cur <- X.new
                Xcount.cur <- Xcount.new
                W.cur <- W.new
                popirm.cur <- popirm.new
                pop_prob.cur <- pop_prob.new
                accepts.k <- accepts.k + 1
            }
            
        }
        
        # save proportion accepted
        accepts[k-1] <- mean(accepts.k)
        
        # draw new parameters
        probs[k] <- update_prob(W = W.cur, p.prior = p.prior)
        
        # new rate parameters 
        params.new <- update_rates(Xcount = Xcount.cur, beta.prior = beta.prior, mu.prior = mu.prior, popsize = popsize)
        Beta[k] <- params.new[1]
        
        Mu[k] <- params.new[2]
        
        Alpha[k] <- params.new[3]
        
        loglik[k] <- calc_loglike(Xcount = Xcount.cur, W = W.cur, irm = popirm.cur, b = Beta[k], m = Mu[k], a = Alpha[k], p = probs[k], initdist = initdist, popsize = popsize)  
        
        if(returnX == TRUE) {
            trajectories[[k]] <- X.cur
        }
        
    }
    
    if(returnX == FALSE){
        results <- list(Beta = Beta, Mu = Mu, probs = probs, loglik = loglik)
        
    } else if(returnX == TRUE){
        results <- list(Beta = Beta, Mu = Mu, probs = probs, loglik = loglik, trajectories=trajectories)
        
    }
    
    return(results)
    
}

