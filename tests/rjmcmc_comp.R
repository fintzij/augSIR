library(doParallel)

SIRres<-SIRsim(popsize = 4, initdist = c(0.7, 0.3, 0), b = 0.5, mu=1, a=0, tmax = 4, censusInterval=0.05, sampprob = 0.25, returnX = TRUE, trim = FALSE)

# get data 
dat <- SIRres$results; colnames(dat) <- c("time", "Observed", "Truth")
dat.m <- melt(dat,id.vars="time")
# 
ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + ylim(c(0,4)) + xlim(0,4) + theme_bw() 


# Simulation settings -----------------------------------------------------

# parameters
b <- 0.5 + runif(1,-0.0005, 0.0005)
m <- 1 + runif(1, -0.005, 0.005)
p <- 0.25 + runif(1,-0.005, 0.005)
initdist <- c(0.7, 0.3, 0)

# simulation objects
niter <- 1000000
tmax <- max(SIRres$results[,1])
popsize <- length(unique(SIRres$trajec[,2]))
accepts <- rep(0, length(niter)); accepts[1] <- 1
accept_probs <- rep(0, length(niter)); accept_probs[1] <- 0
moves <- matrix(0, nrow = niter, ncol = 2); moves[1,] <- c(0,0)
subjID <- 4 #ID of subject to resample

# rjmcmc settings
insert.prob = 1/3; remove.prob = 1/3; shift.prob = 1/3
shift.int <- 0.1

log_likelihood <- rep(0, length(niter))

trajectories <- vector("list", length(niter))


# run the samplers --------------------------------------------------------

cl <- makeCluster(2)
registerDoParallel(cl)

rjmcmc_comp_results <- foreach(s = 1:2, .packages=c('augSIR')) %dopar% {
    
    if(s == 1){ # rjmcmc
        
        # initialize objects
        accepts <- rep(0, length(niter)); accepts[1] <- 1
        accept_ratios <- rep(0, length(niter)); accept_probs[1] <- 0
        moves <- matrix(0, nrow = niter, ncol = 2); moves[1,] <- c(0,0)
        log_likelihood <- rep(0, length(niter))
        trajectories <- vector("list", length(niter))
        
        # observation matrix
        W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$Observed, augmented = 0))
        
        # individual trajectories
        X.cur <- SIRres$trajectory
        
        # count matrix
        Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
        trajectories[[1]] <- Xcount.cur
                
        # update observation matrix
        W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
        
        # save initial log-likelihood
        log_likelihood[1] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = b, m = m, a = 0, p = p, initdist = initdist, popsize = popsize)
        
        for(k in 2:niter){
            
            path.cur <- X.cur[X.cur[,2] == subjID, ]
            path.new <- rjmcmc_draw(path.cur = path.cur, Xcount.cur, j = subjID, initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, tmax = max(W.cur[,1]), b = b, m = m, p = p)
            
            # update bookkeeping objects
            X.new <- X.cur; X.new[X.new[,2] == subjID, ] <- path.new
            Xcount.new <- build_countmat(X = X.new, popsize = popsize)
            W.new <- updateW(W = W.cur, Xcount = Xcount.new)
            
            # calculate acceptance ratio
            rjmcmc.ratio <- rjmcmc_ratio(W.cur = W.cur, W.new = W.new, X.cur = X.cur, X.new = X.new, Xcount.cur = Xcount.cur, Xcount.new = Xcount.new, path.cur = path.cur, path.new = path.new, initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, tmax = tmax, popsize = popsize)
            
            # decide whether to accept. if accept, update current trajectories
            if(rjmcmc.ratio > log(runif(1))){
                
                accepts[k] <- 1
                X.cur <- X.new
                Xcount.cur <- Xcount.new
                W.cur <- W.new
                
            }
            
            # store log-likelihood, move, acceptance prob, and trajectory
            accept_ratios[k] <- rjmcmc.ratio
            moves[k,] <- c(sum(path.cur[,3] != 0), sum(path.new[,3] != 0))
            log_likelihood[k] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = b, m = m, a = 0, p = p, initdist = initdist, popsize = popsize)
            trajectories[[k]] <- Xcount.cur            
            
        }
        
    } else if(s == 2){ # augSIR
        
        # initialize objects
        accepts <- rep(0, length(niter)); accepts[1] <- 1
        accept_probs <- rep(0, length(niter)); accept_probs[1] <- 0
        moves <- matrix(0, nrow = niter, ncol = 2); moves[1,] <- c(0,0)
        log_likelihood <- rep(0, length(niter))
        trajectories <- vector("list", length(niter))
        
        # observation matrix
        W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$Observed, augmented = 0))
        
        # individual trajectories
        X.cur <- SIRres$trajectory
        
        # count matrix
        Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
        trajectories[[1]] <- Xcount.cur
        
        # update observation matrix
        W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
        
        # save initial log-likelihood
        log_likelihood[1] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = b, m = m, a = 0, p = p, initdist = initdist, popsize = popsize)
        
        pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
        
        # build irm matrices
        pathirm.cur <- build_irm(Xcount = Xcount.cur, b = b, m = m, a = 0, popsize = popsize, pop = FALSE)
        patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
        
        if(popsize != pathirm.cur[4,4,dim(pathirm.cur)[3]]){
            
            while(popsize != pathirm.cur[4,4,dim(pathirm.cur)[3]]){
            
            new.numinf <- pathirm.cur[4,4,dim(pathirm.cur)[3]]+1
            
            pathirm.cur <- update_irm(irm = pathirm.cur, new.numinf = new.numinf, b = b, m = m, a = 0, popsize = popsize)
            patheigen.cur <- update_eigen(patheigen = patheigen.cur, pathirm = pathirm.cur)
            
            }
            
        } 
        
        # get .other matrices
        Xother <- X.cur[X.cur[,2]!=subjID,]
        Xcount.other <- build_countmat(X = Xother, popsize = popsize)
        W.other <- updateW(W = W.cur, Xcount = Xcount.other)
        
        for(j in 2:niter){
                        
            # draw new path
            path.new<- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = p, initdist = initdist, tmax = tmax)
            
            # update X, Xcount, W
            X.new <- updateX(X = X.cur, path = path.new, j = subjID)
            Xcount.new <- build_countmat(X = X.new, popsize = popsize)
            W.new <- updateW(W = W.other, Xcount = Xcount.new)
            
            
            pop_prob.new <- pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
            
            path_prob.new <- path_prob(path = path.new, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
            path_prob.cur <- path_prob(path = path.cur, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
            
            a.prob <- accept_prob(pop_prob.new = pop_prob.new, pop_prob.cur = pop_prob.cur, path_prob.cur = path_prob.cur, path_prob.new = path_prob.new)
            
            if(min(a.prob, 0) > log(runif(1))) {
                X.cur <- X.new
                Xcount.cur <- Xcount.new
                W.cur <- W.new
                pop_prob.cur <- pop_prob.new
                accepts[j] <- 1
            }
            
            # store log-likelihood, move, acceptance prob, and trajectory
            accept_ratios[j] <- a.prob
            moves[j,] <- c(sum(path.cur != 0), sum(path.new != 0))
            log_likelihood[j] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = b, m = m, a = 0, p = p, initdist = initdist, popsize = popsize)
            trajectories[[j]] <- Xcount.cur                
        }
    }
    
}

stopCluster(cl)
save(rjmcmc_comp_results, file = "rjmcmc_comp_results.Rdata")