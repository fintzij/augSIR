library(doParallel)

SIRres<-SIRsim(popsize = 4, initdist = c(0.7, 0.3, 0), b = 0.5, mu=1, a=0, tmax = 4, censusInterval=0.05, sampprob = 0.5, returnX = TRUE, trim = FALSE)

# get data 
dat <- SIRres$results; colnames(dat) <- c("time", "Observed", "Truth")
dat.m <- melt(dat,id.vars="time")
# 
ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + ylim(c(0,4)) + xlim(0,4) + theme_bw() 


# Simulation settings -----------------------------------------------------

# parameters
b <- 0.5 + runif(1,-0.0005, 0.0005)
m <- 1 + runif(1, -0.005, 0.005)
samp_prob <- 0.5 + runif(1,-0.005, 0.005)
initdist <- c(0.7, 0.3, 0)

# simulation objects
niter <- 500000
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

# cl <- makeCluster(2)
# registerDoParallel(cl)

rjmcmc_results <- list()

# initialize objects
accepts <- rep(0, length(niter)); accepts[1] <- 1
accept_ratios <- rep(0, length(niter)); accept_probs[1] <- 0
moves <- matrix(0, nrow = niter, ncol = 2); moves[1,] <- c(0,0)
log_likelihood <- rep(0, length(niter))
trajectories <- vector("list", length(niter))
paths <- vector("list", length(niter))
obs_mats <- vector('list', length(niter))

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
log_likelihood[1] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = b, m = m, a = 0, p = samp_prob, initdist = initdist, popsize = popsize)

for(k in 2:niter){
    
    if(k %% 10000 == 0) print(k)
    
    path.cur <- X.cur[X.cur[,2] == subjID, ]
    path.new <- rjmcmc_draw(path.cur = path.cur, Xcount.cur, j = subjID, initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, tmax = max(W.cur[,1]), b = b, m = m, p = samp_prob)
    
    # update bookkeeping objects
    X.new <- X.cur; X.new[X.new[,2] == subjID, ] <- path.new
    Xcount.new <- build_countmat(X = X.new, popsize = popsize)
    W.new <- updateW(W = W.cur, Xcount = Xcount.new)
    
    # calculate acceptance ratio
    rjmcmc.ratio <- rjmcmc_ratio(W.cur = W.cur, W.new = W.new, X.cur = X.cur, X.new = X.new, Xcount.cur = Xcount.cur, Xcount.new = Xcount.new, path.cur = path.cur, path.new = path.new, initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, samp_prob = samp_prob, tmax = tmax, popsize = popsize)
    
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
    log_likelihood[k] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = b, m = m, a = 0, p = samp_prob, initdist = initdist, popsize = popsize)
    trajectories[[k]] <- Xcount.cur 
    paths[[k]] <- path.cur
    obs_mats[[k]] <- W.cur
    
}

rjmcmc_results <- list(trajectories, paths, obs_mats, accepts, accept_ratios, log_likelihood, moves)

save(rjmcmc_results, file = "rjmcmc_results.Rdata")


# augSIR draws ------------------------------------------------------------

# initialize objects
accepts <- rep(0, length(niter)); accepts[1] <- 1
accept_probs <- rep(0, length(niter)); accept_probs[1] <- 0
moves <- matrix(0, nrow = niter, ncol = 2); moves[1,] <- c(0,0)
log_likelihood <- rep(0, length(niter))
trajectories <- vector("list", length(niter))
paths <- vector("list", length(niter))

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
log_likelihood[1] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = b, m = m, a = 0, p = samp_prob, initdist = initdist, popsize = popsize)

pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)

# build irm matrices
pathirm.cur <- build_irm(Xcount = Xcount.cur, b = b, m = m, a = 0, popsize = popsize, pop = FALSE)
patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)

if(popsize > pathirm.cur[4,4,dim(pathirm.cur)[3]]){
    
    while(popsize > pathirm.cur[4,4,dim(pathirm.cur)[3]]){
        
        new.numinf <- pathirm.cur[4,4,dim(pathirm.cur)[3]]+1
        
        pathirm.cur <- update_irm(irm = pathirm.cur, new.numinf = new.numinf, b = b, m = m, a = 0, popsize = popsize)
        patheigen.cur <- update_eigen(patheigen = patheigen.cur, pathirm = pathirm.cur)
        
    }
    
} 

# get .other matrices
Xother <- X.cur[X.cur[,2]!=subjID,]
Xcount.other <- build_countmat(X = Xother, popsize = popsize)
W.other <- updateW(W = W.cur, Xcount = Xcount.other)

paths[[1]] <- X.cur[X.cur[,2] == subjID, c(1,3)]
path.cur <- getpath(X.cur, subjID)

for(j in 2:niter){
    
    if(j %% 10000 == 0) print(j)
    
    # draw new path
    path.new<- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = samp_prob, initdist = initdist, tmax = tmax)
    
    # update X, Xcount, W
    X.new <- updateX(X = X.cur, path = path.new, j = subjID)
    Xcount.new <- build_countmat(X = X.new, popsize = popsize)
    W.new <- updateW(W = W.other, Xcount = Xcount.new)
    
    moves[j,] <- c(sum(X.cur[X.cur[,2] == subjID, 3] != 0), sum(X.new[X.new[,2] == subjID, 3] != 0))
    
    
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
        path.cur <- path.new
    }
    
    # store log-likelihood, move, acceptance prob, and trajectory
    accept_ratios[j] <- a.prob
    log_likelihood[j] <- calc_loglike(Xcount = Xcount.cur, tmax = tmax, W = W.cur, b = b, m = m, a = 0, p = samp_prob, initdist = initdist, popsize = popsize)
    trajectories[[j]] <- Xcount.cur       
    paths[[j]] <- X.cur[X.cur[,2] == subjID, c(1,3)]
    obs_mats[[j]] <- W.cur
}

augSIR_results_rjmcmccomp <- list(trajectories, paths, obs_mats, accept_ratios, log_likelihood, moves)

save(augSIR_results_rjmcmccomp, file = "augSIR_rjmcmccomp.Rdata")


#### Plot results - trajectories
censusInterval <- 0.05

trajectories2 <- list(); observations2 <- list()

# for(k in 1:(length(results2[[6]]))){
for(k in 2:(length(trajectories))){
    if ((k%%500)==0){
        traj <- trajectories[[k]]
        Xobs <- data.frame(time = traj[,1], 
                           infected = traj[,2], 
                           simnum = k)
        trajectories2[[k]] <- Xobs
        
    }
    
}

trajecs <- do.call(rbind,trajectories2)
trajecs$infected <- trajecs$infected + rnorm(nrow(trajecs), sd = 0.05)

truetrajec <- data.frame(build_countmat(SIRres$trajectory, popsize), simnum = 0)

ggplot(data = trajecs, aes(x = time, y = infected, group = simnum)) + geom_step(alpha = 0.1) + geom_step(data = truetrajec, aes(x = time, y = numsick),colour="red", size = 1) +geom_point(data=data.frame(SIRres$results,simnum=0), aes(x=time,y=Observed),size=4,colour="blue", alpha = 0.6) + labs(title = "rjmcmc trajecs") + theme_bw()


# plot average infection status at each time point

path_comp <- data.frame(time = seq(0,4,by=0.05), infected = 0, method = rep(c("augSIR", "rjmcmc"), each = length(seq(0,4,by=0.05))))

augSIR_obsmats <- do.call(rbind,augSIR_results_rjmcmccomp[[3]])
rjmcmc_obsmats <- do.call(rbind, rjmcmc_results[[3]])

for(k in 1:81){
    print(k)
    path_comp$infected[k] <- mean(augSIR_obsmats[augSIR_obsmats[,1] == path_comp$time[k], 3])
    
    path_comp$infected[k + 81] <- mean(rjmcmc_obsmats[rjmcmc_obsmats[,1] == path_comp$time[k], 3])    
    
}

ggplot(path_comp, aes(x=time, y = infected, colour = method)) + geom_line()