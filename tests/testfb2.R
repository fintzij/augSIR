# Testing the path_prob and pop_prob functions for small population size (n = 4, fix 3, sim 1) - CHECKS OUT

# initialize values
Xcount.other <- matrix(c(0,1,2,1,2,1,2,3,0,3,2,0,4,1,0,5,0,0),ncol=3, byrow=T); colnames(Xcount.other) <- c("Time", "Infected", "Susceptible")
tmax <-  10; b <- runif(1); m <- runif(1); a <- 0; initdist <- c(0.95, 0.05, 0); popsize = 4
pathirm <- build_irm(Xcount = Xcount.other, b = b, m = m, a = a, popsize = popsize, pop = FALSE)


# Case 1: Subject not infected by time tmax 
path <- c(0, 0); Xcount.new <- update_Xcount(Xcount.other, path)

#comparison of path probabilities - CHECKS OUT
log(initdist[1]) - sum(b*Xcount.other[1:(dim(Xcount.other)[1] - 1),2]*diff(Xcount.other[,1]))
path_prob(path, Xcount.other, pathirm, initdist, tmax)


#comparison of population probabilities - CHECKS OUT
dmultinom(x = c(Xcount.new[1,3], Xcount.new[1,2], 0), prob = initdist, log = TRUE) + log(b*Xcount.new[1,2]*Xcount.new[1,3]) - (b*Xcount.new[1,2]*Xcount.new[1,3] + m*Xcount.new[1,2]) * diff(Xcount.new[,1])[1] + 
    log(b * Xcount.new[2,2] * Xcount.new[2,3]) - (b*Xcount.new[2,2] * Xcount.new[2,3] + m*Xcount.new[2,2]) * diff(Xcount.new[,1])[2] + 
    log(m * Xcount.new[3,2]) - (m*Xcount.new[3,2] + b*Xcount.new[3,2]*Xcount.new[3,3])*diff(Xcount.new[,1])[3] + 
    log(m * Xcount.new[4,2]) - (m*Xcount.new[4,2] + b*Xcount.new[4,2]*Xcount.new[4,3])*diff(Xcount.new[,1])[4] + 
    log(m * Xcount.new[5,2]) - (m*Xcount.new[5,2] + b*Xcount.new[5,2]*Xcount.new[5,3])*diff(Xcount.new[,1])[5] - 
    (m * Xcount.new[6,2] + b*Xcount.new[6,3]*Xcount.new[6,2])*(tmax - Xcount.new[6, 1])
    

pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)


# case 2: Subject initially infected, but does not recover by tmax 

#comparison of path probabilities
path <- c(0, Inf); Xcount.new <- update_Xcount(Xcount.other, path)

log(initdist[2]) - m*tmax
path_prob(path, Xcount.other, pathirm, initdist, tmax)


#comparison of population probabilities - CHECKS OUT
dmultinom(x = c(Xcount.new[1,3], Xcount.new[1,2], 0), prob = initdist, log = TRUE) + log(b*Xcount.new[1,2]*Xcount.new[1,3]) - (b*Xcount.new[1,2]*Xcount.new[1,3] + m*Xcount.new[1,2]) * diff(Xcount.new[,1])[1] + 
    log(b * Xcount.new[2,2] * Xcount.new[2,3]) - (b*Xcount.new[2,2] * Xcount.new[2,3] + m*Xcount.new[2,2]) * diff(Xcount.new[,1])[2] + 
    log(m * Xcount.new[3,2]) - (m*Xcount.new[3,2] + b*Xcount.new[3,2]*Xcount.new[3,3])*diff(Xcount.new[,1])[3] + 
    log(m * Xcount.new[4,2]) - (m*Xcount.new[4,2] + b*Xcount.new[4,2]*Xcount.new[4,3])*diff(Xcount.new[,1])[4] + 
    log(m * Xcount.new[5,2]) - (m*Xcount.new[5,2] + b*Xcount.new[5,2]*Xcount.new[5,3])*diff(Xcount.new[,1])[5] - 
    (m * Xcount.new[6,2] + b*Xcount.new[6,3]*Xcount.new[6,2])*(tmax - Xcount.new[6, 1])

pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)


# case 3: Subject initially infected and recovers prior to tmax 

# comparison of path probabilities - CHECKS OUT
path <- c(0, runif(1, c(0,tmax))); Xcount.new <- update_Xcount(Xcount.other, path)

log(initdist[2]) + log(m) - m*path[2]
path_prob(path, Xcount.other, pathirm, initdist, tmax)


#comparison of population probailities - CHECKS OUT
dmultinom(x = c(Xcount.new[1,3], Xcount.new[1,2], 0), prob = initdist, log = TRUE) + 
    ifelse(diff(Xcount.new[,2])[1] > 0, log(b*Xcount.new[1,2]*Xcount.new[1,3]), log(m * Xcount.new[1,2])) - (m*Xcount.new[1,2] + b*Xcount.new[1,2]*Xcount.new[1,3])*diff(Xcount.new[,1])[1] + 
    ifelse(diff(Xcount.new[,2])[2] > 0, log(b*Xcount.new[2,2]*Xcount.new[2,3]), log(m * Xcount.new[2,2])) - (m*Xcount.new[2,2] + b*Xcount.new[2,2]*Xcount.new[2,3])*diff(Xcount.new[,1])[2] + 
    ifelse(diff(Xcount.new[,2])[3] > 0, log(b*Xcount.new[3,2]*Xcount.new[3,3]), log(m * Xcount.new[3,2])) - (m*Xcount.new[3,2] + b*Xcount.new[3,2]*Xcount.new[3,3])*diff(Xcount.new[,1])[3] + 
    ifelse(diff(Xcount.new[,2])[4] > 0, log(b*Xcount.new[4,2]*Xcount.new[4,3]), log(m * Xcount.new[4,2])) - (m*Xcount.new[4,2] + b*Xcount.new[4,2]*Xcount.new[4,3])*diff(Xcount.new[,1])[4] + 
    ifelse(diff(Xcount.new[,2])[5] > 0, log(b*Xcount.new[5,2]*Xcount.new[5,3]), log(m * Xcount.new[5,2])) - (m*Xcount.new[5,2] + b*Xcount.new[5,2]*Xcount.new[5,3])*diff(Xcount.new[,1])[5] + 
    ifelse(diff(Xcount.new[,2])[6] > 0, log(b*Xcount.new[6,2]*Xcount.new[6,3]), log(m * Xcount.new[6,2])) - (m*Xcount.new[6,2] + b*Xcount.new[6,2]*Xcount.new[6,3])*diff(Xcount.new[,1])[5] - 
    (m * Xcount.new[7,2] + b*Xcount.new[7,3]*Xcount.new[7,2])*(tmax - Xcount.new[7, 1])

pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)



# case 4: subject not initially infected, becomes infected, and recovery is not observed
path <- c(runif(1, 0, 5), Inf); Xcount.new <- update_Xcount(Xcount.other, path)

#comparison of path probabilities - CHECKS OUT
infec.int <- findInterval(path[1], c(Xcount.other[,1], tmax), all.inside = TRUE); infec.int

log(initdist[1]) - sum(b * Xcount.other[which(1:nrow(Xcount.other) < infec.int), 2] * diff(Xcount.other[,1])[which(1:nrow(Xcount.other)< infec.int)]) + 
    log(b * Xcount.other[infec.int, 2]) - (b * Xcount.other[infec.int, 2] * (path[1] - Xcount.other[infec.int, 1])) - m * (tmax - path[1])

path_prob(path, Xcount.other, pathirm, initdist, tmax)


#comparison of population probabilities - CHECKS OUT
dmultinom(x = c(Xcount.new[1,3], Xcount.new[1,2], 0), prob = initdist, log = TRUE) + 
    ifelse(diff(Xcount.new[,2])[1] > 0, log(b*Xcount.new[1,2]*Xcount.new[1,3]), log(m * Xcount.new[1,2])) - (m*Xcount.new[1,2] + b*Xcount.new[1,2]*Xcount.new[1,3])*diff(Xcount.new[,1])[1] + 
    ifelse(diff(Xcount.new[,2])[2] > 0, log(b*Xcount.new[2,2]*Xcount.new[2,3]), log(m * Xcount.new[2,2])) - (m*Xcount.new[2,2] + b*Xcount.new[2,2]*Xcount.new[2,3])*diff(Xcount.new[,1])[2] + 
    ifelse(diff(Xcount.new[,2])[3] > 0, log(b*Xcount.new[3,2]*Xcount.new[3,3]), log(m * Xcount.new[3,2])) - (m*Xcount.new[3,2] + b*Xcount.new[3,2]*Xcount.new[3,3])*diff(Xcount.new[,1])[3] + 
    ifelse(diff(Xcount.new[,2])[4] > 0, log(b*Xcount.new[4,2]*Xcount.new[4,3]), log(m * Xcount.new[4,2])) - (m*Xcount.new[4,2] + b*Xcount.new[4,2]*Xcount.new[4,3])*diff(Xcount.new[,1])[4] + 
    ifelse(diff(Xcount.new[,2])[5] > 0, log(b*Xcount.new[5,2]*Xcount.new[5,3]), log(m * Xcount.new[5,2])) - (m*Xcount.new[5,2] + b*Xcount.new[5,2]*Xcount.new[5,3])*diff(Xcount.new[,1])[5] + 
    ifelse(diff(Xcount.new[,2])[6] > 0, log(b*Xcount.new[6,2]*Xcount.new[6,3]), log(m * Xcount.new[6,2])) - (m*Xcount.new[6,2] + b*Xcount.new[6,2]*Xcount.new[6,3])*diff(Xcount.new[,1])[6] - 
    (m * Xcount.new[7,2] + b*Xcount.new[7,3]*Xcount.new[7,2])*(tmax - Xcount.new[7, 1])

pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)


# case 5: subject not initially infected, becomes infected, recovers
path <- c(0,0); path[1] <- runif(1,0,5); path[2] <- runif(1, path[1], tmax); Xcount.new <- update_Xcount(Xcount.other, path)

#comparison of path probabilities - CHECKS OUT
infec.int <- findInterval(path[1], c(Xcount.other[,1], tmax), all.inside = TRUE)

log(initdist[1]) - sum(b * Xcount.other[which(1:nrow(Xcount.other) < infec.int), 2] * diff(Xcount.other[,1])[which(1:nrow(Xcount.other)< infec.int)]) + 
    log(b * Xcount.other[infec.int, 2]) - (b * Xcount.other[infec.int, 2] * (path[1] - Xcount.other[infec.int, 1])) + log(m) - m * (path[2] - path[1])

path_prob(path, Xcount.other, pathirm, initdist, tmax)


#comparison of population probabilities - CHECKS OUT
dmultinom(x = c(Xcount.new[1,3], Xcount.new[1,2], 0), prob = initdist, log = TRUE) + 
    ifelse(diff(Xcount.new[,2])[1] > 0, log(b*Xcount.new[1,2]*Xcount.new[1,3]), log(m * Xcount.new[1,2])) - (m*Xcount.new[1,2] + b*Xcount.new[1,2]*Xcount.new[1,3])*diff(Xcount.new[,1])[1] + 
    ifelse(diff(Xcount.new[,2])[2] > 0, log(b*Xcount.new[2,2]*Xcount.new[2,3]), log(m * Xcount.new[2,2])) - (m*Xcount.new[2,2] + b*Xcount.new[2,2]*Xcount.new[2,3])*diff(Xcount.new[,1])[2] + 
    ifelse(diff(Xcount.new[,2])[3] > 0, log(b*Xcount.new[3,2]*Xcount.new[3,3]), log(m * Xcount.new[3,2])) - (m*Xcount.new[3,2] + b*Xcount.new[3,2]*Xcount.new[3,3])*diff(Xcount.new[,1])[3] + 
    ifelse(diff(Xcount.new[,2])[4] > 0, log(b*Xcount.new[4,2]*Xcount.new[4,3]), log(m * Xcount.new[4,2])) - (m*Xcount.new[4,2] + b*Xcount.new[4,2]*Xcount.new[4,3])*diff(Xcount.new[,1])[4] + 
    ifelse(diff(Xcount.new[,2])[5] > 0, log(b*Xcount.new[5,2]*Xcount.new[5,3]), log(m * Xcount.new[5,2])) - (m*Xcount.new[5,2] + b*Xcount.new[5,2]*Xcount.new[5,3])*diff(Xcount.new[,1])[5] + 
    ifelse(diff(Xcount.new[,2])[6] > 0, log(b*Xcount.new[6,2]*Xcount.new[6,3]), log(m * Xcount.new[6,2])) - (m*Xcount.new[6,2] + b*Xcount.new[6,2]*Xcount.new[6,3])*diff(Xcount.new[,1])[6] +
    ifelse(diff(Xcount.new[,2])[7] > 0, log(b*Xcount.new[7,2]*Xcount.new[7,3]), log(m * Xcount.new[7,2])) - (m*Xcount.new[7,2] + b*Xcount.new[7,2]*Xcount.new[7,3])*diff(Xcount.new[,1])[7] -
    (m * Xcount.new[8,2] + b*Xcount.new[8,3]*Xcount.new[8,2])*(tmax - Xcount.new[8, 1])

pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)


# Xcount = Xcount.other; irm = pathirm.cur; irm.eig = patheigen.cur; W = W.other; p = probs[k-1]; initdist = initdist; tmax = tmax
# 
# # # checking the getpath function
# # for(t in 1:50000){
# #     if(t%%500 ==0) print(t)
# #     
# #     Xt <- cbind(W[,1], 0)
# #     
# #     # tpm_seq returns a list whose first element is the sequence of tpms between observation times, and whose second element is a list of tpm subsequences
# #     tpms <- tpm_seq(Xcount = Xcount, obstimes = W[,1], irm.eig = irm.eig)
# #     emits <- emit_seq(W = W, p = p)
# #     
# #     fbmats <- fwd(tpms = tpms[[1]], emits = emits, initdist = initdist)
# #     
# #     Xt[,2] <- bwd(mats = fbmats)
# #     
# #     path.new <- draw_eventtimes(Xt = Xt, Xcount = Xcount, tpm.seqs = tpms[[2]], irm = irm, tmax = tmax)
# #     
# #     X.new <- updateX(X = X.cur, path = path.new, j = subjects[j]); path.new2 <- getpath(X = X.new, j = subjects[j])
# #     
# #     if(any(path.new != path.new2)) break
# # }
# # 
# # 
# # 
# 
# 
# 
# # Test whether simulation for one individual via fb matches simulation via forward sampling
# 
# 
# # Simulation setup --------------------------------------------------------
# 
# # evaluating how well the method reconstructs the true trajectory
# 
# # Simulate data ----------------------------------------------
# 
# SIRres<-SIRsim(popsize = 200, initdist = c(0.9, 0.1, 0), b = 0.005, mu=0.5, a=0, tmax = 25, censusInterval=0.1, sampprob = 0.25, binomsamp = TRUE, returnX = TRUE)
# 
# if(dim(SIRres$results)[1] < 30){
#     while(dim(SIRres$results)[1] < 30){
#         SIRres<-SIRsim(popsize = 200, initdist = c(0.9, 0.1, 0), b = 0.005, mu=0.5, a=0, tmax = 25, censusInterval=0.1, sampprob = 0.25, binomsamp = TRUE, returnX = TRUE)
#         
#     }
# }
# 
# # get data 
# dat <- SIRres$results
# dat.m <- melt(dat,id.vars="time")
# 
# ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + theme_bw()
# 
# 
# sim.settings <- list(popsize = 200,
#                      tmax = max(dat[,1])+5,
#                      niter = 50000,
#                      amplify = 10,
#                      initdist = c(0.9, 0.1, 0))
# 
# inits <- list(beta.init = 0.005 + runif(1,-0.0005, 0.0005),
#               mu.init = 0.5 + runif(1, -0.005, 0.005),
#               alpha.init = 0, 
#               probs.init = 0.25 + runif(1,-0.005, 0.005))
# 
# priors <- list(beta.prior = c(.051, 9.997),
#                mu.prior = c(0.96, 1.96),
#                alpha.prior = NULL,
#                p.prior = c(4.5, 13.35))
# 
# 
# popsize <- sim.settings$popsize # size of the population; 
# tmax <- sim.settings$tmax # maximum time of observation
# amplify <- sim.settings$amplify # amplification parameter for initialization
# niter <- sim.settings$niter # number of iterations in the sampler
# initdist <- sim.settings$initdist # initial distribution for individual infection status
# 
# # vectors for parameters
# Beta <- vector(length=niter); Beta[1] <- inits$beta.init
# Mu <- vector(length = niter); Mu[1] <- inits$mu.init
# Alpha <- vector(length = niter); Alpha[1] <- inits$alpha.init 
# probs <- vector(length = niter); probs[1] <- inits$probs.init
# accepts <- vector(length = niter)
# 
# # vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
# beta.prior <- priors$beta.prior
# mu.prior <- priors$mu.prior
# # alpha.prior <- c(6, 12000)
# p.prior <- priors$p.prior
# 
# # log-likelihood vector
# loglik <- vector(length=niter)
# 
# # list to store trajectories
# trajectories <- list()
# 
# # observation matrix
# W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$Observed, augmented = 0)); 
# 
# # matrix with event times, subject id and event codes. 
# # Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
# X.cur <- SIRres$trajectory
# # X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], p=probs[1], amplify = amplify, tmax=20, popsize = popsize)
# Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
# # update observation matrix
# W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
# 
# if(!checkpossible(X=X.cur, W=W.cur)) {
#     while(!checkpossible(X=X.cur,W=W.cur)){
#         X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], p=probs[1], amplify = amplify, tmax=20, popsize = popsize)
#         Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
#         W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
#     }
# }
# 
# pop_prob.cur <- pop_prob(Xcount = Xcount.cur, b = Beta[1], m = Mu[1], a = Alpha[1], initdist = initdist, popsize = popsize)
# 
# loglik[1] <- calc_loglike(Xcount = Xcount.cur, W = W.cur, b = Beta[1], m = Mu[1], a = Alpha[1], p = probs[1], initdist = initdist, popsize = 200)
# 
# subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
# 
# pathirm.cur <- build_irm(Xcount = Xcount.cur, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = FALSE)
# patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
# 
# accepts.k <- 0 
# 
# Xother <- X.cur[X.cur[,2]!=subjects[j],]
# 
# path.cur <- getpath(X.cur, subjects[j])
# 
# Xcount.other <- get_Xcount_other(Xcount = Xcount.cur, path = path.cur)
# W.other <-get_W_other(W.cur = W.cur, path = path.cur)
# 
# obstimes <- W.cur[,1]
# # Simulation --------------------------------------------------------------
# 
# 
# cl <- makeCluster(3)
# registerDoParallel(cl)
# draws.varinf.sim1 <- foreach(t = 1:niter, .packages='augSIR') %dopar% {
#         
#     path.new <- sim_one_SIR(Xcount = Xcount.other, obstimes = obstimes, b=Beta[1], m=Mu[1], initdist = initdist, tmax = tmax, returnpath=TRUE)
#     pathvec <- ifelse(obstimes < path.new[1], 1, ifelse(obstimes >= path.new[2], 3, 2))
#     
#     Xcount.cur <- update_Xcount(Xcount.other, path.new)
#     pop_prob.cur <- pop_prob(Xcount.cur, b = Beta[1], m = Mu[1], a = 0, initdist = initdist, popsize = popsize)
#         
#     # update X and W matrices with the new marginal sample path
#     W.cur <- updateW(W = W.other, path = path.new)
#     
#     # draw new binomial samples and update the W matrix
#     W.cur[,2] <- rbinom(n = dim(W.cur)[1], size = W.cur[,3], prob = probs[1]) # draw new binomial samples from the update observation matrix
#     
#     W.other <- get_W_other(W.cur = W.cur, path = path.new)
#     
#     fbpath <- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = probs[1], initdist = initdist, tmax = tmax)
#     fbpathvec <- ifelse(obstimes < fbpath[1], 1, ifelse(obstimes >= fbpath[2], 3, 2))
#     
# #     Xcount.new <- update_Xcount(Xcount.other, fbpath)
# #     pop_prob.new <- pop_prob(Xcount.new, b = Beta[1], m = Mu[1], a = 0, initdist = initdist, popsize = popsize)
# #     
# #     path_prob.cur <- path_prob(path = path.new, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
# #     path_prob.new <- path_prob(path = fbpath, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
# #     
# #     a.prob <- accept_prob(pop_prob.new = pop_prob.new, pop_prob.cur = pop_prob.cur, path_prob.cur = path_prob.cur, path_prob.new = path_prob.new)
# #     
# #     accept.fb <- ifelse(min(a.prob, 0) > log(runif(1)), 1, 0)            
# #     
#     # draw the new trajectory
#     list(pathvec, fbpathvec)
# #     list(pathvec, fbpathvec, accept.fb)
#     # draw new trajectory using new binomial samples
# }
# stopCluster(cl)
# 
# 
# method1.draws <- lapply(draws.varinf.sim1, '[[', 1)
# method1.draws <- lapply(method1.draws, function(x) matrix(unlist(x),ncol=dim(W.cur)[1],byrow=TRUE))
# method1.draws <- do.call(rbind,method1.draws)
# method1.draws <- data.frame(cbind(1:dim(method1.draws)[1],method1.draws)); names(method1.draws) <- c("simnum", obstimes)
# 
# method2.draws <- lapply(draws.varinf.sim1, '[[', 2)
# method2.draws <- lapply(method2.draws, function(x) matrix(unlist(x),ncol=dim(W.cur)[1],byrow=TRUE))
# method2.draws <- do.call(rbind,method2.draws)
# method2.draws <- data.frame(cbind(1:dim(method2.draws)[1],method2.draws)); names(method2.draws) <- c("simnum", obstimes)
# # method2.draws <- data.frame(cbind(1:dim(method2.draws)[1],lapply(draws.varinf.sim1,'[[',3),method2.draws)); names(method2.draws) <- c("simnum","accept", obstimes)
# 
# # method2.draws[method2.draws$accept == 0, 3:dim(method2.draws)[2]] <- method1.draws[method2.draws$accept == 0, 2:dim(method1.draws)[2]]
# 
# comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
#                            Time = rep(obstimes,2), 
#                            Susceptible = 0,
#                            Infected = 0, 
#                            Recovered = 0)
# 
# for(r in 1:length(obstimes)){
#     comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
#     comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
#     comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
#     
#     comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
#     comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
#     comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
#     
# }
# 
# summary.melt <- melt(comp.summary, id = c("Method", "Time"))
# ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "Method 1 vs. Method 2 - Full Population")
