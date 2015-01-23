# evaluating how well the method reconstructs the true trajectory

# Simulate data, first for the case where we observe with error ----------------------------------------------

# get data
# SIRres<-SIRsim(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 20, censusInterval=0.25, sampprob = 0.25)
# 
# # get data 
# dat <- SIRres[,c(1:2)]
# dat.m <- melt(SIRres,id.vars="time")
# 
# ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + theme_bw()

SIRres<-SIRsim(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 200, censusInterval=.25, sampprob = 0.25, returnX = TRUE)

if(dim(SIRres$results)[1] < 10){
    while(dim(SIRres$results)[1] < 10){
        SIRres<-SIRsim(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 200, censusInterval=.25, sampprob = 0.25, returnX = TRUE)
        
    }
}

# get data 
dat <- SIRres$results
dat.m <- melt(dat,id.vars="time")

ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + theme_bw()

sim.settings <- list(popsize = 200,
                     tmax = 200,
                     niter = 500,
                     amplify = 10,
                     initdist = c(0.995, 0.005, 0))

inits <- list(beta.init = 0.01 + runif(1,-0.005, 0.005),
              mu.init = 0.5 + runif(1, -0.05, 0.05),
              alpha.init = 0, 
              probs.init = 0.25 + runif(1,-0.05, 0.05))

priors <- list(beta.prior = c(.012, 1.1),
               mu.prior = c(0.96, 1.96),
               alpha.prior = NULL,
               p.prior = c(4.5, 13.35))

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
W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$Observed, augmented = 0)); 

# matrix with event times, subject id and event codes. 
# Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))

X.cur <- SIRres$trajectory
# X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], p=probs[1], amplify = amplify, tmax=20, popsize = popsize)

# update observation matrix
W.cur <- updateW(W.cur,X.cur)

if(!checkpossible(X=X.cur, W=W.cur)) {
    while(!checkpossible(X=X.cur,W=W.cur)){
        X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], p=probs[1], amplify = amplify, tmax=20, popsize = popsize)
        W.cur <- updateW(W.cur,X.cur)
    }
}

trajectories[[1]] <- X.cur

popirm.cur <- buildirm(X.cur, b = Beta[1], m = Mu[1], a = Alpha[1], popsize = popsize, pop = TRUE)
pop_prob.cur <- pop_prob(X.cur, irm = popirm.cur, initdist = initdist, popsize = popsize)

# M-H sampler
for(k in 2:niter){
    # Update trajectories
    print(k)
    subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
    
    pathirm.cur <- buildirm(X = X.cur, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = FALSE)
    patheigen.cur <- irm_decomp(pathirm.cur)
    
    for(j in 1:length(subjects)){
        Xother <- X.cur[X.cur[,2]!=subjects[j],]
        
        path.cur <- getpath(X.cur, subjects[j])
        
        W.other <-updateW(W.cur, Xother)
        Xt <- drawXt(Xother = Xother, irm = pathirm.cur, irm.eig = patheigen.cur, W=W.other, p=probs[k-1], b=Beta[k-1], m=Mu[k-1], a=Alpha[k-1], initdist = initdist)
             
#         durs <- rep(0,200)
#         for(s in 1:200){
#             Xt <- drawXt(Xother = Xother, irm = pathirm.cur, irm.eig = patheigen.cur, W=W.other, p=probs[k-1], b=Beta[k-1], m=Mu[k-1], a=Alpha[k-1], initdist = initdist)
#             path <- drawpath(Xt,Xother, pathirm.cur, tmax)
#             
#             durs[s] <- path[2]-path[1]
#             
#             if(all(path==Inf)) durs[s] <- 0
#             
#             if(is.nan(durs[s])) stop()
#         }
#         
#         durations[j] <- mean(durs[durs<Inf])
#         
#         if(is.na(durations[j])) stop()
        
        path.new<- drawpath(Xt, Xother, pathirm.cur, tmax)
        
        X.new <- updateX(X.cur,path.new,subjects[j]); path.new <- getpath(X.new,subjects[j])
        
        if(max(cumsum(X.new[,3])) == pathirm.cur[4,4,dim(pathirm.cur)[3]]){

            new.numinf <- pathirm.cur[4,4,dim(pathirm.cur)[3]]+1

            pathirm.cur <- update_irm(irm = pathirm.cur, new.numinf = new.numinf, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize)
            patheigen.cur <- update_eigen(patheigen = patheigen.cur, pathirm = pathirm.cur)
        } 
        
        popirm.new <- buildirm(X.new, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = TRUE)
        pop_prob.new <- pop_prob(X.new, irm = popirm.new, initdist = initdist, popsize = popsize)
        
        path_prob.new <- path_prob(path.new, Xother, pathirm.cur, initdist, tmax)
        path_prob.cur <- path_prob(path.cur, Xother, pathirm.cur, initdist, tmax)
        
        a.prob <- accept_prob(pop_prob.new, pop_prob.cur, path_prob.cur, path_prob.new)
        
        if(min(a.prob, 0) > log(runif(1))) {
            X.cur <- X.new
            W.cur <- updateW(W.cur,X.new)
            popirm.cur <- popirm.new
            pop_prob.cur <- pop_prob.new
        }
        
    }
        
    # draw new parameters
    probs[k] <- update_prob(W = W.cur, p.prior = p.prior)
    
    # new rate parameters 
    params.new <- update_rates(X = X.cur, beta.prior = beta.prior, mu.prior = mu.prior, popsize = popsize)
    Beta[k] <- params.new[1]
    
    Mu[k] <- params.new[2]
    
    Alpha[k] <- params.new[3]
    
    loglik[k] <- calc_loglike(X = X.cur, W = W.cur, irm = popirm.cur, b = Beta[k], m = Mu[k], a = Alpha[k], p = probs[k], initdist = initdist, popsize = popsize)  
    
    trajectories[[k]] <- X.cur
    
}

results2 <- list(Beta = Beta, Mu = Mu, probs=probs, loglik = loglik, trajectories = trajectories) 


# Results for the case with error, p=0.2  -----------------------------------------------------------------

censusInterval <- 0.1; p <- 0.25
trajectories2 <- list(); observations2 <- list(); likelihoods <- list()

for(k in 1:(length(results2[[4]]))){
# for(k in 1:(length(trajectories))){
    if ((k%%1)==0){
    traj <- results2[[4]][[k]]
    traj <- trajectories[[k]]
        Xobs <- data.frame(time = unique(traj[,1]), 
                       infected = c(sum(traj[traj[,1]==0,3]),sum(traj[traj[,1]==0,3]) + cumsum(traj[traj[,1]!=0,3])), 
                       simnum = k)
    trajectories2[[k]] <- Xobs
    
    obs <- data.frame(time = seq(0,max(Xobs$time),by=0.25), 
                      truth = 0, 
                      sampled = 0, 
                      simnum = k)
    
    for(j in 1:dim(obs)[1]){
        obs$truth[j] <- sum(traj[traj[,1]<=obs$time[j],3])
    }
    
    obs$sampled <- rbinom(dim(obs)[1], obs$truth, prob = 0.2)
    
    observations2[[k]] <- obs
    }
    
}

trajecs <- do.call(rbind,trajectories2)
samples <- do.call(rbind,observations2)
likelihoods <- do.call(rbind,likelihoods)

truetrajec <- data.frame(time = unique(SIRres$trajectory[,1]), 
                         infected = c(sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]), sum(SIRres$trajectory[SIRres$trajectory[,1]==0,3]) + cumsum(SIRres$trajectory[SIRres$trajectory[,1]!=0,3])),
                         simnum = 0)

trajecs.gg <- ggplot(data = trajecs, aes(x = time, y = infected, group = simnum)) + geom_line(alpha = 0.1) + geom_line(data = truetrajec, aes(x = time, y = infected),colour="red", size = 1) + 
    geom_boxplot(data = samples, aes(x=time, y=sampled, group= time),outlier.shape=NA) + geom_point(data=data.frame(dat,simnum=0), aes(x=time,y=Observed),size=4,colour="blue") + theme_bw()

# # Simulate data, first for the case where we observe without error ----------------------------------------------
# 
# SIRres<-SIRsim(N = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, maxtime = 20,censusInterval=1);
# SIRres = cbind(SIRres,200 - rowSums(SIRres[,2:3]))
# colnames(SIRres)<-c("time","susceptible","infected","recovered")
# 
# # get data and add some noise
# dat <- data.frame(SIRres); dat$infected<-rbinom(n=dim(dat)[1], size=dat$infected, prob=1)
# dat.m <- melt(dat,id.vars="time")
# 
# ggplot(dat.m,aes(x=time,y=value,colour=variable))+geom_point() + theme_bw()
# 
# # initialize simulation settings
# popsize <- 200 # size of the population; 
# tmax <- 20 # maximum time of observation
# 
# niter <- 1000 # number of iterations in the sampler
# 
# initdist <- c(0.995,0.005,0) # initial distribution for individual infection status
# 
# # vectors for parameters
# Beta <- vector(length=niter+1); Beta[1] <- 0.01 + runif(1,-0.0005,0.0005)
# Mu <- vector(length = niter+1); Mu[1] <- 0.5 + runif(1,-0.0005,0.0005)
# Alpha <- vector(length = niter+1); Alpha[1] <- 0 
# probs <- vector(length = niter+1); #probs[1] <- 0.2 + runif(1, -0.1, 0.1)
# probs[1] <-1
# 
# results.noerr <- vector("list",niter)
# # # vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
# # beta.prior <- c(12, 1200)
# # mu.prior <- c(6, 12)
# # alpha.prior <- c(6, 12000)
# # p.prior <- c(9.5, 38)
# 
# # log-likelihood vector
# loglik <- vector(length=niter)
# 
# # observation matrix
# W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$infected, augmented = 0)); 
# if(W.cur[1,2]==0) W.cur[1,2]<-1
# 
# # matrix with event times, subject id and event codes. 
# # Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
# X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
# 
# X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20, popsize = popsize)
# 
# # update observation matrix
# W.cur <- updateW(W.cur,X.cur)
# 
# if(!checkpossible(X=X.cur, W=W.cur)) {
#     while(!checkpossible(X=X.cur,W=W.cur)){
#         X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20, popsize = popsize)
#         W.cur <- updateW(W.cur,X.cur)
#     }
# }
# 
# popirm.cur <- buildirm(X.cur, b = Beta[1], m = Mu[1], a = Alpha[1], pop = TRUE) 
# pop_prob.cur <- pop_prob(X.cur, irm = popirm.cur, initdist = initdist, popsize = popsize)
# 
# # M-H sampler
# for(k in 1:(niter)){
#     print(k)
#     # Update trajectories
#     
#     subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
#     
#     pathirm.cur <- buildirm(X.cur, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = FALSE)
#     patheigen.cur <- irm_decomp(pathirm.cur)
#     
#     for(j in 1:length(subjects)){
#         Xother <- X.cur[X.cur[,2]!=subjects[j],]
#         
#         path.cur <- getpath(X.cur, subjects[j])
#         
#         W.other <-updateW(W.cur, Xother)
#         Xt <- drawXt(Xother = Xother, irm = pathirm.cur, irm.eig = patheigen.cur, W=W.other, p=probs[k-1], b=Beta[k-1], m=Mu[k-1], a=Alpha[k-1], initdist = initdist)
#         
#         path.new<- drawpath(Xt, Xother, pathirm.cur, tmax)
#         
#         X.new <- updateX(X.cur,path.new,subjects[j]); path.new <- getpath(X.new,subjects[j])
#         
#         if(max(cumsum(X.new[,3])) == pathirm.cur[4,4,dim(pathirm.cur)[3]]){
#             new.numinf <- max(cumsum(X.new[,3]))+1
#             pathirm.cur <- update_irm(irm = pathirm.cur, new.numinf = new.numinf, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize)
#             patheigen.cur <- update_eigen(patheigen = patheigen.cur, pathirm = pathirm.cur, ind = new.numinf)
#         } 
#         
#         popirm.new <- buildirm(X.new, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = TRUE)
#         pop_prob.new <- pop_prob(X.new, irm = popirm.new, initdist = initdist, popsize = popsize)
#         
#         path_prob.new <- path_prob(path.new, Xother, pathirm.cur, initdist, tmax)
#         path_prob.cur <- path_prob(path.cur, Xother, pathirm.cur, initdist, tmax)
#         
#         a.prob <- accept_prob(pop_prob.new, pop_prob.cur, path_prob.cur, path_prob.new)
#         
#         if(min(a.prob, 0) > log(runif(1))) {
#             X.cur <- X.new
#             W.cur <- updateW(W.cur,X.cur)
#             popirm.cur <- popirm.new
#             pop_prob.cur <- pop_prob.new
#         }
#     }
#     
#     probs[k] <- rbeta(1,p.prior[1] + sum(W.cur[,2]), p.prior[2] + sum(W.cur[,3]-W.cur[,2]))
#     
#     # draw new rate parameters 
#     Beta[k] <- update_beta(X.cur = X.cur, beta.prior = beta.prior, popsize = popsize)
#     
#     Mu[k] <- update_mu(X.cur = X.cur, mu.prior = mu.prior)
#     
#     #   Alpha[k] <- update_alpha(X.cur = X.cur, alpha.prior = alpha.prior, popsize = popsize)
#     
#     Alpha[k] <- 0
#     
#     loglik[k] <- calc_loglike(X = X.cur, W = W.cur, irm = popirm.cur, b = Beta[k], m = Mu[k], a = Alpha[k], p = probs[k], initdist = initdist, popsize = popsize)  
#     
#     # store results
#     results.noerr[[k]] <- list(trajectories = X.cur, prob = pop_prob(X.new, irm.new))    
# }

# Results for the case without error  -----------------------------------------------------------------
# 
# trajectories <- data.frame(results[[1]][[1]])
# for(j in 2:length(results)){
#     trajectories <- rbind(trajectories, data.frame(results[[j]][[1]]))
# }
# trajectories$simnum <- rep(1:length(results),each=400)
# 
# trajecs <- data.frame(time = c(0,trajectories$time[trajectories$time!=0 & trajectories$simnum==1]),
#                       infected = c(sum(trajectories$event[trajectories$time==0 & trajectories$simnum==1]),
#                                    sum(trajectories$event[trajectories$time==0 & trajectories$simnum==1]) + cumsum(trajectories$event[trajectories$time!=0 & trajectories$simnum==1])),
#                       simnum = 1)
# for(j in 2:length(results)){
#     trajecs <- rbind(trajecs, data.frame(time = c(0,trajectories$time[trajectories$time!=0 & trajectories$simnum==j]),
#                                          infected = c(sum(trajectories$event[trajectories$time==0 & trajectories$simnum==j]),
#                                                       sum(trajectories$event[trajectories$time==0 & trajectories$simnum==j]) + cumsum(trajectories$event[trajectories$time!=0 & trajectories$simnum==j])),
#                                          simnum = j))
# }
# 
# 
# trajecs.gg <- ggplot(data=subset(dat.m,variable=="infected"),aes(x=time,y=value)) + geom_point(size=4,colour="red") + geom_line(data=trajecs,aes(x=time,y=infected, group=simnum),alpha=0.1) + theme_bw()
# 
# dat.infec <-subset(dat.m,variable=="infected")
# trajecs.gg <- ggplot() + geom_line(data=trajecs,aes(x=time,y=infected,group=trajecs$simnum, alpha=0.1),show_guide=FALSE) + geom_point(data=dat.infec, aes(x=time, y=value, colour="red", size=4),show_guide=FALSE)
# 
# 
# # Simulate data, first for the case where we observe with error ----------------------------------------------
# 
# # get data and add some noise
# dat3 <- data.frame(SIRres); dat3$infected<-rbinom(n=dim(dat2)[1], size=dat3$infected, prob=0.8)
# dat.m3 <- melt(dat3,id.vars="time")
# 
# ggplot(dat.m3,aes(x=time,y=value,colour=variable))+geom_point() + theme_bw()
# 
# # initialize simulation settings
# popsize <- 200 # size of the population; 
# tmax <- 20 # maximum time of observation
# niter <- 100 # number of iterations in the sampler
# initdist <- c(0.995,0.005,0) # initial distribution for individual infection status
# 
# # vectors for parameters
# Beta <- vector(length=niter+1); Beta[1] <- 0.01 + runif(1,-0.0005,0.0005)
# Mu <- vector(length = niter+1); Mu[1] <- 0.5 + runif(1,-0.0005,0.0005)
# Alpha <- vector(length = niter+1); Alpha[1] <- 0 
# probs <- vector(length = niter+1); #probs[1] <- 0.2 + runif(1, -0.1, 0.1)
# probs[1] <- 0.8
# 
# results3 <- vector("list",niter)
# # # vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
# # beta.prior <- c(12, 1200)
# # mu.prior <- c(6, 12)
# # alpha.prior <- c(6, 12000)
# # p.prior <- c(9.5, 38)
# 
# # log-likelihood vector
# loglik <- vector(length=niter)
# 
# # observation matrix
# W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$infected, augmented = 0)); 
# if(W.cur[1,2]==0) W.cur[1,2]<-1
# 
# # matrix with event times, subject id and event codes. 
# # Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
# X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
# 
# X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20)
# 
# # update observation matrix
# W.cur <- updateW(W.cur,X.cur)
# 
# if(!checkpossible(X=X.cur, W=W.cur)) {
#     while(!checkpossible(X=X.cur,W=W.cur)){
#         X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20)
#         W.cur <- updateW(W.cur,X.cur)
#     }
# }
# 
# # M-H sampler
# for(k in 1:(niter)){
#     print(k)
#     # Update trajectories
#     
#     subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
#     irm.cur <- buildirm(X.cur, b = Beta[k], m = Mu[k], a = Alpha[k]) 
#     accepts <- rep(0,length(subjects))
#     
#     for(j in 1:length(subjects)){
#         #         print(j)
#         Xother <- X.cur[X.cur[,2]!=subjects[j],]; 
#                     
#         path.cur <- getpath(X.cur, subjects[j])
#         
#         W.other <-updateW(W.cur, Xother)
#         irm.other <- buildirm(Xother, b=Beta[k], m = Mu[k], a=Alpha[k])
#         Xt <- drawXt(Xother = Xother, irm = irm.other, W=W.other, p=probs[k], b=Beta[k], m=Mu[k], a=Alpha[k], initdist = initdist)
#         
#         path.new<- drawpath(Xt, Xother, irm.other, tmax)
#         
#         X.new <- updateX(X.cur,path.new,subjects[j]); path.new <- getpath(X.new,subjects[j])
#         irm.new <- buildirm(X.new, b = Beta[k], m = Mu[k], a = Alpha[k])
#         
#         a.prob <- pop_prob(X.new, irm = irm.new) - pop_prob(X.cur, irm = irm.cur) + 
#             path_prob(path.cur, Xother, irm.other, initdist, tmax) - path_prob(path.new, Xother, irm.other, initdist, tmax)
#         
#         if(min(a.prob, 0) > log(runif(1))) {
#             X.cur <- X.new
#             W.cur <- updateW(W.cur,X.cur)
#             irm.cur <- irm.new
#             accepts[j] <- 1
#         }
#     }
#     
#     # Update observation matrix
#     #     W <- updateW(W, X)
#     
#     # draw new binomial sampling probability parameter
#     #   probs[k+1] <- rbeta(1,p.prior[1] + sum(W[,2]), p.prior[2] + sum(W[,3]-W[,2]))
#     probs[k+1] <- probs[k]
#     
#     # draw new rate parameters 
#     Beta[k+1] <- Beta[k]
#     
#     Mu[k+1] <- Mu[k]
#     
#     #   Alpha[k+1] <- rgamma(1, shape = (alpha.prior[1] + sum(X[,3]==1)),
#     #                        rate = alpha.prior[2] + sum((popsize - cumsum(X[,3]==1))*c(0,X[2:dim(X)[1],1] - X[1:(dim(X)[1]-1),1])*(X[,3]==1)))
#     Alpha[k+1] <- Alpha[k]
#     
#     #     loglik[k] <- calc_loglike(X, W, Beta[k+1], Mu[k+1], Alpha[k+1], probs[k+1])  
#     
#     # store results
#     results3[[k]] <- list(trajectories = X.cur, accepts = mean(accepts), prob = pop_prob(X.new, irm.new))    
# }
# 
# 
# # Results for the case with error, p=0.8  -----------------------------------------------------------------
# 
# trajectories3 <- data.frame(results3[[1]][[1]])
# for(j in 2:length(results3)){
#     trajectories3 <- rbind(trajectories3, data.frame(results3[[j]][[1]]))
# }
# trajectories3$simnum <- rep(1:length(results3),each=400)
# 
# trajecs3 <- data.frame(time = c(0,trajectories3$time[trajectories3$time!=0 & trajectories3$simnum==1]),
#                        infected = c(sum(trajectories3$event[trajectories3$time==0 & trajectories3$simnum==1]),
#                                     sum(trajectories3$event[trajectories3$time==0 & trajectories3$simnum==1]) + cumsum(trajectories3$event[trajectories3$time!=0 & trajectories3$simnum==1])),
#                        simnum = 1)
# for(j in 2:length(results3)){
#     trajecs3 <- rbind(trajecs3, data.frame(time = c(0,trajectories3$time[trajectories3$time!=0 & trajectories3$simnum==j]),
#                                            infected = c(sum(trajectories3$event[trajectories3$time==0 & trajectories3$simnum==j]),
#                                                         sum(trajectories3$event[trajectories3$time==0 & trajectories3$simnum==j]) + cumsum(trajectories3$event[trajectories3$time!=0 & trajectories3$simnum==j])),
#                                            simnum = j))
# }
# 
# trajecs3.gg <- ggplot(data=subset(dat.m3,variable=="infected"),aes(x=time,y=value)) + geom_point(colour="red",size=4) + geom_line(data=trajecs3,aes(x=time,y=infected,group=simnum),alpha=0.2) + theme_bw()
