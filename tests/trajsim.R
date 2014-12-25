# evaluating how well the method reconstructs the true trajectory


# Simulate data, first for the case where we observe without error ----------------------------------------------

SIRres<-SIRsim(N = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, maxtime = 20,censusInterval=1);
SIRres = cbind(SIRres,200 - rowSums(SIRres[,2:3]))
colnames(SIRres)<-c("time","susceptible","infected","recovered")

# get data and add some noise
dat <- data.frame(SIRres); dat$infected<-rbinom(n=dim(dat)[1], size=dat$infected, prob=1)
dat.m <- melt(dat,id.vars="time")

ggplot(dat.m,aes(x=time,y=value,colour=variable))+geom_point() + theme_bw()

# initialize simulation settings
popsize <- 200 # size of the population; 
tmax <- 20 # maximum time of observation

niter <- 1000 # number of iterations in the sampler

initdist <- c(0.995,0.005,0) # initial distribution for individual infection status

# vectors for parameters
Beta <- vector(length=niter+1); Beta[1] <- 0.01 + runif(1,-0.0005,0.0005)
Mu <- vector(length = niter+1); Mu[1] <- 0.5 + runif(1,-0.0005,0.0005)
Alpha <- vector(length = niter+1); Alpha[1] <- 0 
probs <- vector(length = niter+1); #probs[1] <- 0.2 + runif(1, -0.1, 0.1)
probs[1] <-1

results.noerr <- vector("list",niter)
# # vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
# beta.prior <- c(12, 1200)
# mu.prior <- c(6, 12)
# alpha.prior <- c(6, 12000)
# p.prior <- c(9.5, 38)

# log-likelihood vector
loglik <- vector(length=niter)

# observation matrix
W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$infected, augmented = 0)); 
if(W.cur[1,2]==0) W.cur[1,2]<-1

# matrix with event times, subject id and event codes. 
# Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))

X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20)

# update observation matrix
W.cur <- updateW(W.cur,X.cur)

if(!checkpossible(X=X.cur, W=W.cur)) {
    while(!checkpossible(X=X.cur,W=W.cur)){
        X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20)
        W.cur <- updateW(W.cur,X.cur)
    }
}

# M-H sampler
for(k in 1:(niter)){
    print(k)
    # Update trajectories
    
    subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
    irm.cur <- buildirm(X.cur, b = Beta[k], m = Mu[k], a = Alpha[k]) 
    accepts <- rep(0,length(subjects))
    
    for(j in 1:length(subjects)){
#         print(j)
        Xother <- X.cur[X.cur[,2]!=subjects[j],]
        path.cur <- getpath(X.cur, subjects[j])
        
        W.other <-updateW(W.cur, Xother)
        irm.other <- buildirm(Xother, b=Beta[k], m = Mu[k], a=Alpha[k])
        Xt <- drawXt(Xother = Xother, irm = irm.other, W=W.other, p=probs[k], b=Beta[k], m=Mu[k], a=Alpha[k], initdist = initdist)
        
        path.new<- drawpath(Xt, Xother, irm.other, tmax)
        
        X.new <- updateX(X.cur,path.new,subjects[j]); path.new <- getpath(X.new,subjects[j])
        irm.new <- buildirm(X.new, b = Beta[k], m = Mu[k], a = Alpha[k])
        
        a.prob <- pop_prob(X.new, irm = irm.new) - pop_prob(X.cur, irm = irm.cur) + 
            path_prob(path.cur, Xother, irm.other, initdist, tmax) - path_prob(path.new, Xother, irm.other, initdist, tmax)
        
        if(min(a.prob, 0) > log(runif(1))) {
            X.cur <- X.new
            W.cur <- updateW(W.cur,X.cur)
            irm.cur <- irm.new
            accepts[j] <- 1
        }
    }
    
    W.cur
    # Update observation matrix
#     W <- updateW(W, X)
    
    # draw new binomial sampling probability parameter
    #   probs[k+1] <- rbeta(1,p.prior[1] + sum(W[,2]), p.prior[2] + sum(W[,3]-W[,2]))
    probs[k+1] <- probs[k]
    
    # draw new rate parameters 
      Beta[k+1] <- Beta[k]
      
    Mu[k+1] <- Mu[k]
    
    #   Alpha[k+1] <- rgamma(1, shape = (alpha.prior[1] + sum(X[,3]==1)),
    #                        rate = alpha.prior[2] + sum((popsize - cumsum(X[,3]==1))*c(0,X[2:dim(X)[1],1] - X[1:(dim(X)[1]-1),1])*(X[,3]==1)))
    Alpha[k+1] <- Alpha[k]
    
#     loglik[k] <- calc_loglike(X, W, Beta[k+1], Mu[k+1], Alpha[k+1], probs[k+1])  
    
    # store results
    results.noerr[[k]] <- list(trajectories = X.cur, accepts = mean(accepts), prob = pop_prob(X.new, irm.new))    
}


# Results for the case without error  -----------------------------------------------------------------

trajectories <- data.frame(results[[1]][[1]])
for(j in 2:length(results)){
    trajectories <- rbind(trajectories, data.frame(results[[j]][[1]]))
}
trajectories$simnum <- rep(1:length(results),each=400)

trajecs <- data.frame(time = c(0,trajectories$time[trajectories$time!=0 & trajectories$simnum==1]),
                      infected = c(sum(trajectories$event[trajectories$time==0 & trajectories$simnum==1]),
                                   sum(trajectories$event[trajectories$time==0 & trajectories$simnum==1]) + cumsum(trajectories$event[trajectories$time!=0 & trajectories$simnum==1])),
                      simnum = 1)
for(j in 2:length(results)){
    trajecs <- rbind(trajecs, data.frame(time = c(0,trajectories$time[trajectories$time!=0 & trajectories$simnum==j]),
                                         infected = c(sum(trajectories$event[trajectories$time==0 & trajectories$simnum==j]),
                                                      sum(trajectories$event[trajectories$time==0 & trajectories$simnum==j]) + cumsum(trajectories$event[trajectories$time!=0 & trajectories$simnum==j])),
                                         simnum = j))
}


trajecs.gg <- ggplot(data=subset(dat.m,variable=="infected"),aes(x=time,y=value)) + geom_point(size=4,colour="red") + geom_line(data=trajecs,aes(x=time,y=infected, group=simnum),alpha=0.1) + theme_bw()

dat.infec <-subset(dat.m,variable=="infected")
trajecs.gg <- ggplot() + geom_line(data=trajecs,aes(x=time,y=infected,group=trajecs$simnum, alpha=0.1),show_guide=FALSE) + geom_point(data=dat.infec, aes(x=time, y=value, colour="red", size=4),show_guide=FALSE)

# Simulate data, first for the case where we observe with error ----------------------------------------------

# get data
dat2 <- data.frame(SIRres); dat2$infected<-rbinom(n=dim(dat2)[1], size=dat2$infected, prob=0.2)
dat.m2 <- melt(dat2,id.vars="time")

ggplot(dat.m2,aes(x=time,y=value,colour=variable))+geom_point() + theme_bw()

# initialize simulation settings
popsize <- 200 # size of the population; 
tmax <- 20 # maximum time of observation
niter <- 1000 # number of iterations in the sampler
initdist <- c(0.995,0.005,0) # initial distribution for individual infection status

# vectors for parameters
Beta <- vector(length=niter+1); Beta[1] <- 0.01 + runif(1,-0.0005,0.0005)
Mu <- vector(length = niter+1); Mu[1] <- 0.5 + runif(1,-0.0005,0.0005)
Alpha <- vector(length = niter+1); Alpha[1] <- 0 
probs <- vector(length = niter+1); #probs[1] <- 0.2 + runif(1, -0.1, 0.1)
probs[1] <- 0.2
 
results.2 <- vector("list",niter)
# # vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
# beta.prior <- c(12, 1200)
# mu.prior <- c(6, 12)
# alpha.prior <- c(6, 12000)
# p.prior <- c(9.5, 38)

# log-likelihood vector
loglik <- vector(length=niter)

# observation matrix
W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$infected, augmented = 0)); 
if(W.cur[1,2]==0) W.cur[1,2]<-1

# matrix with event times, subject id and event codes. 
# Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))

X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20)

# update observation matrix
W.cur <- updateW(W.cur,X.cur)

if(!checkpossible(X=X.cur, W=W.cur)) {
    while(!checkpossible(X=X.cur,W=W.cur)){
        X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20)
        W.cur <- updateW(W.cur,X.cur)
    }
}

# M-H sampler
for(k in 1:(niter)){
    print(k)
    # Update trajectories
    
    subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
    irm.cur <- buildirm(X.cur, b = Beta[k], m = Mu[k], a = Alpha[k]) 
    accepts <- rep(0,length(subjects))
    
    for(j in 1:length(subjects)){
        #         print(j)
        Xother <- X.cur[X.cur[,2]!=subjects[j],]; 

            
        path.cur <- getpath(X.cur, subjects[j])
        
        W.other <-updateW(W.cur, Xother)
        irm.other <- buildirm(Xother, b=Beta[k], m = Mu[k], a=Alpha[k])
        Xt <- drawXt(Xother = Xother, irm = irm.other, W=W.other, p=probs[k], b=Beta[k], m=Mu[k], a=Alpha[k], initdist = initdist)
        
        path.new<- drawpath(Xt, Xother, irm.other, tmax)
        
        X.new <- updateX(X.cur,path.new,subjects[j]); path.new <- getpath(X.new,subjects[j])
        irm.new <- buildirm(X.new, b = Beta[k], m = Mu[k], a = Alpha[k])
        
        a.prob <- pop_prob(X.new, irm = irm.new) - pop_prob(X.cur, irm = irm.cur) + 
            path_prob(path.cur, Xother, irm.other, initdist, tmax) - path_prob(path.new, Xother, irm.other, initdist, tmax)
        
        if(min(a.prob, 0) > log(runif(1))) {
            X.cur <- X.new
            W.cur <- updateW(W.cur,X.cur)
            irm.cur <- irm.new
            accepts[j] <- 1
        }
    }
    
    # Update observation matrix
    #     W <- updateW(W, X)
    
    # draw new binomial sampling probability parameter
    #   probs[k+1] <- rbeta(1,p.prior[1] + sum(W[,2]), p.prior[2] + sum(W[,3]-W[,2]))
    probs[k+1] <- probs[k]
    
    # draw new rate parameters 
    Beta[k+1] <- Beta[k]
    
    Mu[k+1] <- Mu[k]
    
    #   Alpha[k+1] <- rgamma(1, shape = (alpha.prior[1] + sum(X[,3]==1)),
    #                        rate = alpha.prior[2] + sum((popsize - cumsum(X[,3]==1))*c(0,X[2:dim(X)[1],1] - X[1:(dim(X)[1]-1),1])*(X[,3]==1)))
    Alpha[k+1] <- Alpha[k]
    
    #     loglik[k] <- calc_loglike(X, W, Beta[k+1], Mu[k+1], Alpha[k+1], probs[k+1])  
    
    # store results
    results.2[[k]] <- list(trajectories = X.cur, accepts = mean(accepts), prob = pop_prob(X.new, irm.new))    
}


# Results for the case with error, p=0.2  -----------------------------------------------------------------

trajectories2 <- data.frame(results2[[1]][[1]])
for(j in 2:length(results2)){
    trajectories2 <- rbind(trajectories2, data.frame(results2[[j]][[1]]))
}
trajectories2$simnum <- rep(1:length(results2),each=400)

trajecs2 <- data.frame(time = c(0,trajectories2$time[trajectories2$time!=0 & trajectories2$simnum==1]),
                      infected = c(sum(trajectories2$event[trajectories2$time==0 & trajectories2$simnum==1]),
                                   sum(trajectories2$event[trajectories2$time==0 & trajectories2$simnum==1]) + cumsum(trajectories2$event[trajectories2$time!=0 & trajectories2$simnum==1])),
                      simnum = 1)
for(j in 2:length(results2)){
    trajecs2 <- rbind(trajecs2, data.frame(time = c(0,trajectories2$time[trajectories2$time!=0 & trajectories2$simnum==j]),
                                         infected = c(sum(trajectories2$event[trajectories2$time==0 & trajectories2$simnum==j]),
                                                      sum(trajectories2$event[trajectories2$time==0 & trajectories2$simnum==j]) + cumsum(trajectories2$event[trajectories2$time!=0 & trajectories2$simnum==j])),
                                         simnum = j))
}

trajecs2.gg <- ggplot(data=subset(dat.m2,variable=="infected"),aes(x=time,y=value)) + geom_point(colour="red",size=4) + geom_line(data=trajecs2,aes(x=time,y=infected,group=simnum),alpha=0.2) + theme_bw()



# Simulate data, first for the case where we observe with error ----------------------------------------------

# get data and add some noise
dat3 <- data.frame(SIRres); dat3$infected<-rbinom(n=dim(dat2)[1], size=dat3$infected, prob=0.8)
dat.m3 <- melt(dat3,id.vars="time")

ggplot(dat.m3,aes(x=time,y=value,colour=variable))+geom_point() + theme_bw()

# initialize simulation settings
popsize <- 200 # size of the population; 
tmax <- 20 # maximum time of observation
niter <- 100 # number of iterations in the sampler
initdist <- c(0.995,0.005,0) # initial distribution for individual infection status

# vectors for parameters
Beta <- vector(length=niter+1); Beta[1] <- 0.01 + runif(1,-0.0005,0.0005)
Mu <- vector(length = niter+1); Mu[1] <- 0.5 + runif(1,-0.0005,0.0005)
Alpha <- vector(length = niter+1); Alpha[1] <- 0 
probs <- vector(length = niter+1); #probs[1] <- 0.2 + runif(1, -0.1, 0.1)
probs[1] <- 0.8

results3 <- vector("list",niter)
# # vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
# beta.prior <- c(12, 1200)
# mu.prior <- c(6, 12)
# alpha.prior <- c(6, 12000)
# p.prior <- c(9.5, 38)

# log-likelihood vector
loglik <- vector(length=niter)

# observation matrix
W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$infected, augmented = 0)); 
if(W.cur[1,2]==0) W.cur[1,2]<-1

# matrix with event times, subject id and event codes. 
# Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))

X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20)

# update observation matrix
W.cur <- updateW(W.cur,X.cur)

if(!checkpossible(X=X.cur, W=W.cur)) {
    while(!checkpossible(X=X.cur,W=W.cur)){
        X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20)
        W.cur <- updateW(W.cur,X.cur)
    }
}

# M-H sampler
for(k in 1:(niter)){
    print(k)
    # Update trajectories
    
    subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
    irm.cur <- buildirm(X.cur, b = Beta[k], m = Mu[k], a = Alpha[k]) 
    accepts <- rep(0,length(subjects))
    
    for(j in 1:length(subjects)){
        #         print(j)
        Xother <- X.cur[X.cur[,2]!=subjects[j],]; 
                    
        path.cur <- getpath(X.cur, subjects[j])
        
        W.other <-updateW(W.cur, Xother)
        irm.other <- buildirm(Xother, b=Beta[k], m = Mu[k], a=Alpha[k])
        Xt <- drawXt(Xother = Xother, irm = irm.other, W=W.other, p=probs[k], b=Beta[k], m=Mu[k], a=Alpha[k], initdist = initdist)
        
        path.new<- drawpath(Xt, Xother, irm.other, tmax)
        
        X.new <- updateX(X.cur,path.new,subjects[j]); path.new <- getpath(X.new,subjects[j])
        irm.new <- buildirm(X.new, b = Beta[k], m = Mu[k], a = Alpha[k])
        
        a.prob <- pop_prob(X.new, irm = irm.new) - pop_prob(X.cur, irm = irm.cur) + 
            path_prob(path.cur, Xother, irm.other, initdist, tmax) - path_prob(path.new, Xother, irm.other, initdist, tmax)
        
        if(min(a.prob, 0) > log(runif(1))) {
            X.cur <- X.new
            W.cur <- updateW(W.cur,X.cur)
            irm.cur <- irm.new
            accepts[j] <- 1
        }
    }
    
    # Update observation matrix
    #     W <- updateW(W, X)
    
    # draw new binomial sampling probability parameter
    #   probs[k+1] <- rbeta(1,p.prior[1] + sum(W[,2]), p.prior[2] + sum(W[,3]-W[,2]))
    probs[k+1] <- probs[k]
    
    # draw new rate parameters 
    Beta[k+1] <- Beta[k]
    
    Mu[k+1] <- Mu[k]
    
    #   Alpha[k+1] <- rgamma(1, shape = (alpha.prior[1] + sum(X[,3]==1)),
    #                        rate = alpha.prior[2] + sum((popsize - cumsum(X[,3]==1))*c(0,X[2:dim(X)[1],1] - X[1:(dim(X)[1]-1),1])*(X[,3]==1)))
    Alpha[k+1] <- Alpha[k]
    
    #     loglik[k] <- calc_loglike(X, W, Beta[k+1], Mu[k+1], Alpha[k+1], probs[k+1])  
    
    # store results
    results3[[k]] <- list(trajectories = X.cur, accepts = mean(accepts), prob = pop_prob(X.new, irm.new))    
}


# Results for the case with error, p=0.8  -----------------------------------------------------------------

trajectories3 <- data.frame(results3[[1]][[1]])
for(j in 2:length(results3)){
    trajectories3 <- rbind(trajectories3, data.frame(results3[[j]][[1]]))
}
trajectories3$simnum <- rep(1:length(results3),each=400)

trajecs3 <- data.frame(time = c(0,trajectories3$time[trajectories3$time!=0 & trajectories3$simnum==1]),
                       infected = c(sum(trajectories3$event[trajectories3$time==0 & trajectories3$simnum==1]),
                                    sum(trajectories3$event[trajectories3$time==0 & trajectories3$simnum==1]) + cumsum(trajectories3$event[trajectories3$time!=0 & trajectories3$simnum==1])),
                       simnum = 1)
for(j in 2:length(results3)){
    trajecs3 <- rbind(trajecs3, data.frame(time = c(0,trajectories3$time[trajectories3$time!=0 & trajectories3$simnum==j]),
                                           infected = c(sum(trajectories3$event[trajectories3$time==0 & trajectories3$simnum==j]),
                                                        sum(trajectories3$event[trajectories3$time==0 & trajectories3$simnum==j]) + cumsum(trajectories3$event[trajectories3$time!=0 & trajectories3$simnum==j])),
                                           simnum = j))
}

trajecs3.gg <- ggplot(data=subset(dat.m3,variable=="infected"),aes(x=time,y=value)) + geom_point(colour="red",size=4) + geom_line(data=trajecs3,aes(x=time,y=infected,group=simnum),alpha=0.2) + theme_bw()
