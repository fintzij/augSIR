
# Functions for drawing augmented trajectories ----------------------------

# drawpath draws a path for an endpoint conditioned markov chain
# Xt is the subject level observation matrix with first column a vector of observation times, and second column a vector of infection
# statuses. Xother is the matrix with times of carriage acquisition and clearance for all other subjects. 
# irm is an array of instantaneous rate matrices. tmax is the maximum time of observation.

drawpath <- function(Xt, Xother, irm, tmax){
    path <- c(0,0)
    changetimes <- which(Xt[2:dim(Xt)[1],2] - Xt[1:(dim(Xt)[1]-1),2] !=0)
    
    if(length(changetimes)==0){
        if(all(Xt[,2]==2)){
            path[1] <- 0; t0 <- Xt[dim(Xt)[1],1]; t1 <- Inf
            path[2] <- drawtime(Xother, irm, t0, t1, 2)
            path[2] <- ifelse(path[2]>tmax, Inf, path[2])
            
        } else if(all(Xt[,2]==1)){
            t0 <- Xt[dim(Xt)[1],1]; t1 <- Inf
            path[1] <- drawtime(Xother, irm, t0, t1, 1)
            
            if(path[1]>tmax) {
                path[1] <- path[2] <- Inf
            } else if(path[1] <= tmax){
                path[2] <- drawtime(Xother, irm, path[1], t1, 2)
                path[2] <- ifelse(path[2]>tmax,Inf,path[2])
            }
        }
        
    } else if(length(changetimes)==1){
        if(Xt[changetimes+1,2] - Xt[changetimes,2] ==2){
            t0 <- Xt[changetimes,1]; t1 <- Xt[changetimes+1,1]
            path[1] <- drawtime(Xother, irm, t0, t1, 1)
            path[2] <- ifelse(path[1]==Inf,Inf,drawtime(Xother, irm, path[1], t1, 2))
            
        } else if(Xt[changetimes,2]==1 & (Xt[changetimes+1,2] - Xt[changetimes,2] == 1)){
            t0 <- Xt[changetimes,1]; t1 <- Xt[changetimes+1,1]
            path[1] <- drawtime(Xother, irm, t0, t1, 1)
            
            t0 <- Xt[dim(Xt)[1],1]; t1 <- Inf
            path[2] <- drawtime(Xother, irm, t0, t1, 2)
            path[2] <- ifelse(path[2]>tmax,Inf,path[2])
            
        } else if(Xt[changetimes,2]==2 & (Xt[changetimes+1,2] - Xt[changetimes,2] ==1)){
            path[1] <- 0; t0 <- Xt[changetimes,1]; t1 <- Xt[changetimes+1,1]
            path[2] <- drawtime(Xother, irm, t0, t1, 2)
            
        }
        
    } else if(length(changetimes)==2){
        t0 <- Xt[changetimes[1],1]; t1 <- Xt[changetimes[1]+1,1]
        path[1] <- drawtime(Xother, irm, t0, t1, 1)
        
        t0 <- Xt[changetimes[2],1]; t1 <- Xt[changetimes[2]+1,1]
        path[2] <- drawtime(Xother, irm, t0, t1, 2)
    }  
    
    return(path)
}


# drawtime draws an event time for either carriage acquisition or clearance. Xother is a matrix with the times of acquisition and clearance
# for all other subjects. irm is an array of instantaneous rate matrices. t0 and t1 are the left and right endpoints of the interval. 
# current state is the infection status at t0. 

drawtime <- function(Xother, irm, t0, t1, currentstate){
    times <- unique(Xother[,1])
    
    if(t1 == Inf){
        if(currentstate == 2){
            ind1 <- dim(irm)[3]
            eventtime <- t0 - log(1-runif(1))/irm[2,3,ind1]
            
        } else if(currentstate == 1){
            timeseq <- c(t0, Filter(function(X) X>t0, Xother[,1]));
            if(length(timeseq) == 1) {
                eventtime <- Inf
                
            } else if (length(timeseq) > 1){
                for(k in 1:(length(timeseq)-1)){
                    ind <- sum(times<=timeseq[k])
                    eventtime<- timeseq[k] - log(1-runif(1))/irm[1,2,ind]
                    if(eventtime < timeseq[k+1]) break
                }
                
                if (eventtime>timeseq[length(timeseq)]) {
                    eventtime <- Inf
                }
                
            }
        }
        
    } else if(t1 != Inf){
        if(currentstate == 2){
            eventtime <- t0 - log(1 - runif(1)*(1-exp(-irm[2,3,1] * (t1 - t0))))/irm[2,3,1]
            
        } else if(currentstate == 1){
            if (t0 > max(Xother[,1])){
                eventtime <- Inf
                
            } else if(t0 <= max(Xother[,1])){
                if(t1 < max(Xother[,1])){
                    timeseq <- c(t0, Filter(function(X) X>t0 & X<t1, Xother[,1]), t1)
                } else if(t1 >= max(Xother[,1])){
                    timeseq <- c(t0, Filter(function(X) X>t0, Xother[,1]))
                }
                
                intervalprobs <- rep(0,length(timeseq)-1)
                timediffs <- timeseq[2:length(timeseq)] - timeseq[1:(length(timeseq)-1)]
                indstart <- sum(times<=timeseq[1]); indend <- sum(times<=timeseq[length(timeseq)-1])
                totalprob <- 1 - exp(-sum(irm[1, 2,indstart:indend]*timediffs))
                
                for(k in 1:length(intervalprobs)){
                    ind <- sum(times<=timeseq[k])
                    intervalprobs[k] <- (1-exp(-sum(irm[1,2,indstart:ind]*timediffs[1:k])))/totalprob
                }
                
                interval <- Position(function(x) x>runif(1),intervalprobs) 
                ind <- sum(times <= timeseq[interval])
                
                eventtime<-timeseq[interval]-log(1-runif(1)*(1-exp(-irm[1, 2,ind]*timediffs[interval])))/irm[1, 2,ind]
            }
        }
    }
    
    return(eventtime)
}


# drawXt draws the status of a subject at observation times t1,...,tL conditional on other individuals
# Xother is a matrix with the times of carriage aquisition and clearance for all other subjects. 
# p, is the sampling probability, and b, and m are the epidemiologic parameters. initdist is the initial
# distribution for infection status
drawXt <- function(Xother, irm, W, p, b, m, a, initdist){
    
    Xt <- cbind(W[,1], 0); Xt.fb <- array(0,dim=c(3,3,dim(Xt)[1]-1))
    
#     if (sum(Xother[,3][Xother[,1]==0])<W[1,2]) {
#         initdist <- c(0,1,0)
#     } 
    
    Xt.fb[,,1] <- fwd(Xother, W, irm, distr=initdist, p=p, t0=W[1,1], t1 = W[2,1])
    
    for(k in 2:(dim(Xt.fb)[3])){
        Xt.fb[,,k] <- fwd(Xother, W = W, irm=irm, P.prev=Xt.fb[,,k-1], p=p[1], t0=W[k,1], t1=W[k+1,1])
    }
    
    Xt<-bwd(Xt.fb,Xt)
    return(Xt)
}


# Forward-backward recursion functions ------------------------------------

#fwd produces P_tj from P_t,j-1
fwd <- function(Xother, W, irm, P.prev=NULL, distr=NULL, p, t0, t1){
    P.now <- matrix(0, nrow=3,ncol=3);
    if(is.null(distr)) distr <- colSums(P.prev)
    tpm <- obstpm(Xother,irm,t0,t1)
    
    emit <- dbinom(W[W[,1]==t1,2],sum(Xother[which(Xother[,1]<=t1),3])+(1:3==2),p)
    
    P.now<-outer(distr,emit,FUN="*")*tpm
    if(all(P.now==0)){
        P.now <- diag(1,3)
    } else { 
        P.now <- normalize(P.now)
    }
    return(P.now)
}

# bwd implements the stochastic backward recursion, drawing the infection status at time t from the 
# distribution proportional to the column corresponding to the infection status at time t+1

bwd <- function(Xt.fb,Xt){
    states <- rep(0,dim(Xt)[1]); draws <- runif(dim(Xt)[1])
    initdist<-cumsum(colSums(Xt.fb[,,dim(Xt.fb)[3]]))
    currentstate<- Position(function(X) X>=draws[length(states)],initdist)
    states[length(states)] <- currentstate
    for(k in (dim(Xt)[1] - 1):1){
        dist<-cumsum(normalize(Xt.fb[,currentstate,k]))
        currentstate<- Position(function(X) X>=draws[k],dist)
        states[k] <- currentstate
    }
    Xt[,2] <- states
    return(Xt)
}


# Functions for constructing various matrices -----------------------------

# Normalize is what it sounds like.
normalize <- function(X){
    X/sum(X)
}

# buildirm constructs the instantaneous rate matrices for use in the forward backward algorith
buildirm <- function(X, b, m, a, pop=FALSE){
    Xobs <- rbind(c(0,0,sum(X[X[,1]==0,3])), X[X[,1]!=0,]); times <- unique(Xobs[,1])
    numsusc <- popsize - cumsum(Xobs[,3]==1) - Xobs[1,3]; numinf <- cumsum(Xobs[,3]) #cumulative counts of the numbers of infecteds and recoverds at event times
    
    #   check if the infection has died off. ensure that the number of infecteds is zero forever after.
    if(any(numinf==0) & a==0){
        lastinf <- which(numinf==0)[1]
        numinf[lastinf:length(numinf)] <- 0
    }
    
    irm <- array(0, dim = c(3,3,dim(Xobs)[1]-1))
    if(pop==FALSE){
        irm[1,2,] <- b*numinf[1:(length(numinf)-1)] + a; irm[1,1,] <- -irm[1,2,]
        irm[2,3,] <- m; irm[2,2,] <- -m
        
    } else if(pop==TRUE){
        irm[1,2,] <- (b*numinf[1:(length(numinf)-1)] + a)*numsusc[1:(length(numsusc)-1)]; irm[1,1,] <- -irm[1,2,]
        irm[2,3,] <- m*numinf[1:(length(numinf)-1)]; irm[2,2,] <- irm[2,3,]
    }
    
    return(irm)
}

# buildtpm constructs a transition probability matrix using the matrix exponential
buildtpm <- function(Q, t0, t1){
    Qeig <- eigen(Q)
    Qeig$vectors%*%diag(exp(Qeig$values*(t1-t0)))%*%solve(Qeig$vectors)
}

# obstpm computes the product of transition probability matrices between observation times. takes an array of instantaneous rate matrices,
# an augmented data matrix, two observation times as inputs 
obstpm <- function(Xother, irm, t0, t1){
    timeseq <- c(t0,Filter(function(X) X>t0 & X<t1, Xother[,1]), t1)
    tpm <- diag(1,3); 
    if(0 %in% unique(Xother[,1])){
        times <- unique(Xother[,1])
    } else times <- c(0, unique(Xother[,1]))
    for (k in 1:(length(timeseq)-1)) {
        ind1 <- sum(times<=timeseq[k]); 
        if(ind1 <= dim(irm)[3]){
            irm1 <- irm[,,ind1]
        } else if(ind1 > dim(irm)[3]){
            irm1 <- diag(0,3); irm1[2,3] <- irm[2,3,dim(irm)[3]]; irm1[2,2] <- -irm1[2,3]
        }
        tpm <- tpm %*% buildtpm(irm1,timeseq[k],timeseq[k+1])
    } 
    return(tpm)
}

# numinfected takes a matrix formatted like X or Xother and computes the number of infected individuals at each event time
numinfected <- function(X){
    c(sum(X[,3][X[,1]==0]),cumsum(X[,3][X[,1]!=0]) + sum(X[,3][X[,1]==0]))
}

# initialize_X initializes the matrix of trajectories in the population by dispersing infecteds throughout the observation period
# dat is a matrix with observation times and counts of infecteds, mu is the recovery parameter, p is the sampling probability, amplify
# controls the factor by which the observed number of infecteds is multiplied to get the number of initialized trajectories (up to the population size)
# X is matrix with event times, subject id and event codes. 
# Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
initializeX <- function(W, mu, p, amplify, tmax){
    whichsick <- sample(1:popsize,min(popsize,floor(sum(W[,2])/(amplify*p))))
    
    X <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
    
    for(k in 1:dim(W)[1]){
        if(W[k,2]!=0) {
            if(length(whichsick) < (W[k,2] - sum(X[which(X[,1]<W[k,1]),3]))){
                numneeded <-W[k,2] - sum(X[which(X[,1]<=W[k,1]),3]) - length(whichsick); selected <- rep(0,numneeded)
                ids <- 1:popsize; ids <- ids[!ids %in% whichsick]
                for(r in 1:numneeded){
                    noexcess <- which(W[,3]<=W[,2]); noexcess.times <- W[noexcess,1]
                    keepchoosing <- TRUE
                    while(keepchoosing==TRUE){
                        proposed <- sample(ids, 1)
                        proposed.path <- getpath(X, proposed)
                        if(!any(proposed.path[1] < noexcess.times & proposed.path[2]>noexcess.times)){
                            selected[r] <- proposed
                            ids <- ids[!ids==proposed]
                            X <- updateX(X, c(0,0), proposed)
                            W <- updateW(W, X)
                            keepchoosing <- FALSE
                        } else {
                            ids <- ids[!ids==proposed]
                            keepchoosing<-TRUE
                        }
                    }
                }
                whichsick <- c(whichsick,selected)
            }
            if(sum(X[X[,1]<=W[k,1],3]) <= W[k,2]){
                nowsick <- sample(whichsick,max(0,W[k,2] - sum(X[X[,1]<=W[k,1],3]))); whichsick <- whichsick[-which(whichsick %in% nowsick)]
                
                for(j in 1:length(nowsick)){
                    ind <- which(X[,2] == nowsick[j])
                    if(k == 1){    
                        X[ind,1][1] <- 0; X[ind,3][1]<-1
                        tau <- rexp(1,rate=mu)
                        X[ind,1][2]<-ifelse(tau>tmax,0,tau)
                        X[ind,3][2] <- ifelse(tau>tmax,0,-1)
                    } else{
                        tau <- rexp(1,rate=mu); eventtime <- runif(1,0, tau)
                        
                        X[ind,1][1] <- max(0,W[k,1] - eventtime); X[ind,3][1] <- 1
                        
                        X[ind,1][2]<-ifelse((W[k,1] + (tau-eventtime - min(0,W[k,1]-eventtime)))>tmax,0,W[k,1] + (tau-eventtime- min(0,W[k,1]-eventtime)))
                        X[ind,3][2] <- ifelse(W[k,1] + (tau-eventtime - min(0,W[k,1]-eventtime))>tmax,0,-1)
                    }
                }   
            } else next          
        }   
    }
    
    X<-X[order(X[,1]),]
    
    return(X)
}

# updateW updates an observation matrix. The function takes as arguments an observation matrix with trajectories, X, and a matrix, W, 
# with columns containing observation times, number of infecteds observed, and number of infecteds based on augmented trajectories. 
updateW <- function(W, X){
    for(k in 1:dim(W)[1]){
        W[k,3] <- sum(X[which(X[,1]<=W[k,1]),3])
    }
    return(W)
}

# updateX updates the matix of observations X whose columns are event times, subject id, and event codes. Inputs are the original X matrix,
# a vector of times, and the subject id. The function returns an update X matrix, ordered by event time. 
updateX <- function(X, Xt.path, j){
    if(all(Xt.path==0)){
        X[X[,2]==j,c(1,3)] <- 0
        
    } else if(all(Xt.path==Inf)){
        X[X[,2]==j,c(1,3)] <- 0
        
    } else if(Xt.path[1] != Inf & Xt.path[2] ==Inf){
        X[X[,2]==j,][1,1] <- Xt.path[1]; X[X[,2]==j,][1,3] <- 1
        X[X[,2]==j,][2,1] <- 0 ; X[X[,2]==j,][2,3] <- 0
        
    } else if(all(Xt.path!=Inf)){
        X[X[,2]==j,][1,1] <- Xt.path[1]; X[X[,2]==j,][1,3] <- 1
        X[X[,2]==j,][2,1] <- Xt.path[2] ; X[X[,2]==j,][2,3] <- -1
    }
    
    X <- X[order(X[,1]),]
    
    return(X)
}

# getpath extracts the current path for subject j
getpath <- function(X.cur, j) {
    
    if(all(X.cur[X.cur[,2]==j,3]==0)){ # first the case when no infection is observed
        cur.path <- c(0,0)
        
    } else if(X.cur[X.cur[,2]==j,3][1]==1 & X.cur[X.cur[,2]==j,3][2]==0){ # then the case when an infection but no recovery is observed
        cur.path <- c(X.cur[X.cur[,2]==j,1][1],Inf)
        
    } else if(X.cur[X.cur[,2]==j,3][1] == 1 & X.cur[X.cur[,2]==j,3][2] == -1) { # finally when both an infection and a recovery are observed
        cur.path <- c(X.cur[X.cur[,2]==j,1])
    }
}

# Other functions ---------------------------------------------------------

calc_loglike <- function(X, W, b, m, a=0, p){
#     popsize = length(unique(X[,2])); times <- unique(X[,1]); timediffs <- times[2:length(times)] - times[1:(length(times)-1)]
    Xobs <- rbind(c(0,0,sum(X[X[,1]==0,3])), X[X[,1]!=0,]); irm <- buildirm(X, b, m, a, pop=TRUE)
    
    events <- ifelse(Xobs[Xobs[,1]!=0,3]==1, 1,2); rates <- ifelse(events==1,irm[1,2,], irm[2,3,])
    
    dbinom(sum(W[,2]), sum(W[,3]), prob=p, log=TRUE) + sum(log(rates)) - sum(rates*(Xobs[2:dim(Xobs)[1],1] - Xobs[1:(dim(Xobs)[1]-1),1]))
} 

# pop_prob and path_prob calculate the log-probabilities of the population trajectory and the subject trajectory for use in the M-H ratio

pop_prob <- function(X, irm, initdist){
    Xobs <- rbind(c(0,0,sum(X[X[,1]==0,3])), X[X[,1]!=0,])
    events <- ifelse(Xobs[Xobs[,1]!=0,3]==1, 1,2); rates <- ifelse(events==1,irm[1,2,], irm[2,3,])
    
    sum(log(rates)) - sum(rates*(Xobs[2:dim(Xobs)[1],1] - Xobs[1:(dim(Xobs)[1]-1),1]))
}

path_prob <- function(path, Xother, irm.other, initdist, tmax){
    if(all(path==0)){ # no infection observed
        times <- unique(Xother[,1]); timediffs <- times[2:length(times)] - times[1:(length(times)-1)]
        path.prob <- log(initdist[1])-sum(irm.other[1,2,] * timediffs)
            
    } else if(path[2] > 0 & path[2] < Inf) { # subject is infected and recovery is observed
        
        if(path[1]==0){ # subject is initially infected 
            
            path.prob <- log(initdist[2]) + log(irm.other[2,3,1]) - irm.other[2,3,1]*path[2]
            
        } else if(path[1]!=0){ # subject is not initially infected
            
            times <- unique(Xother[,1]); timediffs <- times[2:length(times)] - times[1:(length(times)-1)]
            
            if(times[2] > path[1]){ # no changes in number of infecteds before subject's infection
                
                path.prob <- log(initdist[1]) + log(irm.other[1,2,1]) - irm.other[1,2,1]*path[1] + log(irm.other[2,3,1]) - irm.other[2,3,1]*(path[2] - path[1])
                
            } else if(times[2] < path[1]){ # if there is at least one change in the number of infecteds before the subject's infection
                
                path.prob <- log(initdist[1]) - sum(irm.other[1,2,which(times < path[1])] * timediffs[which(times < path[1])]) +
                    log(irm.other[1, 2, Position(function(y) y>= path[1], times) - 1]) - irm.other[1,2,Position(function(y) y>= path[1], times) - 1] * (path[1] - times[Position(function(y) y>path[1], times) - 1]) +
                    log(irm.other[2,3,1]) - irm.other[2,3,1]*(path[2] - path[1])
                
            }
        }
        
    } else if(path[2] == Inf) { # subject is infected and no recovery is observed
        
        if(path[1]==0){ # subject is initially infected 
            
            path.prob <- log(initdist[2]) - irm.other[2,3,1]*tmax
            
        } else if(path[1]!=0){ # subject is not initially infected
            
            times <- unique(Xother[,1]); timediffs <- times[2:length(times)] - times[1:(length(times)-1)]
            
            if(times[2] > path[1]){ # no changes in number of infecteds before subject's infection
                
                path.prob <- log(initdist[1]) + log(irm.other[1,2,1]) - irm.other[1,2,1]*path[1] - irm.other[2,3,1]*(tmax - path[1])
                
            } else if(times[2] < path[1]){ # if there is at least one change in the number of infecteds before the subject's infection
                
                path.prob <- log(initdist[1]) - sum(irm.other[1,2,which(times < path[1])] * timediffs[which(times < path[1])]) +
                    log(irm.other[1, 2, Position(function(y) y>= path[1], times) - 1]) - irm.other[1,2,Position(function(y) y>= path[1], times) - 1] * (path[1] - times[Position(function(y) y>path[1], times) - 1]) -
                    irm.other[2,3,1]*(tmax - path[1])
                
            }
        }
    } 
    
    return(path.prob)
}
# checkpossible checks whether removing a subject to resimulate the epidemic would cause an impossible population trajectory. If so, M-H automatically rejects.

checkpossible <- function(X, a=0, W=NULL){
    popseq <-c(sum(X[X[,1]==0,3]), sum(X[X[,1]==0,3]) + cumsum(X[X[,1]!=0,3]))
    if(any(popseq==0) & a==0){
        if(any(which(popseq==0)!=length(popseq)) == TRUE){
            is.possible <-FALSE
        } else if(any(which(popseq==0)!=length(popseq)) == FALSE){
            is.possible <- TRUE
        }
    } else {
        is.possible <- TRUE
    }
    
    if(!is.null(W)){
        W <- updateW(W,X)
        if(any(W[,3]<W[,2])){
            is.possible <- FALSE
        }
    }
    return(is.possible)
}


# augSIR wrapper ----------------------------------------------------------

# augSIR is the wrapper to estimate SIR epidemic parameters via Bayesian data augmentation

augSIR <- function(dat, sim.settings, priors, inits) {
    # initialize simulation settings
    popsize <- sim.settings$popsize # size of the population; 
    tmax <- sim.settings$tmax # maximum time of observation
    niter <- sim.settings$niter # number of iterations in the sampler
    initdist <- sim.settings$initdist # initial distribution for individual infection status
    
    # vectors for parameters
    Beta <- vector(length=niter); Beta[1] <- inits$beta.init
    Mu <- vector(length = niter); Mu[1] <- inits$mu.init
    Alpha <- vector(length = niter); Alpha[1] <- inits$alpha$init 
    probs <- vector(length = niter); probs[1] <- inits$probs.init
    
    # vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
    beta.prior <- priors$beta.prior
    mu.prior <- priors$mu.prior
    # alpha.prior <- c(6, 12000)
    p.prior <- priors$p.prior
    
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
    for(k in 1:(niter-1)){
        # Update trajectories
        
        subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
        irm.cur <- buildirm(X.cur, b = Beta[k], m = Mu[k], a = Alpha[k]) 
        
        for(j in 1:length(subjects)){
            Xother <- X.cur[X.cur[,2]!=subjects[j],]; 
            if(checkpossible(Xother) == TRUE){
                
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
                }
                
            } else next
        }
        
        # Update observation matrix
        W <- updateW(W, X)
        
        # draw new binomial sampling probability parameter
        #   probs[k+1] <- rbeta(1,p.prior[1] + sum(W[,2]), p.prior[2] + sum(W[,3]-W[,2]))
        probs[k+1] <- 1
        
        # draw new rate parameters 
        Beta[k+1] <- rgamma(1, shape = (beta.prior[1] + sum(X[,3]==1)), 
                             rate = beta.prior[2] + sum((popsize - cumsum(X[,3]==1))*(cumsum(X[,3]==1) - cumsum(X[,3]==-1))*c(0,X[2:dim(X)[1],1] - X[1:(dim(X)[1]-1),1])*(X[,3]==1)))
          
        Mu[k+1] <- rgamma(1, shape = mu.prior[1] + sum(X[,3]==-1),
                          rate = mu.prior[2] + sum((cumsum(X[,3]==1) - cumsum(X[,3]==-1))*c(0,X[2:dim(X)[1],1] - X[1:(dim(X)[1]-1)])*(X[,3]==-1)))
        
        #   Alpha[k+1] <- rgamma(1, shape = (alpha.prior[1] + sum(X[,3]==1)),
        #                        rate = alpha.prior[2] + sum((popsize - cumsum(X[,3]==1))*c(0,X[2:dim(X)[1],1] - X[1:(dim(X)[1]-1),1])*(X[,3]==1)))
        Alpha[k+1] <- 0
        
        loglik[k] <- calc_loglike(X, W, Beta[k+1], Mu[k+1], Alpha[k+1], probs[k+1])  
        
        #   trajectories <- data.frame(X); epidemic <- data.frame(W)
        #   trajectories$id <- factor(trajectories$id, levels = unique(as.factor(trajectories$id)))
        #   trajectories$event <- factor(trajectories$event)
        #   trajects <- ggplot(subset(trajectories,time!=0),aes(y=id,x=time,group=id,colour=as.factor(event))) + geom_point()+ geom_line() + geom_vline(xintercept = data.frame(W)$time,alpha=0.2)+
        #     scale_colour_discrete(name="Event", labels=c("Infection Cleared", "Infection Acquired"))
        #   
        #   curve <- ggplot(epidemic,aes(x=time,y=augmented))+geom_line() + theme_bw()
        #   
        #   print(grid.arrange(trajects, curve, ncol=2))
        
    }
    
    results <- list(Beta = Beta, Mu = Mu, loglik = loglik)
    return(results)
    
}
