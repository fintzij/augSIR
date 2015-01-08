
# Functions for drawing augmented trajectories ----------------------------

# drawpath draws a path for an endpoint conditioned markov chain
# Xt is the subject level observation matrix with first column a vector of observation times, and second column a vector of infection
# statuses. Xother is the matrix with times of carriage acquisition and clearance for all other subjects. 
# irm is an array of instantaneous rate matrices. tmax is the maximum time of observation.

drawpath <- function(Xt, Xother, irm, tmax){
    path <- c(0,0)
    changetimes <- which(diff(Xt[,2]) != 0)
    
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
        
    } else if(length(changetimes)==1){ ### all transitions occur in one observation period
        if(Xt[changetimes+1,2] - Xt[changetimes,2] ==2){ ### Subject transitions from healthy to recovered within one observation period
            t0 <- Xt[changetimes,1]; t1 <- Xt[changetimes+1,1]
            path[1] <- drawtime(Xother, irm, t0, t1, 1)
            path[2] <- ifelse(path[1]==Inf,Inf,drawtime(Xother, irm, path[1], t1, 2))
            
        } else if(Xt[changetimes,2]==1 & (Xt[changetimes+1,2] - Xt[changetimes,2] == 1)){ # healthy to infected within one observation period
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
    times <- unique(c(0,Xother[,1]))
    
    if(t1 == Inf){
        if(currentstate == 2){
            eventtime <- t0 - log(1-runif(1))/irm[2,3,1]
            
        } else if(currentstate == 1){
            timeseq <- c(t0, times[times>t0])
            if(length(timeseq) == 1) {
                eventtime <- Inf
                
            } else if (length(timeseq) > 1){
                init.infec <- sum(Xother[,3]==1 & Xother[,1]==0)
                numinf <- c(init.infec, init.infec + cumsum(Xother[Xother[,1]!=0,3]))
                
                for(k in 1:(length(timeseq)-1)){
                    ind <- which(irm[4,4,] == numinf[times <= timeseq[k]][sum(times<=timeseq[k])])
                    
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
            if (t0 > max(times)){
                eventtime <- Inf
                
            } else if(t0 <= max(times)){
                if(t1 < max(Xother[,1])){
                    timeseq <- c(t0, times[times > t0 & times < t1], t1)
                    
                } else if(t1 >= max(times)){
                    timeseq <- c(t0, times[times > t0])
                }
                
                init.infec <- sum(Xother[,3]==1 & Xother[,1]==0)
                numinf <- c(init.infec, init.infec + cumsum(Xother[Xother[,1]!=0,3]))
                
                if(length(timeseq)!=2) {
                    
                    intervalprobs <- rep(0,length(timeseq)-1)
                    timediffs <- diff(timeseq)
                    indstart <- sum(times<=timeseq[1]); indend <- sum(times<=timeseq[length(timeseq)-1])
                    totalprob <- 1 - exp(-sum(irm[1, 2,match(numinf[indstart:indend], irm[4,4,])]*timediffs))
                    
                    for(k in 1:length(intervalprobs)){
                        ind <- sum(times<=timeseq[k])
                        intervalprobs[k] <- (1-exp(-sum(irm[1, 2,match(numinf[indstart:ind], irm[4,4,])]*timediffs[1:k])))/totalprob
                    }
                    
                    interval <- which(intervalprobs > runif(1))[1]
                    ind <- sum(times <= timeseq[interval])
                    
                    eventtime<-timeseq[interval]-log(1-runif(1)*(1-exp(-irm[1, 2, match(numinf[ind], irm[4,4,])]*timediffs[interval])))/irm[1, 2, match(numinf[ind], irm[4,4,])]
                    
                } else if(length(timeseq)==2){
                    timediffs <- t1 - t0; ind <- sum(times <= t0)
                    
                    indstart <- sum(times<=timeseq[1]); indend <- sum(times<=timeseq[length(timeseq)-1])
                    
                    eventtime <- t0 - log(1 - runif(1)*(1 - exp(-irm[1, 2, match(numinf[ind], irm[4,4,])]*(t1 - t0))))/irm[1, 2, match(numinf[ind], irm[4,4,])]
                }
                
            }
        }
    }
    
    return(eventtime)
}


# drawXt draws the status of a subject at observation times t1,...,tL conditional on other individuals
# Xother is a matrix with the times of carriage aquisition and clearance for all other subjects. 
# p, is the sampling probability, and b, and m are the epidemiologic parameters. initdist is the initial
# distribution for infection status

drawXt <- function(Xother, irm, irm.eig, W, p, b, m, a, initdist){
    
    ### ask about this, does it make more sense to use a different initial dist for this case, or should the trajectory not be 
    ### simulated outright.
#     if (sum(Xother[,3][Xother[,1]==0])<W[1,2]) {
#         initdist <- c(0,1,0)
#     } 
    
    Xt <- cbind(W[,1], 0)
    
    Xt.fb <- fwd(Xother, W, irm, irm.eig, initdist, p)
    
    Xt<-bwd(Xt.fb,Xt)
        
    return(Xt)
}

# Forward-backward recursion functions ------------------------------------

#fwd produces P_tj from P_t,j-1
fwd <- function(Xother, W, irm, irm.eig, initdist, p){
    
    Xt.fb <- array(0,dim=c(3,3,dim(W)[1]-1))
    obstimes <- W[,1]
    
    eventtimes <- unique(c(0, Xother[,1]))
    
    init.infec <- sum(Xother[,3]==1 & Xother[,1]==0); numinf <- c(init.infec, init.infec + cumsum(Xother[Xother[,1]!=0,3]))
    
    tpm <- obstpm(numinf = numinf, eventtimes = eventtimes, irm = irm, irm.eig = irm.eig, t0 = obstimes[1], t1 = obstimes[2])
        
    numinf.aug <- numinf[eventtimes<=obstimes[2]][sum(eventtimes<=obstimes[2])]; numinf.obs <- W[2,2]
    
    if(p!=1){
        emit <- dbinom(numinf.obs, numinf.aug +c(0,1,0),p) 
        
    } else if(p==1){
        emit <- dbinom(numinf.obs, min(numinf.aug, numinf.obs) +c(0,1,0),p)
        
    }
    
    Xt.fb[,,1]<-normalize(outer(initdist,emit,FUN="*")*tpm)

    for(k in 2:dim(Xt.fb)[3]){
        tpm <- obstpm(numinf = numinf, eventtimes = eventtimes, irm = irm, irm.eig = irm.eig, t0 = obstimes[k], t1 = obstimes[k+1])
        
        distr <- colSums(Xt.fb[,,k-1])
        
        numinf.aug <- numinf[eventtimes <= obstimes[k+1]][sum(eventtimes <= obstimes[k+1])]; numinf.obs <- W[W[,1]==obstimes[k+1],2]
        
        if(p==1){
            emit <- dbinom(numinf.obs, min(numinf.aug, numinf.obs) +c(0,1,0),p)
            
        } else if(p!=1){
            #### check this with Vladimir. the emission probability in the recovered state may depend on previous and future infection status
            emit <- dbinom(numinf.obs, numinf.aug +c(0,1,0),p)
            
        }
        
        Xt.fb[,,k]<-normalize(outer(distr,emit,FUN="*")*tpm)
    }
    
    return(Xt.fb)
}

# bwd implements the stochastic backward recursion, drawing the infection status at time t from the 
# distribution proportional to the column corresponding to the infection status at time t+1

bwd <- function(Xt.fb,Xt){
    states <- rep(0,dim(Xt)[1]); draws <- runif(dim(Xt)[1])
    
    initdist<-cumsum(colSums(Xt.fb[,,dim(Xt.fb)[3]]))

    states[length(states)] <- which(initdist >= draws[length(states)])[1] 
    
    for(k in (dim(Xt)[1] - 1):1){
        dist<-cumsum(normalize(Xt.fb[,states[k+1],k]))
        states[k] <- which(dist >= draws[k])[1] 
    }
    Xt[,2] <- states
    return(Xt)
}


# Functions for constructing various matrices -----------------------------

# Normalize is what it sounds like.
normalize <- function(X){
    X/sum(X)
}

buildirm <- function(X, b, m, a=0, popsize, pop){
    init.infec <- sum(X[,3]==1 & X[,1]==0)
    
    if(pop == FALSE){
        numinf <- seq(0, init.infec + max(cumsum(X[X[,1]!=0,3]))+1)
        irm <- array(0, dim = c(4,4,length(numinf)))
        irm[1,2,] <- b * numinf + a; irm[1,1,] <- -irm[1,2,]
        irm[2,3,] <- m; irm[2,2,] <- -m
        irm[4,4,] <- numinf
        
    } else if(pop == TRUE){
        numinf <- c(init.infec, init.infec + cumsum(X[X[,1]!=0,3]))
        numsusc <-  popsize - cumsum(c(init.infec,X[X[,1]!=0,3]==1)) 
        
        irm <- array(0, dim = c(3,3, length(numinf) - 1))
        irm[1,2,] <- (b*numinf[1:(length(numinf)-1)] + a)*numsusc[1:(length(numsusc)-1)]; irm[1,1,] <- -irm[1,2,]
        irm[2,3,] <- m*numinf[1:(length(numinf)-1)]; irm[2,2,] <- -irm[2,3,]
    
    }
    
    return(irm)
}

# irm_decomp takes an array of path irms and returns an array of eigenvalues, eigenvectors, and the inverse of the eigenvectors
# irm_decomp returns an array, each element of which is another array containing a diagonal matrix of eigenvalues, a matrix of eigenvectors, and 
# the inverse of the matrix of eigenvectors
irm_decomp <- function(pathirm.cur){
    decomp <- apply(pathirm.cur[c(1:3),c(1:3),], 3, eigen, symmetric = FALSE)
    
    irm.decomp <- array(0,dim = c(3,3,3,length(decomp)))
    
    irm.decomp[,, 1,] <- array(apply(sapply(decomp, '[[', 1), 2, diag), dim = c(3,3,length(decomp)))
    irm.decomp[,, 2,] <- array(sapply(decomp, '[[', 2), dim = c(3,3,length(decomp)))
    irm.decomp[,, 3,] <- array(apply(irm.decomp[,,2,], 3, solve), dim = c(3,3,length(decomp)))
    
    return(irm.decomp)
}

update_irm <- function(irm, new.numinf, b, m , a , popsize ){
    newmat <- diag(0,4)
    newmat[1,2] <- (b * new.numinf + a); newmat[1,1] <- -newmat[1,2]
    newmat[2,3] <- m; newmat[2,2] <- -m
    newmat[4,4] <- new.numinf
    
    irm.new <- array(c(irm,newmat), dim = c(4,4,dim(irm)[3]+1))
    
    return(irm.new)
}

update_eigen <- function(patheigen, pathirm){
    decomp <- eigen(pathirm[c(1:3),c(1:3),dim(pathirm)[3]], symmetric = FALSE)
    
    eigen.new <- array(c(patheigen, 
                             array(c(diag(decomp$values), decomp$vectors, solve(decomp$vectors)), dim = c(3,3,3))), 
                           dim = c(3,3,3,dim(patheigen)[4] + 1))        
    
    return(eigen.new)    
}

# buildtpm constructs a transition probability matrix using the matrix exponential
buildtpm <- function(values, vectors, inv.vecs, t0, t1){
    vectors%*%diag(exp(diag(values)*(t1-t0)))%*%inv.vecs
}

# obstpm computes the product of transition probability matrices between observation times. 

obstpm <- function(numinf, eventtimes, irm, irm.eig, t0, t1){
    timeseq <- c(t0,eventtimes[eventtimes>t0 & eventtimes <t1], t1)
    
    tpm <- diag(1,3) 
    
    for (j in 1:(length(timeseq)-1)) {
        ind <- match(numinf[eventtimes <= timeseq[j]][sum(eventtimes <= timeseq[j])], irm[4,4,])
            
        tpm <- tpm %*% buildtpm(values = irm.eig[,,1,ind],
                                vectors = irm.eig[,,2,ind],
                                inv.vecs = irm.eig[,,3,ind],
                                timeseq[j],timeseq[j+1])
        
    } 
    return(tpm)
}

# initialize_X initializes the matrix of trajectories in the population by dispersing infecteds throughout the observation period
# dat is a matrix with observation times and counts of infecteds, mu is the recovery parameter, p is the sampling probability, amplify
# controls the factor by which the observed number of infecteds is multiplied to get the number of initialized trajectories (up to the population size)
# X is matrix with event times, subject id and event codes. 
# Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
initializeX <- function(W, mu, p, amplify, tmax, popsize){
    whichsick <- sample(1:popsize,min(popsize,floor(sum(W[,2])/p)*amplify))
    
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
    if(all(Xt.path == 0) | all(Xt.path == Inf)){
        X[X[,2]==j,c(1,3)] <- 0
        
    }  else if(Xt.path[1] != Inf & Xt.path[2] ==Inf){
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
getpath <- function(X, j) {
    
    subj <- X[X[,2]==j,]
    
    if(all(subj[,3]==0)){
        path <- c(0,0)
        
    } else if(1 %in% subj[,3] & !(-1 %in% subj[,3])){
        path <- c(subj[subj[,3]==1,1], Inf)
        
    } else if(1 %in% subj[,3] & (-1) %in% subj[,3]){
        path <- subj[,1]
        
    }
    
    return(path)
}

# Other functions ---------------------------------------------------------

calc_loglike <- function(X, W, irm, b, m, a=0, p, initdist, popsize){
    initinfec <- sum(X[,3][X[,1]==0])
    Xobs <- rbind(c(0,0,initinfec), X[X[,1]!=0,]); timediffs <- diff(unique(c(0,X[,1])))
    events <- ifelse(Xobs[Xobs[,1]!=0,3]==1, 1,2); rates <- ifelse(events==1,irm[1,2,], irm[2,3,])
    
    dbinom(sum(W[,2]), sum(W[,3]), prob=p, log=TRUE) + dmultinom(c(popsize - initinfec, initinfec, 0), prob = initdist, log=TRUE) + sum(log(rates)) - sum(rates*timediffs)
} 

# pop_prob and path_prob calculate the log-probabilities of the population trajectory and the subject trajectory for use in the M-H ratio

pop_prob <- function(X, irm, initdist, popsize){

        Xobs <- rbind(c(0,0,sum(X[X[,1]==0,3])), X[X[,1]!=0,])
        init.infec <- sum(X[,3]==1 & X[,1]==0)
        events <- ifelse(Xobs[Xobs[,1]!=0,3]==1, 1,2); rates <- ifelse(events==1,irm[1,2,], irm[2,3,])
        
        # note that it is possible to have an impossible trajectory proposed if the infection is depleted between time 0 and the first observation time, but there are still enough
        # augmented infections by the first observation time (this is the first time for which an emission probability is computed). Therefore, rates are set to max(0, rates)
        
        rates <- pmax(0,rates)
        
        dmultinom(c(popsize - init.infec, init.infec, 0), prob = initdist, log=TRUE) + sum(log(rates)) - sum(rates*diff(Xobs[,1]))

}


path_prob <- function(path, Xother, pathirm, initdist, tmax){
    times <- c(unique(c(0,Xother[,1])), tmax); timediffs <- diff(times)
    init.infec <- sum(Xother[,3]==1 & Xother[,1]==0); numinf <- c(init.infec, init.infec + cumsum(Xother[Xother[,1]!=0,3]))
    
    if(all(path==0)){ # no infection observed
                
        path.prob <- log(initdist[1])-sum(pathirm[1, 2, match(numinf,pathirm[4,4,])] * timediffs)
            
    } else if(path[2] > 0 & path[2] < Inf) { # subject is infected and recovery is observed
        
        if(path[1]==0){ # subject is initially infected 
            
            path.prob <- log(initdist[2]) + log(pathirm[2,3,1]) - pathirm[2,3,1]*path[2]
            
        } else if(path[1]!=0){ # subject is not initially infected
                        
            if(times[2] > path[1]){ # no changes in number of other infecteds before subject's infection
                
                path.prob <- log(initdist[1]) + log(pathirm[1,2, match(init.infec,pathirm[4,4,])]) - pathirm[1,2,match(init.infec, pathirm[4,4,])]*path[1] + 
                    log(pathirm[2,3,1]) - pathirm[2,3,1]*(path[2] - path[1])
                
            } else if(times[2] < path[1]){ # if there is at least one change in the number of infecteds before the subject's infection
                ind1 <- which(times < path[1])[sum(times < path[1])]
                
                path.prob <- log(initdist[1]) - sum(pathirm[1,2, match(numinf[1:(ind1-1)], pathirm[4,4,])] * timediffs[1:(ind1-1)]) + 
                    log(pathirm[1,2,match(numinf[ind1], pathirm[4,4,])]) - pathirm[1,2, match(numinf[ind1], pathirm[4,4,])]*(path[1] - times[ind1]) + 
                    log(pathirm[2,3,1]) - pathirm[2,3,1]*(path[2] - path[1]) 

            }
        }
        
    } else if(path[2] == Inf) { # subject is infected and no recovery is observed
        
        if(path[1]==0){ # subject is initially infected 
            
            path.prob <- log(initdist[2]) - pathirm[2,3,1]*tmax
            
        } else if(path[1]!=0){ # subject is not initially infected
                        
            if(times[2] > path[1]){ # no changes in number of infecteds before subject's infection
                
                path.prob <- log(initdist[1]) + log(pathirm[1,2, match(init.infec,pathirm[4,4,])]) - pathirm[1,2,match(init.infec,pathirm[4,4,])]*path[1] - pathirm[2,3,1]*(tmax - path[1])
                
            } else if(times[2] < path[1]){ # if there is at least one change in the number of infecteds before the subject's infection
                ind1 <- which(times < path[1])[sum(times < path[1])]
                
                path.prob <- log(initdist[1]) - sum(pathirm[1,2, match(numinf[1:(ind1-1)], pathirm[4,4,])] * timediffs[1:(ind1-1)]) + 
                    log(pathirm[1,2,match(numinf[ind1], pathirm[4,4,])]) - pathirm[1,2, match(numinf[ind1], pathirm[4,4,])]*(path[1] - times[ind1])  -
                    pathirm[2,3,1]*(tmax - path[1])
                
            }
        }
    } 
    
    return(path.prob)
}

# Functions to update parameters
update_beta <- function(X.cur, beta.prior, popsize){
    X <- X.cur[X.cur[,1]!=0,]
    
    init.infec <- sum(X.cur[,1] == 0 & X.cur[,3]==1); infections <- X[,3] == 1; 
    numsick <- c(init.infec, init.infec + cumsum(X[,3])); numsusc <- popsize - cumsum(c(init.infec,X[,3]==1)) 
    timediffs <- diff(c(0,X[,1]))
    
    rgamma(1, shape = (beta.prior[1] + sum(infections)), 
           rate = beta.prior[2] + sum(numsick[1:(length(numsick)-1)] * numsusc[1:(length(numsusc) - 1)] * timediffs * infections))
}

update_mu <- function(X.cur, mu.prior, popsize){
    X <- X.cur[X.cur[,1]!=0,]
    
    init.infec <- sum(X.cur[,1] == 0 & X.cur[,3]==1); recoveries <- X[,3] == -1
    numsick <- c(init.infec, init.infec + cumsum(X[,3]))
    timediffs <- diff(c(0,X[,1]), lag = 1)
    
    rgamma(1, shape = (mu.prior[1] + sum(recoveries)), 
           rate = mu.prior[2] + sum(numsick[1:(length(numsick)-1)] * timediffs * recoveries))
    
}

update_alpha <- function(X.cur, alpha.prior, popsize){
    rgamma(1, shape = (alpha.prior[1] + sum(X.cur[,3]==1)),
           rate = alpha.prior[2] + sum((popsize - cumsum(X.cur[,3]==1))*c(0,X.cur[2:dim(X.cur)[1],1] - X.cur[1:(dim(X.cur)[1]-1),1])*(X.cur[,3]==1)))
}

# accept_prob calculates the acceptance ratio for the M-H algorithm
accept_prob <- function(pop_prob.new, pop_prob.cur, path_prob.cur, path_prob.new){
    if(pop_prob.new == -Inf & path_prob.new == -Inf){
        -Inf
    } else{
        pop_prob.new - pop_prob.cur + path_prob.cur - path_prob.new
    }
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

# find.pprior finds the parameters for the beta prior for the binomial sampling probability
find.pprior <- function(mu, sigmasq, inits){
    require(rootSolve)
    
    prior <- function(x) {
        c(F1 = x[1]/(x[1] + x[2]) - mu,
          F2 = x[1]*x[2] / (x[1] + x[2])^2 / (x[1]+ x[2] + 1) - sigmasq)
    }
    
    ss <- multiroot(f = prior, start = c(inits[1],inits[2]))
    return(ss)
}

find.rateprior <- function(mu, sigmasq, inits){
    require(rootSolve)
    
    prior <- function(x){
        c(F1 = x[1]/x[2] - mu,
          F2 = x[1]/x[2]^2 - sigmasq)
    }
    
    ss <- multiroot(f=prior, start = c(inits[1], inits[2]))
    return(ss)
}
# augSIR wrapper ----------------------------------------------------------

# augSIR is the wrapper to estimate SIR epidemic parameters via Bayesian data augmentation

augSIR <- function(dat, sim.settings, priors, inits) {
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
    X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
    
    X.cur <- initializeX(W = W.cur, mu = Mu[1], p=probs[1], amplify = amplify, tmax=20, popsize = popsize)
    
    # update observation matrix
    W.cur <- updateW(W.cur,X.cur)
    
    if(!checkpossible(X=X.cur, W=W.cur)) {
        while(!checkpossible(X=X.cur,W=W.cur)){
            X.cur <- initializeX(W = W.cur, mu = Mu[1], p=probs[1], amplify = amplify, tmax=20, popsize = popsize)
            W.cur <- updateW(W.cur,X.cur)
        }
    }
    
    popirm.cur <- buildirm(X.cur, b = Beta[1], m = Mu[1], a = Alpha[1], popsize = popsize, pop = TRUE)
    pop_prob.cur <- pop_prob(X.cur, irm = popirm.cur, initdist = initdist, popsize = popsize)
    
    # M-H sampler
    for(k in 2:niter){
        # Update trajectories
        print(k)
        subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
        
        pathirm.cur <- buildirm(X.cur, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = FALSE)
        patheigen.cur <- irm_decomp(pathirm.cur)
        
        for(j in 1:length(subjects)){
            Xother <- X.cur[X.cur[,2]!=subjects[j],]
            
            path.cur <- getpath(X.cur, subjects[j])
            
            W.other <-updateW(W.cur, Xother)
            Xt <- drawXt(Xother = Xother, irm = pathirm.cur, irm.eig = patheigen.cur, W=W.other, p=probs[k-1], b=Beta[k-1], m=Mu[k-1], a=Alpha[k-1], initdist = initdist)
            
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
        
        # draw new binomial sampling probability parameter
        probs[k] <- rbeta(1,p.prior[1] + sum(W.cur[,2]), p.prior[2] + sum(W.cur[,3]-W.cur[,2]))
        
        # draw new rate parameters 
        Beta[k] <- update_beta(X.cur = X.cur, beta.prior = beta.prior, popsize = popsize)
        
        Mu[k] <- update_mu(X.cur = X.cur, mu.prior = mu.prior)
        
        #   Alpha[k] <- update_alpha(X.cur = X.cur, alpha.prior = alpha.prior, popsize = popsize)
        
        Alpha[k] <- 0
        
        loglik[k] <- calc_loglike(X = X.cur, W = W.cur, irm = popirm.cur, b = Beta[k], m = Mu[k], a = Alpha[k], p = probs[k], initdist = initdist, popsize = popsize)  
        
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
