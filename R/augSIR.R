
# Functions to draw event times (drawpath, drawtime) -------------------------------------------


# drawpath draws a path for an endpoint conditioned markov chain
# Xt is the subject level observation matrix with first column a vector of observation times, and second column a vector of infection
# statuses. Xcount is the matrix with the population trajectory for all other subjects. 
# irm is an array of instantaneous rate matrices. tmax is the maximum time of observation.

drawpath <- function(Xt, Xcount, irm, tmax){
    path <- c(0,0)
    changetimes <- which(diff(Xt[,2], lag = 1) != 0)
    
    if(length(changetimes)==0){
        if(all(Xt[,2]==2)){
            path[1] <- 0; t0 <- Xt[dim(Xt)[1],1]; t1 <- Inf
            path[2] <- drawtime(Xcount, irm, t0, t1, 2)
            path[2] <- ifelse(path[2]>tmax, Inf, path[2])
            
        } else if(all(Xt[,2]==1)){
            t0 <- Xt[dim(Xt)[1],1]; t1 <- Inf
            path[1] <- drawtime(Xcount, irm, t0, t1, 1)
            
            if(path[1]>tmax) {
                path[1] <- path[2] <- Inf
            } else if(path[1] <= tmax){
                path[2] <- drawtime(Xcount, irm, path[1], t1, 2)
                
                if(path[2] > tmax) path[2] <- Inf
                
            }
            
        }
        
    } else if(length(changetimes)==1){ ### all transitions occur in one observation period
        if((Xt[changetimes+1,2] - Xt[changetimes,2]) ==2){ ### Subject transitions from healthy to recovered within one observation period
            t0 <- Xt[changetimes,1]; t1 <- Xt[changetimes+1,1]
            path[1] <- drawtime(Xcount, irm, t0, t1, 1)
            path[2] <- drawtime(Xcount, irm, path[1], t1, 2)
            
        } else if(Xt[changetimes,2]==1 & ((Xt[changetimes+1,2] - Xt[changetimes,2]) == 1)){ # healthy to infected within one observation period
            t0 <- Xt[changetimes,1]; t1 <- Xt[changetimes+1,1]
            path[1] <- drawtime(Xcount, irm, t0, t1, 1)
            
            t0 <- Xt[dim(Xt)[1],1]; t1 <- Inf
            path[2] <- drawtime(Xcount, irm, t0, t1, 2)
            
            if(path[2] > tmax) path[2] <- Inf
            
        } else if(Xt[changetimes,2]==2 & ((Xt[changetimes+1,2] - Xt[changetimes,2]) ==1)){ #infected to recovered within one observation period
            path[1] <- 0; t0 <- Xt[changetimes,1]; t1 <- Xt[changetimes+1,1]
            path[2] <- drawtime(Xcount, irm, t0, t1, 2)
            
        }
        
    } else if(length(changetimes)==2){
        
        t0 <- Xt[changetimes[1],1]; t1 <- Xt[changetimes[1]+1,1]
        path[1] <- drawtime(Xcount, irm, t0, t1, 1)
        
        t0 <- Xt[changetimes[2],1]; t1 <- Xt[changetimes[2]+1,1]
        path[2] <- drawtime(Xcount, irm, t0, t1, 2)
    }  
    
    return(path)
}


# drawtime draws an event time for either carriage acquisition or clearance. Xother is a matrix with the times of acquisition and clearance
# for all other subjects. irm is an array of instantaneous rate matrices. t0 and t1 are the left and right endpoints of the interval. 
# current state is the infection status at t0. 

drawtime <- function(Xcount, irm, t0, t1, currentstate){
    if(currentstate == 2){
        eventtime <- t0 - log(1 - runif(1)*(1-exp(-irm[2,3,1] * (t1 - t0))))/irm[2,3,1]
        
    } else if(currentstate ==1){
        times <- Xcount[,1]; numinf <- Xcount[,2]
        timeseq <- c(t0, times[times > t0 & times < t1], t1)
        
        if(t1 != Inf){
            if(length(timeseq)!=2) {
                
                intervalprobs <- rep(0,length(timeseq)-1)
                timediffs <- diff(timeseq, lag = 1)
                indstart <- sum(times<=timeseq[1]); indend <- sum(times<=timeseq[length(timeseq)-1])
                totalprob <- 1 - exp(-sum(irm[1, 2,numinf[indstart:indend] + 1]*timediffs))
                
                for(k in 1:length(intervalprobs)){
                    ind <- sum(times<=timeseq[k])
                    intervalprobs[k] <- (1-exp(-sum(irm[1, 2, numinf[indstart:ind] + 1]*timediffs[1:k])))/totalprob
                }
                
                interval <- sample.int(n = length(intervalprobs), size = 1, prob = intervalprobs)
                ind <- sum(times <= timeseq[interval])
                
                if(numinf[ind] != 0) {
                    eventtime <- timeseq[interval]-log(1-runif(1)*(1-exp(-irm[1, 2, numinf[ind] + 1]*timediffs[interval])))/irm[1, 2, numinf[ind] + 1]
                    
                } else if(numinf[ind] == 0) {
                    eventtime <- Inf
                    
                }
                
            } else if(length(timeseq)==2){
                timediffs <- t1 - t0; ind <- sum(times <= t0)
                
                if(numinf[ind] != 0){
                    eventtime <- t0 - log(1 - runif(1)*(1 - exp(-irm[1, 2, numinf[ind] + 1]*(t1 - t0))))/irm[1, 2, numinf[ind] + 1]
                    
                    
                } else if(numinf[ind] == 0){
                    eventtime <- Inf
                    
                }
                
            }
            
        } else if(t1 == Inf){
            for(k in 1:(length(timeseq)-1)){
                ind <- sum(times<=timeseq[k])
                
                eventtime<- timeseq[k] - log(1-runif(1))/irm[1,2, numinf[ind] + 1]
                if(eventtime < timeseq[k+1]) break
            }
            
        }        
        
    }
    
    return(eventtime)
    
}

# Forward-Backward Functions (forward, backward, drawXt) ----------------------------------------------

# forward builds the forward matrices as in Scott (2002), taking as arguments an array of tpms, a matrix of emission probabilities, and an initial distribution over the state at t0
fwd <- function(tpms, emits, initdist){
    fbmats <- array(0, dim = dim(tpms))
    
    pi0 <- normalize(initdist * emits[,1])
    fbmats[,,1] <- normalize(outer(pi0,emits[,2], FUN="*") * tpms[,,1])
    
    for(s in 2:dim(fbmats)[3]){
        fbmats[,,s] <- normalize(outer(colSums(fbmats[,,s-1]), emits[,s+1], FUN="*") * tpms[,,s])
    
    }
    
    return(fbmats)
}

# bwd takes as an argument the array of matrices produced by forward and returns a draw at each observation time. The default argument for the prior on end states is the unit vector
bwd <- function(mats, end.prior = c(1, 1, 1)){
    states <- rep(1, dim(mats)[3]+1)
    
    states[length(states)] <- sample.int(n = 3, size = 1, prob = end.prior%*%mats[,,dim(mats)[3]])
    
    for(s in (length(states)-1):1){
        states[s] <- sample.int(n = 3, size = 1, prob = mats[,states[s+1],s])
        
    }

    return(states)
} 

# fb_marginals calculates the marginal distributions \pi_s(t | theta), taking as arguments an array of fb matrices from forward, and a matrix with columns the initial dist at t0 and emission probs at t0
fb_marginals <- function(fbmats, t0.probs){
    marginals <- matrix(0, dim = c(3, dim(fbmats)[3]+1))
    
    marginals[,1] <- normalize(t0.probs[,1]*t0.probs[,2])
    
    marginals[,2:dim(marginals)] <- apply(fbmats,3, colSums)
}

# drawXt draws the status of a subject at observation times t1,...,tL conditional on other individuals
# Xcount is a matrix with the population trajectory for all other subjects. 
# p, is the sampling probability, and b, and m are the epidemiologic parameters. initdist is the initial
# distribution for infection status
drawXt <- function(Xcount, irm, irm.eig, W, p, initdist){
    
    Xt <- cbind(W[,1], 0)
    
    tpms <- tpm_seq(Xcount = Xcount, obstimes = W[,1], irm.eig = irm.eig)
    emits <- emit_seq(W = W, p = p)
    
    fbmats <- fwd(tpms = tpms, emits = emits, initdist = initdist)
    
    Xt[,2] <- bwd(mats = fbmats)
    
    return(Xt)
}

# Matrix Construction and Update Functions (build_countmat, build_irm, irm_decomp, update_irm, update_eigen, build_tpm, obs_tpm, tpm_seq, emit_seq, updateX, updateW) --------------------------------------

# build_countmat builds the count matrix containing the the numbers of susceptibles and infecteds at event times
build_countmat <- function(X, popsize) {
    initinfec <- sum(X[X[,1]==0,3])
    Xcount <- rbind(c(0, initinfec), X[X[,1]!=0,c(1,3)])
    Xcount[,2] <- cumsum(Xcount[,2])
    
    numsusc <- popsize - cumsum(c(initinfec, X[X[,1]!=0, 3]>0))
    
    Xcount <- cbind(Xcount, numsusc)
    
    return(Xcount)
    
}

# build_irm constructs an array of instantaneous rate matrices, with the bottom right entry denoting the number of infecteds. The top left (n-1)x(n-1) matrix is the irm.
build_irm <- function(Xcount, b, m, a=0, popsize, pop){
    
    if(pop == FALSE){ # construct irm array for an individual
        
        numinf <- 0:(max(Xcount[,2]) + 1)
        irm <- array(0, dim = c(4,4,length(numinf)))
        irm[1,2,] <- b * numinf + a; irm[1,1,] <- -irm[1,2,]
        irm[2,3,] <- m; irm[2,2,] <- -m
        irm[4,4,] <- numinf
        
    } else if(pop == TRUE){ # construct irm array for the population
        
        numinf <- Xcount[,2]
        numsusc <-  Xcount[,3]
        indend <- dim(Xcount)[1]
        
        irm <- array(0, dim = c(3,3, length(numinf) - 1))
        irm[1,2,] <- (b*numinf[1:(indend-1)] + a)*numsusc[1:(indend-1)]; irm[1,1,] <- -irm[1,2,]
        irm[2,3,] <- m*numinf[1:(indend-1)]; irm[2,2,] <- -irm[2,3,]
        
    }
    
    return(irm)
}

# irm_decomp takes an array of path irms and returns an array of eigenvalues, eigenvectors, and the inverse of the eigenvectors
# irm_decomp returns an array, each element of which is another array containing a diagonal matrix of eigenvalues, a matrix of eigenvectors, and 
# the inverse of the matrix of eigenvectors

irm_decomp <- function(pathirm.cur){
    
    irm.decomp <- array(0,dim = c(3,3,3,dim(pathirm.cur)[3]))
    to_reorder <- pathirm.cur[1,1,] > pathirm.cur[2,2,]
    
    for(s in 1:dim(pathirm.cur)[3]){
        
        decomp <- eigen(pathirm.cur[1:3, 1:3, s], symmetric=FALSE)
        
        if(to_reorder[s] == FALSE){
            irm.decomp[,,1,s] <- diag(decomp[[1]])
            irm.decomp[,,2,s] <- decomp[[2]]
            irm.decomp[,,3,s] <- solve(irm.decomp[,,2,s])
            
        } else if(to_reorder[s] == TRUE){
            irm.decomp[,,1,s] <- diag(decomp[[1]][c(2,1,3)])
            irm.decomp[,,2,s] <- decomp[[2]][,c(2,1,3)]
            irm.decomp[,,3,s] <- solve(irm.decomp[,,2,s])
            
        }
    }
    return(irm.decomp)
    
}

# update_irm updates the array of individual level irms, adding a new irm to the end of the array
update_irm <- function(irm, new.numinf, b, m, a, popsize){
    newmat <- diag(0,4)
    newmat[1,2] <- (b * new.numinf + a); newmat[1,1] <- -newmat[1,2]
    newmat[2,3] <- m; newmat[2,2] <- -m
    newmat[4,4] <- new.numinf
    
    irm.new <- array(c(irm,newmat), dim = c(4,4,dim(irm)[3]+1))
    
    return(irm.new)
}

# update_eigen updates the array of eigen decompositions associated with the array of individual level irms, adding a new array at the end
update_eigen <- function(patheigen, pathirm){
    ind <- dim(pathirm)[3]
    
    decomp <- eigen(pathirm[c(1:3),c(1:3),ind], symmetric = FALSE)
    
    if(pathirm[1,1,ind] > pathirm[2,2,ind]){
        decomp$values <- decomp$values[c(2,1,3)]
        decomp$vectors <- decomp$vectors[,c(2,1,3)]
    }
    
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
obstpm <- function(Xcount, irm.eig, t0, t1){
    eventtimes <- Xcount[,1]; numinf <- Xcount[,2]
    timeseq <- c(t0, eventtimes[(eventtimes > t0) & (eventtimes <t1)], t1)
    indstart <- sum(eventtimes <= t0); indend <- indstart + length(timeseq)-2
    inds <- numinf[indstart:indend] + 1
    
    tpm <- diag(1,3) 
    
    for (j in 1:length(inds)) {
        ind <- inds[j]
        tpm <- tpm %*% buildtpm(values = irm.eig[,,1,ind],
                                vectors = irm.eig[,,2,ind],
                                inv.vecs = irm.eig[,,3,ind],
                                t0 = timeseq[j],
                                t1 = timeseq[j+1])
        
    }
    
    return(tpm)
}

# tpm_seq constructs an array of transition probability matrices given a sequence of observation times, a matrix of counts in the form of Xcount, and an array of eigen decompositions for irms
tpm_seq <- function(Xcount, obstimes, irm.eig){
    tpms <- array(0, dim = c(3,3,length(obstimes)-1))
    
    for(s in 1:(length(obstimes)-1)){
        tpms[,,s] <- obstpm(Xcount = Xcount, irm.eig = irm.eig, t0 = obstimes[s], t1 = obstimes[s+1])
        
    }
    
    return(tpms)
    
}

# emit_seq constructs an array of emission matrices given an observation matrix W, the last column of which is the number of augmented trajectories and the binomial sampling probability
emit_seq <- function(W, p){
    emits <- matrix(0, nrow = 3, ncol = dim(W)[1])
    
    for(s in 1:dim(emits)[2]){
        emits[,s] <- dbinom(x = W[s, 2], size = W[s, 3] + c(0,1,0), prob = p)
        
    }
    
    return(emits)
}

# updateW updates an observation matrix. The function takes as arguments an observation matrix with counts, Xcount, and a matrix, W, 
# with columns containing observation times, number of infecteds observed, and number of infecteds based on augmented trajectories. 
updateW <- function(W, Xcount){
    for(k in 1:dim(W)[1]){
        W[k,3] <- Xcount[sum(Xcount[,1]<=W[k,1]),2]
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


# get_Wother obtains the observation matrix excluding subject j
get_W_other <- function(W.cur, path.cur){
    Wother <- W.cur
    
    if(all(path.cur != 0)){
        Wother[,3] <- Wother[,3] - ((Wother[,1] >= path.cur[1]) & (Wother[,1] <= path.cur[2]))
        
    }
    
    return(Wother)
} 

# get_Xcount_other obtains the count matrix excluding subject j
get_Xcount_other <- function(Xcount, path.cur){
    
    if(all(path.cur == 0)){
        Xcount.other <- Xcount
        
    } else if(path.cur[1] == 0 & path.cur[2] != 0){
        Xcount.other <- Xcount[Xcount[,1] != path[2],]
        Xcount.other[,2] <- Xcount.other[,2] - (Xcount.other[,1] <= path[2])
        
    } else if(path.cur)
    
    if(all(path.cur !=0)){
        Xcount.other <- Xcount[!(Xcount %in% path.cur), ]
        
        Xcount.other[,3] <- Xcount.other[,3] - (Xcount.other[,1] < path.cur[1])
        Xcount.other[,2] <- Xcount.other[,2] - ((Xcount.other[,1] >= path.cur[1]) & (Xcount.other[,1] <= path.cur[2]))
    }
    
    return(Xcount.other)    
}

# Metropolis-Hastings Ratio related functions (accept_prob, obs_prob, cacl_loglike, path_prob, pop_prob) ----------------------

# accept_prob calculates the acceptance ratio for the M-H algorithm
accept_prob <- function(pop_prob.new, pop_prob.cur, path_prob.cur, path_prob.new){
    if(pop_prob.new == -Inf & path_prob.new == -Inf){
        -Inf
    } else{
        pop_prob.new - pop_prob.cur + path_prob.cur - path_prob.new
    }
}

# obs_prob calculates the binomial probability of observing the data given counts of the infecteds at observation times
obs_prob <- function(W, p){
    dbinom(x = sum(W[,2]), size = sum(W[,3]), prob = p, log = TRUE)
    
}

# calc_loglik calculates the complete data log-likelihood
calc_loglike <- function(Xcount, W, irm, b, m, a=0, p, initdist, popsize){
    indend <- dim(Xcount)[1]; init.infec <- Xcount[1,2]
    
    events <- diff(Xcount[,2])
    rates <- ifelse(events==1, irm[1,2,], irm[2,3,]); #rates <- pmax(0,rates)
    
    hazards <- irm[1,2,] + irm[2,3,] 
    
    dbinom(sum(W[,2]), sum(W[,3]), prob=p, log=TRUE) + dmultinom(c(popsize - init.infec, init.infec, 0), prob = initdist, log=TRUE) + sum(log(rates)) - sum(hazards*diff(Xcount[,1], lag = 1))
} 


# pop_prob and path_prob calculate the log-probabilities of the population trajectory and the subject trajectory for use in the M-H ratio

pop_prob <- function(Xcount, irm, initdist, popsize){
    indend <- dim(Xcount)[1]; init.infec <- Xcount[1,2]
    
    events <- diff(Xcount[,2])
    rates <- ifelse(events==1, irm[1,2,], irm[2,3,]); #rates <- pmax(0,rates)
    
    hazards <- irm[1,2,] + irm[2,3,]        
    
    dmultinom(c(popsize - init.infec, init.infec, 0), prob = initdist, log=TRUE) + sum(log(rates)) - sum(hazards*diff(Xcount[,1], lag = 1))
    
}


path_prob <- function(path, Xcount, pathirm, initdist, tmax){
    indend <- dim(Xcount)[1]; init.infec <- Xcount[1,2]
    
    times <- Xcount[,1]; timediffs <- diff(times, lag = 1)
    numinf <- Xcount[,2]
    
    if(all(path==0)){ # no infection observed
        
        path.prob <- log(initdist[1])-sum(pathirm[1, 2, numinf[1:(indend - 1)] + 1] * timediffs)
        
    } else if(path[2] > 0 & path[2] < Inf) { # subject is infected and recovery is observed
        
        if(path[1]==0){ # subject is initially infected 
            
            path.prob <- log(initdist[2]) + log(pathirm[2,3,1]) - pathirm[2,3,1]*path[2]
            
        } else if(path[1]!=0){ # subject is not initially infected
            
            if(times[2] > path[1]){ # no changes in number of other infecteds before subject's infection
                
                path.prob <- log(initdist[1]) + log(pathirm[1,2, init.infec + 1]) - pathirm[1,2, init.infec+1]*path[1] + 
                    log(pathirm[2,3,1]) - pathirm[2,3,1]*(path[2] - path[1])
                
            } else if(times[2] < path[1]){ # if there is at least one change in the number of infecteds before the subject's infection
                ind1 <- which(times < path[1])[sum(times < path[1])]
                
                path.prob <- log(initdist[1]) - sum(pathirm[1,2, numinf[1:(ind1-1)] + 1] * timediffs[1:(ind1-1)]) + 
                    log(pathirm[1,2,numinf[ind1] + 1]) - pathirm[1,2, numinf[ind1]+1]*(path[1] - times[ind1]) + 
                    log(pathirm[2,3,1]) - pathirm[2,3,1]*(path[2] - path[1]) 
                
            }
        }
        
    } else if(path[2] == Inf) { # subject is infected and no recovery is observed
        
        if(path[1]==0){ # subject is initially infected 
            
            path.prob <- log(initdist[2]) - pathirm[2,3,1]*tmax
            
        } else if(path[1]!=0){ # subject is not initially infected
            
            if(times[2] > path[1]){ # no changes in number of infecteds before subject's infection
                
                path.prob <- log(initdist[1]) + log(pathirm[1,2, init.infec + 1]) - pathirm[1,2, init.infec + 1]*path[1] - pathirm[2,3,1]*(tmax - path[1])
                
            } else if(times[2] < path[1]){ # if there is at least one change in the number of infecteds before the subject's infection
                ind1 <- which(times < path[1])[sum(times < path[1])]
                
                path.prob <- log(initdist[1]) - sum(pathirm[1,2, numinf[1:(ind1-1)] + 1] * timediffs[1:(ind1-1)]) + 
                    log(pathirm[1,2, numinf[ind1] + 1]) - pathirm[1,2, numinf[ind1] + 1]*(path[1] - times[ind1])  -
                    pathirm[2,3,1]*(tmax - path[1])
                
            }
        }
    } 
    
    return(path.prob)
}


# Functions to update parameters (update_rates, update_prob) ---------------------------------
update_rates <- function(Xcount, beta.prior, mu.prior, alpha.prior = NULL, popsize){
    indend <- dim(Xcount)[1]
    infections <- diff(Xcount[,2], lag = 1)>0; recoveries <- !infections
    
    numsick <- Xcount[1:(indend - 1),2]; numsusc <- Xcount[1:(indend-1),3] 
    timediffs <- diff(Xcount[,1], lag = 1)
    
    beta.new <- rgamma(1, shape = (beta.prior[1] + sum(infections)), 
                       rate = beta.prior[2] + sum(numsick * numsusc * timediffs))
    
    mu.new <- rgamma(1, shape = (mu.prior[1] + sum(recoveries)), 
                     rate = mu.prior[2] + sum(numsick * timediffs))
    
    params.new <- c(beta.new, mu.new,0)
    
    if(!is.null(alpha.prior)){
        alpha.new <- rgamma(1, shape = (alpha.prior[1] + sum(infections)), 
                            rate = alpha.prior[2] + sum(numsusc * timediffs))
        
        params.new[3] <- alpha.new
    }
    
    return(params.new)
}

# update_prob updates the binomial sampling probability parameter
update_prob <- function(W, p.prior){
    rbeta(1,p.prior[1] + sum(W[,2]), p.prior[2] + sum(W[,3]-W[,2]))
    
}


# Auxilliary Functions (normalize, checkpossible, find.pprior, find.rateprior, initializeX, get_path, infec_durs) -----------------------------------------

# Normalize is what it sounds like.
normalize <- function(X){
    X/sum(X)
}

# checkpossible checks whether there is an impossible population trajectory. 
checkpossible <- function(Xcount, a=0, W=NULL){
    popseq <- Xcount[,2]
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
        
        if(any(W[,3]<W[,2])){
            is.possible <- FALSE
        }
    
    }
    return(is.possible)
}

# functions to find parameters of prior distributions
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

# initialize_X initializes the matrix of trajectories in the population by dispersing infecteds throughout the observation period
# dat is a matrix with observation times and counts of infecteds, mu is the recovery parameter, p is the sampling probability, amplify
# controls the factor by which the observed number of infecteds is multiplied to get the number of initialized trajectories (up to the population size)
# X is matrix with event times, subject id and event codes. 
# Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
initializeX <- function(W, b, mu, a, p, amplify, tmax, popsize){
    whichsick <- 1:min(popsize,floor(sum(W[,2])/p)*amplify); totalinfected <- length(whichsick)
    
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
                        durs <- infec_durs(X, ids)
                        proposed <- sample(ids, 1, prob = durs) # preferentially resample individuals with long infection times for reinitialization
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
            if(sum(X[X[,1]<=W[k,1],3]) < W[k,2]){
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
    
    if(sum(infec_durs(X,1:popsize)!=0) < totalinfected){
        
        X<-X[order(X[,1]),]; Xcount <- build_countmat(X = X, popsize = popsize)
        
        infecs_toadd <- which(infec_durs(X,1:totalinfected)==0)
        
        for(s in 1:length(infecs_toadd)){
            indend <- dim(Xcount)[1]
            
            infections <- (diff(Xcount[,2], lag = 1) == 1); recoveries <- !infections
            
            numsick <- Xcount[,2]; numsusc <- Xcount[,3]
            timediffs <- diff(Xcount[,1], lag = 1)
            
            irm <- build_irm(Xcount = Xcount, b = b, m = mu, a = a, popsize = popsize, pop = FALSE)
            
            path <- c(0,0)
            path[1] <- drawtime(Xcount = Xcount, irm = irm, t0=0, t1=max(X[,1]), currentstate = 1)
            path[2] <- ifelse(path[1] == Inf, Inf, drawtime(Xcount = Xcount, irm = irm, t0 = path[1], t1 = Inf, currentstate = 2))
            
            X <- updateX(X, path, infecs_toadd[s])
        }
    }
    
    X<-X[order(X[,1]),]
    
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

# infec_durs retrieves the durations of infections from a matrix X
infec_durs <- function(X,ids){
    durs <- rep(0, length = length(ids))
    for(r in 1:length(durs)){
        durs[r] <- max(X[X[,2]==ids[r],1]) - min(X[X[,2]==ids[r],1])
    }
    
    return(durs)
}


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
    W.cur <- updateW(W = W.cur, X = X.cur)
    
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
            
            Xcount.other <- get_Xcount_other(Xcount = Xcount, popsize = popsize)
            W.other <-get_W_other(W.cur = W.cur, path.cur = path.cur)
            
            Xt <- drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W=W.other, p=probs[k-1], initdist = initdist)
            
            path.new<- drawpath(Xt = Xt, Xcount = Xcount.other, irm = pathirm.cur, tmax = tmax)
            
            X.new <- updateX(X = X.cur, Xt.path = path.new, j = subjects[j]); path.new <- getpath(X = X.new, j = subjects[j])
            Xcount.new <- build_countmat(X = X.new, popsize = popsize)
            W.new <- updateW(W.cur,X.new)
            
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

