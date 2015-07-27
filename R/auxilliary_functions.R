
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

# initializeX initializes the matrix of trajectories in the population by 
# simulating configurations of trajectories until a valid configuration of
# trajectories is acquired
initializeX <- function(W, b, mu, a, samp_prob, initdist, censusInterval, tmax, popsize){
    
    keep.going <- TRUE
    
    while(keep.going == TRUE){
        
        X <- SIRsim(popsize = popsize, initdist = initdist, b = b, mu = mu, a=0, tmax = tmax, censusInterval = censusInterval, sampprob = samp_prob, trim = TRUE, returnX = TRUE)
        
        obsmat <- X$results
        trajecs <- X$trajectory
        
        if(nrow(obsmat) < nrow(W)){
            
            keep.going <- TRUE
            
        } else if(nrow(obsmat) == nrow(W)){
            
            if(all(dbinom(W[,2], obsmat[,3], samp_prob) > 0)) {
                
                init_config <- trajecs
                keep.going <- FALSE
                
            } else {
                keep.going <- TRUE
            }
                
        } else if(nrow(obsmat) > nrow(W)){
            
            obsmat <- obsmat[1:nrow(W), ]
            
            if(all(dbinom(W[,2], obsmat[,3], samp_prob) > 0)) {
                
                init_config <- trajecs
                keep.going <- FALSE
                
            } else {
                keep.going <- TRUE
            }
            
        }
        
    }
    
    return(init_config)
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

# insertRow inserts a row into a data frame
insertRow <- function(df, newrow, j) {
    
    df <- rbind(df, newrow, deparse.level=0)
    
    df <- df[order(c(1:(nrow(df)-1), j-0.5)), ]
    
    return(df)
}



# Matrix consistency checks -----------------------------------------------

# check_subj_pop_mats checks that the count matrix is consistent with the matrix of individual level trajectories

check_subj_pop_mats <- function(X, Xcount){
    
    # in small populations, it is possible for Xcount to be a vector
    if(is.null(dim(Xcount))){
        
        times <- Xcount[1]
        popsize <- length(unique(X[,2]))
        
        inf.vec <- sum(X[X[,1] <= times, 3]) 
        susc.vec <- popsize - sum(X[X[,1] <= times[m], 3] == 1) 
        
        cond <- (inf.vec == Xcount[2]) & (susc.vec == Xcount[3])
        
    } else{
        
        inf.vec <- rep(0, nrow(Xcount)) # vector of count of infecteds
        susc.vec <- rep(0, nrow(Xcount)) # vector of count of susceptibles
        
        times <- unique(Xcount[,1]) # vector of event times
        popsize <- length(unique(X[,2]))
        
        for(m in 1:length(times)){
            
            inf.vec[m] <- sum(X[X[,1] <= times[m], 3]) 
            susc.vec[m] <- popsize - sum(X[X[,1] <= times[m], 3] == 1) 
            
        }
        
        cond <- all(inf.vec == Xcount[,2]) & all(susc.vec == Xcount[,3])
    }
    

    return(cond)
}


# check_count_obs_mats checks that the observation matrix is consistent with the matrix of counts

check_count_obs_mats <- function(Xcount, W){
    
    obs.inf <- rep(0, nrow(W)) # vector of true number of infecteds at observation times
    
    times <- W[,1] # vector of observation times
    
    for(m in 1:length(times)){
        
        obs.inf[m] <- Xcount[sum(Xcount[,1] <= times[m]), 2]
        
    }
    
    if(all(obs.inf == W[,3])) {
        
        return(TRUE)
        
    } else return(FALSE)    
}

