
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
    
    # create vector of ids for subjects to be infected
    whichsick <- 1:min(popsize,floor(sum(W[,2])/p)*amplify); totalinfected <- length(whichsick)
    
    # initialize matrix
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
    
    X[X[,1] > tmax, c(1,3)] <- 0 # censor observations greater than tmax
    
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

# insertRow inserts a row into a data frame
insertRow <- function(df, newrow, j) {
    
    df <- rbind(df, newrow, deparse.level=0)
    
    df <- df[order(c(1:(nrow(df)-1), j-0.5)), ]
    
    return(df)
}