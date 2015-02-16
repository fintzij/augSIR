# SIRsim simulates an SIR epidemic

h1 <- function(a, b, St, It){
  (b*It + a)*St
}

h2 <- function(mu, It){
  mu*It
}

SIRsim <- function(popsize, S0, I0, b, mu, a=0, tmax, censusInterval, sampprob, binomsamp = TRUE, returnX = FALSE) {
    
    if(binomsamp == FALSE) {
        X <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize), observed = rep(0,2*popsize)))
    
        } else if(binomsamp == TRUE) {
        X <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
        
    }
    for(k in seq(1,(2*I0 - 1),by=2)){
        X[k, 3] <- 1 #record an infection, infection time for these cases is zero
        
    }
    
    X <- X[order(X[,1],X[,3]),]
    
    h1t <- h1(a=a,b=b,St=S0,It=I0) 
    h2t <- h2(mu=mu,It=I0)
    
    keep.going <- TRUE; timenow <- 0

    if(I0 != popsize){
        which.susc <- (I0+1):popsize
    } else {
        which.susc <- 0
    }
    
    which.inf <- 1:I0
    infectednow <- I0; susceptiblenow <- S0
    
    iTi <- 0
    
    while(keep.going == TRUE){
        rate <- h1t + h2t
        tau <- rexp(1, rate=rate)
        
        p <- runif(1); probs <- cumsum(c(h1t, h2t))/rate
        iTi <- iTi + tau*infectednow
        
        if (p <= probs[1]) { #infection happens
            timenow <- timenow + tau; infectednow <- infectednow + 1; susceptiblenow <- susceptiblenow - 1
            X[which(X[,2] == which.susc[1])[1],1] <- ifelse(timenow <= tmax, timenow, 0)
            X[which(X[,2] == which.susc[1])[1],3] <- ifelse(timenow <= tmax, 1, 0)
            
            h1t <- h1(a=a, b=b, St=susceptiblenow, It=infectednow); h2t <- h2(mu=mu, It=infectednow)
            which.inf <- c(which.inf, which.susc[1]); which.susc <- which.susc[-1]; 
            X <- X[order(X[,1],X[,3]),]
            
        } else if(p>probs[1]){ # recovery happens
            
            
            timenow <- timenow + tau; infectednow <- infectednow - 1; susceptiblenow <- susceptiblenow
#             durations <- timenow - X[X[,3]!=0 & X[,2] %in% which.inf, 1]; durations <- cumsum(durations)/sum(durations)
#             
#             draw <- runif(1)
#             
#             whorecovers <- Position(function(x) x >= draw, durations)
# 
            whorecovers <- sample(1:length(which.inf),1)

            X[which(X[,2] == which.inf[whorecovers])[1],1] <- ifelse(timenow <= tmax, timenow, 0)
            X[which(X[,2] == which.inf[whorecovers])[1],3] <- ifelse(timenow <= tmax, -1, 0)
            
            h1t <- h1(a=a, b=b, St=susceptiblenow, It=infectednow); h2t <- h2(mu=mu, It=infectednow)
            which.inf <- which.inf[-whorecovers]
            X <- X[order(X[,1],X[,3]),]
        }
        
        if(timenow > tmax | infectednow == 0) keep.going <- FALSE
    }
    
    if(binomsamp == FALSE){
        wasinfected <- unique(X[which(X[,3]==1),2]); observed <- ifelse(runif(length(wasinfected))<= sampprob, 1, 0)
        for(j in 1:length(wasinfected)){
            X[which(X[,2]==wasinfected[j]),4] <- observed[j]
        }
        
        SIRres <- data.frame(time = seq(0, tmax, by = censusInterval),
                             Observed = 0,
                             Truth = 0)
        for(k in 1:dim(SIRres)[1]){
            SIRres$Observed[k] <- sum(X[which(X[,1] <= SIRres$time[k] & X[,4]==1),3])
            SIRres$Truth[k] <- sum(X[which(X[,1] <= SIRres$time[k]),3])
            
        }
        
        
    } else if(binomsamp == TRUE){
        SIRres <- data.frame(time = seq(0, tmax, by = censusInterval),
                             Observed = 0, 
                             Truth = 0)
        for(k in 1:dim(SIRres)[1]){
            SIRres$Observed[k] <- rbinom(n = 1, size = sum(X[which(X[,1] <= SIRres$time[k]),3]), prob = sampprob)
            SIRres$Truth[k] <- sum(X[which(X[,1] <= SIRres$time[k]),3])
            
        }
        
        SIRres$Observed[1] <- max(1,SIRres$Observed[1])
        
    }
    
    if(any(SIRres[,2]==0)){
        ind <- ifelse(any(SIRres[,2]==0 & c(0, diff(SIRres[,2]))==0), which(SIRres[,2]==0 & c(0, diff(SIRres[,2]))==0), length(SIRres[,2]))
        SIRres <- SIRres[1:ind,]
    }
    
    if(returnX == FALSE){
        return(SIRres)
    } else {
        return(list(results = SIRres, trajectory = X, iTi= iTi))
    }
}

# sim_one_SIR simulates a single SIR trajectory for one individual
sim_one_SIR <- function(numsick, eventtimes, obstimes, b, m, initdist, returnpath = FALSE){
    Xt <- rep(1, length(obstimes))
    
    initstate <- sample.int(3,1,prob=initdist)
    
    if(initstate == 2){
        tau <- rexp(1, m)
        Xt[obstimes <= tau] <- 2
        Xt[obstimes > tau] <- 3
        
    } else if(initstate == 1){
        path <- c(0,0)
        cur.time <- 0; ind <- 1
        keep.going <- TRUE
        rate <- b*numsick[1]
        
        while(keep.going == TRUE){
            if(ind < length(eventtimes)){
                tau <- rexp(1, rate)
                if((cur.time + tau) < eventtimes[ind+1]){
                    path[1] <- cur.time + tau
                    cur.time <- cur.time + tau
                    keep.going <- FALSE
                    
                    path[2] <- cur.time + rexp(1, m)
                    
                    Xt[obstimes <= path[1]] <- 1
                    Xt[(obstimes > path[1]) & (obstimes <= path[2])] <- 2
                    Xt[obstimes > path[2]] <- 3
                    
                    
                } else if((cur.time + tau) >= eventtimes[ind + 1]){
                    cur.time <- eventtimes[ind + 1]
                    ind <- ind+1
                    rate <- b*numsick[ind]
                    
                }
            } else if(ind >= length(eventtimes)){
                keep.going <- FALSE
                
                path <- c(Inf, Inf)
                
                Xt[1:length(Xt)] <- 1
                
            }            
            
        }
        
    }
    
    if(returnpath == FALSE) {
        return(Xt)
    } else{
        return(path)
    }
}

# sim_one_tpms simulates infection statuses at a sequence of observation times directly from transition probability matrices
sim_one_tpms <- function(tpms, initdist) {
    states <- vector(0, length = dim(tpms)[3])
    
    states[1] <- sample.int(n = 3, size = 1, prob = initdist)
    
    for(s in 2:length(states)) {
        states[s] <- sample.int(n = 3, size = 1, prob = tpms[k-1,,s])
        
    }
    
    return(states)
}