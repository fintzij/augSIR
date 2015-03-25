# SIRsim simulates an SIR epidemic

infec_rate <- function(a, b, St, It){
  (b*It + a)*St
}

recov_rate <- function(mu, It){
  mu*It
}

SIRsim <- function(popsize, initdist, b, mu, a=0, tmax, censusInterval, sampprob, binomsamp = TRUE, returnX = FALSE) {
    
    # create a matrix to store the individual level trajectories. The first column is the time of infection or recovery, the second column contains subject ids, the 
    # third column records 1 for an infection, -1 for recovery, 0 for no event. The matrix is ordered according to time, then event (only relevant for time=0). 
    
    if(binomsamp == FALSE) {
        # if individuals are sampled individually, there will be another column with an indicator for whether a subject is observed
        X <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize), observed = rep(0,2*popsize)))
    
        } else if(binomsamp == TRUE) {
        # binomial sampling occurs at the end, conditional on the epidemic
        X <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
        
    }
    
    # sample initial number of infected and susceptible individuals according to multinomial distribution
    initcounts <- rmultinom(1, size = popsize, prob = initdist)
    susceptiblenow <- initcounts[1]; infectednow <- initcounts[2]
    
    # if the initial number of infecteds is zero, redraw since there is no need to consider an impossible epidemic
    if(infectednow == 0){
        while(infectednow == 0){
            initcounts <- rmultinom(1, size = popsize, prob = initdist)
            infectednow <- initcounts[2]; susceptiblenow <- initcounts[1]        
        }
    }
    

    X[seq(1,(2*infectednow - 1),by=2), 3] <- 1 #record an infection, infection time for these cases is zero
        
    # infection and recovery rates
    infec.rate <- infec_rate(a=a, b=b, St=susceptiblenow, It=infectednow) 
    recov.rate <- recov_rate(mu=mu, It=infectednow)
    
    keep.going <- TRUE; timenow <- 0

    if(infectednow != popsize){
        which.susc <- (infectednow+1):popsize
    } else {
        which.susc <- 0
    }
    
    which.inf <- 1:infectednow
    
    iTi <- 0
    
    while(keep.going == TRUE){
        
        # the hazard is the sum of the rate of a new infection and the rate or a new recovery
        hazard <- infec.rate + recov.rate
        
        # sample the time of the next infection or recovery
        tau <- rexp(1, rate = hazard)
        
        # conditional on the 
        event <- sample.int(n = 2, size = 1, prob = c(infec.rate, recov.rate))
        
        iTi <- iTi + tau*infectednow
        
        if (event == 1) { #infection happens
            # update the current time, count of infecteds (by adding 1), and count of susceptibles (by subtracting 1)
            timenow <- timenow + tau; infectednow <- infectednow + 1; susceptiblenow <- susceptiblenow - 1
            
            # assign the infection to the next susceptible individual
            X[which(X[,2] == which.susc[1])[1],1] <- timenow # record the infection time
            X[which(X[,2] == which.susc[1])[1],3] <- 1 # record the infection code
            
            # update the infection and recovery rates
            infec.rate <- infec_rate(a=a, b=b, St=susceptiblenow, It=infectednow); recov.rate <- h2(mu=mu, It=infectednow)
            
            # remove the infected individual from the list of susceptibles and add him to the list of infecteds
            which.inf <- c(which.inf, which.susc[1]); which.susc <- which.susc[-1]; 
            
        } else if(event == 2){ # recovery happens
            
            # update time and count of infecteds (subtracting 1). No change to the number of susceptibles.
            timenow <- timenow + tau; infectednow <- infectednow - 1
            
            # choose an individual to recover (discrete uniform over the infected individuals, by the memoryless property of exponentials and the fact that all 
            # individuals recover at the same rate)
            whorecovers <- sample(which.inf, 1)
            
            # record the time and recovery status for the infected individual
            X[which(X[,2] == whorecovers)[2],1] <- timenow #record the recovery time 
            X[which(X[,2] == whorecovers)[2],3] <- -1 #record the recovery if it occurs before tmax
            
            # update infection and recovery rates
            infec.rate <- infec_rate(a=a, b=b, St=susceptiblenow, It=infectednow); recov.rate <- recov_rate(mu=mu, It=infectednow)
            
            #remove recovered individual from list of infecteds
            which.inf <- which.inf[which.inf != whorecovers]
        }
        
        # no more infected individuals, so stop
        if(infectednow == 0) keep.going <- FALSE
    }
    
    # censor individuals outside of observation window
    X[X[,1] > tmax, c(1,3)] <- 0
    X <- X[order(X[,1],X[,3]),]
    
    # individual level binomial sampling
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
        
        
    } else if(binomsamp == TRUE){ #population level binomial sampling (default)
        SIRres <- data.frame(time = seq(0, tmax, by = censusInterval),
                             Observed = 0, 
                             Truth = 0)
        
        # count the true number of infected individuals at each observation time
        for(k in 1:dim(SIRres)[1]){
            SIRres$Truth[k] <- sum(X[which(X[,1] <= SIRres$time[k]),3])
            
        }
        
        #generate binomial samples
        SIRres$Observed <- rbinom(n = length(SIRres$Truth), size = SIRres$Truth, prob = sampprob)
        
        # If first observed count is zero, resample
        if(SIRres$Observed[1] == 0){
            while(SIRres$Observed[1] == 0){
                SIRres$Observed[1] <- rbinom(n = 1, size = SIRres$Truth[1], prob = sampprob)
            }
        }
        
    }
    
    # Censor the observations after two consecutive observations of 0
#     if(any(SIRres[,2]==0)){
#         ind <- ifelse(any(SIRres[,2]==0 & c(0, diff(SIRres[,2]))==0), which(SIRres[,2]==0 & c(0, diff(SIRres[,2]))==0), length(SIRres[,2]))
#         SIRres <- SIRres[1:ind,]
#     }

    # Get rid of the matrix after the infection has died out (no more infecteds)
    if(any(SIRres[,3] == 0)){
        ind <- which(SIRres[,3] == 0)[1]
        SIRres <- SIRres[1:ind,]
    }
    
    if(returnX == FALSE){
        return(SIRres)
    } else {
        return(list(results = SIRres, trajectory = X, iTi= iTi))
    }
}

# sim_one_SIR simulates a single SIR trajectory for one individual
sim_one_SIR <- function(Xcount, obstimes, b, m, initdist, tmax, returnpath = FALSE){
    Xt <- rep(1, length(obstimes))
    eventtimes <- c(Xcount[,1], tmax); numsick <- Xcount[,2]
    
    initstate <- sample.int(n=3, size=1, prob=initdist)
    
    if(initstate == 2){
        tau <- rexp(1, m)
        path <- c(0,tau)
        Xt[obstimes <= tau] <- 2
        Xt[obstimes > tau] <- 3
        
    } else if(initstate == 1){
        path <- c(0,0)
        cur.time <- 0; ind <- 1
        keep.going <- TRUE
        rate <- b*numsick[ind]
        
        while(keep.going == TRUE){
            
            if(cur.time < tmax){
                
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
                    
                    if(numsick[ind] == 0){
                        keep.going <- FALSE
                        path <- c(0,0)
                        Xt[1:length(Xt)] <- 1
                    }
                }
                
            } else if(cur.time >= tmax){
                keep.going <- FALSE
                
                path <- c(0,0)
                
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
    states <- rep(0, dim(tpms)[3] + 1)
    
    states[1] <- sample.int(n = 3, size = 1, prob = initdist)
    
    for(s in 1:dim(tpms)[3]) {
        states[s+1] <- sample.int(n = 3, size = 1, prob = tpms[states[s],,s])
        
    }
    
    return(states)
}