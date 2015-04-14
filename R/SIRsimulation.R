# SIRsim simulates an SIR epidemic

infec_rate <- function(a, b, St, It){
  (b*It + a)*St
}

recov_rate <- function(mu, It){
  mu*It
}

SIRsim <- function(popsize, initdist, b, mu, a=0, tmax, censusInterval, sampprob, trim = TRUE, returnX = FALSE) {
    
    # create a matrix to store the individual level trajectories. The first column is the time of infection or recovery, the second column contains subject ids, the 
    # third column records 1 for an infection, -1 for recovery, 0 for no event. The matrix is ordered according to time, then event (only relevant for time=0). 
    
    if(returnX == TRUE){
        X <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
        
    }
    
    # create matrix to store the simulation results. Observed will contain the binomial samples, time is the time of recovery or infection, truth is the number of infecteds
    # at each event time. 
    
    SIRres <- data.frame(time = seq(0, tmax, by = censusInterval),
                         Observed = 0, 
                         Truth = 0)
    
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
    
    # record the initial number of infecteds
    SIRres$Truth[1] <- infectednow
    
    # if returnX == TRUE, initialize the vectors of IDs for infecteds and susceptibles
    if(returnX == TRUE){
        X[seq(1,(2*infectednow - 1),by=2), 3] <- 1 #record an infection, infection time for these cases is zero
        
        # initialize vector of susceptible IDs
        if(infectednow != popsize){
            which.susc <- (infectednow+1):popsize
        } else {
            which.susc <- 0
        }
        
        # initialize vector of IDs for infecteds
        which.inf <- 1:infectednow
        
    }
        
    # infection and recovery rates
    infec.rate <- infec_rate(a=a, b=b, St=susceptiblenow, It=infectednow) 
    recov.rate <- recov_rate(mu=mu, It=infectednow)
    
    keep.going <- TRUE; timenow <- 0
    
    
    while(keep.going == TRUE){
        
        # the hazard is the sum of the rate of a new infection and the rate or a new recovery
        hazard <- infec.rate + recov.rate
        
        # sample the time of the next infection or recovery
        tau <- rexp(1, rate = hazard)
                
        # conditional on an event happening, sample whether it is an infection or a recovery
        event <- sample.int(n = 2, size = 1, prob = c(infec.rate, recov.rate))
        
        
        if (event == 1) { #infection happens
            
            # update the current time, count of infecteds (by adding 1), and count of susceptibles (by subtracting 1)
            infectednow <- infectednow + 1; susceptiblenow <- susceptiblenow - 1
            
            # update results matrix
            SIRres$Truth[SIRres$time > timenow & SIRres$time <= (timenow + tau)] <- infectednow
            
            # update time
            timenow <- timenow + tau
            
            # update the infection and recovery rates
            infec.rate <- infec_rate(a=a, b=b, St=susceptiblenow, It=infectednow); recov.rate <- recov_rate(mu=mu, It=infectednow)
            
            # if returnX == TRUE, update X
            if(returnX == TRUE){
                # assign the infection to the next susceptible individual
                X[which(X[,2] == which.susc[1])[1],1] <- timenow # record the infection time
                X[which(X[,2] == which.susc[1])[1],3] <- 1 # record the infection code
                
                # remove the infected individual from the list of susceptibles and add him to the list of infecteds
                which.inf <- c(which.inf, which.susc[1]); which.susc <- which.susc[-1]
            }
            
            
        } else if(event == 2){ # recovery happens
            
            # update time and count of infecteds (subtracting 1). No change to the number of susceptibles.
            infectednow <- infectednow - 1
            
            # update results matrix
            SIRres$Truth[SIRres$time > timenow & SIRres$time <= (timenow + tau)] <- infectednow
            
            # update time
            timenow <- timenow + tau
            
            # update infection and recovery rates
            infec.rate <- infec_rate(a=a, b=b, St=susceptiblenow, It=infectednow); recov.rate <- recov_rate(mu=mu, It=infectednow)
            
            # if returnX == TRUE, select an infected individual to recover and update the vector of infecteds
            if(returnX == TRUE){
                # choose an individual to recover (discrete uniform over the infected individuals, by the memoryless property of exponentials and the fact that all 
                # individuals recover at the same rate)
                whorecovers <- sample(which.inf, 1)
                
                # record the time and recovery status for the infected individual
                X[which(X[,2] == whorecovers)[2],1] <- timenow #record the recovery time 
                X[which(X[,2] == whorecovers)[2],3] <- -1 #record the recovery if it occurs before tmax
                
                #remove recovered individual from list of infecteds
                which.inf <- which.inf[which.inf != whorecovers]
            }            

        }
        
        # no more infected individuals so the epidemic has been fully observed and we stop
        if(infectednow == 0) keep.going <- FALSE
    }
    
    # if returnX == TRUE, we may need to censor some observations in the X matrix
    if(returnX == TRUE){
        # censor individuals outside of observation window
        X[X[,1] > tmax, c(1,3)] <- 0
        X <- X[order(X[,1],X[,3]),]
    }
    
    
    #generate binomial samples
    SIRres$Observed <- rbinom(n = length(SIRres$Truth), size = SIRres$Truth, prob = sampprob)
        
    # If first observed count is zero, resample that observation
    if(SIRres$Observed[1] == 0){
        while(SIRres$Observed[1] == 0){
            SIRres$Observed[1] <- rbinom(n = 1, size = SIRres$Truth[1], prob = sampprob)
        }
    }
        

    # Get rid of the matrix after the infection has died out (no more infecteds) if trim == TRUE
    if(trim == TRUE){
        if(any(SIRres[,3] == 0)){
            ind <- which(SIRres[,3] == 0)[1] - 1
            SIRres <- SIRres[1:ind,]
        }
    }

    
    if(returnX == FALSE){
        return(SIRres)
    } else {
        return(list(results = SIRres, trajectory = X))
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