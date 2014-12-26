# SIRsim simulates an SIR epidemic

h1 <- function(alpha, b, St, It){
  (b*It + alpha)*St
}

h2 <- function(mu, It){
  mu*It
}


SIRsim <- function(N, S0, I0, b, mu, a=0, g=0, maxtime, censusInterval){
  
  SIRres <- matrix(c(0,S0,I0),nrow=1); ind=1; 
  h1t <- h1(a,b,S0,I0) 
  h2t <- h2(mu,I0)
  
  timenow <- SIRres[ind,1] ; susceptiblenow <- SIRres[ind,2]; infectednow <- SIRres[ind,3]
  
  while(SIRres[ind,][1] < maxtime & all(SIRres[ind,2:3]!=0) & susceptiblenow!=0 & infectednow!=0){
    if(susceptiblenow==0){
      rate <- h2t
      tau <- rexp(1,rate=rate)
      
      timenow <- timenow + tau; susceptiblenow <- susceptiblenow ; infectednow <- infectednow - 1
      h2t <- h2(mu, infectednow);
      
      if(SIRres[ind,1]<(ind*censusInterval) & timenow>(ind*censusInterval) & timenow<maxtime) {
        SIRres <- rbind(SIRres,c(ind*censusInterval,susceptiblenow,infectednow))
        ind <- ind + 1
      }
    } else if(susceptiblenow!=0){
      rate <- sum(h1t, h2t)
      tau <- rexp(1, rate=rate)
      
      p <- runif(1); probs <- cumsum(c(h1t, h2t)/rate)
      
      if(p < probs[1]){
        timenow <- timenow + tau; susceptiblenow <- susceptiblenow - 1; infectednow <- infectednow +1
        h1t <- h1(a, b, susceptiblenow, infectednow); h2t <- h2(mu, infectednow);
      } else {
        timenow <- timenow + tau; susceptiblenow <- susceptiblenow ; infectednow <- infectednow - 1
        h1t <- h1(a, b, susceptiblenow, infectednow); h2t <- h2(mu, infectednow);
      }
      if(SIRres[ind,1]<(ind*censusInterval) & timenow>(ind*censusInterval) & timenow<maxtime) {
        SIRres <- rbind(SIRres,c(ind*censusInterval,susceptiblenow,infectednow))
        ind <- ind + 1
      }
    }
    if(infectednow==0 & SIRres[ind,3]!=0){
      SIRres <- rbind(SIRres,c((ind*censusInterval),susceptiblenow,infectednow))
    }
  }
  return(SIRres)
}


SIRsim2 <- function(popsize, S0, I0, b, mu, a=0, tmax, censusInterval, prob, returnX = FALSE){
    X <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize), observed = rep(0,2*popsize)))
    
    for(k in seq(1,(2*I0 - 1),by=2)){
        X[k, 3] <- 1 #record an infection
        X[k+1, 1] <- rexp(1, rate=mu) ; X[k+1, 3] <- -1  # record a recovery and the recovery time
        X[k, 4] <- X[k+1, 4] <- ifelse(runif(1) <= prob, 1, 0) #record whether the infection and recovery were observed
        
    }
    
    X <- X[order(X[,1],X[,3]),]
    timenow <- 0 ; susceptiblenow <- S0; infectednow <- I0; whichsusc <- (I0 + 1):popsize
    
    k <- I0 + 1; keep.going <- TRUE
    while(keep.going==TRUE){
        timeseq <- unique(X[,1]); irm <- buildirm(X, b, mu, a, pop=TRUE)
        
        for(j in 1:(length(timeseq)-1)){
            eventtime <- timeseq[j] - log(1-runif(1))/irm[1,2,j]
            if(eventtime < timeseq[j+1]) break
        }
        
        if (eventtime>timeseq[length(timeseq)] | eventtime > tmax) {
            eventtime <- Inf
            keep.going <- FALSE # this condition would indicate that the epidemic has died off
            break
        } else if(eventtime <= timeseq[length(timeseq)]){
            X[which(X[,2]==k)[1],c(1,3)] <- c(eventtime, 1) # record infection time and event code
            X[which(X[,2]==k)[2],c(1,3)] <- c(eventtime - log(1-runif(1))/irm[2,3,1], -1) # record recovery time and event code
            X[which(X[,2]==k), 4] <- ifelse(runif(1) <= prob, 1, 0) # record whether infection was observed
            
            X <- X[order(X[,1],X[,3]),]
            k <- k+1
            
            if(k > popsize){
                keep.going <- FALSE  
            } 
        }
    }
        
    SIRres <- data.frame(time = seq(0, tmax, by = censusInterval),
                         infec.obs = 0,
                         infec.true = 0)
    for(k in 1:dim(SIRres)[1]){
        SIRres$infec.obs[k] <- sum(X[which(X[,1] <= SIRres$time[k] & X[,4]==1),3])
        SIRres$infec.true[k] <- sum(X[which(X[,1] <= SIRres$time[k]),3])
        
    }
    
    if(any(SIRres[,2]==0 & SIRres[,3]==0)){
        SIRres <- SIRres[1:(min(which(SIRres[,2]==0 & SIRres[,3]==0))),]
    }
    
    if(returnX == FALSE){
        return(SIRres)
    } else {
        return(list(results = SIRres, trajectory = X))
    }
}