# SIRsim simulates an SIR epidemic

h1 <- function(a, b, St, It){
  (b*It + a)*St
}

h2 <- function(mu, It){
  mu*It
}

# SIRsim <- function(N, S0, I0, b, mu, a=0, g=0, maxtime, censusInterval){

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


# Deprecated functions ----------------------------------------------------

# SIRsim <- function(N, S0, I0, b, mu, a=0, g=0, maxtime, censusInterval){
#   
#   SIRres <- matrix(c(0,S0,I0),nrow=1); ind=1; 
#   h1t <- h1(a,b,S0,I0) 
#   h2t <- h2(mu,I0)
#   
#   timenow <- SIRres[ind,1] ; susceptiblenow <- SIRres[ind,2]; infectednow <- SIRres[ind,3]
#   
#   while(SIRres[ind,][1] < maxtime & all(SIRres[ind,2:3]!=0) & susceptiblenow!=0 & infectednow!=0){
#     if(susceptiblenow==0){
#       rate <- h2t
#       tau <- rexp(1,rate=rate)
#       
#       timenow <- timenow + tau; susceptiblenow <- susceptiblenow ; infectednow <- infectednow - 1
#       h2t <- h2(mu, infectednow);
#       
#       if(SIRres[ind,1]<(ind*censusInterval) & timenow>(ind*censusInterval) & timenow<maxtime) {
#         SIRres <- rbind(SIRres,c(ind*censusInterval,susceptiblenow,infectednow))
#         ind <- ind + 1
#       }
#     } else if(susceptiblenow!=0){
#       rate <- sum(h1t, h2t)
#       tau <- rexp(1, rate=rate)
#       
#       p <- runif(1); probs <- cumsum(c(h1t, h2t)/rate)
#       
#       if(p < probs[1]){
#         timenow <- timenow + tau; susceptiblenow <- susceptiblenow - 1; infectednow <- infectednow +1
#         h1t <- h1(a, b, susceptiblenow, infectednow); h2t <- h2(mu, infectednow);
#       } else {
#         timenow <- timenow + tau; susceptiblenow <- susceptiblenow ; infectednow <- infectednow - 1
#         h1t <- h1(a, b, susceptiblenow, infectednow); h2t <- h2(mu, infectednow);
#       }
#       if(SIRres[ind,1]<(ind*censusInterval) & timenow>(ind*censusInterval) & timenow<maxtime) {
#         SIRres <- rbind(SIRres,c(ind*censusInterval,susceptiblenow,infectednow))
#         ind <- ind + 1
#       }
#     }
#     if(infectednow==0 & SIRres[ind,3]!=0){
#       SIRres <- rbind(SIRres,c((ind*censusInterval),susceptiblenow,infectednow))
#     }
#   }
#   return(SIRres)
# }

# SIRsim2 depracated, does not work.
# SIRsim2 <- function(popsize, S0, I0, b, mu, a=0, tmax, censusInterval, prob, returnX = FALSE){
#     X <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize), observed = rep(0,2*popsize)))
#     
#     for(k in seq(1,(2*I0 - 1),by=2)){
#         X[k, 3] <- 1 #record an infection
#         X[k+1, 1] <- rexp(1, rate=mu) ; X[k+1, 3] <- -1  # record a recovery and the recovery time
#         X[k, 4] <- X[k+1, 4] <- ifelse(runif(1) <= prob, 1, 0) #record whether the infection and recovery were observed
#         
#     }
#     
#     X <- X[order(X[,1],X[,3]),]
#     timenow <- 0 ; susceptiblenow <- S0; infectednow <- I0; whichsusc <- (I0 + 1):popsize
#     
#     k <- I0 + 1; keep.going <- TRUE
#     while(keep.going==TRUE){
#         timeseq <- unique(X[,1]); irm <- buildirm(X, b, mu, a, pop=TRUE)
#         
#         for(j in 1:(length(timeseq)-1)){
#             eventtime <- timeseq[j] - log(1-runif(1))/irm[1,2,j]
#             if(eventtime < timeseq[j+1]) break
#         }
#         
#         if (eventtime>timeseq[length(timeseq)] | eventtime > tmax) {
#             eventtime <- Inf
#             keep.going <- FALSE # this condition would indicate that the epidemic has died off
#             break
#         } else if(eventtime <= timeseq[length(timeseq)]){
#             X[which(X[,2]==k)[1],c(1,3)] <- c(eventtime, 1) # record infection time and event code
#             X[which(X[,2]==k)[2],c(1,3)] <- c(eventtime - log(1-runif(1))/irm[2,3,1], -1) # record recovery time and event code
#             X[which(X[,2]==k), 4] <- ifelse(runif(1) <= prob, 1, 0) # record whether infection was observed
#             
#             X <- X[order(X[,1],X[,3]),]
#             k <- k+1
#             
#             if(k > popsize){
#                 keep.going <- FALSE  
#             } 
#         }
#     }
#         
#     SIRres <- data.frame(time = seq(0, tmax, by = censusInterval),
#                          infec.obs = 0,
#                          infec.true = 0)
#     for(k in 1:dim(SIRres)[1]){
#         SIRres$infec.obs[k] <- sum(X[which(X[,1] <= SIRres$time[k] & X[,4]==1),3])
#         SIRres$infec.true[k] <- sum(X[which(X[,1] <= SIRres$time[k]),3])
#         
#     }
#     
#     if(any(SIRres[,2]==0 & SIRres[,3]==0)){
#         SIRres <- SIRres[1:(min(which(SIRres[,2]==0 & SIRres[,3]==0))),]
#     }
#     
#     if(returnX == FALSE){
#         return(SIRres)
#     } else {
#         return(list(results = SIRres, trajectory = X))
#     }
# }

