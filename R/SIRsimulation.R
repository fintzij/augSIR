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


