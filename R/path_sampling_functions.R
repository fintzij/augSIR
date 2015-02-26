
# drawtime draws an event time for either carriage acquisition or clearance. Xother is a matrix with the times of acquisition and clearance
# for all other subjects. irm is an array of instantaneous rate matrices. t0 and t1 are the left and right endpoints of the interval. 
# current state is the infection status at t0. 

drawtime <- function(Xcount, irm, t0, t1, currentstate){
    
    if(currentstate == 2){
        eventtime <- t0 - log(1 - runif(1)*(1-exp(-irm[2,3,1] * (t1 - t0))))/irm[2,3,1]
        
    } else if(currentstate ==1){
        times <- Xcount[,1]; numinf <- Xcount[,2]
        
        if(t1 != Inf){
            ind <- sum(times <= t0)
            
            if(numinf[ind] != 0){
                eventtime <- t0 - log(1 - runif(1)*(1 - exp(-irm[1, 2, numinf[ind] + 1]*(t1 - t0))))/irm[1, 2, numinf[ind] + 1]
                
            } else if(numinf[ind] == 0){
                eventtime <- Inf
                
            }
            
        } else if(t1 == Inf){
            
            timeseq <- c(t0, times[times > t0 & times < t1], t1)
            
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


# drawXt draws the path for a subject conditional on the other trajectories. Calls are first made to fwd and bwd in order to draw the infection status at observation times. 
# a call is then made to draw_subpath which samples the infection status at event times, and then samples the event times from the endpoint conditioned markov chain 
draw_path <- function(Xcount, irm, irm.eig, W, p, initdist, tmax){
    
    Xt <- cbind(W[,1], 0)
    
    # tpm_seq returns a list whose first element is the sequence of tpms between observation times, and whose second element is a list of tpm subsequences
    tpms <- tpm_seq(Xcount = Xcount, obstimes = W[,1], irm.eig = irm.eig)
    emits <- emit_seq(W = W, p = p)
    
    fbmats <- fwd(tpms = tpms[[1]], emits = emits, initdist = initdist)
    
    Xt[,2] <- bwd(mats = fbmats)
    
    path <- draw_eventtimes(Xt = Xt, Xcount = Xcount, tpm.seqs = tpms[[2]], irm = irm)    
    
    return(path)
}

# draw_eventtimes takes the matrix Xt and samples the event times appropriately determining whether to call draw_subseq or to draw directly. it returns a path.
draw_eventtimes <- function(Xt, Xcount, tpm.seqs, irm, tmax){
    
    path <- c(0,0)
    
    # first the case where the fb draw is that the subject is infected at all observation times
    if(all(Xt[,2] == 2)){
        
        path[1] <- 0
        
        t0 <- Xt[dim(Xt)[1], 1]; t1 <- Inf
        path[2] <- drawtime(Xcount, irm, t0, t1, 2)
        
        if(path[2] > tmax) path[2] <- Inf
        
        
        # next the case where the fb draw is that the subject is susceptible at all observation times
    } else if(all(Xt[,2] == 1)){
        
        t0 <- Xt[dim(Xt)[1],1]; t1 <- Inf
        path[1] <- drawtime(Xcount, irm, t0, t1, 1)
        
        if(path[1]>tmax) {
            
            path[1] <- path[2] <- Inf
            
        } else if(path[1] <= tmax){
            
            path[2] <- drawtime(Xcount, irm, path[1], t1, 2)
            
            if(path[2] > tmax) path[2] <- Inf
            
        }
        
    } else { # now the case where the fb draw shows that at least one change is observed    
        
        changetimes <- which(diff(Xt[,2], lag = 1) != 0)
        
        # first the case where all three states are observed
        
        if(length(changetimes) == 2){
            
            for(t in 1:2){
                
                tpm.seq <- tpm.seqs[[changetimes[t]]]
                
                t0 <- Xt[changetimes[t], 1]; t1 <- Xt[changetimes[t]+1, 1]
                timeseq <- c(t0, Xcount[Xcount[,1]>t0 & Xcount[,1] < t1,1], t1)
                
                init.state <- t; final.state <- t+1
                states <- draw_subseq(init.state = init.state, final.state = final.state, Xcount = Xcount, tpm.seq = tpm.seq, irm = irm)
                
                ind <- which(diff(states, lag = 1) != 0)
                
                path[t] <- drawtime(Xcount = Xcount, irm = irm, t0 = timeseq[ind], t1 = timeseq[ind + 1], t)
                
            }
            
        } else if(length(changetimes) == 1){
            
            # get the tpm subseq
            tpm.seq <- tpm.seqs[[changetimes]]
            
            # get the sequence of event times, bookended by observation times
            t0 <- Xt[changetimes,1]; t1 <- Xt[changetimes+1,1]
            timeseq <- c(t0, Xcount[Xcount[,1]>t0 & Xcount[,1] < t1,1], t1)
            
            
            if((Xt[changetimes+1,2] - Xt[changetimes,2]) == 2){ ### Subject transitions from healthy to recovered within one observation period
                
                # set initial and final states
                init.state <- 1; final.state <- 3
                
                # draw the states at event times
                states <- draw_subseq(init.state = init.state, final.state = final.state, Xcount = Xcount, tpm.seq = tpm.seq, irm = irm)
                
                # get the indices for when states change
                statediffs <- which(diff(states, lag = 1) != 0)
                
                
                if(length(statediffs) == 2){
                    
                    path[1] <- drawtime(Xcount = Xcount, irm = irm, t0 = timeseq[statediffs[1]], t1 = timeseq[statediffs[1] + 1], 1)
                    
                    path[2] <- drawtime(Xcount = Xcount, irm = irm, t0 = timeseq[statediffs[2]], t1 = timeseq[statediffs[2] + 1], 2)
                    
                    
                } else if(length(statediffs == 1)){
                    
                    ############# This is the case where the process jumps twice in a single interval. Wait for Vladimir's response
                    
                    
                    
                }
                
            } else if(Xt[changetimes,2]==1 & ((Xt[changetimes+1,2] - Xt[changetimes,2]) == 1)){ # healthy to infected within one observation period
                
                # set initial and final states
                init.state <- 1; final.state <- 2
                
                # draw the subsequence of states
                states <- draw_subseq(init.state = init.state, final.state = final.state, Xcount = Xcount, tpm.seq = tpm.seq, irm = irm)
                
                # find the index for when the states transition
                statediffs <- which(diff(states, lag = 1) != 0)               
                
                # sample the event time
                path[1] <- drawtime(Xcount = Xcount, irm = irm, t0 = timeseq[statediffs], t1 = (timeseq[statediffs]+1), 1)
                
                # set t0 and t1 and draw the recovery time
                t0 <- Xt[dim(Xt)[1],1]; t1 <- Inf
                path[2] <- drawtime(Xcount, irm, t0, t1, 2)
                
                if(path[2] > tmax) path[2] <- Inf
                
                
            } else if(Xt[changetimes,2]==2 & ((Xt[changetimes+1,2] - Xt[changetimes,2]) ==1)){ #infected to recovered within one observation period. note that the subject is initially infected
                
                # set infection time to 0
                path[1] <- 0
                
                # set initial and final states
                init.state <- 2; final.state <- 3
                
                # draw the subsequence of states
                states <- draw_subseq(init.state = init.state, final.state = final.state, Xcount = Xcount, tpm.seq = tpm.seq, irm = irm)
                
                # find the index for when the states transition
                statediffs <- which(diff(states, lag = 1) != 0)               
                
                # sample the event time
                path[2] <- drawtime(Xcount = Xcount, irm = irm, t0 = timeseq[statediffs], t1 = (timeseq[statediffs]+1), 2)
                
            }
            
        }
        
    }   
    return(path)    
}

# draw_subseq takes a matrix with the infection status at observation times and returns a sample path by calling draw_times after sampling the infection status at event times
draw_subseq <- function(init.state, final.state, Xcount, tpm.seq, irm){
    eventtimes <- Xcount[,1]    
    
    tpm.subseq <- tpm.seq[,,,1] # the array of tpms for the sequence of event times
    tpm.products <- tpm.seq[,,,2] # the array of products of tpms for the event times
    
    states <- c(init.state, rep(final.state, dim(tpm.subseq)[3]))
    
    if(length(states) >= 3){
        state.probs <- tpm.subseq[init.state, , 1] * tpm.products[ , final.state, 2] / tpm.products[init.state, final.state, 1]
        
        if(any(is.nan(state.probs))) state.probs <- replace(state.probs, is.nan(state.probs), 0)
        
        states[2] <- sample.int(n=3, size = 1, prob = state.probs)
        
        if(length(states) > 3){
            
            for(s in 3:(length(states) - 1)){
                state.probs <- tpm.subseq[states[s-1], , s-1] * tpm.products[ , final.state, s] / tpm.products[states[s-1], final.state, s-1]
                
                if(any(is.nan(state.probs))) state.probs <- replace(state.probs, is.nan(state.probs), 0)
                
                states[s] <- sample.int(n = 3, size = 1, prob = state.probs)
                
            }
            
        }
        
    }
    
    return(states)       
    
}
