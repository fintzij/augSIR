# rjmcmc_draw proposes an insertion, deletion, or shift to the current trajectory
rjmcmc_draw <- function(path.cur, Xcount.cur, j, initdist, shift.int, insert.prob, remove.prob, shift.prob, tmax, b, m , p){
    
    path.new <- path.cur
    
    t_epideath <- ifelse(max(Xcount.cur[,1]) < tmax & Xcount.cur[nrow(Xcount.cur), 2] == 0, max(Xcount.cur[,1]), tmax) # time of epidemic death, if the epidemic is observed to die off    
    
    if(all(path.cur[,3] == 0)){ # neither an infection, nor a recovery are observed in the current path
        
        
        # choose whether to insert both an infection, or an infection and a recovery
        which.insert <- sample.int(2, 1)
        
        if(which.insert == 1){ # with probability 0.5, insert an infection only
            
            # choose whether the subject is initially infected
            init.infec <- sample.int(3, 1, prob = initdist)
            
            if(init.infec == 2){ # with probability initdist[2] subject is initially infected
                
                path.new[1, c(1,3)] <- c(0,1)
                                                
                
            } else if(init.infec == 1){ # with probability 1 - initdist[2] subject is initially susceptible
                
                path.new[1, c(1,3)] <- c(runif(1, 0, t_epideath), 1)
                                
            }
            
            
        } else if (which.insert == 2){ # with probability 0.5 insert both an infection and a recovery
            
            # choose whether the subject is initially infected
            init.infec <- sample.int(3, 1, prob = initdist)
            
            if(init.infec == 2){ # with probability initdist[2] subject is initially infected
                
                path.new[1, c(1,3)] <- c(0,1)
                path.new[2, c(1,3)] <- c(runif(1, 0, tmax), -1)
                
                
            } else if(init.infec == 1){ # with probability 1 - initdist[2] subject is initially susceptible
                
                path.new[1, c(1,3)] <- c(runif(1, 0, t_epideath), 1)
                path.new[2, c(1,3)] <- c(runif(1, path.new[1,1], tmax), -1)
                
            }
            
        }
        
    } else if(all(path.cur[,3] == c(1, 0))){ # an infection is observed in the current path but no recovery
        
        event <- sample.int(3, 1, prob = c(remove.prob, insert.prob, shift.prob)) # select whether a removal, insertion, or shift occurs
        
        if(event == 1) { # remove the infection
            
            path.new[1, c(1,3)] <- 0 # remove the infection
                        
            
        } else if(event == 2) { # insert the recovery
            
            path.new[2, c(1,3)] <- c(runif(1, path.new[1, 1], tmax), -1)
                        
            
        } else if(event == 3) { # shift the infection uniformly within an interval around the current value or to zero if that shift is available
            
            if(path.cur[1,1] == 0){ # if the infection time is 0, we can only shift away from 0
                
                dt <- runif(1, -shift.int, shift.int) 
                                
                path.new[1, c(1, 3)] <- c(path.cur[1,1] + dt , 1)
                
                
            } else if(path.cur[1,1] != 0){
                
                if(0 > (path.cur[1,1] - shift.int)){ # it is possible to shift to 0, so select whether to do so
                                        
                    infec_at_0 <- sample.int(3, 1, prob = initdist) 
                    
                    if(infec_at_0 == 2){# subject is infected at time 0
                        
                        path.new[1, c(1,3)] <- c(0, 1)
                        
                    } else if(infec_at_0 == 1){#subject is susceptible at 0, so draw as usual
                        dt <- runif(1, -shift.int, shift.int)
                        
                        path.new[1, c(1,3)] <- c(path.cur[1,1] + dt, 1)
                        
                    }
                    
                } else if(0 < (path.cur[1,1] - shift.int)){ # draw a new time as usual
                    dt <- runif(1, -shift.int, shift.int)
                    
                    path.new[1, c(1,3)] <- c(path.cur[1,1] + dt, 1)
                    
                }

            }           
            
        }
        
    } else if(all(path.cur[,3] != 0)) {# both an infection or recovery are observed
        
        event <- sample.int(2, 1, prob = c(remove.prob, shift.prob)) # select whether to remove just the recovery, the recovery and the infection, or to shift one of the recovery or infection
        
        event <- ifelse(event == 1, sample.int(2,1), 3)
        
        if(event == 1){# remove the recovery
            
            path.new[2, c(1,3)] <- 0 # remove the recovery
            
            
        } else if(event == 2){ # remove both the recovery and the infection
            
            path.new[, c(1,3)] <- 0
            
            
        } else if(event == 3){ # a shift of either the infection or recovery occurs  
            
            which.shift <- sample.int(2,1) # select whether to shift the infection or the recovery
            
            if(which.shift == 1){ # shift the infection
                
                if(path.cur[1,1] == 0){ # if the infection time is 0, we can only shift away
                    
                    dt <- runif(1, -shift.int, shift.int)
                    
                    path.new[1, c(1, 3)] <- c(path.cur[1,1] + dt , 1)
                    
                    
                } else if(path.cur[1,1] != 0){
                    
                    if(0 > (path.cur[1,1] - shift.int)){ # it is possible to shift to 0, so select whether to do so
                        
                        infec_at_0 <- sample.int(3, 1, prob = initdist) # sample when 
                        
                        if(infec_at_0 == 2){# subject is infected at time 0
                            
                            path.new[1, c(1,3)] <- c(0, 1)
                            
                        } else if(infec_at_0 == 1){#subject is susceptible at 0, so draw as usual
                            dt <- runif(1, -shift.int, shift.int)
                            
                            path.new[1, c(1,3)] <- c(path.cur[1,1] + dt, 1)
                            
                        }
                        
                    } else if(0 < (path.cur[1,1] - shift.int)){ # draw a new time as usual
                        dt <- runif(1, -shift.int, shift.int)
                        
                        path.new[1, c(1,3)] <- c(path.cur[1,1] + dt, 1)
                        
                    }
                    
                }   
                
            } else if(which.shift == 2){ # shift the recovery
                
                dt <- runif(1, -shift.int, shift.int)
                
                path.new[2, c(1, 3)] <- c(path.cur[2, 1] + dt, -1)
                
            }
            
        }
        
        
    }
    
    return(path.new)
    
}

# rjmcmc_ratio computes the acceptance ratio for a proposed move
rjmcmc_ratio <- function(W.cur, W.new, X.cur, X.new, Xcount.cur, Xcount.new, path.cur, path.new, initdist, shift.int, insert.prob, remove.prob, shift.prob, samp_prob, tmax, popsize){
    
    # calculate the time of epidemic death (equal to tmax if the epidemic doesn't die off in the observation window)
    
    t_epideath <- ifelse(max(Xcount.cur[,1]) < tmax & Xcount.cur[nrow(Xcount.cur), 2] == 0, max(Xcount.cur[,1]), tmax)
    
    # determine whether the new path should be automatically rejected
    if(path.new[1,1] > path.new[1,2] | # infection is proposed after recovery
           path.new[1,1] < 0 | # infection is proposed before time 0
           path.new[1,1] > t_epideath | # infection is proposed after the time of epidemic death
           path.new[2,1] > tmax | # recovery is proposed after tmax
           path.new[2,1] < 0){ # recovery is proposed before time 0
        auto_reject <- TRUE
        
    } else{
        auto_reject <- FALSE
        
    }
    
    if(auto_reject == FALSE){
        
        loglike_cur <- dbinom(sum(W.cur[,2]), sum(W.cur[,3]), prob = samp_prob, log = TRUE) + pop_prob(Xcount = Xcount.cur, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
        
        loglike_new <- dbinom(sum(W.new[,2]), sum(W.new[,3]), prob = samp_prob, log = TRUE) + pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
        
        # calculate acceptance ratio for the appropriate move type
        if(all(path.cur[,3] == c(0,0)) & all(path.new[,3] == c(1,0))){ # move from 0 -> 1
            
            rjmcmc.ratio <- min(0, loglike_new - loglike_cur + log(remove.prob) - log(insert.prob) - log(0.5) - log(ifelse(path.new[1,1] == 0, initdist[2], (initdist[1]/t_epideath))))
            
            
        } else if(all(path.cur[,3] == c(0,0)) & all(path.new[,3] == c(1,-1))){ # move from 0 -> 2
            
            rjmcmc.ratio <- min(0, loglike_new - loglike_cur + log(remove.prob) - log(insert.prob) - log(ifelse(path.new[1,1] == 0, initdist[2], (initdist[1]/t_epideath))) - log(1/(tmax - path.new[1,1]))) 
            
            
        } else if(all(path.cur[,3] == c(1,0)) & all(path.new[,3] == c(0,0))){ # move from 1 -> 0
            
            rjmcmc.ratio <- min(0, loglike_new - loglike_cur + log(insert.prob) - log(remove.prob) + log(0.5) + log(ifelse(path.cur[1,1] == 0, initdist[2], initdist[1]/t_epideath)))
            
            
        } else if(all(path.cur[,3] == c(1, 0)) & all(path.new[,3] == c(1, 0))){ # move from 1 -> 1
            
            if(path.new[1,1] == 0){ # new infection time is 0
                
                rjmcmc.ratio <- min(0, loglike_new - loglike_cur - log(initdist[2]))
                
                
            } else if(path.cur[1,1] == 0){ # current infection time is 0
                
                rjmcmc.ratio <- min(0, loglike_new - loglike_cur + log(initdist[2]))
                
                
            } else if(path.new[1,1] != 0 & path.cur[1,1] != 0){
                
                rjmcmc.ratio <- min(0, loglike_new - loglike_cur + log(ifelse(0 > (path.new[1,1] - shift.int), initdist[1], 1)) - log(ifelse(0 > (path.cur[1,1]), initdist[1], 1))) # note that the 1/(2*shift.int) appears in both the numerator and denominator and thus cancels
                
            }
            
            
        } else if(all(path.cur[,3] == c(1, 0)) & all(path.new[,3] == c(1, -1))){ # move from 1 -> 2
            
            rjmcmc.ratio <- min(0, loglike_new - loglike_cur + log(remove.prob) - log(insert.prob) + log(0.5) - log(1/(tmax - path.cur[1,1])))
            
            
        } else if(all(path.cur[,3] == c(1, -1)) & all(path.new[,3] == c(0, 0))){ # move from 2 -> 0
            
            rjmcmc.ratio <- min(0, loglike_new - loglike_cur + log(insert.prob) - log(remove.prob) + log(ifelse(path.cur[1,1]==0, initdist[2], initdist[1]/t_epideath)) + log(1/(tmax - path.cur[1,1])))
            
            
        } else if(all(path.cur[,3] == c(1, -1)) & all(path.new[,3] == c(1, 0))){ # move from 2 -> 1
            
            rjmcmc.ratio <- min(0, loglike_new - loglike_cur + log(insert.prob) - log(remove.prob) - log(0.5) + log(1/(tmax - path.cur[1,1])))
            
            
        } else if(all(path.cur[,3] == c(1, -1)) & all(path.new[,3] == c(1, -1))){ # move from 2 -> 2
            
            if(path.new[1,1] != path.cur[1,1]){ # the infection is shifted
                if(path.new[1,1] == 0){ # new infection time is 0
                    
                    rjmcmc.ratio <- min(0, loglike_new - loglike_cur - log(initdist[2]))
                    
                    
                } else if(path.cur[1,1] == 0){ # current infection time is 0
                    
                    rjmcmc.ratio <- min(0, loglike_new - loglike_cur + log(initdist[2]))
                    
                    
                } else if(path.new[1,1] != 0 & path.cur[1,1] != 0){
                    
                    rjmcmc.ratio <- min(0, loglike_new - loglike_cur + log(ifelse(0 > (path.new[1,1] - shift.int), initdist[1], 1)) - log(ifelse(0 > (path.cur[1,1]), initdist[1], 1)))
                    
                }
                
            } else if(path.new[2,1] != path.cur[2,1]){ # the recovery is shifted
                rjmcmc.ratio <- min(0, loglike_new - loglike.cur)                
                
            }   
            
        }
        
    } else if(auto_reject == TRUE){
        
        rjmcmc.ratio <- -Inf
        
    }    
    
    
    return(rjmcmc.ratio)    
    
}
