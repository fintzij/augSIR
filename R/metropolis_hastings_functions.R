
# Metropolis-Hastings Ratio related functions (accept_prob, obs_prob, cacl_loglike, path_prob, pop_prob) ----------------------

# accept_prob calculates the acceptance ratio for the M-H algorithm
accept_prob <- function(pop_prob.new, pop_prob.cur, path_prob.cur, path_prob.new){
    if(pop_prob.new == -Inf & path_prob.new == -Inf){
        -Inf
    } else{
        pop_prob.new - pop_prob.cur + path_prob.cur - path_prob.new
    }
}

# accept_prob_param_MH calculates the MH acceptance ratio for a new set of parameters
accept_prob_param_MH <- function(log_lik.new, log_lik.cur, prior_probs.new, prior_probs.cur){
    log_lik.new - log_lik.cur + sum(prior_probs.new)- sum(prior_probs.cur)
}

# obs_prob calculates the binomial probability of observing the data given counts of the infecteds at observation times
obs_prob <- function(W, p){
    dbinom(x = sum(W[,2]), size = sum(W[,3]), prob = p, log = TRUE)
    
}

# calc_loglik calculates the complete data log-likelihood
calc_loglike <- function(Xcount, tmax, W,  b, m, a=0, p, initdist, popsize){
    indend <- nrow(Xcount)
        
    numinf <- Xcount[,2]
    numsusc <- Xcount[,3]
    
    events <- diff(Xcount[,2], lag = 1)
    
    infec.rates <- (b * numinf + a) * numsusc
    recov.rates <- m * numinf
    
    hazards <- infec.rates + recov.rates
    
    rates <- ifelse(events==1, infec.rates/numsusc, recov.rates/numinf)
    
    dbinom(sum(W[,2]), sum(W[,3]), prob=p, log=TRUE) + dmultinom(c(Xcount[1,3], Xcount[1,2], 0), prob = initdist, log=TRUE) + sum(log(rates[1:(indend - 1)])) - sum(hazards[1:(indend - 1)]*diff(Xcount[,1], lag = 1)) - hazards[indend]*max(0,tmax - Xcount[indend,1])
    
} 


# pop_prob and path_prob calculate the log-probabilities of the population trajectory and the subject trajectory for use in the M-H ratio

pop_prob <- function(Xcount, tmax, b, m, a = 0, initdist, popsize){
    indend <- dim(Xcount)[1]
    
    numinf <- Xcount[,2]
    numsusc <- Xcount[,3]
    
    events <- diff(Xcount[,2], lag = 1)
    
    infec.rates <- (b * numinf + a) * numsusc
    recov.rates <- m * numinf
  
    hazards <- infec.rates + recov.rates
    
    rates <- ifelse(events==1, infec.rates/numsusc, recov.rates/numinf)
    
    Xcount[1,2]*log(initdist[2]) + Xcount[1,3] * log(initdist[1]) + sum(log(rates[1:(indend - 1)])) - sum(hazards[1:(indend - 1)]*diff(Xcount[,1], lag = 1)) - hazards[indend]*max(0,tmax - Xcount[indend,1])

}

# path_prob <- function(path, fbmats, tpms, Xcount.other, obstimes){
#     
#     # set up vectors for discrete state sequences
#     obstime_seq <- rep(0, dim(fbmats)[3] + 1)
#     eventtime_seq <- rep(0, nrow(Xcount.other))
#     
#     # fill in vectors of discrete state sequences
#     if(all(path == 0)){ #no infection observed
#         obstime_seq[1:length(obstime_seq)] <- 1
#         eventtime_seq[1:length(eventtime_seq)] <- 1
#     } else {
#         obstime_seq[obstimes < path[1]] <- 1
#         obstime_seq[obstimes >= path[1] & obstimes < path[2]] <- 2
#         obstime_seq[obstimes >= path[2]] <- 3
#         
#         eventtime_seq[Xcount.other[,"time"] < path[1]] <- 1
#         eventtime_seq[Xcount.other[,"time"] < path[2] & Xcount.other[,"time"]>= path[1]] <- 2
#         eventtime_seq[Xcount.other[,"time"] >= path[2]] <- 3
#     }
#     
#     HMM_probs <- rep(0, length(obstime_seq))
#     eventtime_probs <- rep(0, length(eventtime_seq))
#     
#     # retrieve sequence of FB probabilities
#     HMM_probs[length(HMM_probs)] <- normalize(colSums(fbmats[,,dim(fbmats)[3]]))[obstime_seq[length(HMM_probs)]]
#     for(k in (length(HMM_probs)-1):1){
#         HMM_probs[k] <- normalize(fbmats[,obstime_seq[k+1],k])[obstime_seq[k]]
#     }
#     
#     # retrieve event time probability and discrete time skeleton probabilities if necessary
#     # first the case where the fb draw is that the state of the subject is unchanged at all observation times
#     if(all(obstime_seq == 1) || all(obstime_seq == 2)){
#         eventtime_probs[1:length(eventtime_probs)] <- 1
#     } else { # at least one transition is observed
#         changetimes <- which(diff(obstime_seq, lag = 1) != 0)
#         
#         if(length(changetimes) == 2){ 
#             
#         }
#     }
# }

path_prob <- function(path, Xcount.other, pathirm, initdist, tmax){
    indend <- dim(Xcount.other)[1]; init.infec <- Xcount.other[1,2]
    
    times <- Xcount.other[,1]; timediffs <- diff(times, lag = 1)
    numinf <- Xcount.other[,2]
    
    if(all(path==0)){ # no infection observed
        
        path.prob <- dmultinom(c(1,0,0), prob = initdist, log = TRUE)-sum(pathirm[1, 2, numinf[1:(indend - 1)] + 1] * timediffs)
        
    } else if(path[2] > 0 & path[2] < Inf) { # subject is infected and recovery is observed
        
        if(path[1]==0){ # subject is initially infected 
            
            path.prob <- dmultinom(c(0,1,0), prob = initdist, log = TRUE) + log(pathirm[2,3,1]) - pathirm[2,3,1]*path[2]
            
        } else if(path[1]!=0){ # subject is not initially infected
            
            if(times[2] > path[1]){ # no changes in number of other infecteds before subject's infection
                
                path.prob <- dmultinom(c(1,0,0), prob = initdist, log = TRUE) + log(pathirm[1,2, init.infec + 1]) - pathirm[1,2, init.infec+1]*path[1] + 
                    log(pathirm[2,3,1]) - pathirm[2,3,1]*(path[2] - path[1])
                
            } else if(times[2] < path[1]){ # if there is at least one change in the number of infecteds before the subject's infection
                ind1 <- sum(times < path[1])
                
                path.prob <- dmultinom(c(1,0,0), prob = initdist, log = TRUE) - # contribution from initial distribution
                    sum(pathirm[1,2, numinf[1:(ind1-1)] + 1] * timediffs[1:(ind1-1)]) + # contribution of escape probabilities before infection
                    log(pathirm[1,2,numinf[ind1] + 1]) - pathirm[1,2, numinf[ind1]+1]*(path[1] - times[ind1]) + # contribution of infection
                    log(pathirm[2,3,1]) - pathirm[2,3,1]*(path[2] - path[1]) # contribution of recovery
                
            }
        }
        
    } else if(path[2] == Inf) { # subject is infected and no recovery is observed
        
        if(path[1]==0){ # subject is initially infected 
            
            path.prob <- dmultinom(c(0,1,0), prob = initdist, log = TRUE) - pathirm[2,3,1]*tmax
            
        } else if(path[1]!=0){ # subject is not initially infected
            
            if(times[2] > path[1]){ # no changes in number of infecteds before subject's infection
                
                path.prob <- dmultinom(c(1,0,0), prob = initdist, log = TRUE) + log(pathirm[1,2, init.infec + 1]) - pathirm[1,2, init.infec + 1]*path[1] - pathirm[2,3,1]*(tmax - path[1])
                
            } else if(times[2] < path[1]){ # if there is at least one change in the number of infecteds before the subject's infection
                ind1 <- which(times < path[1])[sum(times < path[1])]
                
                path.prob <- dmultinom(c(1,0,0), prob = initdist, log = TRUE) - sum(pathirm[1,2, numinf[1:(ind1-1)] + 1] * timediffs[1:(ind1-1)]) + 
                    log(pathirm[1,2, numinf[ind1] + 1]) - pathirm[1,2, numinf[ind1] + 1]*(path[1] - times[ind1])  -
                    pathirm[2,3,1]*(tmax - path[1])
                
            }
        }
    } 
    
    return(path.prob)
}


# Functions to update parameters (update_rates, update_prob) ---------------------------------
update_rates <- function(Xcount, beta.prior, mu.prior, alpha.prior = NULL, initdist.prior = NULL, popsize){
    indend <- dim(Xcount)[1]
    infections <- diff(Xcount[,2], lag = 1)>0; recoveries <- !infections
    
    numsick <- Xcount[1:(indend - 1),2]; numsusc <- Xcount[1:(indend-1),3] 
    timediffs <- diff(Xcount[,1], lag = 1)


    beta.new <- rgamma(1, shape = (beta.prior[1] + sum(infections)), 
                       rate = beta.prior[2] + sum(numsick * numsusc * timediffs))
    
    mu.new <- rgamma(1, shape = (mu.prior[1] + sum(recoveries)), 
                     rate = mu.prior[2] + sum(numsick * timediffs))

    initprob.new <- rbeta(1, shape1 = (initdist.prior[1] + numsick[1]), shape2 = (initdist.prior[2] + popsize - numsick[1]))
    
    initdist.new <- c(1-initprob.new, initprob.new, 0)
    
    params.new <- c(beta.new, mu.new,0, initdist.new)
    
#     if(!is.null(alpha.prior)){
#         alpha.new <- rgamma(1, shape = (alpha.prior[1] + sum(infections)), 
#                             rate = alpha.prior[2] + sum(numsusc * timediffs))
#         
#         params.new[3] <- alpha.new
#     }
    
    return(params.new)
}

update_rates2 <- function(Xcount, beta.prior, mu.prior, alpha.prior = NULL, init.prior = NULL, popsize){
    indend <- dim(Xcount)[1]
    infections <- diff(Xcount[,2], lag = 1)>0; recoveries <- !infections
    
    numsick <- Xcount[1:(indend - 1),2]; numsusc <- Xcount[1:(indend-1),3] 
    timediffs <- diff(Xcount[,1], lag = 1)
    
    suff.stats <- c(numinfec = sum(infections), 
                    numrecov =  sum(recoveries), 
                    beta_suffstat = sum(numsick*numsusc*timediffs),
                    mu_suffstat = sum(numsick*timediffs))
    
    beta.new <- rgamma(1, shape = beta.prior[1]+suff.stats[1], rate = beta.prior[2] + suff.stats[3])
    mu.new <- rgamma(1, shape = mu.prior[1] + suff.stats[2], rate = mu.prior[2] + suff.stats[4])
    
    initprob.new <- rbeta(1, shape1 = (init.prior[1] + numsick[1]), shape2 = (init.prior[2] + popsize - numsick[1]))
    
    initdist.new <- c(1-initprob.new, initprob.new, 0)
    
    params.new <- c(beta.new, mu.new,0, initdist.new)
    
    if(!is.null(alpha.prior)){
        alpha.new <- rgamma(1, shape = (alpha.prior[1] + sum(infections)), 
                            rate = alpha.prior[2] + sum(numsusc * timediffs))
        
        params.new[3] <- alpha.new
    }
    
    return(list(params.new, suff.stats))
}

# update_prob updates the binomial sampling probability parameter
update_prob <- function(W, p.prior){
    rbeta(1,p.prior[1] + sum(W[,2]), p.prior[2] + sum(W[,3]-W[,2]))
    
}


# propose_params generates a vector of new parameters for a M-H update of params
# params is a vector of parameters
# kernel is a function that perturbes the parameters on the estimation scale, it should be symmetric
# toEst is a list of functions that transform the parameters to the estimation scale
# fromEst is a list of functions that transform the parameters from the estimation scale
# which_par is the index in the vector of parameters for which parameter to perturb
propose_params <- function(which_par, params, kernel, toEst, fromEst){
    
    params_est <- sapply(toEst, function(f) f(params))
    
    params_new_est <- kernel(params = params_est, which_par = which_par)
    
    params_new <- sapply(fromEst, function(f) f(params_new_est))
    
    return(list(params_new, params_new_est))
    
}

# evaluate_priors evaluates the prior density for a vector of parameter values
# params is the vector of parameters on the estimation scale
# d_priors is a list of functions to evaluate the prior densities
evaluate_priors <- function(params_est, d_priors){
    
    prior_dens <- mapply(do.call, d_priors, lapply(params_est, list))
    
    return(prior_dens)
    
}

update_prior_prob <- function(which_param, prior_probs, params_est, d_priors){
    
    prior_probs.new <- prior_probs
    
    prior_probs.new[which_param] <- d_priors[[which_param]](params_est[which_param])
    
    return(prior_probs.new)
    
}
