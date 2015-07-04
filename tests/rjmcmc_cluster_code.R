# set working directory and load augSIR files
source("augSIR.R")
source("auxilliary_functions.R")
source("SIRsimulation.R")
source("rjmcmc_functions.R")
source("matrix_build_update.R")
source("metropolis_hastings_functions.R")
source("path_sampling_functions.R")

args <- commandArgs(TRUE)
print(args)

sim_number <- as.numeric(args[1])

popsize = 4; tmax = 4
b <- 0.5 + runif(1,-0.0001, 0.0001)
m <- 1 + runif(1, -0.0001, 0.0001)
samp_prob <- 0.5 
initdist <- c(0.7, 0.3, 0)
insert.prob = 1/3; remove.prob = 1/3; shift.prob = 1/3
shift.int <- 0.1
niter <- 15000
samp_size = 5

keep.going = TRUE

rjmcmc_jointcomp_results <- vector("list", length = samp_size)

for(k in 1:samp_size) {
    
    if(keep.going == TRUE){
        
        print(k)
        
    SIRres<-SIRsim(popsize = 4, initdist = c(0.7, 0.3, 0), b = 0.5, mu=1, a=0, tmax = 4, censusInterval=0.05, sampprob = samp_prob, returnX = TRUE, trim = FALSE)
    
    # make sure the epidemic is interesting (i.e. at least one infection)
    if(max(SIRres$results[,3]) < 2){
        
        while(max(SIRres$results[,3]) < 2){
            
            SIRres<-SIRsim(popsize = 4, initdist = c(0.7, 0.3, 0), b = 0.5, mu=1, a=0, tmax = 4, censusInterval=0.05, sampprob = samp_prob, returnX = TRUE, trim = FALSE)
            
        }
    }
    
    # observation matrix
    W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Truth, augmented = 0))
    
    # individual trajectories
    X.cur <- SIRres$trajectory
    
    # count matrix
    Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
    
    # update observation matrix
    W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
    
    subjects <- rep(1:popsize, each = niter)
    
    
    for(j in 1:length(subjects)){
        
        if(keep.going == TRUE){
        
        path.cur <- X.cur[X.cur[,2] == subjects[j], ]
        path.new <- rjmcmc_draw(path.cur = path.cur, Xcount.cur, j = subjects[j], initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, tmax = max(W.cur[,1]), b = b, m = m, p = samp_prob)
        
        # update bookkeeping objects
        X.new <- X.cur; X.new[X.new[,2] == subjects[j], ] <- path.new
        Xcount.new <- build_countmat(X = X.new, popsize = popsize)
        W.new <- updateW(W = W.cur, Xcount = Xcount.new)
        
        # calculate acceptance ratio
        rjmcmc.ratio <- rjmcmc_ratio(W.cur = W.cur, W.new = W.new, X.cur = X.cur, X.new = X.new, Xcount.cur = Xcount.cur, Xcount.new = Xcount.new, path.cur = path.cur, path.new = path.new, initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, b = b, m = m, samp_prob = samp_prob, tmax = tmax, popsize = popsize)
        if(missing(rjmcmc.ratio)){
            keep.going <- FALSE
            break()           
        }
        # decide whether to accept. if accept, update current trajectories
        if(rjmcmc.ratio > log(runif(1))){
            
            X.cur <- X.new
            Xcount.cur <- Xcount.new
            W.cur <- W.new
            
        }
        }
        
    }
    
    rjmcmc_jointcomp_results[[k]] <- list(Truth = SIRres$trajectory, Observed = SIRres$results, Augmented = Xcount.cur)  
    }
}

assign(paste("rjmcmc_jointcomp_results",sim_number,sep=""),rjmcmc_jointcomp_results)

save(list = paste("rjmcmc_jointcomp_results",sim_number,sep=""), file = paste("augSIR_jointcomp_results_test",sim_number, ".Rdata", sep=""))
