
# Simulation #1 for fixed dataset with no resampling at observation times ---------------

# Procedure
# 1) Simulate dataset using Gillespie - binomial samples taken in increments of 0.5 (coarser grid)
# 2) Remove one individual (say subject 1)
# 3) For k in 1:niter
#    a) Simulate that individual forward (exponential waiting time)
#    b) Redraw the path for the individual, keeping the infection status at observation times from the first step. Note that this step does not involve data, so there is no need to draw a new binomial sample.  


library(doParallel)

for(j in 1:5){
print(j)
# set simulation parameters
niter <- 100000; popsize = 4; tmax = 10
b <- 0.5 + runif(1, -0.0001, 0.0001)
m <- 1 + runif(1, -0.0001, 0.0001)
samp_prob <- 0.5 
initdist <- c(0, 1, 0)
samp_size = 3
obstimes <- seq(0, tmax, by=0.5)

# simlate dataset
SIRres<-SIRsim(popsize = 4, initdist = initdist, b = b, mu=m, a=0, tmax = tmax, censusInterval=0.5, sampprob = 0.5, returnX = TRUE, trim = FALSE)

# make sure the epidemic is interesting (i.e. at least one infection event)
if(max(SIRres$results[,3]) < 2 | SIRres$results[1,2] == 0){
    
    while(max(SIRres$results[,3]) < 2 | SIRres$results[1,2] == 0){
        
        SIRres<-SIRsim(popsize = 4, initdist = initdist, b = b, mu=m, a=0, tmax = tmax, censusInterval=0.5, sampprob = 0.5, returnX = TRUE, trim = FALSE)
        
    }
}

# initialize bookkeeping objects
# observation matrix
W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Observed, augmented = 0))

# individual trajectories
X.cur <- SIRres$trajectory

# count matrix
Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)

# update observation matrix
W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
Xother <- X.cur[X.cur[,2]!=1,]

path.cur <- getpath(X.cur, 1)

# .other objects
W.other <-get_W_other(W.cur = W.cur, path = path.cur)
Xcount.other <- build_countmat(Xother, popsize - 1)

pathirm.cur <- build_irm(Xcount = Xcount.cur, b = b, m = m, a = 0, popsize = popsize, pop = FALSE)
patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)


cl <- makeCluster(4)
registerDoParallel(cl)

strt <- Sys.time()
writeLines(c(""), "log.txt")

sim_one_fixpop_noHMM <- foreach(iter = 1:niter, .packages = "augSIR") %dopar% {
    
    if(iter %% 1000 == 0){
        sink("log.txt", append=TRUE)  
        cat(paste("Starting iteration",iter,"\n"))  
        sink()
    }
    
    # draw new path
    # first via gillespie
    gillespie.path <- c(0, rexp(1,rate = m))
    gillespie.path <- ifelse(gillespie.path > tmax, Inf, gillespie.path)
    
    gillespie.status <- ifelse(obstimes < gillespie.path[2], 2, 3)
    
    Xt <- cbind(obstimes, gillespie.status)
    
    tpms <- tpm_seq(Xcount = Xcount.other, obstimes = obstimes, irm.eig = patheigen.cur)
    
    path.new <- draw_eventtimes(Xt = Xt, Xcount = Xcount.other, tpm.seqs = tpms[[2]], irm = pathirm.cur, tmax = tmax) 
    
    c(ifelse(seq(0, tmax, by =0.05) < gillespie.path[2], 1, 0), ifelse(seq(0, tmax, by =0.05) < path.new[2], 1, 0)) 
}
stopCluster(cl)

print(Sys.time() - strt)

# save(sim_one_fixpop_noHMM, file = paste("sim_one_fixpop_noHMM_",j,".Rdata",sep=""))
# Plots -------------------------------------------------------------------


# plot average infection status at each time point

statuses <- matrix(0, nrow = (2*length(seq(0,tmax, by=0.05))), ncol = length(sim_one_fixpop_noHMM))

for(k in 1:ncol(statuses)){
    
    statuses[,k] <- sim_one_fixpop_noHMM[[k]]
    
}

means <- rowMeans(statuses)

# plots
truepath <- getpath(SIRres$trajectory, 1)

path_comp <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))

path_comp2 <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))

path_comp2[,2] <- path_comp2[,2] - path_comp2[1:(nrow(path_comp[2])/2), 2]


pdf(file = paste("Sim_one_fixpop_noHMM_sim",j,".pdf",sep=""))

print(ggplot(SIRres$results[,c(1,3)], aes(x=time, y=Truth)) + geom_step() + theme_bw() + labs(title = paste("True recovery for subject j = ",signif(truepath[2], digits = 4), sep="")))


print(ggplot(path_comp, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual. No HMM."))

# subtract off gillespie means

print(ggplot(path_comp2, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual. No HMM. Subtracting Gillespie means."))

dev.off()
# Simulation #2 - for fixed dataset --------------------------------------------

# Procedure
# 1) Simulate dataset using Gillespie
# 2) Remove one individual (say subject 1)
# 3) For k in 1:niter
#    a) Simulate that individual forward (exponential waiting time)
#    b) Draw new binomial sample
#    c) Use the binomial sample to redrawthe path for the individual via our method

library(doParallel)

# set simulation parameters
niter <- 100000; popsize = 4; tmax = 10
b <- 0.5 + runif(1, -0.0001, 0.0001)
m <- 1 + runif(1, -0.0001, 0.0001)
samp_prob <- 0.5 
initdist <- c(0, 1, 0)
samp_size = 3
obstimes <- seq(0, tmax, by=0.05)

# simlate dataset
SIRres<-SIRsim(popsize = 4, initdist = initdist, b = b, mu=m, a=0, tmax = tmax, censusInterval=0.05, sampprob = 0.5, returnX = TRUE, trim = FALSE)

# make sure the epidemic is interesting (i.e. at least one infection event)
if(max(SIRres$results[,3]) < 2 | SIRres$results[1,2] == 0){
    
    while(max(SIRres$results[,3]) < 2 | SIRres$results[1,2] == 0){
        
        SIRres<-SIRsim(popsize = 4, initdist = initdist, b = b, mu=m, a=0, tmax = tmax, censusInterval=0.05, sampprob = 0.5, returnX = TRUE, trim = FALSE)
        
    }
}

# initialize bookkeeping objects
# observation matrix
W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Observed, augmented = 0))

# individual trajectories
X.cur <- SIRres$trajectory

# count matrix
Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)

# update observation matrix
W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
Xother <- X.cur[X.cur[,2]!=1,]

path.cur <- getpath(X.cur, 1)

# .other objects
W.other <-get_W_other(W.cur = W.cur, path = path.cur)
Xcount.other <- build_countmat(Xother, popsize - 1)

pathirm.cur <- build_irm(Xcount = Xcount.cur, b = b, m = m, a = 0, popsize = popsize, pop = FALSE)
patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)


cl <- makeCluster(3)
registerDoParallel(cl)

strt <- Sys.time()
writeLines(c(""), "log2.txt")

sim_one_fixpop <- foreach(iter = 1:niter, .packages = "augSIR") %dopar% {
    
    if(iter %% 1000 == 0){
        sink("log2.txt", append=TRUE)  
        cat(paste("Starting iteration",iter,"\n"))  
        sink()
    }
    
    # draw new path
    # first via gillespie
    gillespie.path <- c(0, rexp(1,rate = m))
    gillespie.path <- ifelse(gillespie.path > tmax, Inf, gillespie.path)
    
    # update bookkeeping objects
    Xcount.cur <- update_Xcount(Xcount.other, gillespie.path)

    W.cur <- updateW(W.other, path = gillespie.path)
    W.cur[,2] <- rbinom(nrow(W.cur), W.cur[,3], samp_prob)
    W.other[,2] <- W.cur[,2]
    
    path.new <- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = samp_prob, initdist = initdist, tmax = tmax)
    
    W.new <- updateW(W.other, path = path.new)
    Xcount.new <- update_Xcount(Xcount.other, path.new)
    
    path_prob.cur <- path_prob(gillespie.path, Xcount.other, pathirm.cur, initdist, tmax)
    path_prob.new <- path_prob(path.new, Xcount.other, pathirm.cur, initdist, tmax)
    pop_prob.cur <- pop_prob(Xcount.cur, tmax, b, m, a = 0, initdist, popsize)
    pop_prob.new <- pop_prob(Xcount.new, tmax, b, m, a = 0, initdist, popsize)
    
    a.prob <- accept_prob(pop_prob.new, pop_prob.cur, path_prob.cur, path_prob.new)
    
    if(a.prob >=0 || exp(a.prob)> runif(1)){
        path.new <- path.new
    } else{
        path.new <- gillespie.path
        sink("log2.txt", append = TRUE)
        cat(paste("Path is rejected!"))
        sink()
    }
    
    c(ifelse(obstimes < gillespie.path[2], 1, 0), ifelse(obstimes < path.new[2], 1, 0)) 
}
stopCluster(cl)

print(Sys.time() - strt)

# save(sim_one_fixpop_noHMM, file = paste("sim_one_fixpop_",j,".Rdata",sep=""))


# Plots -------------------------------------------------------------------

statuses <- matrix(0, nrow = (2*length(seq(0,tmax, by=0.05))), ncol = length(sim_one_fixpop))

for(k in 1:ncol(statuses)){
    
    statuses[,k] <- sim_one_fixpop[[k]]
    
}

means <- rowMeans(statuses)

pdf(file = paste("Sim_one_fixpop_sim",j,".pdf",sep=""))

truepath <- getpath(SIRres$trajectory, 1)
print(ggplot(SIRres$results[,c(1,3)], aes(x=time, y=Truth)) + geom_step() + theme_bw() + labs(title = paste("True recovery for subject 1 = ",signif(truepath[2], digits = 4), sep="")))

path_comp <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))

path_comp2 <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))

path_comp2[,2] <- path_comp2[,2] - path_comp2[1:(nrow(path_comp[2])/2), 2]


print(ggplot(path_comp, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual"))
print(ggplot(path_comp2, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual. No HMM. Subtracting Gillespie means."))

dev.off()

}
# Simulation #3 resimulating many different datasets -------------------------

# Procedure: (1) Simulate new population trajectory using Gillespie
#            (2) Draw binomial samples at observation times
#            (3) remove subject 1, then redraw his trajectory using our method 


library(doParallel)

niter <- 50; popsize = 4; tmax = 10
b <- 0.5 + runif(1, -0.0001, 0.0001)
m <- 1 + runif(1, -0.0001, 0.0001)
samp_prob <- 0.5 
initdist <- c(0, 1, 0)
samp_size = 1
obstimes <- seq(0, tmax, by = 0.05)

pathirm.cur <- build_irm(Xcount = Xcount.cur, b = b, m = m, a = 0, popsize = popsize, pop = FALSE)
patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)

cl <- makeCluster(4)
registerDoParallel(cl)

strt <- Sys.time()
writeLines(c(""), "log2.txt")

sim_one_results <- foreach(iter = 1:samp_size, .packages = "augSIR") %do% {
    
    if(iter %% 1000 == 0){
        sink("log2.txt", append=TRUE)  
        cat(paste("Starting iteration",iter,"\n"))  
        sink()
    }
    
    SIRres<-SIRsim(popsize = 4, initdist = initdist, b = b, mu=m, a=0, tmax = tmax, censusInterval=0.05, sampprob = 0.5, returnX = TRUE, trim = FALSE)
    
    # make sure the epidemic is interesting (i.e. at least one infection)
    if(max(SIRres$results[,3]) < 2 | SIRres$results[1,2] == 0){
        
        while(max(SIRres$results[,3]) < 2 | SIRres$results[1,2] == 0){
            
            SIRres<-SIRsim(popsize = 4, initdist = initdist, b = b, mu=m, a=0, tmax = tmax, censusInterval=0.05, sampprob = 0.5, returnX = TRUE, trim = FALSE)
            
        }
    }
    
    # initialize objects
    #     accepts <- rep(0, length(niter)); accepts[1] <- 1
    #     accept_probs <- rep(0, length(niter)); accept_probs[1] <- 0
    #     log_likelihood <- rep(0, length(niter))
    #     trajectories <- vector("list", length(niter))
    #     paths <- vector("list", length(niter))
    
    # observation matrix
    W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Observed, augmented = 0))
    
    # individual trajectories
    X.cur <- SIRres$trajectory
    
    # count matrix
    Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
    #     trajectories[[1]] <- Xcount.cur
    
    # update observation matrix
    W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
    
    pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
    
    # build irm matrices
    #     pathirm.cur <- build_irm(Xcount = Xcount.cur, b = b, m = m, a = 0, popsize = popsize, pop = FALSE)
    #     patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
    
    subjects <- rep(1, niter)
    
    acceptances <- rep(0, length(subjects))
    
    Xother <- X.cur[X.cur[,2]!=subjects[j],]
    Xcount.other <- build_countmat(Xother, popsize - 1)
        
    
    for(j in 1:length(subjects)){
        
#         Xother <- X.cur[X.cur[,2]!=subjects[j],]
        
        path.cur <- getpath(X.cur, subjects[j])
        
        if(!identical(path.cur, Xcount.cur[,1])){
            
#             Xcount.other <- build_countmat(Xother, popsize - 1)
            
            if(nrow(Xcount.other)!=1){
                
                W.other <-get_W_other(W.cur = W.cur, path = path.cur) 
                
                path.new<- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = samp_prob, initdist = initdist, tmax = tmax)
                
                X.new <- updateX(X = X.cur, path = path.new, j = subjects[j])
                Xcount.new <- update_Xcount(Xcount.other = Xcount.other, path = path.new)
                
                W.new <- updateW(W = W.other, path = path.new)
                
                if(max(Xcount.new[,2]) == pathirm.cur[4,4,dim(pathirm.cur)[3]]){
                    
                    new.numinf <- pathirm.cur[4,4,dim(pathirm.cur)[3]]+1
                    
                    pathirm.cur <- update_irm(irm = pathirm.cur, new.numinf = new.numinf, b = b, m = m, a = 0, popsize = popsize)
                    patheigen.cur <- update_eigen(patheigen = patheigen.cur, pathirm = pathirm.cur)
                    
                } 
                
                pop_prob.new <- pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
                
                path_prob.new <- path_prob(path = path.new, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
                path_prob.cur <- path_prob(path = path.cur, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
                
                a.prob <- accept_prob(pop_prob.new = pop_prob.new, pop_prob.cur = pop_prob.cur, path_prob.cur = path_prob.cur, path_prob.new = path_prob.new)
                
                if(min(a.prob, 0) > log(runif(1))) {
                    X.cur <- X.new
                    Xcount.cur <- Xcount.new
                    W.cur <- W.new
                    pop_prob.cur <- pop_prob.new
                    
                }
                
                acceptances[j] <- a.prob
                
            }
        }
        
    }
    #    
#     c(ifelse(obstimes < getpath(SIRres$trajectory, 1)[2], 1, 0), ifelse(obstimes < path.new[2], 1, 0)) 
    
        return(list(Truth = SIRres$trajectory,  Augmented = X.cur, acceptances = acceptances, true_durs = infec_durs(SIRres$trajectory, 1), aug_durs = infec_durs(X.cur, 1)))  
}
stopCluster(cl)

print(Sys.time() - strt)

save(sim_one_results, file = "augSIR_jointcomp_sim_one.Rdata")


# Plots -------------------------------------------------------------------


statuses <- matrix(0, nrow = (2*length(seq(0,tmax, by=0.05))), ncol = length(sim_one_results))

for(k in 1:ncol(statuses)){
    
    statuses[,k] <- sim_one_results[[k]]
    
}

means <- rowMeans(statuses)


path_comp <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))

path_comp2 <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))

path_comp2[,2] <- path_comp2[,2] - path_comp2[1:(nrow(path_comp[2])/2), 2]


print(ggplot(path_comp, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual."))
print(ggplot(path_comp2, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual. Subtracting Gillespie means."))

# more plotting code
# 
# 
# censusInterval <- 0.05; times <- seq(0, tmax, by = censusInterval)
# 
# resultscomp <- list()
# 
# for(k in 1:(length(sim_one_results))){
#     if(k %% 1000 == 0) print(k)
#     
#     traj <- build_countmat(sim_one_results[[k]][[1]], popsize)
#     augSIR_count <- build_countmat(sim_one_results[[k]][[2]], popsize)
#     
#     traj_full <- rep(0, length(times))
#     augSIR_full <- rep(0, length(times))
#     
#     for(j in 1:length(times)){
#         
#         traj_full[j] <- traj[sum(traj[,1] <= times[j]), 2]
#         augSIR_full[j] <- augSIR_count[sum(augSIR_count[,1] <= times[j]), 2]
#         
#     }
#     
#     resultscomp[[k]] <- data.frame(iter = k, cbind(times, traj_full, augSIR_full))
#     
# }
# 
# gillespie_durs <- do.call(rbind,lapply(sim_one_results, "[[", 4))
# augSIR_durs <- do.call(rbind,lapply(sim_one_results, "[[", 5))
# 
# allpaths <- matrix(0, nrow = length(times)*2, ncol = samp_size)
# 
# for(k in 1:samp_size){
#     
#     if(k %% 5000 == 0) print(k)
#     
#     allpaths[,k] <- c(ifelse(times < gillespie_durs[k], 1, 0), ifelse(times < augSIR_durs[k], 1, 0))
#     
# }
# 
# allpaths_comp <- data.frame(time = rep(times, 2), method = rep(c("Gillespie", "augSIR"), each = length(times)), infected = rowMeans(allpaths))
# 
# ggplot(allpaths_comp, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status - resampling one individual")
# 
# 
# # plot average infection status at each time point
# 
# path_comp <- data.frame(time = seq(0,tmax,by=0.05), infected = 0, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))
# 
# 
# cl <- makeCluster(4)
# registerDoParallel(cl)
# 
# means <- foreach(k = 1:length(times)) %dopar% {
#     
#     gillespie_mean <- mean(do.call(rbind,lapply(lapply(resultscomp,"[[",3), "[", k)), na.rm = T)
#     augSIR_mean <- mean(do.call(rbind,lapply(lapply(resultscomp,"[[",4), "[", k)), na.rm = T)
#     
#     c(gillespie_mean, augSIR_mean)
# }
# 
# stopCluster(cl)
# 
# meansbind <- do.call(rbind, means)
# 
# path_comp$infected[1:length(times)] <- meansbind[,1]
# path_comp$infected[(length(times)+1):nrow(path_comp)] <- meansbind[,2]
# 
# ggplot(path_comp, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average number of infecteds - resampling one individual")
# 










# Simulation yielding "incorrect" results ----------------------------------------------------


library(doParallel)

niter <- 10; popsize = 4; tmax = 10
b <- 0.5 + runif(1, -0.0001, 0.0001)
m <- 1 + runif(1, -0.0001, 0.0001)
samp_prob <- 0.5 
initdist <- c(0, 1, 0)
samp_size = 10000


cl <- makeCluster(3)
registerDoParallel(cl)

strt <- Sys.time()
writeLines(c(""), "log2.txt")

sim_one_results <- foreach(iter = 1:samp_size, .packages = "augSIR") %dopar% {
    
    if(iter %% 1000 == 0){
        sink("log2.txt", append=TRUE)  
        cat(paste("Starting iteration",iter,"\n"))  
        sink()
    }
    
    SIRres<-SIRsim(popsize = 4, initdist = initdist, b = b, mu=m, a=0, tmax = tmax, censusInterval=0.05, sampprob = 0.5, returnX = TRUE, trim = FALSE)
    
    # make sure the epidemic is interesting (i.e. at least one infection)
    if(max(SIRres$results[,3]) < 2 | SIRres$results[1,2] == 0){
        
        while(max(SIRres$results[,3]) < 2 | SIRres$results[1,2] == 0){
            
            SIRres<-SIRsim(popsize = 4, initdist = initdist, b = b, mu=m, a=0, tmax = tmax, censusInterval=0.05, sampprob = 0.5, returnX = TRUE, trim = FALSE)
            
        }
    }
    
    # initialize objects
    #     accepts <- rep(0, length(niter)); accepts[1] <- 1
    #     accept_probs <- rep(0, length(niter)); accept_probs[1] <- 0
    #     log_likelihood <- rep(0, length(niter))
    #     trajectories <- vector("list", length(niter))
    #     paths <- vector("list", length(niter))
    
    # observation matrix
    W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Observed, augmented = 0))
    
    # individual trajectories
    X.cur <- SIRres$trajectory
    
    # count matrix
    Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
    #     trajectories[[1]] <- Xcount.cur
    
    # update observation matrix
    W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
    
    pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
    
    # build irm matrices
    pathirm.cur <- build_irm(Xcount = Xcount.cur, b = b, m = m, a = 0, popsize = popsize, pop = FALSE)
    patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
    
    subjects <- 1
    
    acceptances <- rep(0, length(subjects))
    
    for(j in 1:length(subjects)){
        
        Xother <- X.cur[X.cur[,2]!=subjects[j],]
        
        path.cur <- getpath(X.cur, subjects[j])
        
        if(!all(path.cur == Xcount.cur[,1])){
            
            Xcount.other <- build_countmat(Xother, popsize - 1)
            
            if(nrow(Xcount.other)!=1){
                
                W.other <-get_W_other(W.cur = W.cur, path = path.cur)
                
                path.new<- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = samp_prob, initdist = initdist, tmax = tmax)
                
                X.new <- updateX(X = X.cur, path = path.new, j = subjects[j])
                Xcount.new <- update_Xcount(Xcount.other = Xcount.other, path = path.new)
                
                W.new <- updateW(W = W.other, path = path.new)
                
                if(max(Xcount.new[,2]) == pathirm.cur[4,4,dim(pathirm.cur)[3]]){
                    
                    new.numinf <- pathirm.cur[4,4,dim(pathirm.cur)[3]]+1
                    
                    pathirm.cur <- update_irm(irm = pathirm.cur, new.numinf = new.numinf, b = b, m = m, a = 0, popsize = popsize)
                    patheigen.cur <- update_eigen(patheigen = patheigen.cur, pathirm = pathirm.cur)
                    
                } 
                
                
                pop_prob.new <- pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
                
                path_prob.new <- path_prob(path = path.new, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
                path_prob.cur <- path_prob(path = path.cur, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
                
                a.prob <- accept_prob(pop_prob.new = pop_prob.new, pop_prob.cur = pop_prob.cur, path_prob.cur = path_prob.cur, path_prob.new = path_prob.new)
                
                if(min(a.prob, 0) > log(runif(1))) {
                    X.cur <- X.new
                    Xcount.cur <- Xcount.new
                    W.cur <- W.new
                    pop_prob.cur <- pop_prob.new
                    
                }
                
                acceptances[j] <- a.prob
            }
        }
        
    }
    
    return(list(Truth = SIRres$trajectory,  Augmented = X.cur, acceptances = acceptances, true_durs = infec_durs(SIRres$trajectory, 1), aug_durs = infec_durs(X.cur, 1)))  
}
stopCluster(cl)

print(Sys.time() - strt)

save(sim_one_results, file = "augSIR_jointcomp_sim_one_noMH.Rdata")


# Plots -------------------------------------------------------------------


censusInterval <- 0.05; times <- seq(0, tmax, by = censusInterval)

resultscomp <- list()

for(k in 1:(length(sim_one_results))){
    if(k %% 1000 == 0) print(k)
        
    traj <- build_countmat(sim_one_results[[k]][[1]], popsize)
    augSIR_count <- build_countmat(sim_one_results[[k]][[2]], popsize)
    
    traj_full <- rep(0, length(times))
    augSIR_full <- rep(0, length(times))
    
    for(j in 1:length(times)){
        
        traj_full[j] <- traj[sum(traj[,1] <= times[j]), 2]
        augSIR_full[j] <- augSIR_count[sum(augSIR_count[,1] <= times[j]), 2]
        
    }
        
    resultscomp[[k]] <- data.frame(iter = k, cbind(times, traj_full, augSIR_full))

}

gillespie_durs <- do.call(rbind,lapply(sim_one_results, "[[", 4))
augSIR_durs <- do.call(rbind,lapply(sim_one_results, "[[", 5))

allpaths <- matrix(0, nrow = length(times)*2, ncol = samp_size)

for(k in 1:samp_size){
    
    if(k %% 5000 == 0) print(k)
    
    allpaths[,k] <- c(ifelse(times < gillespie_durs[k], 1, 0), ifelse(times < augSIR_durs[k], 1, 0))
    
}
    
allpaths_comp <- data.frame(time = rep(times, 2), method = rep(c("Gillespie", "augSIR"), each = length(times)), infected = rowMeans(allpaths))

ggplot(allpaths_comp, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status - resampling one individual")


# plot average infection status at each time point

path_comp <- data.frame(time = seq(0,tmax,by=0.05), infected = 0, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))


cl <- makeCluster(4)
registerDoParallel(cl)

means <- foreach(k = 1:length(times)) %dopar% {
    
    gillespie_mean <- mean(do.call(rbind,lapply(lapply(resultscomp,"[[",3), "[", k)), na.rm = T)
    augSIR_mean <- mean(do.call(rbind,lapply(lapply(resultscomp,"[[",4), "[", k)), na.rm = T)
    
    c(gillespie_mean, augSIR_mean)
}

stopCluster(cl)

meansbind <- do.call(rbind, means)

path_comp$infected[1:length(times)] <- meansbind[,1]
path_comp$infected[(length(times)+1):nrow(path_comp)] <- meansbind[,2]

ggplot(path_comp, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average number of infecteds - resampling one individual")


