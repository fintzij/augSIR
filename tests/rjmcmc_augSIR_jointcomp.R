library(doParallel)

niter <- 10; popsize = 4; tmax = 10
b <- 0.5 + runif(1, -0.0001, 0.0001)
m <- 1 + runif(1, -0.0001, 0.0001)
samp_prob <- 0.5 
initdist <- c(0, 1, 0)
samp_size = 25000


cl <- makeCluster(4)
registerDoParallel(cl)

strt <- Sys.time()
writeLines(c(""), "log2.txt")

augSIR_jointcomp_results <- foreach(iter = 1:samp_size, .packages = "augSIR") %dopar% {
    
    if(iter %% 10 == 0){
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
    W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Truth, augmented = 0))
    
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
    
    subjects <- rep(1:popsize, each = niter)
    
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
    
    return(list(Truth = SIRres$trajectory, Observed = SIRres$results, Augmented = Xcount.cur, acceptances = acceptances, true_durs = infec_durs(SIRres$trajectory, 1:popsize), aug_durs = infec_durs(X.cur, 1:popsize)))  
}
stopCluster(cl)

print(Sys.time() - strt)

save(augSIR_jointcomp_results, file = "augSIR_jointcomp_results2.Rdata")





# Now repeat for rjmcmc
# 
# # rjmcmc settings
# popsize = 4; tmax = 4
# b <- 0.5 + runif(1,-0.0001, 0.0001)
# m <- 1 + runif(1, -0.0001, 0.0001)
# samp_prob <- 0.5 
# initdist <- c(0.7, 0.3, 0)
# insert.prob = 1/3; remove.prob = 1/3; shift.prob = 1/3
# shift.int <- 0.1
# niter <- 10000
# samp_size = 25000
# 
# cl <- makeCluster(3)
# registerDoParallel(cl)
# 
# strt <- Sys.time()
# writeLines(c(""), "log.txt")
# 
# augSIR_jointcomp_results <- foreach(iter = 1:samp_size, .packages = "augSIR") %dopar% {
#     
#     if(iter %% 5 == 0){
#         sink("log.txt", append=TRUE)  
#         cat(paste("Starting iteration",iter,"\n"))  
#         sink()
#     }
#     
#     SIRres<-SIRsim(popsize = 4, initdist = c(0.7, 0.3, 0), b = 0.5, mu=1, a=0, tmax = 4, censusInterval=0.05, sampprob = samp_prob, returnX = TRUE, trim = FALSE)
#     
#     # make sure the epidemic is interesting (i.e. at least one infection)
#     if(max(SIRres$results[,3]) < 2){
#         
#         while(max(SIRres$results[,3]) < 2){
#             
#             SIRres<-SIRsim(popsize = 4, initdist = c(0.7, 0.3, 0), b = 0.5, mu=1, a=0, tmax = 4, censusInterval=0.05, sampprob = samp_prob, returnX = TRUE, trim = FALSE)
#             
#         }
#     }
#     
#     # initialize objects
#     #     accepts <- rep(0, length(niter)); accepts[1] <- 1
#     #     accept_probs <- rep(0, length(niter)); accept_probs[1] <- 0
#     #     log_likelihood <- rep(0, length(niter))
#     #     trajectories <- vector("list", length(niter))
#     #     paths <- vector("list", length(niter))
#     
#     # observation matrix
#     W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Truth, augmented = 0))
#     
#     # individual trajectories
#     X.cur <- SIRres$trajectory
#     
#     # count matrix
#     Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
#     #     trajectories[[1]] <- Xcount.cur
#     
#     # update observation matrix
#     W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
#     
#     subjects <- rep(1:popsize, each = niter)
#     
#     
#     for(j in 1:length(subjects)){
#         
#         path.cur <- X.cur[X.cur[,2] == subjects[j], ]
#         path.new <- rjmcmc_draw(path.cur = path.cur, Xcount.cur, j = subjects[j], initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, tmax = max(W.cur[,1]), b = b, m = m, p = samp_prob)
#         
#         # update bookkeeping objects
#         X.new <- X.cur; X.new[X.new[,2] == subjects[j], ] <- path.new
#         Xcount.new <- build_countmat(X = X.new, popsize = popsize)
#         W.new <- updateW(W = W.cur, Xcount = Xcount.new)
#         
#         # calculate acceptance ratio
#         rjmcmc.ratio <- rjmcmc_ratio(W.cur = W.cur, W.new = W.new, X.cur = X.cur, X.new = X.new, Xcount.cur = Xcount.cur, Xcount.new = Xcount.new, path.cur = path.cur, path.new = path.new, initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, b = b, m = m, samp_prob = samp_prob, tmax = tmax, popsize = popsize)
#         
#         # decide whether to accept. if accept, update current trajectories
#         if(rjmcmc.ratio > log(runif(1))){
#             
#             X.cur <- X.new
#             Xcount.cur <- Xcount.new
#             W.cur <- W.new
#             
#         }
#         
#     }
#     
#     list(Truth = SIRres$trajectory, Observed = SIRres$results, Augmented = Xcount.cur)  
# }
# stopCluster(cl)
# 
# print(Sys.time() - strt)
# 
# 
# save(augSIR_jointcomp_results, file = "rjmcmc_jointcomp_results.Rdata")



# Plots -------------------------------------------------------------------
censusInterval <- 0.05

resultscomp <- list()

for(k in 1:(length(augSIR_jointcomp_results))){
    if(k %% 1000 == 0) print(k)
    
    traj <- augSIR_jointcomp_results[[k]][[2]]
    augSIR_count <- augSIR_jointcomp_results[[k]][[3]]    
    
    augSIR_full <- rep(0, nrow(traj))
    
    for(j in 1:81){
        
        augSIR_full[j] <- augSIR_count[sum(augSIR_count[,1] <= traj[j,1]), 2]
        
    }
    
    augSIR_samp <- rbinom(n = length(augSIR_full), size = augSIR_full, prob = samp_prob)
    
    resultscomp[[k]] <- data.frame(iter = k, cbind(traj, augSIR_full, augSIR_samp))
}
    



# plot average infection status at each time point

path_comp <- data.frame(time = seq(0,4,by=0.05), infected = 0, method = rep(c("Gillespie Sample","Gillespie Full", "augSIR Full", "augSIR Sample"), each = length(seq(0,4,by=0.05))))


for(k in 1:81){
    print(k)
    
    path_comp[k,2] <- mean(do.call(rbind,lapply(lapply(resultscomp,"[[",3), "[", k)), na.rm = T)
    path_comp[k + 81,2] <- mean(do.call(rbind,lapply(lapply(resultscomp,"[[",4), "[", k)), na.rm = T)
    path_comp[k + 2*81,2] <- mean(do.call(rbind,lapply(lapply(resultscomp,"[[",5), "[", k)),  na.rm = T)
    path_comp[k + 3*81,2] <- mean(do.call(rbind,lapply(lapply(resultscomp,"[[",6), "[", k)), na.rm = T)
        
}

ggplot(path_comp, aes(x=time, y = infected, colour = method)) + geom_line()



# joint likelihood comparison ---------------------------------------------

