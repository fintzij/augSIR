niter <- 2; popsize = 4; tmax = 10
b <- 0.5 + runif(1, -0.0001, 0.0001)
m <- 1 + runif(1, -0.0001, 0.0001)
samp_prob <- 0.5 
initdist <- c(0, 1, 0)
samp_size = 10000
obstimes <- seq(0, tmax, by = 0.05)

results <- vector(mode = "list", length = samp_size)
paths <- matrix(0, nrow = 2*length(obstimes), ncol = samp_size)

for(k in 1:samp_size){
    
    if(k%%1000 == 0) print(k)
    
    # individual trajectories
    X.cur <- SIRres$trajectory
    
    # count matrix
    Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
        
    # observation matrix
    W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Truth, augmented = 0))
    W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
    
    # calculate population trajectory likelihood
    pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
    
    # build irm matrices
    pathirm.cur <- build_irm(Xcount = Xcount.cur, b = b, m = m, a = 0, popsize = popsize, pop = FALSE)
    patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
    
    # vector of subjects to resample (just resample subject 1)
    subjects <- rep(1, niter)
    
    # vector to store objects from kth sample
    results_j <- vector(mode = "list", length = niter)

    for(j in 1:length(subjects)){
        
        # get current path
        path.cur <- getpath(X.cur, subjects[j])
        
        # get .other objects
        Xother <- X.cur[X.cur[,2]!=subjects[j],]
        Xcount.other <- build_countmat(Xother, popsize - 1)
        
        W.other <- get_W_other(W.cur, path.cur)    
        
        
        # draw new path
        path.new <- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = samp_prob, initdist = initdist, tmax = tmax)
        
        
        # get .new objects
        X.new <- updateX(X = X.cur, path = path.new, j = subjects[j])
        Xcount.new <- update_Xcount(Xcount.other = Xcount.other, path = path.new)
        
        W.new <- updateW(W = W.other, path = path.new)
        
        # metropolis hastings
#         pop_prob.new <- pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
#         
#         path_prob.new <- path_prob(path = path.new, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
#         path_prob.cur <- path_prob(path = path.cur, Xcount = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
#         
#         a.prob <- accept_prob(pop_prob.new = pop_prob.new, pop_prob.cur = pop_prob.cur, path_prob.cur = path_prob.cur, path_prob.new = path_prob.new)
#         
#         if(min(a.prob, 0) > log(runif(1))) {
#             X.cur <- X.new
#             Xcount.cur <- Xcount.new
#             W.cur <- W.new
#             pop_prob.cur <- pop_prob.new
#             
#         }
        
        results_j[[j]] <- list(path.cur = path.cur,
                            path.new = path.new,
                            X.cur = X.cur,
                            X.new = X.new,
                            Xcount.cur = Xcount.cur,
                            Xcount.new = Xcount.new,
                            W.cur = W.cur,
                            W.new = W.new,
                            Xother = Xother,
                            Xcount.other = Xcount.other,
                            W.other = W.other)
        
    }

    gillespie.path <- getpath(SIRres$trajectory, 1)
    augSIR.path <- path.new

    paths[,k] <- c(ifelse(obstimes < gillespie.path[2], 1, 0), ifelse(obstimes < augSIR.path[2], 1, 0))
    
    results[[k]] <- results_j
}


# get object lists
X.other.list <- lapply(results[[1]], "[[", 9)
X.other.list <- lapply(results[[1]], "[[", 10)
X.other.list <- lapply(results[[1]], "[[", 11)




# plots -------------------------------------------------------------------

statuses <- matrix(0, nrow = (2*length(seq(0,tmax, by=0.05))), ncol = samp_size)

for(k in 1:ncol(statuses)){
    
    statuses[,k] <- paths[[k]]
    
}

means <- rowMeans(statuses)


path_comp <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))

path_comp2 <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))

path_comp2[,2] <- path_comp2[,2] - path_comp2[1:(nrow(path_comp[2])/2), 2]


print(ggplot(path_comp, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual"))
print(ggplot(path_comp2, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual. No HMM. Subtracting Gillespie means."))

