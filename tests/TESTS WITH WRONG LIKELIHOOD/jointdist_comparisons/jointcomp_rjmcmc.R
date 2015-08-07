# 
# library(ggplot2, lib.loc = "/home/students/fintzij/R_Packages")
# 
# # set working directory and load augSIR files
# source("augSIR.R")
# source("auxilliary_functions.R")
# source("SIRsimulation.R")
# source("rjmcmc_functions.R")
# source("matrix_build_update.R")
# source("metropolis_hastings_functions.R")
# source("path_sampling_functions.R")


# set simulation parameters
niter <- 1; samp_size <- 250000
popsize = 4; tmax = 10
b <- 0.5 + runif(1, -0.0001, 0.0001)
m <- 1 + runif(1, -0.0001, 0.0001)
accepts <- 0
samp_prob <- 0.5 
initdist <- c(0, 1, 0)
obstimes <- seq(0, tmax, by=0.05)
insert.prob = 1/3; remove.prob = 1/3; shift.prob = 1/3
shift.int <- 0.5

# writeLines(c(""), "log1.txt")
statuses <- matrix(0, nrow = (2*length(seq(0,tmax, by=0.05))), ncol = samp_size)

for(k in 1:samp_size){
    if(k%%1000 == 0) print(k)
#     if(k %% 1000 == 0){
#         sink("log1.txt", append=TRUE)  
#         cat(paste("Starting iteration",k,"\n"))  
#         sink()
#     }
    
    # simlate dataset
    SIRres<-SIRsim(popsize = 4, initdist = initdist, b = b, mu=m, a=0, tmax = tmax, censusInterval=0.05, sampprob = 0.5, returnX = TRUE, trim = FALSE)
    
    # observation matrix
    W.cur <- as.matrix(data.frame(time = SIRres$results$time, sampled = SIRres$results$Observed, augmented = 0))
    
    # individual trajectories
    X.cur <- SIRres$trajectory
    
    # count matrix
    Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
    
    # update observation matrix
    W.cur <- updateW(W = W.cur, Xcount = Xcount.cur)
    
    # build irm matrices
    pathirm.cur <- build_irm(Xcount = Xcount.cur, b = b, m = m, a = 0, popsize = popsize, pop = FALSE)
    patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
    
    subjects <- rep(1, niter)
    
    acceptances <- rep(0, length(subjects))
    
    for(j in 1:length(subjects)){
        
        # get current path
        path.cur <- X.cur[X.cur[,2] == subjects[j], ]
        
        # get .other matrices
        Xother <- X.cur[X.cur[,2]!=subjects[j],]
        Xcount.other <- build_countmat(Xother, popsize - 1)
        W.other <-get_W_other(W.cur = W.cur, path = path.cur) 
        
        # draw new path
        path.new <- rjmcmc_draw(path.cur = path.cur, Xcount.cur, j = subjects[j], initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, tmax = max(W.cur[,1]), b = b, m = m, p = samp_prob)
        
        # generate .new objects
        X.new <- X.cur; X.new[X.new[,2] == subjects[j], ] <- path.new
        Xcount.new <- build_countmat(X = X.new, popsize = popsize)
        W.new <- updateW(W = W.cur, Xcount = Xcount.new)
        
        # Metropolis-Hastings (should always accept, but include check to make sure)
        
        # compute log acceptance ratio
        rjmcmc.ratio <- rjmcmc_ratio(W.cur = W.cur, W.new = W.new, X.cur = X.cur, X.new = X.new, Xcount.cur = Xcount.cur, Xcount.new = Xcount.new, path.cur = path.cur, path.new = path.new, initdist = initdist, shift.int = shift.int, insert.prob = insert.prob, remove.prob = remove.prob, shift.prob = shift.prob, b = b, m = m, samp_prob = samp_prob, tmax = tmax, popsize = popsize)
        
        
        if(min(rjmcmc.ratio, 0) > log(runif(1))) {
            X.cur <- X.new
            Xcount.cur <- Xcount.new
            W.cur <- W.new
            # pop_prob.cur <- pop_prob.new
            accepts <- accepts+1
        } else path.new <- path.cur
    }
    
    statuses[,k] <- c(ifelse(obstimes < path.cur[2,1], 1, 0), ifelse(obstimes < path.new[2,1], 1, 0))
}

# plot results
means <- rowMeans(statuses)

path_comp <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "RJMCMC"), each = length(seq(0,tmax,by=0.05))))

path_comp2 <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "RJMCMC"), each = length(seq(0,tmax,by=0.05))))

path_comp2[,2] <- path_comp2[,2] - path_comp2[1:(nrow(path_comp[2])/2), 2]


# png("infecprob_sim1.png", width = 750, height = 550)
print(ggplot(path_comp, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual."))
# dev.off()

# png("infecprob_demeaned_sim1.png", width = 750, height = 550)
print(ggplot(path_comp2, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual. Subtracting Gillespie means."))
# dev.off()
