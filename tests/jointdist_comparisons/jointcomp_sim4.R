
library(ggplot2, lib.loc = "/home/students/fintzij/R_Packages")

# set working directory and load augSIR files
source("augSIR.R")
source("auxilliary_functions.R")
source("SIRsimulation.R")
source("rjmcmc_functions.R")
source("matrix_build_update.R")
source("metropolis_hastings_functions.R")
source("path_sampling_functions.R")


# set simulation parameters
niter <- 1; samp_size <- 100000
popsize = 4; tmax = 10
b <- 0.5 + runif(1, -0.0001, 0.0001)
m <- 1 + runif(1, -0.0001, 0.0001)
samp_prob <- 0.5 
initdist <- c(0, 1, 0)
obstimes <- seq(0, tmax, by=0.05)

writeLines(c(""), "log4.txt")
statuses <- matrix(0, nrow = (2*length(seq(0,tmax, by=0.05))), ncol = samp_size)

for(k in 1:samp_size){
    
    if(k %% 1000 == 0){
        sink("log4.txt", append=TRUE)  
        cat(paste("Starting iteration",k,"\n"))  
        sink()
    }
    
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
    #     W.cur[,2] <- rbinom(nrow(W.cur), W.cur[,2], samp_prob)
    
    
    # build irm matrices
    pathirm.cur <- build_irm(Xcount = Xcount.cur, b = b, m = m, a = 0, popsize = popsize, pop = FALSE)
    patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
    
    subjects <- rep(1:4, niter)
    
    acceptances <- rep(0, length(subjects))
    
    for(j in 1:length(subjects)){
        
        # get current path
        path.cur <- getpath(X.cur, subjects[j])
        
        # get .other matrices
        Xother <- X.cur[X.cur[,2]!=subjects[j],]
        Xcount.other <- build_countmat(Xother, popsize - 1)
        W.other <-get_W_other(W.cur = W.cur, path = path.cur) 
        
        # draw new path
        path.new<- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = samp_prob, initdist = initdist, tmax = tmax)
        
        # generate .new objects
        X.new <- updateX(X = X.cur, path = path.new, j = subjects[j])
        Xcount.new <- update_Xcount(Xcount.other = Xcount.other, path = path.new)
        W.new <- updateW(W = W.other, path = path.new)
        
        # Metropolis-Hastings (should always accept, but include check to make sure)
        # population trajectory likelihoods 
        pop_prob.cur <- pop_prob(Xcount = Xcount.cur, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
        pop_prob.new <- pop_prob(Xcount = Xcount.new, tmax = tmax, b = b, m = m, a = 0, initdist = initdist, popsize = popsize)
        
        # path likelihoods
        path_prob.new <- path_prob(path = path.new, Xcount.other = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
        path_prob.cur <- path_prob(path = path.cur, Xcount.other = Xcount.other, pathirm = pathirm.cur, initdist = initdist, tmax = tmax)
        
        # compute log acceptance ratio
        a.prob <- accept_prob(pop_prob.new = pop_prob.new, pop_prob.cur = pop_prob.cur, path_prob.cur = path_prob.cur, path_prob.new = path_prob.new)
        
        
        if(min(a.prob, 0) > log(runif(1))) {
            X.cur <- X.new
            Xcount.cur <- Xcount.new
            W.cur <- W.new
            pop_prob.cur <- pop_prob.new
            
        }
    }

    statuses[,k] <- c(SIRres$results$Truth, W.cur[,3])
}

# plot results
means <- rowMeans(statuses)

path_comp <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))

path_comp2 <- data.frame(time = rep(seq(0,tmax,by=0.05),2), infected = means, method = rep(c("Gillespie", "augSIR"), each = length(seq(0,tmax,by=0.05))))

path_comp2[,2] <- path_comp2[,2] - path_comp2[1:(nrow(path_comp[2])/2), 2]


png("infecprob_sim4.png", width = 750, height = 550)
print(ggplot(path_comp, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual."))
dev.off()

png("infecprob_demeaned_sim4.png", width = 750, height = 550)
print(ggplot(path_comp2, aes(x=time, y = infected, colour = method)) + geom_line() + labs(title = "Average infection status for one individual. Subtracting Gillespie means."))
dev.off()
