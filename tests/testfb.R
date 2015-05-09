
# Checking emission probabilities -----------------------------------------
# 
# # set objects
# popsize = 200; b = 0.01 + runif(1,-0.005,0.005); m = 0.5 + runif(1,-0.005,0.005); a = 0; tmax = 25; p = 0.25; initdist <- c(0.995, 0.005, 0)
# obstimes <- seq(0,25 ,by=0.1); numsick = c(rep(c(20,75),each=100),rep(20,51))
# eventtimes <- c(0,10,20)
# 
# # set observation matrix
# W.cur <- as.matrix(data.frame(time = seq(0,25,by=0.1), sampled = 0, augmented = c(rep(c(20,75),each=100),rep(20,51))))
# W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
# 
# # set trajectory matrix
# X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
# X.cur[2*(1:75)+1,3] <- 1; X.cur[2*(1:55)+2,3] <- -1
# X.cur[2*(1:55) + 1, 1] <- 10; X.cur[2*(1:55) + 2, 1] <- 20
# 
# 
# X.cur <- X.cur[order(X.cur[,1]),]
# Xcount <- build_countmat(X.cur,popsize); Xcount <- Xcount[c(1,56,111),]
# pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
# patheigen.cur <- irm_decomp(pathirm.cur)
# 
# tpms <- tpm_seq(Xcount, obstimes, patheigen.cur)
# 
# Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- Xcount
# W.other <- updateW(W = W.cur, Xcount = Xcount.other)
# 
# cl <- makeCluster(3)
# registerDoParallel(cl)
# 
# emit.sim<- foreach(t = 1:100000, .packages = 'augSIR') %dopar% {
#     
#     # draw new binomial samples and update the W matrix
#     W.cur[,2] <- rbinom(n = dim(W.cur)[1], size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
#     
#     emit_seq(W.cur,p)    
# 
# }
# stopCluster(cl)
# 
# emit.sim <- simplify2array(emit.sim)
# infec.emits <- rep(0,251)
# for(k in 1:251){
#     infec.emits[k] <- mean(emit.sim[2,k,])
# }
# 


# 2/18/2015 - more diagnostics --------------------------------------------

# set parameters and objects
popsize = 200; b = 0.01 + runif(1,-0.005,0.005); m = 0.5 + runif(1,-0.005,0.005); a = 0; tmax = 25; p = 0.25; initdist <- c(0.995, 0.005, 0)
numinf <- rep(40, length(eventtimes))

# again, with more variation in the rates
obstimes <- seq(0,25 ,by=0.1); numsick = c(rep(c(5,10,20,40,75,40,20,10,5),each=22),rep(5,3)); eventtimes = c(2.7*(0:8))

# set observation matrix
W.cur <- as.matrix(data.frame(time = seq(0,25,by=0.1), sampled = 0, augmented = c(rep(c(5,10,20,40,75,40,20,10,5),each=27),rep(5,8))))
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix

# set trajectory matrix and count matrix
X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:75)+1,3] <- 1; X.cur[2*(6:75)+2,3] <- -1
X.cur[2*(6:10) + 1, 1] <- eventtimes[2]; X.cur[2*(6:10) + 2, 1] <- eventtimes[9]
X.cur[2*(11:20) + 1, 1] <- eventtimes[3]; X.cur[2*(11:20) + 2, 1] <- eventtimes[8]
X.cur[2*(21:40) + 1, 1] <- eventtimes[4]; X.cur[2*(21:40) + 2, 1] <- eventtimes[7]
X.cur[2*(41:75) + 1, 1] <- eventtimes[5]; X.cur[2*(41:75) + 2, 1] <- eventtimes[6]

X.cur <- X.cur[order(X.cur[,1]),]
Xcount <- build_countmat(X.cur,popsize); Xcount <- Xcount[c(1,6,16,36,71,106,126,136,141),]

# initialize irms and tpms
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)

tpms <- tpm_seq(Xcount, obstimes, patheigen.cur)[[1]]

# initialize .other matrices
Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- Xcount
W.other <- updateW(W = W.cur, Xcount = Xcount.other)


cl <- makeCluster(4)
registerDoParallel(cl)
draws.varinf.sim1 <- foreach(t = 1:250000, .packages='augSIR') %dopar% {
        
    pathvec <- sim_one_SIR(Xcount = Xcount, obstimes = obstimes, b=b, m=m, initdist = initdist, tmax = tmax)
    #         pathvec <- sim_one_tpms(tpms = tpms, initdist = initdist)
    
    # update X and W matrices with the new marginal sample path
    W.cur[,3] <- W.cur[,3] + (pathvec == 2) # update the observation matrix with the new trajectory
    
    # draw new binomial samples and update the W matrix
    W.cur[,2] <- rbinom(n = dim(W.cur)[1], size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
    
    W.cur[,3] <- W.cur[,3] - (pathvec == 2)
    
    # draw the new trajectory
    list(pathvec,drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = p, initdist = initdist)[,2])
    # draw new trajectory using new binomial samples
    }
stopCluster(cl)


method1.draws <- lapply(draws.varinf.sim1, '[[', 1)
method1.draws <- lapply(method1.draws, function(x) matrix(unlist(x),ncol=251,byrow=TRUE))
method1.draws <- do.call(rbind,method1.draws)
method1.draws <- data.frame(cbind(1:dim(method1.draws)[1],method1.draws)); names(method1.draws) <- c("simnum", obstimes)
method1.melt <- melt(method1.draws, id = "simnum"); names(method1.melt)[2] <- "time"
method1.melt$time <- as.numeric(method1.melt$time)/10

method2.draws <- lapply(draws.varinf.sim1, '[[', 2)
method2.draws <- lapply(method2.draws, function(x) matrix(unlist(x),ncol=251,byrow=TRUE))
method2.draws <- do.call(rbind,method2.draws)
method2.draws <- data.frame(cbind(1:dim(method2.draws)[1],method2.draws)); names(method2.draws) <- c("simnum", obstimes)
method2.melt <- melt(method2.draws, id = "simnum"); names(method2.melt)[2] <- "time"
method2.melt$time <- as.numeric(method2.melt$time)/10

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "5, 10, 20, 40, 75, 40, 20, 10, 5 Infecteds")

save(draws.varinf.sim1,file="draws.varinf.sim1.Rdata")
rm(method1.draws, method2.draws, draws.varinf.sim1)

######## Now with fewer time varying rates
# 
# # set objects
# popsize = 200; b = 0.01 + runif(1,-0.005,0.005); m = 0.5 + runif(1,-0.005,0.005); a = 0; tmax = 25; p = 0.25; initdist <- c(0.995, 0.005, 0)
# obstimes <- seq(0,25 ,by=0.1); numsick = c(rep(c(20,75),each=100),rep(20,51))
# eventtimes <- c(0,10,20)
# 
# # set observation matrix
# W.cur <- as.matrix(data.frame(time = seq(0,25,by=0.1), sampled = 0, augmented = c(rep(c(20,75),each=100),rep(20,51))))
# W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
# 
# # set trajectory matrix
# X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
# X.cur[2*(1:75)+1,3] <- 1; X.cur[2*(1:55)+2,3] <- -1
# X.cur[2*(1:55) + 1, 1] <- 10; X.cur[2*(1:55) + 2, 1] <- 20
# 
# 
# X.cur <- X.cur[order(X.cur[,1]),]
# Xcount <- build_countmat(X.cur,popsize); Xcount <- Xcount[c(1,56,111),]
# pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
# patheigen.cur <- irm_decomp(pathirm.cur)
# 
# tpms <- tpm_seq(Xcount, obstimes, patheigen.cur)[[1]]
# 
# Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- Xcount
# W.other <- updateW(W = W.cur, Xcount = Xcount.other)
# 
# cl <- makeCluster(3)
# registerDoParallel(cl)
# draws.varinf.sim2 <- foreach(t = 1:250000, .packages = 'augSIR') %dopar% {
#         
#     pathvec.sir <- sim_one_SIR(Xcount = Xcount, obstimes = obstimes, b=b, m=m, initdist = initdist, tmax = tmax)
#     pathvec.tpm <- sim_one_tpms(tpms = tpms, initdist = initdist)
#     
#     # update X and W matrices with the new marginal sample path
#     W.cur[,3] <- W.cur[,3] + (pathvec.tpm == 2) # update the observation matrix with the new trajectory
#     
#     # draw new binomial samples and update the W matrix
#     W.cur[,2] <- rbinom(n = dim(W.cur)[1], size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
#     
#     W.cur[,3] <- W.cur[,3] - (pathvec.tpm == 2)
#     
#     # draw the new trajectory
#     list(pathvec.sir, pathvec.tpm, drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.cur, p = p, initdist = initdist)[,2])
#     # draw new trajectory using new binomial samples
#     }
# 
# stopCluster(cl)
# 
# save(draws.varinf.sim2,file="draws.varinf.sim2.Rdata")
# 
# load("draws.varinf.sim2.Rdata")
# 
# method1.draws <- lapply(draws.varinf.sim2, '[[', 1)
# save(method1.draws,file="method1.sim2.Rdata"); rm(method1.draws)
# 
# method2.draws <- lapply(draws.varinf.sim2, '[[', 2)
# save(method2.draws,file="method2.sim2.Rdata"); rm(method2.draws)
# 
# method3.draws <- lapply(draws.varinf.sim2, '[[', 3)
# save(method3.draws,file="method3.sim2.Rdata"); rm(draws.varinf.sim2)
# 
# 
# # convert model draws to matrices
# load("method1.sim2.Rdata")
# 
# method1.draws.mat <- lapply(method1.draws, function(x) matrix(unlist(x),ncol=251,byrow=TRUE)); rm(method1.draws)
# method1.draws <- do.call(rbind,method1.draws.mat); rm(method1.draws.mat)
# method1.draws.mat <- cbind(1:dim(method1.draws)[1],method1.draws)
# rm(method1.draws); colnames(method1.draws.mat) <- c("simnum", obstimes)
# method1.draws.mat <- as.data.frame(method1.draws.mat)
# save(method1.draws.mat,file="method1.sim2.Rdata"); rm(method1.draws.mat)
# 
# 
# load("method2.sim2.Rdata")
# 
# method2.draws.mat <- lapply(method2.draws, function(x) matrix(unlist(x),ncol=251,byrow=TRUE)); rm(method2.draws)
# method2.draws <- do.call(rbind,method2.draws.mat); rm(method2.draws.mat)
# method2.draws.mat <- cbind(1:dim(method2.draws)[1],method2.draws)
# rm(method2.draws); colnames(method2.draws.mat) <- c("simnum", obstimes)
# method2.draws.mat <- as.data.frame(method2.draws.mat)
# save(method2.draws.mat,file="method2.sim2.Rdata"); rm(method2.draws.mat)
# 
# 
# load("method3.sim2.Rdata")
# 
# method3.draws.mat <- lapply(method3.draws, function(x) matrix(unlist(x),ncol=251,byrow=TRUE)); rm(method3.draws)
# method3.draws <- do.call(rbind,method3.draws.mat); rm(method3.draws.mat)
# method3.draws.mat <- cbind(1:dim(method3.draws)[1],method3.draws)
# rm(method3.draws); colnames(method3.draws.mat) <- c("simnum", obstimes)
# method3.draws.mat <- as.data.frame(method3.draws.mat)
# save(method3.draws.mat,file="method3.sim2.Rdata"); rm(method3.draws.mat)
# 
# 
# load("method1.sim2.Rdata"); load("method2.sim2.Rdata"); load("method3.sim2.Rdata")
# 
# comp.summary <- data.frame(Method = c(rep("Method 1 - SIR", length(obstimes)), rep("Method 1 - TPM", length(obstimes)), rep("Method 2 - FB", length(obstimes))),
#                            Time = rep(obstimes,3), 
#                            Susceptible = 0,
#                            Infected = 0, 
#                            Recovered = 0)
# 
# for(r in 1:length(obstimes)){
#     comp.summary[r,3] <- mean(method1.draws.mat[,r+1]==1)
#     comp.summary[r,4] <- mean(method1.draws.mat[,r+1]==2)
#     comp.summary[r,5] <- mean(method1.draws.mat[,r+1]==3)
#     
#     comp.summary[r + length(obstimes),3] <- mean(method2.draws.mat[,r+1]==1)
#     comp.summary[r + length(obstimes),4] <- mean(method2.draws.mat[,r+1]==2)
#     comp.summary[r + length(obstimes),5] <- mean(method2.draws.mat[,r+1]==3)
#     
#     comp.summary[r + 2*length(obstimes),3] <- mean(method3.draws.mat[,r+1]==1)
#     comp.summary[r + 2*length(obstimes),4] <- mean(method3.draws.mat[,r+1]==2)
#     comp.summary[r + 2*length(obstimes),5] <- mean(method3.draws.mat[,r+1]==3)
#     
# }
# 
# summary.melt <- melt(comp.summary, id = c("Method", "Time"))
# ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_line() + labs(y="Probability", title = "20, 75, 20 Infecteds")
# 
# # rm(method1.draws, method2.draws, draws.varinf.sim2)
# 
# 
# ###### Comparing to theoretical quantities
# # escape probabiliy out of segment 1
# exp(-b*20*10)
# mean(method2.draws[method2.draws[,1]==1,101]==1)
# 
# # escape probability out of segment 2, given escape from segment 1
# mean(method2.draws[method2.draws[,28]==1,55]==1)
# 
# # escape probability out of segment 3, given escape from segment 2
# mean(method2.draws[])
# 



### run simulation - constant rate of 40 infecteds - simulating from tpms

# set objects
popsize = 200; b = 0.01 + runif(1,-0.005,0.005); m = 0.5 + runif(1,-0.005,0.005); a = 0; tmax = 25; p = 0.25; initdist <- c(0.995, 0.005, 0)
obstimes <- seq(0,25 ,by=0.1); numsick <- rep(40, 251)

# set observation matrix
W.cur <- as.matrix(data.frame(time = seq(0,25,by=0.1), sampled = 0, augmented = 40))
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix

# set trajectory matrix
X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:40)+1,3] <- 1
X.cur <- X.cur[order(X.cur[,1]),]

Xcount <- build_countmat(X.cur,popsize)
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)

tpms <- tpm_seq(Xcount, obstimes, patheigen.cur)[[1]]

Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- Xcount
W.other <- updateW(W = W.cur, Xcount = Xcount.other)


cl <- makeCluster(3)
registerDoParallel(cl)
draws.40inf <- foreach(t = 1:250000, .packages='augSIR') %dopar% {
        
        pathvec <- sim_one_tpms(tpms = tpms, initdist = initdist)
        
        # update X and W matrices with the new marginal sample path
        W.cur[,3] <- W.cur[,3] + (pathvec == 2) # update the observation matrix with the new trajectory
        
        # draw new binomial samples and update the W matrix
        W.cur[,2] <- rbinom(n = dim(W.cur)[1], size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
        
        W.cur[,3] <- W.cur[,3] - (pathvec == 2)
        
        # draw fb path
        fbpath <- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.cur, p = p, initdist = initdist, tmax = tmax)
        
        fbpathvec <- ifelse(obstimes <= fbpath[1], 1, ifelse(obstimes > fbpath[2], 3, 2))
        # draw the new trajectory
        list(pathvec, fbpathvec) # draw new trajectory using new binomial samples
    }
stopCluster(cl)

method1.draws <- lapply(draws.40inf, '[[', 1)
method1.draws <- lapply(method1.draws, function(x) matrix(unlist(x),ncol=251,byrow=TRUE))
method1.draws <- do.call(rbind,method1.draws)
method1.draws <- data.frame(cbind(1:dim(method1.draws)[1],method1.draws)); names(method1.draws) <- c("simnum", obstimes)

method2.draws <- lapply(draws.40inf, '[[', 2)
method2.draws <- lapply(method2.draws, function(x) matrix(unlist(x),ncol=251,byrow=TRUE))
method2.draws <- do.call(rbind,method2.draws)
method2.draws <- data.frame(cbind(1:dim(method2.draws)[1],method2.draws)); names(method2.draws) <- c("simnum", obstimes)

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "40 Infecteds")


rm(draws.40inf)


######## Now with time varying rates

# set objects
popsize = 200; b = 0.01 + runif(1,-0.005,0.005); m = 0.5 + runif(1,-0.005,0.005); a = 0; tmax = 25; p = 0.25; initdist <- c(0.995, 0.005, 0)
obstimes <- seq(0,25 ,by=0.1); numsick = c(rep(c(20,75),each=100),rep(20,51))

# set observation matrix
W.cur <- as.matrix(data.frame(time = seq(0,25,by=0.1), sampled = 0, augmented = c(rep(c(20,75),each=100),rep(20,51))))
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix

# set trajectory matrix
X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:75)+1,3] <- 1; X.cur[2*(1:55)+2,3] <- -1
X.cur[2*(1:55) + 1, 1] <- 10; X.cur[2*(1:55) + 2, 1] <- 20


X.cur <- X.cur[order(X.cur[,1]),]
Xcount <- build_countmat(X.cur,popsize); Xcount <- Xcount[c(1,56,111),]
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)

tpms <- tpm_seq(Xcount, obstimes, patheigen.cur)[[1]]

Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- Xcount
W.other <- updateW(W = W.cur, Xcount = Xcount.other)

cl <- makeCluster(3)
registerDoParallel(cl)
draws.varinf <- foreach(t = 1:50000, .packages='augSIR') %dopar% {
        
        pathvec <- sim_one_tpms(tpms = tpms, initdist = initdist)
        
        # update X and W matrices with the new marginal sample path
        W.cur[,3] <- W.cur[,3] + (pathvec == 2) # update the observation matrix with the new trajectory
        
        # draw new binomial samples and update the W matrix
        W.cur[,2] <- rbinom(n = dim(W.cur)[1], size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
        
        W.cur[,3] <- W.cur[,3] - (pathvec == 2)
        
        fbpath <- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.cur, p = p, initdist = initdist, tmax = tmax)
        
        fbpathvec <- ifelse(obstimes <= fbpath[1], 1, ifelse(obstimes > fbpath[2], 3, 2))
        # draw the new trajectory
        list(pathvec, fbpathvec) # draw new trajectory using new binomial samples
    
}
stopCluster(cl)


method1.draws <- lapply(draws.varinf, '[[', 1)
method1.draws <- lapply(method1.draws, function(x) matrix(unlist(x),ncol=251,byrow=TRUE))
method1.draws <- do.call(rbind,method1.draws)
method1.draws <- data.frame(cbind(1:dim(method1.draws)[1],method1.draws)); names(method1.draws) <- c("simnum", obstimes)

method2.draws <- lapply(draws.varinf, '[[', 2)
method2.draws <- lapply(method2.draws, function(x) matrix(unlist(x),ncol=251,byrow=TRUE))
method2.draws <- do.call(rbind,method2.draws)
method2.draws <- data.frame(cbind(1:dim(method2.draws)[1],method2.draws)); names(method2.draws) <- c("simnum", obstimes)

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "20, 75, 20 Infecteds")


rm(draws.varinf)


# again, with more variation in the rates
popsize = 200; b = 0.01 + runif(1,-0.005,0.005); m = 0.5 + runif(1,-0.005,0.005); a = 0; tmax = 25; p = 0.25; initdist <- c(0.995, 0.005, 0)
obstimes <- seq(0,25 ,by=0.1); numsick = c(rep(c(5,10,20,40,75,40,20,10,5),each=22),rep(5,3)); eventtimes = 2.7*(1:8)

# set observation matrix
W.cur <- as.matrix(data.frame(time = seq(0,25,by=0.1), sampled = 0, augmented = c(rep(c(5,10,20,40,75,40,20,10,5),each=27),rep(5,8))))
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix

# set trajectory matrix and count matrix
X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:75)+1,3] <- 1; X.cur[2*(6:75)+2,3] <- -1
X.cur[2*(6:10) + 1, 1] <- eventtimes[1]; X.cur[2*(6:10) + 2, 1] <- eventtimes[8]
X.cur[2*(11:20) + 1, 1] <- eventtimes[2]; X.cur[2*(11:20) + 2, 1] <- eventtimes[7]
X.cur[2*(21:40) + 1, 1] <- eventtimes[3]; X.cur[2*(21:40) + 2, 1] <- eventtimes[6]
X.cur[2*(41:75) + 1, 1] <- eventtimes[4]; X.cur[2*(41:75) + 2, 1] <- eventtimes[5]

X.cur <- X.cur[order(X.cur[,1]),]
Xcount <- build_countmat(X.cur,popsize); Xcount <- Xcount[c(1,6,16,36,71,106,126,136,141),]

# initialize irms and tpms
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)

tpms <- tpm_seq(Xcount, obstimes, patheigen.cur)[[1]]

# initialize .other matrices
Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- Xcount
W.other <- updateW(W = W.cur, Xcount = Xcount.other)


cl <- makeCluster(3)
registerDoParallel(cl)
draws.varinf <- foreach(t = 1:100000, .packages = 'augSIR') %dopar% {
    
    pathvec <- sim_one_tpms(tpms = tpms, initdist = initdist)
    
    # update X and W matrices with the new marginal sample path
    W.cur[,3] <- W.cur[,3] + (pathvec == 2) # update the observation matrix with the new trajectory
    
    # draw new binomial samples and update the W matrix
    W.cur[,2] <- rbinom(n = dim(W.cur)[1], size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
    
    W.cur[,3] <- W.cur[,3] - (pathvec == 2)
    
    fbpath <- draw_path(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.cur, p = p, initdist = initdist, tmax = tmax)
    
    fbpathvec <- ifelse(obstimes <= fbpath[1], 1, ifelse(obstimes > fbpath[2], 3, 2))
    # draw the new trajectory
    list(pathvec, fbpathvec) # draw new trajectory using new binomial samples
}
stopCluster(cl)


method1.draws <- lapply(draws.varinf, '[[', 1)
method1.draws <- lapply(method1.draws, function(x) matrix(unlist(x),ncol=251,byrow=TRUE))
method1.draws <- do.call(rbind,method1.draws)
method1.draws <- data.frame(cbind(1:dim(method1.draws)[1],method1.draws)); names(method1.draws) <- c("simnum", obstimes)

method2.draws <- lapply(draws.varinf, '[[', 2)
method2.draws <- lapply(method2.draws, function(x) matrix(unlist(x),ncol=251,byrow=TRUE))
method2.draws <- do.call(rbind,method2.draws)
method2.draws <- data.frame(cbind(1:dim(method2.draws)[1],method2.draws)); names(method2.draws) <- c("simnum", obstimes)

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "5, 10, 20, 40, 75, 40, 20, 10, 5 Infecteds")


rm(draws.varinf)



# 2/11/2015 - method 1 vs method 2 - one subject  --------

#### First with the number of infecteds constant at 5
obstimes <- seq(0,50,by=0.1)
W.cur <- as.matrix(data.frame(time = seq(0,50,by=0.1), sampled = 0, augmented = 5))
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
numsick = rep(5,400); eventtimes = seq(0,50,length=400)


# cl <- makeCluster(3)
# registerDoParallel(cl)
# method1.5inf <- foreach(t = 1:200000) %dopar% sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = FALSE)
# stopCluster(cl)

X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:5)+1,3] <- 1; X.cur <- X.cur[order(X.cur[,3]),]
Xcount <- build_countmat(X.cur,popsize)
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)

Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- build_countmat(Xother, popsize)
W.other <- updateW(W = W.cur, Xcount = Xcount.other)

cl <- makeCluster(3)
registerDoParallel(cl)
draws.5inf <- foreach(t = 1:10) %do% {
    print(t)
    assign(paste("draws.5inf",t,sep=""), foreach(t = 1:50000, .combine = 'comb', .multicombine=TRUE, .init=list(list(),list())) %dopar% {
    
    path <- sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = TRUE)
    pathvec <- ifelse(obstimes <= path[1], 1, ifelse(obstimes <= path[2] & obstimes > path[1], 2, 3))
    
    # update X and W matrices with the new marginal sample path
    X.cur <- updateX(X = X.cur, Xt.path = path, j=1) # upate the complete data matrix with the new trajectory
    Xcount <- build_countmat(X.cur,popsize)
    W.cur <- updateW(W = W.cur, Xcount = Xcount) # update the observation matrix with the new trajectory
    
    # draw new binomial samples and update the W matrix
    W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
    
    # draw the new trajectory
    list(pathvec,drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = p, initdist = initdist)[,2]) # draw new trajectory using new binomial samples
    })
}
stopCluster(cl)


method1.draws <- do.call(rbind, draws.5inf[[1]])[1,]; method1.draws <- do.call(rbind, method1.draws)
method1.draws <- data.frame(cbind(1:500000,method1.draws)); names(method1.draws) <- c("simnum", obstimes)
method1.melt <- melt(method1.draws, id = "simnum"); names(method1.melt)[2] <- "time"
method1.melt$time <- as.numeric(method1.melt$time)/10

method2.draws <- do.call(rbind,draws.5inf[[2]])[2,]; method2.draws <- do.call(rbind, method2.draws)
method2.draws <- data.frame(cbind(1:500000,method2.draws)); names(method2.draws) <- c("simnum", obstimes)
method2.melt <- melt(method2.draws, id = "simnum"); names(method2.melt)[2] <- "time"
method2.melt$time <- as.numeric(method2.melt$time)/10

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "5 Infecteds")

rm(draws.5inf)

###### Number of infecteds constant at 10

obstimes <- seq(0,30,by=0.1)
W.cur <- as.matrix(data.frame(time = seq(0,30,by=0.1), sampled = 0, augmented = 10))
numsick = rep(10,400); eventtimes = seq(0,30,length=400)
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix

# cl <- makeCluster(3)
# registerDoParallel(cl)
# method1.10inf <- foreach(t = 1:200000) %dopar% sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = FALSE)
# stopCluster(cl)

X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:10)+1,3] <- 1; X.cur <- X.cur[order(X.cur[,3]),]
Xcount <- build_countmat(X.cur,popsize)
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)

Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- build_countmat(Xother, popsize)
W.other <- updateW(W = W.cur, Xcount = Xcount.other)

cl <- makeCluster(3)
registerDoParallel(cl)
draws.10inf <- foreach(t = 1:200000, .combine='comb', .multicombine=TRUE, .init=list(list(), list())) %dopar% {
    
    path <- sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = TRUE)
    pathvec <- ifelse(obstimes <= path[1], 1, ifelse(obstimes <= path[2] & obstimes > path[1], 2, 3))
    
    # update X and W matrices with the new marginal sample path
    X.cur <- updateX(X = X.cur, Xt.path = path, j=1) # upate the complete data matrix with the new trajectory
    Xcount <- build_countmat(X.cur,popsize)
    W.cur <- updateW(W = W.cur, Xcount = Xcount) # update the observation matrix with the new trajectory
    
    # draw new binomial samples and update the W matrix
    W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
    
    # draw the new trajectory
    list(pathvec, drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = p, initdist = initdist)[,2]) # draw new trajectory using new binomial samples
    
}
stopCluster(cl)


method1.draws <- do.call(rbind,draws.10inf[[1]])
method1.draws <- data.frame(cbind(1:200000,method1.draws)); names(method1.draws) <- c("simnum", obstimes)
method1.melt <- melt(method1.draws, id = "simnum"); names(method1.melt)[2] <- "time"
method1.melt$time <- as.numeric(method1.melt$time)/10

method2.draws <- do.call(rbind,draws.10inf[[2]])
method2.draws <- data.frame(cbind(1:200000,method2.draws)); names(method2.draws) <- c("simnum", obstimes)
method2.melt <- melt(method2.draws, id = "simnum"); names(method2.melt)[2] <- "time"
method2.melt$time <- as.numeric(method2.melt$time)/10

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "10 Infecteds")

rm(draws.10inf)
##### Now with 20 infecteds
obstimes <- seq(0,25,by=0.1)
W.cur <- as.matrix(data.frame(time = seq(0,25,by=0.1), sampled = 0, augmented = 20))
numsick = rep(20,400); eventtimes = seq(0,25,length=400)
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix

# cl <- makeCluster(3)
# registerDoParallel(cl)
# method1.20inf <- foreach(t = 1:200000) %dopar% sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = FALSE)
# stopCluster(cl)


X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:20)+1,3] <- 1; X.cur <- X.cur[order(X.cur[,3]),]
Xcount <- build_countmat(X.cur,popsize)
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)

Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- build_countmat(Xother, popsize)
W.other <- updateW(W = W.cur, Xcount = Xcount.other)

cl <- makeCluster(3)
registerDoParallel(cl)
draws.20inf <- foreach(t = 1:200000, .combine='comb', .multicombine=TRUE, .init=list(list(), list())) %dopar% {
    
    path <- sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = TRUE)
    pathvec <- ifelse(obstimes <= path[1], 1, ifelse(obstimes <= path[2] & obstimes > path[1], 2, 3))
    
    # update X and W matrices with the new marginal sample path
    X.cur <- updateX(X = X.cur, Xt.path = path, j=1) # upate the complete data matrix with the new trajectory
    Xcount <- build_countmat(X.cur,popsize)
    W.cur <- updateW(W = W.cur, Xcount = Xcount) # update the observation matrix with the new trajectory
    
    # draw new binomial samples and update the W matrix
    W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
    
    # draw the new trajectory
    list(pathvec, drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = p, initdist = initdist)[,2]) # draw new trajectory using new binomial samples
    
}
stopCluster(cl)


method1.draws <- do.call(rbind,draws.20inf[[1]])
method1.draws <- data.frame(cbind(1:200000,method1.draws)); names(method1.draws) <- c("simnum", obstimes)
method1.melt <- melt(method1.draws, id = "simnum"); names(method1.melt)[2] <- "time"
method1.melt$time <- as.numeric(method1.melt$time)/10

method2.draws <- do.call(rbind,draws.20inf[[2]])
method2.draws <- data.frame(cbind(1:200000,method2.draws)); names(method2.draws) <- c("simnum", obstimes)
method2.melt <- melt(method2.draws, id = "simnum"); names(method2.melt)[2] <- "time"
method2.melt$time <- as.numeric(method2.melt$time)/10

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "20 Infecteds")

rm(draws.20inf)
##### Now with 40 infecteds
obstimes <- seq(0,15,by=0.1)
W.cur <- as.matrix(data.frame(time = seq(0,15,by=0.1), sampled = 0, augmented = 40))
numsick = rep(40,400); eventtimes = seq(0,15,length=400)
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix

# cl <- makeCluster(3)
# registerDoParallel(cl)
# method1.40inf <- foreach(t = 1:200000) %dopar% sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = FALSE)
# stopCluster(cl)


X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:40)+1,3] <- 1; X.cur <- X.cur[order(X.cur[,3]),]
Xcount <- build_countmat(X.cur,popsize)
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)

Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- build_countmat(Xother, popsize)
W.other <- updateW(W = W.cur, Xcount = Xcount.other)


cl <- makeCluster(3)
registerDoParallel(cl)
draws.40inf <- foreach(t = 1:200000, .combine='comb', .multicombine=TRUE, .init=list(list(), list())) %dopar% {
    
    path <- sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = TRUE)
    pathvec <- ifelse(obstimes <= path[1], 1, ifelse(obstimes <= path[2] & obstimes > path[1], 2, 3))
    
    # update X and W matrices with the new marginal sample path
    X.cur <- updateX(X = X.cur, Xt.path = path, j=1) # upate the complete data matrix with the new trajectory
    Xcount <- build_countmat(X.cur,popsize)
    W.cur <- updateW(W = W.cur, Xcount = Xcount) # update the observation matrix with the new trajectory
    
    # draw new binomial samples and update the W matrix
    W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
    
    # draw the new trajectory
    list(pathvec, drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = p, initdist = initdist)[,2]) # draw new trajectory using new binomial samples
    
}
stopCluster(cl)


method1.draws <- do.call(rbind,draws.40inf[[1]])
method1.draws <- data.frame(cbind(1:200000,method1.draws)); names(method1.draws) <- c("simnum", obstimes)
method1.melt <- melt(method1.draws, id = "simnum"); names(method1.melt)[2] <- "time"
method1.melt$time <- as.numeric(method1.melt$time)/10

method2.draws <- do.call(rbind,draws.40inf[[2]])
method2.draws <- data.frame(cbind(1:200000,method2.draws)); names(method2.draws) <- c("simnum", obstimes)
method2.melt <- melt(method2.draws, id = "simnum"); names(method2.melt)[2] <- "time"
method2.melt$time <- as.numeric(method2.melt$time)/10

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "40 Infecteds")


rm(draws.40inf)

##### Now with 75 infecteds
obstimes <- seq(0,10,by=0.1)
W.cur <- as.matrix(data.frame(time = seq(0,10,by=0.1), sampled = 0, augmented = 75))
numsick = rep(75,400); eventtimes = seq(0,10,length=400)
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix

# cl <- makeCluster(3)
# registerDoParallel(cl)
# method1.75inf <- foreach(t = 1:200000) %dopar% sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = FALSE)
# stopCluster(cl)


X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:75)+1,3] <- 1; X.cur <- X.cur[order(X.cur[,3]),]
Xcount <- build_countmat(X.cur,popsize)
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)

Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- build_countmat(Xother, popsize)
W.other <- updateW(W = W.cur, Xcount = Xcount.other)

cl <- makeCluster(4)
registerDoParallel(cl)
draws.75inf <- foreach(t = 1:200000, .combine='comb', .multicombine=TRUE, .init=list(list(), list())) %dopar% {
    
    path <- sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = TRUE)
    pathvec <- ifelse(obstimes <= path[1], 1, ifelse(obstimes <= path[2] & obstimes > path[1], 2, 3))
    
    # update X and W matrices with the new marginal sample path
    X.cur <- updateX(X = X.cur, Xt.path = path, j=1) # upate the complete data matrix with the new trajectory
    Xcount <- build_countmat(X.cur,popsize)
    W.cur <- updateW(W = W.cur, Xcount = Xcount) # update the observation matrix with the new trajectory
    
    # draw new binomial samples and update the W matrix
    W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
    
    # draw the new trajectory
    list(pathvec, drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = p, initdist = initdist)[,2]) # draw new trajectory using new binomial samples
    
}
stopCluster(cl)


method1.draws <- do.call(rbind,draws.75inf[[1]])
method1.draws <- data.frame(cbind(1:200000,method1.draws)); names(method1.draws) <- c("simnum", obstimes)
method1.melt <- melt(method1.draws, id = "simnum"); names(method1.melt)[2] <- "time"
method1.melt$time <- as.numeric(method1.melt$time)/10

method2.draws <- do.call(rbind,draws.75inf[[2]])
method2.draws <- data.frame(cbind(1:200000,method2.draws)); names(method2.draws) <- c("simnum", obstimes)
method2.melt <- melt(method2.draws, id = "simnum"); names(method2.melt)[2] <- "time"
method2.melt$time <- as.numeric(method2.melt$time)/10

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "75 Infecteds")


rm(draws.75inf)

##### Now with 120 infecteds
obstimes <- seq(0,8,by=0.1)
W.cur <- as.matrix(data.frame(time = seq(0,8,by=0.1), sampled = 0, augmented = 120))
numsick = rep(120,400); eventtimes = seq(0,8,length=400)
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix

# cl <- makeCluster(3)
# registerDoParallel(cl)
# method1.120inf <- foreach(t = 1:200000) %dopar% sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = FALSE)
# stopCluster(cl)


X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:120)+1,3] <- 1; X.cur <- X.cur[order(X.cur[,3]),]
Xcount <- build_countmat(X.cur,popsize)
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)

Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- build_countmat(Xother, popsize)
W.other <- updateW(W = W.cur, Xcount = Xcount.other)

cl <- makeCluster(3)
registerDoParallel(cl)
draws.120inf <- foreach(t = 1:200000, .combine='comb', .multicombine=TRUE, .init=list(list(), list())) %dopar% {
    
    path <- sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = TRUE)
    pathvec <- ifelse(obstimes <= path[1], 1, ifelse(obstimes <= path[2] & obstimes > path[1], 2, 3))
    
    # update X and W matrices with the new marginal sample path
    X.cur <- updateX(X = X.cur, Xt.path = path, j=1) # upate the complete data matrix with the new trajectory
    Xcount <- build_countmat(X.cur,popsize)
    W.cur <- updateW(W = W.cur, Xcount = Xcount) # update the observation matrix with the new trajectory
    
    # draw new binomial samples and update the W matrix
    W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
    
    # draw the new trajectory
    list(pathvec, drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = p, initdist = initdist)[,2]) # draw new trajectory using new binomial samples
    
}
stopCluster(cl)


method1.draws <- do.call(rbind,draws.120inf[[1]])
method1.draws <- data.frame(cbind(1:200000,method1.draws)); names(method1.draws) <- c("simnum", obstimes)
method1.melt <- melt(method1.draws, id = "simnum"); names(method1.melt)[2] <- "time"
method1.melt$time <- as.numeric(method1.melt$time)/10

method2.draws <- do.call(rbind,draws.120inf[[2]])
method2.draws <- data.frame(cbind(1:200000,method2.draws)); names(method2.draws) <- c("simnum", obstimes)
method2.melt <- melt(method2.draws, id = "simnum"); names(method2.melt)[2] <- "time"
method2.melt$time <- as.numeric(method2.melt$time)/10

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "120 Infecteds")



rm(draws.120inf)


######## Now with time varying rates

obstimes <- seq(0,25 ,by=0.1)
W.cur <- as.matrix(data.frame(time = seq(0,25,by=0.1), sampled = 0, augmented = c(rep(c(20,75),each=100),rep(20,51))))
numsick = c(rep(c(20,75),each=100),rep(20,51)); 
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix

# cl <- makeCluster(3)
# registerDoParallel(cl)
# method1.120inf <- foreach(t = 1:200000) %dopar% sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = FALSE)
# stopCluster(cl)


X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:75)+1,3] <- 1; X.cur[2*(1:55)+2,3] <- -1
X.cur[2*(1:55) + 1, 1] <- 10; X.cur[2*(1:55) + 2, 1] <- 20


X.cur <- X.cur[order(X.cur[,1]),]
Xcount <- build_countmat(X.cur,popsize); Xcount <- Xcount[c(1,56,111),]
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)
patheigen.cur[,,1,76] <- patheigen.cur[c(2,1,3),c(2,1,3),1,76]
patheigen.cur[,,2,76] <- patheigen.cur[,c(2,1,3),2,76]
patheigen.cur[,,3,76] <- solve(patheigen.cur[,,2,76])

Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- Xcount
W.other <- updateW(W = W.cur, Xcount = Xcount.other)

cl <- makeCluster(3)
registerDoParallel(cl)
draws.varinf <- foreach(t = 1:10) %do% {
    print(t)
    assign(paste("draws.varinf",t,sep=""), foreach(t = 1:50000, .combine = 'comb', .multicombine=TRUE, .init=list(list(),list())) %dopar% {
        
        path <- sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = TRUE)
        pathvec <- ifelse(obstimes <= path[1], 1, ifelse(obstimes <= path[2] & obstimes > path[1], 2, 3))
        
        # update X and W matrices with the new marginal sample path
        X.cur <- updateX(X = X.cur, Xt.path = path, j=1) # upate the complete data matrix with the new trajectory
        Xcount <- build_countmat(X.cur,popsize)
        W.cur <- updateW(W = W.cur, Xcount = Xcount) # update the observation matrix with the new trajectory
        
        # draw new binomial samples and update the W matrix
        W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
        
        # draw the new trajectory
        list(pathvec,drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = p, initdist = initdist)[,2]) # draw new trajectory using new binomial samples
    })
}
stopCluster(cl)


method1.draws <- do.call(rbind,draws.varinf[[1]])[1,]; method1.draws <- do.call(rbind, method1.draws)
method1.draws <- data.frame(cbind(1:500000,method1.draws)); names(method1.draws) <- c("simnum", obstimes)
method1.melt <- melt(method1.draws, id = "simnum"); names(method1.melt)[2] <- "time"
method1.melt$time <- as.numeric(method1.melt$time)/10

method2.draws <- do.call(rbind,draws.varinf[[2]])[1,]; method2.draws <- do.call(rbind, method2.draws)
method2.draws <- data.frame(cbind(1:500000,method2.draws)); names(method2.draws) <- c("simnum", obstimes)
method2.melt <- melt(method2.draws, id = "simnum"); names(method2.melt)[2] <- "time"
method2.melt$time <- as.numeric(method2.melt$time)/10

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "20, 75, 20 Infecteds")


rm(draws.varinf)


# again, with more variation in the rates
obstimes <- seq(0,25 ,by=0.1)
W.cur <- as.matrix(data.frame(time = seq(0,25,by=0.1), sampled = 0, augmented = c(rep(c(5,10,20,40,75,40,20,10,5),each=27),rep(5,8))))
numsick = c(rep(c(5,10,20,40,75,40,20,10,5),each=22),rep(5,3)); eventtimes = 2.7*(1:8)
W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix

# cl <- makeCluster(3)
# registerDoParallel(cl)
# method1.120inf <- foreach(t = 1:200000) %dopar% sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = FALSE)
# stopCluster(cl)


X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))
X.cur[2*(1:75)+1,3] <- 1; X.cur[2*(6:75)+2,3] <- -1
X.cur[2*(6:10) + 1, 1] <- eventtimes[1]; X.cur[2*(6:10) + 2, 1] <- eventtimes[8]
X.cur[2*(11:20) + 1, 1] <- eventtimes[2]; X.cur[2*(11:20) + 2, 1] <- eventtimes[7]
X.cur[2*(21:40) + 1, 1] <- eventtimes[3]; X.cur[2*(21:40) + 2, 1] <- eventtimes[6]
X.cur[2*(41:75) + 1, 1] <- eventtimes[4]; X.cur[2*(41:75) + 2, 1] <- eventtimes[5]

X.cur <- X.cur[order(X.cur[,1]),]
Xcount <- build_countmat(X.cur,popsize); Xcount <- Xcount[c(1,6,16,36,71,106,126,136,141),]
pathirm.cur <- build_irm(Xcount = Xcount, b = b, m = m, a = 0, popsize = popsize, pop=FALSE)
patheigen.cur <- irm_decomp(pathirm.cur)

Xother <- X.cur[X.cur[,2] != 1, ]; Xcount.other <- Xcount
W.other <- updateW(W = W.cur, Xcount = Xcount.other)


cl <- makeCluster(3)
registerDoParallel(cl)
draws.varinf <- foreach(t = 1:10) %do% {
    print(t)
    assign(paste("draws.varinf",t,sep=""), foreach(t = 1:50000, .combine = 'comb', .multicombine=TRUE, .init=list(list(),list())) %dopar% {
        
        path <- sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = TRUE)
        pathvec <- ifelse(obstimes <= path[1], 1, ifelse(obstimes <= path[2] & obstimes > path[1], 2, 3))
        
        # update X and W matrices with the new marginal sample path
        X.cur <- updateX(X = X.cur, Xt.path = path, j=1) # upate the complete data matrix with the new trajectory
        Xcount <- build_countmat(X.cur,popsize)
        W.cur <- updateW(W = W.cur, Xcount = Xcount) # update the observation matrix with the new trajectory
        
        # draw new binomial samples and update the W matrix
        W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
        
        # draw the new trajectory
        list(pathvec,drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = p, initdist = initdist)[,2]) # draw new trajectory using new binomial samples
    })
}
stopCluster(cl)


method1.draws <- do.call(rbind,draws.varinf[[1]])[1,]; method1.draws <- do.call(rbind, method1.draws)
method1.draws <- data.frame(cbind(1:500000,method1.draws)); names(method1.draws) <- c("simnum", obstimes)
method1.melt <- melt(method1.draws, id = "simnum"); names(method1.melt)[2] <- "time"
method1.melt$time <- as.numeric(method1.melt$time)/10

method2.draws <- do.call(rbind,draws.varinf[[2]])[1,]; method2.draws <- do.call(rbind, method2.draws)
method2.draws <- data.frame(cbind(1:500000,method2.draws)); names(method2.draws) <- c("simnum", obstimes)
method2.melt <- melt(method2.draws, id = "simnum"); names(method2.melt)[2] <- "time"
method2.melt$time <- as.numeric(method2.melt$time)/10

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability", title = "5, 10, 20, 40, 75, 40, 20, 10, 5 Infecteds")


rm(draws.varinf)

# 2/10/2015 - Checking against HMM marginal probabilities ------------------
# 
# tMat <- matrix(c(.95, 0, 0, .05, .9, 0, 0, .1, 1), nrow=3)
# 
# eMat <- cbind( rep(1/6,6)            # fair
#                ,c(.1,.1,.3,.3,.1,.1) # biased toward 3 and 4
#                ,c(.3,.1,.1,.1,.1,.3)) # biased toward 1 and 6
# initDist <- c(1/3,1/3,1/3)
# 
# forwd <- function(distr,obs) normalize(( t( tMat ) %*% distr ) * eMat[obs,] )
# backwd <- function(obs,distr) normalize( tMat %*% ( distr * eMat[obs,] ))
# normalize <- function(x) x/sum(x)
# 
# dieRolls <- c(sample(6,10,replace=TRUE), sample(6, 10, replace = TRUE, prob = c(.1,.1,.25,.25,.1,.1)), sample(6, 10, replace = TRUE, prob = c(.3,.1,.1,.1,.1,.3)))
# forwardProbs <- Reduce(forwd,dieRolls,init=initDist,acc=T)
# backwardsProbs <- Reduce(backwd,dieRolls,init=normalize(c(1,1,1)),acc=T,right=T)
# marginalProbs <- apply(mapply(`*`,forwardProbs,backwardsProbs), 2, normalize)
# 
# # Slightly modifying the forward function from my algorithm
# fbmats <- array(0, dim = c(3,3,length(dieRolls) - 1))
#     
# pi0 <- normalize(eMat[dieRolls[1],] * initDist)
# fbmats[,,1] <- normalize(outer(pi0,eMat[dieRolls[2],], FUN="*") * tMat)
#     
# for(s in 2:dim(fbmats)[3]){
#     fbmats[,,s] <- normalize(outer(colSums(fbmats[,,s-1]), eMat[dieRolls[s+1],], FUN="*") * tMat)
#     
# }
#     
# myMarginals <- cbind(normalize(eMat[dieRolls[1],] * initDist), apply(fbmats, 3, colSums))
# 

# # Method 1 vs Method 2 comparison -----------------------------------------
# 
# SIRres<-SIRsim(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=0.5, a=0, tmax = 15, censusInterval=0.25, sampprob = 0.25, binomsamp = TRUE, returnX = TRUE)
# 
# if(dim(SIRres$results)[1] < 10){
#     while(dim(SIRres$results)[1] < 10){
#         SIRres<-SIRsim(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=0.5, a=0, tmax = 15, censusInterval=0.25, sampprob = 0.25, binomsamp = TRUE, returnX = TRUE)
#         
#     }
# }
# 
# # get data 
# dat <- SIRres$results
# dat.m <- melt(dat,id.vars="time")
# 
# ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + theme_bw()
# 
# sim.settings <- list(popsize = 200,
#                      tmax = max(dat[,1] + 1),
#                      niter = 50,
#                      amplify = 10,
#                      initdist = c(0.995, 0.005, 0))
# 
# inits <- list(beta.init = 0.01 + runif(1,-0.0005, 0.0005),
#               mu.init = 0.5 + runif(1, -0.005, 0.005),
#               alpha.init = 0, 
#               probs.init = 0.25 + runif(1,-0.005, 0.005))
# 
# priors <- list(beta.prior = c(.012, 1.1),
#                mu.prior = c(0.96, 1.96),
#                alpha.prior = NULL,
#                p.prior = c(4.5, 13.35))
# 
# popsize <- sim.settings$popsize # size of the population; 
# tmax <- sim.settings$tmax # maximum time of observation
# amplify <- sim.settings$amplify # amplification parameter for initialization
# niter <- sim.settings$niter # number of iterations in the sampler
# initdist <- sim.settings$initdist # initial distribution for individual infection status
# 
# # vectors for parameters
# Beta <- vector(length=niter); Beta[1] <- inits$beta.init
# Mu <- vector(length = niter); Mu[1] <- inits$mu.init
# Alpha <- vector(length = niter); Alpha[1] <- inits$alpha.init 
# probs <- vector(length = niter); probs[1] <- inits$probs.init
# accepts <- vector(length = niter)
# 
# # vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
# beta.prior <- priors$beta.prior
# mu.prior <- priors$mu.prior
# # alpha.prior <- c(6, 12000)
# p.prior <- priors$p.prior
# 
# # log-likelihood vector
# loglik <- vector(length=niter)
# 
# # list to store trajectories
# trajectories <- list()
# 
# # observation matrix
# W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$Observed, augmented = 0)); 
# 
# # matrix with event times, subject id and event codes. 
# # Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
# X.cur <- SIRres$trajectory
# # X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], p=probs[1], amplify = amplify, tmax=20, popsize = popsize)
# Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
# # update observation matrix
# W.cur <- updateW(W = W.cur, X = X.cur)
# 
# if(!checkpossible(X=X.cur, W=W.cur)) {
#     while(!checkpossible(X=X.cur,W=W.cur)){
#         X.cur <- initializeX(W = W.cur, b = Beta[1], mu = Mu[1], a = Alpha[1], p=probs[1], amplify = amplify, tmax=20, popsize = popsize)
#         Xcount.cur <- build_countmat(X = X.cur, popsize = popsize)
#         W.cur <- updateW(W = W.cur, X = X.cur)
#     }
# }
# 
# popirm.cur <- build_irm(Xcount = Xcount.cur, b = Beta[1], m = Mu[1], a = Alpha[1], popsize = popsize, pop = TRUE)
# pop_prob.cur <- pop_prob(Xcount = Xcount.cur, irm = popirm.cur, initdist = initdist, popsize = popsize)
# 
# k=2; j=1
# 
# subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
# 
# pathirm.cur <- build_irm(Xcount = Xcount.cur, b = Beta[k-1], m = Mu[k-1], a = Alpha[k-1], popsize = popsize, pop = FALSE)
# patheigen.cur <- irm_decomp(pathirm.cur = pathirm.cur)
# 
# accepts.k <- 0
# 
# Xother <- X.cur[X.cur[,2]!=subjects[j],]
# Xcount.other <- build_countmat(X = Xother, popsize = popsize)
# 
# path.cur <- getpath(X.cur, subjects[j])
# 
# W.other <-updateW(W = W.cur, Xcount = Xcount.other)
# 
# Xcount = Xcount.other; irm = pathirm.cur; irm.eig = patheigen.cur; W=W.other; p=probs[k-1]; b=Beta[k-1]; m=Mu[k-1]; a=Alpha[k-1]; initdist = initdist
# 
# indend <- dim(Xcount.other)[1]
# 
# numsick <- Xcount.other[,2] 
# eventtimes <- Xcount.other[,1]
# 
# Xt.fb <- array(0,dim=c(3,3,dim(W)[1]-1))
# obstimes <- W[,1]

# Second Simulation -------------------------------------------------------

### Method 1

method1.draws <- matrix(0,nrow = 200000, ncol = length(obstimes))

for(t in 1:200000){
    if(t %% 100 == 0) print(t)
    method1.draws[t,] <- sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = FALSE)
}

### Method 2

method2.draws <- matrix(0,nrow = 200000, ncol = length(obstimes))

for(t in 1:200000){
    if(t %% 100 == 0) print(t)
    
    path <- sim_one_SIR(numsick = numsick, eventtimes = eventtimes, obstimes = obstimes, b = b, m = m, initdist = initdist, returnpath = TRUE)
    
    X.cur <- updateX(X = X.cur, Xt.path = path, j=1) # upate the complete data matrix with the new trajectory
    W.cur <- updateW(W = W.cur, Xcount = Xcount) # update the observation matrix with the new trajectory
    
    W.cur[,2] <- rbinom(n = length(obstimes), size = W.cur[,3], prob = p) # draw new binomial samples from the update observation matrix
    W.other <- updateW(W = W.cur, Xcount=Xcount) # update the observation matrix with the new binomial samples
    
    method2.draws[t,] <- drawXt(Xcount = Xcount.other, irm = pathirm.cur, irm.eig = patheigen.cur, W = W.other, p = p, initdist = initdist)[,2] # draw new trajectory using new binomial samples

}

method1.draws <- data.frame(cbind(1:200000,method1.draws)); names(method1.draws) <- c("simnum", obstimes)
method1.melt <- melt(method1.draws, id = "simnum"); names(method1.melt)[2] <- "time"
method1.melt$time <- as.numeric(method1.melt$time)/10


method2.draws <- data.frame(cbind(1:200000,method2.draws)); names(method2.draws) <- c("simnum", obstimes)
method2.melt <- melt(method2.draws, id = "simnum"); names(method2.melt)[2] <- "time"
method2.melt$time <- as.numeric(method2.melt$time)/10

comp.summary <- data.frame(Method = c(rep("Method 1", length(obstimes)), rep("Method 2", length(obstimes))),
                           Time = rep(obstimes,2), 
                           Susceptible = 0,
                           Infected = 0, 
                           Recovered = 0)

for(r in 1:length(obstimes)){
    comp.summary[r,3] <- mean(method1.draws[,r+1]==1)
    comp.summary[r,4] <- mean(method1.draws[,r+1]==2)
    comp.summary[r,5] <- mean(method1.draws[,r+1]==3)
    
    comp.summary[r + length(obstimes),3] <- mean(method2.draws[,r+1]==1)
    comp.summary[r + length(obstimes),4] <- mean(method2.draws[,r+1]==2)
    comp.summary[r + length(obstimes),5] <- mean(method2.draws[,r+1]==3)
    
}

summary.melt <- melt(comp.summary, id = c("Method", "Time"))
ggplot(summary.melt, aes(x = Time, y = value, colour = variable, linetype = Method)) + geom_path() + labs(y="Probability")

# # First Simulation - 1/31/2015 --------------------------------------------------------
# 
# 
# draws <- matrix(0,nrow = 200000, ncol = length(obstimes))
# 
# for(t in 1:200000){
#     if(t%%200000 == 0) print(t)
#     
#     Xt <- drawXt(Xother = Xother, irm = pathirm.cur, irm.eig = patheigen.cur, W=W.other, p=probs[k-1], initdist = initdist)
#     
#     draws[t,] <- Xt[,2]
#     
# }
# 
# draws <- cbind(1:200000, draws)
# draws <- data.frame(draws); names(draws) <- c("simnum", obstimes)
# draws.melt <- melt(draws,id="simnum"); names(draws.melt)[2] <- "Obs_Number"
# draws.melt$time <- as.numeric(draws.melt$Obs_Number)
# 
# marginals <- data.frame(t(marginals))
# names(marginals) <- c("Susceptible", "Infected", "Recovered")
# marginals <- cbind(obstimes,marginals); marginals$obstimes <- marginals$obstimes*4
# names(marginals)[1] <- "Obs_Number"
# marginals.melt <- melt(marginals, id = "Obs_Number")
# 
# draws.summary <- marginals
# for(r in 1:length(obstimes)){
#     draws.summary[r,2] <- mean(draws[,r+1]==1)
#     draws.summary[r,3] <- mean(draws[,r+1]==2)
#     draws.summary[r,4] <- mean(draws[,r+1]==3)
#     
# }
# draws.summary.melt <- melt(draws.summary, id = "Obs_Number")
# 
# marginals.melt <- cbind("Expected", marginals.melt); names(marginals.melt)[1] <- "which_prob"
# draws.summary.melt <- cbind("Observed", draws.summary.melt);  names(draws.summary.melt)[1] <- "which_prob"
# 
# results.melt <- rbind(marginals.melt, draws.summary.melt)
# results.melt$Obs_Number <- results.melt$Obs_Number/10
# 
# ggplot(results.melt, aes(x = Obs_Number, y = value, colour = variable, linetype = which_prob)) + geom_path() + labs(y="Probability")
