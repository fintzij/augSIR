######################################################################
### eval.R contains commands to test functionality of the package. ###
######################################################################

require(ggplot2)
require(reshape2)

# Simulate data -----------------------------------------------------------

SIRres<-SIRsim(N = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, maxtime = 365,censusInterval=14);
SIRres = cbind(SIRres,200 - rowSums(SIRres[,2:3]))
colnames(SIRres)<-c("time","susceptible","infected","recovered")

# get data 
dat <- data.frame(SIRres); dat$infected<-rbinom(n=dim(dat)[1], size=dat$infected, prob=1)
dat.m <- melt(dat,id.vars="time")

ggplot(dat.m,aes(x=time,y=value,colour=variable))+geom_point() + theme_bw()


# initialize simulation settings
popsize <- 200 # size of the population; 
tmax <- 20 # maximum time of observation
niter <- 2000 # number of iterations in the sampler
initdist <- c(0.995,0.005,0) # initial distribution for individual infection status

# vectors for parameters
Beta <- vector(length=niter); Beta[1] <- 0.01 + runif(1,-0.005,0.005)
Mu <- vector(length = niter); Mu[1] <- 0.5 + runif(1,-0.5,0.5)
Alpha <- vector(length = niter); Alpha[1] <- 0 
probs <- vector(length = niter); #probs[1] <- 0.2 + runif(1, -0.1, 0.1)
probs[1] <-1

# vectors for parameters of distributions for beta, mu, and p. beta and mu have gamma distributions, p has beta distribution.
beta.prior <- c(12, 1200)
mu.prior <- c(6, 12)
# alpha.prior <- c(6, 12000)
# p.prior <- c(9.5, 38)

# log-likelihood vector
loglik <- vector(length=niter)

# observation matrix
W.cur <- as.matrix(data.frame(time = dat$time, sampled = dat$infected, augmented = 0)); 
if(W.cur[1,2]==0) W.cur[1,2]<-1

# matrix with event times, subject id and event codes. 
# Event codes: 1=carriage aquired, -1=carriage infected, 0=event out of time range
X.cur <- as.matrix(data.frame(time=rep(0,popsize*2), id=rep(1:popsize,each=2), event=rep(0,2*popsize)))

X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20)

# update observation matrix
W.cur <- updateW(W.cur,X.cur)

if(!checkpossible(X=X.cur, W=W.cur)) {
    while(!checkpossible(X=X.cur,W=W.cur)){
        X.cur <- initializeX(W.cur, Mu[1], probs[1], 1, tmax=20)
        W.cur <- updateW(W.cur,X.cur)
    }
}



# M-H sampler
for(k in 1:(niter-1)){
  print(k)
  # Update trajectories
  
  subjects <- sample(unique(X.cur[,2]),length(unique(X.cur[,2])),replace=TRUE)
  irm.cur <- buildirm(X.cur, b = Beta[k], m = Mu[k], a = Alpha[k]) 
  
  for(j in 1:length(subjects)){
      print(j)
      Xother <- X.cur[X.cur[,2]!=subjects[j],]; 
      if(checkpossible(Xother) == TRUE){
          
          path.cur <- getpath(X.cur, subjects[j])
          
          W.other <-updateW(W.cur, Xother)
          irm.other <- buildirm(Xother, b=Beta[k], m = Mu[k], a=Alpha[k])
          Xt <- drawXt(Xother = Xother, irm = irm.other, W=W.other, p=probs[k], b=Beta[k], m=Mu[k], a=Alpha[k], initdist = initdist)
          
          path.new<- drawpath(Xt, Xother, irm.other, tmax)
          
          X.new <- updateX(X.cur,path.new,subjects[j]); path.new <- getpath(X.new,subjects[j])
          irm.new <- buildirm(X.new, b = Beta[k], m = Mu[k], a = Alpha[k])
          
          a.prob <- pop_prob(X.new, irm = irm.new) - pop_prob(X.cur, irm = irm.cur) + 
              path_prob(path.cur, Xother, irm.other, initdist, tmax) - path_prob(path.new, Xother, irm.other, initdist, tmax)
          
          if(min(a.prob, 0) > log(runif(1))) {
              X.cur <- X.new
              W.cur <- updateW(W.cur,X.cur)
              irm.cur <- irm.new
          }
          
      } else next
  }
    
  # Update observation matrix
  W <- updateW(W, X)
  
  # draw new binomial sampling probability parameter
#   probs[k+1] <- rbeta(1,p.prior[1] + sum(W[,2]), p.prior[2] + sum(W[,3]-W[,2]))
  probs[k+1] <- 1
  
  # draw new rate parameters 
#   Beta[k+1] <- rgamma(1, shape = (beta.prior[1] + sum(X[,3]==1)), 
#                              rate = beta.prior[2] + sum((popsize - cumsum(X[,3]==1))*(cumsum(X[,3]==1) - cumsum(X[,3]==-1))*c(0,X[2:dim(X)[1],1] - X[1:(dim(X)[1]-1),1])*(X[,3]==1)))
#   
  Mu[k+1] <- rgamma(1, shape = mu.prior[1] + sum(X[,3]==-1),
                           rate = mu.prior[2] + sum((cumsum(X[,3]==1) - cumsum(X[,3]==-1))*c(0,X[2:dim(X)[1],1] - X[1:(dim(X)[1]-1)])*(X[,3]==-1)))
  
#   Alpha[k+1] <- rgamma(1, shape = (alpha.prior[1] + sum(X[,3]==1)),
#                        rate = alpha.prior[2] + sum((popsize - cumsum(X[,3]==1))*c(0,X[2:dim(X)[1],1] - X[1:(dim(X)[1]-1),1])*(X[,3]==1)))
  Alpha[k+1] <- 0

  loglik[k] <- calc_loglike(X, W, Beta[k+1], Mu[k+1], Alpha[k+1], probs[k+1])  
  
#   trajectories <- data.frame(X); epidemic <- data.frame(W)
#   trajectories$id <- factor(trajectories$id, levels = unique(as.factor(trajectories$id)))
#   trajectories$event <- factor(trajectories$event)
#   trajects <- ggplot(subset(trajectories,time!=0),aes(y=id,x=time,group=id,colour=as.factor(event))) + geom_point()+ geom_line() + geom_vline(xintercept = data.frame(W)$time,alpha=0.2)+
#     scale_colour_discrete(name="Event", labels=c("Infection Cleared", "Infection Acquired"))
#   
#   curve <- ggplot(epidemic,aes(x=time,y=augmented))+geom_line() + theme_bw()
#   
#   print(grid.arrange(trajects, curve, ncol=2))
  
}

# Plots -------------------------------------------------------------------

trajectories <- data.frame(X)
trajectories$id <- factor(trajectories$id, levels = unique(as.factor(trajectories$id)))
trajectories$event <- factor(trajectories$event)
ggplot(subset(trajectories,time!=0),aes(y=id,x=time,group=id,colour=as.factor(event))) + geom_point()+ geom_line() + geom_vline(xintercept = data.frame(W)$time,alpha=0.2)
  scale_colour_discrete(name="Event", labels=c("Infection Cleared", "Infection Acquired"))
