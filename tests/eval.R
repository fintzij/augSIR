######################################################################
### eval.R contains commands to test functionality of the package. ###
######################################################################

# Simulate data -----------------------------------------------------------

# SIRres<-SIRsim(N = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, maxtime = 20,censusInterval=1, prob=0.2)
# SIRres = cbind(SIRres,200 - rowSums(SIRres[,2:3]))
# colnames(SIRres)<-c("time","susceptible","infected","recovered")
# 
# # get data 
# dat <- data.frame(SIRres); dat$infected<-rbinom(n=dim(dat)[1], size=dat$infected)
# dat.m <- melt(dat,id.vars="time")

SIRres<-SIRsim2(popsize = 200, S0 = 199, I0 = 1, b = 0.002, mu=.5, a=0, tmax = 15,censusInterval=0.01, prob=0.2)
colnames(SIRres)<-c("time","Observed","Truth")

# get data 
dat <- data.frame(SIRres); dat$Binomial.Count<-rbinom(n=dim(dat)[1], size=dat$Truth, prob = 0.2)
dat.m <- melt(dat,id.vars="time"); dat.m$variable <- factor(dat.m$variable, levels = c("Truth", "Observed", "Binomial.Count"))

ggplot(dat.m,aes(x=time,y=value,colour=variable))+geom_line() + theme_bw()

sim.settings <- list(popsize = 200,
                     tmax = 20,
                     niter = 2000,
                     amplify = 2000,
                     initdist = c(0.995, 0.005, 0))

inits <- list(beta.init = 0.01 + runif(1,-0.005, 0.005),
              mu.init = 0.5 + runif(1, -0.5, 0.5),
              alpha.init = 0, 
              probs.init = 0.2 + runif(1,-0.1, 0.1))

priors <- list(beta.prior = c(6e-4, 0.05),
               mu.prior = c(0.13, 0.24),
               alpha.prior = NULL,
               p.prior = c(0.022, 0.084))

# run sampler

results <- augSIR(dat, sim.settings, priors, inits)
# results.prof <- lineprof(augSIR(dat, sim.settings, priors, inits))


# Plots -------------------------------------------------------------------

trajectories <- data.frame(X)
trajectories$id <- factor(trajectories$id, levels = unique(as.factor(trajectories$id)))
trajectories$event <- factor(trajectories$event)
ggplot(subset(trajectories,time!=0),aes(y=id,x=time,group=id,colour=as.factor(event))) + geom_point()+ geom_line() + geom_vline(xintercept = data.frame(W)$time,alpha=0.2)
  scale_colour_discrete(name="Event", labels=c("Infection Cleared", "Infection Acquired"))
