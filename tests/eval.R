######################################################################
### eval.R contains commands to test functionality of the package. ###
######################################################################

# augSIR simulation -------------------------------------------------------


# simulation
SIRres<-SIRsim(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 20, censusInterval=0.25, sampprob = 0.25)

# get data 
dat <- SIRres[,c(1:2)]
dat.m <- melt(SIRres,id.vars="time")

ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + theme_bw()

sim.settings <- list(popsize = 200,
                     tmax = 20,
                     niter = 10000,
                     amplify = 5,
                     initdist = c(0.995, 0.005, 0))

inits <- list(beta.init = 0.01 + runif(1,-0.005, 0.005),
              mu.init = 0.5 + runif(1, -0.05, 0.05),
              alpha.init = 0, 
              probs.init = 0.2 + runif(1,-0.1, 0.1))

priors <- list(beta.prior = c(.01, 1),
               mu.prior = c(1, 2),
               alpha.prior = NULL,
               p.prior = c(0.022, 0.084))

# run sampler

results <- augSIR(dat, sim.settings, priors, inits)

#profile
Rprof("~/School/UW/Year 3 +/Dissertation/Code/augSIR/tests/profile/augSIRprofile.out")
results <- augSIR(dat, sim.settings, priors, inits)
Rprof()
summaryRprof("~/School/UW/Year 3 +/Dissertation/Code/augSIR/tests/profile/augSIRprofile.out")
# results.prof <- lineprof(augSIR(dat, sim.settings, priors, inits))
plotProfileCallGraph(readProfileData("~/School/UW/Year 3 +/Dissertation/Code/augSIR/tests/profile/augSIRprofile.out"),score = "total")

# Plots -------------------------------------------------------------------

trajectories <- data.frame(X)
trajectories$id <- factor(trajectories$id, levels = unique(as.factor(trajectories$id)))
trajectories$event <- factor(trajectories$event)
ggplot(subset(trajectories,time!=0),aes(y=id,x=time,group=id,colour=as.factor(event))) + geom_point()+ geom_line() + geom_vline(xintercept = data.frame(W)$time,alpha=0.2)
  scale_colour_discrete(name="Event", labels=c("Infection Cleared", "Infection Acquired"))


# Testing various simulation functions-----------------------------------------------------------

# SIRres<-SIRsim(N = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, maxtime = 20,censusInterval=.01)
# SIRres = cbind(SIRres,200 - rowSums(SIRres[,2:3]))
# colnames(SIRres)<-c("time","susceptible","infected","recovered")
# 
# # get data 
# dat <- data.frame(SIRres); dat$infected<-rbinom(n=dim(dat)[1], size=dat$infected, prob=0.2)
# dat.m <- melt(dat,id.vars="time")
# 
# ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + theme_bw()
# 
# SIRres<-SIRsim3(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 20,censusInterval=0.5, prob=0.2)
# 
# # get data 
# dat <- data.frame(SIRres); dat$Binomial.Count<-rbinom(n=dim(dat)[1], size=dat$Truth, prob = 0.2)
# dat.m <- melt(dat,id.vars="time"); dat.m$variable <- factor(dat.m$variable, levels = c("Truth", "Observed", "Binomial.Count"))
# 
# ggplot(dat.m,aes(x=time,y=value,colour=variable))+geom_line() + theme_bw()
# 
# # check sampling method - binomial vs. individual
# 
# SIRsims <- data.frame(simnum=1, SIRres)
# for(k in 2:2000){
#     print(k)
#     SIRres<-SIRsim3(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 20,censusInterval=0.5, prob=0.2)
#     
#     if(max(SIRres$Truth)==1){
#         keep.going <- TRUE
#         while(keep.going == TRUE){
#             SIRres<-SIRsim3(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 20,censusInterval=0.5, prob=0.2)
#             if(max(SIRres$Truth)>1) keep.going <-FALSE
#         }
#     }
#     
#     SIRsims <- rbind(SIRsims,data.frame(simnum = k, SIRres))
# }
# 
# SIRsims$Binomial.Count <- rbinom(n=dim(SIRsims)[1], size = SIRsims$Truth, prob = 0.2)
# SIRsims.m <- melt(SIRsims, id.vars = c("time", "simnum")); SIRsims.m$variable <- factor(SIRsims.m$variable, levels = c("Truth", "Observed", "Binomial.Count"))
# 
# ggplot(SIRsims.m, aes(x = factor(time), y = value, fill = variable)) + geom_boxplot(outlier.shape=NA) + labs(x = "Observation Time", y="Number of Infecteds") 
# 
