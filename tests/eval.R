######################################################################
### eval.R contains commands to test functionality of the package. ###
######################################################################
library(profr)
library(proftools)
# augSIR simulation -------------------------------------------------------

# simulation
SIRres<-SIRsim(popsize = 200, S0 = 199, I0 = 1, b = 0.01, mu=.5, a=0, tmax = 20, censusInterval=0.25, sampprob = 0.25)

# get data 
dat <- SIRres[,c(1:2)]
dat.m <- melt(SIRres,id.vars="time")

ggplot(dat.m, aes(x=time, y=value, colour=variable)) + geom_point() + theme_bw()

sim.settings <- list(popsize = 200,
                     tmax = 20,
                     niter = 500,
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

#   trajectories <- data.frame(X); epidemic <- data.frame(W)
#   trajectories$id <- factor(trajectories$id, levels = unique(as.factor(trajectories$id)))
#   trajectories$event <- factor(trajectories$event)
#   trajects <- ggplot(subset(trajectories,time!=0),aes(y=id,x=time,group=id,colour=as.factor(event))) + geom_point()+ geom_line() + geom_vline(xintercept = data.frame(W)$time,alpha=0.2)+
#     scale_colour_discrete(name="Event", labels=c("Infection Cleared", "Infection Acquired"))
#   
#   curve <- ggplot(epidemic,aes(x=time,y=augmented))+geom_line() + theme_bw()
#   
#   print(grid.arrange(trajects, curve, ncol=2))