library(pomp)

# single step simulator for SIR via Gillespie
# 
# SIR_proc_sim <- function(xt, params, ...){
#     
#     # unpack the rate parameters
#     b <- params["b"]
#     m <- params["m"]
#     
#     # unpack the number of susceptibles and infecteds at time tau
#     St <- xt["St"]
#     It <- xt["It"]
#     
#     # unpack the current time
#     tau <- xt["tau"]
#     
#     # calculate the rates and hazards
#     infec.rate <- b * It * St
#     recov.rate <- m * It
#     hazard <- infec.rate + recov.rate
#     
#     # sample the next event time
#     delta.tau <- rexp(1, rate = hazard)
#     
#     # sample whether an infection or recovery occurs
#     event <- sample.int(2, prob = c(infec.rate, recov.rate))
#     
#     # return new state and time
#     xnew <- c(St = unname(ifelse(event == 1, St - 1, St)), It = unname(ifelse(event == 1, It + 1, It - 1)), tau = unname(tau+delta.tau))
#     
#     return(xnew)
#     
# }
# 


pomp(
    data=read.csv2(text=dat),
    times="time",
    t0=0,
    params=c(
        gamma=24,mu=1/70,iota=0.1,
        beta1=330,beta2=410,beta3=490,
        rho=0.1,
        S.0=0.05,I.0=1e-4,R.0=0.95,
        pop=1000000,
        beta.sd=0
    ),
    rprocess=gillespie.sim(
        rate.fun="_sir_rates",
        PACKAGE="pomp",
        v=cbind(
            birth=c(1,0,0,1,0),
            sdeath=c(-1,0,0,-1,0),
            infection=c(-1,1,0,0,0),
            ideath=c(0,-1,0,-1,0),
            recovery=c(0,-1,1,0,1),
            rdeath=c(0,0,-1,-1,0)
        ),
        d=cbind(
            birth=c(0,0,0,1,0),
            sdeath=c(1,0,0,0,0),
            infection=c(1,1,0,1,0),
            ideath=c(0,1,0,0,0),
            recovery=c(0,1,0,0,0),
            rdeath=c(0,0,1,0,0)
        )
    ),
    skeleton.type="vectorfield",
    skeleton="_sir_ODE",
    measurement.model=reports~binom(size=cases,prob=rho),
    PACKAGE="pomp",
    obsnames = c("reports"),
    statenames=c("S","I","R","N","cases"),
    paramnames=c(
        "gamma","mu","iota",
        "beta1","beta.sd","pop","rho",
        "S.0","I.0","R.0"
    ),
    zeronames=c("cases"),
    comp.names=c("S","I","R"),
    ic.names=c("S.0","I.0","R.0"),
    fromEstimationScale="_sir_par_trans",
    toEstimationScale="_sir_par_untrans",
    nbasis=3L,
    degree=3L,
    period=1.0,
    initializer=function(params, t0, comp.names, ic.names, ...) {
        x0 <- numeric(5)
        names(x0) <- c("S","I","R","N","cases")
        fracs <- params[ic.names]
        x0["N"] <- params["pop"]
        x0[comp.names] <- round(params["pop"]*fracs/sum(fracs))
        x0
    }
) -> gillespie.sir

## originally, the data were created via:
## simulate(po,nsim=1,seed=1165270654L) -> gillespie.sir

c("gillespie.sir")
