library(pomp)

# single step simulator for SIR via Gillespie

SIR_proc_sim <- function(xt, params, ...){
    
    # unpack the rate parameters
    b <- params["b"]
    m <- params["m"]
    
    # unpack the number of susceptibles and infecteds at time tau
    St <- xt["St"]
    It <- xt["It"]
    
    # unpack the current time
    tau <- xt["tau"]
    
    # calculate the rates and hazards
    infec.rate <- b * It * St
    recov.rate <- m * It
    hazard <- infec.rate + recov.rate
    
    # sample the next event time
    delta.tau <- rexp(1, rate = hazard)
    
    # sample whether an infection or recovery occurs
    event <- sample.int(2, prob = c(infec.rate, recov.rate))
    
    # return new state and time
    xnew <- c(St = unname(ifelse(event == 1, St - 1, St)), It = unname(ifelse(event == 1, It + 1, It - 1)), tau = unname(tau+delta.tau))
    
    return(xnew)
    
}



# define sir rates function

sir_rates <- function(j, x, t, params, ...){
    
    
    
}


# define pomp object
sir_pomp <- pomp(
                data = data.frame(
                        time = seq(0,25,by=0.25), 
                        Y = NA),
                times = "time",
                statenames = c("S", "I", "R"),
                rprocess = gillespie.sim(
                    rate.fun = "_sir_rates",
                    PACKAGE = "pomp",
                    v = cbind(
                        infection = c(-1, 1, 0),
                        recovery = c(0, -1, 1)
                        ),
                    d = cbind(
                        infection = c(1, 1, 0),
                        recovery = c(0, 1, 0)
                        )
                    ),
                measurement.model = Y ~ binom(size = I, prob = p),
                comp.names=c("S","I","R"),
                ic.names=c("S.0","I.0","R.0"),
                initializer=function(params, t0, comp.names, ic.names, ...) {
                    x0 <- numeric(5)
                    names(x0) <- c("S","I","R")
                    fracs <- params[ic.names]
                    x0["N"] <- params["pop"]
                    x0[comp.names] <- round(params["pop"]*fracs/sum(fracs))
                    x0
                }
                )