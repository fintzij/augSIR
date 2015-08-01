library(batch)
k=j=0
for (popsize in c(50, 150, 300)) {
    for (censusInterval in c(0.05, 0.2)) {
        for (samp_prob in c(0.1,0.4)) {
            for (resample_prop in c(1)) {
                for (initialization_num in 1:3) {
                    for(R0 in c(4, 8)){
                        k = k+1; j = j+1
                        rbatch("./convcomp_augSIR_reparMH.R", seed = k, popsize = popsize, censusInterval = censusInterval, samp_prob = samp_prob, resample_prop = resample_prop, initialization_num = initialization_num, R0 = R0)
                    }
                }
            }
        }
    }
}

# rbatch("./convcomp_cluster_code.R", seed = 1, popsize = 25, censusInterval = 0.2, samp_prob = 0.5, resample_prop = 1, initialization_num = 1)
# rbatch("./convcomp_rjmcmc_statsclust.R", seed = 2, popsize = 25, censusInterval = 0.2, samp_prob = 0.5, resample_prop = 1, initialization_num = 1)
