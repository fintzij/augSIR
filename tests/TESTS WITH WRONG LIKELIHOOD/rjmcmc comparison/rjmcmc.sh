#!/bin/sh

# set working directory
#cd ~/augSIR_rjmcmc_comp

# qsub -o /dev/null -e /dev/null

# set local variables
set outfile = "~/augSIR_rjmcmc_comp/".$arg1.$$.Rout

# execute commands
R CMD BATCH "--args $arg1" "rjmcmc_cluster_code.R" $outfile

