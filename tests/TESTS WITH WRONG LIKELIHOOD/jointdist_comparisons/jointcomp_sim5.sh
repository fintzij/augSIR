#!/bin/sh

# set working directory
#cd ~/augSIR_rjmcmc_comp

# qsub -o /dev/null -e /dev/null

# set local variables
set outfile = "~/augSIR_jointcomp/".sim5.$$.Rout

# execute commands
R CMD BATCH "jointcomp_sim5.R" $outfile