#!/bin/sh

# set working directory
#cd ~/augSIR_convcomp/

# qsub -o /dev/null -e /dev/null

# set local variables
set outfile = "~/convergence_comp/"rjmcmc.$arg1.$arg2.$arg3.$arg4.$arg5.$$.Rout

# execute commands
R CMD BATCH "--args $arg1 $arg2 $arg3 $arg4 $arg5" "convcomp_rjmcmc_cluster_code.R" $outfile

