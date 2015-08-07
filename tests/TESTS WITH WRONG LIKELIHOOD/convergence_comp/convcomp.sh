#!/bin/sh

# set working directory
cd ~/convergence_comp/

# qsub -o /dev/null -e /dev/null

# set local variables
set outfile = "~/convergence_comp/"augSIR.$arg1._$arg2._$arg3._$arg4._$arg5._$$.Rout

# execute commands
R CMD BATCH "--args $arg1 $arg2 $arg3 $arg4 $arg5" "convcomp_cluster_code.R" $outfile

