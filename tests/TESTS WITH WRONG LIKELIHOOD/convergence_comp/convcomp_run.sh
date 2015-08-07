#!/bin/sh

# set working directory
cd ~/convergence_comp

# run qsubs 
# arg1 is population size - 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000
# arg2 is censusInterval/number of observations - 0.2
# arg3 is binomial sampling probability - 0.5
# arg4 is fraction of subjects to resample - 1
# arg5 is initialization number - 1, 2, 3
qsub -cwd -v arg1=25 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp.sh
qsub -cwd -v arg1=50 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp.sh
qsub -cwd -v arg1=100 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp.sh
qsub -cwd -v arg1=250 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp.sh
qsub -cwd -v arg1=500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp.sh
qsub -cwd -v arg1=1000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp.sh
qsub -cwd -v arg1=2500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp.sh
qsub -cwd -v arg1=5000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp.sh
qsub -cwd -v arg1=10000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp.sh
qsub -cwd -v arg1=25 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp.sh
qsub -cwd -v arg1=50 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp.sh
qsub -cwd -v arg1=100 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp.sh
qsub -cwd -v arg1=250 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp.sh
qsub -cwd -v arg1=500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp.sh
qsub -cwd -v arg1=1000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp.sh
qsub -cwd -v arg1=2500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp.sh
qsub -cwd -v arg1=5000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp.sh
qsub -cwd -v arg1=10000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp.sh
qsub -cwd -v arg1=25 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp.sh
qsub -cwd -v arg1=50 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp.sh
qsub -cwd -v arg1=100 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp.sh
qsub -cwd -v arg1=250 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp.sh
qsub -cwd -v arg1=500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp.sh
qsub -cwd -v arg1=1000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp.sh
qsub -cwd -v arg1=2500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp.sh
qsub -cwd -v arg1=5000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp.sh
qsub -cwd -v arg1=10000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp.sh
qsub -cwd -v arg1=25 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp_rjmcmc.sh
qsub -cwd -v arg1=50 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp_rjmcmc.sh
qsub -cwd -v arg1=100 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp_rjmcmc.sh
qsub -cwd -v arg1=250 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp_rjmcmc.sh
qsub -cwd -v arg1=500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp_rjmcmc.sh
qsub -cwd -v arg1=1000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp_rjmcmc.sh
qsub -cwd -v arg1=2500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp_rjmcmc.sh
qsub -cwd -v arg1=5000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp_rjmcmc.sh
qsub -cwd -v arg1=10000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=1 convcomp_rjmcmc.sh
qsub -cwd -v arg1=25 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp_rjmcmc.sh
qsub -cwd -v arg1=50 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp_rjmcmc.sh
qsub -cwd -v arg1=100 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp_rjmcmc.sh
qsub -cwd -v arg1=250 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp_rjmcmc.sh
qsub -cwd -v arg1=500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp_rjmcmc.sh
qsub -cwd -v arg1=1000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp_rjmcmc.sh
qsub -cwd -v arg1=2500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp_rjmcmc.sh
qsub -cwd -v arg1=5000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp_rjmcmc.sh
qsub -cwd -v arg1=10000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=2 convcomp_rjmcmc.sh
qsub -cwd -v arg1=25 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp_rjmcmc.sh
qsub -cwd -v arg1=50 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp_rjmcmc.sh
qsub -cwd -v arg1=100 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp_rjmcmc.sh
qsub -cwd -v arg1=250 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp_rjmcmc.sh
qsub -cwd -v arg1=500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp_rjmcmc.sh
qsub -cwd -v arg1=1000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp_rjmcmc.sh
qsub -cwd -v arg1=2500 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp_rjmcmc.sh
qsub -cwd -v arg1=5000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp_rjmcmc.sh
qsub -cwd -v arg1=10000 -v arg2=0.1 -v arg3=0.2 -v arg4=0.8 -v arg5=3 convcomp_rjmcmc.sh
