#!/bin/sh

# set working directory
cd ~/augSIR_jointcomp

# run qsubs
qsub -cwd jointcomp_sim1.sh
qsub -cwd jointcomp_sim2.sh
qsub -cwd jointcomp_sim3.sh
qsub -cwd jointcomp_sim4.sh
qsub -cwd jointcomp_sim5.sh
qsub -cwd jointcomp_sim6.sh
