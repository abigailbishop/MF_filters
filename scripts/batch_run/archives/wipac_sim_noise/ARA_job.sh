#!/bin/bash

# load in variables
setup=$1
run=$2
result=$3

# run the reconstruction script
export HDF5_USE_FILE_LOCKING='FALSE'
source /cvmfs/ara.opensciencegrid.org/trunk/centos7/setup.sh
source /home/mkim/analysis/MF_filters/setup.sh
cd /home/mkim/analysis/AraSoft/AraSim/

./AraSim ${setup} ${run} ${result}

