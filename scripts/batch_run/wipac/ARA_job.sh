#!/bin/bash

# load in variables
key=sub_info
station=$1
run=$2
blind_dat=1
condor_run=0
not_override=0
qual_type=1
l2_data=0
no_tqdm=0

# run the reconstruction script
export HDF5_USE_FILE_LOCKING='FALSE'
source /cvmfs/ara.opensciencegrid.org/trunk/centos7/setup.sh
source /home/abishop/ara/a23/MF_filters/setup.sh
cd /home/abishop/ara/a23/MF_filters/scripts/

python3 /home/abishop/ara/a23/MF_filters/scripts/script_executor.py -k ${key} -s ${station} -r ${run} -b ${blind_dat} -c ${condor_run} -n ${not_override} -q ${qual_type} -t ${no_tqdm} -l ${l2_data}

