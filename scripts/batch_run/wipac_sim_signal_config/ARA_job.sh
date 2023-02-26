#!/bin/bash

# load in variables
key=signal
st=$1
run=$2
fla=$3
sim_run=$4
user_path=/misc/disk19/users/
if [ -d "$user_path" ]; then
    echo "There is ${user_path}"
else
    echo "Switch to /data/user/"
    user_path=/data/user/
fi
setup=${user_path}mkim/OMF_filter/ARA0${st}/sim_${key}_setup_full/${key}_F${fla}_A${st}_R${run}.txt
result=${user_path}mkim/OMF_filter/ARA0${st}/sim_${key}_full

if [ -d "$result" ]; then
    echo "There is ${result}"
else
    echo "Make ${result}"
    mkdir ${result}
fi

# run the reconstruction script
export HDF5_USE_FILE_LOCKING='FALSE'
source /cvmfs/ara.opensciencegrid.org/trunk/centos7/setup.sh
source /home/mkim/analysis/MF_filters/setup.sh
cd /home/mkim/analysis/AraSoft/AraSim/

./AraSim ${setup} ${sim_run} $TMPDIR

# at the end, move the results back
mv $TMPDIR/*.root ${result}

