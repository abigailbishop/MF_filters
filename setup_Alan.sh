echo "before :("
echo "ARA_UTIL_INSTALL_DIR = "${ARA_UTIL_INSTALL_DIR}
echo "ARA_ROOT_DIR = "${ARA_ROOT_DIR}
echo "ARA_SIM_DIR = "${ARA_SIM_DIR}
echo "GTRS_DIR = "${GTRS_DIR}
echo "LD_LIBRARY_PATH = "${LD_LIBRARY_PATH}

export HDF5_USE_FILE_LOCKING='FALSE'
source /cvmfs/ara.opensciencegrid.org/trunk/centos7/setup.sh

export RAW_PATH=/data/exp/ARA
export OUTPUT_PATH=/data/user/alansalgo/analysisA23

#export ARA_UTIL_INSTALL_DIR=/cvmfs/ara.opensciencegrid.org/trunk/RHEL_7_x86_64/ara_build
export ARA_UTIL_INSTALL_DIR=/home/mkim/analysis/AraSoft/AraUtil
export ARA_ROOT_DIR=/cvmfs/ara.opensciencegrid.org/trunk/RHEL_7_x86_64/source/AraRoot
export ARA_SIM_DIR=/home/source/AraSim
export GTRS_DIR=
export LD_LIBRARY_PATH=${ARA_UTIL_INSTALL_DIR}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${ARA_SIM_DIR}:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${GTRS_DIR}:$LD_LIBRARY_PATH

echo "after  :)"
echo "ARA_UTIL_INSTALL_DIR = "${ARA_UTIL_INSTALL_DIR}
echo "ARA_ROOT_DIR = "${ARA_ROOT_DIR}
echo "ARA_SIM_DIR = "${ARA_SIM_DIR}
echo "GTRS_DIR = "${GTRS_DIR}
echo "LD_LIBRARY_PATH = "${LD_LIBRARY_PATH}


