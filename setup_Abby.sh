echo "before :("
echo "ARA_UTIL_INSTALL_DIR = "${ARA_UTIL_INSTALL_DIR}
echo "ARA_ROOT_DIR = "${ARA_ROOT_DIR}
echo "ARA_SIM_DIR = "${ARA_SIM_DIR}
echo "GTRS_DIR = "${GTRS_DIR}
echo "LD_LIBRARY_PATH = "${LD_LIBRARY_PATH}

export HDF5_USE_FILE_LOCKING='FALSE'
source /cvmfs/ara.opensciencegrid.org/trunk/centos7/setup.sh

export RAW_PATH=/data/exp/ARA
export OUTPUT_PATH=/data/ana/ARA
export ARA_UTIL_INSTALL_DIR=/home/mkim/analysis/AraSoft/AraUtil
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ARA_UTIL_INSTALL_DIR
export DYLD_LIBRARY_PATH=$DYLA_LIBRARY_PATH:$ARA_UTIL_INSTALL_DIR
export PATH=$PATH:$ARA_UTIL_INSTALL_DIR

export ARA_ROOT_DIR=/home/mkim/analysis/AraSoft/AraRoot
export ARA_SIM_DIR=/home/mkim/analysis/AraSoft/AraSim
export GTRS_DIR=/home/mkim/analysis/AraSoft/GruanToolRs92
export LD_LIBRARY_PATH=${ARA_UTIL_INSTALL_DIR}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${ARA_SIM_DIR}:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${GTRS_DIR}:$LD_LIBRARY_PATH

export MFFILTERS_PATH=/home/abishop/ara/a23/MF_filters

echo "after  :)"
echo "ARA_UTIL_INSTALL_DIR = "${ARA_UTIL_INSTALL_DIR}
echo "ARA_ROOT_DIR = "${ARA_ROOT_DIR}
echo "ARA_SIM_DIR = "${ARA_SIM_DIR}
echo "GTRS_DIR = "${GTRS_DIR}
echo "LD_LIBRARY_PATH = "${LD_LIBRARY_PATH}


