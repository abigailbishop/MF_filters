import numpy as np
import os, sys
from glob import glob
import h5py
from tqdm import tqdm

curr_path = os.getcwd()
sys.path.append(curr_path+'/../')
from tools.ara_run_manager import file_sorter
from tools.ara_utility import size_checker
from tools.ara_known_issue import known_issue_loader

Station = int(sys.argv[1])
Blind = int(sys.argv[2])
if Blind == 1: b_name = '_full'
else: b_name = ''
print(b_name)

# sort
m_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/reco_ele{b_name}/'
ml_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/reco_ele_lite{b_name}/'

r_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/Hist/'
file_name = f'Info_Summary_A{Station}.h5'
hf = h5py.File(r_path + file_name, 'r')
runs = hf['runs'][:]
num_runs = len(runs)
num_evts = hf['num_evts'][:]
del hf, r_path, file_name

ang_num = np.arange(2, dtype = int)
ang_len = len(ang_num)
pol_num = np.arange(2, dtype = int)
pol_len = len(pol_num)
theta = 90 - np.linspace(0.5, 179.5, 179 + 1)
the_len = len(theta)
z = np.sin(np.radians(theta)) * float(41)

for r in tqdm(range(num_runs)):
    
  #if r <10:

    m_name = f'{m_path}reco_ele{b_name}_A{Station}_R{runs[r]}.h5'
    ml_name = f'{ml_path}reco_ele_lite{b_name}_A{Station}_R{runs[r]}.h5'

    hf = h5py.File(ml_name, 'r+')
    mf_list = list(hf)
    try:
        mf_lite_idx = mf_list.index('trig_type') # or coord_cal
    except ValueError:
        mf_lite_idx = -1
    del mf_list

    if mf_lite_idx != -1:
        print(f'{m_name} already has mf_lite! move on!')
        pass
    else:
        try:
            evt_num = np.arange(num_evts[r], dtype = int)

            hf_l = h5py.File(m_name, 'r')
            trig_type = hf_l['trig_type'][::] # evt
            del hf_l

            hf.create_dataset('trig_type', data=trig_type, compression="gzip", compression_opts=9)
        except Exception as e:
            print(e)
            print(f'{m_name}')
            print(f'{ml_name}')
    hf.close()
    
print('done!')





