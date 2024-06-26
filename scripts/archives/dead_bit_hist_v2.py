import numpy as np
import os, sys
import re
from glob import glob
import h5py
from tqdm import tqdm

curr_path = os.getcwd()
sys.path.append(curr_path+'/../')
from tools.run import file_sorter
from tools.run import bin_range_maker
from tools.utility import size_checker
from tools.ara_known_issue import known_issue_loader

Station = int(sys.argv[1])

knwon_issue = known_issue_loader(Station)
bad_runs = knwon_issue.get_knwon_bad_run()

# sort
d_path = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{Station}/dead_bit/*'
print(d_path)
d_list, d_run_tot, d_run_range = file_sorter(d_path)
del d_run_range

# config array
config_arr = []
run_arr = []
dead_bit = []
dda_volt = []
tda_volt = []
atri_volt = []
atri_curr = []
trig_ratio = []
trig_ratio_wo_bv = []
trig_ratio_wo_cal = []
is_sensor = []
num_evts = []
pre_cut = []
evt_type = np.full((3), 0, dtype = int)
rf_evt_type = np.full((3), 0, dtype = int)
bias_rf_evt_type = np.full((3), 0, dtype = int)
bias_only_rf_evt_type = np.full((3), 0, dtype = int)


for r in tqdm(range(len(d_run_tot))):

 #if r <10:

    hf = h5py.File(d_list[r], 'r')
    trig_type = hf['trig_type'][:]
    rf_count = np.count_nonzero(trig_type == 0)
    evt_type[0] += rf_count
    evt_type[1] += np.count_nonzero(trig_type == 1)
    evt_type[2] += np.count_nonzero(trig_type == 2)
    rf_evt_type[0] += rf_count
    qual_cut = hf['pre_qual_cut'][:]

    #debug
    bias_qual_cut = (qual_cut[:,-3] + qual_cut[:,-2]).astype(int)
    bias_rf_evt_type[0] += rf_count
    bias_rf_evt_type[1] += np.count_nonzero(np.logical_and(bias_qual_cut == 0, trig_type == 0))
    bias_rf_evt_type[2] += np.count_nonzero(np.logical_and(bias_qual_cut != 0, trig_type == 0))
    del bias_qual_cut

    bias_only_qual_cut = (qual_cut[:,-3]).astype(int)
    bias_only_rf_evt_type[0] += rf_count
    bias_only_rf_evt_type[1] += np.count_nonzero(np.logical_and(bias_only_qual_cut == 0, trig_type == 0))
    bias_only_rf_evt_type[2] += np.count_nonzero(np.logical_and(bias_only_qual_cut != 0, trig_type == 0))
    del bias_only_qual_cut

    if d_run_tot[r] in bad_runs:
        print('bad run:', d_list[r], d_run_tot[r])
        rf_evt_type[2] += rf_count
        continue
    else:
        rf_evt_type[1] += len(hf['clean_rf_evt'][:]) 
     
        qual_cut_temp = np.copy(qual_cut)
        qual_cut_temp[:, -1] = 0
        qual_cut_sum = np.nansum(qual_cut_temp, axis = 1)
        del qual_cut_temp 

        rf_cut = np.logical_and(qual_cut_sum != 0, trig_type == 0)
        rf_cut_count = np.count_nonzero(rf_cut)
        rf_evt_type[2] += rf_cut_count
    del trig_type, rf_count

    #if rf_evt_type[0] != (rf_evt_type[1] + rf_evt_type[2]):
    #    print(rf_evt_type[0], rf_evt_type[1], rf_evt_type[2])
    #    sys.exit(1)
    
    trig_ratio.append(hf['trig_ratio'][:])
    trig_ratio_wo_bv.append(hf['trig_ratio_wo_bv'][:])
    trig_ratio_wo_cal.append(hf['trig_ratio_wo_cal'][:])

    config = hf['config'][2]
    config_arr.append(config)

    run_arr.append(d_run_tot[r])

    dead_bit.append(hf['dead_bit_hist'][:])

    is_sensor.append(hf['is_sensor'][:])

    dda_volt.append(hf['dda_volt_hist'][:])
    tda_volt.append(hf['tda_volt_hist'][:])
    atri_volt.append(hf['atri_volt_hist'][:])
    atri_curr.append(hf['atri_curr_hist'][:])

    qual_cut_count = np.count_nonzero(qual_cut, axis = 0)
    pre_cut.append(qual_cut_count)

    num_evts_run = len(hf['evt_num'][:])
    num_evts.append(num_evts_run)

    del hf

config_arr = np.asarray(config_arr)
run_arr = np.asarray(run_arr)

print(config_arr.shape)
print(run_arr.shape)

dead_bit = np.asarray(dead_bit)
dda_volt = np.asarray(dda_volt)
tda_volt = np.asarray(tda_volt)
atri_volt = np.asarray(atri_volt)
atri_curr = np.asarray(atri_curr)
trig_ratio = np.asarray(trig_ratio)
trig_ratio_wo_bv = np.asarray(trig_ratio_wo_bv)
trig_ratio_wo_cal = np.asarray(trig_ratio_wo_cal)
is_sensor = np.asarray(is_sensor)
num_evts = np.asarray(num_evts)
pre_cut = np.asarray(pre_cut)

print(dead_bit.shape)
print(dda_volt.shape)
print(tda_volt.shape)
print(atri_volt.shape)
print(atri_curr.shape)
print(trig_ratio.shape)
print(trig_ratio_wo_bv.shape)
print(trig_ratio_wo_cal.shape)
print(is_sensor.shape)
print(num_evts.shape)
print(pre_cut.shape)
print(evt_type)
print(rf_evt_type)
print(bias_rf_evt_type)
print(bias_only_rf_evt_type)

path = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{Station}/Hist/'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)
file_name = f'Dead_Bit_RF_wo_Bad_Runs_v14_A{Station}.h5'
hf = h5py.File(file_name, 'w')

hf.create_dataset('config_arr', data=config_arr, compression="gzip", compression_opts=9)
hf.create_dataset('run_arr', data=run_arr, compression="gzip", compression_opts=9)
hf.create_dataset('dead_bit', data=dead_bit, compression="gzip", compression_opts=9)
hf.create_dataset('dda_volt', data=dda_volt, compression="gzip", compression_opts=9)
hf.create_dataset('tda_volt', data=tda_volt, compression="gzip", compression_opts=9)
hf.create_dataset('atri_volt', data=atri_volt, compression="gzip", compression_opts=9)
hf.create_dataset('atri_curr', data=atri_curr, compression="gzip", compression_opts=9)
hf.create_dataset('trig_ratio', data=trig_ratio, compression="gzip", compression_opts=9)
hf.create_dataset('trig_ratio_wo_bv', data=trig_ratio_wo_bv, compression="gzip", compression_opts=9)
hf.create_dataset('trig_ratio_wo_cal', data=trig_ratio_wo_cal, compression="gzip", compression_opts=9)
hf.create_dataset('is_sensor', data=is_sensor, compression="gzip", compression_opts=9)
hf.create_dataset('num_evts', data=num_evts, compression="gzip", compression_opts=9)
hf.create_dataset('pre_cut', data=pre_cut, compression="gzip", compression_opts=9)
hf.create_dataset('evt_type', data=evt_type, compression="gzip", compression_opts=9)
hf.create_dataset('rf_evt_type', data=rf_evt_type, compression="gzip", compression_opts=9)
hf.create_dataset('bias_rf_evt_type', data=bias_rf_evt_type, compression="gzip", compression_opts=9)
hf.create_dataset('bias_only_rf_evt_type', data=bias_only_rf_evt_type, compression="gzip", compression_opts=9)

hf.close()

print('file is in:',path+file_name)

# quick size check
size_checker(path+file_name)







