import numpy as np
import os, sys
import re
from glob import glob
import h5py
from tqdm import tqdm

curr_path = os.getcwd()
sys.path.append(curr_path+'/../')
from tools.ara_run_manager import file_sorter
from tools.ara_utility import size_checker
from tools.ara_known_issue import known_issue_loader

Station = int(sys.argv[1])

knwon_issue = known_issue_loader(Station)
bad_runs = knwon_issue.get_knwon_bad_run()

# sort
d_path = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{Station}/l1_full/*'
print(d_path)
d_list, d_run_tot, d_run_range = file_sorter(d_path)
del d_run_range

# config array
config_arr = []
config_arr_cut = []
run_arr = []
run_arr_cut = []

l1_hist = []
l1_cut_hist = []

l1_hist2d_max = []
l1_cut_hist2d_max = []

l1_range = np.arange(0,1000)
l1_bins = np.linspace(0,1000,1000+1)
l1_bin_center = (l1_bins[1:] + l1_bins[:-1]) / 2
min_range = np.arange(0, 360)
min_bins = np.linspace(0, 360, 360 + 1)
min_bin_center = (min_bins[1:] + min_bins[:-1]) / 2

l1_bin_len = len(l1_bin_center)
min_bin_len = len(min_bin_center)
l1_hist2d = np.full((min_bin_len, l1_bin_len, 16), 0, dtype = int)
l1_cut_hist2d = np.copy(l1_hist2d)

d_shape = l1_hist2d.shape

def get_2d_max(dat_2d, use_min = False):
    dat_max = np.copy(dat_2d)
    dat_max = dat_max.astype(float)
    dat_max[dat_max != 0] = 1
    dat_max *= l1_bin_center[np.newaxis, :, np.newaxis]
    if use_min:
        dat_max = np.nanmin(dat_max, axis = 1)
    else:
        dat_max = np.nanmax(dat_max, axis = 1)
    return dat_max

for r in tqdm(range(len(d_run_tot))):
    
  #if r <10:

    hf = h5py.File(d_list[r], 'r')
    config = hf['config'][2]
    config_arr.append(config)
    run_arr.append(d_run_tot[r])
    
    l1_unix = hf['unix_time'][:]
    min_unix = np.copy(l1_unix)
    min_unix = min_unix.astype(float)
    min_unix = (min_unix - min_unix[0]) / 60
    trig_ch = hf['trig_ch'][:]
    l1_rate = hf['l1_rate'][:] / 32
    l1_rate = l1_rate[:, trig_ch]
    del trig_ch

    l1_hist_r = np.full((l1_bin_len, 16), 0, dtype = int)
    l1_hist2d_r = np.full(d_shape, 0, dtype = int)
    for a in range(16):
        l1_hist_r[:, a] = np.histogram(l1_rate[:, a], bins = l1_bins)[0].astype(int)
        l1_hist2d_r[:, :, a] = np.histogram2d(min_unix, l1_rate[:, a], bins = (min_bins, l1_bins))[0].astype(int)
    l1_hist2d_max_r = get_2d_max(l1_hist2d_r)
    l1_hist2d_max.append(l1_hist2d_max_r)
    l1_hist.append(l1_hist_r)
    l1_hist2d += l1_hist2d_r

    if d_run_tot[r] in bad_runs:
        #print('bad run:', d_list[r], d_run_tot[r])
        continue

    config_arr_cut.append(config)
    run_arr_cut.append(d_run_tot[r])

    q_path = f'/data/user/mkim/OMF_filter/ARA0{Station}/qual_cut_full/qual_cut_full_A{Station}_R{d_run_tot[r]}.h5'
    hf_q = h5py.File(q_path, 'r')
    unix_time = hf_q['unix_time'][:]
    total_qual_cut = hf_q['total_qual_cut'][:]
    total_qual_cut[:, 17] = 0 #remove unlock unix time
    qual_cut_sum = np.nansum(total_qual_cut, axis = 1)
    l1_unix_bins = np.copy(l1_unix)
    l1_unix_bins = l1_unix_bins.astype(float)
    l1_unix_bins -= 0.5
    l1_unix_bins = np.append(l1_unix_bins, l1_unix_bins[-1]+1)
    unix_clean = np.histogram(unix_time, bins = l1_unix_bins, weights = qual_cut_sum)[0].astype(int)
    
    l1_rate_cut = np.copy(l1_rate)
    l1_rate_cut = l1_rate_cut.astype(float)
    l1_rate_cut[unix_clean != 0] = np.nan

    l1_cut_hist_r = np.full((l1_bin_len, 16), 0, dtype = int)
    l1_cut_hist2d_r = np.full(d_shape, 0, dtype = int)
    for a in range(16):
        l1_cut_hist_r[:, a] = np.histogram(l1_rate_cut[:, a], bins = l1_bins)[0].astype(int)
        l1_cut_hist2d_r[:, :, a] = np.histogram2d(min_unix, l1_rate_cut[:, a], bins = (min_bins, l1_bins))[0].astype(int)
    l1_cut_hist2d_max_r = get_2d_max(l1_cut_hist2d_r)
    l1_cut_hist2d_max.append(l1_cut_hist2d_max_r)
    l1_cut_hist.append(l1_cut_hist_r)
    l1_cut_hist2d += l1_cut_hist2d_r

    del hf, l1_rate, min_unix, l1_unix, l1_unix_bins, unix_time, l1_rate_cut, unix_clean, total_qual_cut, hf_q, q_path

path = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{Station}/Hist/'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

file_name = f'L1_Rate_A{Station}.h5'
hf = h5py.File(file_name, 'w')
hf.create_dataset('config_arr', data=np.asarray(config_arr), compression="gzip", compression_opts=9)
hf.create_dataset('config_arr_cut', data=np.asarray(config_arr_cut), compression="gzip", compression_opts=9)
hf.create_dataset('run_arr', data=np.asarray(run_arr), compression="gzip", compression_opts=9)
hf.create_dataset('run_arr_cut', data=np.asarray(run_arr_cut), compression="gzip", compression_opts=9)
hf.create_dataset('l1_range', data=l1_range, compression="gzip", compression_opts=9)
hf.create_dataset('l1_bins', data=l1_bins, compression="gzip", compression_opts=9)
hf.create_dataset('l1_bin_center', data=l1_bin_center, compression="gzip", compression_opts=9)
hf.create_dataset('min_range', data=min_range, compression="gzip", compression_opts=9)
hf.create_dataset('min_bins', data=min_bins, compression="gzip", compression_opts=9)
hf.create_dataset('min_bin_center', data=min_bin_center, compression="gzip", compression_opts=9)
hf.create_dataset('l1_hist', data=np.asarray(l1_hist), compression="gzip", compression_opts=9)
hf.create_dataset('l1_cut_hist', data=np.asarray(l1_cut_hist), compression="gzip", compression_opts=9)
hf.create_dataset('l1_hist2d', data=l1_hist2d, compression="gzip", compression_opts=9)
hf.create_dataset('l1_cut_hist2d', data=l1_cut_hist2d, compression="gzip", compression_opts=9)
hf.create_dataset('l1_hist2d_max', data=np.asarray(l1_hist2d_max), compression="gzip", compression_opts=9)
hf.create_dataset('l1_cut_hist2d_max', data=np.asarray(l1_cut_hist2d_max), compression="gzip", compression_opts=9)
hf.close()
print('file is in:',path+file_name)
# quick size check
size_checker(path+file_name)






