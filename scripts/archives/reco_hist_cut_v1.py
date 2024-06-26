import numpy as np
import os, sys
from glob import glob
import h5py
from tqdm import tqdm

curr_path = os.getcwd()
sys.path.append(curr_path+'/../')
from tools.ara_run_manager import file_sorter
from tools.ara_utility import size_checker
from tools.ara_run_manager import run_info_loader
from tools.ara_known_issue import known_issue_loader

Station = int(sys.argv[1])
count_i = int(sys.argv[2])
count_f = int(sys.argv[3])
count_ff = count_i + count_f

if Station == 2: num_configs = 7
if Station == 3: num_configs = 9

known_issue = known_issue_loader(Station, verbose = True)
bad_runs = known_issue.get_knwon_bad_run(use_qual = True)
del known_issue

# sort
d_path = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{Station}/reco/*'
d_list, d_run_tot, d_run_range = file_sorter(d_path)
q_path = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{Station}/qual_cut_full/'
b_path = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{Station}/baseline/'
del d_run_range

z_bins = np.linspace(0, 180, 180 + 1)
z_bin_center = (z_bins[1:] + z_bins[:-1]) / 2
z_bin_len = len(z_bin_center)
a_bins = np.linspace(0, 360, 360 + 1)
a_bin_center = (a_bins[1:] + a_bins[:-1]) / 2
a_bin_len = len(a_bin_center)

map_tot = np.full((a_bin_len, z_bin_len, 3, 2, 2, 2, num_configs), 0, dtype = int) # a, z, trig, pol, rad, sol, config
map_cut = np.copy(map_tot)
del z_bin_len, a_bin_len

for r in tqdm(range(len(d_run_tot))):
    
  #if r <10:
  if r >= count_i and r < count_ff:

    try:
        hf = h5py.File(d_list[r], 'r')
    except OSError: 
        print(d_list[r])
        continue

    b_name = f'{b_path}baseline_A{Station}_R{d_run_tot[r]}.h5'
    hf_b = h5py.File(b_name, 'r')
    trig = hf_b['trig_type'][:]
    del b_name, hf_b

    ara_run = run_info_loader(Station, d_run_tot[r])
    g_idx = ara_run.get_config_number() - 1
    del ara_run

    bad_run = d_run_tot[r] in bad_runs

    rf_t = trig == 0
    cal_t = trig == 1
    soft_t = trig == 2
    t_list = [rf_t, cal_t, soft_t]
    coord = hf['coord'][:] # pol, thephi, rad, sol, evt
    evt = hf['evt_num'][:]
    del hf, trig

    q_name = f'{q_path}qual_cut_full_A{Station}_R{d_run_tot[r]}.h5'
    hf_q = h5py.File(q_name, 'r')
    evt_full = hf_q['evt_num'][:]
    qual = hf_q['tot_qual_cut_sum'][:] != 0
    cut = np.in1d(evt, evt_full[qual])
    del q_name, hf_q, qual, evt_full, evt

    coord_cut = np.copy(coord)
    coord_cut[:, :, :, :, cut] = np.nan
    del cut

    for t in range(3):
        for pol in range(2):
            for rad in range(2):
                for sol in range(2):       
                    map_tot[:, :, t, pol, rad, sol, g_idx] += np.histogram2d(coord[pol, 1, rad, sol][t_list[t]], coord[pol, 0, rad, sol][t_list[t]], bins = (a_bins, z_bins))[0].astype(int)
                    if bad_run: continue
                    map_cut[:, :, t, pol, rad, sol, g_idx] += np.histogram2d(coord_cut[pol, 1, rad, sol][t_list[t]], coord_cut[pol, 0, rad, sol][t_list[t]], bins = (a_bins, z_bins))[0].astype(int)
    del g_idx, bad_run, coord, coord_cut, t_list, rf_t, cal_t, soft_t

path = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{Station}/Hist/'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

file_name = f'Reco_Map_New_A{Station}_R{count_i}.h5'
hf = h5py.File(file_name, 'w')
hf.create_dataset('a_bins', data=a_bins, compression="gzip", compression_opts=9)
hf.create_dataset('a_bin_center', data=a_bin_center, compression="gzip", compression_opts=9)
hf.create_dataset('z_bins', data=z_bins, compression="gzip", compression_opts=9)
hf.create_dataset('z_bin_center', data=z_bin_center, compression="gzip", compression_opts=9)
hf.create_dataset('map_tot', data=map_tot, compression="gzip", compression_opts=9)
hf.create_dataset('map_cut', data=map_cut, compression="gzip", compression_opts=9)
hf.close()
print('file is in:',path+file_name)
# quick size check
size_checker(path+file_name)






