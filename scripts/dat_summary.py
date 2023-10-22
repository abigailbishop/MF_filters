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
count_i = int(sys.argv[2])
count_f = int(sys.argv[3])
count_ff = count_i + count_f

if Station == 2: num_configs = 7
if Station == 3: num_configs = 9

# sort
d_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/reco_ele/'
m_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/mf/'

r_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/Hist/'
file_name = f'Info_Summary_A{Station}.h5'
hf = h5py.File(r_path + file_name, 'r')
runs = hf['runs'][:]
configs = hf['configs'][:]
run_ep = hf['run_ep'][:]
evt_ep = hf['evt_ep'][:]
trig_ep = hf['trig_ep'][:]
con_ep = hf['con_ep'][:]
unix_ep = hf['unix_ep'][:]
num_evts = hf['num_evts'][:]
evt_i = int(np.nansum(num_evts[:count_i]))
evt_f = int(np.nansum(num_evts[:count_ff]))
print(evt_i, evt_f)
evt_len = int(evt_f - evt_i)
run_pa = run_ep[evt_i:evt_f]
runs_pa = runs[count_i:count_ff]
num_evts_pa = num_evts[count_i:count_ff]
num_runs = len(runs_pa)
print(evt_len, len(run_pa))
del hf, r_path, file_name

known_issue = known_issue_loader(Station, verbose = True)
bad_runs = known_issue.get_knwon_bad_run(use_qual = True)
b_runs = np.in1d(runs, bad_runs).astype(int)
del known_issue, bad_runs

ang_num = np.arange(2, dtype = int)
ang_len = len(ang_num)
pol_num = np.arange(2, dtype = int)
pol_len = len(pol_num)
sol_num = np.arange(3, dtype = int)
sol_len = len(sol_num)
rad = np.array([41, 170, 300, 450, 600], dtype = float)
rad_len = len(rad)
rad_num = np.arange(rad_len, dtype = int)
theta = 90 - np.linspace(0.5, 179.5, 179 + 1)
the_len = len(theta)

coef_max = np.full((pol_len, evt_len), np.nan, dtype = float) # pol, evt
coord_max = np.full((ang_len + 2, pol_len, evt_len), np.nan, dtype = float) # thepirz, pol, evt
coef_s_max = np.copy(coef_max)
coord_s_max = np.copy(coord_max)
mf_max = np.full((pol_len, evt_len), np.nan, dtype = float) # pols, evts
mf_temp = np.full((ang_len, pol_len, evt_len), np.nan, dtype = float) # thephi, pols, evts 

flat_len = the_len * rad_len * sol_len
theta_ex = np.full((the_len, rad_len, sol_len), np.nan, dtype = float)
theta_ex[:] = theta[:, np.newaxis, np.newaxis]
theta_flat = np.reshape(theta_ex, (flat_len))
rad_ex = np.full((the_len, rad_len, sol_len), np.nan, dtype = float)
rad_ex[:] = rad[np.newaxis, :, np.newaxis]
rad_flat = np.reshape(rad_ex, (flat_len))
z_flat = np.sin(np.radians(theta_flat)) * rad_flat
del theta_ex, rad_ex

sur_ang = np.array([180, 180, 37, 24, 17], dtype = float)
theta_map = np.full((pol_len, the_len, rad_len, sol_len), np.nan, dtype = float)
theta_map[:] = theta[np.newaxis, :, np.newaxis, np.newaxis]
sur_bool = theta_map <= sur_ang[np.newaxis, np.newaxis, :, np.newaxis]
sur_bool_flat = np.reshape(sur_bool, (pol_len, flat_len))
del sur_ang, theta_map, sur_bool

for r in tqdm(range(num_runs)):
    
  #if r <10:

    run_idx = np.in1d(run_pa, runs_pa[r])
    num_evts = num_evts_pa[r]
    evt_num = np.arange(num_evts, dtype = int)

    r_name = f'{d_path}reco_ele_A{Station}_R{runs_pa[r]}.h5'
    hf = h5py.File(r_name, 'r')
    coef_tot = hf['coef'][:] # pol, theta, rad, sol, evt
    coord_tot = hf['coord'][:] # pol, theta, rad, sol, evt
    coef_tot[np.isnan(coef_tot)] = -1
    del hf, r_name
    coef_re = np.reshape(coef_tot, (pol_len, flat_len, -1))
    coord_re = np.reshape(coord_tot, (pol_len, flat_len, -1))
    del coef_tot, coord_tot
    for t in range(2):
        if t == 1:
            coef_re[sur_bool_flat] = -1
            coord_re[sur_bool_flat] = np.nan
        coef_max_idx = np.nanargmax(coef_re, axis = 1)
        coef_max1 = coef_re[pol_num[:, np.newaxis], coef_max_idx, evt_num[np.newaxis, :]] # pol, evt
        neg_idx = coef_max1 < 0
        coef_max1[neg_idx] = np.nan
        coord_max1 = np.full((ang_len + 2, pol_len, num_evts), np.nan, dtype = float) # thepir, pol, evt
        coord_max1[0] = theta_flat[coef_max_idx]
        coord_max1[1] = coord_re[pol_num[:, np.newaxis], coef_max_idx, evt_num[np.newaxis, :]]
        coord_max1[2] = rad_flat[coef_max_idx]
        coord_max1[3] = z_flat[coef_max_idx]
        coord_max1[:, neg_idx] = np.nan
        del coef_max_idx, neg_idx
        if t == 0:
            coef_max[:, run_idx] = coef_max1
            coord_max[:, :, run_idx] = coord_max1
        else:
            coef_s_max[:, run_idx] = coef_max1
            coord_s_max[:, :, run_idx] = coord_max1
        del coef_max1, coord_max1
    del coef_re, coord_re, evt_num, num_evts

    m_name = f'{m_path}mf_A{Station}_R{runs_pa[r]}.h5'
    hf = h5py.File(m_name, 'r')
    mf_t_p = hf['mf_temp'][:, 1:3] # pol, thepi, evt
    mf_max[:, run_idx] = hf['mf_max'][:pol_len]
    mf_temp[:, :, run_idx] = np.transpose(mf_t_p, (1, 0, 2)) # thepi, pol, evt
    del m_name, hf, mf_t_p
    del run_idx
    
print(coef_max.shape)
print(coord_max.shape)
print(mf_max.shape)
print(mf_temp.shape)

path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/Hist/'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

file_name = f'Data_Summary_v16_A{Station}_R{count_i}.h5'
hf = h5py.File(file_name, 'w')
hf.create_dataset('runs', data=runs, compression="gzip", compression_opts=9)
hf.create_dataset('b_runs', data=b_runs, compression="gzip", compression_opts=9)
hf.create_dataset('configs', data=configs, compression="gzip", compression_opts=9)
hf.create_dataset('run_ep', data=run_ep, compression="gzip", compression_opts=9)
hf.create_dataset('evt_ep', data=evt_ep, compression="gzip", compression_opts=9)
hf.create_dataset('trig_ep', data=trig_ep, compression="gzip", compression_opts=9)
hf.create_dataset('con_ep', data=con_ep, compression="gzip", compression_opts=9)
hf.create_dataset('unix_ep', data=unix_ep, compression="gzip", compression_opts=9)
hf.create_dataset('coef_s_max', data=coef_s_max, compression="gzip", compression_opts=9)
hf.create_dataset('coord_s_max', data=coord_s_max, compression="gzip", compression_opts=9)
hf.create_dataset('coef_max', data=coef_max, compression="gzip", compression_opts=9)
hf.create_dataset('coord_max', data=coord_max, compression="gzip", compression_opts=9)
hf.create_dataset('mf_max', data=mf_max, compression="gzip", compression_opts=9)
hf.create_dataset('mf_temp', data=mf_temp, compression="gzip", compression_opts=9)
hf.close()
print('file is in:',path+file_name, size_checker(path+file_name))
print('done!')





