import numpy as np
import os, sys
from glob import glob
import h5py
from tqdm import tqdm

curr_path = os.getcwd()
sys.path.append(curr_path+'/../')
from tools.ara_run_manager import file_sorter
from tools.ara_utility import size_checker

Station = int(sys.argv[1])
Data = str(sys.argv[2])

def main(Station, Data):

    # sort
    d_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/rms_sim/*'
    d_list, d_run_tot, d_run_range, d_len = file_sorter(d_path)
    del d_run_range

    if Station == 2: num_configs = 7
    if Station == 3: num_configs = 9

    r_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/rms_sim_merge/'
    rms_mean = np.full((16, num_configs), np.nan, dtype = float)
    for c in range(num_configs):
        file_name = f'{r_path}rms_A{Station}_R{c + 1}.h5'
        hf = h5py.File(file_name, 'r')
        rms_mean[:, c] = hf['rms_mean'][:]
        print(rms_mean[:, c])
        print(file_name, size_checker(file_name))
        del file_name, hf

    output_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/snr_sim/'
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for r in tqdm(range(len(d_run_tot))):
       
      #if r < 10:
        snr_path = d_list[r].replace('rms', 'snr')
        if os.path.exists(snr_path):
            continue

        hf = h5py.File(d_list[r], 'r')
        config = hf['config'][:]
        configs = int(config[2] - 1)
        p2p = hf['p2p'][:]
        entry_num = hf['entry_num'][:]
        del hf

        snr = p2p / 2 / rms_mean[:, configs][:, np.newaxis]
        del configs

        #snr_path = d_list[r].replace('rms', 'snr')
        hf = h5py.File(snr_path, 'w')
        hf.create_dataset('config', data=config, compression="gzip", compression_opts=9)
        hf.create_dataset('entry_num', data=entry_num, compression="gzip", compression_opts=9)
        hf.create_dataset('snr', data=snr, compression="gzip", compression_opts=9)
        hf.close()
        #print(snr_path, size_checker(snr_path))
        del snr_path, entry_num, snr, p2p

    print('done!')


main(Station, Data)
