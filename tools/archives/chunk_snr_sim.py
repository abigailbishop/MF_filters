import os
import numpy as np
from tqdm import tqdm
import h5py

def snr_sim_collector(Data, Station, Year):

    print('Collecting sim snr starts!')

    from tools.ara_sim_load import ara_root_loader
    from tools.ara_constant import ara_const
    from tools.ara_wf_analyzer import wf_analyzer

    # const. info.
    ara_const = ara_const()
    num_ants = ara_const.USEFUL_CHAN_PER_STATION
    del ara_const

    # data config
    ara_root = ara_root_loader(Data, Station, Year)
    ara_root.get_sub_info(Data, get_angle_info = False)
    num_evts = ara_root.num_evts
    entry_num = ara_root.entry_num
    dt = ara_root.time_step
    wf_len = ara_root.waveform_length
    wf_time = ara_root.wf_time    
    #pnu = ara_root.pnu
    #inu_thrown = ara_root.inu_thrown
    #weight = ara_root.weight
    #probability = ara_root.probability
    #nuflavorint = ara_root.nuflavorint
    #nu_nubar = ara_root.nu_nubar
    #currentint = ara_root.currentint
    #elast_y = ara_root.elast_y
    #posnu = ara_root.posnu
    #nnu = ara_root.nnu

    # wf analyzer
    wf_int = wf_analyzer(dt = dt)
 
    # output array
    rms = np.full((num_ants, num_evts), np.nan, dtype = float)
    p2p = np.copy(rms)
 
    # loop over the events
    for evt in tqdm(range(num_evts)):
      #if evt <100: # debug 

        ara_root.get_entry(evt)
        for ant in range(num_ants):
            wf_v = ara_root.get_rf_ch_wf(ant)[1]
            rms[ant, evt] = np.nanstd(wf_v)
            p2p[ant, evt] = wf_int.get_p2p(wf_v, use_max = True)
            del wf_v
            ara_root.del_TGraph()
    del ara_root, num_ants, num_evts

    rms_mean = np.nanmean(rms, axis = 1)

    signal_key = 'signal_F'
    if Data.find(signal_key) != -1:
        r_idx = Data.find('_R')
        e_idx = Data.find('.txt', r_idx + 2)
        run = int(Data[r_idx + 2:e_idx])
        n_path =  os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{Station}/snr_sim/snr_AraOut.noise_A{Station}_R{run}.txt.run0.h5'
        print('noise_snr_path:', n_path)
        n_hf = h5py.File(n_path, 'r')
        noise_rms_mean = n_hf['rms_mean'][:]
        snr = p2p / 2 / noise_rms_mean[:, np.newaxis]
        del r_idx, e_idx, run, n_path, n_hf, noise_rms_mean
    else:
        snr = p2p / 2 / rms_mean[:, np.newaxis]

    print('Sim snr collecting is done!')

    return {'entry_num':entry_num,
            'dt':dt,
            'wf_time':wf_time,
            #'pnu':pnu,
            #'inu_thrown':inu_thrown,
            #'weight':weight,
            #'probability':probability,
            #'nuflavorint':nuflavorint,
            #'nu_nubar':nu_nubar,
            #'currentint':currentint,
            #'elast_y':elast_y,
            #'posnu':posnu,
            #'nnu':nnu,
            'snr':snr,
            'p2p':p2p,
            'rms':rms,
            'rms_mean':rms_mean} 

