import numpy as np
from tqdm import tqdm
import h5py

def reco_mf_collector(Data, Ped, analyze_blind_dat = False):

    print('Collecting reco mf starts!')

    from tools.ara_data_load import ara_uproot_loader
    from tools.ara_data_load import ara_root_loader
    from tools.ara_constant import ara_const
    from tools.ara_wf_analyzer import wf_analyzer
    from tools.ara_py_interferometers import py_interferometers
    from tools.ara_run_manager import run_info_loader

    # geom. info.
    ara_const = ara_const()
    num_ants = ara_const.USEFUL_CHAN_PER_STATION
    del ara_const

    # data config
    ara_uproot = ara_uproot_loader(Data)
    st = ara_uproot.station_id
    yr = ara_uproot.year
    run = ara_uproot.run
    num_evts = ara_uproot.num_evts
    ara_root = ara_root_loader(Data, Ped, st, yr)
    del ara_uproot

    # pre quality cut
    run_info = run_info_loader(st, run, analyze_blind_dat = analyze_blind_dat)
    daq_dat = run_info.get_result_path(file_type = 'qual_cut', verbose = True)
    daq_hf = h5py.File(daq_dat, 'r')
    daq_qual_cut_sum = daq_hf['daq_qual_cut_sum'][:]
    del daq_dat, daq_hf

    # snr info
    snr_dat = run_info.get_result_path(file_type = 'mf', verbose = True)
    snr_hf = h5py.File(snr_dat, 'r')
    snr = snr_hf['evt_wise_ant'][:]
    del run_info, snr_dat, snr_hf

    # wf analyzer
    wf_int = wf_analyzer(use_time_pad = True, use_band_pass = True, add_double_pad = True)

    # interferometers
    ara_int = py_interferometers(wf_int.pad_len, wf_int.dt, st, yr, run)
    pairs = ara_int.pairs
    v_pairs_len = ara_int.v_pairs_len
    snr_weights = snr[pairs[:, 0]] * snr[pairs[:, 1]]
    snr_v_sum = np.nansum(snr_weights[:v_pairs_len], axis = 0)
    snr_h_sum = np.nansum(snr_weights[v_pairs_len:], axis = 0)
    snr_weights[:v_pairs_len] /= snr_v_sum[np.newaxis, :]
    snr_weights[v_pairs_len:] /= snr_h_sum[np.newaxis, :] 
    del st, yr, run, snr, snr_v_sum, snr_h_sum, v_pairs_len, pairs

    # output array  
    coef = np.full((2, 2, num_evts), np.nan, dtype = float) # pol, rad
    coord = np.full((2, 2, 2, num_evts), np.nan, dtype = float) # thephi, pol, rad

    # loop over the events
    for evt in tqdm(range(num_evts)):
      #if evt <100:        
  
        print(snr_weights[:, evt])
        """
        if daq_qual_cut_sum[evt]:
            continue

        # get entry and wf
        ara_root.get_entry(evt)
        ara_root.get_useful_evt(ara_root.cal_type.kLatestCalib)
        
        # loop over the antennas
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_zero_pad = True, use_band_pass = True)
            del raw_t, raw_v
            ara_root.del_TGraph()
        ara_root.del_usefulEvt()   

        print(snr_weights[:, evt]) 
        coef[:, :, evt], coord[:, :, :, evt] = ara_int.get_sky_map(wf_int.pad_v, weights = snr_weights[:, evt])
        """
    del ara_root, num_evts, num_ants, wf_int, ara_int, daq_qual_cut_sum, snr_weights

    print('Reco mf collecting is done!')

    return {'coef':coef,
            'coord':coord}










