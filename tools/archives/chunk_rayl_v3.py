import numpy as np
from tqdm import tqdm

def rayl_collector(Data, Ped, analyze_blind_dat = False):

    print('Collecting rayl. starts!')

    from tools.ara_data_load import ara_uproot_loader
    from tools.ara_data_load import ara_root_loader
    from tools.ara_constant import ara_const
    from tools.ara_quality_cut import qual_cut_loader
    from tools.ara_quality_cut import cw_qual_cut_loader
    from tools.ara_wf_analyzer import wf_analyzer
    from tools.ara_wf_analyzer import get_rayl_distribution

    # geom. info.
    ara_const = ara_const()
    num_ants = ara_const.USEFUL_CHAN_PER_STATION 
    del ara_const

    # data config
    ara_uproot = ara_uproot_loader(Data)
    ara_uproot.get_sub_info()
    evt_num = ara_uproot.evt_num
    unix_time = ara_uproot.unix_time
    pps_number = ara_uproot.pps_number
    trig_type = ara_uproot.get_trig_type()
    st = ara_uproot.station_id
    run = ara_uproot.run
    ara_root = ara_root_loader(Data, Ped, st, ara_uproot.year)
    del ara_uproot

    # qulity cut
    ara_qual = qual_cut_loader(analyze_blind_dat = analyze_blind_dat, verbose = True)
    total_qual_cut = ara_qual.load_qual_cut_result(st, run)
    clean_rf_entry = ara_qual.get_useful_events(use_entry = True, use_qual = True, trig_idx = 0)
    clean_rf_evt = ara_qual.get_useful_events(use_qual = True, use_rp_ants = True, trig_idx = 0)
    num_clean_rf_evts = np.copy(ara_qual.num_useful_evts)
    clean_rf_rp_ants = np.copy(ara_qual.clean_rp_ants)
    clean_soft_entry = ara_qual.get_useful_events(use_entry = True, use_qual = True, trig_idx = 2)
    clean_soft_evt = ara_qual.get_useful_events(use_qual = True, use_rp_ants = True, trig_idx = 2)
    num_clean_soft_evts = np.copy(ara_qual.num_useful_evts)
    clean_soft_rp_ants = np.copy(ara_qual.clean_rp_ants)

    #clean_rf_rp_ant_sum = np.nansum(clean_rf_rp_ants, axis = 0) == 0
    #clean_soft_rp_ant_sum = np.nansum(clean_soft_rp_ants, axis = 0) == 0
    

    print(f'Number of clean rf event is {num_clean_rf_evts}') 
    print(f'Number of clean soft event is {num_clean_soft_evts}') 
    del ara_qual

    # cw quality cut
    cw_qual = cw_qual_cut_loader(st, run, evt_num, pps_number, verbose = True)
    cw_params = cw_qual.ratio_cut
    del cw_qual, st, run

    # wf analyzer
    wf_int = wf_analyzer(use_time_pad = True, use_freq_pad = True, use_rfft = True, use_band_pass = True, use_cw = True, cw_params = cw_params)
    fft_len = wf_int.pad_fft_len
    freq_range = wf_int.pad_zero_freq 
    del cw_params

    # output 
    clean_rf_rffts = np.full((fft_len, num_ants, num_clean_rf_evts), np.nan, dtype = float)
    clean_soft_rffts = np.full((fft_len, num_ants, num_clean_soft_evts), np.nan, dtype = float)
    print(f'fft array dim.: {clean_rf_rffts.shape}, {clean_soft_rffts.shape}')
    print(f'fft array size: ~{np.round(clean_rf_rffts.nbytes/1024/1024)} MB, ~{np.round(clean_soft_rffts.nbytes/1024/1024)} MB')
    del fft_len

    # loop over the events
    c_evts = [num_clean_rf_evts, num_clean_soft_evts]
    c_entry = [clean_rf_entry, clean_soft_entry]
    r_ant = [clean_rf_rp_ants, clean_soft_rp_ants]
    for trig in range(2): # uhhhhh.....
        num_clean_evts = c_evts[trig]
        clean_entry = c_entry[trig]
        rp_ant = r_ant[trig].astype(bool)
        for evt in tqdm(range(num_clean_evts)):
           #if evt <100:        

            # get entry and wf
            ara_root.get_entry(clean_entry[evt])
            ara_root.get_useful_evt(ara_root.cal_type.kLatestCalib)
        
            # loop over the antennas
            for ant in range(num_ants):
                raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
                wf_int.get_int_wf(raw_t, raw_v, ant, use_zero_pad = True, use_band_pass = True, use_cw = rp_ant[ant, evt])
                del raw_t, raw_v 
                ara_root.del_TGraph()
            ara_root.del_usefulEvt()

            wf_int.get_fft_wf(use_zero_pad = True, use_rfft = True, use_abs = True)
            if trig == 0:
                clean_rf_rffts[:, :, evt] = wf_int.pad_fft
            else:
                clean_soft_rffts[:, :, evt] = wf_int.pad_fft
        del num_clean_evts, clean_entry, rp_ant  
    del num_ants, ara_root, wf_int, c_evts, c_entry, r_ant, clean_rf_entry, clean_soft_entry, num_clean_rf_evts, num_clean_soft_evts
   
    # rayl fit 
    binning = np.array([1000], dtype = int)
    rf_rayl, clean_rf_rfft_2d, clean_rf_bin_edges = get_rayl_distribution(clean_rf_rffts, binning = binning[0])
    soft_rayl, clean_soft_rfft_2d, clean_soft_bin_edges = get_rayl_distribution(clean_soft_rffts, binning = binning[0])
    del clean_rf_rffts, clean_soft_rffts

    print('Rayl. collecting is done!')

    return {'evt_num':evt_num,
            'clean_rf_evt':clean_rf_evt,
            'clean_soft_evt':clean_soft_evt,
            #'clean_rf_rp_ants':clean_rf_rp_ants,
            #'clean_soft_rp_ants':clean_soft_rp_ants,
            'trig_type':trig_type,
            'unix_time':unix_time,
            'pps_number':pps_number,
            'total_qual_cut':total_qual_cut,
            'freq_range':freq_range,
            'clean_rf_bin_edges':clean_rf_bin_edges,
            'clean_soft_bin_edges':clean_soft_bin_edges,
            'clean_rf_rfft_2d':clean_rf_rfft_2d,
            'clean_soft_rfft_2d':clean_soft_rfft_2d,
            'binning':binning,
            'rf_rayl':rf_rayl,
            'soft_rayl':soft_rayl}



