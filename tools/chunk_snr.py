import numpy as np
from tqdm import tqdm
import h5py

def snr_collector(Data, Ped, analyze_blind_dat = False, use_l2 = False, no_tqdm = False):

    print('Collecting snr starts!')
    if use_l2:
        from tools.ara_data_load import ara_l2_loader
    else:
        from tools.ara_data_load import ara_uproot_loader
        from tools.ara_data_load import ara_root_loader
    from tools.ara_constant import ara_const
    from tools.ara_wf_analyzer import wf_analyzer
    from tools.ara_quality_cut import get_bad_events

    # geom. info.
    ara_const = ara_const()
    num_ants = ara_const.USEFUL_CHAN_PER_STATION
    del ara_const

    # data config
    if use_l2:
        ara_root = ara_l2_loader(Data)
        ara_root.get_sub_info()
        num_evts = ara_root.num_evts
        evt_num = ara_root.evt_num
        trig_type = ara_root.trig_type
        st = ara_root.station_id
        run = ara_root.run
    else:
        ara_uproot = ara_uproot_loader(Data)
        ara_uproot.get_sub_info()
        evt_num = ara_uproot.evt_num
        num_evts = ara_uproot.num_evts
        trig_type = ara_uproot.get_trig_type()
        st = ara_uproot.station_id
        run = ara_uproot.run
        ara_root = ara_root_loader(Data, Ped, st, ara_uproot.year)
        del ara_uproot

    # pre quality cut
    daq_qual_cut_sum, tot_qual_cut_sum = get_bad_events(st, run, analyze_blind_dat = analyze_blind_dat, verbose = True, evt_num = evt_num, qual_type = 2)

    # wf analyzer
    wf_int = wf_analyzer(use_time_pad = True, use_band_pass = True, use_cw = True, use_l2 = use_l2, analyze_blind_dat = analyze_blind_dat, st = st, run = run)
    del st, run

    # output array  
    rms = np.full((num_ants, num_evts), np.nan, dtype = float)
    p2p = np.copy(rms)

    # loop over the events
    for evt in tqdm(range(num_evts), disable = no_tqdm):
      #if evt <100:        

        if daq_qual_cut_sum[evt]:
            continue
 
        # get entry and wf
        ara_root.get_entry(evt)
        ara_root.get_useful_evt(ara_root.cal_type.kLatestCalibWithOutTrimFirstBlock)


        
        # loop over the antennas
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_band_pass = True, use_cw = True, use_p2p = True, evt = evt,  dedisperse=True)
            p2p[ant, evt] = wf_int.int_p2p
            
            # split waveform in quarters, save rms as average rms of two lowest rms quarters
            segmented_wf = np.array_split(raw_v[raw_v !=0], 4) 
            rms_values = [np.sqrt(np.mean(segment**2)) for segment in segmented_wf]
            sorted_rms = sorted(rms_values)
            mean_of_two_lowest_rms = np.mean(sorted_rms[:2])
            rms[ant, evt] = mean_of_two_lowest_rms

            del raw_t, raw_v, sorted_rms, mean_of_two_lowest_rms, segmented_wf, rms_values
            ara_root.del_TGraph()
        ara_root.del_usefulEvt()   

        #rms[:, evt] = np.nanstd(wf_int.pad_v, axis = 0)
    del ara_root, num_ants, wf_int

    rms_qual = np.full((num_evts, 4), False, dtype = bool)
    rms_qual[:, 0] = np.logical_and(tot_qual_cut_sum == 0, trig_type == 2)
    rms_qual[:, 1] = np.logical_and(tot_qual_cut_sum == 0, trig_type == 0)
    rms_qual[:, 2] = np.logical_and(daq_qual_cut_sum == 0, trig_type == 2)
    rms_qual[:, 3] = np.logical_and(daq_qual_cut_sum == 0, trig_type == 0)

    clean_idx = np.full((num_evts), False, dtype = bool)
    for q in range(rms_qual.shape[1]):
        if np.count_nonzero(rms_qual[:, q]) > 0:
            clean_idx = rms_qual[:, q]
            break
    rms_qual = np.count_nonzero(rms_qual, axis = 0)
    print(f'rms quality: {rms_qual}')
    rms_copy = np.copy(rms)
    rms_copy[:, ~clean_idx] = np.nan
    rms_mean = np.nanmean(rms_copy, axis = 1)
    snr = p2p / 2 / rms
    
    del num_evts, rms_copy, daq_qual_cut_sum, tot_qual_cut_sum

    print('SNR collecting is done!')

    return {'evt_num':evt_num,
            'trig_type':trig_type,
            'clean_idx':clean_idx,
            'snr':snr,
            'p2p':p2p,
            'rms':rms,
            'rms_mean':rms_mean,
            'rms_qual':rms_qual}








