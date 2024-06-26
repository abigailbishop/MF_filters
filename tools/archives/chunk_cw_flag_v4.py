import os
import h5py
import numpy as np
from tqdm import tqdm

def cw_flag_collector(Data, Ped, analyze_blind_dat = False, use_l2 = False, no_tqdm = False):

    print('Collecting cw flag starts!')

    from tools.ara_data_load import ara_uproot_loader
    from tools.ara_data_load import ara_root_loader
    from tools.ara_constant import ara_const
    from tools.ara_wf_analyzer import wf_analyzer
    from tools.ara_cw_filters import py_phase_variance
    from tools.ara_cw_filters import py_testbed
    from tools.ara_cw_filters import group_bad_frequency
    from tools.ara_quality_cut import get_bad_events
    from tools.ara_known_issue import known_issue_loader
    from tools.ara_utility import size_checker

    # geom. info.
    ara_const = ara_const()
    num_ants = ara_const.USEFUL_CHAN_PER_STATION
    del ara_const

    # data config
    ara_uproot = ara_uproot_loader(Data)
    ara_uproot.get_sub_info()
    trig_type = ara_uproot.get_trig_type()
    evt_num = ara_uproot.evt_num
    entry_num = ara_uproot.entry_num
    num_evts = ara_uproot.num_evts
    st = ara_uproot.station_id
    yr = ara_uproot.year
    run = ara_uproot.run
    ara_root = ara_root_loader(Data, Ped, st, yr)
    del ara_uproot

    # pre quality cut
    daq_qual_cut_sum = get_bad_events(st, run, analyze_blind_dat = analyze_blind_dat, verbose = True, evt_num = evt_num, use_1st = True)[0]

    known_issue = known_issue_loader(st)
    bad_ant = known_issue.get_bad_antenna(run, print_integer = True)
    del known_issue

    # wf analyzer
    wf_int = wf_analyzer(use_time_pad = True, use_freq_pad = True, use_rfft = True)
    freq_range = wf_int.pad_zero_freq

    # cw class
    cw_testbed = py_testbed(st, run, freq_range, analyze_blind_dat = analyze_blind_dat, verbose = True)
    testbed_params = np.array([cw_testbed.dB_cut, cw_testbed.dB_cut_broad, cw_testbed.num_coinc, cw_testbed.freq_range_broad, cw_testbed.freq_range_near])
    cw_phase = py_phase_variance(st, run, freq_range)
    evt_len = cw_phase.evt_len
    phase_params = np.array([cw_phase.sigma_thres, evt_len])

    # output array  
    sigma = []
    phase_idx = []
    testbed_idx = []
    empty = np.full((0), 0, dtype = int)
    empty_float = np.full((0), np.nan, dtype = float)
    clean_entry = entry_num[np.logical_and(~daq_qual_cut_sum, trig_type != 1)]
    del entry_num

    # loop over the events
    evt = 0
    evt_backup = 0
    evt_counts = 0
    pbar = tqdm(total = num_evts, disable = no_tqdm)
    while evt < num_evts:
  
        if evt == evt_backup:
            pbar.update(1)
 
        if daq_qual_cut_sum[evt]:
            if evt == evt_backup:
                sigma.append(empty_float)
                phase_idx.append(empty)
                testbed_idx.append(empty)
            evt_backup += 1
            evt = evt_backup
            continue

        # get entry and wf
        ara_root.get_entry(evt)
        ara_root.get_useful_evt(ara_root.cal_type.kLatestCalibWithOutTrimFirstBlock)
    
        # loop over the antennas
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_zero_pad = True)
            del raw_t, raw_v
            ara_root.del_TGraph()
        ara_root.del_usefulEvt()   

        wf_int.get_fft_wf(use_zero_pad = True, use_rfft = True, use_phase = True, use_abs = True, use_norm = True, use_dBmHz = True)
        rfft_phase = wf_int.pad_phase

        cw_phase.get_phase_differences(rfft_phase, evt_counts % evt_len, trig_type[evt])
        cw_phase.get_bad_phase()
        sigmas = cw_phase.bad_sigma 
        phase_idxs = cw_phase.bad_idx
        if evt != evt_backup:
            sigma[evt] = np.concatenate((sigma[evt], sigmas))
            phase_idx[evt] = np.concatenate((phase_idx[evt], phase_idxs)) 
        else:
            rfft_dbmhz = wf_int.pad_fft
            cw_testbed.get_bad_magnitude(rfft_dbmhz, trig_type[evt])
            testbed_idxs = cw_testbed.bad_idx
            testbed_idx.append(testbed_idxs)
            sigma.append(sigmas)
            phase_idx.append(phase_idxs)
            del rfft_dbmhz
        del rfft_phase
        
        if trig_type[evt] == 1:
            evt_backup += 1
            evt = evt_backup
            continue
        
        time_travel_idx = evt_counts - evt_len + 1
        if time_travel_idx >= 0:
            time_travel_entry = clean_entry[time_travel_idx]
            sigma[time_travel_entry] = np.concatenate((sigma[time_travel_entry], sigmas))
            phase_idx[time_travel_entry] = np.concatenate((phase_idx[time_travel_entry], phase_idxs))
        else:
            time_travel_entry = 0
        evt_counts += 1
          
        time_travel_cal_entry = time_travel_entry - 1
        if time_travel_idx >= 0 and time_travel_cal_entry >= 0 and trig_type[time_travel_cal_entry] == 1:  
            evt = time_travel_cal_entry
        else:
            evt_backup += 1
            evt = evt_backup
        del time_travel_idx, time_travel_entry#, time_travel_cal_entry
    del ara_root, num_ants, wf_int, cw_phase, cw_testbed, evt_len, clean_entry
    pbar.close()

    # to numpy array
    sigma = np.asarray(sigma)
    phase_idx = np.asarray(phase_idx)
    testbed_idx = np.asarray(testbed_idx)

    # group bad frequency
    cw_freq = group_bad_frequency(st, run, freq_range, verbose = True) # constructor for bad frequency grouping function
     
    # output array
    bad_range = []

    # loop over the events
    for evt in tqdm(range(num_evts), disable = no_tqdm):
      #if evt == 0:

        # quality cut
        if daq_qual_cut_sum[evt]:
            bad_range.append(empty_float)
            continue

        bad_range_evt = cw_freq.get_pick_freqs_n_bands(sigma[evt], phase_idx[evt], testbed_idx[evt]).flatten()
        bad_range.append(bad_range_evt)
    del num_evts, daq_qual_cut_sum, cw_freq

    # to numpy array
    bad_range = np.asarray(bad_range)
 
    blind_type = ''
    if analyze_blind_dat:
        blind_type = '_full'
    output_path = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{st}/cw_band{blind_type}/'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    h5_file_name = f'cw_band{blind_type}_A{st}_R{run}.h5'
    hf = h5py.File(f'{output_path}{h5_file_name}', 'w')
    hf.create_dataset('evt_num', data=evt_num, compression="gzip", compression_opts=9)
    try:
        hf.create_dataset('bad_range', data=bad_range, compression="gzip", compression_opts=9)
    except TypeError:
        dt = h5py.vlen_dtype(np.dtype(float))
        hf.create_dataset('bad_range', data=bad_range, dtype = dt, compression="gzip", compression_opts=9)    
    hf.close()
    print(f'output is {output_path}{h5_file_name}.', size_checker(f'{output_path}{h5_file_name}'))
    del st, run, blind_type, output_path, h5_file_name

    print('CW flag collecting is done!')

    return {'evt_num':evt_num,
            'bad_ant':bad_ant,
            'freq_range':freq_range,
            'sigma':sigma,
            'phase_idx':phase_idx,
            'testbed_idx':testbed_idx,
            'bad_range':bad_range,
            'testbed_params':testbed_params,
            'phase_params':phase_params}









