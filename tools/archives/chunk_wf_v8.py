import numpy as np
from tqdm import tqdm
import h5py
import sys
from scipy.signal import medfilt

def wf_collector(Data, Ped, analyze_blind_dat = False, sel_evts = None):

    print('Collecting wf starts!')

    from tools.ara_run_manager import run_info_loader
    from tools.ara_data_load import ara_geom_loader
    from tools.ara_data_load import ara_uproot_loader
    from tools.ara_data_load import ara_root_loader
    from tools.ara_data_load import analog_buffer_info_loader
    from tools.ara_constant import ara_const
    from tools.ara_wf_analyzer import wf_analyzer
    from tools.ara_py_interferometers import py_interferometers
    from tools.ara_py_interferometers import get_products
    from tools.ara_known_issue import known_issue_loader
    from tools.ara_cw_filters import py_phase_variance
    from tools.ara_cw_filters import py_testbed

    # const. info.
    ara_const = ara_const()
    num_ants = ara_const.USEFUL_CHAN_PER_STATION
    num_samps = ara_const.SAMPLES_PER_BLOCK
    num_eles = ara_const.CHANNELS_PER_ATRI
    num_pols = ara_const.POLARIZATION
    del ara_const

    # data config
    ara_uproot = ara_uproot_loader(Data)
    ara_uproot.get_sub_info()
    evt_num = ara_uproot.evt_num
    entry_num = ara_uproot.entry_num
    unix_time = ara_uproot.unix_time
    pps_number = ara_uproot.pps_number
    trig_type = ara_uproot.get_trig_type()
    num_evts = ara_uproot.num_evts
    st = ara_uproot.station_id
    run = ara_uproot.run
    year = ara_uproot.year
    run_info = np.array([st, run, year], dtype = int)
    print('run info:', run_info)
    buffer_info = analog_buffer_info_loader(st, run, year, incl_cable_delay = True)
    buffer_info.get_int_time_info()
    int_num_samples = buffer_info.int_num_samples
    ara_root = ara_root_loader(Data, Ped, st, year)

    # channel mapping 
    ara_geom = ara_geom_loader(st, year, verbose = True)
    rf_ch = np.arange(num_ants, dtype = int)
    ele_ch = ara_geom.get_ele_ch_idx()
    trig_ch = ara_geom.get_trig_ch_idx()
    del ara_geom

    # bad antenna
    known_issue = known_issue_loader(st)
    bad_ant = known_issue.get_bad_antenna(run, print_integer = True)
    del known_issue

    # quality cut
    run_info1 = run_info_loader(st, run, analyze_blind_dat = analyze_blind_dat)
    qual_dat = run_info1.get_result_path(file_type = 'qual_cut', verbose = True, force_blind = True)
    daq_hf = h5py.File(qual_dat, 'r')
    evt_full = daq_hf['evt_num'][:]
    tot_qual_cut = daq_hf['tot_qual_cut'][:]
    daq_qual_cut = daq_hf['daq_qual_cut_sum'][:] != 0
    daq_qual_cut_sum = np.in1d(evt_num, evt_full[daq_qual_cut]).astype(int)
    del qual_dat, daq_hf

    # baseline
    base_dat =  run_info1.get_result_path(file_type = 'baseline', verbose = True)
    base_hf = h5py.File(base_dat, 'r')
    baseline_orig0 = base_hf['baseline'][:]
    baseline_orig = np.full(baseline_orig0.shape, np.nan, dtype = float)
    for a in range(16):
        baseline_orig[:, a] = medfilt(baseline_orig0[:, a], kernel_size = 39)
    del base_dat, base_hf
    print(baseline_orig.shape)

    # weight
    wei_key = 'snr'
    wei_dat = run_info1.get_result_path(file_type = wei_key, verbose = True)
    wei_hf = h5py.File(wei_dat, 'r')
    if wei_key == 'mf':
        wei_ant = wei_hf['evt_wise_ant'][:]
        weights = np.full((num_ants, num_evts), np.nan, dtype = float)
        weights[:8] = wei_ant[0, :8]
        weights[8:] = wei_ant[1, 8:]
        del wei_ant
    else:
        weights = wei_hf['snr'][:]
    del run_info1, wei_key, wei_dat, wei_hf

    print('example event number')
    print('entries, events, trigs')
    for e in range(20):
        print(entry_num[e], evt_num[e], trig_type[e])
    if sel_evts is not None:
        sel_evts = sel_evts.split(',')
        sel_evts = np.asarray(sel_evts).astype(int)
        sel_evt_idx = np.in1d(evt_num, sel_evts)
        sel_entries = entry_num[sel_evt_idx]
        sel_evts = evt_num[sel_evt_idx]
        sel_trig = trig_type[sel_evt_idx]
    else:
        sel_entries = entry_num[:20]
        sel_evts = evt_num[sel_entries]
        sel_trig = trig_type[sel_entries]
    print(f'Selected events are {sel_evts}')
    print(f'Selected entries are {sel_entries}')
    print(f'Selected triggers are {sel_trig}')
    sel_evt_len = len(sel_entries)

    # wf analyzer
    wf_int = wf_analyzer(use_time_pad = True, use_freq_pad = True, use_band_pass = True, use_rfft = True, use_cw = True, baseline = baseline_orig, testbed_dB = 6)
    dt = np.asarray([wf_int.dt])
    print(dt[0])
    pad_time = wf_int.pad_zero_t
    pad_len = wf_int.pad_len
    pad_freq = wf_int.pad_zero_freq
    pad_fft_len = wf_int.pad_fft_len
    df = wf_int.df

    # cw detection
    cw_phase1 = py_phase_variance(st, run, pad_freq, use_debug = True)
    evt_len = cw_phase1.evt_len
    start_evt = int(evt_len - 1)
    phase_var_freq_range = cw_phase1.useful_freq_range
    cw_testbed = py_testbed(st, run, pad_freq, 12, 11, 3, analyze_blind_dat = analyze_blind_dat, verbose = True, use_debug = True)
    baseline = cw_testbed.baseline_copy
    baseline_fft = cw_testbed.baseline_fft_copy
    baseline_fft_medi = cw_testbed.baseline_fft_medi_copy
    testbed_freq_range = cw_testbed.useful_freq_range

    # interferometers
    ara_int = py_interferometers(pad_len, dt[0], st, year, run = run, get_sub_file = True)
    pairs = ara_int.pairs
    v_pairs_len = ara_int.v_pairs_len    
    lags = ara_int.lags
    wei_pairs = get_products(weights, pairs, v_pairs_len) 

    # output array
    blk_max = 48
    samp_pad = blk_max * num_samps
    wf_all = np.full((samp_pad, 2, num_ants, sel_evt_len), np.nan, dtype=float)
    int_wf_all = np.full((wf_int.pad_len, 2, num_ants, sel_evt_len), np.nan, dtype=float)
    bp_wf_all = np.copy(int_wf_all)
    cw_wf_all = np.copy(int_wf_all)
    cw_bp_wf_all = np.copy(int_wf_all)
    ele_wf_all = np.full((samp_pad, 2, num_eles, sel_evt_len), np.nan, dtype=float)
    int_ele_wf_all = np.full((wf_int.pad_len, 2, num_eles, sel_evt_len), np.nan, dtype=float)
    adc_all = np.full((samp_pad, 2, num_ants, sel_evt_len), np.nan, dtype=float)
    ped_all = np.copy(adc_all)
    ele_adc_all = np.full((samp_pad, 2, num_eles, sel_evt_len), np.nan, dtype=float)
    ele_ped_all = np.copy(ele_adc_all)
    freq = np.full((wf_int.pad_fft_len, num_ants, sel_evt_len), np.nan, dtype=float)
    fft = np.copy(freq)
    bp_fft = np.copy(freq)
    cw_fft = np.copy(freq)
    cw_bp_fft = np.copy(freq)
    ele_freq = np.full((wf_int.pad_fft_len, num_eles, sel_evt_len), np.nan, dtype=float)
    ele_fft = np.copy(ele_freq)
    phase = np.copy(freq)
    bp_phase = np.copy(freq)
    cw_phase = np.copy(freq)
    cw_bp_phase = np.copy(freq)
    ele_phase = np.copy(ele_freq)
    bad_pad = 1000
    fft_dB = np.copy(freq)
    fft_dB_tilt = np.full((len(testbed_freq_range), num_ants, sel_evt_len), np.nan, dtype=float)
    baseline_tilt = np.copy(fft_dB_tilt)
    delta_mag = np.copy(fft_dB_tilt)
    testbed_bad_freqs = np.full((len(testbed_freq_range), num_pols, sel_evt_len), np.nan, dtype=float)
    testbed_bad_idx = np.full((bad_pad, sel_evt_len), np.nan, dtype = float)
    phase_variance = np.full((len(phase_var_freq_range), len(pairs), sel_evt_len), np.nan, dtype=float)
    phase_var_median = np.full((len(pairs), sel_evt_len), np.nan, dtype=float)
    phase_var_sigma = np.copy(phase_var_median)
    sigma_variance = np.copy(phase_variance)
    sigma_variance_avg = np.full((len(phase_var_freq_range), num_pols, sel_evt_len), np.nan, dtype=float)
    phase_var_bad_freqs = np.full((len(phase_var_freq_range), sel_evt_len), np.nan, dtype=float)
    phase_var_bad_idx = np.copy(testbed_bad_idx)
    phase_var_bad_sigma = np.copy(testbed_bad_idx)
    blk_est_range = 50
    blk_idx = np.full((blk_est_range, sel_evt_len), np.nan, dtype=float)
    samp_idx = np.full((num_samps, blk_est_range, num_ants, sel_evt_len), np.nan, dtype=float)
    time_arr = np.copy(samp_idx)
    int_time_arr = np.copy(samp_idx)
    num_samps_in_blk = np.full((blk_est_range, num_ants, sel_evt_len), np.nan, dtype=float)
    num_int_samps_in_blk = np.copy(num_samps_in_blk)
    mean_blk = np.full((blk_est_range, num_ants, sel_evt_len), np.nan, dtype=float)
    int_mean_blk = np.copy(mean_blk)
    bp_mean_blk = np.copy(mean_blk)
    cw_mean_blk = np.copy(mean_blk)
    cw_bp_mean_blk = np.copy(mean_blk)
    corr = np.full((ara_int.lag_len, ara_int.pair_len, sel_evt_len), np.nan, dtype = float)
    corr_nonorm = np.copy(corr)
    corr_01 = np.copy(corr)
    bp_corr = np.copy(corr)
    bp_corr_nonorm = np.copy(corr)
    bp_corr_01 = np.copy(corr)    
    cw_corr = np.copy(corr)
    cw_corr_nonorm = np.copy(corr)
    cw_corr_01 = np.copy(corr)
    cw_bp_corr = np.copy(corr)
    cw_bp_corr_nonorm = np.copy(corr)
    cw_bp_corr_01 = np.copy(corr)
    coval = np.full(ara_int.table_shape, np.nan, dtype = float)
    coval = np.repeat(coval[:, :, :, :, :, np.newaxis], sel_evt_len, axis = 5)
    bp_coval = np.copy(coval)
    cw_coval = np.copy(coval)
    cw_bp_coval = np.copy(coval)
    sky_map = np.full((ara_int.table_shape[0], ara_int.table_shape[1], 2, 2, 2, sel_evt_len), np.nan, dtype = float) 
    bp_sky_map = np.copy(sky_map)
    cw_sky_map = np.copy(sky_map)
    cw_bp_sky_map = np.copy(sky_map)

    # cw detection
    evt_counts = 0
    rs_idx = np.logical_and(daq_qual_cut_sum == 0, trig_type != 1)
    rs_entry = entry_num[rs_idx]
    sel_idxs = np.in1d(rs_entry, sel_entries)
    rs_idxs = np.arange(len(rs_entry), dtype = int)
    rs_sel_idxs = rs_idxs[sel_idxs]
    if len(rs_sel_idxs) == 0:
        rs_entry_idx = np.full((0), False, dtype = bool)
    else:
        rs_entry_idx = np.logical_and(rs_idxs > np.nanmin(rs_sel_idxs) - 15, rs_idxs < np.nanmax(rs_sel_idxs) + 1)
    rs_entry = rs_entry[rs_entry_idx]
    rs_evt = evt_num[rs_idx]
    rs_evt = rs_evt[sel_idxs]
    print(f'rs_entry for phase vatriance: {len(rs_entry)}, {rs_entry}')
    for evt in tqdm(range(len(rs_entry))):
        # get entry and wf
        ara_root.get_entry(rs_entry[evt])
        ara_root.get_useful_evt(ara_root.cal_type.kLatestCalib)
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_zero_pad = True)
            ara_root.del_TGraph()
        ara_root.del_usefulEvt()
        wf_int.get_fft_wf(use_zero_pad = True, use_rfft = True, use_phase = True, use_abs = True, use_norm = True, use_dbmHz = True)
        rfft_phase = wf_int.pad_phase
        cw_phase1.get_phase_differences(rfft_phase, evt_counts % evt_len)
        rs_evts = np.where(sel_entries == rs_entry[evt])[0]
        if len(rs_evts) != 0:
            print('CW TestBed!!!!!!!!')
            rs_evts1 = rs_evts[0]
            rfft_dbmhz = wf_int.pad_fft
            fft_dB[:, :, rs_evts1] = np.copy(rfft_dbmhz)
            cw_testbed.get_bad_magnitude(rfft_dbmhz)
            fft_dB_t = cw_testbed.fft_dB_tilt
            baseline_t = cw_testbed.baseline_tilt
            delta_m = cw_testbed.delta_mag_copy
            testbed_bad_f = cw_testbed.bad_freq_pol_copy
            testbed_bad_i = cw_testbed.bad_idx
            bad_len = len(testbed_bad_i)
            fft_dB_tilt[:, :, rs_evts1] = fft_dB_t
            baseline_tilt[:, :, rs_evts1] = baseline_t
            delta_mag[:, :, rs_evts1] = delta_m
            testbed_bad_freqs[:, :, rs_evts1] = testbed_bad_f
            testbed_bad_idx[:bad_len, rs_evts1] = testbed_bad_i
            print(testbed_bad_i)
        if evt_counts < start_evt:
            pass
        else:
            print('CW Phase!!!!!!')
            cw_phase1.get_bad_phase()
            phase_v = cw_phase1.phase_variance_copy
            phase_v_m = cw_phase1.median_copy
            phase_v_s = cw_phase1.sigma_copy
            sigma_v = cw_phase1.sigma_variance_copy
            sigma_v_a = cw_phase1.sigma_variance_avg
            phase_v_bad_f = cw_phase1.bad_bool_copy
            sigmas = cw_phase1.bad_sigma
            phase_idxs = cw_phase1.bad_idx
            print(phase_idxs)

            if st == 2: sigma_thres = 1.5
            if st == 3: sigma_thres = 2

            sigmas_bad = sigmas > sigma_thres
            sigmas = sigmas[sigmas_bad]
            phase_idxs = phase_idxs[sigmas_bad] 

            bad_len = len(phase_idxs)
            phase_variance[:, :, rs_evts1] = phase_v
            phase_var_median[:, rs_evts1] = phase_v_m
            phase_var_sigma[:, rs_evts1] = phase_v_s
            sigma_variance[:, :, rs_evts1] = sigma_v
            sigma_variance_avg[:, :, rs_evts1] = sigma_v_a
            phase_var_bad_freqs[:, rs_evts1] = phase_v_bad_f
            phase_var_bad_idx[:bad_len, rs_evts1] = phase_idxs
            phase_var_bad_sigma[:bad_len, rs_evts1] = sigmas
            print(phase_idxs)
        evt_counts += 1

    bad_idxs = np.append(testbed_bad_idx, phase_var_bad_idx, axis = 0)

    # loop over the events
    for evt in tqdm(range(sel_evt_len)):
    
        # sample info and timing info
        blk_idx_arr, blk_idx_len = ara_uproot.get_block_idx(sel_entries[evt], trim_1st_blk = True)
        blk_idx[:blk_idx_len, evt] = blk_idx_arr
        buffer_info.get_num_samp_in_blk(blk_idx_arr)
        buffer_info.get_num_samp_in_blk(blk_idx_arr, use_int_dat = True)
        num_samps_in_blk[:blk_idx_len, :, evt] = buffer_info.samp_in_blk
        num_int_samps_in_blk[:blk_idx_len, :, evt] = buffer_info.int_samp_in_blk
        samp_idx[:, :blk_idx_len, :, evt] = buffer_info.get_samp_idx(blk_idx_arr)
        time_arr[:, :blk_idx_len, :, evt] = buffer_info.get_time_arr(blk_idx_arr, trim_1st_blk = True)
        int_time_arr[:int_num_samples, :blk_idx_len, :, evt] = buffer_info.get_time_arr(blk_idx_arr, trim_1st_blk = True, use_int_dat = True)
        
        # get entry and wf
        ara_root.get_entry(sel_entries[evt])

        # raw wfs 
        ara_root.get_useful_evt(ara_root.cal_type.kOnlyADCWithOut1stBlock)
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_len = len(raw_t)
            adc_all[:wf_len, 0, ant, evt] = raw_t
            adc_all[:wf_len, 1, ant, evt] = raw_v
            ara_root.del_TGraph()
        for ant in range(num_eles):
            raw_t, raw_v = ara_root.get_ele_ch_wf(ant)
            wf_len = len(raw_t)
            ele_adc_all[:wf_len, 0, ant, evt] = raw_t
            ele_adc_all[:wf_len, 1, ant, evt] = raw_v
            ara_root.del_TGraph()
        ara_root.del_usefulEvt()

        # pedestal
        ara_root.get_useful_evt(ara_root.cal_type.kOnlyPedWithOut1stBlock)
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_len = len(raw_t)
            ped_all[:wf_len, 0, ant, evt] = raw_t
            ped_all[:wf_len, 1, ant, evt] = raw_v
            ara_root.del_TGraph()
        for ant in range(num_eles):
            raw_t, raw_v = ara_root.get_ele_ch_wf(ant)
            wf_len = len(raw_t)
            ele_ped_all[:wf_len, 0, ant, evt] = raw_t
            ele_ped_all[:wf_len, 1, ant, evt] = raw_v
            ara_root.del_TGraph()
        ara_root.del_usefulEvt()

        # calibrated wf, interpolated wf and mean of wf 
        ara_root.get_useful_evt(ara_root.cal_type.kLatestCalib)
        for ant in range(num_ants):        
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_len = len(raw_t)
            wf_all[:wf_len, 0, ant, evt] = raw_t
            wf_all[:wf_len, 1, ant, evt] = raw_v
            mean_blk[:blk_idx_len, ant, evt] = buffer_info.get_mean_blk(ant, raw_v)
            wf_int.get_int_wf(raw_t, raw_v, ant)
            wf_int_len = wf_int.pad_num[ant]
            int_v = wf_int.pad_v[:wf_int_len, ant] 
            int_mean_blk[:blk_idx_len, ant, evt] = buffer_info.get_mean_blk(ant, int_v, use_int_dat = True) 
            ara_root.del_TGraph()
        int_wf_all[:, 0, :, evt] = wf_int.pad_t
        int_wf_all[:, 1, :, evt] = wf_int.pad_v
        wf_int.get_fft_wf(use_rfft = True, use_abs = True, use_norm = True, use_phase = True)
        freq[:, :, evt] = wf_int.pad_freq
        fft[:, :, evt] = wf_int.pad_fft
        phase[:, :, evt] = wf_int.pad_phase
     
        # band-pass wf  
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_band_pass = True)
            wf_int_len = wf_int.pad_num[ant]
            bp_v = wf_int.pad_v[:wf_int_len, ant]
            bp_mean_blk[:blk_idx_len, ant, evt] = buffer_info.get_mean_blk(ant, bp_v, use_int_dat = True)
            ara_root.del_TGraph()
        bp_wf_all[:, 0, :, evt] = wf_int.pad_t
        bp_wf_all[:, 1, :, evt] = wf_int.pad_v
        wf_int.get_fft_wf(use_rfft = True, use_abs = True, use_norm = True, use_phase = True)        
        bp_fft[:, :, evt] = wf_int.pad_fft
        bp_phase[:, :, evt] = wf_int.pad_phase
        
        # cw wf
        bad_idx_evt = bad_idxs[:, evt]
        bad_idx_evt = bad_idx_evt[~np.isnan(bad_idx_evt)]
        bad_idx_evt = bad_idx_evt.astype(int)
        bad_idx_evt = np.unique(bad_idx_evt).astype(int)
        
        df_10 = np.round(0.01 / df).astype(int)

        diff_idx = np.diff(bad_idx_evt)
        diff_idx = diff_idx > (df_10 - 1)
        diff_idx_len = np.count_nonzero(diff_idx) + 1
       
        bad_band = np.full((diff_idx_len, 3), 0, dtype = int)
        bad_band[0, 0] = bad_idx_evt[0]
        bad_band[-1, 1] = bad_idx_evt[-1]
        bad_band[1:, 0] = bad_idx_evt[1:][diff_idx]
        bad_band[:-1, 1] = bad_idx_evt[:-1][diff_idx]
        bad_band[:, 2] = np.round((bad_band[:, 1] + bad_band[:, 0]) / 2)
 
        bad_idx_evt = np.copy(bad_band[:, 2])
        print(bad_band[:, 0])
        print(bad_band[:, 1])
        print(bad_band[:, 2])
        freq_band = np.full((len(bad_idx_evt), 2), np.nan, dtype = float)   
        freq_band[:, 0] = pad_freq[bad_band[:, 0]] - 5
        freq_band[:, 1] = pad_freq[bad_band[:, 1]] + 5

        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_cw = True, bad_idx = bad_idx_evt, band = freq_band)
            wf_int_len = wf_int.pad_num[ant]
            cw_v = wf_int.pad_v[:wf_int_len, ant]
            cw_mean_blk[:blk_idx_len, ant, evt] = buffer_info.get_mean_blk(ant, cw_v, use_int_dat = True)
            ara_root.del_TGraph()
        cw_wf_all[:, 0, :, evt] = wf_int.pad_t
        cw_wf_all[:, 1, :, evt] = wf_int.pad_v
        wf_int.get_fft_wf(use_rfft = True, use_abs = True, use_norm = True, use_phase = True)
        cw_fft[:, :, evt] = wf_int.pad_fft
        cw_phase[:, :, evt] = wf_int.pad_phase

        # cw (and band-passed) wf
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_band_pass = True, use_cw = True, bad_idx = bad_idx_evt, band = freq_band)
            wf_int_len = wf_int.pad_num[ant]
            cw_bp_v = wf_int.pad_v[:wf_int_len, ant]
            cw_bp_mean_blk[:blk_idx_len, ant, evt] = buffer_info.get_mean_blk(ant, cw_bp_v, use_int_dat = True)
            ara_root.del_TGraph()
        cw_bp_wf_all[:, 0, :, evt] = wf_int.pad_t
        cw_bp_wf_all[:, 1, :, evt] = wf_int.pad_v
        wf_int.get_fft_wf(use_rfft = True, use_abs = True, use_norm = True, use_phase = True)
        cw_bp_fft[:, :, evt] = wf_int.pad_fft
        cw_bp_phase[:, :, evt] = wf_int.pad_phase

        # reco w/ interpolated wf
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_zero_pad = True)
            ara_root.del_TGraph()
        ara_int.get_sky_map(wf_int.pad_v, weights = wei_pairs[:, evt], sum_pol = False, return_debug_dat = True)
        corr[:,:,evt] = ara_int.corr
        corr_nonorm[:,:,evt] = ara_int.corr_nonorm
        corr_01[:,:,evt] = ara_int.nor_fac
        coval_evt = ara_int.coval
        coval[:,:,:,:,:,evt] = coval_evt
        sky_map[:,:,:,:,0,evt] = np.nansum(coval_evt[:, :, :, :, :v_pairs_len], axis = 4)
        sky_map[:,:,:,:,1,evt] = np.nansum(coval_evt[:, :, :, :, v_pairs_len:], axis = 4)

        # reco w/band-passed wf
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_zero_pad = True, use_band_pass = True)
            ara_root.del_TGraph()
        ara_int.get_sky_map(wf_int.pad_v, weights = wei_pairs[:, evt], sum_pol = False, return_debug_dat = True)
        bp_corr[:,:,evt] = ara_int.corr
        bp_corr_nonorm[:,:,evt] = ara_int.corr_nonorm
        bp_corr_01[:,:,evt] = ara_int.nor_fac
        bp_coval_evt = ara_int.coval
        bp_coval[:,:,:,:,:,evt] = bp_coval_evt
        bp_sky_map[:,:,:,:,0,evt] = np.nansum(bp_coval_evt[:, :, :, :, :v_pairs_len], axis = 4)
        bp_sky_map[:,:,:,:,1,evt] = np.nansum(bp_coval_evt[:, :, :, :, v_pairs_len:], axis = 4)
        
        
        # reco w/ cw wf
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_zero_pad = True, use_cw = True, bad_idx = bad_idx_evt, band = freq_band)
            ara_root.del_TGraph()
        ara_int.get_sky_map(wf_int.pad_v, weights = wei_pairs[:, evt], sum_pol = False, return_debug_dat = True)
        cw_corr[:,:,evt] = ara_int.corr
        cw_corr_nonorm[:,:,evt] = ara_int.corr_nonorm
        cw_corr_01[:,:,evt] = ara_int.nor_fac
        cw_coval_evt = ara_int.coval
        cw_coval[:,:,:,:,:,evt] = cw_coval_evt
        cw_sky_map[:,:,:,:,0,evt] = np.nansum(cw_coval_evt[:, :, :, :, :v_pairs_len], axis = 4)
        cw_sky_map[:,:,:,:,1,evt] = np.nansum(cw_coval_evt[:, :, :, :, v_pairs_len:], axis = 4)
        
        # reco w/ cw (and band-passed) wf
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_zero_pad = True, use_band_pass = True, use_cw = True, bad_idx = bad_idx_evt, band = freq_band)
            ara_root.del_TGraph()
        ara_int.get_sky_map(wf_int.pad_v, weights = wei_pairs[:, evt], sum_pol = False, return_debug_dat = True)
        cw_bp_corr[:,:,evt] = ara_int.corr
        cw_bp_corr_nonorm[:,:,evt] = ara_int.corr_nonorm
        cw_bp_corr_01[:,:,evt] = ara_int.nor_fac
        cw_bp_coval_evt = ara_int.coval
        cw_bp_coval[:,:,:,:,:,evt] = cw_bp_coval_evt
        cw_bp_sky_map[:,:,:,:,0,evt] = np.nansum(cw_bp_coval_evt[:, :, :, :, :v_pairs_len], axis = 4)
        cw_bp_sky_map[:,:,:,:,1,evt] = np.nansum(cw_bp_coval_evt[:, :, :, :, v_pairs_len:], axis = 4)

    # interpolated all ele chs
    wf_int = wf_analyzer(use_time_pad = True, use_freq_pad = True, use_rfft = True, use_ele_ch = True)
    # loop over the events
    for evt in tqdm(range(sel_evt_len)):
        # get entry and wf
        ara_root.get_entry(sel_entries[evt])
        ara_root.get_useful_evt(ara_root.cal_type.kLatestCalib)
        for ant in range(num_eles):
            raw_t, raw_v = ara_root.get_ele_ch_wf(ant)
            wf_len = len(raw_t)
            ele_wf_all[:wf_len, 0, ant, evt] = raw_t
            ele_wf_all[:wf_len, 1, ant, evt] = raw_v
            wf_int.get_int_wf(raw_t, raw_v, ant)
            ara_root.del_TGraph()
        int_ele_wf_all[:, 0, :, evt] = wf_int.pad_t
        int_ele_wf_all[:, 1, :, evt] = wf_int.pad_v
        wf_int.get_fft_wf(use_rfft = True, use_abs = True, use_norm = True, use_phase = True)
        ele_freq[:, :, evt] = wf_int.pad_freq
        ele_fft[:, :, evt] = wf_int.pad_fft
        ele_phase[:, :, evt] = wf_int.pad_phase
        ara_root.del_usefulEvt()

    print('WF collecting is done!')

    #output
    return {'evt_num':evt_num,
            'entry_num':entry_num,
            'evt_full':evt_full,
            'daq_qual_cut_sum':daq_qual_cut_sum,
            'tot_qual_cut':tot_qual_cut,
            'weights':weights,
            'wei_pairs':wei_pairs,
            'trig_type':trig_type,
            'pps_number':pps_number,
            'unix_time':unix_time,
            'run_info':run_info,
            'rf_ch':rf_ch,
            'ele_ch':ele_ch,
            'trig_ch':trig_ch,
            'bad_ant':bad_ant,
            'sel_entries':sel_entries,
            'sel_evts':sel_evts,
            'sel_trig':sel_trig,
            'dt':dt,
            'pad_time':pad_time,
            'pad_freq':pad_freq,
            'pairs':pairs,
            'lags':lags,
            'wf_all':wf_all,
            'int_wf_all':int_wf_all,
            'bp_wf_all':bp_wf_all,
            'cw_wf_all':cw_wf_all,
            'cw_bp_wf_all':cw_bp_wf_all,
            'ele_wf_all':ele_wf_all,
            'int_ele_wf_all':int_ele_wf_all,
            'adc_all':adc_all,
            'ped_all':ped_all,
            'ele_adc_all':ele_adc_all,
            'ele_ped_all':ele_ped_all,
            'freq':freq,
            'fft':fft,
            'bp_fft':bp_fft,
            'cw_fft':cw_fft,
            'cw_bp_fft':cw_bp_fft,
            'ele_freq':ele_freq,
            'ele_fft':ele_fft,
            'phase':phase,
            'bp_phase':bp_phase,
            'cw_phase':cw_phase,
            'cw_bp_phase':cw_bp_phase,
            'ele_phase':ele_phase,
            'rs_entry':rs_entry,
            'rs_evt':rs_evt,
            'phase_var_freq_range':phase_var_freq_range,
            'baseline':baseline,
            'baseline_fft':baseline_fft,
            'baseline_fft_medi':baseline_fft_medi,
            'testbed_freq_range':testbed_freq_range,
            'fft_dB':fft_dB,
            'fft_dB_tilt':fft_dB_tilt,
            'baseline_tilt':baseline_tilt,
            'delta_mag':delta_mag,
            'testbed_bad_freqs':testbed_bad_freqs,
            'testbed_bad_idx':testbed_bad_idx,
            'phase_variance':phase_variance,
            'phase_var_median':phase_var_median,
            'phase_var_sigma':phase_var_sigma,
            'sigma_variance':sigma_variance,
            'sigma_variance_avg':sigma_variance_avg,
            'phase_var_bad_freqs':phase_var_bad_freqs,
            'phase_var_bad_idx':phase_var_bad_idx,
            'phase_var_bad_sigma':phase_var_bad_sigma,
            'blk_idx':blk_idx,
            'samp_idx':samp_idx,
            'time_arr':time_arr,
            'int_time_arr':int_time_arr,
            'num_samps_in_blk':num_samps_in_blk,
            'num_int_samps_in_blk':num_int_samps_in_blk,
            'mean_blk':mean_blk,
            'int_mean_blk':int_mean_blk,
            'bp_mean_blk':bp_mean_blk,
            'cw_mean_blk':cw_mean_blk,
            'cw_bp_mean_blk':cw_bp_mean_blk,
            'corr':corr,
            'bp_corr':bp_corr,
            'cw_corr':cw_corr,
            'cw_bp_corr':cw_bp_corr,
            'corr_nonorm':corr_nonorm,
            'bp_corr_nonorm':bp_corr_nonorm,
            'cw_corr_nonorm':cw_corr_nonorm,
            'cw_bp_corr_nonorm':cw_bp_corr_nonorm,
            'corr_01':corr_01,
            'bp_corr_01':bp_corr_01,
            'cw_corr_01':cw_corr_01,
            'cw_bp_corr_01':cw_bp_corr_01,
            'coval':coval,
            'bp_coval':bp_coval,
            'cw_coval':cw_coval,
            'cw_bp_coval':cw_bp_coval,
            'sky_map':sky_map,
            'bp_sky_map':bp_sky_map,
            'cw_sky_map':cw_sky_map,
            'cw_bp_sky_map':cw_bp_sky_map}























