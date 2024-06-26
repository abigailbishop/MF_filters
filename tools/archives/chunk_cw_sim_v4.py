import numpy as np
from tqdm import tqdm

def cw_sim_collector(Data, Station, Year):

    print('Collecting cw starts!')

    from tools.ara_sim_load import ara_root_loader
    from tools.ara_constant import ara_const
    from tools.ara_wf_analyzer import wf_analyzer
    from tools.ara_wf_analyzer import hist_loader

    # geom. info.
    ara_const = ara_const()
    num_ants = ara_const.USEFUL_CHAN_PER_STATION 
    del ara_const

    # data config
    ara_root = ara_root_loader(Data, Station, Year)
    num_evts = ara_root.num_evts
    evt_num = np.arange(num_evts, dtype = int)
    entry_num = np.copy(evt_num)
    unix_time = evt_num * 60
    unix_min_bins = np.append(unix_time, unix_time[-1] + 60)

    # wf analyzer
    cw_params = np.full((num_ants), 0.02, dtype = float)
    wf_int = wf_analyzer(use_time_pad = True, use_freq_pad = True, use_rfft = True, use_band_pass = True, use_cw = True, cw_params = cw_params)
    freq_range = wf_int.pad_zero_freq
    freq_bins = np.fft.rfftfreq(200, wf_int.dt)
    amp_range = np.arange(-5, 5, 0.05)
    amp_bins = np.linspace(-5, 5, 200 + 1)
    ara_hist = hist_loader(freq_bins, amp_bins)
    freq_bin_center = ara_hist.bin_x_center
    amp_bin_center = ara_hist.bin_y_center

    # output
    freq_bin_len = len(freq_bin_center)
    amp_bin_len = len(amp_bin_center)
    fft_rf_cut_map = np.full((freq_bin_len, amp_bin_len, num_ants), 0, dtype = int)
    sol_pad = 200
    sub_freq = np.full((sol_pad, num_ants, num_evts), np.nan, dtype = float)
    sub_amp_err = np.copy(sub_freq)
    sub_ratio = np.copy(sub_freq)
    sub_power = np.copy(sub_freq)
    sub_amp_err[0] = 0
    sub_ratio[0] = 0
    del freq_bin_len, amp_bin_len

    # loop over the events
    for evt in tqdm(range(num_evts)):

        # get entry and wf
        ara_root.get_entry(evt)
        
        # loop over the antennas
        for ant in range(num_ants):
            raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
            wf_int.get_int_wf(raw_t, raw_v, ant, use_zero_pad = True, use_band_pass = True, use_cw = True)
            num_sols = wf_int.sin_sub.num_sols + 1
            sub_freq[1:num_sols, ant, evt] = wf_int.sin_sub.sub_freqs
            sub_amp_err[1:num_sols, ant, evt] = wf_int.sin_sub.sub_amp_errs
            sub_ratio[:num_sols, ant, evt] = wf_int.sin_sub.sub_ratios
            sub_power[:num_sols, ant, evt] = wf_int.sin_sub.sub_powers
            del raw_t, raw_v, num_sols 
            ara_root.del_TGraph()

        wf_int.get_fft_wf(use_zero_pad = True, use_rfft = True, use_abs = True)
        fft_evt = np.log10(wf_int.pad_fft)      
        for ant in range(num_ants):
            fft_rf_cut_map[:, :, ant] += np.histogram2d(freq_range, fft_evt[:, ant], bins = (freq_bins, amp_bins))[0].astype(int)        
        del fft_evt
    del ara_root, wf_int 

    sub_sum = np.nansum(sub_power, axis = 0)
    sub_weight = sub_power / sub_sum[np.newaxis, :, :]
    del sub_sum

    print('power')
    power_range = np.arange(0, 2000, 4)
    power_bins = np.linspace(0, 2000, 500 + 1)
    ara_hist = hist_loader(power_bins)
    power_bin_center = ara_hist.bin_x_center
    power_rf_cut_hist = ara_hist.get_1d_hist(sub_power, weight = sub_weight, use_flat = True)
    del ara_hist

    print('ratio')
    ratio_range = np.arange(0, 1, 0.02)
    ratio_bins = np.linspace(0, 1, 50 + 1)
    ara_hist = hist_loader(ratio_bins)
    ratio_bin_center = ara_hist.bin_x_center    
    ratio_rf_cut_hist = ara_hist.get_1d_hist(sub_ratio, weight = sub_weight, use_flat = True)
    del ara_hist

    print('amp error')
    amp_err_range = np.arange(0, 150, 1) 
    amp_err_bins = np.linspace(0, 150, 150 + 1)
    ara_hist = hist_loader(amp_err_bins)
    amp_err_bin_center = ara_hist.bin_x_center
    amp_err_rf_cut_hist = ara_hist.get_1d_hist(sub_amp_err, weight = sub_weight, use_flat = True)
    del ara_hist

    print('2d')
    ara_hist = hist_loader(amp_err_bins, ratio_bins)
    amp_err_ratio_rf_cut_map = ara_hist.get_2d_hist(sub_amp_err, sub_ratio, weight = sub_weight, use_flat = True)
    del ara_hist

    print('time')
    unix_min_range = np.copy(unix_time)
    ara_hist = hist_loader(unix_min_bins, ratio_bins)
    unix_min_bin_center = ara_hist.bin_x_center

    clean_unix = np.copy(unix_time)
    if len(clean_unix) == 0:
        clean_unix_ant = np.full((num_ants, len(clean_unix)), np.nan, dtype = float)
        clean_unix_all = np.full((sol_pad, num_ants, len(clean_unix)), np.nan, dtype = float)
    else:
        clean_unix_ant = np.repeat(clean_unix[np.newaxis, :], num_ants, axis = 0)
        clean_unix_all = np.repeat(clean_unix_ant[np.newaxis, :, :], sol_pad, axis = 0)

    unix_ratio_rf_cut_map = ara_hist.get_2d_hist(clean_unix_all, sub_ratio, weight = sub_weight, use_flat = True)
    unix_ratio_rf_cut_map_max = ara_hist.get_2d_hist_max(unix_ratio_rf_cut_map)
    del ara_hist, sol_pad

    ara_hist = hist_loader(unix_min_bins, power_bins)   
    unix_power_rf_cut_map = ara_hist.get_2d_hist(clean_unix_all, sub_power, weight = sub_weight, use_flat = True)
    unix_power_rf_cut_map_max = ara_hist.get_2d_hist_max(unix_power_rf_cut_map)
    del ara_hist

    ara_hist = hist_loader(unix_min_bins, amp_err_bins)
    unix_amp_err_rf_cut_map = ara_hist.get_2d_hist(clean_unix_all, sub_amp_err, weight = sub_weight, use_flat = True)
    del ara_hist

    ara_hist = hist_loader(unix_min_bins, freq_bins)
    unix_freq_rf_cut_map = ara_hist.get_2d_hist(clean_unix_all, sub_freq, weight = sub_weight, use_flat = True)
    temp_max_idx = np.nanargmax(unix_freq_rf_cut_map, axis = 1)
    unix_freq_rf_cut_map_max = freq_bin_center[temp_max_idx]
    del ara_hist, clean_unix_ant, clean_unix_all, num_ants, temp_max_idx

    print('cw collecting is done!')

    return {'evt_num':evt_num,
            'unix_time':unix_time,
            'clean_unix':clean_unix,
            'sub_freq':sub_freq,
            'sub_amp_err':sub_amp_err,
            'sub_power':sub_power,
            'sub_ratio':sub_ratio,
            'sub_weight':sub_weight,
            'unix_min_range':unix_min_range,
            'unix_min_bins':unix_min_bins,
            'unix_min_bin_center':unix_min_bin_center,
            'freq_range':freq_range,
            'freq_bins':freq_bins,
            'freq_bin_center':freq_bin_center,
            'power_range':power_range,
            'power_bins':power_bins,
            'power_bin_center':power_bin_center,
            'ratio_range':ratio_range,
            'ratio_bins':ratio_bins,
            'ratio_bin_center':ratio_bin_center,
            'amp_range':amp_range,
            'amp_bins':amp_bins,
            'amp_bin_center':amp_bin_center,
            'amp_err_range':amp_err_range,
            'amp_err_bins':amp_err_bins,
            'amp_err_bin_center':amp_err_bin_center,
            'fft_rf_cut_map':fft_rf_cut_map,
            'power_rf_cut_hist':power_rf_cut_hist,
            'ratio_rf_cut_hist':ratio_rf_cut_hist,
            'amp_err_rf_cut_hist':amp_err_rf_cut_hist,
            'amp_err_ratio_rf_cut_map':amp_err_ratio_rf_cut_map,
            'unix_ratio_rf_cut_map':unix_ratio_rf_cut_map,
            'unix_ratio_rf_cut_map_max':unix_ratio_rf_cut_map_max,
            'unix_power_rf_cut_map':unix_power_rf_cut_map,
            'unix_power_rf_cut_map_max':unix_power_rf_cut_map_max,
            'unix_amp_err_rf_cut_map':unix_amp_err_rf_cut_map,
            'unix_freq_rf_cut_map':unix_freq_rf_cut_map,
            'unix_freq_rf_cut_map_max':unix_freq_rf_cut_map_max}



