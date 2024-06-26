import os, sys
import numpy as np
from tqdm import tqdm
import h5py

# custom lib
from tools.ara_constant import ara_const
from tools.ara_run_manager import run_info_loader

ara_const = ara_const()
num_ddas = ara_const.DDA_PER_ATRI
num_blks = ara_const.BLOCKS_PER_DDA
num_ants = ara_const.USEFUL_CHAN_PER_STATION
num_eles = ara_const.CHANNELS_PER_ATRI
num_samps = ara_const.SAMPLES_PER_BLOCK
num_chs = ara_const.RFCHAN_PER_DDA

def quick_qual_check(dat_bool, dat_idx, ser_val):

    bool_len = np.count_nonzero(dat_bool)
    if bool_len > 0:
        print(f'Qcut, {ser_val}:', bool_len, dat_idx[dat_bool])
    del bool_len

class pre_qual_cut_loader:

    def __init__(self, ara_uproot, analyze_blind_dat = False, verbose = False):

        self.st = ara_uproot.station_id
        self.run = ara_uproot.run
        self.evt_num = ara_uproot.evt_num 
        self.num_evts = ara_uproot.num_evts
        self.trig_type = ara_uproot.get_trig_type()
        self.unix_time = ara_uproot.unix_time
        self.pps_number = ara_uproot.pps_number
        self.irs_block_number = ara_uproot.irs_block_number
        self.channel_mask = ara_uproot.channel_mask
        self.blk_len = ara_uproot.read_win//num_ddas
        self.verbose = verbose

        self.run_info = run_info_loader(self.st, self.run, analyze_blind_dat = analyze_blind_dat)
        evt_rate_dat = self.run_info.get_result_path(file_type = 'evt_rate', verbose = self.verbose, force_blind = True)
        self.evt_rate_hf = h5py.File(evt_rate_dat, 'r')
        self.evt_sort = self.evt_rate_hf['evt_num_sort'][:]
        self.unix_sort = self.evt_rate_hf['unix_time_sort'][:]
        self.pps_num = self.evt_rate_hf['pps_number_sort_reset'][:]
        del evt_rate_dat

    def get_daq_structure_errors(self):

        bi_ch_mask = 1 << np.arange(num_chs, dtype = int)
        dda_ch = np.arange(num_ddas, dtype = int)
        dda_idx = (self.channel_mask & 0x300) >> 8

        daq_st_err = np.full((self.num_evts, 5), 0, dtype = int)
        # bad block lengh
        # bad block index
        # block gap
        # bad dda index
        # bad channel mask

        for evt in range(self.num_evts):

            # bad block length
            blk_idx_evt = np.asarray(self.irs_block_number[evt], dtype = int)
            daq_st_err[evt, 0] = len(blk_idx_evt) % num_ddas
            if daq_st_err[evt, 0] != 0:
                continue

            # bad block index
            blk_idx_reshape = np.reshape(blk_idx_evt, (-1, num_ddas))
            daq_st_err[evt, 1] = int(np.any(blk_idx_reshape != blk_idx_reshape[:,0][:, np.newaxis]))
            del blk_idx_evt

            # block gap
            for dda in range(num_ddas):
                blk_idx_dda = blk_idx_reshape[:, dda]
                first_block_idx = blk_idx_dda[0]
                last_block_idx = blk_idx_dda[-1]
                block_diff = len(blk_idx_dda) - 1

                if first_block_idx + block_diff != last_block_idx:
                    if num_blks - first_block_idx + last_block_idx != block_diff:
                        daq_st_err[evt, 2] += 1
                del first_block_idx, last_block_idx, block_diff, blk_idx_dda
            del blk_idx_reshape

            # bad dda channel
            dda_idx_evt = np.asarray(dda_idx[evt], dtype = int)
            dda_idx_reshape = np.reshape(dda_idx_evt, (-1, num_ddas))
            daq_st_err[evt, 3] = int(np.any(dda_idx_reshape != dda_ch[np.newaxis, :]))
            del dda_idx_reshape, dda_idx_evt

            # bad channel mask
            ch_mask_evt = np.asarray(self.channel_mask[evt], dtype = int)
            ch_mask_reshape = np.repeat(ch_mask_evt[:, np.newaxis], num_chs, axis = 1)
            ch_mask_bit = ch_mask_reshape & bi_ch_mask[np.newaxis, :]
            daq_st_err[evt, 4] = int(np.any(ch_mask_bit != bi_ch_mask[np.newaxis, :]))
            del ch_mask_reshape, ch_mask_bit, ch_mask_evt
        del bi_ch_mask, dda_ch, dda_idx

        if self.verbose:
            quick_qual_check(daq_st_err[:, 0] != 0, self.evt_num, 'bad block length events')
            quick_qual_check(daq_st_err[:, 1] != 0, self.evt_num, 'bad block index events')
            quick_qual_check(daq_st_err[:, 2] != 0, self.evt_num, 'block gap events')
            quick_qual_check(daq_st_err[:, 3] != 0, self.evt_num, 'bad dda index events')
            quick_qual_check(daq_st_err[:, 4] != 0, self.evt_num, 'bad channel mask events')

        return daq_st_err

    def get_readout_window_errors(self):

        rf_read_win_len, soft_read_win_len = self.get_read_win_limit()

        read_win_err = np.full((self.num_evts, 4), 0, dtype = int) 
        # single block
        # bad rf readout window
        # bad cal readout window
        # bad soft readout window
    
        read_win_err[:, 0] = (self.blk_len < 2).astype(int)
        read_win_err[:, 1] = np.logical_and(self.blk_len < rf_read_win_len, self.trig_type == 0).astype(int)
        read_win_err[:, 2] = np.logical_and(self.blk_len < rf_read_win_len, self.trig_type == 1).astype(int)
        read_win_err[:, 3] = np.logical_and(self.blk_len < soft_read_win_len, self.trig_type == 2).astype(int)
        del rf_read_win_len, soft_read_win_len

        if self.verbose:
            quick_qual_check(read_win_err[:, 0] != 0, self.evt_num, 'single block events')         
            quick_qual_check(read_win_err[:, 1] != 0, self.evt_num, 'bad rf readout window events')         
            quick_qual_check(read_win_err[:, 2] != 0, self.evt_num, 'bad cal readout window events')         
            quick_qual_check(read_win_err[:, 3] != 0, self.evt_num, 'bad soft readout window events')         

        return read_win_err

    def get_read_win_limit(self):

        if self.st == 2:
            if self.run < 4029:
                rf_readout_limit = 20
            elif self.run > 4028 and self.run < 9749:
                rf_readout_limit = 26
            elif self.run > 9748:
                rf_readout_limit = 28
        elif self.st == 3:
            if self.run < 3104:
                rf_readout_limit = 20
            elif self.run > 3103 and self.run < 10001:
                rf_readout_limit = 26
            elif self.run > 10000:
                rf_readout_limit = 28
            
        if self.st == 2:
            if self.run < 9505:
                soft_readout_limit = 8
            else:
                soft_readout_limit = 12
        elif self.st == 3:
            if self.run < 10001:
                soft_readout_limit = 8
            else:
                soft_readout_limit = 12

        return rf_readout_limit, soft_readout_limit

    def get_bad_unix_time_sequence(self):

        bad_unix_sequence = np.full((self.num_evts), 0, dtype = int)
        
        if self.st == 3 and self.run == 3461: # condamn this run...
            bad_unix_sequence[:] = 1
            if self.verbose:
                quick_qual_check(bad_unix_sequence != 0, self.evt_num, 'bad unix sequence')
            return bad_unix_sequence
        
        n_idxs = np.where(np.diff(self.unix_sort) < 0)[0]
        if len(n_idxs):
            tot_bad_evts = []
            for n_idx in n_idxs:
                n_idx = int(n_idx)
                unix_peak = self.unix_sort[:-1][n_idx]
                bad_idx = np.where(self.unix_sort[n_idx + 1:] < unix_peak + 1)
                bad_evts = self.evt_sort[n_idx + 1:][bad_idx]
                tot_bad_evts.extend(bad_evts)
                del unix_peak, bad_idx

            tot_bad_evts = np.asarray(tot_bad_evts)
            tot_bad_evts = np.unique(tot_bad_evts)
            bad_unix_sequence[:] = np.in1d(self.evt_num, tot_bad_evts).astype(int)
            del tot_bad_evts
        else:
            return bad_unix_sequence
        del n_idxs

        if self.verbose:
            quick_qual_check(bad_unix_sequence != 0, self.evt_num, 'bad unix sequence')

        return bad_unix_sequence

    def get_bad_unix_time_events(self, add_unchecked_unix_time = False):

        from tools.ara_known_issue import known_issue_loader
        ara_known_issue = known_issue_loader(self.st)

        bad_unix_evts = np.full((self.num_evts), 0, dtype = int)
        for evt in range(self.num_evts):
            bad_unix_evts[evt] = ara_known_issue.get_bad_unixtime(self.unix_time[evt])
        
        if add_unchecked_unix_time == True:
            for evt in range(self.num_evts):
                if ara_known_issue.get_unchecked_unixtime(self.unix_time[evt]):
                   bad_unix_evts[evt] = 1 
        del ara_known_issue
        
        if self.verbose:
            quick_qual_check(bad_unix_evts != 0, self.evt_num, 'bad unix time')

        return bad_unix_evts
        
    def get_first_minute_events(self, first_evt_limit = 7):

        unix_time_full = self.evt_rate_hf['unix_time'][:]
        unix_cut = unix_time_full[0] + 60
        first_min_evt_bools = (self.unix_time < unix_cut)
        del unix_cut, unix_time_full

        unix_sort_cut = self.unix_sort[0] + 60
        first_min_evt_sort = self.evt_sort[self.unix_sort < unix_sort_cut]
        first_min_evt_sort_bools = np.in1d(self.evt_num, first_min_evt_sort)
        del unix_sort_cut, first_min_evt_sort

        first_min_evts = np.logical_or(first_min_evt_bools, first_min_evt_sort_bools).astype(int)
        del first_min_evt_bools, first_min_evt_sort_bools

        if self.verbose:
            quick_qual_check(first_min_evts != 0, self.evt_num, f'first minute events')

        return first_min_evts

    def get_bias_voltage_events(self, volt_cut = [3, 3.5]):

        volt_cut = np.asarray(volt_cut, dtype = float)
        bias_volt_evts = np.full((self.num_evts), 0, dtype = int)

        sensor_dat = self.run_info.get_result_path(file_type = 'sensor', verbose = self.verbose, force_blind = True)
        sensor_hf = h5py.File(sensor_dat, 'r')
        sensor_unix = sensor_hf['unix_time'][:]
        if any(np.isnan(sensor_unix)):
            print('There is empty sensorHk file!')
            return bias_volt_evts
        sensor_unix_len = len(sensor_unix)
        if sensor_unix_len == 0:
            print('There is empty sensorHk file!')
            bias_volt_evts[:] = 1
            if self.verbose:
                quick_qual_check(bias_volt_evts != 0, self.evt_num, f'bias voltage events')
            return bias_volt_evts

        dda_volt = sensor_hf['dda_volt'][:]
        del sensor_dat, sensor_hf
        good_dda_bool = np.logical_and(dda_volt > volt_cut[0], dda_volt < volt_cut[1])
        if sensor_unix_len == 1:
            print('There is single sensorHk values!')
            dda_digi_idx = np.array([0], dtype = int)
            good_digi_bool = np.copy(good_dda_bool)
        else:
            dda_digi_idx = np.arange(sensor_unix_len, dtype = int)[1:]
            good_digi_bool = np.logical_and(good_dda_bool[1:], good_dda_bool[:-1])
        del dda_volt, good_dda_bool
 
        unix_digi = np.digitize(self.unix_time, sensor_unix) 
        for dda in range(num_ddas):
            good_digi_idx = dda_digi_idx[good_digi_bool[:, dda]]
            bias_volt_evts += np.in1d(unix_digi, good_digi_idx, invert =True).astype(int)
            del good_digi_idx
        bias_volt_evts[bias_volt_evts != 0] = 1
        del volt_cut, sensor_unix, sensor_unix_len, unix_digi, dda_digi_idx, good_digi_bool

        if self.verbose:
            quick_qual_check(bias_volt_evts != 0, self.evt_num, f'bias voltage events')

        return bias_volt_evts
   
    def get_no_calpulser_events(self, ratio_cut = 0.02, apply_bias_volt = None):
     
        no_cal_evts = np.full((self.num_evts), 0, dtype = int)
        if self.st == 3 and (self.run > 1124 and self.run < 1429):
            return no_cal_evts

        if apply_bias_volt is not None:
            trig_type_evt = self.trig_type[apply_bias_volt == 0]
        else:
            trig_type_evt = self.trig_type
       
        num_evts = len(trig_type_evt)
        if len(trig_type_evt) != 0:
            num_cal_evts = np.count_nonzero(trig_type_evt == 1)
            cal_evt_ratio = num_cal_evts / len(trig_type_evt)
        else:
            cal_evt_ratio = np.nan

        if cal_evt_ratio < ratio_cut:
            no_cal_evts[:] = 1

        if self.verbose:
            quick_qual_check(no_cal_evts != 0, self.evt_num, f'no calpulser events')

        return no_cal_evts

    def get_bad_trig_rate_events(self, rate, lower_cut, upper_cut = None, use_sec = False):

        if upper_cut == None:
            bad_rate_idx = rate < lower_cut
        else:
            bad_rate_idx = np.logical_or(rate < lower_cut, rate > upper_cut)

        if use_sec:
            sec_arr = np.arange(1, dtype = int)
        else:   
            sec_arr = np.arange(60, dtype = int)
        
        bad_sec = self.rate_bins[bad_rate_idx]
        bad_sec = np.repeat(bad_sec[:, np.newaxis], len(sec_arr), axis = 1)
        bad_sec += sec_arr[np.newaxis, :]
        bad_sec = bad_sec.flatten()
        del sec_arr, bad_rate_idx

        bad_pps_idx = np.in1d(self.pps_num, bad_sec)
        bad_evt_sort = self.evt_sort[bad_pps_idx]
        del bad_pps_idx

        return bad_evt_sort

    def get_bad_rate_events(self, use_sec = False):
    
        if use_sec:
            bin_type = 'sec'
            rf_rate_cut = 1
            cal_rate_cut = 1
            cal_upper_cut = 1
            soft_rate_cut = 0
            soft_upper_cut = 2
        else:
            bin_type = 'min'
            rf_rate_cut = 2.8 
            cal_rate_cut = 0.85
            cal_upper_cut = 1.1
            soft_rate_cut = 0.75
            soft_upper_cut = 1.1
            if self.st == 3 and self.run < 10001:
                rf_rate_cut = 4
            if self.st == 3 and self.run > 10000:
                rf_rate_cut = 2
            if self.st == 2 and (self.run > 6499 and self.run < 7175):
                cal_rate_cut = 0.75
                rf_rate_cut = 2.6
            if self.st == 3 and (self.run > 6001 and self.run < 6678):
                cal_rate_cut = 0.8 

        self.rate_bins = (self.evt_rate_hf[f'pps_{bin_type}_bins'][:-1] + 0.5).astype(int) # bin edge to corresponding minute
        rf_evt_rate = self.evt_rate_hf[f'rf_{bin_type}_rate_pps'][:]
        cal_evt_rate = self.evt_rate_hf[f'cal_{bin_type}_rate_pps'][:]
        soft_evt_rate = self.evt_rate_hf[f'soft_{bin_type}_rate_pps'][:]

        bad_rf_sort = self.get_bad_trig_rate_events(rf_evt_rate, rf_rate_cut, use_sec = use_sec)
        bad_cal_sort = self.get_bad_trig_rate_events(cal_evt_rate, cal_rate_cut, cal_upper_cut, use_sec = use_sec)
        bad_soft_sort = self.get_bad_trig_rate_events(soft_evt_rate, soft_rate_cut, soft_upper_cut, use_sec = use_sec)
        del cal_evt_rate, rf_evt_rate, soft_evt_rate

        bad_rate_evts = np.full((self.num_evts, 3), 0, dtype = int)
        bad_rate_evts[:, 0] = np.in1d(self.evt_num, bad_rf_sort).astype(int)
        bad_rate_evts[:, 1] = np.in1d(self.evt_num, bad_cal_sort).astype(int)
        bad_rate_evts[:, 2] = np.in1d(self.evt_num, bad_soft_sort).astype(int)
        del bad_rf_sort, bad_cal_sort, bad_soft_sort, self.rate_bins

        if self.st == 3 and (self.run > 1124 and self.run < 1429):
            bad_rate_evts[:, 1] = 0

        if self.verbose:
            quick_qual_check(bad_rate_evts[:, 0] != 0, self.evt_num, f'bad rf {bin_type} rate events')
            quick_qual_check(bad_rate_evts[:, 1] != 0, self.evt_num, f'bad calpulser {bin_type} rate events')
            quick_qual_check(bad_rate_evts[:, 2] != 0, self.evt_num, f'bad software {bin_type} rate events')

        return bad_rate_evts

    def run_pre_qual_cut(self):

        tot_pre_qual_cut = np.full((self.num_evts, 20), 0, dtype = int)
        tot_pre_qual_cut[:, :5] = self.get_daq_structure_errors()
        tot_pre_qual_cut[:, 5:9] = self.get_readout_window_errors()
        tot_pre_qual_cut[:, 9] = self.get_bad_unix_time_sequence()
        tot_pre_qual_cut[:, 10] = self.get_bad_unix_time_events(add_unchecked_unix_time = True)
        tot_pre_qual_cut[:, 11] = self.get_first_minute_events()
        tot_pre_qual_cut[:, 12] = self.get_bias_voltage_events()
        tot_pre_qual_cut[:, 13] = self.get_no_calpulser_events(apply_bias_volt = tot_pre_qual_cut[:,12])
        tot_pre_qual_cut[:, 14:17] = self.get_bad_rate_events()
        tot_pre_qual_cut[:, 17:] = self.get_bad_rate_events(use_sec = True)

        self.daq_qual_cut_sum = np.nansum(tot_pre_qual_cut[:, :6], axis = 1)
        self.pre_qual_cut_sum = np.nansum(tot_pre_qual_cut, axis = 1)

        if self.verbose:
            quick_qual_check(self.daq_qual_cut_sum != 0, self.evt_num, 'daq error cut!')
            quick_qual_check(self.pre_qual_cut_sum != 0, self.evt_num, 'total pre qual cut!')

        return tot_pre_qual_cut

class post_qual_cut_loader:

    def __init__(self, ara_uproot, ara_root, dt = 0.5, verbose = False):

        #from tools.ara_wf_analyzer import wf_analyzer
        #wf_int = wf_analyzer(dt = dt)
        #self.dt = wf_int.dt
        self.st = ara_uproot.station_id
        self.run = ara_uproot.run
        self.evt_num = ara_uproot.evt_num
        self.num_evts = ara_uproot.num_evts
        self.ara_root = ara_root
        self.verbose = verbose
 
        from tools.ara_known_issue import known_issue_loader
        ara_known_issue = known_issue_loader(self.st)
        self.bad_ant = ara_known_issue.get_bad_antenna(self.run)
        del ara_known_issue, ara_uproot#, wf_int

        # spare
        # cw (testbad, phase, anita)
        # spikey 
     
        self.unlock_cal_evts = np.full((self.num_evts), 0, dtype = int)

    def get_unlocked_calpulser_events(self, raw_v, cal_amp_limit = 2200):

        raw_v_max = np.nanmax(raw_v)
        unlock_cal_flag = int(raw_v_max > cal_amp_limit)
        del raw_v_max

        return unlock_cal_flag

    def run_post_qual_cut(self, evt):

        if self.st == 3 and (self.run > 1124 and self.run < 1429):

            self.ara_root.get_entry(evt)
            self.ara_root.get_useful_evt(self.ara_root.cal_type.kOnlyADCWithOut1stBlockAndBadSamples)
            raw_v = self.ara_root.get_rf_ch_wf(2)[1]   
            self.unlock_cal_evts[evt] = self.get_unlocked_calpulser_events(raw_v)    
            del raw_v
            self.ara_root.del_TGraph()
            self.ara_root.del_usefulEvt()
    
    def get_channel_cerrelation_flag(self, dat, ant_limit = 2, st_limit = 1, apply_bad_ant = False):

        dat_copy = np.copy(dat)

        if apply_bad_ant:
            dat_copy[self.bad_ant != 0] = 0

        flagged_events = np.full((self.num_evts), 0, dtype = int)
        for string in range(num_ddas):
            dat_sum = np.nansum(dat_copy[string::num_ddas], axis = 0)
            flagged_events += (dat_sum > ant_limit).astype(int)
            del dat_sum
        flagged_events = (flagged_events > st_limit).astype(int)

        return flagged_events

    def get_post_qual_cut_value(self):

        return self.unlock_cal_evts

    def get_post_qual_cut(self):

        tot_post_qual_cut = np.full((self.num_evts, 1), 0, dtype = int)
        tot_post_qual_cut[:, 0] = self.unlock_cal_evts

        self.post_qual_cut_sum = np.nansum(tot_post_qual_cut, axis = 1)

        if self.verbose:
            quick_qual_check(tot_post_qual_cut[:, 0] != 0, self.evt_num, 'unlocked calpulser events!')
            quick_qual_check(self.post_qual_cut_sum != 0, self.evt_num, 'total post qual cut!')
        
        return tot_post_qual_cut

class ped_qual_cut_loader:

    def __init__(self, ara_uproot, total_qual_cut, daq_cut_sum, analyze_blind_dat = False, verbose = False):
    
        self.analyze_blind_dat = analyze_blind_dat
        self.verbose = verbose
        self.ara_uproot = ara_uproot
        self.trig_type = self.ara_uproot.get_trig_type()
        self.num_evts = self.ara_uproot.num_evts
        self.evt_num = ara_uproot.evt_num
        self.st = self.ara_uproot.station_id
        self.run = self.ara_uproot.run
        self.total_qual_cut = total_qual_cut
        self.daq_cut_sum = daq_cut_sum
        self.num_qual_type = 3

    def get_clean_events(self):

        clean_evts_qual_type = np.full((self.total_qual_cut.shape[1], self.num_qual_type), 0, dtype = int)
        clean_evts = np.full((self.num_evts, self.num_qual_type), 0, dtype = int)

        # 0~4 daq error
        # 5 single block
        # 6~8 readout window
        # 9 bad unix sequence
        # 10 bad unix time
        # 11 first minute
        # 12 dda voltage
        # 13 bad cal ratio
        # 14 bad rf min rate
        # 15 bad cal min rate
        # 16 bad soft min rate
        # 17 bad rf sec rate
        # 18 bad cal sec rate
        # 19 bad soft sec rate
        # 20 unlock calpulser

        # turn on all cuts
        clean_evts_qual_type[:, 0] = 1
        clean_evts[:, 0] = np.logical_and(np.nansum(self.total_qual_cut, axis = 1) == 0, self.trig_type != 1).astype(int)

        # hardware error only. not use 1) 10 bad unix time, 3) 13 bad cal ratio, and 4) 14 bad rf rate
        qual_type = np.array([0,1,2,3,4,5,6,7,8,9,11,12,15,16,17,18,19,20], dtype = int)
        clean_evts_qual_type[qual_type, 1] = 1
        clean_evts[:, 1] = np.logical_and(np.nansum(self.total_qual_cut[:, qual_type], axis = 1) == 0, self.trig_type != 1).astype(int)
        del qual_type

        # only rf/software
        qual_type = np.array([0,1,2,3,4,5], dtype = int)
        clean_evts_qual_type[qual_type, 2] = 1
        clean_evts[:, 2] = np.logical_and(np.nansum(self.total_qual_cut[:, qual_type], axis = 1) == 0, self.trig_type != 1).astype(int)
        del qual_type
    
        # clean evts for repeder
        clean_num_evts = np.nansum(clean_evts, axis = 0)
        print(f'total uesful events for ped: {clean_num_evts}')

        return clean_evts, clean_evts_qual_type, clean_num_evts

    def get_block_usage(self, clean_evts):

        # ped counter
        block_usage = np.full((num_blks, self.num_qual_type), 0, dtype = int)
        for evt in range(self.num_evts):

            if self.daq_cut_sum[evt] != 0:
                continue

            blk_idx_arr = self.ara_uproot.get_block_idx(evt, trim_1st_blk = True)[0]
            if clean_evts[evt, 0] == 1:
                block_usage[blk_idx_arr, 0] += 1
            if clean_evts[evt, 1] == 1:
                block_usage[blk_idx_arr, 1] += 1
            if clean_evts[evt, 2] == 1:
                block_usage[blk_idx_arr, 2] += 1
            del blk_idx_arr

        low_block_usage = np.any(block_usage < 2, axis = 0).astype(int)
        print(f'low_block_usage flag: {low_block_usage}')

        return block_usage, low_block_usage

    def get_pedestal_qualities(self, clean_evts, block_usage, low_block_usage):

        # select final type
        final_type = np.full((1), len(low_block_usage) - 1, dtype = int)
        ped_counts = np.copy(block_usage[:, -1])
        ped_qualities = np.copy(clean_evts[:, -1])

        for t in range(self.num_qual_type):
            if low_block_usage[t] == 0:
                ped_counts = np.copy(block_usage[:, t])
                ped_qualities = np.copy(clean_evts[:, t])
                final_type[:] = t
                break
        print(f'type {final_type} was chosen for ped!')

        Output = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{self.st}/ped_full/'
        if not os.path.exists(Output):
            os.makedirs(Output) 

        txt_file_name = f'{Output}ped_full_qualities_A{self.st}_R{self.run}.dat'
        np.savetxt(txt_file_name, ped_qualities.astype(int), fmt='%i')
        print(f'output is {txt_file_name}')

        return ped_qualities, ped_counts, final_type

    def get_pedestal_information(self):

        if self.analyze_blind_dat:
            clean_evts, clean_evts_qual_type, clean_num_evts = self.get_clean_events()
            block_usage, low_block_usage = self.get_block_usage(clean_evts)
            ped_qualities, ped_counts, final_type = self.get_pedestal_qualities(clean_evts, block_usage, low_block_usage)
            self.ped_counts = ped_counts
        else:
            clean_evts = np.full((self.num_evts, self.num_qual_type), np.nan, dtype = float)
            clean_evts_qual_type = np.full((self.total_qual_cut.shape[1], self.num_qual_type), np.nan, dtype = float)
            clean_num_evts = np.full((self.num_qual_type), np.nan, dtype = float)
            block_usage = np.full((num_blks, self.num_qual_type), np.nan, dtype = float)
            low_block_usage = np.full((self.num_qual_type), np.nan, dtype = float)
            ped_qualities = np.full((self.num_evts), np.nan, dtype = float)
            ped_counts = np.full((num_blks), np.nan, dtype = float)
            final_type = np.full((1), np.nan, dtype = float)         

        return clean_evts, clean_evts_qual_type, clean_num_evts, block_usage, low_block_usage, ped_qualities, ped_counts, final_type

    def run_ped_qual_cut(self):

        if self.analyze_blind_dat:
            ped_counts = np.copy(self.ped_counts)    
        else:
            run_info = run_info_loader(self.st, self.run, analyze_blind_dat = self.analyze_blind_dat)
            ped_count_dat = run_info.get_result_path(file_type = 'qual_cut', verbose = self.verbose, force_blind = True)
            ped_count_hf = h5py.File(ped_count_dat, 'r')
            ped_counts = ped_count_hf['ped_counts'][:]
            del ped_count_dat, ped_count_hf, run_info
        zero_ped_counts = ped_counts < 1
        ped_blk_counts = ped_counts == 1
        del ped_counts

        ped_qual_cut = np.full((self.num_evts, 2), 0, dtype = int)
        for evt in range(self.num_evts):

            if self.daq_cut_sum[evt] != 0:
                continue

            blk_idx_arr = self.ara_uproot.get_block_idx(evt, trim_1st_blk = True)[0]
            ped_qual_cut[evt, 0] = np.nansum(zero_ped_counts[blk_idx_arr])

            if self.trig_type[evt] == 1:
                continue

            ped_qual_cut[evt, 1] = np.nansum(ped_blk_counts[blk_idx_arr])
            del blk_idx_arr
        del ped_blk_counts, zero_ped_counts

        self.ped_qual_cut_sum = np.nansum(ped_qual_cut, axis = 1)

        if self.verbose:
            quick_qual_check(ped_qual_cut[:, 0] != 0, self.evt_num, f'zero pedestal events')
            quick_qual_check(ped_qual_cut[:, 1] != 0, self.evt_num, f'pedestal block events')
            quick_qual_check(self.ped_qual_cut_sum != 0, self.evt_num, 'total pedestal qual cut!')

        return ped_qual_cut

def get_bad_run(st, run, qual_cut_sum, ped_cut_sum):

    # bad run
    sum_flag = np.all(qual_cut_sum != 0)
    ped_flag = np.any(ped_cut_sum != 0)
    bad_run = np.array([0, 0], dtype = int)

    if sum_flag or ped_flag:
        bad_run[0] = int(sum_flag)
        bad_run[1] = int(ped_flag)
        print(f'A{st} R{run} is bad!!! Bad type:{bad_run}')
        bad_path = f'/home/mkim/analysis/MF_filters/data/qual_runs/qual_run_A{st}.txt'
        bad_run_info = f'{run} {bad_run[0]} {bad_run[1]}\n'
        if os.path.exists(bad_path):
            print(f'There is {bad_path}')
            bad_run_arr = []
            with open(bad_path, 'r') as f:
                for lines in f:
                    run_num = int(lines.split()[0])
                    bad_run_arr.append(run_num)
            bad_run_arr = np.asarray(bad_run_arr, dtype = int)
            if run in bad_run_arr:
                print(f'Run{run} is already in {bad_path}!')
            else:
                print(f'Add run{run} in {bad_path}!')
                with open(bad_path, 'a') as f:
                    f.write(bad_run_info)
            del bad_run_arr
        else:
            print(f'There is NO {bad_path}')
            print(f'Add run{run} in {bad_path}!')
            with open(bad_path, 'w') as f:
                f.write(bad_run_info)
        del bad_path, bad_run_info
    del sum_flag, ped_flag

    return bad_run

class qual_cut_loader:

    def __init__(self, analyze_blind_dat = False, verbose = False):

        self.analyze_blind_dat = analyze_blind_dat
        self.verbose = verbose

    def load_qual_cut_result(self, st, run):

        if self.analyze_blind_dat:
            d_key = 'qual_cut_full'
        else:
            d_key = 'qual_cut'

        d_path = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{st}/'
        d_path += f'{d_key}/'
        d_path += f'{d_key}_A{st}_R{run}.h5'
        qual_file = h5py.File(d_path, 'r')
        if self.verbose:
            print(f'quality cut path:', d_path)

        self.evt_num = qual_file['evt_num'][:]
        self.unix_time = qual_file['unix_time'][:]
        total_qual_cut = qual_file['total_qual_cut'][:]
        self.daq_qual_cut_sum = qual_file['daq_qual_cut_sum'][:]
        self.total_qual_cut_sum = qual_file['total_qual_cut_sum'][:]

        if self.verbose:
            quick_qual_check(self.daq_qual_cut_sum != 0, self.evt_num, 'daq error cut!')
            quick_qual_check(self.total_qual_cut_sum != 0, self.evt_num, 'total qual cut!')
        del d_key, d_path, qual_file

        return total_qual_cut

    """
    def get_qual_cut_class(self, ara_root, ara_uproot, dt = 0.5):

        self.pre_qual = pre_qual_cut_loader(ara_uproot, analyze_blind_dat = self.analyze_blind_dat, verbose = self.verbose)
        self.post_qual = post_qual_cut_loader(ara_uproot, ara_root, dt = dt)

    def get_qual_cut_result(self):

        pre_qual_cut = self.pre_qual.run_pre_qual_cut()
        post_qual_cut = self.post_qual.run_post_qual_cut()
        total_qual_cut = np.append(pre_qual_cut, post_qual_cut, axis = 1)
        del pre_qual_cut, post_qual_cut

        if self.verbose:
            quick_qual_check(np.nansum(total_qual_cut, axis = 1) != 0, self.pre_qual.evt_num, 'total qual cut!')

        return total_qual_cut
    """


















