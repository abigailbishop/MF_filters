import os
import numpy as np
import h5py
from scipy.signal import butter, filtfilt, hilbert, fftconvolve, correlation_lags
from scipy.ndimage import maximum_filter1d
from tqdm import tqdm

# custom lib
from tools.ara_constant import ara_const
from tools.ara_run_manager import run_info_loader
from tools.ara_known_issue import known_issue_loader

ara_const = ara_const()
num_ants = ara_const.USEFUL_CHAN_PER_STATION
num_pols = ara_const.POLARIZATION

class ara_matched_filter:

    def __init__(self, st, run, dt, pad_len, roll_win = 120, get_sub_file = False, use_all_chs = False, use_debug = False, verbose = False):

        self.st = st
        self.run = run
        self.dt = dt
        self.pad_len = pad_len
        self.roll_win_idx = int(roll_win / self.dt) + 1
        self.use_debug = use_debug
        self.use_all_chs = use_all_chs
        self.verbose = verbose
        self.run_info = run_info_loader(self.st, self.run, analyze_blind_dat = True)

        if get_sub_file:
            self.get_zero_pad()
            self.lags = correlation_lags(self.double_pad_len, self.double_pad_len, 'same') * self.dt
            self.lag_len = len(self.lags)
            self.get_psd()
            self.get_template()
            self.normaization()
            if self.use_debug == False:
                del self.temp_rfft
            
            known_issue = known_issue_loader(self.st)
            self.good_chs = known_issue.get_bad_antenna(self.run, good_ant_true = True, print_ant_idx = True)
            self.good_v_len = np.count_nonzero(self.good_chs < 8)
            if self.use_all_chs == False:
                self.temp = self.temp[:, self.good_chs]
                self.zero_pad = self.zero_pad[:, self.good_chs]
            del known_issue
            
            if self.verbose:
                print('sub tools are ready!')
        else:
            self.good_chs = np.arange(num_ants, dtype = int)
            self.lags = correlation_lags(self.pad_len, self.pad_len, 'full') * self.dt
            self.lag_len = len(self.lags)

        if self.verbose:
            print('useful antenna chs for mf:', self.good_chs)

    def get_zero_pad(self):

        self.double_pad_len = self.pad_len * 2
        self.zero_pad = np.full((self.double_pad_len, num_ants), 0, dtype = float)
        self.quater_idx = self.pad_len // 2

        if self.verbose:
            print('pad is on!')

    def get_psd(self):

        ## get psd from rayl sigma
        rayl_dat = self.run_info.get_result_path(file_type = 'rayl', verbose = self.verbose)
        rayl_hf = h5py.File(rayl_dat, 'r')
        soft_rayl = np.nansum(rayl_hf['soft_rayl'][:], axis = 0)
        if self.use_debug:
            self.soft_rayl = np.copy(soft_rayl)
        del rayl_dat, rayl_hf

        ## sigma to mean
        sigma_to_mean = np.sqrt(np.pi / 2)
        self.psd = (soft_rayl * sigma_to_mean) ** 2
        del soft_rayl, sigma_to_mean

        if self.verbose:
            print('psd is on!')

    def get_normalization(self):

        norm_fac = np.abs(self.temp_rfft)**2 / self.psd[:, :, np.newaxis]
        norm_fac = np.sqrt(np.nansum(norm_fac, axis = 0))
        if self.use_debug:
            self.norm_fac = np.copy(norm_fac)

        self.temp = self.temp[::-1] / norm_fac[np.newaxis, :, :]
        del norm_fac

        if self.verbose:
            print('norm is on!')

    def get_template(self, use_sc = False):

        config = self.run_info.get_config_number()
        temp_dat = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{self.st}/temp_sim/temp_AraOut.A{self.st}_R{config}.txt.run0.h5'
        if self.verbose:
            print('template:', temp_dat)
        temp_hf = h5py.File(temp_dat, 'r')
        self.rec_ang = temp_hf['rec_ang'][:]
        self.launch_ang = temp_hf['launch_ang'][:]
        self.view_ang = temp_hf['view_ang'][:]
        self.arr_time = temp_hf['arrival_time'][:]
        diff_idx = -np.round((self.arr_time - np.nanmean(self.arr_time, axis = 0)[:, np.newaxis]) / self.dt).astype(int)

        temp_ori = temp_hf['temp'][:]
        self.num_temps = temp_ori.shape[2]
        self.temp = np.full(temp_ori.shape, 0, dtype = float)
        for t in range(self.num_temps):
            for a in range(num_ants):
                if diff_idx[a, t] > 0:
                    self.temp[diff_idx[a, t]:, :, t] = temp_ori[:-diff_idx[a, t], :, t]
                elif diff_idx[a, t] < 0:
                    self.temp[:diff_idx[a, t], :, t] = temp_ori[-diff_idx[a, t]:, :, t]
                else:
                    pass
        self.temp_rfft = temp_hf['temp_rfft'][:]

        self.corr_sum_fla_shape = (num_pols, self.lag_len * self.num_temps)
        if self.use_debug:
            self.corr_sum_shape = (num_pols, self.lag_len, self.num_temps)
            self.temp_ori = np.copy(temp_ori)
            self.temp_roll = np.copy(self.temp)
        del config, temp_dat, temp_hf, diff_idx, temp_ori

        if self.verbose:
            print('template dim:', self.temp.shape)
            print('template is on!')

    def get_padded_wf(self, pad_v):

        self.zero_pad[:] = 0
        self.zero_pad[self.quater_idx:-self.quater_idx] = pad_v

    def get_mf_wfs(self):

        # fft correlation w/ multiple array at once
        self.corr = fftconvolve(self.temp, self.zero_pad[:, :, np.newaxis], 'same', axes = 0)
        if self.use_debug:
            self.corr_no_hill = np.copy(self.corr)

        self.corr = np.abs(hilbert(self.corr, axis = 0))
        if self.use_debug:
            self.corr_hill = np.copy(self.corr)

    def get_evt_wise_corr(self, corr):

        ## smoothing corr
        corr_roll_max = maximum_filter1d(self.corr, axis = 0, size = self.roll_win_idx, mode='constant')
        if self.use_debug:
            self.corr_roll_max = np.copy(corr_roll_max)

        corr_sum = np.full(self.corr_sum_shape, np.nan, dtype = float)
        corr_sum[0] = np.nansum(corr_roll_max[:, :good_v_len], axis = 1).flatten()
        corr_sum[1] = np.nansum(corr_roll_max[:, good_v_len:], axis = 1).flatten()
        del corr_roll_max
        if self.use_debug:
            self.corr_sum = np.reshape(corr_sum, self.corr_sum_shape)

        self.mf_max = np.nanmax(corr_sum, axis = 1)
        self.mf_temp = np.nanargmax(corr_sum, axis = 1) % self.num_temps
        del corr_sum

    def get_evt_wise_snr(self, pad_v, weights = None):

        ## zero pad
        self.get_padded_wf(pad_v[:, self.good_chs])

        ## matched filter
        self.get_mf_wfs()

        if weights is not None:
           self.corr *= weights

        ## event wise snr
        self.get_evt_wise_corr()
        if self.use_debug == False:
            del self.corr

def get_products(weights, good_chs, good_v_len):

    weights = weights[good_chs]
    wei_v_sum = np.nansum(weights[:good_v_len], axis = 0)
    wei_h_sum = np.nansum(weights[good_v_len:], axis = 0)
    weights[:good_v_len] /= wei_v_sum[np.newaxis, :]
    weights[good_v_len:] /= wei_h_sum[np.newaxis, :]
    del wei_v_sum, wei_h_sum

    return weights












