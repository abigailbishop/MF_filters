import os
import numpy as np
from scipy.signal import fftconvolve
from scipy.signal import correlation_lags
from scipy.signal import hilbert
from tqdm import tqdm
import h5py
from itertools import combinations

# custom lib
from tools.ara_constant import ara_const

ara_const = ara_const()
num_ants = ara_const.USEFUL_CHAN_PER_STATION

class py_interferometers:

    def __init__(self, pad_len, dt, st, yrs, run = None, get_sub_file = False):

        self.dt = dt
        self.st = st
        self.yrs = yrs
        self.run = run
        self.lags = correlation_lags(pad_len, pad_len, 'same') * self.dt
        self.lag_len = len(self.lags)

        if get_sub_file:
            self.get_pair_info()
            self.get_arrival_time_tables()
            self.get_coval_time()
            self.pad_one = np.full((pad_len, num_ants), 1, dtype = float)

    def get_pair_info(self):

        if self.run is not None:
            from tools.ara_known_issue import known_issue_loader
            known_issue = known_issue_loader(self.st)
            good_ant = known_issue.get_bad_antenna(self.run, good_ant_true = True, print_ant_idx = True) 
            del known_issue
        else:
            good_ant = np.arange(num_ants, dtype = int)
        print('useful antenna chs for reco:', good_ant)

        v_pairs = np.asarray(list(combinations(good_ant[good_ant < 8], 2)))
        h_pairs = np.asarray(list(combinations(good_ant[good_ant > 7], 2)))
        self.pairs = np.append(v_pairs, h_pairs, axis = 0)
        self.pair_len = len(self.pairs)
        self.v_pairs_len = len(v_pairs)
        del v_pairs, h_pairs, good_ant
        print('number of pairs:', len(self.pairs))

    def get_arrival_time_tables(self):

        if self.st == 2 or (self.st == 3 and self.yrs <= 1515974400):
            year = 2015
        else:
            year = 2018

        table_path = os.path.expandvars("$OUTPUT_PATH") + f'/OMF_filter/ARA0{self.st}/arr_time_table/'
        table_name = f'arr_time_table_A{self.st}_Y{year}.h5'
        print('arrival time table:', table_path + table_name)        
        del year

        table_hf = h5py.File(table_path + table_name, 'r')
        theta = table_hf['theta_bin'][:] - 90 # zenith to elevation angle
        phi = table_hf['phi_bin'][:]
        radius_arr = table_hf['radius_bin'][:]
        num_ray_sol = table_hf['num_ray_sol'][0]
        arr_table = table_hf['arr_time_table'][:]
        del table_path, table_name, table_hf
 
        self.table = np.full((len(theta), len(phi), len(radius_arr), num_ray_sol, self.pair_len), np.nan, dtype = float)
        table_p1 = np.copy(self.table)
        table_p2 = np.copy(self.table)
        for p in range(self.pair_len):
            p_1st = self.pairs[p, 0]
            p_2nd = self.pairs[p, 1]
            self.table[:, :, :, :, p] = arr_table[:, :, :, p_1st, :] - arr_table[:, :, :, p_2nd, :]
            table_p1[:, :, :, :, p] = arr_table[:, :, :, p_1st, :]
            table_p2[:, :, :, :, p] = arr_table[:, :, :, p_2nd, :]
            del p_1st, p_2nd
        del theta, phi, radius_arr, num_ray_sol, arr_table

        self.bad_arr = np.logical_or(table_p1 < -100, table_p2 < -100)
        self.table_shape = self.table.shape
        print('arr table shape:', self.table_shape)
        del table_p1, table_p2

    def get_coval_time(self):

        self.p0_idx = np.floor((self.table - self.lags[0]) / self.dt).astype(int)
        self.p0_idx[self.p0_idx < 0] = 0
        self.p0_idx[self.p0_idx >= self.lag_len - 1] = self.lag_len - 2

        self.int_factor = (self.table - self.lags[self.p0_idx])/self.dt

    def get_coval_sample(self, corr, sum_pol = False):

        corr_diff = corr[1:] - corr[:-1]

        coval = np.full(self.table_shape, 0, dtype=float)
        for p in range(self.pair_len):
            coval[:, :, :, :, p] = corr_diff[:, p][self.p0_idx[:, :, :, :, p]] * self.int_factor[:, :, :, :, p] + corr[:, p][self.p0_idx[:, :, :, :, p]]
        coval[self.bad_arr] = 0
        del corr_diff

        if sum_pol:
            corr_v_sum = np.nansum(coval[:, :, :, :, :self.v_pairs_len], axis = 4)
            corr_h_sum = np.nansum(coval[:, :, :, :, self.v_pairs_len:], axis = 4)

            corr_v_max = np.nanmax(corr_v_sum, axis = (0,1))
            corr_h_max = np.nanmax(corr_h_sum, axis = (0,1))
            coval = np.asarray([corr_v_max, corr_h_max]) # pol, rad, sol

            coord = np.full((2, 2, 2, 2), np.nan, dtype = float) # pol, thetapi, rad, sol
            for r in range(2):
                for s in range(2):
                    coord[0, :, r, s] = np.asarray(np.where(corr_v_sum[:, :, r, s] == corr_v_max[r, s]), dtype = float).flatten()                 
                    coord[1, :, r, s] = np.asarray(np.where(corr_h_sum[:, :, r, s] == corr_h_max[r, s]), dtype = float).flatten()
            del corr_v_sum, corr_h_sum

            return coval, coord
        else:
            return coval

    def get_cross_correlation(self, pad_v, return_debug_dat = False):

        # fft correlation w/ multiple array at once
        corr = fftconvolve(pad_v[:, self.pairs[:, 0]], pad_v[::-1, self.pairs[:, 1]], 'same', axes = 0)
        if return_debug_dat:
            corr_nonorm = np.copy(corr)

        # normalization factor by wf weight
        nor_fac = fftconvolve(pad_v**2, self.pad_one, 'same', axes = 0)
        nor_fac = np.sqrt(nor_fac[::-1, self.pairs[:, 0]] * nor_fac[:, self.pairs[:, 1]])
        corr /= nor_fac
        corr[np.isnan(corr) | np.isinf(corr)] = 0 # convert x/nan result

        # hilbert
        corr = np.abs(hilbert(corr, axis = 0))

        if return_debug_dat:
            return corr, corr_nonorm, nor_fac
        else:
            del nor_fac
            return corr
    
    def get_sky_map(self, pad_v, weights = None):
        
        # correlation
        corr = self.get_cross_correlation(pad_v)

        if weights is not None:
           corr *= weights

        #coval
        coval, coord = self.get_coval_sample(corr, sum_pol = True)
        del corr

        return coval, coord

def get_products(weights, pairs, v_pairs_len):
   
    wei_pairs = weights[pairs[:, 0]] * weights[pairs[:, 1]]
    wei_v_sum = np.nansum(wei_pairs[:v_pairs_len], axis = 0)
    wei_h_sum = np.nansum(wei_pairs[v_pairs_len:], axis = 0)
    wei_pairs[:v_pairs_len] /= wei_v_sum[np.newaxis, :]
    wei_pairs[v_pairs_len:] /= wei_h_sum[np.newaxis, :]
    del wei_v_sum, wei_h_sum 

    return wei_pairs





















