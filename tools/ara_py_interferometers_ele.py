import os
import numpy as np
from scipy.signal import fftconvolve
from scipy.signal import correlation_lags
from scipy.signal import hilbert
from scipy.ndimage import uniform_filter1d
from tqdm import tqdm
import h5py
from itertools import combinations

# custom lib
from tools.ara_constant import ara_const
from tools.ara_run_manager import get_pair_info

ara_const = ara_const()
num_ants = ara_const.USEFUL_CHAN_PER_STATION
num_pols = ara_const.POLARIZATION

class py_interferometers:

    def __init__(
        self, pad_len, dt, st, 
        run = None, get_sub_file = False, use_debug = False, 
        verbose = False, use_only_max = False
    ):
        """
        pad_len

        Parameters
        ----------
        pad_len : int
            Waveform padding length to use when getting the number of 
              lag/displacement indices for 1D cross correlation.
        dt : float
            Timestep.
        st : int
            Station number.
        run : int
            Run Number
        get_sub_file : bool
            If True, ??
        use_debug : bool
            If True, will printout debugging information.
        verbose : bool
            If True, will printout all calculated information.
        use_only_max : bool
            If True, will only save data for ??
            Appears to be used with the 100% dataset to save space.
        """

        self.verbose = verbose
        self.dt = dt
        self.st = st
        self.run = run
        self.use_debug = use_debug
        self.use_only_max = use_only_max
        if self.verbose and self.use_only_max:
            print('ONLY MAX VALUE!')

        if get_sub_file:
            self.get_zero_pad(pad_len)

            # Get the number of lag/displacement indices for 1D cross correlation
            self.lags = correlation_lags(
                self.double_pad_len, self.double_pad_len, 'same') * self.dt
            self.lag_len = len(self.lags)

            # Collect every pair combination of good antennas
            self.pairs, self.pair_len, self.v_pairs_len = get_pair_info(
                self.st, self.run, verbose = self.verbose)
            self.pair_range = np.arange(self.pair_len, dtype = int)

            self.get_arrival_time_tables()
            self.get_coval_time()
            if self.verbose:
                print('sub tools are ready!')
        else:      
            # Get the number of lag/displacement indices for 1D cross correlation          
            self.lags = correlation_lags(pad_len, pad_len, 'full') * self.dt    
            self.lag_len = len(self.lags)

    def get_zero_pad(self, pad_len):
        """
        Creates a 2D array of zeros, `self.zero_pad`, containing padding that is 
          double the provided padding length for each antenna.

        Parameters
        ----------
        pad_len = int
            Length of zero padding array to create.
        """
    
        self.double_pad_len = pad_len * 2 
        self.double_pad_len_float = float(self.double_pad_len)
        self.zero_pad = np.full((self.double_pad_len, num_ants), 0, dtype = float)
        self.quater_idx = pad_len // 2
        if self.verbose:
            print('pad is on!')

    def get_arrival_time_tables(self):
        """
        Build a table of time delays between every vpol+vpol pair and
          every hpol+hpol pair for all thetas, phis, radii, and rays
        Build table of booleans reporting if arrival time table was calculated
          properly for all radii in addition to all 
          vertically polarized signal for vpol pairs and
          horizontally polarized signal for hpol pairs
          for a given set of thetas and phis. 

        If use_only_max is `True`, will build a table of booleans 
          reporting if events likely came from the suface based on their
          zenith angles for different radii. 

        Attributes
        ----------
        self.table : numpy.ndarray
            Table of time delays for all thetas, phis, radii, rays, and 
              antenna pairs. 
        self.arr_table : numpy.ndarray
            If use_debug = `True`, saves the arrival time table to this 
              class object.
        """

        # Load the arrival time table
        table_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{self.st}/arr_time_table/'
        table_name = f'arr_time_table_A{self.st}_all.h5'
        if self.verbose:
            print('arrival time table:', table_path + table_name)        
        table_hf = h5py.File(table_path + table_name, 'r')

        # Import theta, phi, radii, rays to this class object
        self.theta = 90 - table_hf['theta_bin'][:] # nadir to elevation angle
        self.num_thetas = len(self.theta)
        self.phi = table_hf['phi_bin'][:]
        self.num_phis = len(self.phi)
        self.radius = table_hf['radius_bin'][:]
        if self.verbose:
            print('radius param:', self.radius)
        self.num_rads = len(self.radius)
        self.num_rays = len(table_hf['num_ray_sol'][:])
        self.num_rays += 1 # it is for d and r merged results 

        # Import the flight times and reorganize it so the dimensions are
        #   theta, phi, rad, ray, ant
        arr_table = table_hf['arr_time_table'][:] # theta, phi, rad, ant, ray
        arr_table = np.transpose(arr_table, (0, 1, 2, 4, 3)) # theta, phi, rad, ray, ant 
        if self.use_debug:
            self.arr_table = np.copy(arr_table)
        del table_path, table_name, table_hf

        # Initialize that we want theta and phi calculated 
        #   for v and h polarizations
        self.num_angs = int(2) # theta and phi
        self.num_pols = num_pols # v and h

        # If using light data saving mode, create arrays for event depths and 
        #   surface event zenith identifiers.
        if self.use_only_max:

            # Calculate event depths for calculation in multidimensional array
            self.z_cal = np.sin(np.radians(self.theta)) * float(self.radius[0])

            # Total number of calculations to make: every theta, radius, and ray
            flat_len = self.num_thetas * self.num_rads * self.num_rays

            # Calculate event depths for calculation in a flattened array
            theta_ex = np.full(
                (self.num_thetas, self.num_rads, self.num_rays), 
                np.nan, dtype = float)
            theta_ex[:] = self.theta[:, np.newaxis, np.newaxis]
            self.theta_flat = np.reshape(theta_ex, (flat_len))
            rad_ex = np.full(theta_ex.shape, np.nan, dtype = float)
            rad_ex[:] = self.radius[np.newaxis, :, np.newaxis]
            self.rad_flat = np.reshape(rad_ex, (flat_len))
            self.z_flat = np.sin(np.radians(self.theta_flat)) * self.rad_flat
            del theta_ex, rad_ex

            # Create array of bools that will tell you if an event likely
            #   came from the surface for a given pol, theta, rad, and ray
            sur_ang = np.array([180, 180, 37, 24, 17], dtype = float)
            theta_map = np.full(
                (self.num_pols, self.num_thetas, self.num_rads, self.num_rays), 
                np.nan, dtype = float)
            theta_map[:] = self.theta[np.newaxis, :, np.newaxis, np.newaxis]
            sur_bool = theta_map <= sur_ang[np.newaxis, np.newaxis, :, np.newaxis]
            self.sur_bool_flat = np.reshape(sur_bool, (self.num_pols, flat_len))
            del sur_ang, theta_map, sur_bool, flat_len

        # Calculate the time difference between each pair of good antennas
        table_p1 = arr_table[:, :, :, :, self.pairs[:, 0]]
        table_p2 = arr_table[:, :, :, :, self.pairs[:, 1]] 
        self.table = table_p1 - table_p2 # theta, phi, rad, ray, pair
        del arr_table

        # Set all time differences to `nan` if arrival time is too negative 
        #   as there must've been an error in ray tracing
        minus_arr = np.logical_or(table_p1 < -100, table_p2 < -100) 
        self.table[minus_arr] = np.nan

        # For channel pairs that did not have ray solutions, set values in 
        #   table to `nan`
        partial_ray = np.any(np.isnan(self.table), axis = 4)
        partial_ray_ex = np.repeat(
            partial_ray[:, :, :, :, np.newaxis], self.pair_len, axis = 4) # match table shape
        self.table[partial_ray_ex] = np.nan

        # Cleanup
        del table_p1, table_p2, minus_arr, partial_ray, partial_ray_ex

        # Save more metadata to the class object
        self.table_shape = self.table.shape # (theta, phi, ray, rad, pair)
        self.sky_map_shape = (self.num_pols, self.num_thetas, self.num_phis, 
                              self.num_rads, self.num_rays)
        self.results_shape = (self.num_pols, self.num_thetas, 
                              self.num_rads, self.num_rays)

        # Initialize array indicating failed/bad reconstructions
        self.bad_coval = np.isnan(self.table)
        self.bad_sky = np.full(self.sky_map_shape, False, dtype = bool)

        # For vertically polarized signal along direct rays and refracted rays,
        #   save whether all arrival timetable calculations for all 
        #   vpol pairs are valid
        self.bad_sky[0, :, :, :, :2] = np.all(
            np.isnan(self.table[:, :, :, :, :self.v_pairs_len]), 
            axis = 4
        )
        
        # For horizontally polarized signal along direct and refracted rays,
        #   save whether all arrival timetable calculations for all
        #   hpol pairs are valid
        self.bad_sky[1, :, :, :, :2] = np.all(
            np.isnan(self.table[:, :, :, :, self.v_pairs_len:]), 
            axis = 4
        ) 
        
        # For vertically polarized signal in the direct+refracted solution,
        #   check that solutions at all radii are valid
        self.bad_sky[0, :, :, :, 2] = np.all(
            self.bad_sky[0, :, :, :, :2], 
            axis = 3
        ) 
        
        # For horizontally polarized signal in the direct+refracted solution,
        #   check that solutions at all radii are valid
        self.bad_sky[1, :, :, :, 2] = np.all(
            self.bad_sky[1, :, :, :, :2], 
            axis = 3
        )
        
        if self.verbose:
            print('arr table shape:', self.table_shape)

        # for advanced indexing
        self.pol_range = np.arange(self.num_pols, dtype = int)
        self.rad_range = np.arange(self.num_rads, dtype = int)
        self.ray_range = np.arange(self.num_rays, dtype = int)
        self.theta_range = np.arange(self.num_thetas, dtype = int)
        self.phi_range = np.arange(self.num_phis, dtype = int)
   
    def get_coval_time(self):
        """
        Get correlation value timing?
        """

        # Convert all `nan` values in the time difference table to a giant number
        self.p0_idx = np.floor((self.table - self.lags[0]) / self.dt).astype(int)

        # For all data where the lag is greater than the time difference, set 
        #   value to 0 (replaces some of the giant values, supposedly)
        self.p0_idx[self.p0_idx < 0] = 0 # it also replace giant vaues

        # Replaces more of the giant values, supposedly
        self.p0_idx[self.p0_idx >= self.lag_len - 1] = self.lag_len - 2

        # Save the ( time delays minus the lags ) divided by the time step
        self.int_factor = (self.table - self.lags[self.p0_idx])/self.dt

        if self.verbose:
            print('coval time is on!')

    def get_coval_sample(self):
        """
        Get correlation value for the given sample.
        """

        ## coval (theta, phi, rad, ray, pair)
        # Get the correlation values for each pair of antennas
        coval = ( np.diff(self.corr, axis = 0)[self.p0_idx, self.pair_range] 
                  * self.int_factor + self.corr[self.p0_idx, self.pair_range]  )

        # Make sure coval array is set to `nan` where
        #   time difference table has `nan` values
        coval[self.bad_coval] = np.nan 
        if self.use_debug:
            self.coval = np.copy(coval) # individual pairs sky map

        # Initialize sky map with shape (n_pol, n_theta, n_phi, n_r, n_rays)
        sky_map = np.full(self.sky_map_shape, np.nan, dtype = float) 

        # Add all non-nan correlation values for all vpol pairs 
        #   for vertically polarized signal in direct rays and refracted rays.
        sky_map[0, :, :, :, :2] = np.nansum( 
            coval[:, :, :, :, :self.v_pairs_len], 
            axis = 4
        )

        # Add all non-nan correlation values for all hpol pairs
        #   for horizontally polarized signal in direct rays and refracted rays.
        sky_map[1, :, :, :, :2] = np.nansum(
            coval[:, :, :, :, self.v_pairs_len:], 
            axis = 4
        )

        # For all locations that were identified to be "bad", set the sum
        #   to `nan` (apparently sets sum of `nan`s to zero by default?)
        sky_map[self.bad_sky] = np.nan 

        # Average the sum of correlation values from direct and refracted rays 
        #   over the radii and save this as the Direct+Refracted correlation
        sky_map[0, :, :, :, 2] = np.nanmean(
            sky_map[0, :, :, :, :2], 
            axis = 3
        ) # Vertically polarized signal
        sky_map[1, :, :, :, 2] = np.nanmean(
            sky_map[1, :, :, :, :2], 
            axis = 3
        ) # Horizontally polarized signal

        if self.use_debug:
            self.sky_map = np.copy(sky_map)
        del coval
  
        # Switch all `nan` values we set in the sky map to `-1`
        sky_map[np.isnan(sky_map)] = -1

        # For each polarization, theta, radius, and ray:
        #   Get the index of phi that has the highest correlation coefficient
        # Shape: (# of pols, # of thetas, # of rs, # of rays)
        coef_phi_max_idx = np.nanargmax(sky_map, axis = 2) 

        # For each polarization, theta, radius, and ray:
        #   Get the value for phi with the greatest correlation value
        # Shape: (# of pols, # of thetas, # of rs, # of rays)
        self.coord_max_ele = self.phi[coef_phi_max_idx] 

        # For each polarization, theta, radius, and ray:
        #   Save the maximum correlation coefficient value corresponding 
        #     to the azimuthal angle phi that we saved.
        # Skymap shape: (n_pol, n_theta, n_phi, n_r, n_rays)
        # coef_phi_max_idx shape: (n_pol, n_theta, n_r, n_rays)
        # Shape: (# of pols, # of thetas, # of rs, # of rays)
        self.coef_max_ele = sky_map[
            self.pol_range[:, np.newaxis, np.newaxis, np.newaxis], 
            self.theta_range[np.newaxis, :, np.newaxis, np.newaxis], 
            coef_phi_max_idx, 
            self.rad_range[np.newaxis, np.newaxis, :, np.newaxis], 
            self.ray_range[np.newaxis, np.newaxis, np.newaxis, :]
        ] 
        
        # If we're only analyzing the coordinates with maximum correlation
        #   (Usually the case when analyzing 100% data)
        if self.use_only_max:
            coef = self.coef_max_ele[:, :, 0, 0] # pol, theta, (rad), (ray)
            coef_idx = np.nanargmax(coef, axis = 1)
            self.coef_cal = coef[self.pol_range, coef_idx] # pol
            neg_idx = self.coef_cal < 0
            self.coef_cal[neg_idx] = np.nan
            del coef

            self.coord_cal = np.full( (self.num_angs + 1, self.num_pols), 
                                      np.nan, dtype = float) # thepiz, pol
            self.coord_cal[0] = self.theta[coef_idx]
            self.coord_cal[1] = self.coord_max_ele[:, :, 0, 0][self.pol_range, coef_idx] # pol, theta, (rad), (ray)
            self.coord_cal[2] = self.z_cal[coef_idx]
            self.coord_cal[:, neg_idx] = np.nan
            del coef_idx, neg_idx

            coef_re = np.reshape(self.coef_max_ele, (self.num_pols,  -1))
            coord_re = np.reshape(self.coord_max_ele, (self.num_pols, -1))
            self.coef_max = np.full( (self.num_pols), 
                                     np.nan, dtype = float ) # pol, evt
            self.coord_max = np.full( (self.num_angs + 2, self.num_pols), 
                                      np.nan, dtype = float ) # thepirz, pol
            self.coef_s_max = np.copy(self.coef_max)
            self.coord_s_max = np.copy(self.coord_max)
            for t in range(2):
                if t == 1:
                    coef_re[self.sur_bool_flat] = -1
                    coord_re[self.sur_bool_flat] = np.nan
                coef_max_idx = np.nanargmax(coef_re, axis = 1)
                coef_max1 = coef_re[self.pol_range, coef_max_idx] # pol
                neg_idx = coef_max1 < 0
                coef_max1[neg_idx] = np.nan
                coord_max1 = np.full( (self.num_angs + 2, self.num_pols), 
                                      np.nan, dtype = float ) # thepir, pol
                coord_max1[0] = self.theta_flat[coef_max_idx]
                coord_max1[1] = coord_re[self.pol_range, coef_max_idx]
                coord_max1[2] = self.rad_flat[coef_max_idx]
                coord_max1[3] = self.z_flat[coef_max_idx]
                coord_max1[:, neg_idx] = np.nan
                del coef_max_idx, neg_idx
                if t == 0:
                    self.coef_max[:] = coef_max1
                    self.coord_max[:] = coord_max1
                else:
                    self.coef_s_max[:] = coef_max1
                    self.coord_s_max[:] = coord_max1
                del coef_max1, coord_max1
            del coef_re, coord_re
        else:
            # Set all correlation coefficients and coordinates to `nan`
            #   if the max correlation coefficient for that event is negative
            neg_idx = self.coef_max_ele < 0
            self.coord_max_ele[neg_idx] = np.nan
            self.coef_max_ele[neg_idx] = np.nan

            # Save the indices corresponding to the phi with the max correlation
            #   coefficient for all polarizations, thetas, radii, and rays
            #   if debugging. Set the indices of all failed correlation
            #   coefficient calculations to -1.
            if self.use_debug:
                self.coef_phi_max_idx = np.copy(coef_phi_max_idx)
                self.coef_phi_max_idx[neg_idx] = -1

            del neg_idx

        del sky_map, coef_phi_max_idx

    def get_padded_wf(self):
        """
        Save a zero padded version of `self.pad_v` to `self.zero_pad`.
        """

        self.zero_pad[:] = 0
        self.zero_pad[self.quater_idx:-self.quater_idx] = self.pad_v

    def get_cross_correlation(self):
        """
        Cross correlate waveforms by hilbert enveloping the
          the FFT Convolution of two waveforms scaled by waveform weight.

        Attributes
        ----------
        self.corr : numpy.ndarray
            2D array of shape (n_samples, n_pairs) containing the hilbert 
              enveloped and weighted cross correlation for each pair
              of antennas.
        """

        # fft correlation between padded waveforms for each antenna pair
        self.corr = fftconvolve(
            self.zero_pad[:, self.pairs[:, 0]], 
            self.zero_pad[::-1, self.pairs[:, 1]], 
            'same', axes = 0
        )
        if self.use_debug:
            self.corr_nonorm = np.copy(self.corr)

        # Calculate the normalization factor by wf weight
        nor_fac = uniform_filter1d(
            self.zero_pad**2, 
            size = self.double_pad_len, mode = 'constant', axis = 0
        ) * self.double_pad_len_float
        nor_fac = np.sqrt( nor_fac[::-1, self.pairs[:, 0]] 
                           * nor_fac[:, self.pairs[:, 1]])
        
        # Scale the fft correlations by the normalization factor
        self.corr /= nor_fac
        self.corr[np.isnan(self.corr) | np.isinf(self.corr)] = 0 # convert x/nan result
        if self.use_debug:
            self.nor_fac = nor_fac
        else:
            del nor_fac

        # Calculate hilbert envelope over the cross correlated waveforms
        self.corr = np.abs(hilbert(self.corr, axis = 0))
    
    def get_sky_map(self, pad_v, weights = None):
        """
        For the provided list of waveforms in a detector in `pad_v`, pad the
          waveforms, cross correlate them, and get coval.

        Parameters
        ----------
        pad_v : numpy.ndarray
            2D array of waveforms in a station with shape (n_samples, n_ants).
        weights : number, array
            Weights to muliply against the cross correlated waveforms.
        """

        # Zero pad the provided waveforms, stored in self.zero_pad with shape
        #   (n_samples, n_pairs)
        self.pad_v = pad_v
        self.get_padded_wf()
        del self.pad_v        

        # Cross correlate the provided waveforms
        self.get_cross_correlation()

        # Incorporate weights into the cross correlation if provided
        if weights is not None:
            self.corr *= weights

        # coval
        self.get_coval_sample()
        if self.use_debug == False:
            del self.corr

def get_products(weights, pairs, v_pairs_len):
    """
    Return the product of weight divided by sum of weights for each 
      polarization. So for a vpol pair, multiply the weights of each antenna
      together then divide by the sum of all vpol weights.

    Parameters
    ----------
    weights : numpy.ndarray
        Contains the weights array for each antenna
        Shape: (n_ants, len(weights))
    pairs : numpy.ndarray
        Contains the index of each antenna in a pair.
        Shape: (n_pairs, 2 (one for each antenna in the pair) )
    v_pairs_len : int
        Number of vpol pairs. These will be the first listed in the list of 
          products of weights for each pair. Hpol pairs come after.
    """
   
    # Multiply the weights of each antenna in a pair together
    wei_pairs = weights[pairs[:, 0]] * weights[pairs[:, 1]]

    # Get the sum of weights over every vpol pair and every hpol pair
    # Shape: (n_polarizaitons, n_pairs)
    # FYI, the sum of hpol weights for a vpol pair will be set to `nan`
    #   same goes for vpol weights for hpol pairs
    wei_pol = np.full((num_pols, weights.shape[1]), np.nan, dtype = float)
    wei_pol[0] = np.nansum(wei_pairs[:v_pairs_len], axis = 0) # save vpol sum of weights
    wei_pol[1] = np.nansum(wei_pairs[v_pairs_len:], axis = 0) # save hpol sum of weights

    # Divide product of weights for each antenna pair by the sum of 
    #   weights of the same polarization.
    wei_pairs[:v_pairs_len] /= wei_pol[0][np.newaxis, :] # vpol pair
    wei_pairs[v_pairs_len:] /= wei_pol[1][np.newaxis, :] # hpol pairs
    del wei_pol

    return wei_pairs





















