##
# @file e1_effective_area_calculator.py
#
# @section Created on 09/14/2022, mkim@icecube.wisc.edu
#
# @brief This is designed to calculate effective area (trigger level) from AraSim output

import os
import numpy as np
import h5py
import click # 'pip3 install click' will make you very happy
import uproot # 'pip3 install uproot' will also make you very happy
from glob import glob
from tqdm import tqdm

@click.command()
@click.option('-d', '--data_path', type = str, help = 'ex) /data/user/mkim/OMF_filter/ARA02/sim_signal/')
@click.option('-o', '--output_path', type = str, help = 'ex) /home/mkim/')
def main(data_path, output_path):
    """! This is designed to calculate effective area (trigger level) from AraSim output
        This script will only work if ARASim settings are 
        EXPONENT=1, E^-1 spectrum
        ONLY_PASSED_EVENTS=1, AraSim throws events until the number of events that pass the trigger is equal to NNU_PASSED
        INTERACTION_MODE=0, Aeff mode. spherical volume

    @param data_path  string. data path that contains multiple sim outputs
    @param output_path. string. output path
    """
   
    ## check the number of sim output
    tot_file_path = glob(f'{data_path}*')
    print(f'total number of sim outputs: {len(tot_file_path)}')

    ## load output and save information
    pnu = np.full((len(tot_file_path) * 2), np.nan, dtype = float) # Neutrino energy
    cos_angle = np.copy(pnu) # Neutrino zenith angle
    probability = np.copy(pnu) # interaction probability 
    inu_thrown = 0 # total thrown event

    for run in tqdm(range(len(tot_file_path))):
      #if run < 10: # for debug 

        file = uproot.open(tot_file_path[run])
        if run == 0:
            radius = np.asarray(file['AraTree/settings/POSNU_RADIUS'], dtype = int)[0] # radius of spherical volume [m]
        ara_tree_2 = file['AraTree2/event']
        pnu[run*2:run*2+2] = np.asarray(ara_tree_2['pnu'], dtype = float)
        cos_angle[run*2:run*2+2] = np.asarray(ara_tree_2['Nu_Interaction/Nu_Interaction.nnu.theta'], dtype = float)
        probability[run*2:run*2+2] = np.asarray(ara_tree_2['Nu_Interaction/Nu_Interaction.probability'], dtype = float) 
        inu_thrown += np.asarray(ara_tree_2['inu_thrown'], dtype = int)[-1]
        del file, ara_tree_2
    del tot_file_path

    pnu /= 1e9 # eV to GeV
    cos_angle = np.cos(cos_angle) # degree to radian

    ## IceCube oneweight calculation for E-1 spectrum
    ## related document: https://www.dropbox.com/s/yfl9z55zy4p6njl/effective_area_V1.pdf?dl=0
    solid_angle = 4 * np.pi
    area = np.pi * (radius**2)    
    log_emax = 1e12 # 12 GeV
    log_emin = 1e7 # 7 Gev
    one_weight = probability * pnu * (np.log(log_emax) - np.log(log_emin)) * solid_angle * area

    ## bin space for effective area
    energy_bins = np.logspace(np.log10(log_emin), np.log10(log_emax), 40 + 1)
    cos_bins = np.linspace(-1, 1, 100 + 1)

    ## effective area as function of energy [m^2]
    aeff_1d = np.histogram(pnu, weights = one_weight, bins = energy_bins)[0] 
    aeff_1d /= inu_thrown * np.diff(energy_bins) * solid_angle

    ## effective area as function of energy and cosine angle [m^2]
    aeff_2d = np.histogram2d(pnu, cos_angle, weights = one_weight, bins=(energy_bins, cos_bins))[0]
    aeff_2d /= inu_thrown * np.diff(energy_bins)[:, np.newaxis] * np.diff(cos_bins)[np.newaxis, :] * solid_angle

    ## create output path
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    hf_file_name = f'{output_path}'
    hf = h5py.File(hf_file_name, 'w')
    hf.create_dataset('energy_bins', data=energy_bins, compression="gzip", compression_opts=9)
    hf.create_dataset('cos_bins', data=cos_bins, compression="gzip", compression_opts=9)
    hf.create_dataset('aeff_1d', data=aeff_1d, compression="gzip", compression_opts=9)
    hf.create_dataset('aeff_2d', data=aeff_2d, compression="gzip", compression_opts=9)
    hf.close()

    print(f'Output is {hf_file_name}')
    print('Done!')

if __name__ == "__main__":

    main()















