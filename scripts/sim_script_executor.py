import numpy as np
import os, sys
import h5py
import click
from importlib import import_module
from tqdm import tqdm

# custom lib
curr_path = os.getcwd()
sys.path.append(curr_path+'/../')
from tools.ara_utility import size_checker
from tools.ara_run_manager import get_path_info_v2
from tools.ara_run_manager import get_file_name

@click.command()
@click.option('-k', '--key', type = str)
@click.option('-s', '--station', type = int)
@click.option('-y', '--year', default = 2015, type = int)
@click.option('-d', '--data', type = str)
@click.option('-a', '--act_evt', default = None)
@click.option('-n', '--not_override', default = False, type = bool)
def script_loader(key, station, year, data, act_evt, not_override):

    if not_override:
        output = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{station}/{key}_sim/'
        data_name = get_file_name(data)
        if key == 'cw_flag_signal':
            key_p = key[:-7]
            h5_file_name = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{station}/{key_p}_sim/{key_p}_{data_name}'
        else:
            h5_file_name = f'{output}{key}_{data_name}'
        h5_file_name_out = h5_file_name + '.h5'
        if os.path.exists(h5_file_name_out):
            print(f'{h5_file_name_out} is already there!!')
            return    

    print(data, size_checker(data))
    sim_type = get_path_info_v2(data, 'AraOut.', '_')
    config = int(get_path_info_v2(data, '_R', '.txt'))
    energy = int(get_path_info_v2(data, '_E', '_F'))
    flavor = int(get_path_info_v2(data, '_F', '_A'))
    sim_run = int(get_path_info_v2(data, 'txt.run', '.root'))
    if config < 6:
        year = 2015
    else:
        year = 2018
    print('St:', station, 'Type:', sim_type, 'Flavor:', flavor, 'Config:', config, 'Year:', year, 'Sim Run:', sim_run, 'Energy:', energy)
    del sim_type

    # run the chunk code
    module = import_module(f'tools.chunk_{key}_sim')
    method = getattr(module, f'{key}_sim_collector')
    if key == 'wf':
        results = method(data, station, year, act_evt)
    else:
        results = method(data, station, year)
    del module, method

    # create output dir
    if key == 'cw_flag_signal':
        key = key[:-7]
    output = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{station}/{key}_sim/'
    print(f'Output path check:{output}')
    if not os.path.exists(output):
        os.makedirs(output)

    data_name = get_file_name(data)
    h5_file_name = f'{output}{key}_{data_name}'
    if key == 'wf' and act_evt is not None:
        act_evt = act_evt.split(',')
        act_evt = np.asarray(act_evt).astype(int)
        if len(act_evt) == 1:
            h5_file_name += f'_E{act_evt[0]}'
        else:
            h5_file_name += f'_E{act_evt[0]}_to_E{act_evt[-1]}'
        del act_evt
    h5_file_name_out = h5_file_name + '.h5'
    hf = h5py.File(h5_file_name_out, 'w')

    #saving result
    hf.create_dataset('config', data=np.array([station, sim_run, config, year, flavor, energy]), compression="gzip", compression_opts=9)
    for r in results:
        print(r, results[r].shape)
        try:
            hf.create_dataset(r, data=results[r], compression="gzip", compression_opts=9)
        except TypeError:
            dt = h5py.vlen_dtype(np.dtype(float))
            hf.create_dataset(r, data=results[r], dtype = dt, compression="gzip", compression_opts=9)
    hf.close()
    print(f'output is {h5_file_name_out}.', size_checker(h5_file_name_out))

if __name__ == "__main__":

    script_loader()













    
