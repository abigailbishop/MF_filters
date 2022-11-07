##
# @file config_collector.py
#
# @section Created on 08/03/2022, mkim@icecube.wisc.edu
#
# @brief This is designed to collect data from configFile.run[xxxxxx].dat, runStart.run[xxxxxx].dat, and runStop.run[xxxxxx].dat

import numpy as np
import os, sys
from glob import glob
from tqdm import tqdm
from datetime import datetime
import h5py
import click

def context_finder(txt_read, key, end_key = ';', num_vals = 1):
    """! function for scrap information between two string variables. It will find the string between the variables and convert into numpy array

    @param txt_read  
    @param key  string
    @param end_key  string
    @param num_vals  integer
    @return val integer
    """
    
    ## check whether there is a same config but commented out one
    ## if there is, delete it before begin the search
    old_key = '//'+key
    new_txt_read = txt_read.replace(old_key, '')    

    key_idx = new_txt_read.find(key) # find the index of the key
    if key_idx != -1: # check whether key is in the txt_read or not
        key_idx += len(key)
        end_key_idx = new_txt_read.find(end_key, key_idx) # find the end_key index after key_idx

        if num_vals != 1: # multiple elements
            val = np.asarray(new_txt_read[key_idx:end_key_idx].split(","), dtype = int)
        else:
            val = int(new_txt_read[key_idx:end_key_idx])
    else:
        val = np.full((num_vals), np.nan, dtype = float) # it there is no key in the txt_read, output numpy nan

    return val

def scrap_run_dat_info(run_list, run_key = 'Run:', time_key = 'Time:', run_end_key = 'Message'):
    """! function for scraping the run number, and unix time from run dat file and convert into UTC date time

    @param run_list  string. dat file path
    @param run_key  string  
    @param time_key  string
    @param run_end_key  string
    @return run_num  integer
    @return unix_time integer
    @return date_time integer
    """

    ## open dat file and scrap the information
    with open(run_list,'r') as run_file:
        run_read = run_file.read()

        ## scrap the information from context
        ## check whether output of context_finder has a information or not
        ## some of files dont have a all the information...
        
        run_num = context_finder(run_read, run_key, time_key)
        if not np.isfinite(run_num):
            print(f'{run_list} doesnt has run number!! Use file name for identifying run number!')
            run_num = int(run_list[-10:-4]) # I dont like scrapping info from file name... but this is best...
        
        unix_time = context_finder(run_read, time_key, run_end_key)
        if not np.isfinite(unix_time):
            print(f'{run_list} doesnt has unix time!! Output date time as a NAN...')
            date_time = np.nan
        else:
            date_time = datetime.utcfromtimestamp(unix_time)
            date_time = date_time.strftime('%Y%m%d%H%M%S')
            date_time = int(date_time) # save date time as an interger for convenience
    
    return run_num, unix_time, date_time

@click.command()
@click.option('-s', '--station', type = int, help = 'any ara station id ex) 2 or 3')
@click.option('-o', '--output', type = str, help = 'ex) /home/mkim/')
@click.option('-b', '--blinded', default = False, type = bool, help = 'ex) 0: unblinded or 1: blinded')
@click.option('-p', '--pad_len', default = 22000, type = int, help = 'ex) 22000 runs')
@click.option('-i', '--year_start', default = 2013, type = int, help = 'ex) 2013')
@click.option('-f', '--year_end', default = 2020, type = int, help = 'ex) 2020')
def main(station, output, blinded, pad_len, year_start, year_end):
    """! main function for collecting configuration
    This is designed to collect data from configFile.run[xxxxxx].dat, runStart.run[xxxxxx].dat, and runStop.run[xxxxxx].dat
    This script needs 1) h5df, 2) tqdm, and 3) click package
    Please intstall aobve three package before run the code
    'pip install <package>'  
 
    @param station  integer. station id ex) 1 ~ 5. not for PA at this moment
    @param output  string. desired path for saving output
    """

    print('Collectiong file paths from each year!')
    ## all the file name variables
    config_dat_str = 'configFile.run[0-9]*.dat'
    run_start_str = 'runStart.run[0-9]*.dat'
    run_stop_str = 'runStop.run[0-9]*.dat'
    data_path = '/data/exp/ARA/'

    ## list for all file path
    config_list = [] # saving all configFile.run[xxxxxx].dat paths
    run_start_list = [] # saving for all runStart.run[xxxxxx].dat paths
    run_stop_list = [] # saving for all runStop.run[xxxxxx].dat paths
    run_list = [] # saving for all run paths
    evt_list = [] # saving for all event paths
    ped_list = [] # saving for all pedestal paths

    ## Seach the config files by loop over all L1 data path (2013 ~ 2019) and save it into lists
    if blinded:
        type_path_2013 = 'full2013Data'
        type_path = 'blinded'
    else:
        type_path_2013 = 'burnSample1in10'
        type_path = 'unblinded'
    for y in tqdm(range(year_start, year_end + 1)):
        yrs = int(y)
        if yrs == 2013: # 2013 has different file path 
            #g_path = f'{yrs}/filtered/unzippedTGZFiles/ARA0{station}/run_[0-9]*/logs/' # alternate path
            c_path = f'{yrs}/raw/ARA0{station}-SPS-ARA/run_[0-9]*/logs/'
            r_path = f'{yrs}/filtered/{type_path_2013}/ARA0{station}/root/run*'
            p_path = f'/data/user/mkim/ARA_2013_Ped/ARA0{station}/' # use this nicely unpacked pedestal for scrapping
            ped_list += glob(f'{p_path}pedestalValues*')
            paths = [c_path] # event and pestal path. At this moment, not sure where is config for pedestal...
        else:
            c_path = f'{yrs}/{type_path}/L1/ARA0{station}/*/run*/' # event path
            r_path = c_path[:-1] # run path
            p_path = f'{yrs}/calibration/pedestals/ARA0{station}/' # pestal path
            ped_list += glob(f'{data_path}{p_path}pedestalValues*')
            paths = [c_path, p_path] # event and pestal path
        run_list += glob(f'{data_path}{r_path}')
        evt_list += glob(f'{data_path}{r_path}/event[0-9]*')
        for p in range(len(paths)):
            config_list += glob(f'{data_path}{paths[p]}{config_dat_str}')
            run_start_list += glob(f'{data_path}{paths[p]}{run_start_str}')
            run_stop_list += glob(f'{data_path}{paths[p]}{run_stop_str}')

    run_len = len(run_list)
    evt_len = len(evt_list)
    ped_len = len(ped_list)
    config_len = len(config_list)
    run_start_len = len(run_start_list)
    run_stop_len = len(run_stop_list)
    print('Total # of run paths:', run_len)
    print('Total # of event files:', evt_len)
    print('Total # of pedestal files:', ped_len)
    print('Total # of config files (event + pedestal):', config_len)
    print('Total # of runstart files (event + pedestal):', run_start_len)
    print('Total # of runstop files (event + pedestal):', run_stop_len)

    print('Collecting information. 6 for loop')
    ## giant numpy array pad for storing all the configuration from each run's config files
    ## if there is no run config for corresponding run number. array element would be Nan
    countings = np.array([run_len, evt_len, ped_len, config_len, run_start_len, run_stop_len], dtype = int) # total number of run, event, pedestal, config, runstart, and runstop counting
    run_num = np.arange(pad_len, dtype = int) 
    run_num_len = len(run_num) 
    run_type = np.full((run_num_len), np.nan, dtype = float) # whether run number is corresponding to event or pedestal. 0: pedestal, 1: event, 2: only run path but no event file, and nan: There are no event or pedetal run for corresponding run number
    unix_time = np.full((2, run_num_len), np.nan, dtype = float) # array for unix time stored in runStart.run[xxxxxx].dat and runStop.run[xxxxxx].dat
    date_time = np.copy(unix_time) # array for UTC date time converted from unix time
    rf_block_num = np.full((run_num_len), np.nan, dtype = float) # number of RF blocks conifg stored in each run config file. we assuming 1 bkock = 20 ns
    soft_block_num = np.copy(rf_block_num) # number of Software blocks config stroed in each run config file
    trig_win_num = np.copy(rf_block_num) # number of trigger window length stired in each run config file
    delay_enable  = np.copy(rf_block_num) # trigger delay enable. 1: enable, 0: not enable
    delay_num = np.full((16, run_num_len), np.nan, dtype = float) # number of trigger delay config stored in each run config file
    masked_ant = np.full((20, run_num_len), np.nan, dtype = float) # which antenn is masked. 1: unmasked, 0: masked. This mapping is following trigger channel mapping
    scaler_goal = np.copy(delay_num) # L1 servo goal value. unit is Hz. I guess...
    calpulser_info = np.full((6, run_num_len), np.nan, dtype = float) # info for calpulser configuration. antennaIceA#I1=, antennaIceB#I1=, opIceA#I1=, opIceB#I1=0, attIceA#I1=, and attIceB#I1=
    ## related document for calpulser configuration... https://aradocs.wipac.wisc.edu/0005/000502/007/icecalsoftware.pdf

    ## 5 each 'for loop' for collecting information. Since not all the runs have a config, runstart, or runstop, we need to scrap the information from each list      
    ## variables for searching config information
    masked_ant_key = 'enableL1Trigger#I20='
    rf_blk_key = 'numRF0TriggerBlocks#1='
    soft_blk_key = 'numSoftTriggerBlocks#1='
    trig_win_key = 'triggerWindowSize#1='
    delay_enable_key = 'enableTriggerDelays#I1='
    delay_key = 'triggerDelays#I16='
    goal_key = 'scalerGoalValues#I16='
    calpulser_key = ['antennaIceA#I1=', 'antennaIceB#I1=', 'opIceA#I1=', 'opIceB#I1=', 'attIceA#I1=', 'attIceB#I1=']
    cal_key_len = len(calpulser_key)
    evt_key = 'event'
    evt_key_len = len(evt_key)
    run_key = 'run'
    run_key_len = len(run_key)
    root_key = '.root'
    dat_key = '.dat'

    ## 1st event files
    for evts in tqdm(range(evt_len)):
        run_idx = evt_list[evts].find(evt_key) + evt_key_len          
        format_idx = evt_list[evts].find(root_key, run_idx)           
        val = int(evt_list[evts][run_idx:format_idx])
        run_type[val] = 1        
        
    ## 2nd pedestal files
    for peds in tqdm(range(ped_len)):
        run_idx = ped_list[peds].find(run_key) + run_key_len
        format_idx = ped_list[peds].find(dat_key, run_idx)
        val = int(ped_list[peds][run_idx:format_idx])
        run_type[val] = 0

    ## 3rd check the run that has no event file...
    for runs in tqdm(range(run_len)):
        run_idx = run_list[runs].find(run_key) + run_key_len
        val = int(run_list[runs][run_idx:])
        if np.isnan(run_type[val]):
            run_type[val] = 2
            print(f'{run_list[runs]} has no event file!!')

    ## 4th run start files
    for runs in tqdm(range(run_start_len)):
        run_start, unix_start, date_start = scrap_run_dat_info(run_start_list[runs])
        unix_time[0, run_start] = unix_start
        date_time[0, run_start] = date_start
 
    ## 5th run stop files. I believe if the run is suddenly terminated, that run usually dont have a runstop files 
    for runs in tqdm(range(run_stop_len)):
        run_stop, unix_stop, date_stop = scrap_run_dat_info(run_stop_list[runs])
        unix_time[1, run_stop] = unix_stop
        date_time[1, run_stop] = date_stop
 
    ## 6th config files
    for runs in tqdm(range(config_len)):
        run_config = int(config_list[runs][-10:-4]) # scrap the run number from config file path
        
        ## open the each config file and scrap the information by context_finder function
        with open(config_list[runs],'r') as config_file:
            config_read = config_file.read()

            rf_block_num[run_config] = context_finder(config_read, rf_blk_key)
            soft_block_num[run_config] = context_finder(config_read, soft_blk_key)
            delay_enable[run_config] = context_finder(config_read, delay_enable_key)
            masked_ant[:, run_config] = context_finder(config_read, masked_ant_key, num_vals = 20)
            scaler_goal[:, run_config] = context_finder(config_read, goal_key, num_vals = 16)
            for cal in range(cal_key_len):
                calpulser_info[cal, run_config] = context_finder(config_read, calpulser_key[cal])
            
            ## sometime ara4/5 doesn't have a ';' at the end of values. use diffeent key for scrapping infromation
            try:
                trig_win_num[run_config] = context_finder(config_read, trig_win_key)
            except ValueError:
                trig_win_num[run_config] = context_finder(config_read, trig_win_key, end_key = '//')
            try: 
                delay_num[:, run_config] = context_finder(config_read, delay_key, num_vals = 16)
            except ValueError:
                delay_num[:, run_config] = context_finder(config_read, delay_key, end_key = '<', num_vals = 16)

    print('Saving information')
    ## create output dir
    if not os.path.exists(output):
        os.makedirs(output)

    ## save into h5 file
    h5_file_name=f'{output}Config_A{station}_{type_path}_{year_start}_to_{year_end}.h5'
    hf = h5py.File(h5_file_name, 'w')
    hf.create_dataset('run_num', data=run_num, compression="gzip", compression_opts=9)
    hf.create_dataset('run_type', data=run_type, compression="gzip", compression_opts=9)
    hf.create_dataset('unix_time', data=unix_time, compression="gzip", compression_opts=9)
    hf.create_dataset('date_time', data=date_time, compression="gzip", compression_opts=9)
    hf.create_dataset('masked_ant', data=masked_ant, compression="gzip", compression_opts=9)
    hf.create_dataset('rf_block_num', data=rf_block_num, compression="gzip", compression_opts=9)
    hf.create_dataset('soft_block_num', data=soft_block_num, compression="gzip", compression_opts=9)
    hf.create_dataset('trig_win_num', data=trig_win_num, compression="gzip", compression_opts=9)
    hf.create_dataset('delay_enable', data=delay_enable, compression="gzip", compression_opts=9)
    hf.create_dataset('delay_num', data=delay_num, compression="gzip", compression_opts=9)
    hf.create_dataset('scaler_goal', data=scaler_goal, compression="gzip", compression_opts=9)
    hf.create_dataset('calpulser_info', data=calpulser_info, compression="gzip", compression_opts=9)
    hf.create_dataset('countings', data=countings, compression="gzip", compression_opts=9)
    hf.close() 

    print(f'output is {h5_file_name}')
    print('Done!')

if __name__ == "__main__":

    ## excute main function
    main()
