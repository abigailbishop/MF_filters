"""
For a provided station and blind vs unblinded data, compile the following data 
  from sub_info files: 
    For each run:  
        run number, config number, and number of events
    For each event in each run: 
        run number, config number, config number, unix time, trigger type

Arguments
---------
Station : int
    Station number
Blind : int
    Whether to use blind or unblinded data. 
      For 100% data, provide `1`. For 10% data, provide `0`.

Outputs
-------
hf : h5py.File
    File titled `Info_Summary_*.h5` available in 
      `$OUTPUT_PATH/ARA0{Station}/Hist/`
"""

import numpy as np
import os, sys
import h5py
from tqdm import tqdm

curr_path = os.getcwd()
sys.path.append(curr_path+'/../')
from tools.ara_run_manager import file_sorter
from tools.ara_utility import size_checker

Station = int(sys.argv[1])
Blind = int(sys.argv[2])

if Station == 2: num_configs = 7
if Station == 3: num_configs = 9

# Get list of runs available in sub_info directory (where event data was collected)
if Blind == 1:
    burn_key = '_full'
    burn_name = '_Full'
else: 
    burn_key = '_burn'
    burn_name = ''
d_path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/sub_info{burn_key}/*'
d_list, d_run_tot, d_run_range, d_len = file_sorter(d_path)
del d_run_range

# Create objects to track data for each run available in `d_path`
runs = np.copy(d_run_tot)                  # Run numbers (sorted already)
configs = np.full((d_len), 0, dtype = int) # Configuration number
num_evts = np.copy(configs)                # Number of events

# Anticipated number of events in all runs for any station. 
# May need to be increased as time passes and ARA stations collect more data. 
pad_len = 120_000_000
if Blind == 1: 
    pad_len = 140_000_000
    pad_len *= 10

# Create objects to track data for each event
# Assumes a massively high number of events. Will sift out entries that 
#   were not filled with actual data later.
run_ep = np.full((pad_len), 0, dtype = int) # Run numbers event lives in
evt_ep = np.copy(run_ep)                    # Event number of this event
trig_ep = np.copy(run_ep)                   # Trigger type of this event
con_ep = np.copy(run_ep)                    # Configuration number for this event
unix_ep = np.copy(run_ep)                   # Unix time for this event

# Loop over files for all runs available in `d_path` 
#   and extract data for each event
count_b = 0 # used to track start index to save event data to in aforementioned arrays
count_a = 0 # used to track last index to save event data to in aforementioned arrays
print(f"Analyzing {d_len} runs: ", end="", flush=True)
for r in tqdm(range(len(d_run_tot))):
    if r%25 == 0: print(r, end="", flush=True)
    else: print(">", end="", flush=True)
    
  #if r > 960:

    # Extract data from the file
    hf = h5py.File(d_list[r], 'r')
    configs[r] = hf['config'][2]
    evt = hf['evt_num'][:]
    trig_type = hf['trig_type'][:]
    unix_time = hf['unix_time'][:]
    num_evts_r = len(evt)
    num_evts[r] = num_evts_r

    # Save values extracted from the `hf` h5 file
    count_a += num_evts_r
    if count_a > run_ep.size: 
        raise IndexError(f"Arrays to save data to are too small. "
                         f"Update `pad_len` variable to be larger.")
    run_ep[count_b:count_a] = d_run_tot[r]
    con_ep[count_b:count_a] = configs[r]
    evt_ep[count_b:count_a] = evt
    trig_ep[count_b:count_a] = trig_type
    unix_ep[count_b:count_a] = unix_time
    count_b += num_evts_r
    del hf, evt, trig_type, unix_time, num_evts_r
print()

# Count full number of anticipated events available and trim all 
#   oversized arrays to just include this many events 
tot_runs = int(np.nansum(num_evts))
print(tot_runs)
run_ep = run_ep[:tot_runs]
evt_ep = evt_ep[:tot_runs]
trig_ep = trig_ep[:tot_runs]
con_ep = con_ep[:tot_runs]
unix_ep = unix_ep[:tot_runs]
print(runs.shape)
print(run_ep.shape)

# Identify/create the directory that will store outputs from summary makers
path = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/Hist/'
if not os.path.exists(path):
    os.makedirs(path)
os.chdir(path)

# Save information from each event to an h5 file
file_name = f'Info_Summary{burn_name}_A{Station}.h5'
hf = h5py.File(file_name, 'w')
hf.create_dataset(
    'runs', data=runs, compression="gzip", compression_opts=9
) # Run numbers (size = n_runs)
hf.create_dataset(
    'num_evts', data=num_evts, compression="gzip", compression_opts=9
) # Number of events in each run (size = n_runs)
hf.create_dataset(
    'configs', data=configs, compression="gzip", compression_opts=9
) # Livetime configuration number for each run (size = n_runs)
hf.create_dataset(
    'run_ep', data=run_ep, compression="gzip", compression_opts=9
) # Run number for each event in full event list (all events in all runs from d_path)\
# (size = n_events)
hf.create_dataset(
    'evt_ep', data=evt_ep, compression="gzip", compression_opts=9
) # Event number for each event in full event list (size = n_events)
hf.create_dataset(
    'trig_ep', data=trig_ep, compression="gzip", compression_opts=9
) # Trigger type for each event in full event list (size = n_events)
hf.create_dataset(
    'con_ep', data=con_ep, compression="gzip", compression_opts=9
) # Configuration number for each event in full event list (size = n_events)
hf.create_dataset(
    'unix_ep', data=unix_ep, compression="gzip", compression_opts=9
) # Unix time for each event in the full event list (size = n_events)
hf.close()
print('file is in:',path+file_name, size_checker(path+file_name))
print('done!')