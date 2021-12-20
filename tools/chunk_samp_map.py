import os, sys
import numpy as np
from tqdm import tqdm

def samp_map_collector_dat(Data, Ped, Station, Year, clean_info = None):

    print('Collecting wf starts!')

    from tools.ara_data_load import ara_uproot_loader
    from tools.ara_data_load import ara_root_loader
    from tools.ara_data_load import analog_buffer_info_loader
    from tools.constant import ara_const

    # geom. info.
    ara_const = ara_const()
    num_Ants = ara_const.USEFUL_CHAN_PER_STATION
    num_Buffers = ara_const.SAMPLES_PER_DDA
    num_Bits = ara_const.BUFFER_BIT_RANGE
    num_Strs = ara_const.DDA_PER_ATRI
    del ara_const

    # data config
    buffer_info = analog_buffer_info_loader(Station, Year, incl_cable_delay = True)
    ara_root = ara_root_loader(Data, Ped, Station, Year)
    ara_uproot = ara_uproot_loader(Data)
    ara_uproot.get_sub_info()

    # quality cut results
    if clean_info is not None:
        clean_evt, clean_entry, clean_st = clean_info.values()
    else:
        clean_evt = ara_uproot.evt_num
        clean_entry = ara_uproot.entry_num
        clean_st = np.full((num_Strs, len(clean_evt)), 0, dtype = int)

    # output array
    samp_map = np.full((num_Buffers, num_Bits, num_Ants), 0, dtype = int)

    # loop over the events
    for evt in tqdm(range(len(clean_evt))):
      #if clean_entry[evt] < 5000:
        # get entry and wf
        ara_root.get_entry(clean_entry[evt])
        ara_root.get_useful_evt()

        # sample index
        blk_idx_arr = ara_uproot.get_block_idx(clean_entry[evt], trim_1st_blk = True)[0]
        samp_idx = buffer_info.get_samp_idx(blk_idx_arr, ch_shape = True)   
        del blk_idx_arr

        # loop over the antennas
        for ant in range(num_Ants):

            if clean_st[ant%num_Strs, evt] != 0:
                print('not good string!', ant%num_Strs)
                continue       
 
            # stack in sample map
            samp_idx_ant = samp_idx[:,ant][~np.isnan(samp_idx[:,ant])].astype(int)
            raw_v = ara_root.get_rf_ch_wf(ant)[1].astype(int)
            samp_map[samp_idx_ant, raw_v, ant] += 1
            del samp_idx_ant, raw_v
        del samp_idx
    del ara_root, ara_uproot, clean_entry, num_Strs
  
    samp_medi = np.full((num_Buffers, num_Ants), np.nan, dtype = float)
    buffer_bit_range = np.arange(num_Bits)
    bin_width = np.diff(buffer_bit_range)[0]
    for sam in tqdm(range(num_Buffers)):
        for ant in range(num_Ants): 

            samp_medi[sam, ant] = buffer_info.get_median_from_hist(samp_map[sam, :, ant], bin_center = buffer_bit_range, bin_width = bin_width)

    buffer_sample_range = np.arange(num_Buffers)
    del num_Ants, num_Buffers, num_Bits, bin_width, buffer_info 

    print('WF collecting is done!')

    return {'buffer_bit_range':buffer_bit_range,
            'buffer_sample_range':buffer_sample_range,
            'samp_map':samp_map, 
            'samp_medi':samp_medi, 
            'clean_evt':clean_evt}








