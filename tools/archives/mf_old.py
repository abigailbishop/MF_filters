import os, sys
import numpy as np
from scipy.signal import hilbert
from tqdm import tqdm

# custom lib
from tools.ara_root import useful_evt_maker
#from tools.ara_root import useful_evt_maker_del
from tools.wf import station_pad

def off_pad_maker(pad_len, dt):

    off_pad = np.arange(0,pad_len,1) * dt

    return off_pad, len(off_pad), off_pad[0], off_pad[-1]

def Band_Square(freq, band_amp, num_temp = 1):
    
    #Bothpol
    band_amp[(freq>=430e6) & (freq<=480e6)]=1e-50
    band_amp[(freq<=-430e6) & (freq>=-480e6)]=1e-50
    
    if num_temp >= 1:
        band_amp[(freq>=-150e6) & (freq<=150e6),:8]=1e-50
        band_amp[(freq>=-217e6) & (freq<=217e6),8:]=1e-50
    
    if num_temp >= 2:
        band_amp[(freq>=680e6) | (freq<=-680e6),:8,1]=1e-50
        band_amp[(freq>=630e6) | (freq<=-630e6),8:,1]=1e-50

    if num_temp >= 3:
        band_amp[(freq>=665e6) | (freq<=-665e6),:8,2]=1e-50
        band_amp[(freq>=530e6) | (freq<=-530e6),8:,2]=1e-50

    #Both theta
    band_amp[(freq>=680e6) | (freq<=-680e6),:,0]=1e-50   
 
    return(band_amp)

def Matched_filter_SNR(dat, temp, psd, dt, df, hil = False):

    snr = temp * dat.conjugate()
    snr /= psd
    snr = np.abs(2*np.fft.ifft(snr, axis = 0) / dt)

    if hil == True:
        snr = hilbert(snr.real(), axis = 0)
    else:
        pass

    snr /= np.sqrt(np.abs(2 * np.nansum(temp * temp.conjugate() / psd, axis = 0) * df))[np.newaxis,:,np.newaxis]

    print('snr making is done!')

    return snr
        
def OMF_v2(psd, rf_f_w, rf_v, temp_v, ndf, t_pad_l):
    
    # matched filtering
    conjugation = temp_v * rf_v.conjugate()
    conjugation /= psd
    #conjugation = hilbert(np.abs(2*np.fft.ifft(conjugation,axis=0)).real, axis = 0)
    #conjugation = np.abs(2*np.fft.ifft(conjugation,axis=0))
    conjugation = 2*np.fft.ifft(conjugation,axis=0) * ndf
    conjugation = conjugation[::-1,:,:]
    #conjugation[conjugation<0] = 0
    
    nor_fac_sq = np.abs(2*np.nansum(temp_v * temp_v.conjugate() / psd,axis=0) * rf_f_w)
    nor_fac_sq = np.repeat(nor_fac_sq[np.newaxis,:,:], t_pad_l, axis=0) 

    return conjugation, nor_fac_sq

def soft_psd_maker(R, evtTree, rawEvt, num_evts, cal, q # ara root
                    , num_Ants, bad_ant_i, t_width_ns, ndf # known config
                    , t_pad_len, time_pad_i, time_pad_f # time
                    , f # freq
                    , n_theta # theta
                    , DMode, Station = None, Run = None, Output = None, sel_evt = None): # argv

    print('PSD making starts!')
    # custom lib
    from tools.fft import psd_maker
    from tools.array import arr_2d

    # array for psd
    psd = arr_2d(t_pad_len, num_Ants, 0, float)
    #psd = arr_2d(t_pad_len, num_Ants, 0, complex)
    num_psd = 0

    # loop over the events
    for event in tqdm(range(num_evts)):

        # make a useful event
        usefulEvent = useful_evt_maker(R, evtTree, rawEvt, event, cal)

        # trigger filter
        #if rawEvt.isSoftwareTrigger() == 0 and rawEvt.isCalpulserEvent() == 0:
        if rawEvt.isSoftwareTrigger() == 1:

            # quality cut
            if q.isGoodEvent(usefulEvent) == 1:

                # make padded wf and interpolated wf length
                ant_arr, int_time_len = station_pad(usefulEvent, num_Ants, t_width_ns, t_pad_len, time_pad_i, time_pad_f)

                # make psd
                psd += psd_maker(ant_arr, ndf, t_pad_len, int_time_len)
                num_psd += 1 

                if event == sel_evt and DMode == 'debug':

                    from tools.debug import psd_indi_debug
                    ant_arr_copy, ant_arr_fft, ant_arr_fft_band, indi_psd, indi_psd_band = psd_indi_debug(Station, Run, Output, event
                            , t_pad_len, time_pad_i, time_pad_f, t_width_ns, ndf
                            , f
                            , n_theta
                            , ant_arr, int_time_len)
                    selected_evt = np.copy(sel_evt)

                elif sel_evt is None and DMode == 'debug' and num_psd == 10:

                    print(f'event# {event} ({num_psd}th soft event) is selected for the debug')

                    from tools.debug import psd_indi_debug
                    ant_arr_copy, ant_arr_fft, ant_arr_fft_band, indi_psd, indi_psd_band = psd_indi_debug(Station, Run, Output, event
                            , t_pad_len, time_pad_i, time_pad_f, t_width_ns, ndf
                            , f
                            , n_theta
                            , ant_arr, int_time_len)
                    selected_evt = np.copy(event)

                else:
                    pass

                del ant_arr, int_time_len

        #R.~UsefulAtriStationEvent()
        #useful_evt_maker_del(R, usefulEvent) 
        #usefulEvent.Delete()
        del usefulEvent
        #usefulEvent.Delete()

    # averaging
    psd /= num_psd
    print('The number of soft triggered events is', num_psd)

    if DMode == 'debug':

        from tools.debug import psd_debug
        psd_copy = psd_debug(Station, Run, Output
                , f
                , n_theta
                , psd, num_psd)

    else:
        pass

    # band pass filter
    psd_wo_band = np.copy(psd)
    psd = Band_Square(f, np.repeat(psd[:,:,np.newaxis], n_theta, axis=2))

    # remove bad antenna
    try:
        psd[:, bad_ant_i, :] = np.nan
    except IndexError:
        pass

    print('PSD making is done!')

    if DMode == 'debug':
        del psd_wo_band
        return psd, psd_copy, num_psd, ant_arr_copy, ant_arr_fft, ant_arr_fft_band, indi_psd, indi_psd_band, selected_evt

    elif DMode == 'normal':
    #else:
        del num_psd
        return psd, psd_wo_band

    else:
 
        print('DMode is not set. Choose 1) normal or 2) debug w/ events')
        sys.exit(1)

def evt_snr_maker(R, evtTree, rawEvt, num_evts, cal, q # ara root
                    , num_Ants, bad_ant_i, t_width_ns, ndf # known config
                    , t_pad_l, time_pad_i, time_pad_f # time
                    , f, f_w # freq
                    , n_psd # psd
                    , temp_v, n_theta, theta_w, peak_i # temp
                    , DMode, Station, Year, CPath, Run = None, Output = None, sel_evt = None): # argv

    print('Event-wise SNR making starts!')
    from scipy.ndimage import maximum_filter1d

    #custom lib
    from tools.ara_root import trig_checker
    from tools.arr_table import table_loader
    from tools.array import arr_3d

    # table_loader
    peak_w = 100
    if DMode == 'debug':
        mov_i, pad_t_l, p_len_front, p_len_end, ps_len_i, mov_t, pad_t = table_loader(CPath, Station, Year, theta_w, peak_w)
    else:
        mov_i, pad_t_l, p_len_front, p_len_end, ps_len_i = table_loader(CPath, Station, Year, theta_w, peak_w)[:-2]

    # tale remove range
    tale_i_front= -1*peak_i + p_len_front
    tale_i_end= peak_i + p_len_end

    # array for snr
    snr_wf = arr_3d(pad_t_l, num_Ants, n_theta, 0, float)
    snr_wf_nor = np.copy(snr_wf)
    snr_wf_01 = arr_3d(pad_t_l, num_Ants, n_theta, 0, int)

    # array for 2d map
    snr_wf_2d_v = arr_3d(mov_i.shape[0], mov_i.shape[2], mov_i.shape[3], 0, float)
    snr_wf_2d_h = np.copy(snr_wf_2d_v)
    snr_wf_2d_01_v = np.copy(snr_wf_2d_v)
    snr_wf_2d_01_h = np.copy(snr_wf_2d_v)

    #array for event-wise snr
    if DMode == 'debug':
        evt_snr = []
    else:
        pass
    evt_snr_v = []
    evt_snr_h = []
    evt_num = []
    trigger = []

    # half antenna number
    half_ant = int(num_Ants/2)

    # loop over the events
    for event in tqdm(range(num_evts)):

        # make a useful event
        usefulEvent = useful_evt_maker(R, evtTree, rawEvt, event, cal)

        # trigger filter
        #if rawEvt.isSoftwareTrigger() == 0:

        # quality cut
        if q.isGoodEvent(usefulEvent) == 1:

                # make padded wf with interpolated wf to fft
                ant_arr = station_pad(usefulEvent, num_Ants, t_width_ns, t_pad_l, time_pad_i, time_pad_f)[0]

                # OMF
                snr_wf[:] = 0
                snr_wf_nor[:] = 0
                snr_wf[p_len_front:-p_len_end], snr_wf_nor[p_len_front:-p_len_end] = OMF_v2(n_psd, f_w, Band_Square(f, np.repeat((np.fft.fft(ant_arr, axis=0)/ndf)[:,:,np.newaxis], n_theta, axis=2)), temp_v, ndf, t_pad_l)

                if event == sel_evt and DMode == 'debug':

                    trig_index = trig_checker(rawEvt)
                    if trig_index == 0:
                        trig_type = 'RF'
                    elif trig_index == 1:
                        trig_type = 'Cal'
                    elif trig_index == 2:
                        trig_type = 'Soft'

                    print(f'Evt#{event} is selected!')
                    print('Trigger type is',trig_type) 

                    from tools.debug import evt_snr_indi_debug_0
                    ant_arr_copy, ant_arr_fft, ant_arr_fft_band, snr_wf_copy = evt_snr_indi_debug_0(Station, Run, Output, event, trig_type
                            , time_pad_i, time_pad_f, t_width_ns, ndf
                            , f
                            , n_theta, pad_t
                            , ant_arr, snr_wf, snr_wf_nor)

                else:
                    pass

                # remove tale step 1
                ant_arr[ant_arr != 0] = 1
                snr_wf_01[:] = 0
                for n_t in range(n_theta):
                    snr_wf_01[tale_i_front[0][n_t]:-tale_i_end[0][n_t], :half_ant, n_t] = ant_arr[:, :half_ant]
                    snr_wf_01[tale_i_front[1][n_t]:-tale_i_end[1][n_t], half_ant:, n_t] = ant_arr[:, half_ant:]
                del ant_arr

                # remove tale step 2
                snr_wf[snr_wf_01 == 0] = 0
                snr_wf_nor[snr_wf_01 == 0] = 0

                if event == sel_evt and DMode == 'debug':

                    from tools.debug import evt_snr_indi_debug_1
                    snr_wf_copy_1 = evt_snr_indi_debug_1(Station, Run, Output, event, trig_type
                            , pad_t
                            , snr_wf, snr_wf_nor)
                
                else:
                    pass

                # picking max
                snr_wf = maximum_filter1d(np.abs(snr_wf), size=ps_len_i, axis=0, mode='constant')
                snr_wf_nor = maximum_filter1d(snr_wf_nor, size=ps_len_i, axis=0, mode='constant')

                if event == sel_evt and DMode == 'debug':

                    from tools.debug import evt_snr_indi_debug_2
                    snr_wf_r_max_copy = evt_snr_indi_debug_2(Station, Run, Output, event, trig_type
                            , pad_t, peak_w
                            , snr_wf, snr_wf_nor)

                else:
                    pass
                
                # remove bad antenna
                try:
                    snr_wf[:, bad_ant_i, :] = 0
                    snr_wf_nor[:, bad_ant_i, :] = 0
                except IndexError:
                    pass

                # 2d map array
                snr_wf_2d_v[:] = 0
                snr_wf_2d_h[:] = 0
                snr_wf_2d_01_v[:] = 0
                snr_wf_2d_01_h[:] = 0

                # stack snr into 2d map
                for half in range(half_ant):

                    snr_wf_2d_v[:,2] += snr_wf[:,half,0][mov_i[:,half,2,:]]
                    snr_wf_2d_v[:,3] += snr_wf[:,half,1][mov_i[:,half,3,:]]
                    snr_wf_2d_v[:,1] += snr_wf[:,half,1][mov_i[:,half,1,:]]
                    snr_wf_2d_v[:,0] += snr_wf[:,half,2][mov_i[:,half,0,:]]
                    snr_wf_2d_v[:,4] += snr_wf[:,half,2][mov_i[:,half,4,:]]

                    snr_wf_2d_01_v[:,2] += snr_wf_nor[:,half,0][mov_i[:,half,2,:]]
                    snr_wf_2d_01_v[:,3] += snr_wf_nor[:,half,1][mov_i[:,half,3,:]]
                    snr_wf_2d_01_v[:,1] += snr_wf_nor[:,half,1][mov_i[:,half,1,:]]
                    snr_wf_2d_01_v[:,0] += snr_wf_nor[:,half,2][mov_i[:,half,0,:]]
                    snr_wf_2d_01_v[:,4] += snr_wf_nor[:,half,2][mov_i[:,half,4,:]]

                    snr_wf_2d_h[:,2] += snr_wf[:,half+half_ant,0][mov_i[:,half+half_ant,2,:]]
                    snr_wf_2d_h[:,3] += snr_wf[:,half+half_ant,1][mov_i[:,half+half_ant,3,:]]
                    snr_wf_2d_h[:,1] += snr_wf[:,half+half_ant,1][mov_i[:,half+half_ant,1,:]]
                    snr_wf_2d_h[:,0] += snr_wf[:,half+half_ant,2][mov_i[:,half+half_ant,0,:]]
                    snr_wf_2d_h[:,4] += snr_wf[:,half+half_ant,2][mov_i[:,half+half_ant,4,:]]

                    snr_wf_2d_01_h[:,2] += snr_wf_nor[:,half+half_ant,0][mov_i[:,half+half_ant,2,:]]
                    snr_wf_2d_01_h[:,3] += snr_wf_nor[:,half+half_ant,1][mov_i[:,half+half_ant,3,:]]
                    snr_wf_2d_01_h[:,1] += snr_wf_nor[:,half+half_ant,1][mov_i[:,half+half_ant,1,:]]
                    snr_wf_2d_01_h[:,0] += snr_wf_nor[:,half+half_ant,2][mov_i[:,half+half_ant,0,:]]
                    snr_wf_2d_01_h[:,4] += snr_wf_nor[:,half+half_ant,2][mov_i[:,half+half_ant,4,:]]

                if DMode == 'debug':
                    evt_snr.append(np.nanmax(np.sqrt((snr_wf_2d_v + snr_wf_2d_h)**2 / (snr_wf_2d_01_v + snr_wf_2d_01_h))))
                else:
                    pass
                evt_snr_v.append(np.nanmax(np.sqrt(snr_wf_2d_v**2 / snr_wf_2d_01_v)))
                evt_snr_h.append(np.nanmax(np.sqrt(snr_wf_2d_h**2 / snr_wf_2d_01_h)))
                evt_num.append(event)
                trigger.append(trig_checker(rawEvt))

                if event == sel_evt and DMode == 'debug':

                    from tools.debug import evt_snr_indi_debug_3
                    v_match, h_match, v_sum_match, h_sum_match, v_sum_match_01, h_sum_match_01, v_avg_match, h_avg_match, evt_snr_v_sky_2d, evt_snr_h_sky_2d, nadir_range, phi_range, v_opt_angle, h_opt_angle = evt_snr_indi_debug_3(Station, Run, Output, event, trig_type
                            , num_Ants, half_ant
                            , mov_i, mov_t, theta_w, peak_w
                            , snr_wf, snr_wf_nor)

                else:
                    pass

        #R.~UsefulAtriStationEvent()
        del usefulEvent

    del mov_i, pad_t_l, p_len_front, p_len_end, ps_len_i, tale_i_front, tale_i_end, snr_wf, snr_wf_nor, snr_wf_01, snr_wf_2d_v, snr_wf_2d_01_v, snr_wf_2d_h, snr_wf_2d_01_h, half_ant

    if DMode == 'debug':
        evt_snr = np.asarray(evt_snr)
    else:
        pass
    evt_snr_v = np.asarray(evt_snr_v)
    evt_snr_h = np.asarray(evt_snr_h)
    evt_num = np.asarray(evt_num)
    trigger = np.asarray(trigger)
 
    if DMode == 'debug':

        from tools.debug import evt_snr_debug
        evt_snr_debug(Station, Run, Output, trigger, evt_snr_v, evt_snr_h)

    else:
        pass

    print('Event-wise SNR making is done!')

    if DMode == 'debug':
   
        del trig_type 
        return evt_snr, evt_snr_v, evt_snr_h, evt_num, trigger, trig_index, ant_arr_copy, ant_arr_fft, ant_arr_fft_band, snr_wf_copy, snr_wf_copy_1, snr_wf_r_max_copy, v_match, h_match, v_sum_match, h_sum_match, v_sum_match_01, h_sum_match_01, v_avg_match, h_avg_match, evt_snr_v_sky_2d, evt_snr_h_sky_2d, nadir_range, phi_range, np.arange(time_pad_i, time_pad_f+t_width_ns, t_width_ns), mov_t, pad_t, peak_w, v_opt_angle, h_opt_angle

    elif DMode == 'normal':
    #else:
        del peak_w
        return evt_snr_v, evt_snr_h, evt_num, trigger
                
    else:

        print('DMode is not set. Choose 1) normal or 2) debug w/ events')
        sys.exit(1)
"""
def indi_snr_maker(Station, Run, Output # argv
                    , R, evtTree, rawEvt, num_evts, cal, q # ara root
                    , num_Ants, bad_ant_i, t_width_ns # known config
                    , t_pad_l, time_pad_i, time_pad_f # time
                    , f, f_w # freq
                    , n_psd # psd
                    , temp_v, n_theta): # temp

    print('Indi SNR making starts!')

    import h5py

    # create output file
    os.chdir(Output)
    h5_file_name='Indi_SNR_A'+str(Station)+'_R'+str(Run)+'.h5'
    hf = h5py.File(h5_file_name, 'w')

    trigger = []

    # loop over the events
    for event in tqdm(range(num_evts)):

        # make a useful event
        usefulEvent = useful_evt_maker(R, evtTree, rawEvt, event, cal)

        # trigger filter
        if rawEvt.isSoftwareTrigger() == 0:

            # quality cut
            if q.isGoodEvent(usefulEvent) == 1:

                # make padded wf with interpolated wf to fft
                ant_arr = station_pad(usefulEvent, num_Ants, t_width_ns, t_pad_l, time_pad_i, time_pad_f)

                # OMF
                snr_wf = OMF(n_psd, f_w, Band_Square(f, np.repeat(np.fft.fft(ant_arr, axis=0)[:,:,np.newaxis], n_theta, axis=2)), temp_v)
                del ant_arr

                # remove bad antenna
                snr_wf[:, bad_ant_i, :] = np.nan

                #saving result
                hf.create_dataset(str(event), data=snr_wf, compression="gzip", compression_opts=9)
                del snr_wf

                # trigger check
                trigger.append(trig_checker(rawEvt))

        del usefulEvent

    #saving result
    hf.create_dataset('Trigger', data=np.asarray(trigger), compression="gzip", compression_opts=9)

    hf.close()
    del hf

    print('output is',Output+h5_file_name)
    del Output, h5_file_name

    print('Indi SNR making is done!')

def indi_snr_maker_debug(Station, Run, Output, sel_evt # argv
                    , R, evtTree, rawEvt, num_evts, cal, q # ara root
                    , num_Ants, bad_ant_i, t_width_ns # known config
                    , t_pad_l, time_pad_i, time_pad_f # time
                    , f, f_w # freq
                    , n_psd # psd
                    , temp_v, n_theta): # temp

    print('Indi SNR making starts!')

    import h5py

    # create output file
    os.chdir(Output)
    h5_file_name='Indi_SNR_A'+str(Station)+'_R'+str(Run)+'_debug.h5'
    hf = h5py.File(h5_file_name, 'w')

    g2 = hf.create_group('SNR')

    trigger = []

    # loop over the events
    for event in tqdm(range(num_evts)):

        # make a useful event
        usefulEvent = useful_evt_maker(R, evtTree, rawEvt, event, cal)

        # trigger filter
        if rawEvt.isSoftwareTrigger() == 0:

            # quality cut
            if q.isGoodEvent(usefulEvent) == 1:

                # make padded wf with interpolated wf to fft
                ant_arr = station_pad(usefulEvent, num_Ants, t_width_ns, t_pad_l, time_pad_i, time_pad_f)

                if event == sel_evt:

                    trig_index = trig_checker(rawEvt)
                    if trig_index == 0:
                        trig_type = 'RF'
                    elif trig_index == 1:
                        trig_type = 'Cal'

                    print('Evt#'+str(event)+' is selected!')
                    print('Trigger type is',trig_type)

                    # wf plot
                    from tools.plot import plot_16
                    ant_arr_copy = np.copy(ant_arr)

                    plot_16(r'Time [ $ns$ ]',r'Amplitude [ $V$ ]',trig_type+' WF, A'+str(Station)+', Run'+str(Run)+', Evt'+str(event)
                                ,np.arange(time_pad_i, time_pad_f+t_width_ns, t_width_ns),ant_arr_copy
                                ,np.round(np.nanmax(np.abs(ant_arr_copy),axis=0),2)
                                ,time_pad_i,time_pad_f
                                ,Output,trig_type+'_WF_A'+str(Station)+'_Run'+str(Run)+'_Evt'+str(event)+'.png'
                                ,'Event# '+str(event)+' '+trig_type+' WF plot was generated!')

                    # fft plot
                    ant_arr_fft = np.fft.fft(ant_arr_copy, axis=0)
                    ant_arr_fft_band = Band_Square_debug(f, np.repeat(ant_arr_fft[:,:,np.newaxis], n_theta, axis=2))

                    from tools.plot import plot_16_log_theta

                    plot_16_log_theta(r'Frequency [ $GHz$ ]',r'Amplitude [ $V$ ]', trig_type+' FFT, A'+str(Station)+', Run'+str(Run)+', Evt'+str(event)
                                ,f/1e9,ant_arr_fft
                                ,f/1e9,ant_arr_fft_band[:,:,0]
                                ,f/1e9,ant_arr_fft_band[:,:,1]
                                ,f/1e9,ant_arr_fft_band[:,:,2]
                                ,1e-4,1e2
                                ,Output,trig_type+'_FFT_A'+str(Station)+'_Run'+str(Run)+'_Evt'+str(event)+'.png'
                                ,'Event# '+str(event)+' '+trig_type+' FFT plot was generated!')

                else:
                    pass

                # OMF
                snr_wf = OMF(n_psd, f_w, Band_Square(f, np.repeat(np.fft.fft(ant_arr, axis=0)[:,:,np.newaxis], n_theta, axis=2)), temp_v)
                del ant_arr

                if event == sel_evt:

                    #snr plot
                    from tools.plot import plot_16_3
                    snr_wf_copy = np.copy(snr_wf)

                    plot_16_3(r'Offset Time [ $ns$ ]',r'SNR [ $V/RMS$ ]', trig_type+' SNR, A'+str(Station)+', Run'+str(Run)+', Evt'+str(event)+', w/ tale'
                                ,np.arange(time_pad_i, time_pad_f+t_width_ns, t_width_ns),snr_wf_copy[:,:,0],np.round(np.nanmax(snr_wf_copy[:,:,0],axis=0),2)
                                ,np.arange(time_pad_i, time_pad_f+t_width_ns, t_width_ns),snr_wf_copy[:,:,1],np.round(np.nanmax(snr_wf_copy[:,:,1],axis=0),2)
                                ,np.arange(time_pad_i, time_pad_f+t_width_ns, t_width_ns),snr_wf_copy[:,:,2],np.round(np.nanmax(snr_wf_copy[:,:,2],axis=0),2)
                                ,Output,trig_type+'_SNR_A'+str(Station)+'_Run'+str(Run)+'_Evt'+str(event)+'_w_tale.png'
                                ,'Event# '+str(event)+' '+trig_type+' SNR w/ tale plot was generated!')

                else:
                    pass

                # remove bad antenna
                snr_wf[:, bad_ant_i, :] = np.nan

                #saving result
                g2.create_dataset(str(event), data=snr_wf, compression="gzip", compression_opts=9)
                del snr_wf

                # trigger check
                trigger.append(trig_checker(rawEvt))

        del usefulEvent

    #saving result
    g2.create_dataset('Trigger', data=np.asarray(trigger), compression="gzip", compression_opts=9)
    del g2

    g3 = hf.create_group('SNR_indi')
    g3.create_dataset('Trig_index', data=np.array([trig_index]), compression="gzip", compression_opts=9)
    g3.create_dataset('Freq', data=f, compression="gzip", compression_opts=9)
    g3.create_dataset('Time', data=np.arange(time_pad_i, time_pad_f+t_width_ns, t_width_ns), compression="gzip", compression_opts=9)
    g3.create_dataset('WF_evt'+str(sel_evt), data=ant_arr_copy, compression="gzip", compression_opts=9)
    g3.create_dataset('FFT_evt'+str(sel_evt), data=ant_arr_fft, compression="gzip", compression_opts=9)
    g3.create_dataset('FFT_band_evt'+str(sel_evt), data=ant_arr_fft_band, compression="gzip", compression_opts=9)
    g3.create_dataset('SNR_evt'+str(sel_evt), data=snr_wf_copy, compression="gzip", compression_opts=9)    
    del g3, trig_index, ant_arr_copy, ant_arr_fft, ant_arr_fft_band, snr_wf_copy, trig_type

    hf.close()
    del hf

    print('output is',Output+h5_file_name)
    del Output, h5_file_name

    print('Indi SNR making is done!')
"""
"""
def psd_maker(vol, dt, i_time_len):

    vol = np.fft.fft(vol, axis=0)
    vol *= dt
    vol *= vol.conjugate()
    vol *= 2
    vol /= dt
    vol /= i_time_len[np.newaxis, :]

    return vol
"""
"""
def OMF(psd, rf_f_w, rf_v, temp_v):

    # matched filtering
    conjugation = temp_v * rf_v.conjugate()
    conjugation /= psd
    conjugation = np.abs(hilbert(2*np.fft.ifft(conjugation,axis=0).real, axis = 0))
    #conjugation = np.abs(2*np.fft.ifft(conjugation,axis=0))
    conjugation /= np.sqrt(np.abs(2*np.nansum(temp_v * temp_v.conjugate() / psd,axis=0) * rf_f_w))[np.newaxis, :, :]
    conjugation = conjugation[::-1,:,:]
    conjugation[conjugation<0] = 0

    return conjugation
"""

