import h5py
import numpy as np
import os
import matplotlib.pyplot as plt

def get_skymap(
    station, run, event,
    use_l2 = False, condor_run = False, blind_dat = True, no_tqdm = False,
    file_type = 'event', verbose = False, return_none = False, return_dat_only = False
):
    """
    For a given station, run, and event, use the MF_filters framework to return
      the skymap of correlation values. 

    Parameters
    ----------
    station : int
        Station number
    run : int
        Run number
    event : int
        Event number

    Returns
    -------
    sky_map : numpy.nd_array
        The skymap of correlation coefficients with shape (antenna polarization,
        zenith angle, azimuth angle, reconstruction radius, ray solution)
    """

    if use_l2:
        from tools.ara_data_load import ara_l2_loader
    else:
        from tools.ara_data_load import ara_uproot_loader
        from tools.ara_data_load import ara_root_loader
    from tools.ara_constant import ara_const
    from tools.ara_wf_analyzer import wf_analyzer
    from tools.ara_py_interferometers_ele import py_interferometers
    from tools.ara_py_interferometers_ele import get_products
    from tools.ara_run_manager import run_info_loader, condor_info_loader
    from tools.ara_known_issue import known_issue_loader
    from tools.ara_quality_cut import get_bad_events

    # geom. info.
    ara_const = ara_const()
    num_ants = ara_const.USEFUL_CHAN_PER_STATION
    num_pols = ara_const.POLARIZATION
    del ara_const

    run_info = run_info_loader(station, run, analyze_blind_dat = blind_dat)
    Data, Ped = run_info.get_data_ped_path(
        file_type = file_type, return_none = return_none,
        verbose = verbose, return_dat_only = return_dat_only, l2_data = use_l2)
    station, run, Config, Year, Month, Date = run_info.get_data_info()

    condor_info = condor_info_loader(use_condor = condor_run, verbose = True)
    Data = condor_info.get_target_to_condor_path(Data)
    Ped = condor_info.get_target_to_condor_path(Ped)

    # data config
    if use_l2:
        ara_root = ara_l2_loader(Data)
        num_evts = ara_root.num_evts
        evt_num = ara_root.evt_num
        daq_qual_cut_sum = ara_root.daq_cut
        trig_type = ara_root.trig_type
        st = ara_root.station_id
        run = ara_root.run
    else:
        ara_uproot = ara_uproot_loader(Data)
        evt_num = ara_uproot.evt_num
        num_evts = ara_uproot.num_evts
        trig_type = ara_uproot.get_trig_type()
        st = ara_uproot.station_id
        yr = ara_uproot.year
        run = ara_uproot.run
        ara_root = ara_root_loader(Data, Ped, st, yr)
        del ara_uproot, yr

    # pre quality cut
    if use_l2 == False:
        daq_qual_cut_sum = get_bad_events(
            st, run, 
            analyze_blind_dat = blind_dat, 
            verbose = True, 
            evt_num = evt_num, 
            qual_type = 2
        )[0]

    known_issue = known_issue_loader(st)
    bad_ant = known_issue.get_bad_antenna(run, print_integer = True)
    del known_issue

    # sub info
    run_info = run_info_loader(st, run, analyze_blind_dat = blind_dat)
    wei_dat = run_info.get_result_path(file_type = 'snr', verbose = True)
    wei_hf = h5py.File(wei_dat, 'r')
    weights = wei_hf['snr'][:]
    evt_num_b = np.full((10), -1, dtype = int)
    if blind_dat:
        print('BURN!!!!!')
        reco_dat = run_info.get_result_path(
            file_type = 'reco_ele_lite', 
            verbose = True, return_none = True, force_unblind = True
        )
    del wei_dat, wei_hf, run_info

    # wf analyzer
    wf_int = wf_analyzer(
        use_time_pad = True, use_band_pass = True, use_cw = True, 
        verbose = True, use_l2 = use_l2, analyze_blind_dat = blind_dat, 
        st = st, run = run
    )

    # interferometers
    ara_int = py_interferometers(
        wf_int.pad_len, wf_int.dt, 
        st, run = run, 
        get_sub_file = True, verbose = True, use_only_max = blind_dat
    )
    only_max = ara_int.use_only_max
    num_angs = ara_int.num_angs
    radius = ara_int.radius
    theta = ara_int.theta
    re_shape = ara_int.results_shape
    wei_pairs = get_products(weights, ara_int.pairs, ara_int.v_pairs_len)
    del st, run, weights

    # output array 
    coef = np.full(
        (re_shape[0], re_shape[1], re_shape[2], re_shape[3], num_evts), 
        np.nan, dtype = float
    ) # pol, theta, rad, ray, evt
    coord = np.copy(coef) # pol, theta, rad, ray, evt
    del re_shape  

    # get entry and wf
    ara_root.get_entry(event)
    ara_root.get_useful_evt(ara_root.cal_type.kLatestCalibWithOutTrimFirstBlock)

    # loop over the antennas
    for ant in range(num_ants):
        raw_t, raw_v = ara_root.get_rf_ch_wf(ant)
        wf_int.get_int_wf(
            raw_t, raw_v, ant, 
            use_zero_pad = True, use_band_pass = True, 
            evt = event
        )
        del raw_t, raw_v
        ara_root.del_TGraph()
    ara_root.del_usefulEvt()   

    sky_map = ara_int.get_sky_map(wf_int.pad_v, weights = wei_pairs[:, event], return_skymap=True)

    return sky_map


def plot_skymap(
    sky_map, radius_index, 
    polarization=0, solution=0, norm=None, 
    cmap='BuPu', cbar_min=None, cbar_max=None,
    title=None, save_name=None
):
    """
    Given a skymap of correlation values, plot the skymaps for the requested 
      reconstruction radius, antenna polarization, and ray solution.
    """

    fig, ax = plt.subplots(figsize=(8,4))
    plotter = sky_map[polarization, :, :, radius_index, solution]
    extent = [-179.5, 179.5, -89.5, 89.5]
    mappable = ax.imshow(
        plotter, extent=extent, cmap=cmap, 
        origin='upper', aspect='auto', norm=norm)
    
    # Add Labels
    ax.set_xlabel("Azimuth [deg]")
    ax.set_ylabel("Zenith [deg]")
    if title==None: 
        ax.set_title(f"Polarization {polarization}, Solution {solution}")
    else: 
        ax.set_title(title)

    # Build colorbar
    if cbar_min !=None: mappable.norm.vmin = cbar_min
    if cbar_max !=None: mappable.norm.vmax = cbar_max
    plt.colorbar(mappable=mappable, label="Correlation Coefficient")

    # Cleanup and save plot (if requested)
    plt.tight_layout()
    if save_name!=None: plt.savefig(save_name, dpi=400)

    return fig, ax


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('station', type=int)
    parser.add_argument('run', type=int)
    parser.add_argument('event', type=int)
    parser.add_argument('radius_index', type=int)
    parser.add_argument('-p', '--plotdir', type=str, default=None)
    parser.add_argument('-d', '--datadir', type=str, default="reco_ele_lite")
    parser.add_argument('-t', '--title',   type=str, default=None)
    parser.add_argument('--polarization', type=int, default=0)
    parser.add_argument('--solution', type=int, default=0)
    parser.add_argument('--cbar_norm', type=str, default=None)
    parser.add_argument('--cbar_min', type=float, default=None)
    parser.add_argument('--cbar_max', type=float, default=None)
    args = parser.parse_args()

    if args.plotdir == None: 
        plotdir = f"/data/ana/ARA/ARA0{args.station}/plots/reco_maps/"
    else:
        plotdir = args.plotdir
    if not os.path.exists(plotdir):
        raise ValueError(f"Requested plot directory does not exist: {plotdir}")
    save_name = f"{plotdir}/skymap_A{args.station}_R{args.run}_E{args.event}.png"
    print(f"Will save plot to: {save_name}")

    if args.title == None:
        title = f"ARA0{args.station}, Run {args.run}, Event {args.event}, Polarization {args.polarization}, Solution {args.solution}"
    else: 
        title = args.title

    # Load and plot the reco file
    # reco_file = h5py.File(reco_file_path)
    # plot_event_zen_phi(
    #     reco_file, args.event, args.radius_index,
    #     pol=args.polarization, sol=args.solution, title=title,
    #     save_name=save_name,
    #     norm=args.cbar_norm, cbar_min=args.cbar_min, cbar_max=args.cbar_max
    # )
    # reco_file.close()

    plot_skymap(
        get_skymap(args.station, args.run, args.event), args.radius_index, 
        save_name=save_name, title=title,
        polarization=args.polarization, solution=args.solution,
        norm=args.cbar_norm, cbar_min=args.cbar_min, cbar_max=args.cbar_max
    )

    return


if __name__ == "__main__":
    main()