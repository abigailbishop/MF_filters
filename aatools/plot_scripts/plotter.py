"""
Script: plot_snr_coef.py
Authors: Alan Salcedo Gomez, Abigail Bishop

Description:
This script reads data from HDF5 files with paths specified in a text file, fills a histogram
with the SNR (Signal-to-Noise Ratio) data, and plots the histogram using ROOT. 
It saves the plot as an image.

Usage:
python plot_snr_coef.py <file_in> <targetFile> <plot_title>
- file_in: Path to the text file containing paths to HDF5 files.
- targetFile: Name of the output plot file.
- plot_title: Title of the plot.
Example: python plot_snr_coef.py snr_paths_ARA02.txt plot_1d_snr_all_A2.png "A2 All Triggers (Oct. 2020 - May 2021)"
"""

import sys
import os
import numpy as np
import h5py
import ROOT
from ROOT import TCanvas, TH1D
import matplotlib.pyplot as plt
from tqdm import tqdm as tq
import re

def int_to_shorthand(number): 
    """
    Take in an integer, and return its shorthand as a string. Ex: 1 -> 1st
    """
    number = int(number)
    if number == 1: 
        return "1st"
    elif number == 2: 
        return "2nd"
    elif number == 3: 
        return "3rd"
    else: 
        return f"{number}th"
    
def trig_type_str(trig_type, shorthand=False):
    """
    From the integer version of trigger type, return string representation of
      trigger type
    """
    trig_type = int(trig_type)
    if trig_type == 0: 
        if shorthand: return "RF"
        else: return "Normal Triggers"
    elif trig_type == 1:
        if shorthand: return "CP"
        else: return "Calibration Pulsers"
    elif trig_type == 2:
        if shorthand: return "Soft"
        else: return "Software Triggers"
    else: 
        raise ValueError(f"Provided trigger type of {trig_type} is invalid.")

def get_stats_text(data):
    mean_value = np.nanmean(data)
    std_deviation = np.nanstd(data)
    stats_text = f"Mean: {mean_value:.2f}\nStd Dev: {std_deviation:.2f}"
    stats_text += f"\nEvents: {np.count_nonzero(~np.isnan(data))}"
    return stats_text

def Distribution_1D(file_in, targetFile, plot_title, ana_variable, trigger_type):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(1)

    hist_1d = TH1D(f"{ana_variable} (1D)", plot_title, 100, 0, 15)
    hist_1d.GetXaxis().SetTitle(f"{ana_variable}")
    hist_1d.GetYaxis().SetTitle("Counts")

    # Set fill color to blue
    hist_1d.SetFillColor(ROOT.kBlue)

    # Set title font weight to bold
    hist_1d.SetTitleFont(43, "t")

    total_runs = len(file_in)
    pbar = tq(total=total_runs)
    for lines in file_in:
        if not os.path.exists(lines.strip()):
            print(f"Skipping missing file {lines}")
            continue
        run_number = int(re.search(r'R(\d+)', lines).group(1))
        with h5py.File(lines.strip(), 'r') as hf:
            freq_range = np.array(hf.get(f'{ana_variable}')[0,:])
            trig = np.array(hf.get('trig_type'))
            for i in range(len(trig)):
                if trig[i] == trigger_type:
                    hist_1d.Fill(freq_range[i])
        pbar.update()

    c1 = TCanvas("c1", "c1", 40, 10, 1400, 900)
    hist_1d.Draw()
    c1.Print(targetFile, 'png')

def Distribution_1D_plt(
    file_in, targetFile, 
    plot_title, ana_variable, trigger_type
):
    total_runs = len(file_in)
    max_entries_per_run = 150000
    data_to_plot = np.full(total_runs*max_entries_per_run, np.nan)
    current_index = 0
    for lines in file_in:
        if not os.path.exists(lines.strip()):
            print(f"Skipping missing file {lines}")
            continue
        file = h5py.File(lines.strip(), "r")
        data = np.nanmax( 
            np.array(file.get(f'{ana_variable}')),
            axis=0
        )
        trig = np.array(file.get('trig_type'))
        for i in range(len(trig)):
            if trig[i] == trigger_type:
                data_to_plot[current_index] = data[i]
                current_index += 1
        file.close()
        del file, data

    fig, ax = plt.subplots()
    ax.hist(data_to_plot[:current_index], bins=20)
    ax.set_xlabel(f"{ana_variable}")
    ax.set_ylabel("Counts")
    ax.set_title(plot_title)

    # Calculate statistics and put in a box on the plot
    stats_text = get_stats_text(data_to_plot)
    plt.text(0.95, 0.95, stats_text, ha='right', va='top', 
             transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', alpha=0.5, edgecolor='black'))
    
    plt.tight_layout()
    plt.savefig(targetFile)
    
    return fig, ax

def plot_ant_stats(
    files,  save_name, 
    plot_title, analysis_variable, trigger_type, 
    index_from_highest = 3, 
):
    total_runs = len(files)
    max_entries_per_run = 150000
    data_to_plot = np.full(total_runs*max_entries_per_run, np.nan)
    current_index = 0
    for lines in files:

        if not os.path.exists(lines.strip()):
            print(f"Skipping missing file {lines}")
            continue

        file = h5py.File(lines.strip(), "r")

        # Get third highest RPR for each file
        data = np.array(file['snr'])
        data = np.sort(data, axis=0)[-index_from_highest]

        # Check trigger type and save to array
        trig = np.array(file.get('trig_type'))
        print("data", data.shape, "trig", trig.shape)
        for i in range(len(trig)):
            if trig[i] == trigger_type:
                data_to_plot[current_index] = data[i]
                current_index += 1
        file.close()
        del file, data

    if len(files) == 1: 
        xlabel  = f"{int_to_shorthand(index_from_highest)} Max {analysis_variable}"
        xlabel += f" in Run {get_run_number_from_file(files[0])}"
    else: 
        xlabel = f"{int_to_shorthand(index_from_highest)} Max {analysis_variable}"

    plot_title = f"{plot_title} - {trig_type_str(trigger_type)}"

    fig, ax = plt.subplots()
    ax.hist(data_to_plot[:current_index], bins=20)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")
    ax.set_title(plot_title)

    # Calculate statistics and put in a box on the plot
    stats_text = get_stats_text(data_to_plot)
    plt.text(0.95, 0.95, stats_text, ha='right', va='top', 
             transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', alpha=0.5, edgecolor='black'))

    plt.tight_layout()
    plt.savefig(save_name)
    
    return fig, ax

def hist2d(
    x_array, y_array, 
    cmap="ocean_r", figsize=(5,4),
    bins=10, xyranges=None, 
    weights=None, density=None, norm=None,
    x_label=None, y_label=None, title=None,
    cbar_min=None, cbar_max=None,
    save_name=None,
):
    """
    Plots a 2D histogram without all the classic clutter
    """
    
    fig, ax = plt.subplots(figsize=figsize)
    H, xbins, ybins = np.histogram2d(
        x_array, y_array, bins=bins, range=xyranges,
        weights=weights, density=density
    )
    extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
    mappable = ax.imshow(H.T, extent=extent, cmap=cmap, 
                         origin='lower', aspect='auto', norm=norm)
    
    if x_label  !=None: ax.set_xlabel(x_label)
    if y_label  !=None: ax.set_ylabel(y_label)
    if title    !=None: ax.set_title(title)
    if cbar_min !=None: mappable.norm.vmin = cbar_min
    if cbar_max !=None: mappable.norm.vmax = cbar_max
    if xyranges !=None: 
        if np.array(xyranges).ndim == 2: 
            ax.set_xlim(*xyranges[0]); ax.set_ylim(*xyranges[1])
        else: 
            ax.set_xlim(*xyranges); ax.set_ylim(*xyranges)
    
    if isinstance(weights, np.ndarray):
        plt.colorbar(mappable=mappable, label="Weighted Number of Events")
    elif isinstance(weights, list):
        plt.colorbar(mappable=mappable, label="Weighted Number of Events")
    else:
        plt.colorbar(mappable=mappable, label="Number of Events")
    plt.tight_layout()
    
    if save_name!=None: plt.savefig(save_name, dpi=300)
        
    return fig, ax

def plot_zen_phi(
    reco_files, trigger_type, save_name, cmap="BuPu",
    pol=0, sol=0, radius_index=0, bins=10, figsize=(5,4), plot_title=""
):

    # Dataset {2, 180, 5, 3, 7866}
    
    total_runs = len(reco_files)
    max_entries_per_run = 150000
    zeniths = np.full(total_runs*max_entries_per_run, np.nan)
    phis = np.full(total_runs*max_entries_per_run, np.nan)
    current_index = 0

    for file in reco_files:

        if not os.path.exists(file):
            print(f"Skipping missing file {file}")
            continue

        file = h5py.File(file.strip(), "r")

        # Currently not in use (didn't work lol)
        # Picked the event based on the reconstruction radius with the best result
        if radius_index == None: 
            # Only use the best radius
            coefs = file['coef'][pol, :, :, sol, :]
            coords = file['coord'][pol, :, :, sol, :]   
            coef_max_r_idx = np.nanargmax(coefs, axis=2)  
            coefs = np.array(coefs)[
                np.arange(coefs.shape[0])[:, np.newaxis, np.newaxis],
                coef_max_r_idx,
                np.arange(coefs.shape[2])[np.newaxis, np.newaxis, :],
            ]
            coords = np.array(coords)[
                np.arange(coords.shape[0])[:, np.newaxis, np.newaxis],
                coef_max_r_idx,
                np.arange(coords.shape[2])[np.newaxis, np.newaxis, :],
            ]
        else: 
            coefs = file['coef'][pol, :, radius_index, sol, :]
            coords = file['coord'][pol, :, radius_index, sol, :]

        coef_max_idx = np.nanargmax(coefs, axis=0)
        best_coefs = np.array(coefs)[
            coef_max_idx,
            np.arange(coefs.shape[1])[np.newaxis, :],
        ]

        best_thetas = np.array(file['theta'])[coef_max_idx]
        best_phis = np.array(coords)[
            coef_max_idx,
            np.arange(coords.shape[1])[np.newaxis, :],
        ]

        # Check trigger type and save to array
        trig = np.array(file.get('trig_type'))
        for i in range(len(trig)):
            if trig[i] == trigger_type:
                zeniths[current_index] = best_thetas[i]
                phis[current_index] = best_phis[0,i]
                current_index += 1
        file.close()
        del file

    if plot_title != "":
        plot_title = f"{plot_title} - {trig_type_str(trigger_type)}"
    else:
        plot_title = trig_type_str(trigger_type)

    fig, ax = hist2d(
        phis[:current_index], zeniths[:current_index], 
        bins=bins, save_name=save_name, figsize=figsize,
        cmap = "Purples",
        x_label="Reconstructed Azimuthal Angle [deg]",
        y_label="Reconstructed Elevation Angle [deg]",
        title=plot_title
    )

    return fig, ax

def plot_rad_zen_coef(
    reco_file, event, 
    pol=0, sol=0, save_name=None,
    cmap="Purples", cbar_min=None, cbar_max=None, 
):
    """"
    For an event index in the provided h5 reconstruction file, plot the 
      correlation coefficients calculated for the various radii and zeniths.

    Parameters
    ----------
    reco_file : h5py._hl.files.File
        The h5 file (already loaded into python with h5py.File()) with data
          to be plotted.
    event : int
        The index in the `reco_file` corresponding to the event to be plotted.
    pol : int
        The index of the polarization that was reconstructed.
    sol : int
        The index of the ray solution that was reconstructed.
    save_name : str
        Path to where you would like to save the created plot.
    cmap : str or matplotlib.colors.Colormap
        Colormap with which you'd like the heatmap to be plotted.
    cbar_min : float
        Value corresponding to the minimum value on the color bar.
    cbar_max : float
        Value corresponding to the maximum value on the colorbar.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure that was plotted.
    ax : matplotlib.axes._axes.Axes
        The axes of the figure that was plotted.
    """

    # Parse data from the provided reco_file
    radii = reco_file['radius']
    max_radius = max(radii)
    zeniths = reco_file['theta']
    coefs = reco_file['coef'][pol, :, :, sol, event]

    # Size (in meters) of each bin created along the radius axis
    radius_bin_width = 30

    # Initialize the array of values we'll be plotting 
    # Will plot data for all zeniths
    # Will plot radii from 0 to to max_radius + (3*radius_bin_width)
    plotter = np.full( 
        ( coefs.shape[0],                    # Number of zenith bins to plot
         (max_radius//radius_bin_width)+3 ), # Number of radius bins to plot 
         # Increment radius bins by 3 to put padding on the right of the plot
        np.nan # Initialize to values of `nan` which will appear blank on plot
    )

    # Add correlation coefficients to the plotting array
    for radius_index, radius in enumerate(radii):
        plotter_index = radius//radius_bin_width # ex: 170//20 = 8
        plotter[:,plotter_index] = coefs[:,radius_index]

    # Calculate horizontal and vertical extents of the plotting array
    # Shift radius extent by -radius_bin_width/2 to center the bins 
    extent = [
        -radius_bin_width/2, # Lower radius
        (radius_bin_width*plotter.shape[1])-radius_bin_width/2, # Upper radius 
        zeniths[-1], # Upper zenith
        zeniths[0], # Lower zenith
    ]
 
    # Plot the plotter array
    fig, ax = plt.subplots()   
    mappable = ax.imshow(plotter, extent=extent, cmap=cmap, 
                         origin="upper", aspect='auto') 
    
    # Add Labels
    ax.set_xlabel("Reconstruction Radii [m]")
    ax.set_ylabel("Reconstruction Zeniths [deg]")
    ax.set_title(f"Event {event}, Polarization {pol}, Solution {sol}")

    # Build colorbar
    if cbar_min !=None: mappable.norm.vmin = cbar_min
    if cbar_max !=None: mappable.norm.vmax = cbar_max
    plt.colorbar(mappable=mappable, label="Correlation Coefficient")

    # Add reconstruction radii to the plot
    for radius in radii:
        plot_position = (radius+radius_bin_width, -45)
        ax.annotate(f"Reconstruction Radius: {radius:.0f} m", 
                    plot_position, xytext=plot_position, alpha=0.25,
                    horizontalalignment="center", rotation="vertical")

    # Cleanup and save plot (if requested)
    plt.tight_layout()
    if save_name!=None: plt.savefig(save_name, dpi=300)
        
    return fig, ax

def get_run_number_from_file(file_name):
    split_file_name = file_name.split("_R")[1]
    run_number = ""
    for char in split_file_name: 
        if char.isnumeric():  
            run_number += char
        else: 
            return int(run_number)
    return int(run_number)

def get_filelist(data_dir, first_run, last_run):
    dir_contents = sorted(os.listdir(data_dir+"/"))
    filelist = []
    for this_file in dir_contents:
        run_number = get_run_number_from_file(this_file)
        if (run_number>=first_run and run_number<=last_run):
            filelist.append(data_dir+"/"+this_file) 
    return filelist

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plot_dir', type=str,
        help="Directory in which to store plots")
    parser.add_argument('data_dir', type=str,
        help="Directory with output files")
    parser.add_argument('plot_title', type=str,
        help="Plot title")
    parser.add_argument('-k', '--keys', type=str, nargs="*", required=True,
        help="key corresponding to the script called by script_executor.py")
    parser.add_argument('-a', '--ana_variables', type=str, nargs="*", required=True,
        help="key in the h5files corresponding to the data we want to plot")
    parser.add_argument('--trig_type', '-t', type=int,
        help="Trigger type for the events to plot")
    parser.add_argument('--run_range', type=int, nargs=2, required=True,
        help="Runs to plot will follow range(args.run_range[0], args.run_range[1]+1)")
    args = parser.parse_args()

    plot_title = args.plot_title

    for key, ana_variable in zip(args.keys, args.ana_variables):
         file_in = get_filelist(args.data_dir+"/"+key, args.run_range[0], args.run_range[1]+1)
         targetFile = args.plot_dir+"/"+key+"_1D_Distribution_"+str(args.trig_type)+"_"+ana_variable+".png"
         Distribution_1D(file_in, targetFile, plot_title, ana_variable, args.trig_type)

