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

    total_runs = sum(1 for _ in open(file_in))
    pbar = tq(total=total_runs)
    for lines in file_in:
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

