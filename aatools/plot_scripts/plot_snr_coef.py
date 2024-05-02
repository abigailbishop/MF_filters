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
import numpy as np
import h5py
import ROOT
from ROOT import TCanvas, TH1D
import matplotlib.pyplot as plt
from tqdm import tqdm as tq
import re

def main(file_in, targetFile, plot_title):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(1)

    hist_1d = TH1D("SNR (1D)", plot_title, 100, 0, 15)
    hist_1d.GetXaxis().SetTitle("SNR")
    hist_1d.GetYaxis().SetTitle("Counts")

    # Set fill color to blue
    hist_1d.SetFillColor(ROOT.kBlue)

    # Set title font weight to bold
    hist_1d.SetTitleFont(43, "t")

    total_runs = sum(1 for _ in open(file_in))
    pbar = tq(total=total_runs)
    with open(file_in, 'r') as file_open:
        for lines in file_open.readlines():
            run_number = int(re.search(r'R(\d+)', lines).group(1))
            with h5py.File(lines.strip(), 'r') as hf:
                freq_range = np.array(hf.get('snr')[0,:])
                trig = np.array(hf.get('trig_type'))
                for i in range(len(trig)):
                    #if trig[i] == 1:
                    hist_1d.Fill(freq_range[i])
            pbar.update()

    c1 = TCanvas("c1", "c1", 40, 10, 1400, 900)
    hist_1d.Draw()
    c1.Print(targetFile, 'png')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <file_in> <targetFile> <plot_title>")
        sys.exit(1)
    file_in = sys.argv[1]
    targetFile = sys.argv[2]
    plot_title = sys.argv[3]
    main(file_in, targetFile, plot_title)

