import numpy as np
import h5py
import ROOT
from ROOT import TCanvas, TH1D
import matplotlib.pyplot as plt
from tqdm import tqdm as tq
import re

def draw_horizontal_line(ax, x_min, x_max, y):
    # Create an array of x values for the line
    x_values = np.linspace(x_min, x_max, 100)
    # Create an array of y values for the line (all equal to y)
    y_values = np.full_like(x_values, y)
    # Plot the horizontal line
    ax.plot(x_values, y_values, 'r-', linewidth=1)

def main():
    file_in = "/data/user/abishop/ara/a23/a2/step2/ARA02/snr_full/snr_paths.txt"
    targetFile = "plot_1d_snr.png"
#    file_in1 = "/data/ana/ARA/ARA04/reco/paths_reco.txt"
#    targetFile1 = "plot_1d_reco.png"
    fileType = "png"
    c1 = TCanvas("c1", "c1", 40, 10, 1400, 900)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(1)

    hist_1d = TH1D("SNR (1D)", "SNR (1D)", 100, 0, 15)
    hist_1d.GetXaxis().SetTitle("SNR")
    hist_1d.GetYaxis().SetTitle("Counts")
    # Create a 1D histogram
    # Create a 1D histogram

    total_runs = sum(1 for _ in open(file_in))
    pbar = tq(total=total_runs)
    file_open = open(file_in, 'r')
    for lines in file_open.readlines():
        run_number = int(re.search(r'R(\d+)', lines).group(1))
        hf = h5py.File(str(lines.strip()), 'r')
        freq_range = np.array(hf.get('snr')[0,:])
        trig = np.array(hf.get('trig_type'))
        for i in range(len(trig)):
            #if trig[i] == 1:
               hist_1d.Fill(freq_range[i])
        pbar.update()
    hist_1d.Draw()
    plt.axvline(x=0.9, color='red', linewidth=1)
    c1.Print(targetFile, fileType)
    file_open.close()
#    c2 = TCanvas("c2", "c2", 40, 10, 1400, 900)
#    hist = TH1D("Correlation (1D)", "Correlation (1D)", 100, 0, 1)
#    hist.GetXaxis().SetTitle("Corr coef")
#    hist.GetYaxis().SetTitle("Counts")

#    file_open1 = open(file_in1, 'r')
#    for lines1 in file_open1.readlines():
#        hf1 = h5py.File(str(lines1.strip()), 'r')
#        coef = np.array(hf1.get('coef')[0,0,0,:])
#        trig = np.array(hf1.get('trig_type'))
#        for i in range(len(trig)):
#            if trig[i] == 1:
#               hist.Fill(coef[i])
#        pbar.update() 

    # Create the plot for the 1D histogram
#    hist.Draw()
#    c2.Print(targetFile1,fileType)
#    file_open1.close()     
#    pbar.close()

if __name__ == "__main__":
    main()

