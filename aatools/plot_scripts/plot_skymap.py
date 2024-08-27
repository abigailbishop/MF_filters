from plotter import plot_event_zen_phi
import argparse
import h5py
import os

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

# Look for the data directory and reco file
datadir = f"/data/ana/ARA/ARA0{args.station}/{args.datadir}"
if not os.path.exists(datadir):
    raise ValueError(f"Requested data directory doesn't exist: {datadir}")
reco_file_path = f"{datadir}/{args.datadir}_A{args.station}_R{args.run}.h5"
if not os.path.exists(reco_file_path):
    raise ValueError(f"Reco file doesn't exist for requested run: {reco_file_path}")

if args.plotdir == None: 
    plotdir = f"/data/ana/ARA/ARA0{args.station}/plots/reco_maps/"
else:
    plotdir = args.plotdir
if not os.path.exists(plotdir):
    raise ValueError(f"Requested plot directory does not exist: {plotdir}")

if args.title == None:
    title = f"ARA0{args.station}, Run {args.run}, Event {args.event}, Polarization {args.polarization}, Solution {args.solution}"
else: 
    title = args.title

# Load and plot the reco file
reco_file = h5py.File(reco_file_path)
plot_event_zen_phi(
    reco_file, args.event, args.radius_index,
    pol=args.polarization, sol=args.solution, title=title,
    save_name=f"{plotdir}/skymap_A{args.station}_R{args.run}_E{args.event}.png",
    norm=args.cbar_norm, cbar_min=args.cbar_min, cbar_max=args.cbar_max
)
reco_file.close()