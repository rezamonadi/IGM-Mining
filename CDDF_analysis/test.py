import calc_cddf
import os.path as path
import numpy as np
import h5py
import matplotlib
from scipy.io import loadmat
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import calc_cddf
import make_plots

calc_cddf.compute_all_snrs(raw_file="../data/dr12/processed/preloaded_qsos_C13_full.mat", processed_file="pltData.mat", save_file="snrs_qsos_dr12.mat")
# ff = h5py.File('snrs_qsos_dr12.mat', 'r')
print(ff.keys())
print(ff['snrs'][:100])