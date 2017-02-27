"""
Plot liquid water weighted average pH (2D plot of the whole cloudy field)

"""
# This Python file uses the following encoding: utf-8
import sys
import subprocess
import os

import numpy as np
import math
import h5py as h5

# libcloud bindings to python (to have access to library constants) 
sys.path.insert(0, "../../../../../build/bindings/python/")
from libcloudphxx import common as cm

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

for case in ('case_base', 'case1a', 'case3', 'case4', 'case5', 'case6'):
    # path to build directory with data and plots
    dir_path = '../../../build/thesis_2017/thesis_chap_6/'

    # read in the data
    data = h5.File(dir_path + 'out_' + case + '/timestep0000011800.h5', 'r')
    mesh = h5.File(dir_path + 'out_' + case + '/const.h5', 'r')

    # left and right edges of bins for dry and wet radius
    wet_edges = np.fromiter( (1e-6 * 10**(-3 + i * .1) for i in xrange(56)), dtype="float")
    num_wet   = wet_edges.size - 1

    # calculate average pH
    # (weighted with volume of cloud droplets)
    tmp    = data["rw_rng000_mom3"][:]
    x_grid = mesh["X"][:-1,:-1]
    y_grid = mesh["Y"][:-1,:-1]

    nom    = np.zeros(tmp.shape)
    den    = np.zeros(tmp.shape)
    pH     = np.zeros(tmp.shape)

    for i in range(num_wet-2): # first two bins are for total cloud and rain water conc
        r3_name = "rw_rng" + str(i+2).zfill(3) + "_mom3"
        H_name =  "chem_H_rw_rng" + str(i).zfill(3) + "_mom0"

        tmp_r3  = data[r3_name][:]
        tmp_nH  = data[H_name][:] / cm.M_H

        idx = np.where(tmp_r3 > 0)

        nom[idx] += (tmp_nH[idx] / (4./3 * math.pi * tmp_r3[idx] * 1e3)) * tmp_r3[idx]
        den[idx] += tmp_r3[idx]                             # to liters

    idx = np.where(den > 0)
    pH[idx]  = -1 * np.log10(nom[idx] / den[idx])

    print case
    print pH.max()
    print pH.min()

    #plot settings
    fig = plt.figure(figsize=(17,13))
    plt.rcParams.update({'font.size': 30})

    ax=fig.add_subplot(111)
   
    cmap = plt.get_cmap('Set1')
    cmap.set_under('White')
    vmin = 3
    vmax = 5.8
    cplt = ax.pcolormesh(x_grid, y_grid, pH, cmap = cmap, vmin=vmin, vmax=vmax)     
    cbar = fig.colorbar(cplt)

    ax.set_xlabel('X [km]'); 
    ax.set_ylabel('Z [km]'); 
    ax.set_title('water weighted average pH')

    ax.set_xlim([0, 75])
    ax.set_ylim([0, 75])

    ax.set_xticks([0, 15, 30, 45, 60, 75])
    ax.set_yticks([0, 15, 30, 45, 60, 75])
    ax.set_xticklabels(["0", "0.3", "0.6", "0.9", "1.2", "1.5"])
    ax.set_yticklabels(["0", "0.3", "0.6", "0.9", "1.2", "1.5"])

    plt.grid()

    # save fig to eps
    # I'm not sure why, but Python tends to produce very heavy figures
    # saving them as eps and then converting to pdf solves this problem...
    output_dir = dir_path + '/plots_of_pH/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plt.savefig(output_dir + case + "_ph_profile.eps")
