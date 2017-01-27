# This Python file uses the following encoding: utf-8
import sys

from scipy.io import netcdf
import numpy as np
import math
import subprocess
import colormaps as cmaps


import h5py as h5

sys.path.insert(0, "../../../../../build/bindings/python/")
from libcloudphxx import common as cm

sys.path.insert(0, "../../../../../../parcel/")
import functions as fn

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt


for case in ('case_base', 'case3'):

    data = h5.File(case + '/out_hall_pinsky_stratocumulus/timestep0000011800.h5', 'r')
    mesh = h5.File(case + '/out_hall_pinsky_stratocumulus/const.h5', 'r')

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

    fig = plt.figure()
    ax=fig.add_subplot(111)

    cmap = plt.get_cmap('coolwarm')
    cmap.set_under('LightGray')

    #vmin = 1
    #vmax = 5

    print case
    print pH.max()
    print pH.min()

    cplt = ax.pcolormesh(x_grid, y_grid, pH, cmap = cmap, vmin=3.8, vmax=5)#cmap=cmaps.viridis)
    cbar = fig.colorbar(cplt)

    ax.set_xlabel('x'); 
    ax.set_ylabel('y'); 
    ax.set_title('water weighted average pH')
    ax.set_xlim([0, 75])
    ax.set_ylim([0, 75])

    plt.savefig(case + "_ph_profile.pdf")

#plot settings
#fig = plt.figure(figsize=(28,13))
#plt.rcParams.update({'font.size': 30})

#ax.set_yticks([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600])
#ax.set_yticklabels(y_labels2)
#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")
#ax.set_xlim([3.6, 5])
#ax.set_xticks([3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5])

#plt.grid()
#ax.set_xlabel('average pH')
#ax.set_ylabel('height above the ground [km]')
#plt.plot(pH, t, "b", lw=4.)
