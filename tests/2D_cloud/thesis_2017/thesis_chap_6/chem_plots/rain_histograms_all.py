import numpy   as np
import h5py    as h5
import Gnuplot as gp
import math    as mt
import sys
import os

# path to build directory with data
dir_path = '../../../build/tests/chem_sandbox/'

# all the collision kernels
kernels = { "hall_pinsky_stratocumulus",
            "hall", 
            "hall_davis_no_waals", 
            "vohl_davis_no_waals", 
            "onishi_hall", 
            "onishi_hall_davis_no_waals"}
                                         #   dry, cld,  rin
cover = {"hall_pinsky_stratocumulus" :   [0.,   0., 0.],
         "hall":                         [0.,   0., 0.], 
         "hall_davis_no_waals":          [0.,   0., 0.],  
         "vohl_davis_no_waals" :         [0.,   0., 0.], 
         "onishi_hall" :                 [0.,   0., 0.], 
         "onishi_hall_davis_no_waals" :  [0.,   0., 0.]}

# only grid-cells with rain water mixing ratio greater than cutoff will be shown
cutoff_cld = 0.01
cutoff_rin = 0.01 #g/kg

# left and right edges of bins for dry and wet radius
dry_edges = np.fromiter( (1e-6 * 10**(-3 + i * .1) for i in xrange(41)), dtype="float")
wet_edges = np.fromiter( (1e-6 * 10**(-3 + i * .1) for i in xrange(56)), dtype="float")
num_dry   = dry_edges.size - 1
num_wet   = wet_edges.size - 1

for kernel in kernels:

    # helper arrays for storing initial and cutoff dry radii as well as final and cutoff wet radii
    n_dry_ini = np.zeros(num_dry)
    n_dry_dry = np.zeros(num_dry)
    n_dry_cld = np.zeros(num_dry)
    n_dry_rin = np.zeros(num_dry)

    n_wet_all = np.zeros(num_wet)
    n_wet_dry = np.zeros(num_wet)
    n_wet_cld = np.zeros(num_wet)
    n_wet_rin = np.zeros(num_wet)

    n_run = 0.

    for run in {"44", "9", "13", "6", "42"}:
        # open hdf5 files with data
        h5f_ini       = h5.File(dir_path + 'out_case2_' + kernel + '_' + run + '/timestep0000000000.h5', 'r')
        h5f_fin       = h5.File(dir_path + 'out_case2_' + kernel + '_' + run + '/timestep0000011800.h5', 'r')

        # choose indexes of grid-cell with r_c and r_r > cutoff
        cld_mixr_fin = h5f_fin['rw_rng000_mom3'][:] * 4./3 * 3.14 * 1e3 * 1e3 #g/kg
        rin_mixr_fin = h5f_fin['rw_rng001_mom3'][:] * 4./3 * 3.14 * 1e3 * 1e3 #g/kg
        
        idx_dry = np.where(cld_mixr_fin <= cutoff_cld)
        idx_cld = np.where(cld_mixr_fin >  cutoff_cld)
        idx_rin = np.where(rin_mixr_fin >= cutoff_rin)

        # count dry cloudy and rainy grid cells
        cover[kernel][0] += idx_dry[0].size
        cover[kernel][1] += idx_cld[0].size
        cover[kernel][2] += idx_rin[0].size

        n_run += 1.
   
        for i in range(num_dry-1): # first bin id for total conc
            name = "rd_rng" + str(i+1).zfill(3) + "_mom0"
            tmp_ini  = 1e-6 * h5f_ini[name][:]
            tmp_fin  = 1e-6 * h5f_fin[name][:]

            n_dry_ini[i] += tmp_ini.sum()          / 76. / 76.       / (mt.log(dry_edges[i+1], 10) - mt.log(dry_edges[i], 10))
   
            n_dry_dry[i] += tmp_fin[idx_dry].sum() / idx_dry[0].size / (mt.log(dry_edges[i+1], 10) - mt.log(dry_edges[i], 10))
            n_dry_cld[i] += tmp_fin[idx_cld].sum() / idx_cld[0].size / (mt.log(dry_edges[i+1], 10) - mt.log(dry_edges[i], 10))
            n_dry_rin[i] += tmp_fin[idx_rin].sum() / idx_rin[0].size / (mt.log(dry_edges[i+1], 10) - mt.log(dry_edges[i], 10))
   
        for i in range(num_wet-2): # first two bins are for total cloud and rain water conc
            name = "rw_rng" + str(i+2).zfill(3) + "_mom0"
            tmp_fin  = 1e-6 * h5f_fin[name][:]

            n_wet_all[i]  += tmp_fin.sum()         / 76. / 76.       / (mt.log(wet_edges[i+1], 10) - mt.log(wet_edges[i], 10))

            n_wet_dry[i] += tmp_fin[idx_dry].sum() / idx_dry[0].size / (mt.log(wet_edges[i+1], 10) - mt.log(wet_edges[i], 10))
            n_wet_cld[i] += tmp_fin[idx_cld].sum() / idx_cld[0].size / (mt.log(wet_edges[i+1], 10) - mt.log(wet_edges[i], 10))
            n_wet_rin[i] += tmp_fin[idx_rin].sum() / idx_rin[0].size / (mt.log(wet_edges[i+1], 10) - mt.log(wet_edges[i], 10))
    
        # close hdf5 files
        h5f_ini.close()
        h5f_fin.close()

    # average over all runs
    n_dry_ini  /= n_run
    n_dry_dry  /= n_run
    n_dry_cld  /= n_run
    n_dry_rin  /= n_run

    n_wet_all  /= n_run
    n_wet_dry  /= n_run
    n_wet_cld  /= n_run
    n_wet_rin  /= n_run

    cover[kernel][0] = cover[kernel][0] / n_run / 76 / 76 * 100 
    cover[kernel][1] = cover[kernel][1] / n_run / 76 / 76 * 100 
    cover[kernel][2] = cover[kernel][2] / n_run / 76 / 76 * 100 

    print kernel, "  ", cover[kernel][0]
    print kernel, "  ", cover[kernel][1]
    print kernel, "  ", cover[kernel][2]

    # create directory for output (if it doesn't exist)
    output_dir = dir_path + 'plots_of_size_distr/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # plotting
    ymin = 0.001
    ymax = 1000
    xmin = 0.001
    xmax = 500

    g = gp.Gnuplot()
    g('reset')
    g('set term svg dynamic enhanced font "Verdana, 14"')
    g('set output "' + output_dir + 'case2_' + kernel + '_size_distr.svg" ')
    g('set logscale xy')
    g('set key samplen 1.2')
    g('set xtics rotate by 65 right (.01, .1, 1, 10, 100)')
    g('set ytics (.001, .01, .1, 1, 10, 100, 1000)')
    g('set xlabel "particle radius [um]"')
    g('set ylabel "dN/dlog_{10}(r) [mg^{-1} per log_{10}(size interval)]"')
    g('set xrange [' +  str(xmin) + ':' + str(xmax) + ']')
    g('set yrange [' +  str(ymin) + ':' + str(ymax) + ']')
    g('set grid')
    g('set nokey')
 
    plot_rd_ini  = gp.PlotItems.Data(dry_edges[:-1] * 1e6 , n_dry_ini,  with_="steps lw 4 lc rgb 'black'")
    #plot_rd_dry  = gp.PlotItems.Data(dry_edges[:-1] * 1e6 , n_dry_dry,  with_="steps lw 4 lc rgb 'green'")
    plot_rd_cld  = gp.PlotItems.Data(dry_edges[:-1] * 1e6 , n_dry_cld,  with_="steps lw 4 lc rgb 'green'")
    plot_rd_rin  = gp.PlotItems.Data(dry_edges[:-1] * 1e6 , n_dry_rin,  with_="steps lw 4 lc rgb 'red'")

    #plot_rw_all  = gp.PlotItems.Data(wet_edges[:-1] * 1e6 , n_wet_all,  with_="steps lw 4 lc rgb 'purple'")
    #plot_rw_dry  = gp.PlotItems.Data(wet_edges[:-1] * 1e6 , n_wet_dry,  with_="steps lw 4 lc rgb 'brown'")
    plot_rw_cld  = gp.PlotItems.Data(wet_edges[:-1] * 1e6 , n_wet_cld,  with_="steps lw 4 lc rgb 'magenta'")
    plot_rw_rin  = gp.PlotItems.Data(wet_edges[:-1] * 1e6 , n_wet_rin,  with_="steps lw 4 lc rgb 'blue'")
    
    g.plot(plot_rd_ini, plot_rd_cld, plot_rd_rin, plot_rw_cld, plot_rw_rin)
