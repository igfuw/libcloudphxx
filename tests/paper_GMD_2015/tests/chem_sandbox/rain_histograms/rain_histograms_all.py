import numpy   as np
import h5py    as h5
import Gnuplot as gp
import math    as mt

# all the collision kernels
kernels = {"out_hall", "out_hall_davis_no_waals", 
           "out_hall_pinsky_stratocumulus", 
            "out_vohl_davis_no_waals", 
            "out_onishi_hall", "out_onishi_hall_davis_no_waals"}

rain_cover = {"out_hall": 0, "out_hall_davis_no_waals": 0,  
              "out_hall_pinsky_stratocumulus" : 0, 
              "out_vohl_davis_no_waals" : 0, 
              "out_onishi_hall" : 0, "out_onishi_hall_davis_no_waals" : 0}

# only grid-cells with rain water mixing ratio greater than cutoff will be shown
cutoff = 0.01 #g/kg

# left and right edges of bins for dry and wet radius
dry_edges = np.fromiter( (1e-6 * 10**(-3 + i * .1) for i in xrange(41)), dtype="float")
wet_edges = np.fromiter( (1e-6 * 10**(-3 + i * .1) for i in xrange(56)), dtype="float")
num_dry   = dry_edges.size - 1
num_wet   = wet_edges.size - 1

for kernel in kernels:

    # helper arrays for storing initial and cutoff dry radii as well as final and cutoff wet radii
    n_dry_ini  = np.zeros(num_dry)
    n_dry_fin  = np.zeros(num_dry)
    n_wet_all  = np.zeros(num_wet)
    n_wet_rain = np.zeros(num_wet)

    n_run = 0.

    for run in {"seed_9", "seed_13", "seed_30", "seed_42", "seed_44"}:
        # open hdf5 files with data
        h5f_ini       = h5.File('data_for_rain_histograms/' + run + '/' + kernel + '/timestep0000000000.h5', 'r')
        h5f_fin       = h5.File('data_for_rain_histograms/' + run + '/' + kernel + '/timestep0000011800.h5', 'r')
        rain_mixr_fin = h5f_fin['rw_rng001_mom3'][:] * 4./3 * 3.14 * 1e3 * 1e3 #g/kg
    
        idx  = np.where(rain_mixr_fin >= cutoff)
        idx2 = np.where(rain_mixr_fin < cutoff)
        rain_cover[kernel] += idx[0].size
        n_run += 1.
        #print "for " + kernel + " the % of grid-cells with rain water mixing ratio > " + str(cutoff) + " is " +  str(idx[0].size/76./76.*100)
   
        for i in range(num_dry-1): # first bin id for total conc
            name = "rd_rng" + str(i+1).zfill(3) + "_mom0"
            tmp_ini  = 1e-6 * h5f_ini[name][:]
            tmp_fin  = 1e-6 * h5f_fin[name][:]
    
            n_dry_ini[i] += tmp_ini.sum() / 76. / 76. / (mt.log(dry_edges[i+1], 10) - mt.log(dry_edges[i], 10))
            n_dry_fin[i] += tmp_fin[idx].sum() / idx[0].size / (mt.log(dry_edges[i+1], 10) - mt.log(dry_edges[i], 10))
    
        for i in range(num_wet-2): # first two bins are for total cloud and rain water conc
            name = "rw_rng" + str(i+2).zfill(3) + "_mom0"
            tmp_fin  = 1e-6 * h5f_fin[name][:]
    
            n_wet_all[i]  += tmp_fin.sum() / 76. / 76. / (mt.log(wet_edges[i+1], 10) - mt.log(wet_edges[i], 10))
            n_wet_rain[i] += tmp_fin[idx].sum() / idx[0].size / (mt.log(wet_edges[i+1], 10) - mt.log(wet_edges[i], 10))
    
        # close hdf5 files
        h5f_ini.close()
        h5f_fin.close()

    # average over all runs
    n_dry_ini  /= n_run
    n_dry_fin  /= n_run
    n_wet_all  /= n_run
    n_wet_rain /= n_run
    rain_cover[kernel] = rain_cover[kernel] / n_run / 76 / 76 * 100 
    print kernel, "  ", rain_cover[kernel]

    # plotting
    ymin = 0.001
    ymax = 1000
    xmin = 0.001
    xmax = 500

    g = gp.Gnuplot()
    g('reset')
    g('set term svg dynamic enhanced font "Verdana, 14"')
    g('set output "' + kernel + '_rain_distr.svg" ')
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
 
    plot_rd_ini  = gp.PlotItems.Data(dry_edges[:-1] * 1e6 , n_dry_ini,  with_="steps lw 4 lc rgb 'coral'")
    plot_rd_fin  = gp.PlotItems.Data(dry_edges[:-1] * 1e6 , n_dry_fin,  with_="steps lw 4 lc rgb 'brown'")

    plot_rw_all  = gp.PlotItems.Data(wet_edges[:-1] * 1e6 , n_wet_all,  with_="steps lw 4 lc rgb 'purple'")
    plot_rw_rain = gp.PlotItems.Data(wet_edges[:-1] * 1e6 , n_wet_rain, with_="steps lw 4 lc rgb 'blue'")
    
    g.plot(plot_rd_ini, plot_rd_fin, plot_rw_rain, plot_rw_all)
