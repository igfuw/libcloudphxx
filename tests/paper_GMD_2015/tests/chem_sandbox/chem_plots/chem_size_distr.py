import numpy   as np
import h5py    as h5
import Gnuplot as gp
import math    as mt

# all the collision kernels
kernels = {"out_hall_pinsky_stratocumulus"}

ub_cover = {"out_hall_pinsky_stratocumulus" : 0}
ua_cover = {"out_hall_pinsky_stratocumulus" : 0}
db_cover = {"out_hall_pinsky_stratocumulus" : 0}
da_cover = {"out_hall_pinsky_stratocumulus" : 0}

# left and right edges of bins for dry and wet radius
dry_edges = np.fromiter( (1e-6 * 10**(-3 + i * .1) for i in xrange(41)), dtype="float")
wet_edges = np.fromiter( (1e-6 * 10**(-3 + i * .1) for i in xrange(56)), dtype="float")
num_dry   = dry_edges.size - 1
num_wet   = wet_edges.size - 1

for kernel in kernels:

    # helper arrays for storing initial and cutoff dry radii as well as final and cutoff wet radii
    n_dry_ub = np.zeros(num_dry)
    n_dry_ua = np.zeros(num_dry)
    n_dry_db = np.zeros(num_dry)
    n_dry_da = np.zeros(num_dry)

    n_wet_ub = np.zeros(num_wet)
    n_wet_ua = np.zeros(num_wet)
    n_wet_db = np.zeros(num_wet)
    n_wet_da = np.zeros(num_wet)

    for run in {"seed_44"}:
        # open hdf5 files with data
        h5f_cst   = h5.File('case_base/const.h5', 'r')
        h5f_ini   = h5.File('case_base/timestep0000000000.h5', 'r')
        h5f_fin   = h5.File('case_base/timestep0000011800.h5', 'r')

        x_grid = h5f_cst['X'][:-1][:-1]
        y_grid = h5f_cst['Y'][:-1][:-1]
        
        idx_ub = np.where((x_grid > 1 ) & (x_grid < 37) & (y_grid > 15) & (y_grid < 45))
        idx_ua = np.where((x_grid > 1 ) & (x_grid < 37) & (y_grid > 45) & (y_grid < 70))
        idx_db = np.where((x_grid > 39) & (x_grid < 76) & (y_grid > 15) & (y_grid < 45))
        idx_da = np.where((x_grid > 39) & (x_grid < 76) & (y_grid > 45) & (y_grid < 70))

        print idx_ub
        #print idx_ub.shape

        # count dry cloudy and rainy grid cells
        ub_cover[kernel] += idx_ub[0].size
        ua_cover[kernel] += idx_ua[0].size
        db_cover[kernel] += idx_db[0].size
        da_cover[kernel] += idx_da[0].size

        #print "for " + kernel + " the % of grid-cells with rain water mixing ratio > " + str(cutoff) + " is " +  str(idx[0].size/76./76.*100)
   
        for i in range(num_dry-1): # first bin id for total conc
            name = "rd_rng" + str(i+1).zfill(3) + "_mom0"
            tmp_ini  = 1e-6 * h5f_ini[name][:]
            tmp_fin  = 1e-6 * h5f_fin[name][:]

            n_dry_ub[i] += tmp_fin[idx_ub].sum() / idx_ub[0].size / (mt.log(dry_edges[i+1], 10) - mt.log(dry_edges[i], 10))
            n_dry_ua[i] += tmp_fin[idx_ua].sum() / idx_ua[0].size / (mt.log(dry_edges[i+1], 10) - mt.log(dry_edges[i], 10))
            n_dry_db[i] += tmp_fin[idx_db].sum() / idx_db[0].size / (mt.log(dry_edges[i+1], 10) - mt.log(dry_edges[i], 10))
            n_dry_da[i] += tmp_fin[idx_da].sum() / idx_da[0].size / (mt.log(dry_edges[i+1], 10) - mt.log(dry_edges[i], 10))
   
        for i in range(num_wet-2): # first two bins are for total cloud and rain water conc
            name = "rw_rng" + str(i+2).zfill(3) + "_mom0"
            tmp_fin  = 1e-6 * h5f_fin[name][:]

            n_wet_ub[i] += tmp_fin[idx_ub].sum() / idx_ub[0].size / (mt.log(wet_edges[i+1], 10) - mt.log(wet_edges[i], 10))
            n_wet_ua[i] += tmp_fin[idx_ua].sum() / idx_ua[0].size / (mt.log(wet_edges[i+1], 10) - mt.log(wet_edges[i], 10))
            n_wet_db[i] += tmp_fin[idx_db].sum() / idx_db[0].size / (mt.log(wet_edges[i+1], 10) - mt.log(wet_edges[i], 10))
            n_wet_da[i] += tmp_fin[idx_da].sum() / idx_da[0].size / (mt.log(wet_edges[i+1], 10) - mt.log(wet_edges[i], 10))
    
        # close hdf5 files
        h5f_ini.close()
        h5f_fin.close()

    ub_cover[kernel] = ub_cover[kernel] / 76. / 76 * 100 
    ua_cover[kernel] = ua_cover[kernel] / 76. / 76 * 100 
    db_cover[kernel] = db_cover[kernel] / 76. / 76 * 100 
    da_cover[kernel] = da_cover[kernel] / 76. / 76 * 100 

    print kernel, "  ", ub_cover[kernel]
    print kernel, "  ", ua_cover[kernel]
    print kernel, "  ", db_cover[kernel]
    print kernel, "  ", da_cover[kernel]

    # plotting
    ymin = 0.001
    ymax = 1000
    xmin = 0.001
    xmax = 5

    g = gp.Gnuplot()
    g('reset')
    g('set term svg dynamic enhanced font "Verdana, 14"')
    g('set output "' + kernel + '_aerosol_ud.svg" ')
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
    
    plot_rd_ub  = gp.PlotItems.Data(dry_edges[:-1] * 1e6 , n_dry_ub,  with_="steps lw 4 lc rgb 'black'")
    plot_rd_ua  = gp.PlotItems.Data(dry_edges[:-1] * 1e6 , n_dry_ua,  with_="steps lw 4 lc rgb 'green'")
    plot_rd_db  = gp.PlotItems.Data(dry_edges[:-1] * 1e6 , n_dry_db,  with_="steps lw 4 lc rgb 'red'")
    plot_rd_da  = gp.PlotItems.Data(dry_edges[:-1] * 1e6 , n_dry_da,  with_="steps lw 4 lc rgb 'blue'")

    g.plot(plot_rd_ub, plot_rd_ua, plot_rd_db, plot_rd_da)

    # plotting
    ymin = 0.001
    ymax = 1000
    xmin = 0.001
    xmax = 500

    g = gp.Gnuplot()
    g('reset')
    g('set term svg dynamic enhanced font "Verdana, 14"')
    g('set output "' + kernel + '_rain_ud.svg" ')
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

    plot_rw_ub  = gp.PlotItems.Data(wet_edges[:-1] * 1e6 , n_wet_ub,  with_="steps lw 4 lc rgb 'black'")
    plot_rw_ua  = gp.PlotItems.Data(wet_edges[:-1] * 1e6 , n_wet_ua,  with_="steps lw 4 lc rgb 'green'")
    plot_rw_db  = gp.PlotItems.Data(wet_edges[:-1] * 1e6 , n_wet_db,  with_="steps lw 4 lc rgb 'red'")
    plot_rw_da  = gp.PlotItems.Data(wet_edges[:-1] * 1e6 , n_wet_da,  with_="steps lw 4 lc rgb 'blue'")
    
    g.plot(plot_rw_ub, plot_rw_ua, plot_rw_db, plot_rw_da) 
