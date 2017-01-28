"""
Prints how effective (in terms of the depleted moles / ug of dry air) was:
    oxidation in total (looking at the depletes S_IV)
    oxidation by H2O2  (looking at the depleted H2O2)
    oxidation by O3    (looking at the depleted O3)
"""
import numpy   as np
import h5py    as h5
import sys

# libcloudph++ bindings to python (to have access to library constants)
sys.path.insert(0, "../../../../../build/bindings/python/")
from libcloudphxx import common as cm

for case in ('case_base', 'case3', 'case4', 'case5', 'case6'):

    # open hdf5 files with data
    h5f_ini = h5.File('data/' + case + '/out_hall_pinsky_stratocumulus/timestep0000000000.h5', 'r')
    h5f_end = h5.File('data/' + case + '/out_hall_pinsky_stratocumulus/timestep0000011800.h5', 'r')
    
    # helper dict for chem names and molar mass
                 #name   gas molar mass   aqueous molar mass    label in hdf5     ini    spn   end
    help_dict = {
                'H2O2' : [cm.M_H2O2,      cm.M_H2O2,            'H2O2',            0,     0,    0],
                'O3'   : [cm.M_O3,        cm.M_O3,              'O3',              0,     0,    0],
                'SO2'  : [cm.M_SO2,       cm.M_SO2_H2O,         'S_IV',            0,     0,    0]
                }
    
    # calulate the initial and final number of moles per micro gram of dry air
    # and store them in dict
    # (additionally one could also check state after spinup)
    for key, val in help_dict.iteritems():
     
        name1 = key + "g"                   # gas phase chem species
        name2 = "chem_" + val[2] + "_aq"    # aq  phase chem species

        # choose indexes of grid-cell with r_c > 0.01
        cld_mixr_fin = h5f_end['rw_rng000_mom3'][:] * 4./3 * 3.14 * 1e3 * 1e3 #g/kg
        idx_cld      = np.where(cld_mixr_fin >  0)
        size_cld     = idx_cld[0].shape[0]
        size_all     = 76. * 76

        # total concentration [moles/ug of dry air]
        total_ini_conc = h5f_ini[name1][:] / val[0] + h5f_ini[name2][:] / val[1]
        total_end_conc = h5f_end[name1][:] / val[0] + h5f_end[name2][:] / val[1]
    
        # zoom in on the cloudy parts and calculate average over them
        ini = (total_ini_conc[:]).sum() / size_all * 1e9 
        end = (total_end_conc[:]).sum() / size_all * 1e9 
    
        val[3] = ini
        val[5] = end
 
    print " "
    print case
    print "-------------- % difference % --------------------------------"
    
    print "delta O3 / init O3     ", (help_dict['O3'][3] - help_dict['O3'][5]) / help_dict['O3'][3] * 100
    print "delta H2O2 / init H2O2 ", (help_dict['H2O2'][3] - help_dict['H2O2'][5]) / help_dict['H2O2'][3] * 100
    print "delta SO2 / init SO2   ", (help_dict['SO2'][3] - help_dict['SO2'][5]) / help_dict['SO2'][3] * 100

