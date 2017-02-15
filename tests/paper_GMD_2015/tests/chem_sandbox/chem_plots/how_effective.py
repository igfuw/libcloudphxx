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

# path to build directory with data
dir_path = '../../../build/tests/chem_sandbox/'

for case in ['case_base', 'case3', 'case4', 'case5', 'case6', 'case1a']:

    # open hdf5 files with data
    h5f_ini = h5.File(dir_path + 'out_' + case + '/timestep0000000000.h5', 'r')
    h5f_end = h5.File(dir_path + 'out_' + case + '/timestep0000011800.h5', 'r')
    
    # helper dict for chem names and molar mass
                 #name   gas molar mass   aqueous molar mass    label in hdf5     ini    spn   end
    help_dict = {
                'H2O2' : [cm.M_H2O2,      cm.M_H2O2,            'H2O2',            0,     0,    0],
                'O3'   : [cm.M_O3,        cm.M_O3,              'O3',              0,     0,    0],
                'SO2'  : [cm.M_SO2,       cm.M_SO2_H2O,         'S_IV',            0,     0,    0],
                'H2SO4': [0,              cm.M_H2SO4,           'S_VI',            0,     0,    0]
                }
    
    # calulate the initial and final number of moles per micro gram of dry air
    # and store them in dict
    # (additionally one could also check state after spinup)
    for key, val in help_dict.iteritems():
     
        name1 = key + "g"                   # gas phase chem species
        size_all     = 76. * 76

        if key == 'H2SO4':
            name2 = 'chem_S_VI_aq'

            # total concentration [moles/ug of dry air]
            ini = (h5f_ini[name2][:] / val[1]).sum() / size_all * 1e9
            end = (h5f_end[name2][:] / val[1]).sum() / size_all * 1e9

        else:
            name2 = "chem_" + val[2] + "_aq"    # aq  phase chem species

            # total concentration [moles/ug of dry air]
            ini = (h5f_ini[name1][:].sum() / val[0] + h5f_ini[name2][:].sum() / val[1]) / size_all * 1e9
            end = (h5f_end[name1][:].sum() / val[0] + h5f_end[name2][:].sum() / val[1]) / size_all * 1e9
    
        val[3] = ini
        val[5] = end

    print " "
    print "---------------------- " + case + " ------------------------------"
   
    dn_S6   = (help_dict['H2SO4'][5] - help_dict['H2SO4'][3]) 
    dn_S4   = (help_dict['SO2'][3]   - help_dict['SO2'][5]) 
    dn_O3   = (help_dict['O3'][3]    - help_dict['O3'][5]) 
    dn_H2O2 = (help_dict['H2O2'][3]  - help_dict['H2O2'][5]) 
 
    print "init. SO2 is depleted by ", (help_dict['SO2'][3] - help_dict['SO2'][5])   / help_dict['SO2'][3]  * 100, "%"
    print "oxidation by O3          ", dn_O3   / dn_S4 * 100, "%" , " vs ", dn_O3   / dn_S6 * 100, "%"
    print "oxidation by H2O2        ", dn_H2O2 / dn_S4 * 100, "%" , " vs ", dn_H2O2 / dn_S6 * 100, "%"
 
    print " "
    print "initial mass of S6       ", help_dict['H2SO4'][3] * val[1] , "ug/kg of dry air"
    #print "fin S6   ", help_dict['H2SO4'][5] * val[1] , "ug/kg of dry air"
    print "increase of mass of S6   ", help_dict['H2SO4'][5] * val[1] - help_dict['H2SO4'][3] * val[1], "ug/kg of dry air"
    print "final S6 is increased by ", (help_dict['H2SO4'][5] - help_dict['H2SO4'][3]) / help_dict['H2SO4'][5] * 100, "%"


