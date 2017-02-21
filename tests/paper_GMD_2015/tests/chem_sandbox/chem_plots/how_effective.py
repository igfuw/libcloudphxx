"""
Prints how effective (in terms of the depleted moles / ug of dry air) was:
    oxidation in total (looking at the depletes S_IV)
    oxidation by H2O2  (looking at the depleted H2O2)
    oxidation by O3    (looking at the depleted O3)
"""
import numpy   as np
import h5py    as h5
import sys
import math

# libcloudph++ bindings to python (to have access to library constants)
sys.path.insert(0, "../../../../../build/bindings/python/")
from libcloudphxx import common as cm

# path to build directory with data
dir_path = '../../../build/tests/chem_sandbox/'

for case in ['case_base', 'case3', 'case4', 'case5', 'case6']:

    # open hdf5 files with data
    h5f_ini = h5.File(dir_path + 'out_' + case + '/timestep0000010000.h5', 'r')
    h5f_end = h5.File(dir_path + 'out_' + case + '/timestep0000011800.h5', 'r')
    h5f_cst = h5.File(dir_path + 'out_' + case + '/const.h5', 'r')
    
    # dry air density
    rho_d   = h5f_cst["G"][:]
    # volume of grid cells
    dv = np.ones(shape=(76,76))
    dv[0,:]   = 0.5
    dv[-1,:]  = 0.5
    dv[:,0]   = 0.5
    dv[:,-1]  = 0.5
    dv[0,0]   = 0.25
    dv[-1,-1] = 0.25
    dv[0,-1]  = 0.25
    dv[-1,0]  = 0.25
    dv *= 20.
   
    # total domain volume 
    total_vol = dv.sum()

    # helper dict for chem names and molar mass
                 #name   gas molar mass   aqueous molar mass    label in hdf5      spn   end
    help_dict = {
                'H2O2' : [cm.M_H2O2,      cm.M_H2O2,            'H2O2',            0,    0],
                'O3'   : [cm.M_O3,        cm.M_O3,              'O3',              0,    0],
                'SO2'  : [cm.M_SO2,       cm.M_SO2_H2O,         'S_IV',            0,    0],
                'H2SO4': [0,              cm.M_H2SO4,           'S_VI',            0,    0]
                }
    
    # calulate the initial and final number of moles and store them in dict
    # (additionally one could also check state after spinup)
    for key, val in help_dict.iteritems():
     
        name1 = key + "g"                   # gas phase chem species

        if key == 'H2SO4':
            name2 = 'chem_S_VI_aq'

            # total moles
            ini = (h5f_ini[name2][:] * rho_d * dv).sum() / val[1] 
            end = (h5f_end[name2][:] * rho_d * dv).sum() / val[1]

        else:
            name2 = "chem_" + val[2] + "_aq"    # aq  phase chem species

            # total moles 
            ini = (h5f_ini[name1][:] * dv * rho_d).sum() / val[0] + (h5f_ini[name2][:] * dv * rho_d).sum() / val[1]
            end = (h5f_end[name1][:] * dv * rho_d).sum() / val[0] + (h5f_end[name2][:] * dv * rho_d).sum() / val[1]
    
        val[3] = ini
        val[4] = end

    dn_S6   = (help_dict['H2SO4'][4] - help_dict['H2SO4'][3]) 
    dn_S4   = (help_dict['SO2'][3]   - help_dict['SO2'][4]) 
    dn_O3   = (help_dict['O3'][3]    - help_dict['O3'][4]) 
    dn_H2O2 = (help_dict['H2O2'][3]  - help_dict['H2O2'][4]) 
 
    ini_dry_mass = 4./3 * math.pi * (h5f_ini["rd_rng000_mom3"][:] * dv * rho_d).sum() * 1800 / total_vol * 1e9 # ug/m3 dry air
    fin_dry_mass = 4./3 * math.pi * (h5f_end["rd_rng000_mom3"][:] * dv * rho_d).sum() * 1800 / total_vol * 1e9 # ug/m3 dry air
    ini_dry_mass_chem = help_dict['H2SO4'][3] * (cm.M_NH3 + cm.M_H2SO4) / total_vol * 1e9                      # ug/m3 dry air


    print " "
    print "---------------------- " + case + " ------------------------------"
    print "init. SO2 is depleted by ", (help_dict['SO2'][3] - help_dict['SO2'][4])   / help_dict['SO2'][3]  * 100, "%"
    print "oxidation by O3          ", dn_O3   / dn_S6 * 100, "%"
    print "oxidation by H2O2        ", dn_H2O2 / dn_S6 * 100, "%"
    print " "
    #print "initial mass of S6       ", help_dict['H2SO4'][3] * val[1] / total_vol * 1e9, "ug/m3 of dry air"
    #print "final mass of S6         ", help_dict['H2SO4'][4] * val[1] / total_vol * 1e9, "ug/m3 of dry air"
    #print " "
    #print "increase of mass of S6   ", (help_dict['H2SO4'][4] - help_dict['H2SO4'][3]) * val[1] / total_vol * 1e9, "ug/m3 of dry air"
    #print "final S6 is increased by ", (help_dict['H2SO4'][4] - help_dict['H2SO4'][3]) / help_dict['H2SO4'][4] * 100, "%"
    print "ini particulate mass (r dry)  ", ini_dry_mass, "ug/m3 dry air"
    print "ini particulate mass (chem)   ", ini_dry_mass_chem, "ug/m3 dry air"
    print "fin particulate mass (r_dry)  ", fin_dry_mass, "ug/m3 dry air"
    print "delta particulate mass        ", fin_dry_mass - ini_dry_mass, "ug/m3 dry air"
    print " "
    print "init mass [% final mass]      ", ini_dry_mass / fin_dry_mass * 100, "%"
    print "created mass [% final mass]   ", (fin_dry_mass - ini_dry_mass) / fin_dry_mass * 100, "%"

