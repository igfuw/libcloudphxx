"""
Prints how many moles were depleted vs gained during the oxidation reaction
"""
import numpy   as np
import h5py    as h5
import sys
import math

# libcloudph++ bindings to python (to have access to library constants)
sys.path.insert(0, "../../../../../build/bindings/python/")
from libcloudphxx import common as cm

# path to directory with data
dir_path = ''

for case in ['case_base', 'case3', 'case4', 'case5', 'case6']:

    # open hdf5 files with data
    h5f_ini = h5.File(dir_path + 'out_' + case + '/timestep0000010000.h5', 'r') # model after spinup
    h5f_end = h5.File(dir_path + 'out_' + case + '/timestep0000011800.h5', 'r') # model at the end of simulation
    h5f_cst = h5.File(dir_path + 'out_' + case + '/const.h5', 'r')              # constant rhod 

    # dry air density                ... see the comment in "how_effective.py" test
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

    #helper dict for chem names and molar mass
                 #name   gas molar mass   aqueous molar mass    label in hdf5     ini    spn   end
    help_dict = {
                'H2O2' : [cm.M_H2O2,      cm.M_H2O2,            'H2O2',            0,     0,    0],
                'O3'   : [cm.M_O3,        cm.M_O3,              'O3',              0,     0,    0],
                'SO2'  : [cm.M_SO2,       cm.M_SO2_H2O,         'S_IV',            0,     0,    0],
                'CO2'  : [cm.M_CO2,       cm.M_CO2_H2O,         'C_IV',            0,     0,    0],
                'NH3'  : [cm.M_NH3,       cm.M_NH3_H2O,         'N_III',           0,     0,    0],
                'HNO3' : [cm.M_HNO3,      cm.M_HNO3,            'N_V',             0,     0,    0],
                'H2SO4': [0,              cm.M_H2SO4,           'S_VI',            0,     0,    0]
                }
    
    # calulate the initial (after spin-up) and final number of moles and store them in dict
    for key, val in help_dict.iteritems():
        name1 = key + "g"
        if key == 'H2SO4':
            name2 = 'chem_S_VI_aq'
            ini_tab = h5f_ini[name2][:]
            end_tab = h5f_end[name2][:]

            # total moles
            ini = (ini_tab * rho_d * dv).sum() / val[1]
            end = (end_tab * rho_d * dv).sum() / val[1]
        else:
            name2 = "chem_" + val[2] + "_aq"
            # total moles
            ini = ((h5f_ini[name1][:] * rho_d * dv) / val[0] + (h5f_ini[name2][:] * rho_d * dv) / val[1]).sum()
            end = ((h5f_end[name1][:] * rho_d * dv) / val[0] + (h5f_end[name2][:] * rho_d * dv) / val[1]).sum()
 
        val[3] = ini
        val[5] = end

    # how many moles of H2O2, O3 and S_IV were depleted?
    depleted_O3   = help_dict['O3'][3]    - help_dict['O3'][5] 
    depleted_H2O2 = help_dict['H2O2'][3]  - help_dict['H2O2'][5] 
    depleted_SO2  = help_dict['SO2'][3]   - help_dict['SO2'][5] 
    # how many moles of S_VI were created?
    gained_S6     = help_dict['H2SO4'][5] - help_dict['H2SO4'][3]
    
    relative_error_1 = (depleted_O3 + depleted_H2O2 - gained_S6) / gained_S6 
    relative_error_2 = (depleted_SO2 - gained_S6) / gained_S6 

    print " "
    print "----------------- " + case +  " --------------------------"
    print "  "
    print "depleted O3 and H2O2          = ", depleted_O3 + depleted_H2O2
    print "depleted S4                   = ", depleted_SO2
    print "gained S6                     = ", gained_S6
    print "  "
    print "error in O3 + H2O2 in realtion to delta S6 = ", relative_error_1 * 100, " %"
    print "error in S_IV      in relation to delta S6 = ", relative_error_2 * 100, " %"
