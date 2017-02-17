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
dir_path = '../../../build/tests/chem_sandbox/'

#for case in ['case_base', 'case3', 'case4', 'case5', 'case6']:
for case in [#'case6',\
             #'case6_no_coll_no_sedi',\
             #'case6_no_sedi',\
             #'case6_no_coll',\
             #'case6_yes_coll_yes_sedi_no_react',\
             #'case6_yes_coll_no_sedi_no_react',\
             #'case6_yes_coll_yes_sedi_yes_react_no_cloud',\
             #'case6_puddle'#,\
             #'case6_puddle_no_react'
             #'test6_no_coll_no_sedi_no_react',
             #'test6_no_coll_yes_sedi_no_react',
             "test6_no_coll_no_sedi_no_chem",
             "test6_no_coll_yes_sedi_no_chem",
             "test6_no_coll_yes_sedi_no_chem_pred_corr"
             #'test6_yes_coll_yes_sedi_no_react',
            ]:

    # open hdf5 files with data
    h5f_ini = h5.File(dir_path + 'out_' + case + '/timestep0000000000.h5', 'r')
    #h5f_end = h5.File(dir_path + 'out_' + case + '/timestep0000011800.h5', 'r')
    h5f_end = h5.File(dir_path + 'out_' + case + '/timestep0000000600.h5', 'r')

    # helper dict for chem names and molar mass
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
    
    # calulate the initial and final number of moles per micro gram of dry air
    # and store them in dict
    # (additionally one could also check state after spinup)
    #for key, val in help_dict.iteritems():
    #    name1 = key + "g"
    #    if key == 'H2SO4':
    #        name2 = 'chem_S_VI_aq'
    #        # moles/ug of dry air
    #        ini = (h5f_ini[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
    #        end = (h5f_end[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
    #    else:
    #        name2 = "chem_" + val[2] + "_aq"
    #        # moles/ug of dry air
    #        ini = (h5f_ini[name1][:] / val[0] + h5f_ini[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
    #        end = (h5f_end[name1][:] / val[0] + h5f_end[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
    #    val[3] = ini
    #    val[5] = end
 
    print " "
    print "----------------- " + case +  " --------------------------"
 
    ini_dry_vol = h5f_ini["rd_rng000_mom3"][:].sum() / 76. / 76. * 1e9 * 4./3 * math.pi * 1e3
    end_dry_vol = h5f_end["rd_rng000_mom3"][:].sum() / 76. / 76. * 1e9 * 4./3 * math.pi * 1e3

    out_H2SO4 = 44
    out_dry_vol = 44

    if (case == "test6_no_coll_yes_sedi_no_chem"):
        out_dry_vol = 2.4276329210759506e-11                         / 20 / 20 / 76 / 76 * 1e9 * 1e3
    elif (case == "test6_no_coll_no_sedi_no_chem"):
        out_dry_vol = 2.8202570487301237e-12                         / 20 / 20 / 76 / 76 * 1e9 * 1e3
    elif (case == "test6_no_coll_yes_sedi_no_chem_pred_corr"):
        out_dry_vol = 2.4276329210759506e-11                         / 20 / 20 / 76 / 76 * 1e9 * 1e3
 
    print " " 
    #print "ini  S6   = ", help_dict['H2SO4'][3], "mol/ug"
    #print "end  S6   = ", help_dict['H2SO4'][5], "mol/ug"
    #print "lost S6   = ", out_H2SO4,             "mol/ug"
    #print " " 
    print "ini  dry  = ", ini_dry_vol,           "dm3/ug"
    print "end  dry  = ", end_dry_vol,           "dm3/ug"
    print "lost dry  = ", out_dry_vol,           "dm3/ug"

    #if (case == "case6_puddle"): 
    #    out_H2SO4 = 5.1103535270158318e-07  / cm.M_H2SO4   / 20 / 20 / 76 / 76 * 1e9
    #elif (case == "case6_puddle_no_react"):
    #    out_H2SO4 = 5.1103535270158318e-07  / cm.M_H2SO4   / 20 / 20 / 76 / 76 * 1e9
    #if (case == "out_test6"):
    #    out_H2SO4 = 3.7237722239069626e-08 / cm.M_H2SO4   / 20 / 20 / 76 / 76 * 1e9
    #if (case == "test6_no_coll_no_sedi_no_react"):
    #    out_H2SO4 = 4.3260288684621173e-09            / cm.M_H2SO4   / 20 / 20 / 76 / 76 * 1e9
    #    out_dry_vol = 2.8202570487301237e-12                         / 20 / 20 / 76 / 76 * 1e9 * 1e3
    #elif (case == "test6_yes_coll_yes_sedi_no_react"):
    #    out_H2SO4 = 3.7237722239069626e-08            / cm.M_H2SO4   / 20 / 20 / 76 / 76 * 1e9
    #    out_dry_vol = 2.4276329210759506e-11                         / 20 / 20 / 76 / 76 * 1e9 * 1e3
    #elif (case == "test6_no_coll_yes_sedi_no_react"):
    #    out_H2SO4 = 3.7237722239069626e-08            / cm.M_H2SO4   / 20 / 20 / 76 / 76 * 1e9
    #    out_dry_vol =  2.4276329210759506e-11                        / 20 / 20 / 76 / 76 * 1e9 * 1e3

    #depleted_O3   = help_dict['O3'][3]    - help_dict['O3'][5] 
    #depleted_H2O2 = help_dict['H2O2'][3]  - help_dict['H2O2'][5] 
    #depleted_SO2  = help_dict['SO2'][3]   - help_dict['SO2'][5] 
    #gained_S6     = help_dict['H2SO4'][5] - help_dict['H2SO4'][3] + out_H2SO4 
    #
    #relative_error_1 = (depleted_O3 + depleted_H2O2 - gained_S6) / gained_S6 
    #relative_error_2 = (depleted_SO2 - gained_S6) / gained_S6 
 
    #print "  "
    #print "depleted O3 and H2O2          = ", depleted_O3 + depleted_H2O2
    #print "depleted S4                   = ", depleted_SO2
    #print "gained S6                     = ", gained_S6
    #print "lost S6 ", out_H2SO4
    #print "error in realtion to delta S6 = ", relative_error_1 * 100, " %"
    #print "error in relation to delta S6 = ", relative_error_2 * 100, " %"
