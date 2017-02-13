"""
Prints how many moles were depleted vs gained during the oxidation reaction
"""
import numpy   as np
import h5py    as h5
import sys

# libcloudph++ bindings to python (to have access to library constants)
sys.path.insert(0, "../../../../../build/bindings/python/")
from libcloudphxx import common as cm

#for case in ['case_base', 'case_base_rk', 'case3', 'case4', 'case4_no_O3', 'case5']:
for case in ['moles']:

    # open hdf5 files with data
    h5f_ini = h5.File('data/' + case + '/out_hall_pinsky_stratocumulus/timestep0000000000.h5', 'r')
    #h5f_spn = h5.File(case'/out_hall_pinsky_stratocumulus/timestep0000010000.h5', 'r')
    h5f_end = h5.File('data/' + case + '/out_hall_pinsky_stratocumulus/timestep0000011800.h5', 'r')
    
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
    for key, val in help_dict.iteritems():
     
        name1 = key + "g"
    
        if key == 'H2SO4':
            name2 = 'chem_S_VI_aq'
    
            # moles/ug of dry air
            ini = (h5f_ini[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
            #spn = (h5f_spn[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
            end = (h5f_end[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
    
        else:
            name2 = "chem_" + val[2] + "_aq"
    
            # moles/ug of dry air
            ini = (h5f_ini[name1][:] / val[0] + h5f_ini[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
            #spn = (h5f_spn[name1][:] / val[0] + h5f_spn[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
            end = (h5f_end[name1][:] / val[0] + h5f_end[name2][:] / val[1]).sum() / 76. / 76. * 1e9 
    
        val[3] = ini
        #val[4] = spn
        val[5] = end
 
#print " "
#print "-------------- init vs end of spin-up ----------------------"
#
#for key in ['CO2', 'NH3', 'HNO3', 'SO2', 'O3', 'H2O2', 'H2SO4']:
#    
#    relative_error = abs(help_dict[key][4] - help_dict[key][3]) / help_dict[key][3]
#    print key , " relative error ", relative_error * 100, " %"
#
    print " "
    print case
 
    print " "
    print "-------------- init vs end  --------------------------------"
    
    for key in ['CO2', 'NH3', 'HNO3']:
        
        relative_error = abs(help_dict[key][5] - help_dict[key][3]) / help_dict[key][3]
        print key , " relative error ", relative_error * 100, " %"
        
    print " "
    print "-------------- test react  --------------------------------"
        
    depleted_O3   = help_dict['O3'][3]    - help_dict['O3'][5] 
    depleted_H2O2 = help_dict['H2O2'][3]  - help_dict['H2O2'][5] 
    depleted_SO2  = help_dict['SO2'][3]   - help_dict['SO2'][5] 
    gained_S6     = help_dict['H2SO4'][5] - help_dict['H2SO4'][3] 
    
    relative_error_1 = (depleted_O3 + depleted_H2O2 - gained_S6) / gained_S6 
    relative_error_2 = (depleted_SO2 - gained_S6) / gained_S6 
    
    print " "
    print "error in realtion to delta S6 = ", relative_error_1 * 100, " %"
    print "depleted O3 and H2O2          = ", depleted_O3 + depleted_H2O2
    print "gained S6                     = ", gained_S6
    print " "
    print "error in relation to delta S6 =  ", relative_error_2 * 100, " %"
    print "depleted S4                   = ", depleted_SO2
    print "gained S6                     = ", gained_S6
    print " "
    print "ini S6   = ",  help_dict['H2SO4'][3]
    print "end S6   = ",   help_dict['H2SO4'][5]

