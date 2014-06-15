import sys
sys.path.append(".")
import pytest
import inspect
import pdb

from numpy import array as arr_t
import analytic_blk_1m_pytest as anpy
from libcloudphxx  import common

"""it checks consistency between thermodynamical equations used in the libcloudph and my python tests""" 

@pytest.mark.parametrize("arg", [
    ({"T" : 273.,   "press" : 1.e5,  "rv" : 0.}), ({"T" : 273.,   "press" : 1.e5,  "rv" : 1.e-2}),
    ({"T" : 283.15, "press" : 9.e4,  "rv" : 0.}), ({"T" : 283.15, "press" : 9.e4,  "rv" : 1.e-2})
    ])
def test_thetadry(arg, epsilon = 0.01):
    rho_d = anpy.density_dry(**arg)
    th_d = anpy.pottemp_dry(**arg)
    temp_md = common.T(th_d, rho_d)
    assert abs(arg["T"] - temp_md) <= epsilon * arg["T"]


@pytest.mark.parametrize("arg", [
    ({"T" : 273.,   "press" : 1.e5,  "rv" : 0.}), ({"T" : 273.,   "press" : 1.e5,  "rv" : 1.e-2}),
    ({"T" : 283.15, "press" : 9.e4,  "rv" : 0.}), ({"T" : 283.15, "press" : 9.e4,  "rv" : 1.e-2})
        ])
def test_pressure(arg, epsilon = 0.01):
    rho_d = anpy.density_dry(**arg)
    press_md = common.p(rho_d, arg["rv"], arg["T"])
    assert abs(arg["press"] - press_md) <= epsilon * arg["press"]
                    

@pytest.mark.parametrize("arg", [{"T" : 273.}, {"T" : 283.}, {"T" : 293.}])
def test_pvs(arg, epsilon = 0.03):
    pvs_md = common.p_vs(arg["T"])
    pvs_anpy = anpy.press_sat(arg["T"])
    assert abs(pvs_anpy - pvs_md) <= epsilon * pvs_anpy
                    
