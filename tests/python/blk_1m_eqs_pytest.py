import sys
sys.path.append(".")
import pytest
import inspect
import pdb

from numpy import array as arr_t
import analytic_blk_1m_pytest as anpy

#TODO from libcloudphxx_blk_1m_pytest import adj_cellwise

def theta_model_test(temp,pres):
    return 273.

@pytest.mark.parametrize("arg, expected", [
        ({"T" : arr_t([273.]),   "p" : arr_t([1.e5])}, {"theta" : arr_t([0.])}),
        ({"T" : arr_t([283.15]), "p" : arr_t([9.e4])}, {"theta" : arr_t([0.])})
    ])
def test_theta(arg, expected, epsilon = 0.05):
    #TODO
    theta_md = theta_model_test(arg["T"], arg["p"])
    theta_anpy = anpy.pot_temp(**arg)
    assert abs(theta_anpy - theta_md) <= epsilon * theta_anpy


