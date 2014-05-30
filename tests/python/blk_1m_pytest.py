import sys
sys.path.append(".")
sys.path.append("../../../tests/python/wrf_microphys/kessler/")
#sys.path.append("/Users/dorota/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages")
#sys.path.append("/Users/dorota/libcloudphxx/build/tests/python")
import pytest
import inspect
import pdb

from numpy import array as arr_t

from libcloudphxx_blk_1m_pytest import adj_cellwise
#uncomment if you want wrf kessler; should be run from different directory TODO!
#from wrf_blk_1m_pytest import adj_cellwise

# typical values a example
press_0 = arr_t([900.e2  ])
th_0   = arr_t([291.8])
T_0    = arr_t([283.15])
rv_0   = arr_t([8.e-3])
rc_0   = arr_t([5.e-4])
rr_0   = arr_t([0.  ])
dt_0   = 1

def condensation(press = None, T = None,
                 rv = None, rc = None, rr = None, dt=dt_0):

    print "\n na poczatku srawdzam, czy ma qv", rv
    #pdb.set_trace()
    T = th if T!=None else T_0.copy()
    rv = rv if rv!=None else rv_0.copy()
    rc = rc if rc!=None else rc_0.copy()
    rr = rr if rr!=None else rr_0.copy()
    press = press if press!=None else press_0.copy()
    print "\n In condensation. Who is calling..?", inspect.stack()[1][3]
    rv, rc = adj_cellwise(press, T, rv, rc, rr, dt)
    return rv, rc


#@pytest.mark.skipif
@pytest.mark.parametrize("arg", [
    {'T':arr_t([255.])},    pytest.mark.xfail({'T':arr_t([500. ])}),
    {'rv':arr_t([-1.e-5])}, pytest.mark.xfail({'rv':arr_t([0.1 ])}),
    {'rc':arr_t([-1.e-5])}, pytest.mark.xfail({'rc':arr_t([0.01])}),
    {'rr':arr_t([-1.e-5])}, pytest.mark.xfail({'rr':arr_t([0.01])})
    ])
def test_exeptions_wrongvalue(arg):
    print "\n jestem w test_exeption", arg
    with pytest.raises(Exception) as excinfo:
        condensation(**arg)
    #the line below can give you information about the exception 
    #print "exception info:", excinfo.value
    

#TODO: wypisac znane outputy, moze tez theta?
#TODO: polaczyc z plikiem analytic_blk_1m_pytest??
#@pytest.mark.skipif
@pytest.mark.parametrize("arg, expected", [
        ({"rv" :  arr_t([0.]),     "rc" :  arr_t([0.])},
         {"rv" :  arr_t([0.]),     "rc" :  arr_t([0.])}), # no water
        ({"rv" :  arr_t([7.e-3]),  "rc" :  arr_t([0.])},
         {"rv" :  arr_t([7.e-3]),  "rc" :  arr_t([0.])}), # no cl water and subsat.
        ({"rv" :  arr_t([10.e-3]), "rc" :  arr_t([0.])},
         {"rv" :  arr_t([9.44e-3]), "rc" :  arr_t([.56e-3])}), # no cl water and supersat.
        ({"rv" :  arr_t([5.e-3]),  "rc" :  arr_t([1.e-3])},
         {"rv" :  arr_t([6.e-3]),  "rc" :  arr_t([0.])}), # subsat. leads to coplete evap.
        ({"rv" :  arr_t([8.e-3]),  "rc" :  arr_t([1.e-3])},
         {"rv" :  arr_t([8.26e-3]), "rc" :  arr_t([0.74e-3])}), # subsat. leads to some evap.
        ({"rv" :  arr_t([9.e-3]),  "rc" :  arr_t([1.e-3])},
         {"rv" :  arr_t([8.85e-3]), "rc" :  arr_t([1.15e-3])}), # supersat. leads to cond.
    ])
#TODO zastanowic sie nad epsilonem
def test_expected_output_evapcond(arg, expected, epsilon = 0.13):
    print "\n w test_expected value przed", arg
    rv, rc = condensation(**arg)
    #print "rv, rc po", rv, rc
    for key, value in expected.items():
        print "\n key, valuu, eval(key)", key, value, eval(key)
        assert abs(eval(key) - value) <= epsilon * abs(value)

