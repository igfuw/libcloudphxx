import sys
sys.path.append(".")
#sys.path.append("/Users/dorota/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages")
#sys.path.append("/Users/dorota/libcloudphxx/build/tests/python")
import pytest

from numpy import array as arr_t

from libcloudphxx import blk_1m

#TODO: mysle, ze moznaby tu dodac params?
# czy te cond, cevp to jakies atrybuty? nie mozna tego krocej zapisac?
# zrezygnowalam chwilowo z fixture

def opts_cr(cond = True, cevp = True, revp = True, conv = True,
            accr = True, sedi = False):
    opts = blk_1m.opts_t()
    opts.cond = cond
    opts.cevp = cevp
    opts.revp = revp
    opts.conv = conv
    opts.accr = accr
    opts.sedi = sedi
    return opts

# typical values a example
rhod_0 = arr_t([1.1  ])
th_0   = arr_t([291.8])
rv_0   = arr_t([8.e-3])
rc_0   = arr_t([5.e-4])
rr_0   = arr_t([0.  ])
dt_0   = 1

#ta f-cja jest tylko po to, aby byly keword arg., konieczna?
def adj_cellwise(opts, rhod = rhod_0, th = th_0,
                 rv = rv_0, rc = rc_0, rr = rr_0, dt = dt_0):
    blk_1m.adj_cellwise(opts, rhod, th, rv, rc, rr, dt)
    return rv

@pytest.mark.parametrize("arg", [
    {'th':arr_t([200])},    pytest.mark.xfail({'th':arr_t([500 ])}),
    pytest.mark.xfail({'rv':arr_t([-1.e-5])}), pytest.mark.xfail({'rv':arr_t([0.1 ])}),
    #pytest.mark.xfail({'rc':arr_t([-1.e-5])}), pytest.mark.xfail({'rc':arr_t([0.01])}),
    #pytest.mark.xfail({'rr':arr_t([-1.e-5])}), pytest.mark.xfail({'rr':arr_t([0.01])})
    ])
def test_exeptions_wrongvalue(arg):
    opts = opts_cr()
    with pytest.raises(Exception):
        adj_cellwise(opts, **arg) 

#TODO: wypisac znane outputy, moze tez theta?
@pytest.mark.skipif
@pytest.mark.parametrize("arg, expected", [
    ({"rv" :  arr_t([0.]),     "rc" :  arr_t([0.])},
     {"rv" :  arr_t([0.]),     "rc" :  arr_t([0.])}), # no water
    ({"rv" :  arr_t([7.e-3]),  "rc" :  arr_t([0.])},
     {"rv" :  arr_t([7.e-3]),  "rc" :  arr_t([0.])}), # no cl water and subsat.
    ({"rv" :  arr_t([10.e-3]), "rc" :  arr_t([0.])},
     {"rv" :  arr_t([8.6e-3]), "rc" :  arr_t([1.4 e-3])}), # no cl water and supersat.
    ({"rv" :  arr_t([5.e-3]),  "rc" :  arr_t([1.e-3])},
     {"rv" :  arr_t([6.e-3]),  "rc" :  arr_t([0.])}), # subsat. leads to coplete evap.
    ({"rv" :  arr_t([8.e-3]),  "rc" :  arr_t([1.e-3])},
     {"rv" :  arr_t([8.6e-3]), "rc" :  arr_t([0.4e-3])}), # subsat. leads to some evap.
    ({"rv" :  arr_t([9.e-3]),  "rc" :  arr_t([1.e-3])},
     {"rv" :  arr_t([8.6e-3]), "rc" :  arr_t([1.4e-3])}), # supersat. leads to cond.
    ])
#TODO zastanowic sie nad epsilonem
def test_expected_output_evapcond(arg, expected, epsilon = 0.1):
    opts = opts_cr(conv = False, accr = False )
    rv = adj_cellwise(opts, **arg)
    for key, value in expected.items():
        print key, value, eval(key)
        assert abs(eval(key) - value) <= epsilon * abs(value)

