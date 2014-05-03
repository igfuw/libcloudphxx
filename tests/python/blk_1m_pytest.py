import sys
sys.path.append(".")
sys.path.append("/Users/dorota/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages")
sys.path.append("/Users/dorota/libcloudphxx/build/tests/python")
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
rhod_0 = arr_t([1.  ])
th_0   = arr_t([300.])
rv_0   = arr_t([8.e-3])
rc_0   = arr_t([5.e-4])
rr_0   = arr_t([0.  ])
dt_0   = 1

#ta f-cja jest tylko po to, aby byly keword arg., konieczna?
def adj_cellwise(opts, rhod = rhod_0, th = th0,
                 rv = rv_0, rc = rc_0, rr = rr_0, dt = dt_0):
    return  blk_1m.adj_cellwise(opts, rhod, th, rv, rc, rr, dt)

@pytest.mark.parametrize("arg", [
    {'th':arr_t([200])},    pytest.mark.xfail({'th':arr_t([500 ])}),
    {'rv':arr_t([-1.e-5])}, pytest.mark.xfail({'rv':arr_t([0.1 ])}),
    {'rc':arr_t([-1.e-5])}, pytest.mark.xfail({'rc':arr_t([0.01])}),
    {'rr':arr_t([-1.e-5])}, pytest.mark.xfail({'rr':arr_t([0.01])})
    ])
def test_exeptions_wrongvalue(arg):
    opts = opts_cr()
    with pytest.raises(Exception):
        adj_cellwise(opts, **arg) 

#TODO: wypisac znane outputy
@pytest.mark.parametrize("arg, expected", [
    ({}, ),
    ])
#TODO zastanowic sie nad epsilonem
def test_expected_output(arg, expected, epsilon = 0.05):
    opts = opts_cr()
    assert abs(adj_cellwise(opts, **arg) - expected) <= epsilon * abs(expected)



#dot_rc = arr_t([0.])
#dot_rr = arr_t([0.])
#blk_1m.rhs_cellwise(opts, dot_rc, dot_rr, rc, rr)
#print "Python dot_rc =", dot_rc
