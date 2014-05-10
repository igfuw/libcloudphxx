import sys
sys.path.append(".")
#sys.path.append("/Users/dorota/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages")
#sys.path.append("/Users/dorota/libcloudphxx/build/tests/python")
import pytest
import inspect
import pdb

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
def adj_cellwise(opts, rhod = rhod_0.copy(), th = th_0.copy(),
                 rv = rv_0.copy(), rc = rc_0.copy(), rr = rr_0.copy(), dt = dt_0):
    print "\n In adj_cellwise. Who is calling..?", inspect.stack()[1][3]
    print "th, rv, rc initial ", th, rv, rc
    print "th_0, rv_0, rc_0 initial", th_0, rv_0, rc_0
    #pdb.set_trace()
    blk_1m.adj_cellwise(opts, rhod, th, rv, rc, rr, dt)
    print "th, rv, rc", th, rv, rc
    print "th_0, rv_0, rc_0", th_0, rv_0, rc_0
    #pdb.set_trace()
    return rv, rc

#@pytest.mark.skipif
@pytest.mark.parametrize("arg", [
    {'th':arr_t([255.])},    pytest.mark.xfail({'th':arr_t([500. ])}),
    {'rv':arr_t([-1.e-5])}, pytest.mark.xfail({'rv':arr_t([0.1 ])}),
    {'rc':arr_t([-1.e-5])}, pytest.mark.xfail({'rc':arr_t([0.01])}),
    {'rr':arr_t([-1.e-5])}, pytest.mark.xfail({'rr':arr_t([0.01])})
    ])
def test_exeptions_wrongvalue(arg):
    print "\n jestem w test_exeption", arg
    opts = opts_cr()
    with pytest.raises(Exception) as excinfo:
#        pdb.set_trace()
        adj_cellwise(opts, **arg)
    #the line below can give you information about the exception 
    #print "exception info:", excinfo.value
    

#TODO: wypisac znane outputy, moze tez theta?
@pytest.mark.skipif
@pytest.mark.parametrize("arg, expected", [
#ok        ({"rv" :  arr_t([0.]),     "rc" :  arr_t([0.])},
#ok         {"rv" :  arr_t([0.]),     "rc" :  arr_t([0.])}), # no water
#ok        ({"rv" :  arr_t([7.e-3]),  "rc" :  arr_t([0.])},
#ok         {"rv" :  arr_t([7.e-3]),  "rc" :  arr_t([0.])}), # no cl water and subsat.
#        ({"rv" :  arr_t([10.e-3]), "rc" :  arr_t([0.])},
#         {"rv" :  arr_t([8.6e-3]), "rc" :  arr_t([1.4e-3])}), # no cl water and supersat.
#ok        ({"rv" :  arr_t([5.e-3]),  "rc" :  arr_t([1.e-3])},
#ok         {"rv" :  arr_t([6.e-3]),  "rc" :  arr_t([0.])}), # subsat. leads to coplete evap.
        ({"rv" :  arr_t([8.e-3]),  "rc" :  arr_t([1.e-3])},
         {"rv" :  arr_t([8.6e-3]), "rc" :  arr_t([0.4e-3])}), # subsat. leads to some evap.
#        ({"rv" :  arr_t([9.e-3]),  "rc" :  arr_t([1.e-3])},
#         {"rv" :  arr_t([8.6e-3]), "rc" :  arr_t([1.4e-3])}), # supersat. leads to cond.
    ])
#TODO zastanowic sie nad epsilonem
def test_expected_output_evapcond(arg, expected, epsilon = 0.1):
    opts = opts_cr(conv = False, accr = False )
    rv, rc = adj_cellwise(opts, **arg)
    for key, value in expected.items():
        #pdb.set_trace()
        print "\n key, valuu, eval(key)", key, value, eval(key)
        assert abs(eval(key) - value) <= epsilon * abs(value)

