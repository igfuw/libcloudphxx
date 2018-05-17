import sys
sys.path.insert(0, "../../bindings/python/")

from numpy import array as arr_t # ndarray dtype default to float64, while array's is int64!

from libcloudphxx import blk_1m
from libcloudphxx import common

opts = blk_1m.opts_t()
r_eps_def = opts.r_eps
nwtrph_iters_def = opts.nwtrph_iters

def supersaturation(T, p, rv):
    r_vs = common.r_vs(T, p)
    ss = 100. * (rv / r_vs - 1)
    return ss

def initial_state(init_sup_sat):
    rhod = arr_t([1.  ])
    th   = arr_t([300.])
    if init_sup_sat:
        rv   = arr_t([0.02])
    else:
        rv   = arr_t([0.002])
    rc   = arr_t([0.015])
    rr   = arr_t([0.  ])
    dt   = 1
    
    T = common.T(th[0], rhod[0])
    p = common.p(rhod[0], rv[0], T)
    ss = supersaturation(T, p, rv[0])
    print "initial supersaturation", ss

    return rhod, th, rv, rc, rr, dt

def test_adj_cellwise(init_sup_sat, r_eps = r_eps_def):
    opts.r_eps = r_eps
    print "[standard adj_cellwise]"
    rhod, th, rv, rc, rr, dt = initial_state(init_sup_sat)
    blk_1m.adj_cellwise(opts, rhod, th, rv, rc, rr, dt)
    
    T = common.T(th[0], rhod[0])
    p = common.p(rhod[0], rv[0], T)
    ss = supersaturation(T, p, rv[0])
    print "final supersaturation", ss
    return ss

def test_adj_cellwise_constp(init_sup_sat, r_eps = r_eps_def):
    opts.r_eps = r_eps
    print "[constp adj_cellwise]"
    rhod, th, rv, rc, rr, dt = initial_state(init_sup_sat)

    # define pressure consistent with adj_cellwise to compare results
    p   = arr_t([common.p(rhod[0], rv[0], common.T(th[0], rhod[0]))])
    p_d = arr_t([p[0] - common.p_v(p[0], rv[0])])

    blk_1m.adj_cellwise_constp(opts, rhod, p, p_d, th, rv, rc, rr, dt)
    
    T = common.exner(p_d[0]) * th[0]
    ss = supersaturation(T, p[0], rv[0])
    print "final supersaturation", ss
    return ss

def test_adj_cellwise_nwtrph(init_sup_sat, nwtrph_iters = nwtrph_iters_def):
    opts.nwtrph_iters = nwtrph_iters
    print "[nwtrph adj_cellwise]"
    rhod, th, rv, rc, rr, dt = initial_state(init_sup_sat)

    # define pressure consistent with adj_cellwise to compare results
    p   = arr_t([common.p(rhod[0], rv[0], common.T(th[0], rhod[0]))])
    p_d = arr_t([p[0] - common.p_v(p[0], rv[0])])

    blk_1m.adj_cellwise_nwtrph(opts, p, p_d, th, rv, rc, dt)
    
    T = common.exner(p_d[0]) * th[0]
    ss = supersaturation(T, p[0], rv[0])
    print "final supersaturation", ss
    return ss

eps = {
        # supersaturation
        True  : { 'org'      : 3e-2, 'constp'      : 3e-2, 'nwtrph'      : 1e-3,
                  'org_prec' : 6e-4, 'constp_prec' : 6e-4, 'nwtrph_prec' : 1e-12 },
        # subsaturation
        False : { 'org'      : 0.5 , 'constp'      : 0.5 , 'nwtrph'      : 1e-4,
                  'org_prec' : 2e-3, 'constp_prec' : 2e-3, 'nwtrph_prec' : 1e-13 }
      }

for init_sup_sat in [True, False]:
    # default precision
    assert(abs(test_adj_cellwise(init_sup_sat)) < eps[init_sup_sat]['org'])
    assert(abs(test_adj_cellwise_constp(init_sup_sat)) < eps[init_sup_sat]['constp'])
    assert(abs(test_adj_cellwise_nwtrph(init_sup_sat)) < eps[init_sup_sat]['nwtrph'])
    
    # cranked up precision
    assert(abs(test_adj_cellwise(init_sup_sat, 1e-7)) < eps[init_sup_sat]['org_prec'])
    assert(abs(test_adj_cellwise_constp(init_sup_sat, 1e-7)) < eps[init_sup_sat]['constp_prec'])
    assert(abs(test_adj_cellwise_nwtrph(init_sup_sat, 5)) < eps[init_sup_sat]['nwtrph_prec'])
