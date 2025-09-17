import sys
sys.path.insert(0, "../../bindings/python/")
from timeit import default_timer as timer


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
    th_dry = arr_t([300.]) 
    if init_sup_sat:
        rv   = arr_t([0.02])
    else:
        rv   = arr_t([0.002])
    rc   = arr_t([0.015])
    rr   = arr_t([0.  ])
    dt   = 1

    T = common.T(th_dry[0], rhod[0])
    p = arr_t([common.p(rhod[0], rv[0], T)])
    ss = supersaturation(T, p[0], rv[0])
    print("initial supersaturation", ss, "th_dry", th_dry[0], "th_std", common.th_dry2std(th_dry[0], rv[0]), "rv", rv[0], "T", T, "p", p[0], "rc", rc[0])

    return rhod, th_dry, rv, rc, rr, p, dt

def test_adj_cellwise(init_sup_sat, _opts, name, r_eps = r_eps_def, nwtrph_iters = nwtrph_iters_def):
    _opts.r_eps = r_eps
    _opts.nwtrph_iters = nwtrph_iters
    rhod, th_dry, rv, rc, rr, p, dt = initial_state(init_sup_sat)
    th_std = arr_t([common.th_dry2std(th_dry[0], rv[0])]) 

    # only two combinations allowed
    if _opts.th_dry and not _opts.const_p:
      start = timer()
      blk_1m.adj_cellwise(_opts, rhod, p, th_dry, rv, rc, rr, dt)
      end = timer()
      T = common.T(th_dry[0], rhod[0])
      p = arr_t([common.p(rhod[0], rv[0], T)])
      th_std = arr_t([common.th_dry2std(th_dry[0], rv[0])]) 
    elif not _opts.th_dry and _opts.const_p:
      start = timer()
      blk_1m.adj_cellwise(_opts, rhod, p, th_std, rv, rc, rr, dt)
      end = timer()
      T = common.exner(p[0]) * th_std[0]
      th_dry = arr_t([common.th_std2dry(th_std[0], rv[0])]) 

    ss = supersaturation(T, p[0], rv[0])
    print("final supersaturation", ss, "th_dry", th_dry[0], "th_std", th_std[0], "r_v", rv[0], "T", T, "p", p[0], "rc", rc[0])
    print(f"Execution time: {(end - start) * 1e6:.0e} microseconds")
    return ss

# def test_adj_cellwise(init_sup_sat, r_eps = r_eps_def):
#     opts.r_eps = r_eps
#     print("[standard adj_cellwise]")
#     rhod, th, rv, rc, rr, dt = initial_state(init_sup_sat)
#     blk_1m.adj_cellwise(opts, rhod, th, rv, rc, rr, dt)
#     
#     T = common.T(th[0], rhod[0])
#     p = common.p(rhod[0], rv[0], T)
#     ss = supersaturation(T, p, rv[0])
#     print("final supersaturation", ss, th[0], rv[0])
#     return ss
# 
# def test_adj_cellwise_constp(init_sup_sat, r_eps = r_eps_def):
#     opts.r_eps = r_eps
#     print("[constp adj_cellwise]")
#     rhod, th, rv, rc, rr, dt = initial_state(init_sup_sat)
# 
#     # define pressure consistent with adj_cellwise to compare results
#     p   = arr_t([common.p(rhod[0], rv[0], common.T(th[0], rhod[0]))])
# 
#     #constp requires th_std input
#     th_std = arr_t([common.th_dry2std(th[0], rv[0])])
#     blk_1m.adj_cellwise_constp(opts, rhod, p, th_std, rv, rc, rr, dt)
#     
#     T = common.exner(p[0]) * th_std[0]
#     ss = supersaturation(T, p[0], rv[0])
#     print("final supersaturation", ss, th_std[0], rv[0])
#     return ss
# 
# def test_adj_cellwise_nwtrph(init_sup_sat, nwtrph_iters = nwtrph_iters_def):
#     opts.nwtrph_iters = nwtrph_iters
#     print("[nwtrph adj_cellwise]")
#     rhod, th, rv, rc, rr, dt = initial_state(init_sup_sat)
# 
#     # define pressure consistent with adj_cellwise to compare results
#     p   = arr_t([common.p(rhod[0], rv[0], common.T(th[0], rhod[0]))])
# 
#     #nwtrph requires th_std input
#     th_std = arr_t([common.th_dry2std(th[0], rv[0])])
#     blk_1m.adj_cellwise_nwtrph(opts, p, th_std, rv, rc, dt)
#    
#     T = common.exner(p[0]) * th_std[0]
#     ss = supersaturation(T, p[0], rv[0])
#     print("final supersaturation", ss, th_std[0], rv[0])
#     return ss


test_cases = {
  "RK4, variable pressure and dry theta" : [0, 1, 0], # adj_nwtrph, th_dry, const_p
  "RK4, constant pressure and 'standard' theta" : [0, 0, 1],
  "Newton-Raphson, variable pressure and dry theta" : [1, 1, 0],
  "Newton-Raphson, constant pressure and 'standard' theta" : [1, 0, 1],
}

# precision
eps = {
        # supersaturation
        True  : {
          "RK4, variable pressure and dry theta" : 3e-2,
          "RK4, constant pressure and 'standard' theta" : 3e-2,
          "Newton-Raphson, variable pressure and dry theta" : 3, # !!
          "Newton-Raphson, constant pressure and 'standard' theta" : 1,
        },
        False  : {
          "RK4, variable pressure and dry theta" : 0.5,
          "RK4, constant pressure and 'standard' theta" : 0.5,
          "Newton-Raphson, variable pressure and dry theta" : 0.8, # !!
          "Newton-Raphson, constant pressure and 'standard' theta" : 5e-3,
        }
      }
# increased precision
eps_crank = {
        # supersaturation
        True  : {
          "RK4, variable pressure and dry theta" : 6e-4,
          "RK4, constant pressure and 'standard' theta" : 6e-4,
          "Newton-Raphson, variable pressure and dry theta" : 0.8, # !!
          "Newton-Raphson, constant pressure and 'standard' theta" : 9e-4,
        },
        False  : {
          "RK4, variable pressure and dry theta" : 2e-3,
          "RK4, constant pressure and 'standard' theta" : 2e-3,
          "Newton-Raphson, variable pressure and dry theta" : 0.8, # !!
          "Newton-Raphson, constant pressure and 'standard' theta" : 2e-6,
        }
      }

#eps = {
#        # supersaturation
#        True  : { 'org'      : 3e-2, 'constp'      : 0.1 , 'nwtrph'      : 1e-2,
#                  'org_prec' : 6e-4, 'constp_prec' : 6e-4, 'nwtrph_prec' : 9e-4 },
#        # subsaturation
#        False : { 'org'      : 0.5 , 'constp'      : 0.5 , 'nwtrph'      : 5e-3,
#                  'org_prec' : 2e-3, 'constp_prec' : 2e-3, 'nwtrph_prec' : 2e-6 }
#      }



for name, opt in test_cases.items():
  opts.adj_nwtrph = opt[0]
  opts.th_dry = opt[1]
  opts.const_p = opt[2]

  for init_sup_sat in [True, False]:
    # default precision
    if init_sup_sat:
      print("[ " + name + " default precision ] initial supersaturation")
    else:
      print("[ " + name + " default precision ] initial subsaturation")
    print("maximum final supersaturation allowed (absolute value): ", eps[init_sup_sat][name])
    assert(abs(test_adj_cellwise(init_sup_sat, opts, name)) < eps[init_sup_sat][name])

   # cranked up precision
    if init_sup_sat:
      print("[ " + name + " increased precision ] initial supersaturation")
    else:
      print("[ " + name + " increased precision ] initial subsaturation")
    print("maximum final supersaturation allowed (absolute value): ", eps_crank[init_sup_sat][name])
    if opts.adj_nwtrph:
      ss = test_adj_cellwise(init_sup_sat, opts, name, nwtrph_iters = 5)
    else:
      ss = test_adj_cellwise(init_sup_sat, opts, name, r_eps = 1e-7)
    assert(abs(ss) < eps_crank[init_sup_sat][name])

    print("")
