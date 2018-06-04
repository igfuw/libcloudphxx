import sys
sys.path.insert(0, "../../bindings/python/")
sys.path.insert(0, "../../../build/bindings/python/")

from numpy import array as arr_t # ndarray dtype default to float64, while array's is int64!
from numpy import arange
from numpy import isclose
from numpy import frombuffer
from math import exp, log, sqrt, pi
import timeit

from libcloudphxx import lgrngn
from libcloudphxx import common

# wrapper for timing excecution time
def wrapper(func, opts, th, rv, rhod):
    def wrapped():
        return func(opts, th, rv, rhod)
    return wrapped

def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev  = 1.4
  n_tot  = 60e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);


def lognormal2(lnr):
  mean_r = .4e-6 / 2
  stdev  = 1.2
  n_tot  = 20e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

opts = lgrngn.opts_t()

opts_init = lgrngn.opts_init_t()
kappa = .61
opts_init.dry_distros = {kappa:lognormal, kappa:lognormal2}
opts_init.coal_switch = False
opts_init.sedi_switch = False
opts_init.RH_max = 1.01
opts_init.dt = 1
opts_init.sd_conc = int(1e3)
opts_init.n_sd_max = opts_init.sd_conc

backend = lgrngn.backend_t.serial

opts.adve = False
opts.sedi = False
opts.cond = True
opts.coal = False
opts.chem = False
opts.RH_max = 1.01

#expected theta and rv after condensation:
exp_th = { True : 309.357, # constp
           False: 307.798}   # varp
exp_rv = { True : 1.64e-2, # constp
           False: 1.7e-2}  # varp


def supersaturation(prtcls):
    prtcls.diag_RH()
    return (frombuffer(prtcls.outbuf())[0] - 1) * 100

def mean_r(prtcls):
    prtcls.diag_all()
    prtcls.diag_wet_mom(1)
    mom1 = frombuffer(prtcls.outbuf())[0]
    prtcls.diag_wet_mom(0)
    mom0 = frombuffer(prtcls.outbuf())[0]
    return mom1 / mom0 * 1e6

def act_conc(prtcls):
    prtcls.diag_rw_ge_rc()
    prtcls.diag_wet_mom(0)
    return (frombuffer(prtcls.outbuf())[0]) * 1e6

def initial_state():
    # a little below saturation
    rhod = arr_t([1.1  ])
    th   = arr_t([305.])
    rv   = arr_t([0.0085])
    dt   = 1
    
    T = common.T(th[0], rhod[0])
    p = arr_t([common.p(rhod[0], rv[0], T)])

    return rhod, th, rv, p, dt

def test(RH_formula, step_count, substep_count, exact_substep, constp):
    print "[RH_formula = ", RH_formula,"]"
    print "step_count = ", step_count, " substep_count = ", substep_count, "exact substepping = ", exact_substep, "constp = ", constp

    opts_init.sstp_cond=substep_count
    opts_init.exact_sstp_cond=exact_substep
    opts_init.RH_formula = lgrngn.RH_formula_t.pv_cc

    rhod, th, rv, p, dt = initial_state()
    prtcls = lgrngn.factory(backend, opts_init)
    if constp == False:
      prtcls.init(th, rv, rhod)
    else:
      prtcls.init(th, rv, rhod, p)
    ss = supersaturation(prtcls)
    print "initial supersaturation", ss

    # go to supersaturated air (same conditions as in sat_adj_blk_1m test), density changes to test density substepping too
    rhod[0] = 1.1
    th[0]   = 300.
    rv[0]   = 0.02
    rv_init = rv.copy()
    th_init = th.copy()

    exectime = 0
    for step in arange(step_count):
      wrapped = wrapper(prtcls.step_sync, opts, th, rv, rhod)
      exectime += timeit.timeit(wrapped, number=1)
      prtcls.step_async(opts)
      print step, supersaturation(prtcls), th[0], rv[0], mean_r(prtcls), act_conc(prtcls)
    
    ss_post_cond = supersaturation(prtcls)
    print "supersaturation after condensation", ss_post_cond, th[0], rv[0], mean_r(prtcls), act_conc(prtcls)

#    assert(abs(th[0] - exp_th[constp]) < 1e-5 * exp_th[constp])
#    assert(abs(rv[0] - exp_rv[constp]) < 1e-3 * exp_rv[constp])
    rv_diff = rv_init.copy() - rv[0].copy()
    th_diff = th_init.copy() - th[0].copy()
  
    # change to subsaturated air - test evaporation
    rhod[0] = 1.1
    th[0]   = 305
    rv[0]   = 0.0085
    rv_init = rv.copy()
    th_init = th.copy()

    for step in arange(step_count):
      wrapped = wrapper(prtcls.step_sync, opts, th, rv, rhod)
      exectime += timeit.timeit(wrapped, number=1)
      prtcls.step_async(opts)
      print step, supersaturation(prtcls), th[0], rv[0], mean_r(prtcls), act_conc(prtcls)

    ss_post_evap = supersaturation(prtcls)
    print "supersaturation after evaporation", ss_post_evap, th[0], rv[0], mean_r(prtcls), act_conc(prtcls)
    print 'execution time: ', exectime
    
    return ss_post_cond, th[0] - th_init[0] - th_diff[0], rv[0] - rv_init[0] - rv_diff[0]


for constp in [False, True]:
  for exact_sstp in [False]:#, [True, False]:
    for RH_formula in [lgrngn.RH_formula_t.pv_cc]:#, lgrngn.RH_formula_t.rv_cc, lgrngn.RH_formula_t.pv_tet, lgrngn.RH_formula_t.rv_tet]:
      ss, th_diff_1  , rv_diff = test(RH_formula, 400, 1, exact_sstp, constp) 
      print ss, th_diff_1  , rv_diff
#      assert(abs(ss) < 4.5e-3)
#      assert(abs(rv_diff) < 1e-9)
#      assert(abs(th_diff_1) < 4.2e-2)

      ss, th_diff_10 , rv_diff = test(RH_formula, 400, 10, exact_sstp, constp)
      print ss, th_diff_10 , rv_diff
#      assert(abs(ss) < 4.5e-3)
#      assert(abs(rv_diff) < 1e-9)
#      assert(abs(th_diff_10) < 4.2e-3)

      ss, th_diff_100, rv_diff = test(RH_formula, 400, 100, exact_sstp, constp)
      print ss, th_diff_100, rv_diff
#      assert(abs(ss) < 4.5e-3)
#      assert(abs(rv_diff) < 1e-9)
#      assert(abs(th_diff_100) < 4.2e-4)




