import sys
sys.path.insert(0, "../../bindings/python/")

from numpy import array as arr_t # ndarray dtype default to float64, while array's is int64!
from numpy import arange
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


opts = lgrngn.opts_t()

opts_init = lgrngn.opts_init_t()
kappa = .61
opts_init.dry_distros = {kappa:lognormal}
opts_init.coal_switch = False
opts_init.sedi_switch = False
opts_init.RH_max = 0.999 # to comply with the assert(RH<1) at init
opts_init.dt = 1
opts_init.sd_conc = int(1e2)
opts_init.n_sd_max = opts_init.sd_conc

backend = lgrngn.backend_t.serial

opts.adve = False
opts.sedi = False
opts.cond = True
opts.coal = False
opts.chem = False

#expected theta and rv after condensation:
exp_th = { True : 306.9, # constp
           False: 307.78}  # varp
exp_rv = { True : 1.628e-2, # constp
           False: 1.7e-2}  # varp

def supersaturation(prtcls):
    prtcls.diag_RH()
    return (frombuffer(prtcls.outbuf())[0] - 1) * 100

def temperature(prtcls):
    prtcls.diag_temperature()
    return frombuffer(prtcls.outbuf())[0]

def pressure(prtcls):
    prtcls.diag_pressure()
    return frombuffer(prtcls.outbuf())[0]

def initial_state():
    rhod = arr_t([1.  ])
    th   = arr_t([300.])
    rv   = arr_t([0.02])

    T = common.T(th[0], rhod[0])
    p = arr_t([common.p(rhod[0], rv[0], T)])

    return rhod, th, rv, p

def test(RH_formula, step_count, substep_count, exact_substep, constp):
    print("[RH_formula = ", RH_formula,"]")
    print("step_count = ", step_count, " substep_count = ", substep_count, "exact substepping = ", exact_substep, "constp = ", constp)

    opts_init.sstp_cond=substep_count
    opts_init.exact_sstp_cond=exact_substep
    opts_init.RH_formula = RH_formula

    rhod, th, rv, p = initial_state()
    rv_init = rv.copy()

    # in constp mode, th_std is expected instead of th_dry
    if constp == True:
      # dry/std conversions assume p = rhod (Rd + rv * Rv) T
      # which in general is not true in constp, but is true at init so we use it here
      th[0] = common.th_dry2std(th[0], rv[0])

    th_init = th.copy()
    prtcls = lgrngn.factory(backend, opts_init)
    if constp == False:
      prtcls.init(th, rv, rhod)
    else:
      prtcls.init(th, rv, rhod, p)
    ss = supersaturation(prtcls)
    print("initial supersaturation", ss)

    exectime = 0
    # first step without condesnation just to see diag output
    opts.cond = False
    for step in arange(step_count):
      wrapped = wrapper(prtcls.step_sync, opts, th, rv, rhod)
      exectime += timeit.timeit(wrapped, number=1)
      prtcls.step_async(opts)
      opts.cond = True
  #    print step, supersaturation(prtcls), temperature(prtcls), pressure(prtcls), th[0], rv[0]

    ss_post_cond = supersaturation(prtcls)
    print("supersaturation after condensation", ss_post_cond, th[0], rv[0])

    assert(abs(th[0] - exp_th[constp]) < 1e-4 * exp_th[constp])
    assert(abs(rv[0] - exp_rv[constp]) < 1e-3 * exp_rv[constp])
    rv_diff = rv_init.copy() - rv[0].copy()
  
    # change to subsaturated air - test evaporation
    rv[0]   = 0.002
    rv_init = rv.copy()

    for step in arange(step_count):
      wrapped = wrapper(prtcls.step_sync, opts, th, rv, rhod)
      exectime += timeit.timeit(wrapped, number=1)
      prtcls.step_async(opts)
#      print step, supersaturation(prtcls), temperature(prtcls), pressure(prtcls), th[0], rv[0]

    ss_post_evap = supersaturation(prtcls)
    print("supersaturation after evaporation", ss_post_evap, th[0], rv[0])
    print('execution time: ', exectime)
    
    return ss_post_cond, th[0] - th_init[0], rv[0] - rv_init[0] - rv_diff[0]


for constp in [False, True]:
  for exact_sstp in [False, True]:
    for RH_formula in [lgrngn.RH_formula_t.pv_cc, lgrngn.RH_formula_t.rv_cc, lgrngn.RH_formula_t.pv_tet, lgrngn.RH_formula_t.rv_tet]:
      ss, th_diff_1  , rv_diff = test(RH_formula, 40, 1, exact_sstp, constp)
      print(ss, th_diff_1  , rv_diff)
      assert(abs(ss) < 4.5e-3)
      assert(abs(rv_diff) < 1e-9)

      ss, th_diff_10 , rv_diff = test(RH_formula, 40, 10, exact_sstp, constp)
      print(ss, th_diff_10 , rv_diff)
      assert(abs(ss) < 4.5e-3)
      assert(abs(rv_diff) < 1e-9)

      ss, th_diff_100, rv_diff = test(RH_formula, 40, 100, exact_sstp, constp)
      print(ss, th_diff_100, rv_diff)
      assert(abs(ss) < 4.5e-3)
      assert(abs(rv_diff) < 1e-9)

      if constp == False:
        assert(abs(th_diff_1) < 4.2e-2)
        assert(abs(th_diff_10) < 4.2e-3)
        assert(abs(th_diff_100) < 4.2e-4)
      else :
        # TODO: why with constant pressure the error doesn't scale so well?
        #       is there a systematic error caused by the fact that with constant pressure,
        #       pressure doesnt agree with T, rv and rhod?
        assert(abs(th_diff_1) < 1.1e-1)
        assert(abs(th_diff_10) < 7.4e-2)
        assert(abs(th_diff_100) < 7.3e-2)



