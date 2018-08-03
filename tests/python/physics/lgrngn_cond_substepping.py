# test condensation with th, rv, rhod substepping
# and low RH_max to see how substepping affects activation
# it's with GCCNs - large differences between Tetens and Clau-Clap in this case!
#                   also makes substepping important!

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


def lognormal2(lnr):
  mean_r = 4e-6 / 2
  stdev  = 1.2
  n_tot  = 10e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

opts = lgrngn.opts_t()

opts_init = lgrngn.opts_init_t()
opts_init.dry_distros = {.61:lognormal, 1.28:lognormal2} # normal mode + GCCNs
opts_init.coal_switch = False
opts_init.sedi_switch = False
opts_init.RH_max = 1.0001
opts_init.dt = 1
opts_init.sd_conc = int(1e3)
opts_init.n_sd_max = opts_init.sd_conc

backend = lgrngn.backend_t.serial

opts.adve = False
opts.sedi = False
opts.cond = True
opts.coal = False
opts.chem = False
#opts.RH_max = 1.005

#expected theta and rv after condensation:
exp_th = { True : 298.884, # constp
           False: 300.115}   # varp
exp_rv = { True : 9.05e-3, # constp
           False: 9.06e-3}  # varp
exp_rv_diff = { True : 1e-6, # constp
                False: 5e-9}   # varp

# range of supersat after condesation
ss_max = -0.71
ss_min = -0.96


def supersaturation(prtcls):
    prtcls.diag_RH()
    return (frombuffer(prtcls.outbuf())[0] - 1) * 100

def mean_r(prtcls):
    prtcls.diag_wet_rng(0.5e-6, 1)
    prtcls.diag_wet_mom(1)
    mom1 = frombuffer(prtcls.outbuf())[0]
    prtcls.diag_wet_mom(0)
    mom0 = frombuffer(prtcls.outbuf())[0]
    return mom1 / mom0 * 1e6

def second_r(prtcls):
    prtcls.diag_wet_rng(0.5e-6, 1)
    prtcls.diag_wet_mom(2)
    mom2 = frombuffer(prtcls.outbuf())[0]
    prtcls.diag_wet_mom(0)
    mom0 = frombuffer(prtcls.outbuf())[0]
    return mom2 / mom0

def third_r(prtcls):
    prtcls.diag_wet_rng(0.5e-6, 1)
    prtcls.diag_wet_mom(3)
    mom3 = frombuffer(prtcls.outbuf())[0]
    prtcls.diag_wet_mom(0)
    mom0 = frombuffer(prtcls.outbuf())[0]
    return mom3 / mom0

def act_conc(prtcls):
    prtcls.diag_wet_rng(0.5e-6, 1)
    prtcls.diag_wet_mom(0)
    return (frombuffer(prtcls.outbuf())[0] / 1e3) # per gram

def gccn_conc(prtcls):
    prtcls.diag_dry_rng(0.5e-6, 1)
    prtcls.diag_wet_mom(0)
    return (frombuffer(prtcls.outbuf())[0] / 1e3) # per gram

def initial_state():
    # a little below saturation
    rhod = arr_t([1.1  ])
    th   = arr_t([305.])
    rv   = arr_t([0.0085])

    T = common.T(th[0], rhod[0])
    p = arr_t([common.p(rhod[0], rv[0], T)])

    return rhod, th, rv, p

def supersat_state():
    rhod = arr_t([1.  ])
    th   = arr_t([300.])
    rv   = arr_t([0.0091])

    T = common.T(th[0], rhod[0])
    p = arr_t([common.p(rhod[0], rv[0], T)])

    return rhod, th, rv, p

def test(RH_formula, step_count, substep_count, exact_substep, constp):
    print "[RH_formula = ", RH_formula,"]"
    print "step_count = ", step_count, " substep_count = ", substep_count, "exact substepping = ", exact_substep, "constp = ", constp

    opts_init.sstp_cond=substep_count
    opts_init.exact_sstp_cond=exact_substep
    opts_init.RH_formula = RH_formula

    rhod, th, rv, p = initial_state()
    rhod_ss, th_ss, rv_ss, p_ss = supersat_state()

    # in constp mode, th_std is expected instead of th_dry
    if constp == True:
      # dry/std conversions assume p = rhod (Rd + rv * Rv) T
      # which in general is not true in constp, but is true at init so we use it here
      th[0] = common.th_dry2std(th[0], rv[0])
      th_ss[0] = common.th_dry2std(th_ss[0], rv_ss[0])

    prtcls = lgrngn.factory(backend, opts_init)

    if constp == False:
      prtcls.init(th, rv, rhod)
    else:
      prtcls.init(th, rv, rhod, p_ss)

    # go to supersaturated air, density changes to test density substepping too
    rhod[0] = rhod_ss
    th[0]   = th_ss
    rv[0]   = rv_ss
    rv_init = rv.copy()
    th_init = th.copy()

    exectime = 0
    opts.cond = 0
    for step in arange(step_count):
      wrapped = wrapper(prtcls.step_sync, opts, th, rv, rhod)
      exectime += timeit.timeit(wrapped, number=1)
      prtcls.step_async(opts)
      if step == 9:
        # some parameters are analyzed after 10 steps, before small CCNs evaporate
        act_conc_post_cond = act_conc(prtcls)
        mean_r_post_cond = mean_r(prtcls)
        second_r_post_cond = second_r(prtcls)
        third_r_post_cond = third_r(prtcls)
      if step == 0:
        print "initial supersaturation", supersaturation(prtcls)
        opts.cond = 1
#      print step, supersaturation(prtcls), th[0], rv[0], mean_r(prtcls), second_r(prtcls), third_r(prtcls), act_conc(prtcls)
    
    ss_post_cond = supersaturation(prtcls)
    print "supersaturation after condensation", ss_post_cond, th[0], rv[0], mean_r(prtcls), second_r(prtcls), third_r(prtcls), act_conc(prtcls)

    assert(abs(th[0] - exp_th[constp]) < 1e-4 * exp_th[constp])
    assert(abs(rv[0] - exp_rv[constp]) < 1e-3 * exp_rv[constp])
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
 #     print step, supersaturation(prtcls), th[0], rv[0], mean_r(prtcls), second_r(prtcls), third_r(prtcls), act_conc(prtcls)

    ss_post_evap = supersaturation(prtcls)
    print "supersaturation after evaporation", ss_post_evap, th[0], rv[0], mean_r(prtcls), second_r(prtcls), third_r(prtcls), act_conc(prtcls), gccn_conc(prtcls)
    # after evaporation, only larger mode particles should have r > 0.5 microns
    assert(act_conc(prtcls) == gccn_conc(prtcls))
    print 'execution time: ', exectime
    
    return ss_post_cond, th[0] - th_init[0] - th_diff[0], rv[0] - rv_init[0] - rv_diff[0], act_conc_post_cond, mean_r_post_cond, second_r_post_cond, third_r_post_cond


for constp in [False, True]:
  for exact_sstp in [False, True]:
    for RH_formula in [lgrngn.RH_formula_t.pv_cc, lgrngn.RH_formula_t.rv_cc, lgrngn.RH_formula_t.pv_tet, lgrngn.RH_formula_t.rv_tet]:

      ss, th_diff_1  , rv_diff, act, mr, sr, tr = test(RH_formula, 100, 1, exact_sstp, constp)
      print ss, th_diff_1  , rv_diff, act, mr, sr, tr
      assert(ss_min < ss < ss_max) # GCCNs condensate even at ss<0
      assert(abs(rv_diff) < exp_rv_diff[constp])
      assert(abs(th_diff_1) < 7e-3)

      # expected concentration of droplets with r>0.5um
      exp_act ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 40809,
                          lgrngn.RH_formula_t.rv_cc  : 40809,
                          lgrngn.RH_formula_t.pv_tet : 40809,
                          lgrngn.RH_formula_t.rv_tet : 8128,
                       },
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 40809,
                          lgrngn.RH_formula_t.rv_cc  : 40809,
                          lgrngn.RH_formula_t.pv_tet : 40809,
                          lgrngn.RH_formula_t.rv_tet : 8125,
                        }
               }
      # expected mean radius of droplets with r>0.5um
      exp_mr ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 3.141,
                          lgrngn.RH_formula_t.rv_cc  : 3.156,
                          lgrngn.RH_formula_t.pv_tet : 2.767,
                          lgrngn.RH_formula_t.rv_tet : 8.170,
                       },
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 2.905,
                          lgrngn.RH_formula_t.rv_cc  : 2.92,
                          lgrngn.RH_formula_t.pv_tet : 2.417,
                          lgrngn.RH_formula_t.rv_tet : 8.099,
                        }
               }

      # expected second moment of droplets with r>0.5um
      exp_sr ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 1.662e-11,
                          lgrngn.RH_formula_t.rv_cc  : 1.668e-11,
                          lgrngn.RH_formula_t.pv_tet : 1.528e-11,
                          lgrngn.RH_formula_t.rv_tet : 6.760e-11,
                       }, 
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 1.559e-11,
                          lgrngn.RH_formula_t.rv_cc  : 1.564e-11,
                          lgrngn.RH_formula_t.pv_tet : 1.421e-11,
                          lgrngn.RH_formula_t.rv_tet : 6.641e-11,
                        } 
               }

      # expected third moment of droplets with r>0.5um
      exp_tr ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 1.226e-16,
                          lgrngn.RH_formula_t.rv_cc  : 1.228e-16,
                          lgrngn.RH_formula_t.pv_tet : 1.179e-16,
                          lgrngn.RH_formula_t.rv_tet : 5.662e-16,
                       }, 
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 1.173e-16,
                          lgrngn.RH_formula_t.rv_cc  : 1.174e-16,
                          lgrngn.RH_formula_t.pv_tet : 1.131e-16,
                          lgrngn.RH_formula_t.rv_tet : 5.513e-16,
                        } 
               }

      assert(abs(act - exp_act[constp][RH_formula]) < 1e-3 * exp_act[constp][RH_formula])
      assert(abs(mr - exp_mr[constp][RH_formula]) < 1e-3 * exp_mr[constp][RH_formula])
      assert(abs(sr - exp_sr[constp][RH_formula]) < 1e-3 * exp_sr[constp][RH_formula])
      assert(abs(tr - exp_tr[constp][RH_formula]) < 2e-3 * exp_tr[constp][RH_formula])


      ss, th_diff_10  , rv_diff, act, mr, sr, tr = test(RH_formula, 100, 10, exact_sstp, constp)
      print ss, th_diff_10  , rv_diff, act, mr, sr, tr
      assert(ss_min < ss < ss_max) # GCCNs condensate even at ss<0
      assert(abs(rv_diff) < exp_rv_diff[constp])
      assert(abs(th_diff_10) < 5e-3)


      exp_act ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 16561,
                          lgrngn.RH_formula_t.rv_cc  : 17036,
                          lgrngn.RH_formula_t.pv_tet : 9485,
                          lgrngn.RH_formula_t.rv_tet : 8125,
                       }, 
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 14931,
                          lgrngn.RH_formula_t.rv_cc  : 15502,
                          lgrngn.RH_formula_t.pv_tet : 8131,
                          lgrngn.RH_formula_t.rv_tet : 8125,
                        } 
               }

      exp_mr ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 4.726,
                          lgrngn.RH_formula_t.rv_cc  : 4.642,
                          lgrngn.RH_formula_t.pv_tet : 7.200,
                          lgrngn.RH_formula_t.rv_tet : 8.212,
                       }, 
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 4.878,
                          lgrngn.RH_formula_t.rv_cc  : 4.743,
                          lgrngn.RH_formula_t.pv_tet : 8.191,
                          lgrngn.RH_formula_t.rv_tet : 8.133,
                        } 
               }

      exp_sr ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 3.526e-11,
                          lgrngn.RH_formula_t.rv_cc  : 3.435e-11,
                          lgrngn.RH_formula_t.pv_tet : 5.961e-11,
                          lgrngn.RH_formula_t.rv_tet : 6.828e-11,
                       }, 
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 3.78e-11,
                          lgrngn.RH_formula_t.rv_cc  : 3.646e-11,
                          lgrngn.RH_formula_t.pv_tet : 6.8e-11,
                          lgrngn.RH_formula_t.rv_tet : 6.699e-11,
                        } 
               }

      exp_tr ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 2.948e-16,
                          lgrngn.RH_formula_t.rv_cc  : 2.867e-16,
                          lgrngn.RH_formula_t.pv_tet : 5.054e-16,
                          lgrngn.RH_formula_t.rv_tet : 5.749e-16,
                       }, 
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 3.166e-16,
                          lgrngn.RH_formula_t.rv_cc  : 3.051e-16,
                          lgrngn.RH_formula_t.pv_tet : 5.72e-16,
                          lgrngn.RH_formula_t.rv_tet : 5.585e-16,
                        } 
               }

      assert(abs(act - exp_act[constp][RH_formula]) < 1.5e-2 * exp_act[constp][RH_formula])
      assert(abs(mr - exp_mr[constp][RH_formula]) < 1.5e-2 * exp_mr[constp][RH_formula])
      assert(abs(sr - exp_sr[constp][RH_formula]) < 1.5e-2 * exp_sr[constp][RH_formula])
      assert(abs(tr - exp_tr[constp][RH_formula]) < 1.5e-2 * exp_tr[constp][RH_formula])

      ss, th_diff_100  , rv_diff, act, mr, sr, tr = test(RH_formula, 100, 100, exact_sstp, constp)
      print ss, th_diff_100  , rv_diff, act, mr, sr, tr
      assert(ss_min < ss < ss_max) # GCCNs condensate even at ss<0
      assert(abs(rv_diff) < exp_rv_diff[constp])
      assert(abs(th_diff_100) < 5e-3)

      exp_act ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 16097,
                          lgrngn.RH_formula_t.rv_cc  : 16404,
                          lgrngn.RH_formula_t.pv_tet : 8979,
                          lgrngn.RH_formula_t.rv_tet : 8125,
                       }, 
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 13040,
                          lgrngn.RH_formula_t.rv_cc  : 13751,
                          lgrngn.RH_formula_t.pv_tet : 8123,
                          lgrngn.RH_formula_t.rv_tet : 8125,
                        } 
               }

      exp_mr ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 4.792,
                          lgrngn.RH_formula_t.rv_cc  : 4.741,
                          lgrngn.RH_formula_t.pv_tet : 7.572,
                          lgrngn.RH_formula_t.rv_tet : 8.222,
                       }, 
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 5.438,
                          lgrngn.RH_formula_t.rv_cc  : 5.203,
                          lgrngn.RH_formula_t.pv_tet : 8.2095,
                          lgrngn.RH_formula_t.rv_tet : 8.145,
                        } 
               }

      exp_sr ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 3.622e-11,
                          lgrngn.RH_formula_t.rv_cc  : 3.560e-11,
                          lgrngn.RH_formula_t.pv_tet : 6.307e-11,
                          lgrngn.RH_formula_t.rv_tet : 6.844e-11,
                       }, 
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 4.324e-11,
                          lgrngn.RH_formula_t.rv_cc  : 4.105e-11,
                          lgrngn.RH_formula_t.pv_tet : 6.8239e-11,
                          lgrngn.RH_formula_t.rv_tet : 6.717e-11,
                        } 
               }

      exp_tr ={ True: {        # constp
                          lgrngn.RH_formula_t.pv_cc  : 3.041e-16,
                          lgrngn.RH_formula_t.rv_cc  : 2.985e-16,
                          lgrngn.RH_formula_t.pv_tet : 5.357e-16,
                          lgrngn.RH_formula_t.rv_tet : 5.768e-16,
                       }, 
                 False: {        # varp
                          lgrngn.RH_formula_t.pv_cc  : 3.64e-16,
                          lgrngn.RH_formula_t.rv_cc  : 3.452e-16,
                          lgrngn.RH_formula_t.pv_tet : 5.7412e-16,
                          lgrngn.RH_formula_t.rv_tet : 5.609e-16,
                        } 
               }

      # reduced precision due to differences in results between Linux and OSX
      assert(abs(act - exp_act[constp][RH_formula]) < 1.5e-2 * exp_act[constp][RH_formula])
      assert(abs(mr - exp_mr[constp][RH_formula]) < 1.5e-2 * exp_mr[constp][RH_formula])
      assert(abs(sr - exp_sr[constp][RH_formula]) < 1.5e-2 * exp_sr[constp][RH_formula])
      assert(abs(tr - exp_tr[constp][RH_formula]) < 1.5e-2 * exp_tr[constp][RH_formula])




