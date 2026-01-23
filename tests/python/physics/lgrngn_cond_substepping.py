# test condensation with th, rv, rhod substepping
# and low RH_max to see how substepping affects activation
# it's with GCCNs - large differences between Tetens and Clau-Clap in this case!
#                   also makes substepping important!
#
# Usage:
#   python lgrngn_cond_substepping.py              # Run tests and save results
#   python lgrngn_cond_substepping.py --save-ref   # Also save results as reference data

import sys
import os
sys.path.insert(0, "../../../build/bindings/python/")
sys.path.insert(0, "../../bindings/python/")

from numpy import array as arr_t # ndarray dtype default to float64, while array's is int64!
from numpy import arange
from numpy import frombuffer
from math import exp, log, sqrt, pi
import timeit

from libcloudphxx import lgrngn
from libcloudphxx import common

import pandas as pd

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
opts_init.dry_distros = {(.61,0.):lognormal, (1.28,0.):lognormal2} # normal mode + GCCNs
opts_init.coal_switch = False
opts_init.sedi_switch = False
# opts_init.RH_max = 1.0001
opts_init.RH_max = 0.95
opts_init.ice_switch = False
opts_init.dt = 1
opts_init.sd_conc = int(1e3)
opts_init.n_sd_max = opts_init.sd_conc

opts_init.rc2_T = 10 # results are the same for 0C to 100C
opts_init.sstp_cond_adapt_drw2_eps = 1e-3 #1e-3
opts_init.sstp_cond_adapt_drw2_max = 2 #2


# backend = lgrngn.backend_t.CUDA
backend = lgrngn.backend_t.OpenMP
# backend = lgrngn.backend_t.serial

opts.adve = False
opts.sedi = False
opts.cond = True
opts.coal = False
opts.chem = False
opts.RH_max = 1.005
# opts.RH_max = 0.9
opts.ice_nucl = False
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

def test(RH_formula, step_count, substep_count, exact_substep, constp, mixing, adaptive, sstp_cond_act):
    print("[RH_formula = ", RH_formula,"]", flush=True)
    print("step_count = ", step_count, " substep_count = ", substep_count, "exact substepping = ", exact_substep, "constp = ", constp, "mixing per substep = ", mixing, "adaptive substep no. = ", adaptive, "sstp_cond_act = ", sstp_cond_act, flush=True)

    opts_init.sstp_cond=substep_count
    opts_init.exact_sstp_cond=exact_substep
    opts_init.RH_formula = RH_formula
    opts_init.sstp_cond_mix = mixing
    opts_init.adaptive_sstp_cond=adaptive
    opts_init.sstp_cond_act=sstp_cond_act

    # opts_init.rd_min=1e-10
    # opts_init.rd_max=1e-5

    rhod, th, rv, p = initial_state()
    rhod_ss, th_ss, rv_ss, p_ss = supersat_state()

    # in constp mode, th_std is expected instead of th_dry
    if constp == True:
      # dry/std conversions assume p = rhod (Rd + rv * Rv) T
      # which in general is not true in constp, but is true at init so we use it here
      th[0] = common.th_dry2std(th[0], rv[0])
      th_ss[0] = common.th_dry2std(th_ss[0], rv_ss[0])
      opts_init.const_p = True
      opts_init.th_dry = False
    else:
      opts_init.const_p = False
      opts_init.th_dry = True

    prtcls = lgrngn.factory(backend, opts_init)

    if constp == False:
      prtcls.init(th, rv, rhod)
    else:
      prtcls.init(th, rv, rhod, p_ss)

    # go to supersaturated air, density changes to test density substepping too
    rhod[0] = rhod_ss[0]
    th[0]   = th_ss[0]
    rv[0]   = rv_ss[0]
    rv_init = rv.copy()
    th_init = th.copy()

    exectime = 0
    opts.cond = 0
    for step in arange(step_count):
      # print("step ", step, flush=True)
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
        print("initial supersaturation", supersaturation(prtcls))
        opts.cond = 1
    
    ss_post_cond = supersaturation(prtcls)
    print("supersaturation after condensation", ss_post_cond, th[0], rv[0], mean_r(prtcls), second_r(prtcls), third_r(prtcls), act_conc(prtcls))

    th_post_cond = th[0]
    rv_post_cond = rv[0]

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

    ss_post_evap = supersaturation(prtcls)
    print("supersaturation after evaporation", ss_post_evap, th[0], rv[0], mean_r(prtcls), second_r(prtcls), third_r(prtcls), act_conc(prtcls), gccn_conc(prtcls))
    print('execution time: ', exectime)

    results = {
      'ss': ss_post_cond,
      'th_diff': th[0] - th_init[0] - th_diff[0],
      'rv_diff': rv[0] - rv_init[0] - rv_diff[0],
      'act': act_conc_post_cond,
      'mr': mean_r_post_cond,
      'sr': second_r_post_cond,
      'tr': third_r_post_cond,
      'exectime': exectime,     
      'act_post_evap': act_conc(prtcls),
      'gccn_post_evap': gccn_conc(prtcls),
      'th_post_cond': th_post_cond,
      'rv_post_cond': rv_post_cond,
    }
    
    return results

records = []

for adaptive in [True, False]: # adaptive condensation substepping?
  for mixing in [False, True]: # communicate changes in rv an theta between SDs after each substep?
    for constp in [False, True]:
      for exact_sstp in [False, True]:
        # for RH_formula in [lgrngn.RH_formula_t.pv_cc]:
        for RH_formula in [lgrngn.RH_formula_t.pv_cc, lgrngn.RH_formula_t.rv_cc, lgrngn.RH_formula_t.pv_tet, lgrngn.RH_formula_t.rv_tet]:
          for sstp_cond in [1, 2, 3, 4, 6, 8, 32]:
            for sstp_cond_act in [1, 8]:
              if(mixing == False and exact_sstp == False):
                continue # mixing can be turned off only with exact substepping
              if(exact_sstp == False and adaptive == True):
                continue # adaptive substepping requires exact substepping
              if(sstp_cond_act > 1 and adaptive == False):
                continue # sstp_cond_act > 1 requires adaptive substepping
              if(adaptive == True and exact_sstp == True and mixing == True):
                continue # adaptive perparticle doesnt work with mixing yet
              
              results = test(RH_formula, 100, sstp_cond, exact_sstp, constp, mixing, adaptive, sstp_cond_act)
              print(results)
              results['mixing'] = mixing
              results['adaptive'] = adaptive
              results['constp'] = constp
              results['exact_sstp'] = exact_sstp
              results['RH_formula'] = RH_formula
              results['sstp_cond'] = sstp_cond
              results['sstp_cond_act'] = sstp_cond_act
              results['sstp_cond_adapt_drw2_eps'] = opts_init.sstp_cond_adapt_drw2_eps
              results['sstp_cond_adapt_drw2_max'] = opts_init.sstp_cond_adapt_drw2_max
              records.append(results)

# save results to a CSV file for refdata comparison and for plotting
df = pd.DataFrame(records)
df['sd_conc'] = opts_init.sd_conc  # Add the column to all rows at once
df['RH_max'] = opts_init.RH_max
df['dt'] = 1
os.makedirs("test_results", exist_ok=True)
df.to_csv("test_results/lgrngn_cond_substepping_results.csv", index=False)

# Optionally save as reference data
if '--save-ref' in sys.argv:
  print("\nSaving results as reference data...")
# Get the directory where this script is located
  script_dir = os.path.dirname(os.path.abspath(__file__))    
  refdata_dir = os.path.join(script_dir, "refdata")
  os.makedirs(refdata_dir, exist_ok=True)
  refdata_file = os.path.join(refdata_dir, "lgrngn_cond_substepping_refdata.csv")
  df.to_csv(refdata_file, index=False)
  print("Reference data saved to: test_results/lgrngn_cond_substepping_refdata.csv")
  print("Future runs can be compared against this reference using:")
  print("  python lgrngn_cond_substepping_test.py")

