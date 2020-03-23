import sys
sys.path.insert(0, "../../bindings/python/")

from numpy import array as arr_t # ndarray dtype default to float64, while array's is int64!
from numpy import arange
from numpy import frombuffer
from math import exp, log, sqrt, pi

from libcloudphxx import lgrngn
from libcloudphxx import common

def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev  = 1.4
  n_tot  = 60e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

opts = lgrngn.opts_t()

opts_init = lgrngn.opts_init_t()
opts_init.dry_distros = {.61:lognormal, 1.28:lognormal}
opts_init.coal_switch = False
opts_init.sedi_switch = False
opts_init.RH_max = 0.999 # to comply with the assert(RH<1) at init
opts_init.dt = 0.1
opts_init.sd_conc = int(1e2)
opts_init.n_sd_max = opts_init.sd_conc
opts_init.diag_incloud_time = True

backend = lgrngn.backend_t.serial

opts.adve = False
opts.sedi = False
opts.cond = True
opts.coal = False
opts.chem = False

def supersaturation(prtcls):
    prtcls.diag_RH()
    return (frombuffer(prtcls.outbuf())[0] - 1) * 100

def temperature(prtcls):
    prtcls.diag_temperature()
    return frombuffer(prtcls.outbuf())[0]

def pressure(prtcls):
    prtcls.diag_pressure()
    return frombuffer(prtcls.outbuf())[0]

def diag_incloud_time(prtcls):
    prtcls.diag_incloud_time_mom(1)
    m1 = frombuffer(prtcls.outbuf())[0]
    prtcls.diag_incloud_time_mom(0)
    m0 = frombuffer(prtcls.outbuf())[0]
    return m1/m0

def initial_state():
    rhod = arr_t([1.  ])
    th   = arr_t([300.])
    rv   = arr_t([0.009 - 0.00005])
    return rhod, th, rv

rhod, th, rv = initial_state()

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

#ss = supersaturation(prtcls)
#print "initial supersaturation", ss

# first step without condesnation just to see diag output
opts.cond = True
for step in arange(400):
  rv[0] += 0.00001 * opts_init.dt
  prtcls.sync_in(th, rv, rhod)

#  ss_post_add = supersaturation(prtcls)
#  print "supersaturation after increaseing rv", ss_post_add

  prtcls.step_cond(opts, th, rv)
  prtcls.step_async(opts)

#  ss_post_cond = supersaturation(prtcls)
#  print "supersaturation after condensation", ss_post_cond

prtcls.diag_all()
all_mean_incloud_time = diag_incloud_time(prtcls)
print "all: mean incloud time:", all_mean_incloud_time
prtcls.diag_dry_rng(0, 0.02e-6)
rdlt2_mean_incloud_time = diag_incloud_time(prtcls)
print "rd < 2um: mean incloud time:", rdlt2_mean_incloud_time
prtcls.diag_dry_rng(0.02e-6, 1)
rdgt2_mean_incloud_time = diag_incloud_time(prtcls)
print "rd > 2um: mean incloud time:", rdgt2_mean_incloud_time
prtcls.diag_dry_rng(0.02e-6, 1)
prtcls.diag_kappa_rng_cons(1, 10)
rdgt2_kpagt1_mean_incloud_time = diag_incloud_time(prtcls)
print "rd > 2um and kappa > 1: mean incloud time:", rdgt2_kpagt1_mean_incloud_time
prtcls.diag_dry_rng(0.02e-6, 1)
prtcls.diag_kappa_rng_cons(0, 1)
rdgt2_kpalt1_mean_incloud_time = diag_incloud_time(prtcls)
print "rd > 2um and kappa < 1: mean incloud time:", rdgt2_kpalt1_mean_incloud_time

assert(rdlt2_mean_incloud_time < all_mean_incloud_time)
assert(all_mean_incloud_time < rdgt2_mean_incloud_time)
assert(rdgt2_mean_incloud_time < rdgt2_kpagt1_mean_incloud_time)
assert(rdgt2_kpalt1_mean_incloud_time < rdgt2_mean_incloud_time)

