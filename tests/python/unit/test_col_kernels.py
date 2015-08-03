import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn
from math import exp, log, sqrt, pi
import numpy as np

rhod = 1. * np.ones((1,))
th = 300. * np.ones((1,))
rv = 0.01 * np.ones((1,))

def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev = 1.4
  n_tot = 60e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

kappa = .61

count = 0
for kernel in [lgrngn.kernel_t.geometric, lgrngn.kernel_t.geometric, lgrngn.kernel_t.long,  lgrngn.kernel_t.hall, lgrngn.kernel_t.hall_davis_no_waals, lgrngn.kernel_t.golovin]:
  opts_init = lgrngn.opts_init_t()
  opts_init.dt = 1
  opts_init.dry_distros = {kappa:lognormal}
  opts_init.sd_conc_mean = 50.
  opts_init.terminal_velocity=lgrngn.vt_t.beard
  opts_init.kernel = kernel
  opts_init.kernel_parameters = np.array([])
  if(kernel == lgrngn.kernel_t.golovin):
    opts_init.kernel_parameters = np.array([1.])
  if(kernel == lgrngn.kernel_t.geometric):
    if(count == 0):
    # geometric
      count+=1
    else:
    # geometric with multiplier
      opts_init.kernel_parameters = np.array([10.])

  try:
    prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
  except:
    prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)

  prtcls.init(th, rv, rhod)

  Opts = lgrngn.opts_t()
  Opts.adve = False
  Opts.sedi = False
  Opts.cond = False
  Opts.coal = True
  Opts.chem_dsl = False
  Opts.chem_dsc = False
  Opts.chem_rct = False

  prtcls.step_sync(Opts,th,rv,rhod)
  prtcls.step_async(Opts)
