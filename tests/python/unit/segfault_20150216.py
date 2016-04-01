import sys, boost.mpi
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn
from math import exp, log, sqrt, pi
import numpy as np 

opts_init = lgrngn.opts_init_t()
opts_init.coal_switch = False
opts_init.sedi_switch = False

opts_init.dt = 1

def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev = 1.4
  n_tot = 60e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);
kappa = .61

opts_init.dry_distros = {kappa:lognormal}
opts_init.sd_conc = 50
opts_init.n_sd_max = 50

try:
  prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init) # the segfault actually was detected on CUDA, but we don't have yet a mechanism do detect if CUDA hardware is present (e.g. on Travis)
except:
  prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)

rhod = 1. * np.ones((1,))
th = 300. * np.ones((1,))
rv = 0.01 * np.ones((1,))

prtcls.init(th, rv, rhod)
