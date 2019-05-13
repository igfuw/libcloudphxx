import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn

import numpy as np 
from math import exp, log, sqrt, pi
from time import time

def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev  = 1.4
  n_tot  = 60e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

Opts_init = lgrngn.opts_init_t()
kappa = .61
Opts_init.dry_distros = {kappa:lognormal}
Opts_init.coal_switch = False
Opts_init.sedi_switch = True
Opts_init.terminal_velocity = lgrngn.vt_t.beard76

Opts_init.dt = 0.01

Opts_init.nz = 6
Opts_init.nx = 1
Opts_init.dz = 1
Opts_init.dx = 1
Opts_init.z1 = Opts_init.nz * Opts_init.dz
Opts_init.x1 = Opts_init.nx * Opts_init.dx

Opts_init.rng_seed = int(time())
Opts_init.sd_conc = 1000
Opts_init.n_sd_max = Opts_init.sd_conc * (Opts_init.nx * Opts_init.nz)
Opts_init.div_LS = 1. # 1/s large-scale subsidence

Backend = lgrngn.backend_t.serial

Opts = lgrngn.opts_t()
Opts.adve = False
Opts.sedi = True
Opts.cond = False
Opts.coal = False
Opts.chem = False
Opts.rcyc = False

Rhod =   1. * np.ones((Opts_init.nx, Opts_init.nz))
Th   = 300. * np.ones((Opts_init.nx, Opts_init.nz))
Rv   = 0.01 * np.ones((Opts_init.nx, Opts_init.nz))

prtcls = lgrngn.factory(Backend, Opts_init)
prtcls.init(Th, Rv, Rhod)

prtcls.diag_all()
prtcls.diag_sd_conc()
tab_in = np.copy(np.frombuffer(prtcls.outbuf()).reshape(Opts_init.nx, Opts_init.nz))
print "at init \n", tab_in

for it in range(100):
  prtcls.step_sync(Opts, Th, Rv, Rhod)
  prtcls.step_async(Opts)

prtcls.diag_all()
prtcls.diag_sd_conc()
tab_out = np.copy(np.frombuffer(prtcls.outbuf()).reshape(Opts_init.nx, Opts_init.nz))
print "after 1s \n", tab_out

# SD position z(t) = z_0 exp(-div_LS t)
e = 2.7183
res_01 = Opts_init.sd_conc * e  # expected number of SDs in cells 0 and 1 at t=1s
res_2 = Opts_init.sd_conc * (6. - 2. * e) # expected number of SDs in cell 2 at t=1s

print "analytical result in cells 0 and 1: ",res_01
print "analytical result in cell 2: ",res_2

assert(tab_out[0][5] == 0.)
assert(tab_out[0][4] == 0.)
assert(tab_out[0][3] == 0.)
tolerance = 1. / sqrt(Opts_init.sd_conc)
print "relative tolerance: ", tolerance
assert np.isclose(res_2, tab_out[0][2], atol=0., rtol=5*tolerance)
assert np.isclose(res_01, tab_out[0][1], atol=0., rtol=tolerance)
assert np.isclose(res_01, tab_out[0][0], atol=0., rtol=tolerance)

