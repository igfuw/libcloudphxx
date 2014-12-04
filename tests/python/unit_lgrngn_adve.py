import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn

import numpy as np 

from math import exp, log, sqrt, pi

def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev  = 1.4
  n_tot  = 60e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

opts_init = lgrngn.opts_init_t()
kappa = .61
opts_init.dry_distros = {kappa:lognormal}

opts_init.dt = 1

opts_init.nz = 2
opts_init.nx = 3
opts_init.dz = 1
opts_init.dx = 1
opts_init.z1 = opts_init.nz * opts_init.dz
opts_init.x1 = opts_init.nx * opts_init.dx

opts_init.sd_conc_mean = 11. / (opts_init.nx * opts_init.nz)

backend = lgrngn.backend_t.serial

opts = lgrngn.opts_t()
opts.adve = True
opts.sedi = False
opts.cond = False
opts.coal = False
opts.chem = False

rhod = 1. * np.ones((opts_init.nx, opts_init.nz))
th   = 300. * np.ones((opts_init.nx, opts_init.nz))
rv   = 0.01 * np.ones((opts_init.nx, opts_init.nz))

# a 2x2 grid
#
#      |      |
#
#  ->  *  ->  *  ->   
# 
#      |      |
#      
#  ->  *  ->  *  ->
# 
#      |      |

rhod_Cx = 1. * np.ones((opts_init.nx + 1, opts_init.nz))
rhod_Cz = 0. * np.ones((opts_init.nx, opts_init.nz + 1))

print rhod.shape, rhod_Cx.shape

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod, rhod_Cx, rhod_Cz)

prtcls.step_sync(opts, th, rv, rhod)

#prtcls.diag_wet_rng(0,1)
prtcls.diag_sd_conc()
print np.frombuffer(prtcls.outbuf()).reshape(opts_init.nx, opts_init.nz)

prtcls.step_async(opts)
prtcls.step_sync(opts, th, rv, rhod)

#prtcls.diag_wet_rng(0,1)
prtcls.diag_sd_conc()
print np.frombuffer(prtcls.outbuf()).reshape(opts_init.nx, opts_init.nz)
