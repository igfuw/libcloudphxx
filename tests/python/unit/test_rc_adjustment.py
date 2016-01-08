#test if after initialization we have approximately the same water content in each cell

import sys 
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn
import numpy as np
import time
from math import pi

# initial exponential distribution in droplet volume
# as a function of ln(r)
def expvolumelnr(lnr):
  r=np.exp(lnr)
  return n_zero * 3.*np.power(r,3)/np.power(r_zero,3)*np.exp(- np.power((r/r_zero),3));

#initial conditions, ca. 2g / m^3
r_zero = 15e-6
n_zero = 1.42e8

opts_init = lgrngn.opts_init_t()
opts_init.coal_switch=0
opts_init.sedi_switch=0
opts_init.cond_switch=0
opts_init.adve_switch=0
opts_init.dt = 1
opts_init.dx = 1
opts_init.dz = 3
opts_init.dy = 2
opts_init.nx = 2 
opts_init.nz = 3 
opts_init.ny = 4 
opts_init.x1 = opts_init.dx * opts_init.nx
opts_init.z1 = opts_init.dz * opts_init.nz
opts_init.y1 = opts_init.dy * opts_init.ny
opts_init.rng_seed = int(time.time())

th   = 300 * np.ones((opts_init.nx, opts_init.ny, opts_init.nz))  
rv   = 0.01 * np.ones((opts_init.nx, opts_init.ny, opts_init.nz))  
rc_adj   = 0.001 * np.ones((opts_init.nx, opts_init.ny, opts_init.nz))  
rhod = 1. * np.ones((opts_init.nx, opts_init.ny, opts_init.nz)) + .1 * np.mgrid[1:1+opts_init.nx, 1:1+opts_init.ny, 1:1+opts_init.nz][1] # different densities, hence different water content

kappa = 1e-6

opts_init.dry_distros = {kappa:expvolumelnr}

opts_init.sd_conc = 64
opts_init.n_sd_max = opts_init.sd_conc * opts_init.nx * opts_init.ny * opts_init.nz

try:
  prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
except:
  prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)

prtcls.init(th, rv, rhod)

prtcls.diag_all()
prtcls.diag_wet_mom(3) # gives specific moment (divided by rhod)
water_pre =1e3* 4/3*pi* np.frombuffer(prtcls.outbuf())#.mean() # dropping a constant
print water_pre

prtcls.step_rc_adjust(rc_adj)

prtcls.diag_wet_mom(3) # gives specific moment (divided by rhod)
water_post = 1e3*4/3*pi*np.frombuffer(prtcls.outbuf())#.mean() # dropping a constant
water_post -= rc_adj[0][0][0] * np.ones((opts_init.nx * opts_init.ny * opts_init.nz)) # 
print water_post
if(not np.allclose(water_pre, water_post,rtol=0,atol=1e-8)):
  raise Exception("Error in adding rc adjustment")
