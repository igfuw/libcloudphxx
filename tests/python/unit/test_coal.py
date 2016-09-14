# test if coalescence conserves total rd3 rw3 and kpa*rd3

import sys 
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn
import numpy as np

#total time of simulation
simulation_time = 200

#initial conditions, ca. 1g / m^3
r_zero = 30.084e-6
n_zero = pow(2,23)

# initial exponential distribution in droplet volume, cf. Shima 2009
# as a function of ln(r)
def expvolumelnr(lnr):
  r=np.exp(lnr)
  return n_zero * 3.*np.power(r,3)/np.power(r_zero,3)*np.exp(- np.power((r/r_zero),3));

opts_init = lgrngn.opts_init_t()
opts_init.dt = simulation_time
opts_init.sstp_coal = simulation_time

rhod = 1. * np.ones((1,))
th = 300. * np.ones((1,))
rv = 0.01 * np.ones((1,))

kappa1 = 0.1
kappa2 = 0.9

opts_init.dry_distros = {kappa1:expvolumelnr, kappa2:expvolumelnr}

opts_init.sd_conc = pow(2,14)
opts_init.n_sd_max = pow(2,14)

opts_init.kernel = lgrngn.kernel_t.geometric
opts_init.terminal_velocity = lgrngn.vt_t.beard77fast

try:
  prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
except:
  prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)

prtcls.init(th, rv, rhod)

opts = lgrngn.opts_t()
opts.adve = False
opts.sedi = False
opts.cond = False
opts.coal = True
opts.chem = False
opts.rcyc = True

def totrd3():
  prtcls.diag_all()
  prtcls.diag_dry_mom(3)
  return np.frombuffer(prtcls.outbuf())[0]

def totrw3():
  prtcls.diag_all()
  prtcls.diag_wet_mom(3)
  return np.frombuffer(prtcls.outbuf())[0]

bins = pow(10, -6 + np.arange(200)/50.)

#get sum(kappa*rd^3)
def kappa_rd3_sum():
  ret=0
  for i in range(bins.size - 1) :
    prtcls.diag_dry_rng(bins[i], bins[i+1])
    prtcls.diag_kappa(1)
    ret += np.frombuffer(prtcls.outbuf())[0] * pow((bins[i+1] + bins[i]) / 2, 3)
  return ret


def diag(arg):
  arg[0] = totrd3()
  arg[1] = totrw3()
  arg[2] = kappa_rd3_sum()
  
res_init = np.zeros(3)
diag(res_init)

#simulation loop
prtcls.step_sync(opts, th, rv, rhod)
prtcls.step_async(opts)

res_final = np.zeros(3)
diag(res_final)

eps = 1e-10
assert np.isclose(res_final[0], res_init[0], atol=0., rtol=eps),\
  "total dry volume is not conserved during coalescence"
assert np.isclose(res_final[1], res_init[1], atol=0., rtol=eps),\
  "total wet volume is not conserved during coalescence"
eps = 1e-2 # looser condition due to discretization into bins
assert np.isclose(res_final[2], res_init[2], atol=0., rtol=eps),\
  "total kappa*rd^3 is not conserved during coalescence"

