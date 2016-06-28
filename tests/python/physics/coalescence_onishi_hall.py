# test coalescence algorithm for Onishi kernel with Hall efficiencies
# by comparing speedup due to turbulence, see Onishi (2015) in JAS

import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn

import numpy as np
import time

def spherevol(r):
  return 4./3.*pow(r,3)*np.pi;

def diag_rain_volume():
  prtcls.diag_wet_rng(40e-6, 1)
  prtcls.diag_wet_mom(3)
  return 4./3.*np.pi*np.frombuffer(prtcls.outbuf()).mean() * rhod.mean()

#Onishi case
r_zero = 15e-6
n_zero = 1.42e8

epsilon=0.04 # turbulent dissipation rate
re_lambda=66 # Taylor-microscale Reynolds number

n_runs = 50

def expvolumelnr(lnr):
  r=np.exp(lnr)
  return n_zero * 3.*np.power(r,3)/np.power(r_zero,3)*np.exp(- np.power((r/r_zero),3));

opts_init = lgrngn.opts_init_t()

opts_init.dt = 1.
opts_init.terminal_velocity = lgrngn.vt_t.beard77fast
opts_init.sd_conc = 1024
opts_init.n_sd_max = 1024

t10_arr = np.zeros((2, n_runs))
v_zero = spherevol(r_zero)

th   = 300 * np.ones((1,))       
rv   = 1. * np.ones((1,))
rhod  = 1.22419 * np.ones((1,))

opts_init.dry_distros = {0.:expvolumelnr}

Opts = lgrngn.opts_t()
Opts.adve = False
Opts.sedi = False
Opts.cond = False
Opts.coal = True
Opts.chem = False
  
for kernel_t in range(2):
  if(kernel_t == 0):
    opts_init.kernel = lgrngn.kernel_t.hall
  else:
    opts_init.kernel = lgrngn.kernel_t.onishi_hall
    opts_init.kernel_parameters = np.array([epsilon, re_lambda])

  for zz in range(n_runs):
    opts_init.rng_seed = int(time.time())
    
    try:
      prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
    except:
      prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
    
    prtcls.init(th, rv, rhod)
  
    prtcls.diag_all()
    prtcls.diag_wet_mom(3)
    total_volume=4./3.*np.pi*np.frombuffer(prtcls.outbuf()).mean() * rhod.mean()
    
    rain_volume = 0.
    t=0
    t10=0 
    while(t10 == 0):
      prtcls.step_sync(Opts,th,rv,rhod)
      prtcls.step_async(Opts)
      rain_volume = diag_rain_volume()
      t+=opts_init.dt
      if(rain_volume > total_volume/10.):
        t10=t
        t10_arr[kernel_t][zz] = t10

hall_mean = t10_arr[0].mean()
onishi_mean = t10_arr[1].mean()

print('Hall: Time to turn 10% of water into rain drops (r>40um) '+str(t10_arr[0].mean()) + ' +/- ' +str(np.std(t10_arr[0])))
print('Onishi: Time to turn 10% of water into rain drops (r>40um) '+str(t10_arr[1].mean()) + ' +/- ' +str(np.std(t10_arr[1])))
if(hall_mean / onishi_mean > 1.62 or hall_mean / onishi_mean < 1.22):
  raise Exception("Onishi turbulent enhancement not corresct, ratio of t10% times: "+str(hall_mean / onishi_mean))

