import sys
#try:
#  import boost.mpi
#except:
#  pass

sys.path.insert(0, "../../bindings/python/")
sys.path.insert(0, "../../../build/bindings/python/")

from libcloudphxx import lgrngn

from numpy import array as arr_t, frombuffer

from math import exp, log, sqrt, pi
from time import time

def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev  = 1.4
  n_tot  = 60e6 
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

#def lognormal_rlx(lnr):
#  mean_r = .10e-6 / 2
#  stdev  = 1.4
#  n_tot  = 60e4
#  return n_tot * exp(
#    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
#  ) / log(stdev) / sqrt(2*pi);

def test(opts_init):
  opts_init.supstp_rlx = 50
  opts_init.rng_seed = int(time())
#  opts_init.rng_seed_init = int(time())
  print('rng_seed:', opts_init.rng_seed)
#  print('rng_seed_init:', opts_init.rng_seed_init)
  opts_init.dt = 1
  opts_init.nx = 2;
  opts_init.nz = 2;
  opts_init.dx=1.;
  opts_init.dz=1.;
  opts_init.x0=0.;
  opts_init.z0=0.;
  opts_init.x1=opts_init.nx * opts_init.dx;
  opts_init.z1=opts_init.nz * opts_init.dz;
  opts_init.aerosol_independent_of_rhod=1;
  
  opts_init.chem_switch = 0;
  opts_init.coal_switch = 0;
  opts_init.adve_switch = 0;
  opts_init.cond_switch = 0;
  opts_init.sedi_switch = 0;
  opts_init.src_switch = 0;
  opts_init.rlx_switch = 1;

  opts_init.src_z0 = 0;
  opts_init.src_z1 = opts_init.dz; #create aerosol only in the lower cells
  opts_init.src_x0 = 0;
  opts_init.src_x1 = opts_init.dx*opts_init.nx;
  
  opts = lgrngn.opts_t()
  
  opts.adve = 0;
  opts.chem = 0;
  opts.sedi = 0;
  opts.coal = 0;
  opts.cond = 0;
  
  rhod = arr_t([[  1.,    1.  ],[   1.,     1.  ]])
  th   = arr_t([[300.,  300.  ],[ 300.,   300.  ]])
  rv   = arr_t([[   .01,   .01],[    .01,    .01]])
  
  try:
    prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
  except:
    prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
  
  prtcls.init(th, rv, rhod)
  
  # 100 steps during which number of droplets should be doubled in two calls to src
  opts.rlx = 1
  for i in range(100):
    prtcls.step_sync(opts,th,rv,rhod)
    prtcls.step_async(opts)
  
  prtcls.diag_all()
  prtcls.diag_sd_conc()
  sd_conc = frombuffer(prtcls.outbuf()).copy()
  
  prtcls.diag_all()
  prtcls.diag_wet_mom(0)
  wet_mom0 = frombuffer(prtcls.outbuf()).copy()
  
  prtcls.diag_all()
  prtcls.diag_wet_mom(1)
  wet_mom1 = frombuffer(prtcls.outbuf()).copy()

  return sd_conc, wet_mom0, wet_mom1

# test source with dry_distros
kappa = .61
opts_init = lgrngn.opts_init_t()
opts_init.dry_distros = {kappa:lognormal}
opts_init.rlx_dry_distros = {kappa: [lognormal, [0,2],[0,2*opts_init.dz]]}
opts_init.sd_conc = 2#100#24
opts_init.rlx_bins = 2# 512
opts_init.n_sd_max = int((opts_init.sd_conc * 2 + opts_init.rlx_bins * 2) * 2) # assuming nx=nz=2

print(' --- dry_distros rlx ---')

sd_conc, wet_mom0, wet_mom1 = test(opts_init)

print('diag_sd_conc', sd_conc)
#if not((sd_conc[0] == 1164 or sd_conc[0] == 1165) and (sd_conc[2] == 1164 or sd_conc[2] == 1165)):
#  raise Exception("wrong amount of SDs were added")
#if not(sd_conc[1] == 1024 and sd_conc[3] == 1024):
#  raise Exception("SDs were added in wrong cells")

print(('wet mom0', wet_mom0))
if (abs( 2 - (wet_mom0[0] + wet_mom0[2]) / (wet_mom0[1] + wet_mom0[3]) ) > 0.015):
  raise Exception("incorrect multiplicity after source")

print(('wet mom1', wet_mom1))
if (abs( (7.84 / 2.12) - (wet_mom1[0] + wet_mom1[2]) / (wet_mom1[1] + wet_mom1[3]) ) > 0.015):
  raise Exception("incorrect radius after source")

