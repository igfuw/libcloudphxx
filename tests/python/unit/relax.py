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

def lognormal_rlx(lnr):
  mean_r = .04e-6 / 2
  stdev  = 1.4
  n_tot  = 120e6 
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
  opts_init.supstp_rlx = 2
  opts_init.rng_seed = int(time())
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

  opts_init.y0=0.;
  opts_init.y1=1.;
  
  opts_init.chem_switch = 0;
  opts_init.coal_switch = 0;
  opts_init.adve_switch = 0;
  opts_init.cond_switch = 0;
  opts_init.sedi_switch = 0;
  opts_init.src_switch = 0;
  opts_init.rlx_switch = 1;
  
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
  
  # 2 steps during which relaxation should be done once
  opts.rlx = 1
  for i in range(2):
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
opts_init.rlx_dry_distros = {kappa: [lognormal_rlx, [0,2],[0,opts_init.dz]]}
opts_init.sd_conc = 1024
opts_init.rlx_bins = 1024
opts_init.rlx_timescale = 4 # whole simulation time is 2, so this means we should get half of the droplets added




print(' --- dry_distros rlx sd_per_bin = 1 ---')

opts_init.rlx_sd_per_bin = 1
opts_init.n_sd_max = int((opts_init.sd_conc * 2 + opts_init.rlx_bins * opts_init.rlx_sd_per_bin * 2) * 2) # assuming nx=nz=2

sd_conc, wet_mom0, wet_mom1 = test(opts_init)

print('diag_sd_conc', sd_conc)
# relaxation should add SD only in the lower cells
# there is some randomness in the number of SD that should be added
exp_sd_min = 1024 + opts_init.rlx_sd_per_bin * 400
exp_sd_max = 1024 + opts_init.rlx_sd_per_bin * 600
if sd_conc[0] < exp_sd_min or sd_conc[0] > exp_sd_max or sd_conc[2] < exp_sd_min or sd_conc[2] > exp_sd_max:
  raise Exception("wrong amount of SDs were added")
if not(sd_conc[1] == 1024 and sd_conc[3] == 1024):
  raise Exception("SDs were added in wrong cells")

print(('wet mom0', wet_mom0))
print((wet_mom0[0] + wet_mom0[2]) / (wet_mom0[1] + wet_mom0[3])) # relax n_stp is two times bigger
if abs((wet_mom0[0] + wet_mom0[2]) / (wet_mom0[1] + wet_mom0[3]) - 1.5)  > 0.01:  # we expect 1.5, because n_rlx is two tmes bigger, but relaxation time scale is twice the simulation time
  raise Exception("incorrect multiplicity after source")

print(('wet mom1', wet_mom1))
print((wet_mom1[0] + wet_mom1[2]) / (wet_mom1[1] + wet_mom1[3])) # relax n_stp is two times bigger
if abs((wet_mom1[0] + wet_mom1[2]) / (wet_mom1[1] + wet_mom1[3]) - 1.5)  > 0.01:
  raise Exception("incorrect radius after source")




print(' --- dry_distros rlx sd_per_bin = 10 ---')

opts_init.rlx_sd_per_bin = 10
opts_init.n_sd_max = int((opts_init.sd_conc * 2 + opts_init.rlx_bins * opts_init.rlx_sd_per_bin * 2) * 2) # assuming nx=nz=2

sd_conc, wet_mom0, wet_mom1 = test(opts_init)

print('diag_sd_conc', sd_conc)
# relaxation should add SD only in the lower cells
# there is some randomness in the number of SD that should be added
exp_sd_min = 1024 + opts_init.rlx_sd_per_bin * 400
exp_sd_max = 1024 + opts_init.rlx_sd_per_bin * 600
if sd_conc[0] < exp_sd_min or sd_conc[0] > exp_sd_max or sd_conc[2] < exp_sd_min or sd_conc[2] > exp_sd_max:
  raise Exception("wrong amount of SDs were added")
if not(sd_conc[1] == 1024 and sd_conc[3] == 1024):
  raise Exception("SDs were added in wrong cells")

print(('wet mom0', wet_mom0))
print((wet_mom0[0] + wet_mom0[2]) / (wet_mom0[1] + wet_mom0[3])) # relax n_stp is two times bigger
if abs((wet_mom0[0] + wet_mom0[2]) / (wet_mom0[1] + wet_mom0[3]) - 1.5)  > 0.01:
  raise Exception("incorrect multiplicity after source")

print(('wet mom1', wet_mom1))
print((wet_mom1[0] + wet_mom1[2]) / (wet_mom1[1] + wet_mom1[3])) # relax n_stp is two times bigger
if abs((wet_mom1[0] + wet_mom1[2]) / (wet_mom1[1] + wet_mom1[3]) - 1.5)  > 0.01:
  raise Exception("incorrect radius after source")

