import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn

from numpy import array as arr_t, frombuffer, repeat, zeros, float64, ones, roll, pi

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
opts_init.coal_switch=0
opts_init.sedi_switch=0
opts_init.dt = 1
opts_init.sd_conc = 64
opts_init.n_sd_max = 512
opts_init.rng_seed = 396
spinup = 100

backend = lgrngn.backend_t.serial

opts = lgrngn.opts_t()
opts.sedi=0
opts.coal=0
opts.cond=1

# 1D (periodic horizontal domain)
rhod = arr_t([  1.,   1.])
C    = arr_t([   1.,   1.,  1.])

opts_init.nx = 2
opts_init.dx = 1
opts_init.x1 = opts_init.nx * opts_init.dx

print 'spinup: ' + str(spinup)

for sstp_cond in [1,2,20]:
  print ''
  print 'sstp_cond = ' + str(sstp_cond)
  opts.adve=0
  opts_init.sstp_cond = sstp_cond
  prtcls = lgrngn.factory(backend, opts_init)
  th   = arr_t([300., 300.])
  rv   = arr_t([   .009,  .1]) # first cell subsaturated, second cell supersaturated
  prtcls.init(th, rv, rhod, C)

  #equilibrium wet moment post spinup
  prtcls.diag_all()
  prtcls.diag_wet_mom(3);
  wet_post_init = frombuffer(prtcls.outbuf()).copy()
  print 'wet post init: ' + str(wet_post_init)
  print 'water post init: ' + str(1000. * 4./3. * pi * wet_post_init + rv)
  
  #spinup to get equilibrium
#  opts.RH_max = 1.1
  for it in range(spinup):
    prtcls.step_sync(opts, th, rv)
    prtcls.step_async(opts)
  
#  opts.RH_max = 44
  #equilibrium wet moment post spinup
  prtcls.diag_all()
  prtcls.diag_wet_mom(3);
  wet_post_spin = frombuffer(prtcls.outbuf()).copy()
  print 'wet post spin: ' + str(wet_post_spin)
  print 'water post spin: ' + str(1000. * 4./3. * pi * wet_post_spin + rv)
  
  #one more cond to make sure we reached equilibrium
  print th
  print rv
  prtcls.step_sync(opts, th, rv)
  print th
  print rv
  
  #equilibrium wet moment post spinup + one cond
  prtcls.diag_all()
  prtcls.diag_wet_mom(3);
  wet_post_spinpp = frombuffer(prtcls.outbuf()).copy()
  print 'wet post spin and one cond: ' + str(wet_post_spinpp) 
  print 'prcntg chagne: ' + str((wet_post_spinpp - wet_post_spin) / wet_post_spin)
  
  #advect SDs
  opts.adve=1
  prtcls.step_async(opts)
  
  #equilibrium wet moment post adve
  prtcls.diag_all()
  prtcls.diag_wet_mom(3);
  wet_post_adve = frombuffer(prtcls.outbuf()).copy()
  print 'wet post adve: ' + str(wet_post_adve) 
  wet_post_adve_roll = roll(wet_post_adve,1).copy()
  print 'wet post adve: ' + str(wet_post_adve_roll)
  print 'prcntg chagne (rolled): ' + str((wet_post_adve_roll - wet_post_spinpp) / wet_post_spinpp)
  
  #advect rv
  tmp = rv.copy()
  rv[0] = tmp[1]
  rv[1] = tmp[0]
  #advect th
  tmp = th.copy()
  th[0] = tmp[1]
  th[1] = tmp[0]

  #condensation with advected SDs and rv
  print rv
  prtcls.step_sync(opts, th, rv)
  
  #equilibrium wet mom post adve and cond
  prtcls.diag_all()
  prtcls.diag_wet_mom(3);
  wet_post_adve_cond = frombuffer(prtcls.outbuf()).copy()
  print 'wet post adve and cond: ' + str(wet_post_adve_cond) 
  print 'prcntg chagne: ' + str((wet_post_adve_cond - wet_post_adve) / wet_post_adve)
