import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn

from numpy import array as arr_t, frombuffer, repeat, zeros, float64, ones, roll, pi, allclose, copy

from math import exp, log, sqrt, pi

def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev  = 1.4
  n_tot  = 60e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

def test(turb_cond):
  print 'turb_cond = ', turb_cond
  opts_init = lgrngn.opts_init_t()
  kappa = .61
  opts_init.dry_distros = {kappa:lognormal}
  opts_init.coal_switch=0
  opts_init.sedi_switch=0
  opts_init.dt = 1
  opts_init.sd_conc = 64
  opts_init.n_sd_max = 512
  opts_init.rng_seed = 396
  opts_init.exact_sstp_cond = True # test would fail with per-cell sstp logic
  opts_init.turb_cond_switch = turb_cond
  spinup = 20
  
  backend = lgrngn.backend_t.serial
  
  opts = lgrngn.opts_t()
  opts.sedi=0
  opts.coal=0
  opts.cond=1
  opts.turb_cond=turb_cond
  
  opts_init.nx = 2
  opts_init.dx = 1
  opts_init.x1 = opts_init.nx * opts_init.dx
  
  opts_init.nz = 1
  opts_init.dz = 1
  opts_init.z1 = opts_init.nz * opts_init.dz

  opts_init.SGS_mix_len = opts_init.dz * ones(opts_init.nz) 
  
  rhod =   1. * ones((opts_init.nx, opts_init.nz))
  Cx =   1. * ones((opts_init.nx+1+4, opts_init.nz))
  Cz =   0. * ones((opts_init.nx+4, opts_init.nz+1))
  eps =   1e-4 * ones((opts_init.nx, opts_init.nz))
  
  for sstp_cond in [1,2,5]:
    print 'sstp_cond = ' + str(sstp_cond)
    opts.adve=0
    opts_init.sstp_cond = sstp_cond
    prtcls = lgrngn.factory(backend, opts_init)
    th =   300 * ones((opts_init.nx, opts_init.nz))
    rv =   .0025 * ones((opts_init.nx, opts_init.nz))
    rv[1,0] = 0.0095
    prtcls.init(th, rv, rhod, Cx=Cx, Cz=Cz)
  
    #equilibrium wet moment post spinup
    prtcls.diag_all()
    prtcls.diag_wet_mom(3);
    wet_post_init = copy(frombuffer(prtcls.outbuf()).reshape(opts_init.nx, opts_init.nz))
    water_post_init = 1000. * 4./3. * pi * wet_post_init + rv
    
    #spinup to get equilibrium
    if(turb_cond):
      prtcls.step_sync(opts, th, rv, diss_rate=eps)
    else:
      prtcls.step_sync(opts, th, rv)
    for it in range(spinup):
      prtcls.step_async(opts)
      if(turb_cond):
        prtcls.step_sync(opts, th, rv, diss_rate=eps)
      else:
        prtcls.step_sync(opts, th, rv)
    
    #equilibrium wet moment post spinup
    prtcls.diag_all()
    prtcls.diag_wet_mom(3);
    wet_post_spin = copy(frombuffer(prtcls.outbuf()).reshape(opts_init.nx, opts_init.nz))
    water_post_spin = 1000. * 4./3. * pi * wet_post_spin + rv
    assert allclose(water_post_spin, water_post_init, atol=0, rtol=1e-10) #some discrepancy due to water density
    
    #advect SDs
    opts.adve=1
    prtcls.step_async(opts)
    
    #equilibrium wet moment post adve
    prtcls.diag_all()
    prtcls.diag_wet_mom(3);
    wet_post_adve = copy(frombuffer(prtcls.outbuf()).reshape(opts_init.nx, opts_init.nz))
    wet_post_adve_roll = roll(wet_post_adve,1).copy()
    assert all(wet_post_adve_roll == wet_post_spin)
    
    #advect rv
    tmp = rv.copy()
    rv[0,0] = tmp[1,0]
    rv[1,0] = tmp[0,0]
    #advect th
    tmp = th.copy()
    th[0,0] = tmp[1,0]
    th[1,0] = tmp[0,0]
  
    #condensation with advected SDs and rv
    if(turb_cond):
      prtcls.step_sync(opts, th, rv, diss_rate=eps)
    else:
      prtcls.step_sync(opts, th, rv)
    
    #wet mom post adve and cond
    prtcls.diag_all()
    prtcls.diag_wet_mom(3);
    wet_post_adve_cond =  copy(frombuffer(prtcls.outbuf()).reshape(opts_init.nx, opts_init.nz))
    assert allclose(wet_post_adve, wet_post_adve_cond, atol=0, rtol=3e-2)

test(False)
test(True)
