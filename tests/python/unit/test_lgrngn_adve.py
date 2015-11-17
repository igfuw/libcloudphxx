import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn

import numpy as np 
from math import exp, log, sqrt, pi

import pytest

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
Opts_init.sedi_switch = False

Opts_init.dt = 1

Opts_init.nz = 5
Opts_init.nx = 6
Opts_init.dz = 1
Opts_init.dx = 1
Opts_init.z1 = Opts_init.nz * Opts_init.dz
Opts_init.x1 = Opts_init.nx * Opts_init.dx

Opts_init.sd_conc = 50 / (Opts_init.nx * Opts_init.nz)

Backend = lgrngn.backend_t.serial

Opts = lgrngn.opts_t()
Opts.adve = True
Opts.sedi = False
Opts.cond = False
Opts.coal = False
Opts.chem = False

Rhod =   1. * np.ones((Opts_init.nx, Opts_init.nz))
Th   = 300. * np.ones((Opts_init.nx, Opts_init.nz))
Rv   = 0.01 * np.ones((Opts_init.nx, Opts_init.nz))

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



def advection_1step(Cx_arg, Cz_arg, backend=Backend, opts_init=Opts_init, opts=Opts, 
                    rhod=Rhod, th=Th, rv=Rv):
  prtcls = lgrngn.factory(backend, opts_init)
  Cx = Cx_arg * np.ones((opts_init.nx + 1, opts_init.nz))
  Cz = Cz_arg * np.ones((opts_init.nx, opts_init.nz + 1))
  prtcls.init(th, rv, rhod, Cx, Cz)

  prtcls.step_sync(opts, th, rv, rhod)

  #prtcls.diag_wet_rng(0,1)
  prtcls.diag_sd_conc()
  tab_in = np.copy(np.frombuffer(prtcls.outbuf()).reshape(opts_init.nx, opts_init.nz))
  print "tab_in \n", tab_in
  
  prtcls.step_async(opts)
  prtcls.step_sync(opts, th, rv, rhod)
  
  #prtcls.diag_wet_rng(0,1)
  prtcls.diag_sd_conc()
  tab_out = np.frombuffer(prtcls.outbuf()).reshape(opts_init.nx, opts_init.nz)
  print "tab_out \n", tab_out, np.roll(tab_out, -1, 0)
  return tab_in, tab_out


@pytest.mark.parametrize("Cx, Cz, roll_st, roll_ax", [
                          (1., 0., -1, 0), (-1., 0., 1, 0),
                          pytest.mark.xfail((0., 1., -1, 1)), 
                          pytest.mark.xfail((0., -1., 1, 1))
                          ])
def test_advection(Cx, Cz, roll_st, roll_ax):
  tab_in, tab_out = advection_1step(Cx, Cz)
  print "w tescie \n", tab_in, "\n", tab_out
  assert (tab_in == np.roll(tab_out, roll_st, roll_ax)).all()
