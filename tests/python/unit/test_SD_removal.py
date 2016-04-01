import sys, boost.mpi
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn
from math import exp, log, sqrt, pi
import numpy as np

_Chem_g_id = { 
  "SO2_g"  : lgrngn.chem_species_t.SO2, 
  "H2O2_g" : lgrngn.chem_species_t.H2O2, 
  "O3_g"   : lgrngn.chem_species_t.O3,
  "HNO3_g" : lgrngn.chem_species_t.HNO3,
  "NH3_g"  : lgrngn.chem_species_t.NH3,
  "CO2_g"  : lgrngn.chem_species_t.CO2
}


opts_init = lgrngn.opts_init_t()
opts_init.dt = pow(2,15)
opts_init.sstp_coal = 1 

rhod = 1. * np.ones((1,))
th = 300. * np.ones((1,))
rv = 0.01 * np.ones((1,))

ambient_chem = dict((v, 1e-10 * np.ones((1,))) for k,v in _Chem_g_id.iteritems())

def expvolumelnr(lnr):  
  r_zero = 30.531e-6
  n_zero = pow(2,8)
  r=np.exp(lnr)
  return n_zero * 3.*np.power(r,3)/np.power(r_zero,3)*np.exp(- np.power((r/r_zero),3));

kappa = .01 

opts_init.dry_distros = {kappa:expvolumelnr}

opts_init.sd_conc = 64
opts_init.n_sd_max = 64
opts_init.chem_switch = True

opts_init.kernel = lgrngn.kernel_t.geometric
opts_init.terminal_velocity = lgrngn.vt_t.beard76
try:
  prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
except:
  prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)

prtcls.init(th, rv, rhod, ambient_chem = ambient_chem)

Opts = lgrngn.opts_t()
Opts.adve = False
Opts.sedi = False
Opts.cond = False
Opts.coal = True
Opts.chem_dsl = True
Opts.chem_dsc = True
Opts.chem_rct = False

for i in range(900):
  prtcls.step_sync(Opts,th,rv,rhod, ambient_chem=ambient_chem)
  prtcls.step_async(Opts)

prtcls.diag_sd_conc()
tmp = np.frombuffer(prtcls.outbuf())[0]
print 'final no of SDs: ', tmp
if(tmp >= 10 or tmp == 0):
  raise Exception("wrong amount of SDs were removed")
