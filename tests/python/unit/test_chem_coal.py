import sys 
sys.path.insert(0, "../../bindings/python/")
sys.path.insert(0, "/home/pdziekan/praca/libcloudphxx/build/bindings/python/")

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

def expvolumelnr(lnr):  
  r_zero = 30.531e-6
  n_zero = pow(2,8)
  r=np.exp(lnr)
  return n_zero * 3.*np.power(r,3)/np.power(r_zero,3)*np.exp(- np.power((r/r_zero),3));

kappa = 0.

opts_init.dry_distros = {kappa:expvolumelnr}

opts_init.sd_conc = 64
opts_init.n_sd_max = 64

opts_init.chem_switch = True
opts_init.chem_rho = 1.8e-3
ambient_chem = dict((v, np.ones((1,)) ) for k,v in _Chem_g_id.iteritems())

opts_init.kernel = lgrngn.kernel_t.geometric
opts_init.terminal_velocity = lgrngn.vt_t.beard
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
Opts.chem_dsl = False
Opts.chem_dsc = False
Opts.chem_rct = False

prtcls.diag_all()
prtcls.diag_chem(lgrngn.chem_species_t.NH4);
NH4_init = np.frombuffer(prtcls.outbuf())[0]
prtcls.diag_chem(lgrngn.chem_species_t.H);
H_init = np.frombuffer(prtcls.outbuf())[0]
prtcls.diag_chem(lgrngn.chem_species_t.S_VI);
S_VI_init = np.frombuffer(prtcls.outbuf())[0]

for i in range(300):
  prtcls.step_sync(Opts,th,rv,rhod, ambient_chem=ambient_chem)
  prtcls.step_async(Opts)

prtcls.diag_all()
prtcls.diag_chem(lgrngn.chem_species_t.NH4);
NH4_final = np.frombuffer(prtcls.outbuf())[0]
prtcls.diag_chem(lgrngn.chem_species_t.H);
H_final = np.frombuffer(prtcls.outbuf())[0]
prtcls.diag_chem(lgrngn.chem_species_t.S_VI);
S_VI_final = np.frombuffer(prtcls.outbuf())[0]

eps=1e-10

assert np.isclose(NH4_init, NH4_final, atol=0., rtol=eps),\
  "total amount of NH4 is not conserved during coalescence"

assert np.isclose(H_init, H_final, atol=0., rtol=eps),\
  "total amount of H is not conserved during coalescence"

assert np.isclose(S_VI_init, S_VI_final, atol=0., rtol=eps),\
  "total amount of S_VI is not conserved during coalescence"
