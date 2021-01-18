import sys
#try:
#  import boost.mpi
#except:
#  pass

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

ambient_chem = dict((v, 1e-10 * np.ones((1,))) for k,v in _Chem_g_id.items())

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
opts_init.sedi_switch = False

opts_init.kernel = lgrngn.kernel_t.geometric
for vt_t in [lgrngn.vt_t.beard76, lgrngn.vt_t.beard77, lgrngn.vt_t.beard77fast, lgrngn.vt_t.khvorostyanov_spherical, lgrngn.vt_t.khvorostyanov_nonspherical]:
  print("vt_t: ", vt_t)
  opts_init.terminal_velocity = vt_t
  #try:
  #  prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
  #except:
  #  prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
  prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
  
  prtcls.init(th, rv, rhod, ambient_chem = ambient_chem)
  
  Opts = lgrngn.opts_t()
  Opts.adve = False
  Opts.sedi = False
  Opts.cond = False
  Opts.coal = True
  Opts.chem_dsl = True
  Opts.chem_dsc = False
  Opts.chem_rct = False
  Opts.rcyc = True
  
  print('before loop')
  for i in range(900):
#    print 'start of step no ', i
    prtcls.diag_all()
    prtcls.diag_sd_conc()
    sd_conc = np.frombuffer(prtcls.outbuf())[0]
 #   print 'sd conc pre ', sd_conc
    prtcls.step_sync(Opts,th,rv,rhod, ambient_chem=ambient_chem)
 #   print 'post step sync ', i
    prtcls.step_async(Opts)
  #  print 'post step async ', i
  
  prtcls.diag_all()
  prtcls.diag_sd_conc()
  sd_conc = np.frombuffer(prtcls.outbuf())[0]
  print('final no of SDs: ', sd_conc)
  if(sd_conc > 10 or sd_conc == 0):
    raise Exception("wrong amount of SDs were removed")
  
  prtcls.diag_all()
  prtcls.diag_wet_mom(0)
  prtcls_no = np.frombuffer(prtcls.outbuf())[0]
  print('final no of particles: ', prtcls_no)
  if(sd_conc != prtcls_no):
    raise Exception("with rcyc on, droplets were removed while some others had multiplicity > 1")
