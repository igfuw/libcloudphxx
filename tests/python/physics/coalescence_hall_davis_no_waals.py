# test coalescence algorithm for geometric kernel with Hall efficiencies for drops with r>30um and Davis&Rogers efficiencies for smaller ones
# by comparing mass density function with results of EFM modeling

import sys 
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn
from math import exp, sqrt
import numpy as np

#total time of simulation
simulation_time = 1800

#root mean square deviation
def RMSD(a1, a2):
  nonempty = 0 
  tot = 0 
  for i in range(a1.size):
    if(a1[i] > 0 or a2[i] > 0): 
      tot+=pow(a1[i] - a2[i], 2)
      nonempty+=1
  return np.sqrt(tot/nonempty)

#initial conditions, ca. 1g / m^3
r_zero = 30.084e-6                      # can't be greater due to rd_max_init = 1e-5
n_zero = 1.25 * pow(2,23)

# initial exponential distribution in droplet volume, cf. Shima 2009
# as a function of ln(r)
def expvolumelnr(lnr):
  r=np.exp(lnr)
  return n_zero * 3.*np.power(r,3)/np.power(r_zero,3)*np.exp(- np.power((r/r_zero),3));

opts_init = lgrngn.opts_init_t()
opts_init.dt = simulation_time
opts_init.sstp_coal = simulation_time

opts_init.dx = 100
opts_init.dz = 1
opts_init.nx = 1
opts_init.nz = 1
opts_init.x1 = 100
opts_init.z1 = 1
opts_init.sedi_switch = False

rhod =   1. * np.ones((opts_init.nx, opts_init.nz))
th   = 300. * np.ones((opts_init.nx, opts_init.nz))
rv   = 0.01 * np.ones((opts_init.nx, opts_init.nz))

kappa = 1e-10

opts_init.dry_distros = {kappa:expvolumelnr}

opts_init.sd_conc = pow(2,14)
opts_init.n_sd_max = opts_init.sd_conc
opts_init.kernel = lgrngn.kernel_t.hall_davis_no_waals

opts = lgrngn.opts_t()
opts.adve = False
opts.sedi = False
opts.cond = False
opts.coal = True
opts.chem = False

bins = 6. *  pow(10, -6 + np.arange(150)/50.)

#get mass density function
def diag(arg):
  for i in range(arg.size) :
    prtcls.diag_all()
    prtcls.diag_wet_mass_dens( (bins[i] + bins[i+1])/2. ,0.62 ) #sigma0 = 0.62 like in Shima (2009)
    arg[i]= np.frombuffer(prtcls.outbuf()).mean() * 1000 # * 1000 to get grams 

results = np.zeros(bins.size-1)

# results using Bott's method for Beard76 term vel
bott1800 = np.array([1.38199e-05, 1.57908e-05, 1.79305e-05, 2.02468e-05, 2.27597e-05, 2.55066e-05, 2.85295e-05, 3.18587e-05, 3.55041e-05, 3.94589e-05, 4.3714e-05, 4.82583e-05, 5.3107e-05, 5.82765e-05, 6.37843e-05, 6.96515e-05, 7.59039e-05, 8.25861e-05, 8.97817e-05, 9.7613e-05, 0.000106259, 0.000115913, 0.000126733, 0.000138828, 0.000152281, 0.000167097, 0.000183154, 0.000200313, 0.000218507, 0.000237468, 0.000256677, 0.000274845, 0.000290957, 0.000304515, 0.000314301, 0.000317331, 0.000311517, 0.000297198, 0.000275879, 0.000250005, 0.000222306, 0.000201047, 0.000190701, 0.000181325, 0.000172454, 0.000164238, 0.000159034, 0.000156689, 0.000155278, 0.000154838, 0.000154788, 0.000155249, 0.000156931, 0.000160285, 0.000166335, 0.00017442, 0.000184645, 0.00019718, 0.000211995, 0.000229172, 0.000248871, 0.000272454, 0.000300346, 0.000331661, 0.000366821, 0.000406322, 0.000450999, 0.00050156, 0.000558714, 0.00062337, 0.000696966, 0.00078046, 0.0008754, 0.000983683, 0.00110756, 0.00124968, 0.00141387, 0.00160394, 0.00182415, 0.00208052, 0.00238027, 0.00273237, 0.00314791, 0.00364083, 0.00422857, 0.00493568, 0.00578798, 0.00682106, 0.00808022, 0.00962269, 0.0115211, 0.0138674, 0.0167778, 0.0203984, 0.0249108, 0.030539, 0.0375565, 0.0461041, 0.0563682, 0.0687794, 0.0835501, 0.100858, 0.12079, 0.137515, 0.156138, 0.183594, 0.214335, 0.248523, 0.286324, 0.327896, 0.373386, 0.422884, 0.476408, 0.533889, 0.595081, 0.659549, 0.726619, 0.795314, 0.864311, 0.931858, 0.995727, 1.05315, 1.10078, 1.13471, 1.15053, 1.14362, 1.10956, 1.04496, 0.948523, 0.822462, 0.673665, 0.514086, 0.359412, 0.225616, 0.124174, 0.058276, 0.0225813, 0.00696037, 0.00163505, 0.000278816, 3.26663e-05, 2.53925e-06, 1.27047e-07, 3.82823e-09, 6.48661e-11, 5.83051e-13, 2.62537e-15, 5.62474e-18, 0])


for vt_eq in [lgrngn.vt_t.beard76, lgrngn.vt_t.beard77, lgrngn.vt_t.beard77fast]:
  opts_init.terminal_velocity = vt_eq

  try:
    prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
  except:
    prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
  
  prtcls.init(th, rv, rhod)
  
  #simulation loop
  prtcls.step_sync(opts, th, rv, rhod)
  prtcls.step_async(opts)
      
  diag(results)
  
  rmsd = RMSD(results,bott1800)
  
  print 'RMSD = ' + str(rmsd);
  
  if(rmsd > 6e-2):
    raise Exception("Simulation result does not agree with analytic prediction")
