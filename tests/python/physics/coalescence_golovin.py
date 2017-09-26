#roughly test coalescence algorithm by comparing mass density function with analytic prediction of Golovin

import sys 
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn
from math import exp, log, sqrt, pi
from scipy import special
import numpy as np

#total time of simulation
simulation_time = 800

def spherevol(r):
  return 4./3.*pow(r,3)*np.pi;

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
r_zero = 30.084e-6
n_zero = pow(2,23)
v_zero = spherevol(r_zero)

#Golovin kernel parameter
b=1500.

# initial exponential distribution in droplet volume, cf. Shima 2009
# as a function of ln(r)
def expvolumelnr(lnr):
  r=np.exp(lnr)
  return n_zero * 3.*np.power(r,3)/np.power(r_zero,3)*np.exp(- np.power((r/r_zero),3));

# analytic Golovin result for expvolume initial conditions
# returns number density as a function of volume
# cf. Scott et al 1967, eq. 2.7
def golovin(v,t,n0,v0,b): 
  x = v/v0
  T = b*n0*v0*t
  tau = 1-np.exp(-T)
  bessel = special.iv(1,2*x*np.sqrt(tau))
  result = 0.
  if(not np.isinf(bessel)): 
    result =   n0/v0 *  bessel * (1-tau) * np.exp(-x*(tau+1)) / x / np.sqrt(tau)
  if (np.isnan(result)):
    result = 0.
  return result 

opts_init = lgrngn.opts_init_t()
opts_init.dt = simulation_time
opts_init.sstp_coal = simulation_time

rhod = 1. * np.ones((1,))
th = 300. * np.ones((1,))
rv = 0.01 * np.ones((1,))

kappa = 1e-10

opts_init.dry_distros = {kappa:expvolumelnr}

opts_init.kernel = lgrngn.kernel_t.golovin
opts_init.terminal_velocity = lgrngn.vt_t.beard77
opts_init.kernel_parameters = np.array([b]);


opts = lgrngn.opts_t()
opts.adve = False
opts.sedi = False
opts.cond = False
opts.coal = True
opts.chem = False

bins = pow(10, -6 + np.arange(150)/50.)

#get mass density function
def diag(arg):
  for i in range(arg.size) :
    prtcls.diag_all()
    prtcls.diag_wet_mass_dens( (bins[i] + bins[i+1])/2. ,0.62 ) #sigma0 = 0.62 like in Shima (2009)
    arg[i]= np.frombuffer(prtcls.outbuf()).mean()

#get total numer of particles
def partno():
  prtcls.diag_all()
  prtcls.diag_wet_mom(0)
  return np.frombuffer(prtcls.outbuf())[0]

#get Golovin mass density prediction
def calc_golovin(res,t,n0,v0,b):
  for i in range(res.size) :
    vol = spherevol((bins[i]+bins[i+1])/2.)
    res[i] = golovin(vol,t,n0,v0,b)
    res[i] *= vol * vol * 3000.  #mass density function = 3 * dens_water * volume^2 * number_density_function(volume)

results = np.zeros(bins.size-1)
golovin_results = np.zeros(bins.size-1)

#loop to test sd_conc and const_multi options
for i in range(0,2):
  if(i==0):
    opts_init.sd_conc = pow(2,14)
    opts_init.n_sd_max = pow(2,14)
  else:
    opts_init.sd_conc = 0
    opts_init.sd_const_multi = 1000
    opts_init.n_sd_max = int(float(n_zero) / opts_init.sd_const_multi + 10)

  try:
    prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
  except:
    prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
  
  prtcls.init(th, rv, rhod)
  init_number_of_particles = partno()
  
  #simulation loop
  prtcls.step_sync(opts, th, rv, rhod)
  prtcls.step_async(opts)
      
  diag(results)
  calc_golovin(golovin_results,simulation_time,init_number_of_particles,v_zero,b)
  rmsd = RMSD(results,golovin_results)
  
  if(i==0):
    print 'sd_conc RMSD = ' + str(rmsd);
    limit = 1e-5;
  else:
    print 'const_multi RMSD = ' + str(rmsd);
    limit = 2e-5; # constant multiplicity doesn't represent tails of the distribution so well and mass densty function depends on large tail?

  if(rmsd > limit):
    raise Exception("Simulation result does not agree with analytic prediction")
