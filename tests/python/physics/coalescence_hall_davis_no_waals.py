#roughly test coalescence algorithm by comparing mass density function with analytic prediction of Golovin

import sys 
sys.path.insert(0, "../../bindings/python/")
sys.path.insert(0,'/home/pracownicy/pdziekan/praca/code/libcloudphxx/build/bindings/python/')

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
r_zero = 3.0531e-6                      # can't be greater due to rd_max_init = 1e-5
n_zero = 1.25 * pow(2,23)

# initial exponential distribution in droplet volume, cf. Shima 2009
# as a function of ln(r)
def expvolumelnr(lnr):
  r=np.exp(lnr)
  return n_zero * 3.*np.power(r,3)/np.power(r_zero,3)*np.exp(- np.power((r/r_zero),3));

opts_init = lgrngn.opts_init_t()
opts_init.dt = 1
opts_init.sstp_coal = 1

rhod = 1. * np.ones((1,))
th = 300. * np.ones((1,))
rv = 0.01 * np.ones((1,))

kappa = 50. #unrealistic, but we want initial wet radii of the order of 30 um so that coalescence takes place 

opts_init.dry_distros = {kappa:expvolumelnr}

opts_init.sd_conc_mean = pow(2,14)

opts_init.kernel = lgrngn.kernel_t.hall_davis_no_waals

try:
  prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
except:
  prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)

prtcls.init(th, rv, rhod)

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
    arg[i]= np.frombuffer(prtcls.outbuf()).mean()

results = np.zeros(bins.size-1)

#simulation loop
for t in range(int((simulation_time)/opts_init.dt)):
  prtcls.step_sync(opts, th, rv, rhod)
  prtcls.step_async(opts)
    
diag(results)
rmsd = RMSD(results,golovin_results)

print 'RMSD = ' + str(rmsd);

if(rmsd > 5.7e-6):
  raise Exception("Simulation result does not agree with analytic prediction")
