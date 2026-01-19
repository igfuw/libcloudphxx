from parcel import parcel
from rhs_blk_2m import rhs_blk_2m
from rhs_lgrngn import rhs_lgrngn
from libcloudphxx.common import th_std2dry, th_dry2std
from libcloudphxx.common import p_vs
from libcloudphxx.common import eps, p_1000, R_d, c_pd

from libcloudphxx.lgrngn import chem_species_t

from numpy import array as arr_t
from math import exp, log, sqrt, pi


# initial parameters
T    = arr_t([282.2])
p    = arr_t([95000.]) 
p_v  = arr_t([0.95 * p_vs(T[0])])
p_d  = p - p_v 
r_v  = eps * p_v / p_d
#th_d = arr_t([th_std2dry(300., r_v[0])])
th_d = T * pow(p_1000 / p_d[0], R_d / c_pd) 
w    = 0.5
dt   = .1
nt   = int(600 / w / dt) # 600 metres

# blk_2m-specific parameter
# TODO: spectrum

# lgrngn-specific parameters
sd_conc = 44
def lognormal(lnr):
  mean_r = .08e-6 / 2
  stdev = 2
  n_tot = 566e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);
kappa = .61
rd_insol = 0.
# rho = 1.8 # TODO

chem_gas = {
  chem_species_t.SO2  : 200e-12,
  chem_species_t.O3   : 50e-9,
  chem_species_t.H2O2 : 500e-12
}

# running all three
for rhs in [
  rhs_blk_2m(dt), 
  rhs_lgrngn(dt, sd_conc, {(kappa, rd_insol):lognormal}, chem_gas)
]:
  parcel(p_d, th_d, r_v, w, dt, nt, rhs)
  print(p_d, arr_t([th_dry2std(th_d[0], r_v[0])]), r_v)

