import sys
sys.path.append(".")

from libcloudphxx import lgrngn

from numpy import array as arr_t, frombuffer

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
print opts_init.nx
print opts_init.ny
print opts_init.nz
print opts_init.dx
print opts_init.dy
print opts_init.dz
print opts_init.kernel 
assert opts_init.kernel == lgrngn.kernel_t.geometric
print opts_init.sd_conc_mean

# intentionally overwritten - only serial is always available
backend = lgrngn.backend_t.OpenMP
backend = lgrngn.backend_t.CUDA
backend = lgrngn.backend_t.serial

prtcls = lgrngn.factory(backend, opts_init)

opts = lgrngn.opts_t()
print opts.adve
print opts.sedi
print opts.cond
print opts.coal
#print opts.chem # TODO
print opts.RH_max
print opts.sstp_cond
print opts.sstp_coal
#print opts.sstp_chem TODO

opts.cond = False #TODO: it doesn't work yet

rhod = arr_t([  1.])
th   = arr_t([300.])
rv   = arr_t([  0.01])

prtcls.init(th, rv, rhod)
prtcls.step_sync(opts, th, rv)
prtcls.step_async(opts)
prtcls.diag_dry_rng(0.,1.)
prtcls.diag_wet_rng(0.,1.)
prtcls.diag_dry_mom(1)
prtcls.diag_wet_mom(1)
#prtcls.diag_chem(constituent)? TODO
prtcls.diag_sd_conc()

print frombuffer(prtcls.outbuf())
