import sys
sys.path.append(".")

from libcloudphxx import lgrngn

def lognormal(x):
  return 0

opts_init = lgrngn.opts_init_t()
kappa = 0
spctr = lognormal
opts_init.dry_distros = {kappa:spctr};
opts_init.nx = 0
opts_init.ny = 0
opts_init.nz = 0
opts_init.dx = 1
opts_init.dy = 1
opts_init.dz = 1
opts_init.kernel = lgrngn.kernel_t.geometric
opts_init.sd_conc_mean = 64

# intentionally overwritten - only serial is olways available
backend = lgrngn.backend_t.OpenMP
backend = lgrngn.backend_t.CUDA
backend = lgrngn.backend_t.serial

prtcls = lgrngn.factory(backend, opts_init)

opts = lgrngn.opts_t()
opts.adve = True
opts.sedi = True
opts.cond = True
opts.coal = True
opts.RH_max = 44
opts.sstp_cond = 1
opts.sstp_coal = 1

#rhod = 
#th = 
#rv = 

#prtcls.init(th, rv, rhod)
#prtcls.step_sync(opts, th, rv)
#prtcls.step_async(opts)
prtcls.diag_dry_rng(0.,1.)
prtcls.diag_wet_rng(0.,1.)
#prtcls.diag_dry_mom(1)
#prtcls.diag_wet_mom(1)
#prtcls.diag_sd_conc()
#print prtcls.outbuf()
