import sys
sys.path.append(".")

from libcloudphxx import lgrngn

opts = lgrngn.opts_t()

opts_init = lgrngn.opts_init_t()

backend = lgrngn.backend_t.serial
#backend = lgrngn.backend_t.OpenMP
#backend = lgrngn.backend_t.CUDA

prtcls = lgrngn.factory(backend, opts_init)

#prtcls.step_async(opts)
#print prtcls.outbuf()
