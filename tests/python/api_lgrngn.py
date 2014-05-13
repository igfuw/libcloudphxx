import sys
sys.path.append(".")

from libcloudphxx import lgrngn

opts = lgrngn.opts_t()

prtcls = lgrngn.factory(lgrngn.backend_t.serial, lgrngn.opts_init_t())
#prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, lgrngn.opts_init_t())
#prtcls = lgrngn.factory(lgrngn.backend_t.CUDA,   lgrngn.opts_init_t())
