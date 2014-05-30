from libcloudphxx import lgrngn

class rhs_lgrngn:

  def __init__(self, sd_conc):
    opts_init = opts_init = lgrngn.opts_init_t()
    backend = lgrngn.backend_t.serial
    self.prtcls = lgrngn.factory(backend, opts_init)

  def __call__(self, rhod, th_d, r_v, dot_th, dot_rv):
    dot_th += 0
    dot_rv += 0
