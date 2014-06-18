from libcloudphxx import lgrngn

class rhs_lgrngn:

  def __init__(self, dt, sd_conc, dry_distros):
    opts_init = lgrngn.opts_init_t()
    opts_init.dry_distros = dry_distros

    backend = lgrngn.backend_t.serial
    self.prtcls = lgrngn.factory(backend, opts_init)
    self.opts = lgrngn.opts_t()

  def init(self, th_d, r_v, rhod):
    self.prtcls.init(th_d, r_v, rhod)

  def __call__(self, rhod, th_d, r_v, dot_th, dot_rv):
    th_d_copy = th_d.copy()
    r_v_copy = r_v.copy()
    self.prtcls.step_sync(self.opts, th_d_copy, r_v_copy)
    self.prtcls.step_async(self.opts)
    dot_th += th_d_copy - th_d
    dot_rv += r_v_copy - r_v
