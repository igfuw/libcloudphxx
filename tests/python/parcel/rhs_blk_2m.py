from libcloudphxx import blk_2m

from numpy import array as arr_t

class rhs_blk_2m:

  def init(self, th_d, r_v, rhod):
    pass

  def __init__(self, dt):
    self.dt = dt
    self.opts = blk_2m.opts_t()

# TODO!
#    self.opts.dry_distros = [ 
#      {"mean_rd":.04e-6 / 2, "sdev_rd":1.4, "N_stp":60e6, "chem_b":.55},
#      {"mean_rd":.15e-6 / 2, "sdev_rd":1.6, "N_stp":40e6, "chem_b":.55}
#    ]

    self.rc = arr_t([0.])
    self.nc = arr_t([0.])
    self.rr = arr_t([0.])
    self.nr = arr_t([0.])

  def __call__(self, rhod, th, rv, dot_th, dot_rv):
    dot_rc = arr_t([0.])
    dot_nc = arr_t([0.])
    dot_rr = arr_t([0.])
    dot_nr = arr_t([0.])
    blk_2m.rhs_cellwise(self.opts, 
      dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
      rhod,
      th, rv, 
      self.rc, self.nc, self.rr, self.nr,
      self.dt
    );
    self.rc += dot_rc * self.dt
    self.nc += dot_nc * self.dt
    self.rr += dot_rr * self.dt
    self.nr += dot_nr * self.dt
