#from libcloupdhxx.common.earth import g
#from libcloudphxx.common.moist_air import R_d, c_pd
#from libcloudphxx.common.theta_std import p_1000

R_d = p_0 = c_pd = g = 44

from numpy import array as arr_t

# p_d, th_d, r_v should contain initial values
#                and       
def parcel(p_d, th_d, r_v, w, dt, nt, rhs):

  # perfect gas for for dry air
  def rhod(p_d, th_d):
    def T(p_d, th_d):
      return th_d * pow(p_d / p_0, R_d / c_pd)
    return p_d / R_d / T(p_d, th_d)

  # Euler-like integration
  for _ in range(nt):
    # first, adjusting thr pressure using hydrostatic law
    p_d += dt * (-g * rhod(p_d, th_d) * w)
    # computing rhs for th and rv
    dot_th = arr_t([0])
    dot_rv = arr_t([0])
    rhs(rhod(p_d, th_d), th_d, r_v, dot_th, dot_rv)
    # applying the rhs
    th_d += dt * dot_th
    r_v  += dt * dot_rv
