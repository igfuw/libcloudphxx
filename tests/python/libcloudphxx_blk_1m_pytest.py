import sys
sys.path.append(".")

from numpy import array as arr_t
from constants_pytest import Rd, Rv, cp, p0
from libcloudphxx import blk_1m, common
import analytic_blk_1m_pytest as anpy

def opts_cr(cond = True, cevp = True, revp = True, conv = True,
            accr = True, sedi = False):
    opts = blk_1m.opts_t()
    opts.cond = cond
    opts.cevp = cevp
    opts.revp = revp
    opts.conv = conv
    opts.accr = accr
    opts.sedi = sedi
    return opts

#ta f-cja jest tylko po to, aby byly keword arg., konieczna?
def adj_cellwise(press, T, rv, rc, rr, dt):
    opts = opts_cr()
    rho_d = arr_t(anpy.density_dry(rv, press, T))
    th_d = arr_t(anpy.pottemp_dry(rv, press, T))
    blk_1m.adj_cellwise(opts, rho_d, th_d, rv, rc, rr, dt)
    return rv, rc

