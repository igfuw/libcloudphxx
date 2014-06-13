import sys
sys.path.append(".")

from numpy import array as arr_t
from constants_pytest import Rd, Rv, cp, p0
from libcloudphxx import blk_1m, common

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


def density_dry(rv, press, T):
    rho_d = press / T / Rd / (1. + rv * Rv/Rd)
    return arr_t([rho_d])

def pottemp_dry(rv, press, T):
    theta_d = T * (p0 / press * (1. + rv * Rv/Rd))**(Rd/cp) 
    return arr_t([theta_d]) 

#ta f-cja jest tylko po to, aby byly keword arg., konieczna?
def adj_cellwise(press, T, rv, rc, rr, dt):
    opts = opts_cr()
    rho_d = density_dry(rv, press, T)
    th_d = pottemp_dry(rv, press, T)
    blk_1m.adj_cellwise(opts, rho_d, th_d, rv, rc, rr, dt)
    return rv, rc

print common.p_vs(283.15)
