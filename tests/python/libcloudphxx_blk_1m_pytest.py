import sys
sys.path.append(".")

from numpy import array as arr_t

from libcloudphxx import blk_1m

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


def density_dry(press, th):
    #TODO
    return arr_t([1.107])

def pottemp_dry(th):
    #TODO
    return th

#ta f-cja jest tylko po to, aby byly keword arg., konieczna?
def adj_cellwise(press, th, rv, rc, rr, dt):
    opts = opts_cr()
    rhod = density_dry(press, th)
    thd = pottemp_dry(th)
    blk_1m.adj_cellwise(opts, rhod, thd, rv, rc, rr, dt)
    return rv, rc
