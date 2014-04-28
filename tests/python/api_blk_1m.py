import sys
sys.path.append(".")

from numpy import array as arr_t

from libcloudphxx import blk_1m

opts = blk_1m.opts_t()
opts.cond = True
opts.cevp = True
opts.revp = True
opts.conv = True
opts.accr = True
opts.sedi = False


rhod = arr_t([1.  ])
th   = arr_t([300.])
rv   = arr_t([0.  ])
rc   = arr_t([0.01])
rr   = arr_t([0.  ])
dt   = 1

blk_1m.adj_cellwise(opts, rhod, th, rv, rc, rr, dt)
print "Python: th=", th

dot_rc = arr_t([0.])
dot_rr = arr_t([0.])
blk_1m.rhs_cellwise(opts, dot_rc, dot_rr, rc, rr)
print "Python dot_rc =", dot_rc
