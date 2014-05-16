import sys
sys.path.append(".")

from numpy import array as arr_t

from libcloudphxx import blk_2m

opts = blk_2m.opts_t()

opts.acti = True
opts.cond = True
opts.acnv = True
opts.accr = True
opts.sedi = False


rhod = arr_t([1.  ])
th   = arr_t([300.])
rv   = arr_t([0.  ])
rc   = arr_t([0.01])
nc   = arr_t([1e-3])
rr   = arr_t([0.  ])
nr   = arr_t([0.  ])
dt   = 1

dot_th = arr_t([0.])
dot_rv = arr_t([0.])
dot_rc = arr_t([0.])
dot_nc = arr_t([0.])
dot_rr = arr_t([0.])
dot_nr = arr_t([0.])

th_old = th.copy()
blk_2m.rhs_cellwise(
  opts, 
  dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr, 
  rhod, 
  th, rv, rc, nc, rr, nr, 
  dt
)
assert th == th_old # th should be left untouched
assert dot_th != 0 # but it's derivative should be non-zero

dz = 1
dot_rr_old = dot_rr.copy()
dot_nr_old = dot_nr.copy()
flux = blk_2m.rhs_columnwise(opts, dot_rr, dot_nr, rhod, rr, nr, dt, dz) 
#TODO assert flux == 0
assert dot_rr == dot_rr_old and dot_nr == dot_nr_old # no rain water -> no precip

