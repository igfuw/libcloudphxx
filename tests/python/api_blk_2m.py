import sys
sys.path.append("../../bindings/python/")

from numpy import array as arr_t

from libcloudphxx import blk_2m

opts = blk_2m.opts_t()

print "acti =", opts.acti 
print "cond =", opts.cond 
print "acnv =", opts.acnv 
print "accr =", opts.accr
print "sedi =", opts.sedi
print "RH_max =", opts.RH_max

opts.dry_distros = [
  {"mean_rd":.04e-6 / 2, "sdev_rd":1.4, "N_stp":60e6, "chem_b":.55},
  {"mean_rd":.15e-6 / 2, "sdev_rd":1.6, "N_stp":40e6, "chem_b":.55}
]

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
assert flux == 0
assert dot_rr == dot_rr_old and dot_nr == dot_nr_old # no rain water -> no precip

