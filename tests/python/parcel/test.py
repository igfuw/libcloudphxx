from parcel import parcel
from rhs_blk_2m import rhs_blk_2m

from numpy import array as arr_t

# initial parameters
p_d = arr_t([100000.])
th_d = arr_t([300.])
r_v = arr_t([0.01])
w = 1.
dt = .1
nt = 100

rhs = rhs_blk_2m(dt)

parcel(p_d, th_d, r_v, w, dt, nt, rhs)

print p_d, th_d, r_v
