from parcel import parcel
from rhs_blk_2m import rhs_blk_2m
from libcloudphxx.common import th_std2dry, th_dry2std

from numpy import array as arr_t

# initial parameters
p_d = arr_t([100000.]) # TODO: p
r_v = arr_t([0.01])
th_d = arr_t([th_std2dry(300., r_v[0])])
w = 1.
dt = .1
nt = 100

rhs = rhs_blk_2m(dt)

parcel(p_d, th_d, r_v, w, dt, nt, rhs)

print p_d, arr_t([th_dry2std(th_d[0], r_v[0])]), r_v
