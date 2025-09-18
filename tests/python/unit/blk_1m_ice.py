import sys
sys.path.insert(0, "../../bindings/python/")

from numpy import array as arr_t # ndarray dtype default to float64, while array's is int64!
import numpy as np

from libcloudphxx import blk_1m

opts = blk_1m.opts_t()
opts.adj_nwtrph = True
opts.const_p = True
opts.th_dry = False
ice_t = blk_1m.ice_t()

rhod = arr_t([1.  ])
p    = arr_t([1.e5])
dt   = 0.1
dz   = 1
th   = arr_t([263.])
rv   = arr_t([0.1 ])
rc   = arr_t([0.01])
rr   = arr_t([0.01])
ria  = arr_t([0.])
rib  = arr_t([0.])

for _ in range(100):
    dot_th = arr_t([0.])
    dot_rc = arr_t([0.])
    dot_rr = arr_t([0.])
    dot_rv = arr_t([0.])
    dot_ria = arr_t([0.])
    dot_rib = arr_t([0.])
    blk_1m.rhs_cellwise_ice(opts, dot_th, dot_rv, dot_rc, dot_rr, dot_ria, dot_rib, rhod, p, th, rv, rc, rr, ria, rib, dt)
    flux_rain = blk_1m.rhs_columnwise(opts, dot_rr, rhod, rr, dz)
    flux_iceA = blk_1m.rhs_columnwise_ice(opts, dot_ria, rhod, ria, dz, ice_t.iceA)
    flux_iceB = blk_1m.rhs_columnwise_ice(opts, dot_rib, rhod, rib, dz, ice_t.iceB)
    th += dt * dot_th
    rv += dt * dot_rv
    rc += dt * dot_rc
    rr += dt * dot_rr
    ria += dt * dot_ria
    rib += dt * dot_rib

print('dot_rv='+str(dot_rv))
print('dot_rc='+str(dot_rc))
print('dot_rr='+str(dot_rr))
print('dot_ria='+str(dot_ria))
print('dot_rib='+str(dot_rib))
print('th='+str(th))
print('rv='+str(rv))
print('rc='+str(rc))
print('rr='+str(rr))
print('ria='+str(ria))
print('rib='+str(rib))
print('rain_flux='+str(flux_rain))
print('iceA_flux='+str(flux_iceA))
print('iceB_flux='+str(flux_iceB))

assert np.isnan(rv) == False
assert np.isnan(rc) == False
assert np.isnan(rr) == False
assert np.isnan(ria) == False
assert np.isnan(rib) == False

assert rv >= 0
assert rc >= 0
assert rr >= 0
assert ria >= 0
assert rib >= 0

assert np.isclose(ria, 2.7564e-05, rtol=1e-5, atol=1e-8)
assert np.isclose(rib, 3.2808e-06, rtol=1e-5, atol=1e-8)
