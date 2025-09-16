import sys
sys.path.insert(0, "../../bindings/python/")

from numpy import array as arr_t # ndarray dtype default to float64, while array's is int64!

from libcloudphxx import blk_1m

opts = blk_1m.opts_t()
ice_t = blk_1m.ice_t()
print("cond =", opts.cond)
print("cevp =", opts.cevp)
print("revp =", opts.revp)
print("conv =", opts.conv)
print("accr =", opts.accr)
print("sedi =", opts.sedi)
print("homA1 =", opts.homA1)
print("homA2 =", opts.homA2)
print("hetA =", opts.hetA)
print("hetB =", opts.hetB)
print("depA =", opts.depA)
print("depB =", opts.depB)
print("rimA =", opts.rimA)
print("rimB =", opts.rimB)
print("melA =", opts.melA)
print("melB =", opts.melB)
print("r_c0 =", opts.r_c0)
print("r_eps =", opts.r_eps)
print("adj_nwtrph =", opts.adj_nwtrph)
print("th_dry =", opts.th_dry)
print("const_p =", opts.const_p)

test_cases = { \
  "RK4, variable pressure and dry theta" : [0, 1, 0], # adj_nwtrph, th_dry, const_p 
  "RK4, constant pressure and 'standard' theta" : [0, 0, 1],
  "Newton-Raphson, variable pressure and dry theta" : [1, 1, 0],
  "Newton-Raphson, constant pressure and 'standard' theta" : [1, 0, 1],
}

# test saturation adjustment
rhod = arr_t([1.  ])
p    = arr_t([1.e5])
th   = arr_t([300.])
rc   = arr_t([0.01])
rv   = arr_t([0.  ])
rr   = arr_t([0.  ])
dt   = 1
dz   = 1

def test_sat_adj(opts, name):
  print("Testing saturation adjustment with " + name)
  th_new = th.copy()
  rv_new = rv.copy()
  rc_new = rc.copy()
  rr_new = rr.copy()
  blk_1m.adj_cellwise(opts, rhod, p, th_new, rv_new, rc_new, rr_new, dt)
  assert th != th_new # some water should have evaporated
  assert rv != rv_new
  assert rc != rc_new
  assert rr == rr_new

for name, opt in test_cases.items():
  opts.adj_nwtrph = opt[0]
  opts.th_dry = opt[1]
  opts.const_p = opt[2]
  test_sat_adj(opts, name)

# test RHS cellwise
dot_th = arr_t([0.])
dot_rv = arr_t([0.])
dot_rc = arr_t([0.])
dot_rr = arr_t([0.])

def test_rhs_cell(opts, name):
  print("Testing RHS cellwise with " + name)
  dot_th_new = dot_th.copy()
  dot_rv_new = dot_rv.copy()
  dot_rc_new = dot_rc.copy()
  dot_rr_new = dot_rr.copy()
  if opts.adj_nwtrph:
    rr = arr_t([0.01])
    blk_1m.rhs_cellwise_revap(opts, dot_th_new, dot_rv_new, dot_rc_new, dot_rr_new, rhod, p, th, rv, rc, rr, dt)
    assert dot_th_new != 0 # some rain should have evaporated
    assert dot_rv_new != 0
  else:
    rr = arr_t([0.00])
    blk_1m.rhs_cellwise(opts, dot_rc_new, dot_rr_new, rc, rr)
  assert dot_rc_new != 0 # some water should have coalesced
  assert dot_rr_new != 0

for name, opt in test_cases.items():
  opts.adj_nwtrph = opt[0]
  opts.th_dry = opt[1]
  opts.const_p = opt[2]
  test_rhs_cell(opts, name)

# test RHS columnwise
print("Testing RHS columnwise")
rr   = arr_t([0.])
dot_rr_old = dot_rr.copy()
flux = blk_1m.rhs_columnwise(opts, dot_rr, rhod, rr, dz)
assert flux == 0
assert dot_rr == dot_rr_old # no rain water -> no precip

#test ice, TODO: separate file/test?
th   = arr_t([230.])  #testing ice physics
ria = arr_t([0.1])
rib = arr_t([0.1])
dot_rc = arr_t([0.])
dot_rr = arr_t([0.])
dot_rv = arr_t([0.])
dot_ria = arr_t([0.])
dot_rib = arr_t([0.])

def test_rhs_cell_ice(opts, name):
  print("Testing RHS cellwise ice with " + name)
  dot_th_new = dot_th.copy()
  dot_rv_new = dot_rv.copy()
  dot_rc_new = dot_rc.copy()
  dot_rr_new = dot_rr.copy()
  dot_ria_new = dot_ria.copy()
  dot_rib_new = dot_rib.copy()
  blk_1m.rhs_cellwise_ice(opts, dot_th_new, dot_rv_new, dot_rc_new, dot_rr_new, dot_ria_new, dot_rib_new, rhod, p, th, rv, rc, rr, ria, rib, dt)
  assert dot_ria != 0
  assert dot_rib != 0

for name, opt in test_cases.items():
  opts.adj_nwtrph = opt[0]
  opts.th_dry = opt[1]
  opts.const_p = opt[2]
  test_rhs_cell_ice(opts, name)

#testing sedimentation of ice
flux_iceA = blk_1m.rhs_columnwise_ice(opts, dot_ria, rhod, ria, dz, ice_t.iceA)
flux_iceB = blk_1m.rhs_columnwise_ice(opts, dot_rib, rhod, rib, dz, ice_t.iceB)
assert flux_iceA != 0
assert flux_iceB != 0
