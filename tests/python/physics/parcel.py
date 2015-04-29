import sys 
sys.path.insert(0, "../../../build/bindings/python/")

from libcloudphxx import common
from scipy.io import netcdf

def micro(state):
  state["th_d"] += 0
  state["r_v"] += 0

def save_init(outfile, state):
  f = netcdf.netcdf_file(outfile, 'w')
  f.createDimension('t', None)
  
  units = {"z" : "m", "t" : "s", "r_v" : "kg/kg", "th_d" : "K", "rhod" : "kg/m3", "p" : "Pa"}

  for var in state:
    f.createVariable(var, 'd', ('t',))
    f.variables[var].unit = units[var]
  save_record(f, state)
  return f

def save_record(f, state):
  rec = f.variables['t'].shape[0]
  f.variables['t'][rec] = state["t"]
  for var, val in state.iteritems():
    f.variables[var][rec] = val

def save_attrs(f, dictnr):
  for var, val in dictnr.iteritems():
    setattr(f, var, val)

def parcel(dt=.1, z_max=100, w=1, T_0=300, p_0=101300, r_0=.1, outfile="test.nc", outfreq=100):
  th_0 = T_0 * (common.p_1000 / p_0)**(common.R_d / common.c_pd)
  nt = int(z_max / (w * dt))
  state = {
    "t" : 0, "z" : 0,
    "r_v" : r_0, "p" : p_0,
    "th_d" : common.th_std2dry(th_0, r_0), 
    "rhod" : common.rhod(p_0, th_0, r_0)
  }
  with save_init(outfile, state) as f:
    for it in range(1,nt+1):
      # diagnostics
      # the reasons to use analytic solution:
      # - independent of dt
      # - same as in 2D kinematic model
      state["z"] += w * dt
      state["t"] = it * dt
      state["p"] = common.p_hydro(state["z"], th_0, r_0, 0, p_0)
      state["rhod"] = common.rhod(state["p"], th_0, r_0)

      # microphysics
      micro(state)

      # output
      if (it % outfreq == 0): save_record(f, state)

      # stats

    #save_stats(f, )
    import inspect
    args, _, _, _ = inspect.getargvalues(inspect.currentframe())
    save_attrs(f, dict(zip(args, [locals()[k] for k in args])))

parcel()
