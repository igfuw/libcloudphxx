import sys 
sys.path.insert(0, "../../../build/bindings/python/")

from libcloudphxx import common, lgrngn
from scipy.io import netcdf
import inspect, numpy as np
import pdb

def micro_init(opts, state):
  # sanity check
  if (state["RH"] > 1): raise Exception("Please supply initial T,p,r_v below supersaturation")

  # using nested function to get access to opts
  def lognormal(lnr):
    from math import exp, log, sqrt, pi
    return opts["n_tot"] * exp(
      -(lnr - log(opts["mean_r"]))**2 / 2 / log(opts["stdev"])**2
    ) / log(opts["stdev"]) / sqrt(2*pi);

  # lagrangian scheme options
  opts_init = lgrngn.opts_init_t()  
  for opt in ["dt", "sd_conc_mean"]:  
    setattr(opts_init, opt, opts[opt])
  opts_init.dry_distros = {opts["kappa"]:lognormal}
  opts_init.kernel = lgrngn.kernel_t.geometric #TODO: will not be needed soon (libcloud PR #89)

  # initialitation
  micro = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
  micro.init(state["th_d"], state["r_v"], state["rhod"])
  return micro

def micro_step(micro, state, info):
  opts = lgrngn.opts_t()
  micro.step_sync(opts, state["th_d"], state["r_v"], state["rhod"]) 

def stats(state, info):
  state["T"] = np.array([common.T(state["th_d"][0], state["rhod"][0])])
  state["RH"] = state["p"] * state["r_v"] / (state["r_v"] + common.eps) / common.p_vs(state["T"][0])
  info["RH_max"] = max(info["RH_max"], state["RH"])

def output_init(outfile, state):
  f = netcdf.netcdf_file(outfile, 'w')
  f.createDimension('t', None)
  
  units = {"z" : "m", "t" : "s", "r_v" : "kg/kg", "th_d" : "K", "rhod" : "kg/m3", "p" : "Pa", "T" : "K", "RH" : "1"}

  for var in state:
    f.createVariable(var, 'd', ('t',))
    f.variables[var].unit = units[var]
  return f

def output_save(f, state):
  rec = f.variables['t'].shape[0]
  f.variables['t'][rec] = state["t"]
  for var, val in state.iteritems():
    f.variables[var][rec] = val

def save_attrs(f, dictnr):
  for var, val in dictnr.iteritems():
    setattr(f, var, val)

def parcel(dt=.1, z_max=2000, w=1, T_0=300, p_0=101300, r_0=.022, outfile="test.nc", outfreq=100, sd_conc_mean=64, kappa=.5,
  mean_r = .04e-6 / 2, stdev  = 1.4, n_tot  = 60e6
):
  # packing function arguments into "opts" dictionary
  args, _, _, _ = inspect.getargvalues(inspect.currentframe())
  opts = dict(zip(args, [locals()[k] for k in args]))

  th_0 = T_0 * (common.p_1000 / p_0)**(common.R_d / common.c_pd)
  nt = int(z_max / (w * dt))
  state = {
    "t" : 0, "z" : 0,
    "r_v" : np.array([r_0]), "p" : p_0,
    "th_d" : np.array([common.th_std2dry(th_0, r_0)]), 
    "rhod" : np.array([common.rhod(p_0, th_0, r_0)]),
    "T" : None, "RH" : None
  }
  info = { "RH_max" : 0 }
  with output_init(outfile, state) as f:
    # t=0 : init & save
    stats(state, info)
    micro = micro_init(opts, state)
    output_save(f, state)

    # timestepping
    for it in range(1,nt+1):
      # diagnostics
      # the reasons to use analytic solution:
      # - independent of dt
      # - same as in 2D kinematic model
      state["z"] += w * dt
      state["t"] = it * dt
      state["p"] = common.p_hydro(state["z"], th_0, r_0, 0, p_0)
      state["rhod"][0] = common.rhod(state["p"], th_0, r_0)

      # microphysics
      micro_step(micro, state, info)
      stats(state, info)

      # TODO: only if user wants to stop @ RH_max
      #if (state["RH"] < info["RH_max"]): break

      # output
      if (it % outfreq == 0): output_save(f, state)

    save_attrs(f, info)
    save_attrs(f, opts)

parcel()
