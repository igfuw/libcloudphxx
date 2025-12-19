import sys

sys.path.insert(0, "../../bindings/python/")
sys.path.insert(0, "../../../build/bindings/python/")

from libcloudphxx import lgrngn, common
import numpy as np
from math import exp, log, sqrt, pi

n_tot  = 60e6
def lognormal(lnr):
    mean_r = .04e-6 / 2
    stdev  = 1.4
    return n_tot * exp(
        -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
    ) / log(stdev) / sqrt(2*pi);

opts_init = lgrngn.opts_init_t()
backend = lgrngn.backend_t.serial
opts = lgrngn.opts_t()

kappa = .61
rd_insol = 0.5e-6
opts_init.dry_distros = {(kappa, rd_insol): lognormal}
opts_init.dt = 0.1
opts_init.sd_conc = 100
opts_init.n_sd_max = 1000
opts_init.RH_max = 0.95

opts_init.ice_switch = True
opts_init.coal_switch = False
opts_init.sedi_switch = False
opts.ice_nucl = True
opts.cond = True
opts.coal = False
opts.adve = False
opts.sedi = False


for time_dep_switch in [True, False]:
    print("time dependent ice nucleation = ", time_dep_switch)
    p = 80000.
    T = 243.
    RH = 1.
    rv = np.array([RH * common.r_vs(T, p)])
    th = np.array([T / common.exner(p)])
    rhod = np.array([common.rhod(p, th[0], rv[0])])

    opts_init.time_dep_ice_nucl = time_dep_switch
    prtcls = lgrngn.factory(backend, opts_init)
    prtcls.init(th, rv, rhod)
    for _ in range(500):
        prtcls.step_sync(opts, th, rv, rhod)
        prtcls.step_async(opts)
    prtcls.diag_all()
    prtcls.diag_ice_mix_ratio()
    ri = np.frombuffer(prtcls.outbuf())
    print("ice mixing ratio ", ri[0] * 1e3, "g/kg")
    print("water vapor mixing ratio ", rv[0] * 1e3, "g/kg")
    assert np.isnan(ri[0]) == False
    assert np.isnan(rv[0]) == False
    assert rv[0] >= 0
    assert ri[0] >= 0