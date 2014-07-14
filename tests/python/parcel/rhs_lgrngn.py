from libcloudphxx import lgrngn
from numpy import frombuffer, arange, ndarray
from io import open
from os import mkdir

class rhs_lgrngn:

  # ctor
  def __init__(self, dt, sd_conc, dry_distros, chem_gas):
    opts_init = lgrngn.opts_init_t()
    opts_init.sd_conc_mean = sd_conc
    opts_init.dry_distros = dry_distros

    backend = lgrngn.backend_t.serial
    self.prtcls = lgrngn.factory(backend, opts_init)
    self.opts = lgrngn.opts_t()
    self.opts.chem_gas = chem_gas

    # TODO: what's below should not be here...
    # TODO: outfreq
    self.outdir = 'aqq'
    self.bins_dry = 1e-6 * pow(10, -3 + arange(0, 40) * .1)
    self.bins_wet = 1e-6 * pow(10, -3 + arange(0, 25) * .2)

  # t=0 stuff
  def init(self, rhod, th_d, r_v):
    self.prtcls.init(th_d, r_v, rhod)
    try:
      mkdir(self.outdir)
    except OSError:
      pass
    self.out_snd = open(self.outdir + "/sounding.txt", mode='w')
    self.out_snd.write(u"#rhod [kg/m3]\tth_d [K] (theta dry!)\tr_v [kg/kg] (mixing ratio)\tM0 [TODO]\tM1 [TODO]\tM2 [TODO]\tM3 [TODO]\tS_VI [kg/kg]\n")
    self.out_dry = open(self.outdir + "/spec_dry.txt", mode='w')
    self.out_dry.write(u"#r_d [m] (left bin edge)\tn [kg-3] (per mass of dry air)\n")
    self.out_wet = open(self.outdir + "/spec_wet.txt", mode='w')
    self.out_wet.write(u"#r_d [m] (left bin edge)\tn [kg-3] (per mass of dry air)\n")

  def step(self, rhod, th_d, r_v, dot_th, dot_rv):
    th_d_copy = th_d.copy()
    r_v_copy = r_v.copy()
    self.prtcls.step_sync(self.opts, th_d_copy, r_v_copy)
    self.prtcls.step_async(self.opts)
    dot_th += th_d_copy - th_d
    dot_rv += r_v_copy - r_v

  def diag(self, rhod, th_d, r_v, t):
    # helper for gnuplot-readable output
    def save(out, xx, yy):
      for x, y in zip(xx, yy):
	out.write(u"%g\t%g\n" % (x, y))
      out.write(u"\n\n") # gnuplot treats "\n\n" as dataset separator (plot ... index n)

    # outputting spec_dry.txt
    bins = ndarray((self.bins_dry.size - 1,))
    for i in range(0, bins.size) :
      self.prtcls.diag_dry_rng(self.bins_dry[i], self.bins_dry[i+1])
      self.prtcls.diag_dry_mom(0)
      bins[i] = frombuffer(self.prtcls.outbuf())
    save(self.out_dry, self.bins_dry[0:-1], bins)

    # outputting spec_wet.txt
    bins = ndarray((self.bins_wet.size - 1,))
    for i in range(0, bins.size) :
      self.prtcls.diag_wet_rng(self.bins_wet[i], self.bins_wet[i+1])
      self.prtcls.diag_wet_mom(0)
      bins[i] = frombuffer(self.prtcls.outbuf())
    save(self.out_wet, self.bins_wet[0:-1], bins)

    # outputting sounding.txt
    self.out_snd.write(u"%g" % (rhod))
    self.out_snd.write(u"\t%g" % (th_d))
    self.out_snd.write(u"\t%g" % (r_v))

    ## cloud water
    self.prtcls.diag_wet_rng(1e-6, 25e-6) # 0 ... 1 m
    for k in range(0,4):
      self.prtcls.diag_wet_mom(k)
      self.out_snd.write(u"\t%g" % (frombuffer(self.prtcls.outbuf())))

    ## chem stuff
    self.prtcls.diag_wet_rng(0,1) # 0 ... 1 m
    self.prtcls.diag_chem(lgrngn.chem_species_t.S_VI)
    self.out_snd.write(u"\t%g" % (frombuffer(self.prtcls.outbuf())))
   
    self.out_snd.write(u"\n")
