from libcloudphxx import lgrngn
from numpy import frombuffer, arange, ndarray

class rhs_lgrngn:

  def __init__(self, dt, sd_conc, dry_distros):
    opts_init = lgrngn.opts_init_t()
    opts_init.dry_distros = dry_distros

    backend = lgrngn.backend_t.serial
    self.prtcls = lgrngn.factory(backend, opts_init)
    self.opts = lgrngn.opts_t()

    # TODO: what's below should not be here...
    self.bins_dry = arange(0, 40)
    self.bins_dry = 1e-6 * pow(10, -3 + self.bins_dry * .1)

    self.bins_wet = arange(0, 25)
    self.bins_wet = 1e-6 * pow(10, -3 + self.bins_wet * .2)

  def init(self, th_d, r_v, rhod):
    self.prtcls.init(th_d, r_v, rhod)

  def __call__(self, rhod, th_d, r_v, dot_th, dot_rv):
    # diagnostics (TODO: of course should not be here!)
    self.prtcls.diag_wet_rng(0,1) # 0 ... 1 m
    self.prtcls.diag_dry_mom(3)
    print "3rd moment  =", frombuffer(self.prtcls.outbuf()) 
    self.prtcls.diag_chem(lgrngn.chem_aq.S_VI)
    print "sulfur / kg =", frombuffer(self.prtcls.outbuf())
    print "sulfur / m3 =", frombuffer(self.prtcls.outbuf()) * rhod
    
    bins = ndarray((self.bins_dry.size - 1,))
    for i in arange(0, self.bins_dry.size - 1) :
      self.prtcls.diag_dry_rng(self.bins_dry[i], self.bins_dry[i+1])
      self.prtcls.diag_dry_mom(0)
      bins[i] = frombuffer(self.prtcls.outbuf())

    #import subprocess
    #gnuplot = subprocess.Popen(["/usr/bin/gnuplot"], stdin=subprocess.PIPE)
    #gnuplot.stdin.write("set term dumb 79 25\n")
    #gnuplot.stdin.write("set logscale xy\n")
    #gnuplot.stdin.write("plot '-' with histeps \n")
    #for i,j in zip(self.bins_dry[0:-1], bins):
    #  gnuplot.stdin.write("%f %f\n" % (i,j))
    #gnuplot.stdin.write("e\n")
    #gnuplot.stdin.flush()

    # the timestep
    th_d_copy = th_d.copy()
    r_v_copy = r_v.copy()
    self.prtcls.step_sync(self.opts, th_d_copy, r_v_copy)
    self.prtcls.step_async(self.opts)
    dot_th += th_d_copy - th_d
    dot_rv += r_v_copy - r_v

