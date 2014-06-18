import sys
sys.path.append(".")

from libcloudphxx import common

print "R_d =", common.R_d
print "c_pd =", common.c_pd
print "g =", common.g
print "p_1000 =", common.p_1000

th = 300
rv = .01

print common.th_dry2std(th, rv)	
assert common.th_std2dry(common.th_dry2std(th, rv), rv) == th

assert common.p_vs(273.16) == 611.73
