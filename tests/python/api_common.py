import sys
sys.path.append("../../bindings/python/")

#<listing-1>
from libcloudphxx import common
print "common.p_vs(273.16)=", common.p_vs(273.16)
assert abs(common.p_vs(273.16) - 611.73) < .001
#</listing-1>

print "R_d =", common.R_d
print "c_pd =", common.c_pd
print "g =", common.g
print "p_1000 =", common.p_1000

th = 300
rv = .01

print common.th_dry2std(th, rv)	
assert common.th_std2dry(common.th_dry2std(th, rv), rv) == th

rd3 = (.2e-6)**3
assert common.rw3_cr(rd3, .5, 300) > rd3
