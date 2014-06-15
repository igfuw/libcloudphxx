import numpy as np
import math

from constants_pytest import es0, T0, Rv, Rd, cp, L, p0

#in Pa (4.23 C&W)
def press_sat(T, es_ref=es0, T_ref=T0):
    return es_ref * math.exp(-L/Rv * (1./T - 1./T_ref))

# in 1 (4.37 C&W)
def mixrat_sat(T, p, es_ref=es0, T_ref=T0):
    es = press_sat(T)
    return Rd/Rv * es / (p-es)

def pot_temp(T,p):
    return T * (p0 / p)**(Rd/cp)

def density_dry(rv, press, T):
    return press / T / Rd / (1. + rv * Rv/Rd)


def pottemp_dry(rv, press, T):
    return T * (p0 / press * (1. + rv * Rv/Rd))**(Rd/cp)

                

#linearization of: r_v - \delta r_v = r_vs(\theta + \delta \theta) and
# \delta \theta = L/c_p * \theta/T 
def delta_r(rv, T, p, es_ref=es0, T_ref=T0):
    rvs = mixrat_sat(T,p)
    return (rv - rvs) / (1. + rvs / T**2 * L**2 / Rv / cp)
    
if __name__ == "__main__":
    print mixrat_sat(283.15, 900.e2)
    print pot_temp(283.15, 900.e2)
    print press_sat(283.15, 1227., 283.15)
    print delta_r(10.e-3, 283.15, 900.e2)
    print delta_r(8.e-3, 283.15, 900.e2)
    print delta_r(9.e-3, 283.15, 900.e2)
