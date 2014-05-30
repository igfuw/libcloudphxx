import numpy as np
from cffi_kessler import kessler
from stale import Rd, cp #TODO ujednolicic stale
import pdb

#TODO
Rv = 461.4
p0 = 1000.e2

nx = 1
ny = 1
nz = 1

dz8w = np.ones((nx,nz,ny), "float32") * 20.
z = np.ones((nx,nz,ny), "float32") * 700.
rainnc = np.zeros((nx,ny), "float32")
rainncv = np.zeros((nx,ny), "float32")

def exner(press):
    return np.array((press/1.e5)**(Rd/cp), "float32")

def density(rv, press, T):
    #TODO check, if it's full rho
    R_tot = (Rd + rv*Rv) / (1. + rv)
    rho = press / R_tot / T
    return np.array(rho, "float32")

def pottemp(press, T):
     theta = T * (p0 / press)**(Rd/cp) 
     return np.array(theta, "float32")

def adj_cellwise(press_in, T_in, qv_in, qc_in, qr_in, dt):
#    pdb.set_trace()
    pii = exner(press_in)
    rho = density(qv_in, press_in, T_in)
    th = pottemp(press_in, T_in)
    qv = np.array(qv_in, "float32")
    qc = np.array(qc_in, "float32")
    qr = np.array(qr_in, "float32")
    print "qv_przed", qv, qc, pii, rho
    kessler(nx, ny, nz, dt,
            th, qv, qc, qr, rho, pii, dz8w, z,
            rainnc, rainncv)
    print "qv po", qv, qc, qr
    return qv, qc

#adj_cellwise(np.array([9.e4]), t_np, qv_np, qc_np, qr_np, 1)
