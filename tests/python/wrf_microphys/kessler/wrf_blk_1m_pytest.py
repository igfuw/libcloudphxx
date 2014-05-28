import numpy as np
from cffi_kessler import kessler
from stale import Rd, cp
import pdb

nx = 1
ny = 1
nz = 1

dz8w = np.ones((nx,nz,ny), "float32") * 20.
z = np.ones((nx,nz,ny), "float32") * 700.
rainnc = np.zeros((nx,ny), "float32")
rainncv = np.zeros((nx,ny), "float32")

def exner(press):
    return np.array((press/1.e5)**(Rd/cp), "float32")

def density(press, th):
    #TODO check, if it's full rho
    return np.ones((nx,nz,ny), "float32") * 1.

def adj_cellwise(press_in, th_in, qv_in, qc_in, qr_in, dt):
#    pdb.set_trace()
    pii = exner(press_in)
    rho = density(press_in, th_in)
    th = np.array(th_in, "float32")
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
