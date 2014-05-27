import sys
sys.path.append(".")
import numpy as np
from cffi_kessler import kessler

nx = 1
ny = 1
nz = 1
dt_in = 1.

t_np = np.ones((nx,nz,ny), "float32") * 291.8
qv_np = np.ones((nx,nz,ny), "float32") * 10.e-3
qc_np = np.ones((nx,nz,ny), "float32") * 0.e-3
qr_np = np.zeros((nx,nz,ny), "float32")
rho_np = np.ones((nx,nz,ny), "float32")
pii_np = np.ones((nx,nz,ny), "float32") * .97
dz8w_np = np.ones((nx,nz,ny), "float32") * 20.
z_np = np.ones((nx,nz,ny), "float32") * 700.
rainnc_np = np.zeros((nx,ny), "float32")
rainncv_np = np.zeros((nx,ny), "float32")

print "qv_przed", qv_np, qc_np
kessler(nx, ny, nz, dt_in,
            t_np, qv_np, qc_np, qr_np, rho_np, pii_np, dz8w_np, z_np,
            rainnc_np, rainncv_np)
print "qv po", qv_np, qc_np, qr_np
