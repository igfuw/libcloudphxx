from cffi import FFI
import numpy as np
from stale import xlv, cp, EP2, SVP1, SVP2, SVP3, SVPT0, rhowater
# TODO: should check all above!

ffi = FFI()

def kessler(nx, ny, nz, dt_in,
            t_np, qv_np, qc_np, qr_np, rho_np, pii_np, dz8w_np, z_np,
            rainnc_np, rainncv_np):

    ffi.cdef("void c_kessler(float t[], float qv[], float qc[], float qr[], float rho[], float pii[], float dt_in, float z[], float xlv, float cp, float EP2, float SVP1, float SVP2, float SVP3, float SVPT0, float rhowater, float dz8w[], float RAINNC[], float RAINNCV[], int ids, int ide, int jds, int jde, int kds, int kde, int ims, int ime, int jms, int jme, int kms, int kme, int its, int ite, int jts, int jte, int kts, int kt);")

    lib = ffi.dlopen('libkessler.so')

    #TODO: add checking if the t_np etc. are float32! 
    t = ffi.cast("float*", t_np.__array_interface__['data'][0])
    qv = ffi.cast("float*", qv_np.__array_interface__['data'][0])
    qc = ffi.cast("float*", qc_np.__array_interface__['data'][0])
    qr = ffi.cast("float*", qr_np.__array_interface__['data'][0])
    rho = ffi.cast("float*", rho_np.__array_interface__['data'][0]) 
    pii = ffi.cast("float*", pii_np.__array_interface__['data'][0])
    dz8w = ffi.cast("float*", dz8w_np.__array_interface__['data'][0])
    z = ffi.cast("float*", z_np.__array_interface__['data'][0])
    RAINNC = ffi.cast("float*", rainnc_np.__array_interface__['data'][0]) 
    RAINNCV = ffi.cast("float*", rainncv_np.__array_interface__['data'][0])

#its etc. is not actually used; not sure why its etc. is used but can (should?) be the same as ids etc.
    [ims, ime, ids, ide, its, ite] = [1, nx] * 3

    [jms, jme, jds, jde, jts, jte] = [1, ny] * 3

    [kms, kme, kds, kde, kts, kte] = [1, nz] * 3

    lib.c_kessler(t, qv, qc, qr, rho, pii                 
                  ,dt_in, z, xlv, cp                        
                  ,EP2,SVP1,SVP2,SVP3,SVPT0,rhowater        
                  ,dz8w                                     
                  ,RAINNC, RAINNCV                          
                  ,ids,ide, jds,jde, kds,kde                 
                  ,ims,ime, jms,jme, kms,kme                 
                  ,its,ite, jts,jte, kts,kte                 
                  )
