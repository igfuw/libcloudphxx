import matplotlib.pyplot as plt
from scipy.io import netcdf
import numpy as np
import pdb

f_out = netcdf.netcdf_file("test.nc", "r")
z = f_out.variables["z"][:]
RH = f_out.variables["RH"][:]
print type(f_out.bins)
print "RH_max", f_out.RH_max, RH.max()

RH_max_l = [f_out.RH_max] * z.shape[0]

plt.figure(1, figsize=(6,6))
ax = plt.subplot(1,1,1)
ax.plot(RH, z, RH_max_l, z)
plt.show()
