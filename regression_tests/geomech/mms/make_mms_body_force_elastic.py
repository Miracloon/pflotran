import math
import h5py
import numpy as np

# Structured grid from mms_elastic_3d.in
nx = 4
ny = 4
nz = 4
xmin, ymin, zmin = 0.0, 0.0, 0.0
xmax, ymax, zmax = 1.0, 1.0, 1.0

nxv = nx + 1
nyv = ny + 1
nzv = nz + 1
nxyv = nxv * nyv
nnode = nxv * nyv * nzv

xv = np.linspace(xmin, xmax, nxv)
yv = np.linspace(ymin, ymax, nyv)
zv = np.linspace(zmin, zmax, nzv)

# Elastic properties from input
E = 1.0e9
nu = 0.25
rho = 1000.0

lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
mu = E / (2.0 * (1.0 + nu))

pi = math.pi

body_x = np.zeros(nnode, dtype=np.float64)
body_y = np.zeros(nnode, dtype=np.float64)
body_z = np.zeros(nnode, dtype=np.float64)

# Natural ordering: vid = ix + iy*nxv + iz*nxyv (0-based)
idx = 0
for iz in range(nzv):
    z = zv[iz]
    sinz = math.sin(pi * z)
    cosz = math.cos(pi * z)
    for iy in range(nyv):
        y = yv[iy]
        siny = math.sin(pi * y)
        cosy = math.cos(pi * y)
        for ix in range(nxv):
            x = xv[ix]
            sinx = math.sin(pi * x)
            cosx = math.cos(pi * x)

            f = sinx * siny * sinz
            lap_f = -3.0 * (pi ** 2) * f

            ddiv_dx = (pi ** 2) * (-f + cosx * cosy * sinz + cosx * siny * cosz)
            ddiv_dy = (pi ** 2) * (-f + cosx * cosy * sinz + sinx * cosy * cosz)
            ddiv_dz = (pi ** 2) * (-f + cosx * siny * cosz + sinx * cosy * cosz)

            body_x[idx] = (-(lam + mu) * ddiv_dx - mu * lap_f) / rho
            body_y[idx] = (-(lam + mu) * ddiv_dy - mu * lap_f) / rho
            body_z[idx] = (-(lam + mu) * ddiv_dz - mu * lap_f) / rho

            idx += 1

cell_ids = np.arange(1, nnode + 1, dtype=np.int32)

with h5py.File("mms_body_force_elastic.h5", "w") as h5:
    h5.create_dataset("Cell Ids", data=cell_ids)
    h5.create_dataset("Body_force_x", data=body_x)
    h5.create_dataset("Body_force_y", data=body_y)
    h5.create_dataset("Body_force_z", data=body_z)
