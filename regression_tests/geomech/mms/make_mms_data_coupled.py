"""Generate MMS data for coupled flow-geomechanics verification.

Manufactured solution on [0,1]^3 (steady state):

  Pressure:     p = p0 + Ap * sin(pi*x)*sin(pi*y)*sin(pi*z)
  Displacement: u = Au * [sin(pi*x)*sin(pi*y)*sin(pi*z),
                          sin(pi*x)*sin(pi*y)*sin(pi*z),
                          sin(pi*x)*sin(pi*y)*sin(pi*z)]

  (all three displacement components identical, vanish on all faces)

Flow (Richards, steady state, constant rho_l, S_l=1, constant k, mu_l):
  Source Q per cell = rho_l * (k/mu_l) * 3*pi^2 * Ap * sin(pi*xc)*sin(pi*yc)*sin(pi*zc) * V_cell
  (in kg/s, for HETEROGENEOUS_MASS_RATE)

Geomechanics (quasi-static, linear isotropic elasticity + Biot coupling):
  Body force (acceleration, m/s^2) = (-div(sigma) + beta*grad(dp)) / rho_bulk

  where:
    div(sigma) = standard stress divergence from linear isotropic elasticity
    dp = Ap * sin(pi*x)*sin(pi*y)*sin(pi*z) (pressure increment = p - p_init)
    beta = Biot coefficient

  The Biot term is needed because PFLOTRAN's weak form is:
    K*u = integral(rho*b*N dV) + integral(beta*dp*dN/dx dV)
  In strong form: -div(sigma) = rho*b - beta*grad(dp)
  So: rho*b = -div(sigma) + beta*grad(dp)
"""

import argparse
import math

import h5py
import numpy as np


def generate_body_force(nx, ny, nz, E, nu, rho_s, rho_f, phi, Au, Ap, beta,
                        filename):
    """Generate geomechanics body force HDF5 (acceleration, m/s^2) at vertices.

    Body force = (-div(sigma(u_mms)) + beta*grad(dp_mms)) / rho_bulk

    where u_mms = Au*(f, f, f) with f = sin(pi*x)*sin(pi*y)*sin(pi*z)
    and dp_mms = Ap * sin(pi*x)*sin(pi*y)*sin(pi*z).
    """
    nxv, nyv, nzv = nx + 1, ny + 1, nz + 1
    nnode = nxv * nyv * nzv

    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    mu = E / (2.0 * (1.0 + nu))
    rho_bulk = phi * rho_f + (1.0 - phi) * rho_s

    xv = np.linspace(0.0, 1.0, nxv)
    yv = np.linspace(0.0, 1.0, nyv)
    zv = np.linspace(0.0, 1.0, nzv)

    pi = math.pi

    body_x = np.zeros(nnode, dtype=np.float64)
    body_y = np.zeros(nnode, dtype=np.float64)
    body_z = np.zeros(nnode, dtype=np.float64)

    idx = 0
    for iz in range(nzv):
        for iy in range(nyv):
            for ix in range(nxv):
                xc = xv[ix]
                yc = yv[iy]
                zc = zv[iz]

                # f = sin(pi*x)*sin(pi*y)*sin(pi*z) and its derivatives
                sx = math.sin(pi * xc)
                sy = math.sin(pi * yc)
                sz = math.sin(pi * zc)
                cx = math.cos(pi * xc)
                cy = math.cos(pi * yc)
                cz = math.cos(pi * zc)

                f = sx * sy * sz

                # Second derivatives of f
                s_xx = -(pi**2) * f  # d^2f/dx^2
                s_yy = -(pi**2) * f  # d^2f/dy^2
                s_zz = -(pi**2) * f  # d^2f/dz^2
                s_xy = (pi**2) * cx * cy * sz  # d^2f/dxdy
                s_xz = (pi**2) * cx * sy * cz  # d^2f/dxdz
                s_yz = (pi**2) * sx * cy * cz  # d^2f/dydz

                # Stress divergence components (same formula as mms_elastic_3d):
                # div(sigma)_x = (lam+2*mu)*f_xx + (lam+mu)*f_xy +
                #                (lam+mu)*f_xz + mu*f_yy + mu*f_zz
                div_x = ((lam + 2.0 * mu) * s_xx +
                         (lam + mu) * s_xy +
                         (lam + mu) * s_xz +
                         mu * s_yy + mu * s_zz)
                div_y = (mu * s_xx +
                         (lam + mu) * s_xy +
                         (lam + 2.0 * mu) * s_yy +
                         (lam + mu) * s_yz +
                         mu * s_zz)
                div_z = (mu * s_xx +
                         mu * s_yy +
                         (lam + mu) * s_xz +
                         (lam + mu) * s_yz +
                         (lam + 2.0 * mu) * s_zz)

                # Biot pressure gradient:
                # dp = Ap * sin(pi*x)*sin(pi*y)*sin(pi*z)
                # grad(dp)_x = Ap * pi * cos(pi*x)*sin(pi*y)*sin(pi*z)
                # grad(dp)_y = Ap * pi * sin(pi*x)*cos(pi*y)*sin(pi*z)
                # grad(dp)_z = Ap * pi * sin(pi*x)*sin(pi*y)*cos(pi*z)
                grad_dp_x = Ap * pi * cx * sy * sz
                grad_dp_y = Ap * pi * sx * cy * sz
                grad_dp_z = Ap * pi * sx * sy * cz

                # Body force = (-div(sigma) + beta*grad(dp)) / rho_bulk
                body_x[idx] = (-div_x + beta * grad_dp_x) / rho_bulk
                body_y[idx] = (-div_y + beta * grad_dp_y) / rho_bulk
                body_z[idx] = (-div_z + beta * grad_dp_z) / rho_bulk

                idx += 1

    cell_ids = np.arange(1, nnode + 1, dtype=np.int32)

    with h5py.File(filename, "w") as h5:
        h5.create_dataset("Cell Ids", data=cell_ids)
        h5.create_dataset("Body_force_x", data=body_x)
        h5.create_dataset("Body_force_y", data=body_y)
        h5.create_dataset("Body_force_z", data=body_z)


def generate_flow_source(nx, ny, nz, k, mu_l, rho_l, Ap, filename):
    """Generate flow source term HDF5 (HETEROGENEOUS_MASS_RATE, kg/s per cell).

    The MMS source for steady-state Richards with constant properties:
      Q_density = rho_l * (k/mu_l) * 3*pi^2 * Ap * sin(pi*xc)*sin(pi*yc)*sin(pi*zc)
    Per cell:
      Q_cell = Q_density * V_cell   [kg/s]
    """
    ncells = nx * ny * nz
    dx, dy, dz = 1.0 / nx, 1.0 / ny, 1.0 / nz
    V_cell = dx * dy * dz

    # Cell center coordinates
    xc = np.linspace(dx / 2, 1.0 - dx / 2, nx)
    yc = np.linspace(dy / 2, 1.0 - dy / 2, ny)
    zc = np.linspace(dz / 2, 1.0 - dz / 2, nz)

    pi = math.pi
    Kd = k / mu_l  # hydraulic diffusivity m^2/(Pa.s)
    coef = rho_l * Kd * 3.0 * pi**2 * Ap * V_cell

    # Source rate per cell (kg/s)
    source_rate = np.zeros(ncells, dtype=np.float64)

    # Natural ordering: cell_id = ix + iy*nx + iz*nx*ny (0-based), then +1
    idx = 0
    for iz in range(nz):
        sz = math.sin(pi * zc[iz])
        for iy in range(ny):
            sy = math.sin(pi * yc[iy])
            for ix in range(nx):
                sx = math.sin(pi * xc[ix])
                source_rate[idx] = coef * sx * sy * sz
                idx += 1

    # DATASET MAPPED format:
    #   Data group: 1D float64 array of source rates
    #   Map group:  2D int32 array.
    #
    #   HDF5/C shape must be (N, 2) so that Fortran reads it as (2, N),
    #   where mapping(1,:) = data_indices (1-based) and
    #         mapping(2,:) = cell_ids (1-based).
    cell_ids = np.arange(1, ncells + 1, dtype=np.int32)
    data_indices = np.arange(1, ncells + 1, dtype=np.int32)
    # Stack as columns: shape (N, 2) in C order → (2, N) in Fortran
    map_array = np.stack([data_indices, cell_ids], axis=1).astype(np.int32)

    with h5py.File(filename, "w") as h5:
        # Data group
        grp_data = h5.create_group("Flow_Source")
        grp_data.create_dataset("Data", data=source_rate)

        # Map group
        grp_map = h5.create_group("Flow_Source_Map")
        grp_map.create_dataset("Data", data=map_array)


def main():
    parser = argparse.ArgumentParser(
        description="Generate MMS data for coupled flow-geomechanics test."
    )
    parser.add_argument("--nx", type=int, default=4, help="Number of cells in x")
    parser.add_argument("--ny", type=int, default=4, help="Number of cells in y")
    parser.add_argument("--nz", type=int, default=4, help="Number of cells in z")
    args = parser.parse_args()

    nx, ny, nz = args.nx, args.ny, args.nz

    # ---- Material / MMS parameters ----
    # Elastic properties
    E = 1.0e9        # Young's modulus [Pa]
    nu = 0.25        # Poisson's ratio
    rho_s = 2500.0   # Rock density [kg/m^3]
    rho_f = 997.16   # Fluid density [kg/m^3]
    phi = 0.2        # Porosity

    # Biot coupling
    beta = 1.0       # Biot coefficient

    # Flow properties
    k = 1.0e-12      # Permeability [m^2]
    mu_l = 1.0e-3    # Fluid viscosity [Pa.s]
    rho_l = rho_f    # Fluid density (constant) [kg/m^3]

    # MMS amplitudes
    Au = 1.0         # Displacement amplitude [m]
    Ap = 1000.0      # Pressure perturbation amplitude [Pa]

    # ---- Generate files ----
    generate_body_force(nx, ny, nz, E, nu, rho_s, rho_f, phi, Au, Ap, beta,
                        "mms_body_force_coupled.h5")
    generate_flow_source(nx, ny, nz, k, mu_l, rho_l, Ap,
                         "mms_flow_source_coupled.h5")

    print(f"Generated MMS data for {nx}x{ny}x{nz} grid:")
    print(f"  mms_body_force.h5  ({(nx+1)*(ny+1)*(nz+1)} vertices)")
    print(f"  mms_flow_source.h5 ({nx*ny*nz} cells)")


if __name__ == "__main__":
    main()
