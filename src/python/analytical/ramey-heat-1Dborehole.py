#==============================================================================
# Authors: 
#     Piyoosh Jaysaval, PNNL (piyoosh.jaysaval@pnnl.gov)
# Date: February 2026
#
# Description: Script to solve Ramey (1962) analytical solution for temperature
#              with to depth and time in borehole fluid and heat flow.
#
# Implementation follows:
#              1. Chen et al. (2021): Numerical investigation on the capacity
#                  and efficiency of a deep enhanced U-tube borehole heat
#                  exchanger system for building heating, Renewable Energy
#              2. Ramey et al. (1962) Wellbore heat transmission,
#                  J. Petrol. Technol., 14, 427-435.
#              3. https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/
#                            Data/Parabolic/T/BHE_1P/pipe_flow_ebhe.py
#==============================================================================

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Union, Dict

import numpy as np
import matplotlib.pyplot as plt

Number = Union[int, float]
ArrayLike = Union[np.ndarray, Number]


@dataclass
class RameyInputs:
    # Temperatures (°C)
    T_d: float = 55.0   # formation temperature
    T_i: float = 20.0   # inlet fluid temperature

    # Geometry (m)
    L: float = 30.0
    r_pi: float = 0.12913     # inner radius of pipe
    t_pi: float = 0.00587     # pipe wall thickness
    r_b: float = 0.14         # borehole radius

    # Flow rate (m^3/s)
    q: float = 1.99e-4        # Area x v = pi*r_pi^2 * 0.0038 m/s --> 1.99e-4

    # Fluid properties
    rho_f: float = 1000.0     # kg/m^3
    c_p_f: float = 4190.0     # J/kg/K
    mu_f: float = 1.14e-3     # Pa*s
    lambda_f: float = 0.59    # W/m/K

    # Rock properties
    rho_re: float = 1800.0    # kg/m^3
    c_p_re: float = 1778.0    # J/kg/K
    lambda_re: float = 2.78018   # W/m/K

    # Grout + pipe conductivities (W/m/K)
    lambda_g: float = 0.73
    lambda_pi: float = 1.3

    # Numerical
    eps: float = 1e-12

    # Choice of constant Nu model
    # "laminar" -> Nu = 3.66 (Ramey-style fully-developed)
    # "ogs_like" -> OGS-style correlation but with NO z dependence
    nu_mode: str = "ogs_like"   # or "laminar"
    Nu_laminar: float = 3.66    # could also use 4.364 if desired

    # If we *want* to keep OGS entrance-length factor but as constant:
    # set z_ref_for_Nu (e.g., L) and we evaluate (r_pi/z_ref)^(2/3) once.
    use_entrance_factor: bool = False
    z_ref_for_Nu: float = 30.0


@dataclass
class RameyModel:
    params: RameyInputs = field(default_factory=RameyInputs)
    _U_cache: Optional[float] = field(default=None, init=False, repr=False)
    _Nu_cache: Optional[float] = field(default=None, init=False, repr=False)

    def update(self, **kwargs: Dict[str, Number]) -> None:
        for k, v in kwargs.items():
            if not hasattr(self.params, k):
                raise ValueError(f"Unknown parameter: {k}")
            setattr(self.params, k, float(v))
        self._U_cache = None
        self._Nu_cache = None

    # -------------------------------
    # Dimensionless groups Eq. A.10
    # -------------------------------
    def prandtl(self) -> float:
        p = self.params
        return (p.mu_f * p.c_p_f) / max(p.lambda_f, p.eps)

    def reynolds(self) -> float:
        p = self.params
        v = p.q / (np.pi * p.r_pi**2)
        D = 2.0 * p.r_pi
        return (p.rho_f * v * D) / max(p.mu_f, p.eps)

    # -------------------------------
    # Ramey time function Eq. A.3-A.4
    # -------------------------------
    def t_dimless(self, t_s: ArrayLike) -> np.ndarray:
        p = self.params
        t = np.asarray(t_s, dtype=float)
        denom = max(p.rho_re * p.c_p_re * p.r_b**2, p.eps)
        return (p.lambda_re * t) / denom

    def f_time(self, t_s: ArrayLike) -> np.ndarray:
        p = self.params
        tD = np.maximum(self.t_dimless(t_s), p.eps)

        f = np.empty_like(tD)
        mask = tD > 1.5
        f[mask] = (0.4063 + 0.5 * np.log(tD[mask])) * (1.0 + 0.6 / tD[mask])

        rt = np.sqrt(tD[~mask])
        f[~mask] = 1.1281 * rt * (1.0 - 0.3 * rt)
        return f

    # -------------------------------
    # Constant Nu, h, U Eq. A.7 - A.9
    # -------------------------------
    def nusselt_const(self, force_recompute: bool = False) -> float:
        """
        Returns a SINGLE scalar Nu.
        """
        p = self.params
        if (self._Nu_cache is not None) and (not force_recompute):
            return float(self._Nu_cache)

        Re = self.reynolds()
        Pr = self.prandtl()

        if p.nu_mode.lower() == "laminar":
            Nu = float(p.Nu_laminar)

        elif p.nu_mode.lower() == "ogs_like":
            # OGS notebook uses a blended correlation but includes a
            # z-dependent entrance factor.
            # Here we REMOVE that dependence. Optionally,
            # we can evaluate that factor once at z_ref.
            xi = (1.8 * np.log(max(Re, 1.0)) - 1.5) ** (-2)

            base = (
                (xi / 8.0 * Re * Pr)
                / (1.0 + 12.7 * np.sqrt(xi / 8.0) * (Pr ** (2.0 / 3.0) - 1.0))
            )

            if p.use_entrance_factor:
                zref = max(p.z_ref_for_Nu, 1e-10)
                entrance = 1.0 + (p.r_pi / zref) ** (2.0 / 3.0)
            else:
                entrance = 1.0

            Nu_turb = base * entrance
            Nu_lam = 4.364

            if Re < 2300.0:
                Nu = float(Nu_lam)
            elif Re > 10000.0:
                Nu = float(Nu_turb)
            else:
                gamma = (Re - 2300.0) / (10000.0 - 2300.0)
                gamma = float(np.clip(gamma, 0.0, 1.0))

                # “fixed” transitional blend consistent with their structure
                Nu_at_1e4 = (
                    (0.0308 / 8.0 * 1.0e4 * Pr)
                    / (1.0 + 12.7 * np.sqrt(0.0308 / 8.0)
                       * (Pr ** (2.0 / 3.0) - 1.0))
                ) * entrance

                Nu = float((1.0 - gamma) * Nu_lam + gamma * Nu_at_1e4)
        else:
            raise ValueError("nu_mode must be 'laminar' or 'ogs_like'")

        self._Nu_cache = Nu
        return Nu

    def h_inner_const(self) -> float:
        # Chen et al. (2021) Eq. A.6
        p = self.params
        Nu = self.nusselt_const()
        return (p.lambda_f * Nu) / max(2.0 * p.r_pi, p.eps)

    def U_const(self, force_recompute: bool = False) -> float:
        """
        Single scalar U based on single scalar Nu.
        Chen et al. (2021) Eq. A.5
        """
        p = self.params
        if (self._U_cache is not None) and (not force_recompute):
            return float(self._U_cache)

        h = max(self.h_inner_const(), p.eps)

        r_i = p.r_pi
        r_o = p.r_pi + p.t_pi

        term_conv = r_o / (max(r_i, p.eps) * h)
        term_cond = r_o * (
            np.log(r_o / max(r_i, p.eps)) / max(p.lambda_pi, p.eps)
            + np.log(max(p.r_b, p.eps) / max(r_o, p.eps))
            / max(p.lambda_g, p.eps)
        )

        U = 1.0 / max(term_conv + term_cond, p.eps)
        self._U_cache = float(U)
        return float(U)

    # -------------------
    # Ramey X(t) and T(z,t)
    # -------------------
    def X(self, t_s: ArrayLike, U: Optional[float] = None) -> np.ndarray:
        """
        X(t) = (q rho_f cp_f) * (lambda_re + r_pi U f(t))
                                               / (2 pi r_pi U lambda_re)
        Chen et al. (2021) Eq. A.2
        """
        p = self.params
        t = np.asarray(t_s, dtype=float)
        f = self.f_time(t)

        if U is None:
            U = self.U_const()

        num = (p.q * p.rho_f * p.c_p_f) * (p.lambda_re + p.r_pi * U * f)
        den = (2.0 * np.pi * p.r_pi * U * p.lambda_re)
        return num / np.maximum(den, p.eps)

    def temperature(self, z: ArrayLike, t_s: ArrayLike,
                    U: Optional[float] = None) -> np.ndarray:
        """
        T(z,t) = T_d + (T_i - T_d) * exp(-z / X(t))
        Chen et al. (2021) Eq. A.1
        """
        p = self.params
        z = np.asarray(z, dtype=float)
        t = np.asarray(t_s, dtype=float)

        Xt = self.X(t, U=U)  # (Nt,)
        zz = z.reshape((-1, 1)) if z.ndim == 1 else z
        Xtt = Xt.reshape((1, -1))
        return p.T_d + (p.T_i - p.T_d) * np.exp(-zz / np.maximum(Xtt, p.eps))

    # -------------------
    # Plotting
    # -------------------
    def plot_outlet_vs_time(self, t_s: np.ndarray) -> None:
        p = self.params
        U0 = self.U_const()
        Tout = self.temperature(p.L, t_s, U=U0).reshape(-1)
        t_days = t_s / 86400.0

        plt.figure()
        plt.plot(t_days, Tout,'o-',color='red')
        plt.xlabel("Time (days)")
        plt.ylabel("Outlet temperature (°C)")
        plt.title("Ramey analytical: outlet temperature vs time")
        plt.grid(True, which="both", alpha=0.3)
        plt.tight_layout()

    def plot_profiles_vs_depth(self, z: np.ndarray, t_s: np.ndarray,
                               times_to_show_s: np.ndarray) -> None:
        U0 = self.U_const()
        T = self.temperature(z, t_s, U=U0)

        plt.figure()
        for ts in times_to_show_s:
            j = int(np.argmin(np.abs(t_s - ts)))
            plt.plot(z, T[:, j], label=f"t={t_s[j]/86400:.2f} d")
        plt.xlabel("Depth Z (m)")
        plt.ylabel("Temperature (degree C)")
        plt.title("Ramey analytical: temperature vs depth")
        plt.grid(True, which="both", alpha=0.3)
        plt.legend()
        plt.tight_layout()

    def plot_heatmap_T_of_z_t(self, z: np.ndarray, t_s: np.ndarray) -> None:
        U0 = self.U_const()
        T = self.temperature(z, t_s, U=U0)
        t_days = t_s / 86400.0

        plt.figure()
        plt.imshow(
            T,
            aspect="auto",
            origin="lower",
            extent=[t_days.min(), t_days.max(), z.min(), z.max()],
        )
        plt.xlabel("Time (days)")
        plt.ylabel("Depth Z (m)")
        plt.title("Ramey analytical: T(Z,t) heatmap")
        plt.colorbar(label="Temperature (°C)")
        plt.tight_layout()

    def plot_contours_T_of_z_t(self, z: np.ndarray, t_s: np.ndarray,
                               n_levels: int = 12) -> None:
        U0 = self.U_const()
        T = self.temperature(z, t_s, U=U0)
        t_days = t_s / 86400.0
        TT, ZZ = np.meshgrid(t_days, z)

        plt.figure()
        cs = plt.contour(TT, ZZ, T, levels=n_levels)
        plt.clabel(cs, inline=True, fontsize=8)
        plt.xlabel("Time (days)")
        plt.ylabel("Depth Z (m)")
        plt.title("Ramey analytical: T(Z,t) contours")
        plt.grid(True, which="both", alpha=0.2)
        plt.tight_layout()


if __name__ == "__main__":
    model = RameyModel()

    # Depth and time grids
    Z = np.insert(np.linspace(0.3, model.params.L, 100), 0, 1e-10)
    t_end = 30 * 86400  # 30-day solution
    t_s = np.arange(0, t_end + 1, 43200)  # every 12 hours

    Nu0 = model.nusselt_const()
    U0 = model.U_const()
    print(f"Using constant Nu = {Nu0:.6g}")
    print(f"Using constant U  = {U0:.6g} W/m^2/K")

    model.plot_outlet_vs_time(t_s)
    model.plot_profiles_vs_depth(
        Z, t_s,
        times_to_show_s=np.array([0, 1*86400, 3*86400, 5*86400, 10*86400,
                                  30*86400])
    )
    model.plot_heatmap_T_of_z_t(Z, t_s)
    model.plot_contours_T_of_z_t(Z, t_s, n_levels=12)

    plt.show()
