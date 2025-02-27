#
# comparison.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# comparison.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Compare CCL and NumCosmo results."""

import timeit
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
import pyccl
import numcosmo_py.cosmology as ncc
from numcosmo_py.plotting.tools import latex_float
from numcosmo_py.helper import npa_to_seq


class CompareFunc1d:
    """Class for functions comparison."""

    def __init__(
        self,
        x: npt.NDArray[np.float64],
        y1: npt.NDArray[np.float64],
        y2: npt.NDArray[np.float64],
        *,
        model: str = "unnamed",
        name1: str = "1",
        name2: str = "2",
        x_symbol: str = "x",
        y_symbol: str = "y",
        x_unit: str | None = None,
        y_unit: str | None = None,
    ) -> None:
        """Compare CCL and NumCosmo results for a function."""
        assert x.shape == y1.shape
        assert x.shape == y2.shape
        self.x = x
        self.y1 = y1
        self.y2 = y2
        self.model = model
        self.name1 = name1
        self.name2 = name2
        self.x_symbol = x_symbol
        self.y_symbol = y_symbol
        self.x_unit = x_unit
        self.y_unit = y_unit
        self.diff = np.zeros_like(y1)

        non_zero_indices = y1 != 0.0
        zero_indices = ~non_zero_indices
        self.diff[non_zero_indices] = y2[non_zero_indices] / y1[non_zero_indices] - 1.0
        self.diff[zero_indices] = y2[zero_indices] - y1[zero_indices]
        self.abs_diff = np.abs(self.diff)

    @property
    def rel_diff_min(self) -> float:
        """Return the minimum relative difference."""
        return np.min(self.abs_diff)

    @property
    def rel_diff_max(self) -> float:
        """Return the maximum relative difference."""
        return np.max(self.abs_diff)

    @property
    def rel_diff_mean(self) -> float:
        """Return the mean relative difference."""
        return np.mean(self.abs_diff)

    @property
    def rel_diff_std(self) -> float:
        """Return the standard deviation of the relative difference."""
        return np.std(self.abs_diff)

    @property
    def rel_diff_median(self) -> float:
        """Return the median relative difference."""
        return np.median(self.abs_diff)

    @property
    def x_label(self) -> str:
        """Return the label for the x axis."""
        if self.x_unit:
            return f"${self.x_symbol}$ [{self.x_unit}]"
        else:
            return f"${self.x_symbol}$"

    @property
    def y_label(self) -> str:
        """Return the label for the y axis."""
        if self.y_unit:
            return f"${self.y_symbol}$ [{self.y_unit}]"
        else:
            return f"${self.y_symbol}$"

    def plot(
        self,
        axs: list[plt.Axes],
        *,
        xscale: str = "linear",
        yscale: str = "log",
        color: str = "black",
        lw: float = 0.8,
    ) -> None:
        """Plot the comparison."""
        assert len(axs) == 2
        assert isinstance(axs[0], plt.Axes)
        assert isinstance(axs[1], plt.Axes)
        y1_name, y2_name = (
            f"{self.y_symbol}^\\mathrm{{{self.name1}}}",
            f"{self.y_symbol}^\\mathrm{{{self.name2}}}",
        )
        axs[0].plot(
            self.x, self.y1, label=self.model, lw=lw, color=color, linestyle="-"
        )
        axs[0].plot(self.x, self.y2, lw=lw, color=color, linestyle="--")
        if np.sum(self.abs_diff) > 0.0:
            axs[1].plot(self.x, self.abs_diff, lw=lw, color=color)
        axs[1].set_xscale(xscale)
        axs[1].set_yscale(yscale)
        axs[0].legend(loc="best")
        axs[0].set_ylabel(self.y_label)
        axs[1].set_xlabel(self.x_label)
        axs[1].set_ylabel(f"${y1_name}/{y2_name}-1$")

    def summary_row(self, precision: int = 2, convert_g: bool = True) -> list[str]:
        """Return a summary row."""
        rel_diff_min = latex_float(
            self.rel_diff_min, precision=precision, convert_g=convert_g
        )
        rel_diff_max = latex_float(
            self.rel_diff_max, precision=precision, convert_g=convert_g
        )
        rel_diff_mean = latex_float(
            self.rel_diff_mean, precision=precision, convert_g=convert_g
        )
        rel_diff_std = latex_float(
            self.rel_diff_std, precision=precision, convert_g=convert_g
        )
        return [
            self.model,
            self.y_label,
            f"${rel_diff_min}$",
            f"${rel_diff_max}$",
            f"${rel_diff_mean}$",
            f"${rel_diff_std}$",
        ]

    @classmethod
    def table_header(cls) -> list[str]:
        """Return the table header."""
        return [
            "Model",
            "Quantity",
            "$\\Delta P_{\\text{min}}/P$",
            "$\\Delta P_{\\text{max}}/P$",
            "$\\overline{\\Delta P / P}$",
            "$\\sigma_{\\Delta P / P}$",
        ]


def compute_times(func, repeat=5, number=100):
    """Compute times for a function."""
    times = timeit.repeat(func, repeat=repeat, number=number)

    return [np.mean(times) / number, np.std(times) / number]


def compare_Omega_m(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
) -> CompareFunc1d:
    """Compare Omega_m from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    a = 1.0 / (1.0 + z)

    ccl_Om = np.array(pyccl.omega_x(ccl_cosmo, a, "matter"))
    # NumCosmo defines Omega_m including only non-radiation part of the neutrino density
    # (rho_mnu - 3 p_mnu), as it is done in CLASS, however, CCL uses the full neutrino
    # density. Thus, since CCL does not compute p_mnu, we recompute Omega_m from
    # NumCosmo in the same way as CCL does.
    nc_Oc = np.array([cosmo.E2Omega_c(z_i) / cosmo.E2(z_i) for z_i in z])
    nc_Ob = np.array([cosmo.E2Omega_b(z_i) / cosmo.E2(z_i) for z_i in z])
    nc_Omnu = np.array([cosmo.E2Omega_mnu(z_i) / cosmo.E2(z_i) for z_i in z])
    nc_Om = nc_Oc + nc_Ob + nc_Omnu

    return CompareFunc1d(
        x=z, y1=ccl_Om, y2=nc_Om, model=model, x_symbol="z", y_symbol=r"\Omega_m"
    )


def compare_Omega_g(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare Omega_g from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    a = 1.0 / (1.0 + z)

    ccl_Og = np.array(pyccl.omega_x(ccl_cosmo, a, "radiation"))
    nc_Og = np.array([cosmo.E2Omega_g(z_i) / cosmo.E2(z_i) for z_i in z])

    return CompareFunc1d(
        x=z, y1=ccl_Og, y2=nc_Og, model=model, x_symbol="z", y_symbol=r"\Omega_g"
    )


def compare_Omega_k(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare Omega_k from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    a = 1.0 / (1.0 + z)

    ccl_Ok = np.array(pyccl.omega_x(ccl_cosmo, a, "curvature"))
    nc_Ok = np.array([cosmo.E2Omega_k(z_i) / cosmo.E2(z_i) for z_i in z])

    return CompareFunc1d(
        x=z, y1=ccl_Ok, y2=nc_Ok, model=model, x_symbol="z", y_symbol=r"\Omega_k"
    )


def compare_Omega_nu(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare Omega_nu from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    a = 1.0 / (1.0 + z)

    ccl_Onu = np.array(pyccl.omega_x(ccl_cosmo, a, "neutrinos_rel"))
    nc_Onu = np.array([cosmo.E2Omega_nu(z_i) / cosmo.E2(z_i) for z_i in z])

    return CompareFunc1d(
        x=z, y1=ccl_Onu, y2=nc_Onu, model=model, x_symbol="z", y_symbol=r"\Omega_{\nu}"
    )


def compare_Omega_mnu(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare Omega_mnu from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    a = 1.0 / (1.0 + z)

    ccl_Omnu = np.array(pyccl.omega_x(ccl_cosmo, a, "neutrinos_massive"))
    nc_Omnu = np.array([cosmo.E2Omega_mnu(z_i) / cosmo.E2(z_i) for z_i in z])

    return CompareFunc1d(
        x=z,
        y1=ccl_Omnu,
        y2=nc_Omnu,
        model=model,
        x_symbol="z",
        y_symbol=r"\Omega_{m\nu}",
    )


def compare_Omega_de(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare Omega_de from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    a = 1.0 / (1.0 + z)

    ccl_Ode = np.array(pyccl.omega_x(ccl_cosmo, a, "dark_energy"))
    nc_Ode = np.array([cosmo.E2Omega_de(z_i) / cosmo.E2(z_i) for z_i in z])

    return CompareFunc1d(
        x=z,
        y1=ccl_Ode,
        y2=nc_Ode,
        model=model,
        x_symbol="z",
        y_symbol=r"\Omega_\mathrm{DE}",
    )


def compare_Hubble(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare H(z)/H_0 from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    a = 1.0 / (1.0 + z)

    ccl_E = np.array(pyccl.h_over_h0(ccl_cosmo, a))
    nc_E = np.array([cosmo.E(z_i) for z_i in z])

    return CompareFunc1d(
        x=z, y1=ccl_E, y2=nc_E, model=model, x_symbol="z", y_symbol=r"H(z)/H_0"
    )


# Distances comparison


def compare_distance_comoving(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare comoving distance from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    a = 1.0 / (1.0 + z)

    RH_Mpc = cosmo.RH_Mpc()

    ccl_D = np.array(pyccl.comoving_radial_distance(ccl_cosmo, a))
    nc_D = (
        np.array(dist.comoving_array(cosmo, npa_to_seq(z)), dtype=np.float64) * RH_Mpc
    )

    return CompareFunc1d(
        x=z, y1=ccl_D, y2=nc_D, model=model, x_symbol="z", y_symbol=r"D_c", y_unit="Mpc"
    )


def compare_distance_transverse(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare transverse distance from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    a = 1.0 / (1.0 + z)

    RH_Mpc = cosmo.RH_Mpc()

    ccl_D = np.array(pyccl.comoving_angular_distance(ccl_cosmo, a))
    nc_D = (
        np.array(dist.transverse_array(cosmo, npa_to_seq(z)), dtype=np.float64) * RH_Mpc
    )

    return CompareFunc1d(
        x=z, y1=ccl_D, y2=nc_D, model=model, x_symbol="z", y_symbol=r"D_t", y_unit="Mpc"
    )


def compare_distance_angular_diameter(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare angular diameter distance from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    a = 1.0 / (1.0 + z)

    RH_Mpc = cosmo.RH_Mpc()

    ccl_D = np.array(pyccl.angular_diameter_distance(ccl_cosmo, a))
    nc_D = (
        np.array(dist.angular_diameter_array(cosmo, npa_to_seq(z)), dtype=np.float64)
        * RH_Mpc
    )

    return CompareFunc1d(
        x=z, y1=ccl_D, y2=nc_D, model=model, x_symbol="z", y_symbol=r"D_A", y_unit="Mpc"
    )


def compare_distance_luminosity(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare luminosity distance from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    a = 1.0 / (1.0 + z)

    RH_Mpc = cosmo.RH_Mpc()

    ccl_D = np.array(pyccl.luminosity_distance(ccl_cosmo, a))
    nc_D = (
        np.array(dist.luminosity_array(cosmo, npa_to_seq(z)), dtype=np.float64) * RH_Mpc
    )

    return CompareFunc1d(
        x=z, y1=ccl_D, y2=nc_D, model=model, x_symbol="z", y_symbol=r"D_L", y_unit="Mpc"
    )


def compare_distance_modulus(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare modulus distance from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    a = 1.0 / (1.0 + z)

    RH_Mpc = cosmo.RH_Mpc()

    ccl_D = np.array(pyccl.distance_modulus(ccl_cosmo, a))
    nc_D = np.array(dist.dmodulus_array(cosmo, npa_to_seq(z))) + 5 * np.log10(RH_Mpc)

    return CompareFunc1d(
        x=z, y1=ccl_D, y2=nc_D, model=model, x_symbol="z", y_symbol=r"\mu"
    )


def compare_distance_lookback_time(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare lookback time from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    a = 1.0 / (1.0 + z)

    t_H = (
        pyccl.physical_constants.MPC_TO_METER
        / 1.0e14
        / pyccl.physical_constants.YEAR
        / ccl_cosmo["h"]
    )
    ccl_D = np.array(pyccl.lookback_time(ccl_cosmo, a))
    nc_D = np.array(
        [dist.lookback_time(cosmo, z_i) * t_H for z_i in z], dtype=np.float64
    )

    return CompareFunc1d(
        x=z, y1=ccl_D, y2=nc_D, model=model, x_symbol="z", y_symbol=r"t_L", y_unit="Gyr"
    )


def compare_distance_comoving_volume(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare comoving volume from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist
    a = 1.0 / (1.0 + z)

    RH_Mpc = cosmo.RH_Mpc()

    ccl_cve = a**2 * pyccl.comoving_volume_element(ccl_cosmo, a)
    nc_cve = np.array(
        [dist.comoving_volume_element(cosmo, z_i) * RH_Mpc**3 for z_i in z]
    )

    return CompareFunc1d(
        x=z,
        y1=ccl_cve,
        y2=nc_cve,
        model=model,
        x_symbol="z",
        y_symbol=r"\mathrm{d}V/\mathrm{d}z\mathrm{d}\Omega",
        y_unit="Mpc${}^3$",
    )


# Growth factor


def compare_growth_factor(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare growth factor from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    ps_lin = cosmology.ps_ml
    a = 1.0 / (1.0 + z)

    gf = ps_lin.peek_gf()
    gf.prepare(cosmo)

    ccl_D_z = pyccl.growth_factor(ccl_cosmo, a)
    nc_D_z = np.array([gf.eval(cosmo, z_i) for z_i in z])

    return CompareFunc1d(
        x=z, y1=ccl_D_z, y2=nc_D_z, model=model, x_symbol="z", y_symbol=r"D"
    )


def compare_growth_rate(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    z: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare growth rate from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    ps_lin = cosmology.ps_ml
    a = 1.0 / (1.0 + z)

    gf = ps_lin.peek_gf()
    gf.prepare(cosmo)

    ccl_f_z = pyccl.growth_rate(ccl_cosmo, a)
    nc_f_z = np.array(
        [-(1.0 + z_i) / gf.eval(cosmo, z_i) * gf.eval_deriv(cosmo, z_i) for z_i in z]
    )

    return CompareFunc1d(
        x=z, y1=ccl_f_z, y2=nc_f_z, model=model, x_symbol="z", y_symbol=r"f"
    )


# Scale factor as a function of comoving distance


def compare_scale_factor(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    chi: np.ndarray,
    *,
    model: str = "unnamed",
):
    """Compare scale factor from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist

    RH_Mpc = cosmo.RH_Mpc()
    ccl_a = pyccl.scale_factor_of_chi(ccl_cosmo, chi)
    nc_a = np.array([1.0 / (1.0 + dist.inv_comoving(cosmo, d / RH_Mpc)) for d in chi])

    return CompareFunc1d(
        x=chi,
        y1=ccl_a,
        y2=nc_a,
        model=model,
        x_symbol=r"\chi",
        y_symbol=r"a",
        x_unit="Mpc",
    )


# Weaklensing related


def compare_Sigma_crit(
    ccl_cosmo: pyccl.Cosmology,
    cosmology: ncc.Cosmology,
    zs: np.ndarray,
    zl: float,
    *,
    model: str = "unnamed",
):
    """Compare critical surface mass density from CCL and NumCosmo."""
    cosmo = cosmology.cosmo
    dist = cosmology.dist

    a_s = 1.0 / (1.0 + zs)
    a_l = 1.0 / (1.0 + zl)

    sigmaC_ccl = pyccl.sigma_critical(ccl_cosmo, a_lens=a_l, a_source=a_s)
    sigmaC_nc = np.array([dist.sigma_critical(cosmo, zsi, zl) for zsi in zs])

    return CompareFunc1d(
        x=zs,
        y1=sigmaC_ccl,
        y2=sigmaC_nc,
        model=model,
        x_symbol="z_s",
        y_symbol=r"\Sigma_\mathrm{crit}",
    )
