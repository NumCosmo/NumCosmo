#
# cosmology.py
#
# Tue Aug 6 16:11:50 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# cosmology.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""NumCosmo cosmology class."""

from typing import TypedDict
from enum import Enum
from . import Ncm, Nc


class ParameterDesc(TypedDict, total=False):
    """Parameter description."""

    name: str
    symbol: str
    scale: float
    lower_bound: float
    upper_bound: float
    abstol: float
    fit: bool
    value: float


class HIPrimModel(str, Enum):
    """Planck 18 primordial model."""

    ATAN = "atan"
    BPL = "broken-power-law"
    EXPC = "exponential-c"
    POWER_LAW = "power-law"
    SBPL = "smooth-broken-power-law"
    TWO_FLUIDS = "two-fluids"


class Cosmology:
    """NumCosmo cosmology class.

    Class to handle cosmological calculations using NumCosmo. The class is a wrapper
    around the NumCosmo cosmology, distance, and power-spectrum classes.
    """

    def __init__(
        self,
        *,
        cosmo: Nc.HICosmo,
        dist: Nc.Distance,
        ps_ml: Nc.PowspecML | None = None,
        ps_mnl: Nc.PowspecMNL | None = None,
        psf_tophat: Ncm.PowspecFilter | None = None,
        compute_inv_comoving: bool = True,
    ) -> None:
        """Initialize the NumCosmo cosmology class."""
        self.cosmo = cosmo
        self.dist = dist
        self._ps_ml = ps_ml
        self._ps_mnl = ps_mnl
        self._psf_tophat = psf_tophat
        self.dist.compute_inv_comoving(compute_inv_comoving)
        self.recomb = Nc.RecombSeager()

        self._mset = Ncm.MSet.new_array([cosmo])

        self.prepare()

    # Here we set a constructor that will be used to create a default cosmology
    @classmethod
    def default(
        cls,
        dist_max_z: float = 10.0,
        halofit_max_z: float = 5.0,
        halofit_reltol: float = 1.0e-7,
    ) -> "Cosmology":
        """Create a default cosmology."""
        cosmo = Nc.HICosmoDEXcdm()
        cosmo.omega_x2omega_k()
        cosmo["Omegak"] = 0.0
        cosmo.add_submodel(Nc.HIPrimPowerLaw.new())
        cosmo.add_submodel(Nc.HIReionCamb.new())
        dist = Nc.Distance.new(dist_max_z)
        ps_ml = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH())
        ps_mnl = Nc.PowspecMNLHaloFit.new(ps_ml, halofit_max_z, halofit_reltol)
        psf = Ncm.PowspecFilter.new(ps_ml, Ncm.PowspecFilterType.TOPHAT)
        return cls(cosmo=cosmo, dist=dist, ps_ml=ps_ml, ps_mnl=ps_mnl, psf_tophat=psf)

    @classmethod
    def default_minimal(cls, dist_max_z: float = 10.0) -> "Cosmology":
        """Create a minimal default cosmology."""
        cosmo = Nc.HICosmoDEXcdm()
        cosmo.omega_x2omega_k()
        cosmo["Omegak"] = 0.0
        cosmo.add_submodel(Nc.HIPrimPowerLaw.new())
        cosmo.add_submodel(Nc.HIReionCamb.new())
        dist = Nc.Distance.new(dist_max_z)
        return cls(cosmo=cosmo, dist=dist)

    @property
    def ps_ml(self) -> Nc.PowspecML:
        """Return the linear matter power spectrum."""
        if self._ps_ml is None:
            raise AttributeError("Linear matter power spectrum not set.")
        return self._ps_ml

    @property
    def ps_mnl(self) -> Nc.PowspecMNL:
        """Return the non-linear matter power spectrum."""
        if self._ps_mnl is None:
            raise AttributeError("Non-linear matter power spectrum not set.")
        return self._ps_mnl

    @property
    def psf_tophat(self) -> Ncm.PowspecFilter:
        """Return the top-hat power spectrum filter."""
        if self._psf_tophat is None:
            raise AttributeError("Top-hat power spectrum filter not set.")
        return self._psf_tophat

    @property
    def mset(self) -> Ncm.MSet:
        """Return the NumCosmo model set."""
        return self._mset

    def prepare(self) -> None:
        """Prepare the cosmology for calculations."""
        self.dist.prepare_if_needed(self.cosmo)
        self.recomb.prepare_if_needed(self.cosmo)
        if self._ps_ml is not None:
            self._ps_ml.prepare_if_needed(self.cosmo)
        if self._ps_mnl is not None:
            self._ps_mnl.prepare_if_needed(self.cosmo)
        if self._psf_tophat is not None:
            self._psf_tophat.prepare_if_needed(self.cosmo)


def create_cosmo(
    massive_nu: bool = False,
    prim_model: HIPrimModel = HIPrimModel.POWER_LAW,
) -> Nc.HICosmo:
    """Create a cosmology for CMB experiments."""
    if massive_nu:
        cosmo = Nc.HICosmoDEXcdm(massnu_length=1)
    else:
        cosmo = Nc.HICosmoDEXcdm()

    cosmo.params_set_default_ftype()
    cosmo.cmb_params()
    cosmo["H0"] = 70.0
    cosmo["omegab"] = 0.022
    cosmo["omegac"] = 0.12

    if massive_nu:
        cosmo["ENnu"] = 2.0328
        cosmo["massnu_0"] = 0.06
        cosmo.param_set_desc("massnu_0", {"fit": True})

    cosmo.param_set_desc("H0", {"fit": True})
    cosmo.param_set_desc("omegac", {"fit": True})
    cosmo.param_set_desc("omegab", {"fit": True})
    cosmo.param_set_desc("Omegak", {"fit": False})
    cosmo.param_set_desc("w", {"fit": False})

    match prim_model:
        case HIPrimModel.ATAN:
            prim = Nc.HIPrimAtan.new()
        case HIPrimModel.BPL:
            prim = Nc.HIPrimBPL.new()
        case HIPrimModel.EXPC:
            prim = Nc.HIPrimExpc.new()
        case HIPrimModel.POWER_LAW:
            prim = Nc.HIPrimPowerLaw.new()
            prim.param_set_desc("ln10e10ASA", {"fit": True})
            prim.param_set_desc("n_SA", {"fit": True})
        case HIPrimModel.SBPL:
            prim = Nc.HIPrimSBPL.new()
        case HIPrimModel.TWO_FLUIDS:
            prim = Nc.HIPrimTwoFluids(use_default_calib=True)
            prim.param_set_desc("ln10e10ASA", {"fit": True})
            prim.param_set_desc("lnk0", {"fit": True})
            prim.param_set_desc("lnw", {"fit": True})
        case _:
            raise ValueError(f"Invalid primordial model: {prim_model}")

    reion = Nc.HIReionCamb.new()
    reion.param_set_desc("z_re", {"fit": True})

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    return cosmo
