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

from . import Ncm, Nc


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
        psf: Ncm.PowspecFilter | None = None,
        compute_inv_comoving: bool = True,
    ) -> None:
        """Initialize the NumCosmo cosmology class."""
        self.cosmo = cosmo
        self.dist = dist
        self._ps_ml = ps_ml
        self._ps_mnl = ps_mnl
        self._psf = psf
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
        return cls(cosmo=cosmo, dist=dist, ps_ml=ps_ml, ps_mnl=ps_mnl, psf=psf)

    @property
    def ps_ml(self) -> Nc.PowspecML:
        """Return the linear matter power spectrum."""
        if self._ps_ml is None:
            raise ValueError("Linear matter power spectrum not set.")
        return self._ps_ml

    @property
    def ps_mnl(self) -> Nc.PowspecMNL:
        """Return the non-linear matter power spectrum."""
        if self._ps_mnl is None:
            raise ValueError("Non-linear matter power spectrum not set.")
        return self._ps_mnl

    @property
    def psf(self) -> Ncm.PowspecFilter:
        """Return the power spectrum filter."""
        if self._psf is None:
            raise ValueError("Power spectrum filter not set.")
        return self._psf

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
        if self._psf is not None:
            self._psf.prepare_if_needed(self.cosmo)
