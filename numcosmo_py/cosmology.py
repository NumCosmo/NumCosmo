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

from . import Nc


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
        compute_inv_comoving: bool = True,
    ) -> None:
        """Initialize the NumCosmo cosmology class."""
        self.cosmo = cosmo
        self.dist = dist
        self._ps_ml = ps_ml
        self._ps_mnl = ps_mnl
        self.dist.compute_inv_comoving(compute_inv_comoving)
        self.recomb = Nc.RecombSeager()

        self.prepare()

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

    def prepare(self) -> None:
        """Prepare the cosmology for calculations."""
        self.dist.prepare_if_needed(self.cosmo)
        self.recomb.prepare_if_needed(self.cosmo)
        if self._ps_ml is not None:
            self._ps_ml.prepare_if_needed(self.cosmo)
        if self._ps_mnl is not None:
            self._ps_mnl.prepare_if_needed(self.cosmo)
