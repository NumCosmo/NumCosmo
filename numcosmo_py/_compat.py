#
# _compat.py
#
# Sat Oct 3 10:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# _compat.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Deprecated NumCosmo 0.27 Python API, kept until 1.0.

Each shim emits a DeprecationWarning and maps the 0.27 call onto the current
API. The whole module is removed in 1.0.
"""

import warnings
from typing import cast

from . import nc as Nc
from . import ncm as Ncm

_model_add_submodel = Ncm.Model.add_submodel
_hicosmo_param_set_by_name = Nc.HICosmo.param_set_by_name
_reion_camb_set_z_from_tau = Nc.HIReionCamb.set_z_from_tau


def _add_submodel(self: Ncm.Model, submodel: Ncm.Model) -> None:
    """Attach submodel after construction (deprecated 0.27 API)."""
    warnings.warn(
        "Attaching a submodel after construction is deprecated and is dropped in "
        "NumCosmo 1.0; pass it to the constructor instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    _model_add_submodel(self, submodel)


def _param_set_by_name(self: Nc.HICosmo, name: str, val: float) -> None:
    """Set a parameter by name, accepting the 0.27 cosmology parameter Yp."""
    if name != "Yp":
        _hicosmo_param_set_by_name(self, name, val)
        return

    warnings.warn(
        "Yp is no longer a cosmology parameter and is dropped in NumCosmo 1.0; "
        "put an NcBBNParametrized in the cosmology's bbn slot and set Yp there.",
        DeprecationWarning,
        stacklevel=2,
    )
    bbn = self.peek_bbn()
    # In 0.27 a fixed Yp was ignored: Yp came from BBN, as with the default
    # NcBBNParthenope now.
    if isinstance(bbn, Nc.BBNParametrized):
        bbn.param_set_by_name("Yp", val)


def _set_z_from_tau(self: Nc.HIReionCamb, *args: Nc.HICosmo | float) -> None:
    """Set z_re from tau, accepting the 0.27 signature (cosmo, tau)."""
    if len(args) == 1:
        _reion_camb_set_z_from_tau(self, cast(float, args[0]))
        return
    if len(args) != 2:
        raise TypeError(f"set_z_from_tau takes 1 or 2 arguments, {len(args)} given.")

    cosmo = cast(Nc.HICosmo, args[0])
    tau = float(cast(float, args[1]))
    warnings.warn(
        "set_z_from_tau(cosmo, tau) is deprecated and is dropped in NumCosmo 1.0; "
        "attach the reionization model to the cosmology and call set_z_from_tau(tau).",
        DeprecationWarning,
        stacklevel=2,
    )
    host = self.peek_host()
    if host is not None and host is not cosmo:
        raise ValueError(
            "set_z_from_tau: cosmo is not the host cosmology of this reionization model."
        )
    if host is not None or isinstance(self.peek_reparam(), Nc.HIReionCambReparamTau):
        _reion_camb_set_z_from_tau(self, tau)
        return

    # Not attached yet, as 0.27 allowed: compute with the cosmology given.
    z_re = self.calc_z_from_tau(cosmo, tau)
    self.orig_param_set(Nc.HIReionCambSParams.HII_HEII_Z, z_re)


Ncm.Model.add_submodel = _add_submodel  # type: ignore[method-assign]
Nc.HICosmo.param_set_by_name = _param_set_by_name  # type: ignore[method-assign,assignment]
Nc.HIReionCamb.set_z_from_tau = _set_z_from_tau  # type: ignore[method-assign,assignment]
