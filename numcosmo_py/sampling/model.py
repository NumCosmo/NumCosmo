#
# model.py
#
# Min Apr 3 11:30:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# model.py
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

"""A simplified interface for the NumCosmo models."""

import sys

from numcosmo_py import Ncm, GObject

GENERIC_MODEL_NAME = "NcmModelGeneric"

MODEL_BUILDER = Ncm.ModelBuilder.new(
    Ncm.Model, GENERIC_MODEL_NAME, "NumCosmo generic model"
)
MODEL_BUILDER.add_vparam(
    1,
    r"\theta",
    r"theta",
    -sys.float_info.max,
    +sys.float_info.max,
    1.0,
    0.0,
    0.0,
    Ncm.ParamType.FREE,
)

NcmTypeModelGeneric = MODEL_BUILDER.create()  # pylint:disable=invalid-name
GObject.new(NcmTypeModelGeneric)
NcmModelGeneric = NcmTypeModelGeneric.pytype  # pylint:disable=invalid-name
GObject.type_register(NcmModelGeneric)
del MODEL_BUILDER


def build_mset(ndim: int) -> Ncm.MSet:
    """Build a new MSet object."""

    if ndim < 1:
        raise ValueError("ndim must be >= 1")

    model = NcmModelGeneric(theta_length=ndim)

    mset = Ncm.MSet()
    mset.set(model)
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    return mset


def get_generic_model(mset: Ncm.MSet) -> Ncm.Model:
    """Get the generic model from an MSet object."""

    mid = mset.get_id_by_ns(GENERIC_MODEL_NAME)
    model = mset.peek(mid)

    if model is None:
        raise ValueError("No generic model found in MSet")

    return model
