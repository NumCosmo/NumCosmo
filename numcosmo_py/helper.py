"""Helper functions for numcosmo_py."""

from typing import Sequence, Type
import numpy as np
import numpy.typing as npt


from . import Ncm, GObject


def npa_to_seq(npa: npt.NDArray[np.float64]) -> Sequence[float]:
    """Convert a NumPy array to a Python sequence."""
    return np.ravel(npa).tolist()


def register_model_class(mb: Ncm.ModelBuilder) -> Type:
    """Register a model class."""
    NcmTypeModelGeneric = mb.create()  # pylint:disable=invalid-name
    # We need to create a new instance to register the type
    GObject.new(NcmTypeModelGeneric)
    NcmModelGeneric = NcmTypeModelGeneric.pytype  # pylint:disable=invalid-name
    GObject.type_register(NcmModelGeneric)

    return NcmModelGeneric
