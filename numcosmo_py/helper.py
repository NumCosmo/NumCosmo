"""Helper functions for numcosmo_py."""

from typing import Sequence, Type, TypeVar, cast
import numpy as np
import numpy.typing as npt


from . import Ncm, GObject

# Type variable for generic object duplication
T = TypeVar("T")


def npa_to_seq(npa: npt.NDArray[np.float64]) -> Sequence[float]:
    """Convert a NumPy array to a Python sequence."""
    return np.ravel(npa).tolist()


def duplicate_via_serialization(obj: T, ser: Ncm.Serialize | None = None) -> T:
    """Duplicate an object via serialization with type preservation.

    This utility function wraps Ncm.Serialize.dup_obj() with proper type hints,
    allowing the type checker to understand that the returned object has the
    same type as the input. This is particularly useful when duplicating
    NumCosmo objects for testing or when creating independent copies.

    Parameters
    ----------
    obj : T
        The object to duplicate via serialization
    ser : Ncm.Serialize | None, optional
        Serializer to use. If None, creates a new one with CLEAN_DUP option.
        Default is None.

    Returns
    -------
    T
        A duplicate of the input object with the same type

    Examples
    --------
    >>> from numcosmo_py import Nc
    >>> from numcosmo_py.helper import duplicate_via_serialization
    >>> kernel = Nc.XcorKernelCMBLensing(...)
    >>> kernel_dup = duplicate_via_serialization(kernel)
    >>> # No need for cast, type is preserved!
    >>> assert type(kernel_dup) == type(kernel)
    >>> assert kernel_dup is not kernel
    """
    if ser is None:
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    # The cast tells the type checker that dup_obj returns the same type
    # even though the GObject introspection doesn't provide this info
    assert isinstance(obj, GObject.GObject)
    return cast(T, ser.dup_obj(obj))


def register_model_class(mb: Ncm.ModelBuilder) -> Type:
    """Register a model class."""
    NcmTypeModelGeneric = mb.create()  # pylint:disable=invalid-name
    # We need to create a new instance to register the type
    GObject.new(NcmTypeModelGeneric)
    NcmModelGeneric = NcmTypeModelGeneric.pytype  # pylint:disable=invalid-name
    GObject.type_register(NcmModelGeneric)

    return NcmModelGeneric
