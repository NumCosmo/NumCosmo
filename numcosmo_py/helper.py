"""Helper functions for numcosmo_py."""

from typing import Sequence
import numpy as np
import numpy.typing as npt


def npa_to_seq(npa: npt.NDArray[np.float64]) -> Sequence[float]:
    """Convert a NumPy array to a Python sequence."""
    return np.ravel(npa).tolist()
