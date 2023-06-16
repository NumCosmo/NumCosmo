"""Module for NumCosmo Python bindings."""

# The hack below is necessary to make the NumCosmo Python bindings work.
# This allows the use of our stubs and it also makes pylint and mypy happy.

import sys
import gi

gi.require_version("NumCosmo", "1.0")
gi.require_version("NumCosmoMath", "1.0")

# pylint:disable=wrong-import-position,unused-import,wildcard-import,unused-wildcard-import
from gi.repository import NumCosmo
from gi.repository.NumCosmo import *  # type: ignore

sys.modules[__name__] = NumCosmo
