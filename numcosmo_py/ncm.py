"""Module for NumCosmoMath Python bindings."""

# The hack below is necessary to make the NumCosmoMath Python bindings work.
# This allows the use of our stubs and it also makes pylint and mypy happy.

import sys
import gi

gi.require_version("NumCosmoMath", "1.0")

# pylint:disable=wrong-import-position,unused-import,wildcard-import,unused-wildcard-import
from gi.repository import NumCosmoMath
from gi.repository.NumCosmoMath import *  # type: ignore

sys.modules[__name__] = NumCosmoMath
