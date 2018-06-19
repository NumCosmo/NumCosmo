#!/usr/bin/env python

import math
from scipy.stats import norm
import numpy as np

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from py_sline_model import PySLineModel

#
# Creating a new MSetFunc
#
class PyTestFunc (Ncm.MSetFunc1):
  def __init__ (self):
    Ncm.MSetFunc.__init__ (self, dimension  = 1, nvariables = 0)
    self.symbol = "m+b"
    self.name   = "m_plus_b"

  def do_eval1 (self, mset, x):
    mid = mset.get_id_by_ns ("NcPySLineModel")
    slm = mset.peek (mid)

    res = slm.f_x (1.0)

    return [res]
  
#
# Register our new MSetFunc
#
GObject.type_register (PyTestFunc)
