#!/usr/bin/env python

import gi
import math
from scipy.stats import norm
import numpy as np

gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

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

  def do_eval1 (self, mset, x):
    mid = mset.get_id_by_ns ("NcPySLineModel")
    slm = mset.peek (mid)

    res = slm.props.b + slm.props.m

    return [res]
  
#
# Register our new MSetFunc
#
GObject.type_register (PyTestFunc)
