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
# Creating a new class implementing our object Ncm.Data
#
class PySLineData (Ncm.Data):

  #
  # Creating three properties, len (length), dof (degrees of freedom) and data (NcmMatrix containing 
  # the data).
  #
  len  = GObject.Property (type = GObject.TYPE_UINT, default = 0, flags = GObject.PARAM_READWRITE)
  dof  = GObject.Property (type = GObject.TYPE_UINT, default = 0, flags = GObject.PARAM_READWRITE)
  data = GObject.Property (type = Ncm.Matrix, flags = GObject.PARAM_READWRITE)

  #
  # The contructor assigns some default values and calls the father's constructor.
  # 
  def __init__ (self, len = 600):
    Ncm.Data.__init__ (self)
    self.len  = len
    self.dof  = self.len - 2
    self.data = Ncm.Matrix.new (self.len, 3)

  #
  # Implements the virtual method get_length.
  #
  def do_get_length (self):
    return self.len

  #
  # Implements the virtual method get_dof.
  #
  def do_get_dof (self):
    return self.dof
    
  #
  # Implements the virtual method `begin'.
  # This method usually do some groundwork in the data 
  # before the actual calculations. For example, if the likelihood
  # involves the decomposition of a constant matrix, it can be done
  # during `begin' once and then used afterwards.
  #
  def do_begin (self):
    return

  #
  # Implements the virtual method `prepare'.
  # This method should do all the necessary calculations using mset
  # to be able to calculate the likelihood afterwards.
  #
  def do_prepare (self, mset):
    self.dof  = self.len - mset.fparams_len ()
    return
  
  #
  # Implements the virtual method `resample'.
  # This method creates a new sample from the fiducial model.
  # It is necessary to the MC analysis but can be skipped if
  # doing only MCMC.
  #
  def do_resample (self, mset, rng):
    mid = mset.get_id_by_ns ("NcPySLineModel")
    slm = mset.peek (mid)
        
    my_data = []

    for i in range (self.len):
      x = rng.uniform_gen (0.0, 10.0)
      s = x * rng.uniform_gen (0.4, 0.5)
      v = rng.gaussian_gen (slm.f_x (x), s)

      my_data.append ([x, v, s])

    my_array = np.array (my_data).reshape (3 * self.len)
    self.data = Ncm.Matrix.new_array (my_array, 3)
    self.set_init (True)

  #
  # Implements the virtual method `m2lnL'.
  # This method should calculate the value of the likelihood for
  # the model set `mset'.
  # 
  def do_m2lnL_val (self, mset):
    mid = mset.get_id_by_ns ("NcPySLineModel")
    slm = mset.peek (mid)
    
    m2lnL = 0.0
    for i in range (self.len):
      x = self.data.get (i, 0)
      v = self.data.get (i, 1)
      s = self.data.get (i, 2)
      m2lnL += ((slm.f_x (x) - v) / s)**2
    return m2lnL
  
#
# Register our new Python class PySLineData
#
GObject.type_register (PySLineData)
