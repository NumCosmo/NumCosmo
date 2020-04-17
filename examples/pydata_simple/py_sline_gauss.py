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
class PySLineGauss (Ncm.DataGaussCov):
  #
  # We need one vector property to save the independent variables x
  #
  xv = GObject.Property (type = Ncm.Vector, flags = GObject.PARAM_READWRITE)

  #
  # The contructor assigns some default values and calls the father's constructor.
  # 
  def __init__ (self, len = 600):
    Ncm.DataGaussCov.__init__ (self, n_points = len)

    if self.np > 0:
      self.xv = Ncm.Vector.new (self.np)
    else:
      self.xv = None
    
    self.cov_init = False

    #
    # Initializing to sane values
    #
    if self.np > 0:
      self.cov.set_identity ()
      self.xv.set_zero ()

  #
  # Implements the virtual method get_length.
  #
  def do_get_length (self):
    return self.np

  #
  # Implements the virtual method get_dof.
  #
  def do_get_dof (self):
    return self.np
    
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
    return

  #
  # Implements the virtual method `mean_func'.
  # This method should compute the theoretical mean for the gaussian
  # distribution.
  #
  def do_mean_func (self, mset, vp):
    mid = mset.get_id_by_ns ("NcPySLineModel")
    slm = mset.peek (mid)
    
    for i in range (self.np):
      x = self.xv.get (i)
      vp.set (i, slm.f_x (x))

    return

  def create_random_cov (self, slm, rng):    
    ya = [slm.f_x (x) for x in self.xv.dup_array ()]
    yv = Ncm.Vector.new_array (ya)
			
    self.cov.fill_rand_cov2 (yv, 0.5, 2.0, 15.0, rng)
	
    # Remove comment to print the covariance matrix
    # self.cov.log_vals ("# COV: ", "% 10.3g")
#
# Register our new Python class PySLineGauss
#
GObject.type_register (PySLineGauss)
