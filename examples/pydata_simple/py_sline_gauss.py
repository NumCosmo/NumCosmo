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
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init ()

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
    self.dof = self.np
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
    self.dof = self.np - mset.fparams_len ()
    return

  #
  # Implements the virtual method `mean_func'.
  # This method should compute the theoretical mean for the gaussian
  # distribution.
  #
  def do_mean_func (self, mset, vp):
    mid = mset.get_id_by_ns ("NcPySLineModel")
    slm = mset.peek (mid)

    #print self.np, slm.props.b, slm.props.m, vp.len ()
    
    for i in range (self.np):
      x = self.xv.get (i)
      #print slm.props.m, slm.props.b, x, slm.props.m * x + slm.props.b
      vp.set (i, slm.props.m * x + slm.props.b)

    return

  #
  # Implements the virtual method `cov_func'.
  # This method should compute the covariance matrix for the gaussian
  # distribution.
  # 
  # It has an internal flag to assert if the covariance matrix was 
  # already initialized.
  #
  def do_cov_func (self, mset, cov):
    if not self.cov_init:
      f = 0.9
      d = 1.0e-1
      
      cov.set_zero ()
    
      for i in range (self.np):
        x       = self.xv.get (i)
        sigma_i = math.fabs (x * f) + d
        
        cov.set (i, i, + sigma_i)
        
        if i + 1 < self.np:
          xp1       = self.xv.get (i + 1)
          sigma_ip1 = math.fabs (xp1 * f) + d
          
          cov.set (i,     i + 1, + i * math.sqrt (sigma_i * sigma_ip1) / self.np)
          cov.set (i + 1, i,     + i * math.sqrt (sigma_i * sigma_ip1) / self.np)
      
      #
      # sets cov = cov * cov
      # to guarantee that it is posdef 
      #
      tmp1 = cov.dup ()
      cov.dsymm('U', 1.0, tmp1, tmp1, 0.0);

      # Remove to print the covariance matrix
      # cov.log_vals ("# COV: ", "% 15.8g")
      
      self.cov_init = True      
      
      return True
    else:
      return False

#
# Register our new Python class PySLineGauss
#
GObject.type_register (PySLineGauss)
