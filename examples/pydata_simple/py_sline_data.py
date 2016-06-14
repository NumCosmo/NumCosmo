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
# Creating a new class implementing our object NcPySLineData
#
class PySLineData (Ncm.Data):
  len  = GObject.Property (type = GObject.TYPE_UINT, default = 0, flags = GObject.PARAM_READWRITE)
  dof  = GObject.Property (type = GObject.TYPE_UINT, default = 0, flags = GObject.PARAM_READWRITE)
  data = GObject.Property (type = Ncm.Matrix, flags = GObject.PARAM_READWRITE)

  def __init__ (self):
    Ncm.Data.__init__ (self)
    self.len  = 600
    self.dof  = self.len - 2
    self.data = Ncm.Matrix.new (self.len, 3)

  def do_get_length (self):
    return self.len

  def do_get_dof (self):
    return self.dof
    
  def do_begin (self):
    return

  def do_prepare (self, mset):
    self.dof  = self.len - mset.fparams_len ()
    return
  
  def do_resample (self, mset, rng):
    mid = mset.get_id_by_ns ("NcPySLineModel")
    slm = mset.peek (mid)
        
    my_data = []

    for i in range (self.len):
      x = rng.uniform_gen (0.0, 10.0)
      s = x * rng.uniform_gen (0.4, 0.5)
      v = rng.gaussian_gen (slm.props.m * x + slm.props.b, s)

      my_data.append ([x, v, s])

    my_array = np.array (my_data).reshape (3 * self.len)
    self.data = Ncm.Matrix.new_array (my_array, 3)
    self.set_init (True)

  def do_m2lnL_val (self, mset):
    mid = mset.get_id_by_ns ("NcPySLineModel")
    slm = mset.peek (mid)
    
    m2lnL = 0.0
    for i in range (self.len):
      m2lnL += ((slm.props.m * self.data.get (i, 0) + slm.props.b - self.data.get (i, 1)) / self.data.get (i, 2))**2
    return m2lnL
  
#
# Register our new Python class PySLineData
#
GObject.type_register (PySLineData)
