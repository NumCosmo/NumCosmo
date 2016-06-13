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
  def __init__ (self):
    Ncm.Data.__init__ (self)
    self.len  = 2000
    self.dof  = self.len - 2
    self.data = []

  def do_get_length (self):
    return self.len

  def do_get_dof (self):
    return self.dof
    
  def do_begin (self):
    return

  def do_prepare (self, mset):
    self.dof  = self.len - mset.fparams_len ()
    return
  
  def do_m2lnL_val (self, mset):
    mid = mset.get_id_by_ns ("NcPySLineModel")
    slm = mset.peek (mid)
    
    m2lnL = 0.0
    for d in self.data:
      m2lnL += ((slm.props.m * d[0] + slm.props.b - d[1]) / d[2])**2
    return m2lnL
    
  def simulate_data (self, sline):
    assert (GObject.type_is_a (sline, PySLineModel))

    self.data = []
    
    for i in range (self.len):
      x = np.random.uniform (0.0, 10.0)
      s = x * np.random.uniform (0.4, 0.5)
      v = np.random.normal (sline.props.m * x + sline.props.b, s)

      self.data.append ([x, v, s])
    self.set_init (True)
  
#
# Register our new Python class PySLineData
#
GObject.type_register (PySLineData)
