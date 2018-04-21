#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import math
from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init ()

#
# New ModelBuilder object, defines a new model NcHIPrimExample implementing
# the Nc.HIPrim abstract class.
# 
mb = Ncm.ModelBuilder.new (Nc.HIPrim, "NcHIPrimExample", "A example primordial model")

#
# New parameter A_s to describe the spectrum amplitud (it is usually better
# to work with ln(10^10As))
#
mb.add_sparam ("A_s", "As", 0.0, 1.0, 0.1, 0.0, 1.0e-9, Ncm.ParamType.FREE)

#
# New parameter n_s to describe the spectral index
#
mb.add_sparam ("n_s", "ns", 0.5, 1.5, 0.1, 0.0, 0.96,   Ncm.ParamType.FREE)

#
# New parameters a, b and c to describe the spectrum modification
#
mb.add_sparam ("a", "a", 0.0,   1.0, 0.01, 0.0,   0.5, Ncm.ParamType.FREE)
mb.add_sparam ("b", "b", 0.0, 1.0e4, 0.10, 0.0, 100.0, Ncm.ParamType.FREE)
mb.add_sparam ("c", "c", 0.0,   6.0, 0.10, 0.0,   0.0, Ncm.ParamType.FREE)

#
# Creates a new GObject, it is not a Python object yet!
#
GNcHIPrimExample = mb.create ()

#
# Workaraound to make the pygobject "see" the new object
#
GObject.new (GNcHIPrimExample)

#
# Gets the Python version of the object (.pytype) and register
# it as a PyGObject object.
#
NcHIPrimExample = GNcHIPrimExample.pytype
GObject.type_register (NcHIPrimExample)

#
# Creating a new class implementing our object NcHIPrimExample
#
class PyHIPrimExample (NcHIPrimExample):
  #
  # Defining some property which is not part of the model paramers.
  # All model parameters must be defined by the ModelBuilder.
  #
  some_property = GObject.Property (type = str)
  
  #
  # Implementing the adiabatic spectrum function
  #
  def do_lnSA_powspec_lnk (self, lnk):
    lnk0 = self.get_lnk_pivot ()
    lnka = lnk - lnk0
    ka = math.exp (lnka)
    As = self.props.As
    ns = self.props.ns
    a  = self.props.a
    b  = self.props.b
    c  = self.props.c
    a2 = a * a
    
    return (ns - 1.0) * lnka + math.log (As) + math.log1p (a2 * math.cos (b * ka + c)**2) 

#
# Register our new Python class PyHIPrimExample
#
GObject.type_register(PyHIPrimExample)
