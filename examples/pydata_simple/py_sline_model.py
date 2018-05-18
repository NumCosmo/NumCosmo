#!/usr/bin/env python

import gi
import math

gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
# New ModelBuilder object, defines a new model NcHIPrimExample implementing
# the base Ncm.Model abstract class.
# 
mb = Ncm.ModelBuilder.new (Ncm.Model, "NcPySLineModel", "A simple python data example model")

#
# New parameter m to describe the slope
#
mb.add_sparam ("m", "m", 0.0, 5.0, 0.1, 0.0, 2.0, Ncm.ParamType.FREE)

#
# New parameter b to describe the intercept
#
mb.add_sparam ("b", "b", -10.0, 10.0, 0.1, 0.0, 1.0, Ncm.ParamType.FREE)

#
# Creates a new GObject, it is not a Python object yet!
#
GNcPySLineModel = mb.create ()

#
# Workaraound to make the pygobject "see" the new object
#
GObject.new (GNcPySLineModel)

#
# Gets the Python version of the object (.pytype) and register
# it as a PyGObject object.
#
NcPySLineModel = GNcPySLineModel.pytype
GObject.type_register (NcPySLineModel)

#
# Creating a new class implementing our object NcPySLineModel
#
class PySLineModel (NcPySLineModel):
  #
  # Defining some property which is not part of the model paramers.
  # All model parameters must be defined by the ModelBuilder.
  #
  some_property = GObject.Property (type = str)

  #
  # Calling the father's constructor
  #  
  def __init__ (self):
    NcPySLineModel.__init__ (self)
  
#
# Register our new Python class PyNcPySLineModel
#
GObject.type_register (PySLineModel)
