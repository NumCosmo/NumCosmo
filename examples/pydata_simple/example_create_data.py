#!/usr/bin/python2

import sys
import gi
import math
import numpy as np
import matplotlib.pyplot as plt
import os.path

gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from py_sline_model import PySLineModel
from py_sline_data import PySLineData
from py_sline_gauss import PySLineGauss

#
# Instantiating a new empty SLine data object.
#
sld = None
if (len (sys.argv) != 2) or (sys.argv[1] != '--plain' and sys.argv[1] != '--gauss'):
  print "usage: example_create_data.py --plain or --gauss"
  sys.exit (-1)
elif sys.argv[1] == '--plain':
  sld = PySLineData (len = 5000)
else: 
  sld = PySLineGauss (len = 5000)
  sld.xv.set_array (np.linspace (0.0, 10.0, sld.get_size ()))
  
#
# Instantiating a new SLine model object and setting
# some values for its parameters.
#
slm = PySLineModel ()
slm.props.m = 0.987
slm.props.b = 0.123

#
# New Model set object including slm with parameters
# set as free.
#
mset = Ncm.MSet.empty_new ()
mset.set (slm)
mset.param_set_all_ftype (Ncm.ParamType.FREE)
mset.prepare_fparam_map ()

#
# Creating a new Serialization object, with a data
# file does not exists, generate a new sample using
# mset as fiducial model and save it to data_file.
# 
# If data_file already exists, reload sld from it.
#
ser = Ncm.Serialize.new (0)
data_file = "example_data.obj"
rng = Ncm.RNG.new ()
rng.set_random_seed (False)
sld.resample (mset, rng)
ser.to_binfile (sld, data_file)

