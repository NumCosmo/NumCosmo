#!/usr/bin/env python

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import os.path

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
from py_sline_data import PySLineData
from py_sline_gauss import PySLineGauss

import matplotlib.pyplot as plt

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init()

#
# Creating a new NcmRNG object with seed 123 and default algorithm
#
rng = Ncm.RNG.seeded_new(None, 123)

#
# Instantiating a new SLine model object and setting
# some values for its parameters.
#
slm = PySLineModel()
slm.props.alpha = 1.0
slm.props.a = 0.5

#
# Instantiating a new empty SLine data object.
#
sld = None
if (len(sys.argv) != 2) or (sys.argv[1] != '--plain' and sys.argv[1] != '--gauss'):
    print("usage: example_create_data.py --plain or --gauss")
    sys.exit(-1)
elif sys.argv[1] == '--plain':
    sld = PySLineData(len=50)
else:
    sld = PySLineGauss(len=50)
    sld.xv.set_array(np.linspace(0.0, 10.0, sld.get_size()))
    sld.create_random_cov(slm, rng)

#
# New Model set object including slm with parameters
# set as free.
#
mset = Ncm.MSet.empty_new()
mset.set(slm)
mset.param_set_all_ftype(Ncm.ParamType.FREE)
mset.prepare_fparam_map()

#
# Creating a new Serialization object, with a data
# file does not exists, generate a new sample using
# mset as fiducial model and save it to data_file.
# 
# If data_file already exists, reload sld from it.
#
ser = Ncm.Serialize.new(0)
data_file = "example_data.obj"
sld.resample(mset, rng)
ser.to_binfile(sld, data_file)

#
# Plotting the created data points
# 

diag = sld.y.dup()

sld.cov.get_diag(diag)

xa = sld.xv.dup_array()
plt.errorbar(sld.xv.dup_array(), sld.y.dup_array(), yerr=diag.dup_array(), fmt='o', label=r'$y\pm\sigma$')
plt.plot (xa, [slm.f_x(x) for x in xa], lw=1.5, label="y(x)")

plt.savefig("example_data_f.pdf")

plt.show()
