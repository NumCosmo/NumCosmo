#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import scipy.stats as ss
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time
import math
import sys
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
#  New homogeneous and isotropic cosmological model NcHICosmoQGRW
#
cosmo = Nc.HICosmo.new_from_name(Nc.HICosmo, "NcHICosmoSFB")
csq1d = Nc.HIPertAdiab.new ()

ki    = 1.0e0
kf    = 1.0e9
k_a   = np.geomspace (ki, kf, 4)

csq1d.set_k (kf)
(Found2, etafa) = csq1d.find_adiab_time_limit (cosmo, +1.0e-25, +1.0e15, 1.0e1)

print (Found2, etafa)

