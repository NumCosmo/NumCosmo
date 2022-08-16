#!/usr/bin/env python3

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

ki    = 1.0e5
kf    = 1.0e6
k_a   = np.geomspace (ki, kf, 10)

csq1d.set_k (kf)
(Found2, etafa) = csq1d.find_adiab_time_limit (cosmo, t0=20.0, t1=100.0, reltol=1.0e1)

print (f"{Found2}, {etafa}")

csq1d.set_k (ki)
(Found1, etaia) = csq1d.find_adiab_time_limit (cosmo, t0=-100.0, t1=20.0, reltol=1.0e1)

print (f"{Found1}, {etaia}")
fig = plt.figure (dpi = 120)

max_etaf = -200.0
min_etai = 200.0

for k in k_a:
  csq1d.set_k (k)
  csq1d.set_reltol (1.0e-5)

#  (Found1, etai)  = csq1d.find_adiab_time_limit (cosmo, -100.0, 20.0, 1.0e-2)
  #(Found2, etafa) = csq1d.find_adiab_time_limit (cosmo, -20.0, 100.0, 1.0e0)
  #etaf = etafa 
  #csq1d.set_ti (etai)
  #csq1d.set_tf (etaf)
  #csq1d.set_init_cond_adiab (cosmo, etai)
  #csq1d.prepare(cosmo)
 # min_etai = min (etai, min_etai)
#  max_etaf = max (etaf, max_etaf)
#Check if the derivatives are right with sympy and then check why is init_cond_adiab crashing for some momentum. Find a good interval that does not crack also the get time array.
#  eta_a, eta_s = csq1d.get_time_array ()
  print(csq1d.do_eval_powspec_factor(csq1d, cosmo))
  print(csq1d.do_eval_H(csq1d, cosmo, 2, 2))
  print(csq1d.do_eval_x(csq1d, cosmo, 2, 2))