#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import math
import numpy as np
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

gsf = Nc.GalaxySelfunc.new (2)

gsf.load_from_txts ("nc_galaxy_selfunc", None)

ser = Ncm.Serialize.new (0)

ser.to_file (gsf, "gsf_test.obj")

gsf2 = ser.from_file ("gsf_test.obj")

for i in range (2):
  zmin  = gsf.get_zmin (i)
  zmax  = gsf.get_zmax (i)
  zmin2 = gsf2.get_zmin (i)
  zmax2 = gsf2.get_zmax (i)
  
  print ("# SHELL %d (% 22.15g % 22.15g) (% 22.15g % 22.15g)" % (i, zmin, zmax, zmin2, zmax2))

  for z in np.linspace (zmin, zmax, 100):
    print ("% 22.15g % 22.15g % .5e" % (gsf.eval (i, z), gsf2.eval (i, z), gsf2.eval (i, z) / gsf.eval (i, z) - 1.0))

