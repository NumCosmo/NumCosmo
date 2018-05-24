#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import timeit
import sys
import time
import math
import numpy as np
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
sys.argv = Ncm.cfg_init_full (sys.argv)

rng = Ncm.RNG.new (None)
rng.set_random_seed (True)

mj  = Ncm.MPIJobTest.new ()
mj.set_rand_vector (12, rng)

ser = Ncm.Serialize.new (0)
mj.init_all_slaves (ser)

a = []
b = []
for t in np.arange (12.0):
  a.append (Ncm.Vector.new_array ([t]))
  b.append (Ncm.Vector.new (1))

#mj.run_array (a, b)
t1 = timeit.Timer('mj.run_array (a, b)', "from __main__ import mj, a, b")
print ("Timming: ", t1.timeit (1))

for a_i, b_i in zip (a, b):
  i = int (a_i.get (0))
  print ("%.3d % 22.15g % 22.15g %e" % (i, b_i.get (0), mj.props.vector.get (i), b_i.get (0) / mj.props.vector.get (i) - 1.0))
  
  

#time.sleep (600)


