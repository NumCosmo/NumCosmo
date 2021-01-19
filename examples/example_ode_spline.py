#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import ctypes
from math import *
from gi.repository import NumCosmoMath as Ncm
from gi.repository import NumCosmo as Nc
from gi.repository import GObject

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

class TestClass (Ncm.Model):
  def __call__ (self, *args):
    return args[0]


aas = TestClass ()

def test (y, x, data):
  return y

test.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_char_p]
test.restype = ctypes.c_double

s = Ncm.SplineCubicNotaknot.new ()

os = Ncm.OdeSpline.new (s, test)
os.set_reltol (1.0e-3)

os.props.xi = 0.0
os.props.xf = 5.0
os.props.yi = 1.0

nhaca = [1,2,3,4]

os.prepare (id (nhaca))

ss = os.peek_spline()

for i in range (ss.len):
    print ("%d % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g" % (i, ss.xv.get (i), ss.yv.get (i), ss.b.get (i), ss.c.get(i), ss.d.get(i)))

#for i in range (100):
#  x = 1.0 / 99.0 * i
#  expx = exp (x)
#  odex = ss.eval (x)
#  print (x, expx, odex, fabs ((expx - odex) / expx))
