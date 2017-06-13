#!/usr/bin/python2

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

class PyHOAATest (Ncm.HOAA):
  def __init__ (self):
    Ncm.HOAA.__init__ (self, opt = Ncm.HOAAOpt.DLNMNU_ONLY)

  def do_eval_nu (self, model, t, k):
    return k

  def do_eval_m (self, model, t, k):
    return 1.0 / (t * t)

  def do_eval_mnu (self, model, t, k):
    return k / (t * t)

  def do_eval_dlnmnu (self, model, t, k):
    return -2.0 / t

  def do_eval_system (self, model, t, k):
    return k, -2.0 / t, 0.0

  def do_nsing (self, model, k):
    return 1

  def do_get_sing_info (self, model, k, sing):
    return 0.0, -1.0, +1.0, Ncm.HOAASingType.INF

  def do_eval_sing_mnu (self, model, t, k, sing):
    return k / (t * t)

  def do_eval_sing_dlnmnu (self, model, t, k, sing):
    return -2.0 / t

  def do_eval_sing_system (self, model, t, k, sing):
    return k, -2.0 / t, 0.0

def sol_q (k, t):
  a = k * t
  return (-t / k**2) * (math.cos (a) - math.sin (a) / a)

def sol_p (k, t):
  a = k * t
  return math.sin (a) / a

hoaa = PyHOAATest ()

k = 1.0
hoaa.set_ti (-1.0e10)
hoaa.set_tf (+1.0e10)
hoaa.set_k (k)
hoaa.set_reltol (1.0e-14)

ti  = - 10.0

S1  = sol_q (k, ti)
PS1 = sol_p (k, ti)

hoaa.prepare ()

(t0, t1) = hoaa.get_t0_t1 ()

(Aq, Av) = hoaa.eval_solution (None, ti, S1, PS1)

print "# ", t0, t1
print "# ", Aq, Av

ta = np.linspace (-42.0, -39.00, 100000)

for t in ta:
  (q, v, Pq, Pv) = hoaa.eval_QV (None, t)
  (epsilon, gamma, sin_thetab, cos_thetab) = hoaa.eval_AA (None, t)
  
  S  = Aq * q + Av * v
  PS = Aq * Pq + Av * Pv
  
  mnu = hoaa.eval_mnu (None, t, k)
  nu  = hoaa.eval_nu (None, t, k)
  
  I = 0.5 * (mnu * q**2 + Pq**2 / mnu)
  J = 0.5 * (mnu * v**2 + Pv**2 / mnu)

  print t, S, sol_q (k, t), PS, sol_p (k, t), I, J, math.sqrt (I * J), math.cosh (epsilon)

