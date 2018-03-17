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
    return t * t

  def do_eval_mnu (self, model, t, k):
    return t * t * k

  def do_eval_dlnmnu (self, model, t, k):
    return 2.0 / t

  def do_eval_system (self, model, t, k):
    return k, 2.0 / t, 0.0

  def do_nsing (self, model, k):
    return 1

  def do_get_sing_info (self, model, k, sing):
    return 0.0, -1.0, +1.0, Ncm.HOAASingType.ZERO    

  def do_eval_sing_mnu (self, model, t, k, sing):
    return t * t * k

  def do_eval_sing_dlnmnu (self, model, t, k, sing):
    return 2.0 / t

  def do_eval_sing_system (self, model, t, k, sing):
    return k, 2.0 / t, 0.0


def sol_q (k, t):
  a = k * t
  return math.sin (a) / a

def sol_p (k, t):
  a = k * t
  return t * (math.cos (a) - math.sin (a) / a)

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

print ("# ", t0, t1)
print ("# ", Aq, Av)

ta = np.linspace (-5.0e-1, +5.0e-1, 1000000)

for t in ta:
  (q, v, Pq, Pv) = hoaa.eval_QV (None, t)
  (upsilon, gamma, qbar, pbar) = hoaa.eval_AA (None, t)
  
  S  = Aq * q + Av * v
  PS = Aq * Pq + Av * Pv
  
  mnu   = hoaa.eval_mnu (None, t, k)
  nu    = hoaa.eval_nu (None, t, k)
  lnmnu = math.log (mnu)
  
  I = 0.5 * (mnu * q**2 + Pq**2 / mnu)
  J = 0.5 * (mnu * v**2 + Pv**2 / mnu)

  print (t, S, sol_q (k, t), PS, sol_p (k, t), I, J, math.sqrt (I * J), upsilon, gamma, qbar, pbar, (qbar**2 + pbar**2) / math.hypot (upsilon, 1.0 / math.cosh (lnmnu)) - 1.0, Aq * q / (Av * v) + 1.0)

