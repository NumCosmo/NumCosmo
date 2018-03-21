#!/usr/bin/python2

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import math
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

def _blit_draw(self, artists, bg_cache):
    # Handles blitted drawing, which renders only the artists given instead
    # of the entire figure.
    updated_ax = []
    for a in artists:
        # If we haven't cached the background for this axes object, do
        # so now. This might not always be reliable, but it's an attempt
        # to automate the process.
        if a.axes not in bg_cache:
            # bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.bbox)
            # change here
            bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.figure.bbox)
        a.axes.draw_artist(a)
        updated_ax.append(a.axes)

    # After rendering all the needed artists, blit each axes individually.
    for ax in set(updated_ax):
        # and here
        # ax.figure.canvas.blit(ax.bbox)
        ax.figure.canvas.blit(ax.figure.bbox)

# MONKEY PATCH!!
matplotlib.animation.Animation._blit_draw = _blit_draw

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

p = Ncm.QMProp.new ()
p.props.abstol     = 1.0e-230
p.props.reltol     = 1.0e-7
p.props.nknots     = 101
p.props.noboundary = False
p.set_property ("lambda", 0.0)

psi0 = Ncm.QMPropGauss.new (0.0, 1.0, 1.0, -1.0)
#psi0 = Ncm.QMPropExp.new (3.0, 2.0, -1.0)

#print (psi0.eval (1.0))
tstep = 1.0e-3
xf    = 10.0
xfp   = 10.0
p.set_init_cond_gauss (psi0, 1.0e-2, xf)
#p.set_init_cond_exp (psi0, 0.0, xf)

x          = np.linspace (0.0, xfp, 10000) #p.get_knots ()
psi        = p.eval_psi (x)
max_Re_psi = max (np.abs (psi[0::2]))
max_Im_psi = max (np.abs (psi[1::2]))
yb         = max (max_Re_psi, max_Im_psi)

fig = plt.figure()

ax = plt.axes (xlim = (x[0], x[-1]), ylim = (-1.5e0 * yb, 1.5e0 * yb))
#ax = plt.axes (xlim = (x[0], x[-1]), ylim = (-0.1, 0.1))
ax.grid ()

ttl = ax.text (.1, 1.005, '', transform = ax.transAxes)

N = 4
lines = [plt.plot([], [])[0] for _ in range(N)]

lines.append (ttl)

def init():    
  ttl.set_text('')
  lines[0].set_data ([], [])
  lines[1].set_data ([], [])
  lines[2].set_data ([], [])
  lines[3].set_data ([], [])
  return lines

def animate(i):
  tf = tstep * i
  #p.evolve_spec (tf)
  p.evolve (tf)

  psi   = np.array (p.eval_psi (x))
  rho   = np.array (p.eval_rho (x))
  dS    = np.array (p.eval_dS (x))
  rho_s = p.peek_rho_s ()
  
  lines[0].set_data (x, np.sqrt (np.abs (rho)))
  lines[1].set_data (x, psi[0::2])
  lines[2].set_data (x, psi[1::2])
  lines[3].set_data (x, dS / 10.0)

  ttl.set_text ("t = % .15f, norma = % .15f" % (tf, rho_s.eval_integ (0.0, xf)))

  return lines

anim = animation.FuncAnimation (fig, animate, np.arange (0, 200000), init_func = init, interval = 1, blit = True, repeat = False)

#mywriter = animation.FFMpegWriter()
#anim.save ('mymovie.mp4', writer = mywriter)

plt.show ()

