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
p.props.abstol     = 1.0e-30
p.props.reltol     = 1.0e-9
p.props.nknots     = 1000
p.props.noboundary = False
p.set_property ("lambda", 0.0)

psi0 = Ncm.QMPropGauss.new (0.0, 1.0, 1.0, -1.0)
#psi0 = Ncm.QMPropExp.new (3.0, 2.0, -1.0)

#print psi0.eval (1.0)
tstep = 1.1e-6
xf    = 60.0
xfp   = 20.0
p.set_init_cond_gauss (psi0, 0.0, xf)
#p.set_init_cond_exp (psi0, 0.0, xf)

x             = np.linspace (0.0, xfp, 2000) #p.get_knots ()
#x          = p.get_knots ()
Re_psi_s      = p.get_Re_psi ()
Re_psi_s_eval = np.vectorize (Re_psi_s.eval)
Im_psi_s      = p.get_Im_psi ()
Im_psi_s_eval = np.vectorize (Im_psi_s.eval)

max_Re_psi    = max (np.abs (Re_psi_s_eval (x)))
max_Im_psi    = max (np.abs (Im_psi_s_eval (x)))
yb            = max (max_Re_psi, max_Im_psi)

#psi   = np.array (p.get_psi ())

fig = plt.figure()

ax = plt.axes (xlim = (x[0], x[-1]), ylim = (-1.5e0 * yb, 1.5e0 * yb))
ax.grid ()
#ax = plt.axes (xlim = (x[0], x[-1]), ylim = (-700.0, 2.0 * max (psi)))
#ax = plt.axes (xlim = (x[0], x[-1]), ylim = (-xf, xf))

#ax.set_title("nhoca")
ttl = ax.text (.1, 1.005, '', transform = ax.transAxes)

#print psi

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
  p.evolve_spec (tf)
  #psi = np.array (p.get_psi ())
  rho_s      = p.get_rho ()
  dS_s       = p.get_dS ()
  Re_psi_s   = p.get_Re_psi ()
  Im_psi_s   = p.get_Im_psi ()

  Re_psi_s_eval = np.vectorize (Re_psi_s.eval)
  Im_psi_s_eval = np.vectorize (Im_psi_s.eval)
  rho_s_eval    = np.vectorize (rho_s.eval)
  dS_s_eval     = np.vectorize (dS_s.eval)
  
  #lines[0].set_data (x, psi[::2]) 
  #lines[1].set_data (x, psi[1::2]) 
  #lines[2].set_data (x, rho_s_eval (x)) 
  #lines[0].set_data (x, np.log (np.abs(rho_s_eval (x))) / 700.0)
  lines[0].set_data (x, np.sqrt (np.abs (rho_s_eval (x))))
  lines[1].set_data (x, Re_psi_s_eval (x))
  lines[2].set_data (x, Im_psi_s_eval (x))
  lines[3].set_data (x, dS_s_eval (x) / 10.0)
  #ax.set_title("nhoca %d" % i)
  ttl.set_text ("t = % .15f, norma = % .15f" % (tf, rho_s.eval_integ (0.0, xf)))

  return lines

anim = animation.FuncAnimation (fig, animate, np.arange (0, 200000), init_func = init, interval = 1, blit = True, repeat = False)

#mywriter = animation.FFMpegWriter()
#anim.save ('mymovie.mp4', writer = mywriter)

plt.show ()
#plt.savefig ("anim.mpeg")

'''
p.props.np = 1000
p.gauss_ini (0.0, 10.0, 1.0, 1.0)
p.propto (1.0, 0.05)

def psi(x,t):                                     
  r = p.propto (x, t)
  return r[0] + 1j * r[1]

psiv = np.vectorize (lambda x, t: abs (psi(x, t)))

fig, ax = plt.subplots()

x     = np.linspace (0, 60, 1000)
line, = plt.plot (x, psiv(x, 0.01))

def animate(i):
  line.set_ydata (psiv (x, 0.02 * (i + 1)))  # update the data
  return line,

# Init only required for blitting to give a clean slate.
def init():
    line.set_ydata (np.ma.array (x, mask = True))
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(0, 100), init_func=init, interval=5, blit=True)

plt.show ()
'''