#!/usr/bin/env python

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
l1                 = 500.0
H0                 = -1.0e-1
p.props.abstol     = 0.0
p.props.reltol     = 1.0e-11
p.props.nknots     = 151
p.props.noboundary = False
p.set_property ("lambda", l1)
s1 = int (p.props.nknots / 10)
n1 = ceil (p.props.nknots / s1)
offset = 5.0

print ("# ", n1, " ", s1)

psi0 = Ncm.QMPropGauss.new (offset + 10.0, 1.0, 1.0, H0)
#psi0 = Ncm.QMPropExp.new (3.0, 2.0, -1.0)

#print (psi0.eval (1.0))
sim   = True
tstep = 5.0e-3
tf    = 5.0
xf    = offset + 15.01
xfp   = offset + 20.0
xi    = offset
p.set_init_cond_gauss (psi0, offset + 5.0, xf)
#p.set_init_cond_exp (psi0, 0.0, xf)

rho_s      = p.peek_rho_s ()
q          = rho_s.get_xv ().dup_array ()
x          = np.linspace (offset - 1.0, xfp, 10000) #p.get_knots ()
psi        = p.eval_psi (x)
max_Re_psi = max (np.abs (psi[0::2]))
max_Im_psi = max (np.abs (psi[1::2]))
yb         = max (max_Re_psi, max_Im_psi)

#fig = plt.figure()

fig, [ax, ax2] = plt.subplots (1, 2, figsize = (16, 12))

#ax = plt.axes (xlim = (x[0], x[-1]), ylim = (-1.5e0 * yb, 1.5e0 * yb))
#ax = plt.axes (xlim = (x[0], x[-1]), ylim = (-6.1, 6.1))
ax.set_xlim (xi, x[-1])
ax.set_ylim (-0.6, 0.8)
ax.grid ()

ax2.set_xlim (-0.1, tf)
ax2.set_ylim (xi, xfp)

ax.set_xlabel (r'$a(t)$')
ax2.set_xlabel (r'$t$')
ax2.set_ylabel (r'$a(t)$')

N = 4
ttl = ax.text (.1, 1.005, '', transform = ax.transAxes)
lines = []
lines.append (ax.plot ([], [], label=r'$\sqrt{\psi^*\psi}$')[0])
lines.append (ax.plot ([], [], label=r'$\mathrm{Re}(\psi)$')[0])
lines.append (ax.plot ([], [], label=r'$\mathrm{Im}(\psi)$')[0])
lines.append (ax.plot ([], [], label=r'$\partial_aS$')[0])
ax.legend (loc='best')

fig.tight_layout()
#lines = [ax.plot([], [])[0] for _ in range(N)]
lines.append (ttl)

lines.append (ax2.plot ([], [], 'bo', label=r'Bohm')[0])
lines.append (ax2.plot ([], [], 'ro', label=r'$\langle a(t)\rangle$')[0])
ax2.legend (loc='best')

ta   = [0.0]
traj = []
x_t  = []
y_t  = []

for i in range (n1):
  lines.append (ax2.plot ([], [])[0])
  qi = [q[i * s1]]
  traj.append (qi)

lines.append (ax2.plot ([], [])[0])
traj.append ([p.eval_int_xrho ()])

def init():    
  lines[0].set_data ([], [])
  lines[1].set_data ([], [])
  lines[2].set_data ([], [])
  lines[3].set_data ([], [])
  lines[4].set_text (None)
  lines[5].set_data ([], [])
  lines[6].set_data ([], [])
  for i in range (n1):
    lines[i + 7].set_data ([], [])

  return lines

if not sim:
  for i in np.arange (0, ceil (tf / tstep)):
    tf = tstep * i
    p.evolve (tf)

    rho_s = p.peek_rho_s ()
    q = rho_s.get_xv ().dup_array ()
    x = np.linspace (q[0], q[-1], 10000)
  
    psi   = np.array (p.eval_psi (x))
    rho   = np.array (p.eval_rho (x))
    dS    = np.array (p.eval_dS (x))
  
    y = []
  
    qm = p.eval_int_xrho ()
  
    y.append (np.sqrt (rho))
    y.append (psi[0::2])
    y.append (psi[1::2])
    y.append (dS / xf)
    y.append (q[::s1])
    y.append (p.eval_int_rho ())
    y.append (qm)

    x_t.append (x)
    y_t.append (y)

    ta.append (tf)
    for i in range (n1):
      traj[i].append (q[i * s1])
    traj[n1].append (qm)

def animate(i):
  tf = tstep * i
  if sim:
    #p.evolve_spec (tf)
    p.evolve (tf)

    rho_s = p.peek_rho_s ()
    q = rho_s.get_xv ().dup_array ()
    x = np.linspace (q[0], q[-1], 10000)
  
    psi   = np.array (p.eval_psi (x))
    rho   = np.array (p.eval_rho (x))
    dS    = np.array (p.eval_dS (x))
    
    lines[0].set_data (x, np.sqrt (rho))
    lines[1].set_data (x, psi[0::2])
    lines[2].set_data (x, psi[1::2])
    lines[3].set_data (x, dS / xf)
  
    qm = p.eval_int_xrho ()
  
    nq = len (q[::s1])
    tfa = [tf] * nq
    lines[5].set_data (tfa, q[::s1])
    lines[6].set_data ([tf], [qm])

    ta.append (tf)
    for i in range (n1):
      traj[i].append (q[i * s1])
      lines[i + 7].set_data (ta, traj[i])
    
    traj[n1].append (qm)
    lines[n1 + 7].set_data (ta, traj[n1])

    ttl.set_text ("t = % .15f, norma = % .15f" % (tf, p.eval_int_rho ()))
  
  else:
  
    lines[0].set_data (x_t[i], y_t[i][0])
    lines[1].set_data (x_t[i], y_t[i][1])
    lines[2].set_data (x_t[i], y_t[i][2])
    lines[3].set_data (x_t[i], y_t[i][3])
  
    lines[5].set_data (ta[i], y_t[i][4])
    lines[6].set_data (ta[i], y_t[i][6])
    
    for j in range (n1):
      lines[j + 7].set_data (ta[0:i], traj[j][0:i])
    
    lines[n1 + 7].set_data (ta[0:i], traj[n1][0:i])

    ttl.set_text ("t = % .15f, norma = % .15f" % (ta[i], y_t[i][5]))

  return lines

anim = animation.FuncAnimation (fig, animate, np.arange (0, int (tf / tstep)), init_func = init, interval = 1, blit = True, repeat = True)

#mywriter = animation.FFMpegWriter(fps = 24)
#anim.save ('qm3_evol_lambda_%f_H0_%f.mp4' % (l1, H0), writer = mywriter)

plt.show ()

