#!/usr/bin/env python

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import sys
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

p                  = Nc.HIQG1D.new ()
l1                 = 0.0
H0                 = -15.0e-1
p.props.abstol     = 0.0
p.props.reltol     = 1.0e-7
p.props.nknots     = 200
p.props.noboundary = False
p.set_property ("lambda", l1)
offset = 5.0
center = (10.0 + offset) * 0.0 + 0.0

print ("# lambda = % 22.15g, basis a = % 22.15g, acs a = % 22.15g nu = % 22.15g" % (p.get_lambda (), p.get_basis_a (), p.get_acs_a (), p.get_nu ()))

psi0 = Nc.HIQG1DGauss.new (center, 1.0001, 1.0, H0)
#psi0 = Nc.HIQG1DExp.new (p.get_acs_a (), 2.0, H0)
#psi0 = Nc.HIQG1DExp.new (3.0, 2.0, H0)

#print (psi0.eval (1.0))
npp   = 1000
sim   = True
tstep = 5.0e-3
tf    = 4.0
xf    = center + 50.0
xfp   = center + 50.0
xi    = (center - 7.0) * 0.0 + 0.0
dSdiv = 10.0

p.set_init_cond_gauss (psi0, xi, xf)
#p.set_init_cond_exp (psi0, xi, xf)
p.prepare ()
n1 = p.nBohm ()

q          = p.peek_knots ().dup_array ()
xa         = np.linspace (xi, xfp, npp)
psi        = [p.eval_psi (x) for x in xa]
max_Re_psi = max (np.abs (psi[:][0]))
max_Im_psi = max (np.abs (psi[:][1]))
yb         = max (max_Re_psi, max_Im_psi)

fig, [ax, ax2] = plt.subplots (1, 2, figsize = (16, 12))

ax.set_xlim (xi, xfp)
ax.set_ylim (-0.8, 0.8)
ax.grid ()

ax2.set_xlim (0.0, 10.0)
ax2.set_ylim (1.0e-2, 1e2)
ax2.set_xscale ('symlog', linthreshy=0.1)
ax2.set_yscale ('log')

ax.set_xlabel (r'$a(t)$')
ax2.set_xlabel (r'$t$')
ax2.set_ylabel (r'$a(t)$')

N = 4
ttl = ax.text (.1, 1.005, '', transform = ax.transAxes)
lines = []
lines.append (ax.plot ([], [], label=r'$\sqrt{\psi^*\psi}$', animated=True)[0])
lines.append (ax.plot ([], [], label=r'$\mathrm{Re}(\psi)$', animated=True)[0])
lines.append (ax.plot ([], [], label=r'$\mathrm{Im}(\psi)$', animated=True)[0])
lines.append (ax.plot ([], [], label=r'$\partial_aS$', animated=True)[0])
ax.legend (loc='best')

fig.tight_layout()
lines.append (ttl)

lines.append (ax2.plot ([], [], 'bo', label=r'Bohm', animated=True)[0])
lines.append (ax2.plot ([], [], 'ro', label=r'$\langle a(t)\rangle$', animated=True)[0])
ax2.legend (loc='best')

ta   = [0.0]
traj = []
x_t  = []
y_t  = []

for i in range (n1):
  lines.append (ax2.plot ([], [])[0])
  qi = [p.Bohm (i)]
  traj.append (qi)

lines.append (ax2.plot ([], [])[0])
traj.append ([p.int_xrho_0_inf ()])

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

  for p in lines:
    p.set_visible (False)

  return lines

if not sim:
  for i in np.arange (0, ceil (tf / tstep)):
    tf = tstep * i
    p.evol (tf)

    q  = p.peek_knots ().dup_array ()
    xa = np.linspace (xi, xfp, npp)

    psi = np.array ([p.eval_psi (x) for x in xa])
    rho = [np.sum (psi_i**2) for psi_i in psi]
    dS  = np.array ([0.0] * len (xa))
    #dS  = [p.eval_dS (x) for x in xa]
    y = []
  
    qm = p.int_xrho_0_inf ()
  
    y.append (np.sqrt (rho))
    y.append (psi[:][0])
    y.append (psi[:][1])
    y.append (dS / dSdiv)
    y.append (q[::s1])
    y.append (p.int_rho_0_inf ())
    y.append (qm)

    x_t.append (xa)
    y_t.append (y)

    ta.append (tf)
    for i in range (n1):
      traj[i].append (p.Bohm (i))
    traj[n1].append (qm)

def animate(i):
  tf = tstep * i

  if (i == 1):
    for pl in lines:
      pl.set_visible (True)

  if sim:
    p.evol (tf)
    
    q  = p.peek_knots ().dup_array ()
    xa = np.linspace (xi, xfp, npp)
    #print ("% 22.15g % 22.15g % 22.15g % 22.15g" % (q[0], q[1], q[2], q[3]))
  
    psi = np.array ([p.eval_psi (x) for x in xa])
    rho = np.array ([np.sum (psi_i**2) for psi_i in psi])
    dS  = np.array ([0.0] * len (xa))
    #dS  = np.array ([p.eval_dS (x) for x in xa])

    #print (dS / xa)
    
    lines[0].set_data (xa, np.sqrt (rho))
    lines[1].set_data (xa, psi[:,0])
    lines[2].set_data (xa, psi[:,1])
    lines[3].set_data (xa, dS / dSdiv)
    
    #print (dS)
  
    qm = p.int_xrho_0_inf ()
  
    tfa = [tf] * n1
    lines[5].set_data (tfa, [p.Bohm (i) for i in range (n1)])
    lines[6].set_data ([tf], [qm])

    ta.append (tf)
    for i in range (n1):
      traj[i].append (p.Bohm (i))
      lines[i + 7].set_data (ta, traj[i])
    
    traj[n1].append (qm)
    lines[n1 + 7].set_data (ta, traj[n1])

    ttl.set_text ("t = % .15f, norma = % .15f" % (tf, p.int_rho_0_inf ()))
  
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

