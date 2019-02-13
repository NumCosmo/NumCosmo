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

font = {'size'   : 20}

matplotlib.rc('font', **font)

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
l1                 = 1.0
H0                 = -15.0e-1
p.props.abstol     = 0.0
p.props.reltol     = 1.0e-7
p.props.nknots     = 400
p.props.noboundary = False
p.set_property ("lambda", l1)
offset = 5.0
center = (10.0 + offset) * 0.0 + 0.0

print ("# lambda = % 22.15g, basis a = % 22.15g, acs a = % 22.15g nu = % 22.15g mu = % 22.15g" % (p.get_lambda (), p.get_basis_a (), p.get_acs_a (), p.get_nu (), p.get_mu ()))

#psi0 = Nc.HIQG1DGauss.new (center, 1.0, 1.0, H0)
psi0 = Nc.HIQG1DExp.new (p.get_acs_a (), 2.0, H0)
#psi0 = Nc.HIQG1DExp.new (3.0, 2.0, H0)

#print (psi0.eval (1.0))
npp   = 1000
sim   = True
tstep = 1.0e-3
tf    = 4.5
xf    = center + 100.0
xfp   = center + 20.0
xi    = (center - 7.0) * 0.0 + 0.0
dSdiv = 10.0

#p.set_init_cond_gauss (psi0, xi, xf)
p.set_init_cond_exp (psi0, xi, xf)
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
ax.set_ylim (-1.0, 1.0)
ax.grid ()

#ax2.set_xlim (0.0, 10.0)
#ax2.set_ylim (1.0e-2, 1e2)
ax2.set_xlim (-0.1, 12.0)
ax2.set_ylim (-8.0, 8.0)
#ax2.set_xscale ('symlog', linthreshy=0.1)
#ax2.set_yscale ('log')

ax.set_xlabel (r'$V$')
ax2.set_xlabel (r'$V$')
ax2.set_ylabel (r'$P_V$')
ax2.grid ()

N = 4
ttl = ax.text (.1, 1.005, '', transform = ax.transAxes)
lines = []
lines.append (ax.plot ([], [], label=r'$\sqrt{\psi^*\psi}$', animated=True)[0])
lines.append (ax.plot ([], [], label=r'$\mathrm{Re}(\psi)$', animated=True)[0])
lines.append (ax.plot ([], [], label=r'$\mathrm{Im}(\psi)$', animated=True)[0])
lines.append (ax.plot ([], [], animated=True)[0])
#lines.append (ax.plot ([], [], label=r'$\partial_aS$', animated=True)[0])
ax.legend (loc='best')

fig.tight_layout()
lines.append (ttl)

lines.append (ax2.plot ([], [], 'bo', label=r'Bohm', animated=True)[0])
#lines.append (ax2.plot ([], [], 'ro', label=r'$\langle a(t)\rangle$', animated=True)[0])
lines.append (ax2.plot ([], [], 'ro', label=r'Semi-classical', animated=True)[0])
ax2.legend (loc='best')

ta    = [0.0]
traj  = []
trajP = []
x_t   = []
y_t   = []

Pini = H0
Vini = 3.0
Hini = ((Pini * Vini)**2 + l1) / Vini**2
t0   = -0.5 * Pini * Vini / Hini 

def a_sc(t):
  return math.sqrt (4.0 * Hini * (t - t0)**2 + l1 / Hini)
def p_sc(t):
  return 2.0 * Hini * (t - t0) / a_sc (t)

for i in range (n1):
  lines.append (ax2.plot ([], [])[0])
  qi = [p.Bohm (i)]
  pi = [p.Bohm_p (i)]
  traj.append (qi)
  trajP.append (pi)

lines.append (ax2.plot ([], [])[0])
#traj.append ([p.int_xrho_0_inf ()])
traj.append ([a_sc (0.0)])
trajP.append ([p_sc (0.0)])

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
  
#    qm = p.int_xrho_0_inf ()
    
    tfa = [tf] * n1
    lines[5].set_data ([p.Bohm (i) for i in range (n1)], [p.Bohm_p (i) for i in range (n1)])
    lines[6].set_data ([a_sc (tf)], [p_sc (tf)])
#    lines[5].set_data (tfa, [p.Bohm (i) for i in range (n1)])
#    lines[6].set_data ([tf], [qm])

    ta.append (tf)
    for i in range (n1):
      traj[i].append (p.Bohm (i))
      trajP[i].append (p.Bohm_p (i))
#      lines[i + 7].set_data (ta, traj[i])
      lines[i + 7].set_data (traj[i], trajP[i])
    
#    traj[n1].append (qm)
#    lines[n1 + 7].set_data (ta, traj[n1])
#    lines[n1 + 7].set_data (traj[n1], trajP[n1])

    traj[n1].append (a_sc (tf))
    trajP[n1].append (p_sc (tf))
#    lines[n1 + 7].set_data (ta, traj[n1])
    lines[n1 + 7].set_data (traj[n1], trajP[n1])

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

anim = animation.FuncAnimation (fig, animate, np.arange (0, int (tf / tstep)), init_func = init, interval = 1, blit = True, repeat = False)

#mywriter = animation.FFMpegWriter(fps = 24)
#anim.save ('cmp1_lambda_%f_H0_%f.mp4' % (l1, H0), writer = mywriter)

plt.show ()

