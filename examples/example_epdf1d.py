#!/usr/bin/python2

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import matplotlib as mpl
import matplotlib.pyplot as plot
from scipy.stats import norm
from scipy.stats import beta
import numpy as np
import math
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()
np.random.seed (0)

#
# n = number of points to reconstruct the distribution
# Sampling from a beta distribution
#
cut_l  = -0.4
cut_u  =  0.4
peak1  =  0.5
sigma1 =  0.2
peak2  = -0.5
sigma2 =  0.2

def true_cdf (x):
  return 0.5 * (norm.cdf (x, peak1, sigma1) + norm.cdf (x, peak2, sigma2))

cdf_cut = true_cdf (cut_l) + 1.0 - true_cdf (cut_u)
rbnorm = 1.0 - cdf_cut

def true_inv_cdf (u):
  return norm.ppf (u, 0.0, 1.0)
  
def true_p (x):
  return 0.5 * (norm.pdf (x, peak1, sigma1) + norm.pdf (x, peak2, sigma2)) / rbnorm

n = 800000
s = np.concatenate ((np.random.normal (peak1, sigma1, n), np.random.normal (peak2, sigma2, n)), axis = 0)

sa = []
for si in s:
  if si >= cut_l and si <= cut_u:
    sa.append (si)
s = sa
n = len (s)

print "# Number of points = %u" % (n)

#
# Creating a new Ncm.StatsDist1dEPDF object with
# Maximum of 200 saved samples (exceeding points are joined with 
# their nearest points), standard deviation scale of 0.1 factor
# and minimum distance scale of 0.001.
# 
epdf     = Ncm.StatsDist1dEPDF.new_full (2000, Ncm.StatsDist1dEPDFBw.AUTO, 0.01, 0.001)
epdf_rot = Ncm.StatsDist1dEPDF.new_full (2000, Ncm.StatsDist1dEPDFBw.ROT, 0.01,  0.001)

#
# Adding the points to the epdf object.
#
for i in range (n):
  epdf.add_obs (s[i])
  epdf_rot.add_obs (s[i])

#
# Preparing the object from the given sample.
#
epdf.prepare ()
epdf_rot.prepare ()

#
# Plotting the results.
#
p_a = []
p_rot_a = []
pdf_a = []
pdf_rot_a = []
x_a = []
inv_pdf_a = []
inv_pdf_rot_a = []
u_a = []

for i in range (1000):
  x = epdf.xi + (epdf.xf - epdf.xi) / 999.0 * i;
  u = 1.0 / 1000.0 * i
  
  x_a.append (x)
  u_a.append (u)
  
  p_a.append (epdf.eval_p (x))
  p_rot_a.append (epdf_rot.eval_p (x))
  
  pdf_a.append (epdf.eval_pdf (x))
  pdf_rot_a.append (epdf_rot.eval_pdf (x))

  inv_pdf_a.append (epdf.eval_inv_pdf (u))
  inv_pdf_rot_a.append (epdf_rot.eval_inv_pdf (u))

  #a1 = epdf.eval_inv_pdf (u)
  #a2 = beta.ppf (u, a = 2.0, b = 5.0)
  #print u, a1, a2, math.fabs ((a1 - a2) / a2)


#
# Plotting the cumulative distribution.
#
fig = plot.subplot ()
plot.title ("PDF")
fig.plot (x_a, p_a, label = "autobw")
fig.plot (x_a, p_rot_a, label = "rotbw")
fig.plot (x_a, true_p (x_a), label = "true dist")

fig.legend(loc = "upper right")

plot.savefig ("epdf1d_pdf.pdf")
plot.clf ()

#
# Plotting the cumulative distribution.
#
fig = plot.subplot ()
plot.title ("PDF diff")
fig.plot (x_a, np.abs (np.array ((p_a - true_p (x_a)) / true_p (x_a))), label = "autobw")
fig.plot (x_a, np.abs (np.array ((p_rot_a - true_p (x_a)) / true_p (x_a))), label = "rotbw")
fig.set_ylim ([1.0e-6, 1.0e1])
fig.grid ()

fig.legend(loc = "upper right")
fig.set_yscale ("log")

plot.savefig ("epdf1d_pdf_diff.pdf")
plot.clf ()

#
# Plotting the cumulative distribution.
#
fig = plot.subplot ()
plot.title ("CDF")
fig.plot (x_a, pdf_a, label = "autobw")
fig.plot (x_a, pdf_rot_a, label = "rotbw")
fig.plot (x_a, (true_cdf (x_a) - cdf_cut) / rbnorm, label = "true dist")

fig.legend(loc = "upper right")

plot.savefig ("epdf1d_cdf.pdf")
plot.clf ()

#
# Plotting the inverse cumulative distribution.
#
fig = plot.subplot ()
plot.title ("Inverse CDF")
fig.plot (u_a, inv_pdf_a, label = "autobw")
fig.plot (u_a, inv_pdf_rot_a, label = "rotbw")
fig.plot (u_a, beta.ppf (u_a, a = 2.0, b = 5.0), label = "true dist")

fig.legend(loc = "upper right")

plot.savefig ("epdf1d_invcdf.pdf")
plot.clf ()

