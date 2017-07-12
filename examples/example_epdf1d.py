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
# Sampling from a Gaussian distribution
#
cut_l  = -0.4
cut_u  =  0.4
peak1  =  0.5
sigma1 =  0.2
peak2  = -0.5
sigma2 =  0.2

# Cumulative distribution function
def true_cdf (x):
  return 0.5 * (norm.cdf (x, peak1, sigma1) + norm.cdf (x, peak2, sigma2))


cdf_cut_l = true_cdf (cut_l)
cdf_cut_u = 1.0 - true_cdf (cut_u)
cdf_cut   = cdf_cut_l + cdf_cut_u
rbnorm = 1.0 - cdf_cut

# Probability density function  
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

#
# Plotting the probability density function.
#
fig = plot.subplot ()
plot.title ("PDF")
fig.plot (x_a, p_a, label = "auto-bw")
fig.plot (x_a, p_rot_a, label = "RoT-bw")
fig.plot (x_a, true_p (x_a), label = "true dist")

fig.legend(loc = "upper center")

plot.savefig ("epdf1d_pdf.svg")
plot.clf ()

#
# Plotting the relative difference of the reconstructed distributions and the true one.
#
fig = plot.subplot ()
plot.title ("PDF relative difference with respect to the true distribution")
fig.plot (x_a, np.abs (np.array ((p_a - true_p (x_a)) / true_p (x_a))), label = "auto-bw")
fig.plot (x_a, np.abs (np.array ((p_rot_a - true_p (x_a)) / true_p (x_a))), label = "RoT-bw")
fig.set_ylim ([1.0e-6, 1.0e1])
fig.grid ()

fig.legend(loc = "upper right")
fig.set_yscale ("log")

plot.savefig ("epdf1d_pdf_diff.svg")
plot.clf ()

#
# Plotting the cumulative distribution.
#
fig = plot.subplot ()
plot.title ("CDF")
fig.plot (x_a, pdf_a, label = "auto-bw")
fig.plot (x_a, pdf_rot_a, label = "RoT-bw")
fig.plot (x_a, (true_cdf (x_a) - cdf_cut_l) / rbnorm, label = "true dist")

fig.legend(loc = "upper left")

plot.savefig ("epdf1d_cdf.svg")
plot.clf ()

#
# Plotting the relative difference of the reconstructed cumulative distributions and the true one.
#
fig = plot.subplot ()
plot.title ("CDF relative difference with respect to the true distribution")
fig.plot (x_a, np.abs (pdf_a     / ( (true_cdf (x_a) - cdf_cut_l) / rbnorm ) - 1.0), label = "auto-bw")
fig.plot (x_a, np.abs (pdf_rot_a / ( (true_cdf (x_a) - cdf_cut_l) / rbnorm ) - 1.0), label = "RoT-bw")
fig.grid ()

fig.legend(loc = "upper right")
fig.set_yscale ("log")

plot.savefig ("epdf1d_cdf_diff.svg")
plot.clf ()

#
# Plotting the inverse cumulative distribution.
#
fig = plot.subplot ()
plot.title ("Inverse CDF")
fig.plot (u_a, inv_pdf_a, label = "auto-bw")
fig.plot (u_a, inv_pdf_rot_a, label = "RoT-bw")

fig.legend(loc = "upper left")

plot.savefig ("epdf1d_invcdf.svg")
plot.clf ()

