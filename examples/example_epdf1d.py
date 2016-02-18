#!/usr/bin/python2

import gi
gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

import matplotlib as mpl
import matplotlib.pyplot as plot
from scipy.stats import norm
from scipy.stats import beta
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm
import numpy as np
import math

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

#
# n = number of points to reconstruct the distribution
# Sampling from a beta distribution
#
n = 2000
s = np.random.beta (2.0, 5.0, n)

#
# Creating a new Ncm.StatsDist1dEPDF object with
# Maximum of 200 saved samples (exceeding points are joined with 
# their nearest points), standard deviation scale of 0.1 factor
# and minimum distance scale of 0.001.
# 
epdf = Ncm.StatsDist1dEPDF.new (200, 0.1, 0.001)

#
# Adding the points to the epdf object.
#
for i in range (n):
  epdf.add_obs (s[i])

#
# Preparing the object from the given sample.
#
epdf.prepare ()

#
# Plotting the results.
#
pdf_a = []
x_a = []
inv_pdf_a = []
u_a = []

for i in range (1000):
  x = epdf.xi + (epdf.xf - epdf.xi) / 999.0 * i;
  u = 1.0 / 1000.0 * i
  
  x_a.append (x)
  u_a.append (u)
  
  pdf_a.append (epdf.eval_pdf (x))
  inv_pdf_a.append (epdf.eval_inv_pdf (u))
  a1 = epdf.eval_inv_pdf (u)
  a2 = beta.ppf (u, a = 2.0, b = 5.0)
  print u, a1, a2, math.fabs ((a1 - a2) / a2)

fig = plot.subplot ()

#
# Plotting the cumulative distribution.
#
plot.title ("CDF")
fig.plot (x_a, pdf_a)
fig.plot (x_a, beta.cdf (x_a, a = 2.0, b = 5.0))

plot.savefig ("epdf1d_cdf.png")
plot.clf ()

#
# Plotting the inverse cumulative distribution.
#
fig = plot.subplot ()
plot.title ("Inverse CDF")
fig.plot (u_a, inv_pdf_a)
fig.plot (u_a, beta.ppf (u_a, a = 2.0, b = 5.0))

plot.savefig ("epdf1d_invcdf.png")
plot.clf ()

