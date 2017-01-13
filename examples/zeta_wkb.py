#!/usr/bin/python2

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import time
import cProfile
import math
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

#
#  New homogeneous and isotropic cosmological model NcHICosmoQGRW
#
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoQGRW")

cosmo.props.w      = 1.0e-16
cosmo.props.Omegar = 1.0e-5
cosmo.props.Omegaw = 1.0 - 1.0e-5

pert = Nc.HIPertAdiab.new ()
pert.props.reltol  = 1.0e-10
pert.set_mode_k (1.0);

#wkb_prec = pert.reltol
wkb_prec = 1.0e-10

alpha_minf = -cosmo.abs_alpha (1.0e-36)
alphaf     = cosmo.abs_alpha (1.0e25)
alphaf_wkb = pert.wkb_maxtime (cosmo, alpha_minf, -alphaf)
alphai     = pert.wkb_maxtime_prec (cosmo, Nc.HIPertWKBCmp.POTENTIAL, wkb_prec, alpha_minf, -alphaf)
alphai     = -cosmo.abs_alpha (cosmo.x_alpha (alphai) * 1.0)

print "# Required precision ", pert.reltol
print "# wkb ini = % 7.5e, ewkb ini = % 7.5e, wkb end = % 7.5e, end = % 7.5e" % (cosmo.x_alpha (alpha_minf), cosmo.x_alpha (alphai), cosmo.x_alpha (alphaf_wkb), cosmo.x_alpha (alphaf))

pert.prepare_wkb (cosmo, wkb_prec, alpha_minf, -alphaf)
pert.set_stiff_solver (True)

print "# WKB prepared"

pert.set_init_cond_wkb (cosmo, alphai)

for i in range (10000):
  alpha = alphai + (alphaf - alphai) / 10000.0 * (i + 1)
  pert.evolve (cosmo, alpha)
  (alphas, Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta) = pert.get_values ()
  (wkb_Re_zeta, wkb_Im_zeta, wkb_Re_Pzeta, wkb_Im_Pzeta) = pert.wkb_zeta_Pzeta (cosmo, alphas)
  #print alphas, Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta, wkb_Re_zeta, wkb_Im_zeta, wkb_Re_Pzeta, wkb_Im_Pzeta
  abs_zeta = math.hypot (Re_zeta, Im_zeta)
  abs_wkb_zeta = math.hypot (wkb_Re_zeta, wkb_Im_zeta)
  print "%f %e %g %g %e %e" % (alphas, cosmo.x_alpha (alphas), abs_zeta, abs_wkb_zeta, abs ((abs_zeta - abs_wkb_zeta)/abs_zeta), abs ((Re_zeta - wkb_Re_zeta)/Re_zeta))


