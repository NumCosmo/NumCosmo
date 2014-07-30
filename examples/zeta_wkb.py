#!/usr/bin/python2

from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm
import time

#
#  Initializing the library objects, this must be called before 
#  any other library function.
#
Ncm.cfg_init ()

#
#  New homogeneous and isotropic cosmological model NcHICosmoQGRW
#
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoQGRW")

cosmo.props.w      = 1.0e-12
cosmo.props.Omegar = 1.0e-5
cosmo.props.Omegaw = 1.0 - 1.0e-5

pert = Nc.HIPertAdiab.new ()
pert.props.reltol  = 1.0e-5
pert.set_mode_k (1.0);

alpha_minf = -cosmo.abs_alpha (1.0e-30)
alphaf     = cosmo.abs_alpha (1.0e25)
alphaf_wkb = pert.wkb_maxtime (cosmo, alpha_minf, -alphaf)
alphai     = pert.wkb_maxtime_prec (cosmo, alpha_minf, -alphaf)

print "# ", cosmo.x_alpha (alpha_minf), cosmo.x_alpha (alphai), alphaf_wkb, cosmo.x_alpha (alphaf)

pert.prepare_wkb (cosmo, alpha_minf, alphai)
pert.prepare_ode (cosmo, alphai, alphaf)

pert.set_stiff_solver (True)

pert.set_init_cond_wkb (cosmo, alphai)
  
for i in range (100000):
  alpha = alphai + (alphaf - alphai) / 100000.0 * (i + 1)
  pert.evolve (cosmo, alpha)
  (alphas, Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta) = pert.get_values ()
  (wkb_Re_zeta, wkb_Im_zeta, wkb_Re_Pzeta, wkb_Im_Pzeta) = pert.ode_zeta_Pzeta (cosmo, alphas)
  print alphas, Re_zeta, Im_zeta, Re_Pzeta, Im_Pzeta, wkb_Re_zeta, wkb_Im_zeta, wkb_Re_Pzeta, wkb_Im_Pzeta

