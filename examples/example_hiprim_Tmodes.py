#!/usr/bin/python2

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import math
import numpy as np
import matplotlib.pyplot as plt

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init ()

#
# Script parameters
#
# Maximum multipole
lmax = 2500

#
# Creating a new instance of HIPrimPowerLaw
#
prim = Nc.HIPrimPowerLaw.new ()
prim.props.T_SA_ratio = 2.0e-1

#
#  New CLASS backend precision object
#  Lets also increase k_per_decade_primordial since we are
#  dealing with a modified spectrum.
#
cbe_prec = Nc.CBEPrecision.new ()
cbe_prec.props.k_per_decade_primordial = 50.0

#
#  New CLASS backend object
#
cbe = Nc.CBE.prec_new (cbe_prec)

#
#  New CLASS backend object
#
Bcbe = Nc.HIPertBoltzmannCBE.full_new (cbe)
Bcbe.set_TT_lmax (lmax)
# Setting which CMB data to use
Bcbe.set_target_Cls (Nc.DataCMBDataType.TT)
# Setting if the lensed Cl's are going to be used or not.
Bcbe.set_lensed_Cls (True)
# Setting if the tensor contribution is going to be used or not. 
Bcbe.set_tensor (True)

#
#  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm
#
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm")
cosmo.omega_x2omega_k ()
cosmo.param_set_by_name ("Omegak", 0.0)

#
#  New homogeneous and isotropic reionization object
#
reion = Nc.HIReionCamb.new ()

#
# Adding submodels to the main cosmological model.
#
cosmo.add_submodel (reion)
cosmo.add_submodel (prim)

#
# Preparing the Class backend object
#
Bcbe.prepare (cosmo)

Cls1 = Ncm.Vector.new (lmax + 1)
Cls2 = Ncm.Vector.new (lmax + 1)

Bcbe.get_TT_Cls (Cls1)

prim.props.T_SA_ratio = 1.0e-15
Bcbe.prepare (cosmo)
Bcbe.get_TT_Cls (Cls2)

Cls1_a = Cls1.dup_array ()
Cls2_a = Cls2.dup_array ()

Cls1_a = np.array (Cls1_a[2:])
Cls2_a = np.array (Cls2_a[2:])

ell = np.array (range (2, lmax + 1))

Cls1_a = ell * (ell + 1.0) * Cls1_a
Cls2_a = ell * (ell + 1.0) * Cls2_a

#
#  Ploting the TT angular power spcetrum
#

plt.title (r'With and without tensor contribution to $C_\ell$')
plt.xscale('log')
plt.plot (ell, Cls1_a, 'r', label="with-T")
plt.plot (ell, Cls2_a, 'b--', label="without-T")

plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_\ell$')
plt.legend(loc=2)

plt.savefig ("hiprim_tensor_Cls.svg")

plt.clf ()
