#!/usr/bin/python2

import gi
import math
import numpy as np
import matplotlib.pyplot as plt

gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from py_hiprim_example import PyHIPrimExample

#
# Script params
#
lmax = 2500

#
# Creating a new instance of PyHIPrimExample
#
prim = PyHIPrimExample ()

print "# As        = ", prim.props.As
print "# P (k = 1) = ", prim.SA_powspec_k (1.0)
print "# (a, b, c) = ( ", prim.props.a, ", ", prim.props.b, ", ", prim.props.c, " )"

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
#  Set precision parameters
#
#ser.set_property_from_key_file (cbe_prec, "cl_ref.pre")

#
#  New CLASS backend object
#
Bcbe = Nc.HIPertBoltzmannCBE.full_new (cbe)
Bcbe.set_TT_lmax (lmax)
Bcbe.set_target_Cls (Nc.DataCMBDataType.TT)
Bcbe.set_lensed_Cls (True)

#
# Makes sure that Bcbe contain the default values
#
#Bcbe.props.precision.assert_default ()

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
#reion.z_to_tau (cosmo)

#
# Preparing the Class backend object
#
Bcbe.prepare (prim, reion, cosmo)

Cls1 = Ncm.Vector.new (lmax + 1)
Cls2 = Ncm.Vector.new (lmax + 1)

Bcbe.get_TT_Cls (Cls1)

prim.props.a = 0
Bcbe.prepare (prim, reion, cosmo)
Bcbe.get_TT_Cls (Cls2)

Cls1_a = Cls1.dup_array ()
Cls2_a = Cls2.dup_array ()

Cls1_a = np.array (Cls1_a[2:])
Cls2_a = np.array (Cls2_a[2:])

ell = np.array (range (2, lmax + 1))

Cls1_a = ell * (ell + 1.0) * Cls1_a
Cls2_a = ell * (ell + 1.0) * Cls2_a

#
#  Ploting ionization history.
#

plt.title (r'Modified and non-modified $C_\ell$')
plt.xscale('log')
plt.plot (ell, Cls1_a, 'r', label="Modified")
plt.plot (ell, Cls2_a, 'b--', label="Non-modified")

plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_\ell$')
plt.legend(loc=2)

plt.savefig ("hiprim_Cls.png")

plt.clf ()
