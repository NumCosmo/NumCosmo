#!/usr/bin/python2

import gi
import math
gi.require_version('NumCosmo', '1.0')
gi.require_version('NumCosmoMath', '1.0')

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
Ncm.cfg_init ()

#
# New ModelBuilder object, defines a new model NcHIPrimExample implementing
# the Nc.HIPrim abstract class.
# 
mb = Ncm.ModelBuilder.new (Nc.HIPrim, "NcHIPrimExample", "A example primordial model")

#
# New parameter A_s to describe the spectrum amplitud (it is usually better
# to work with ln(10^10As))
#
mb.add_sparam ("A_s", "As", 0.0, 1.0, 0.1, 0.0, 1.0e-9, Ncm.ParamType.FREE)

#
# New parameter n_s to describe the spectral index
#
mb.add_sparam ("n_s", "ns", 0.5, 1.5, 0.1, 0.0, 0.96,   Ncm.ParamType.FREE)

#
# Creates a new GObject, it is not a Python object yet!
#
GNcHIPrimExample = mb.create ()

#
# Gets the Python version of the object (.pytype)
#
NcHIPrimExample = GNcHIPrimExample.pytype

#
# Workaraound to make the pygobject "see" the new object
#
GObject.new (GNcHIPrimExample)
GObject.type_register (NcHIPrimExample)

#
# Creating a new class implementing our object NcHIPrimExample
#
class PyHIPrimExample (NcHIPrimExample):
  #
  # Defining some property which is not part of the model paramers.
  # All model parameters must be defined by the ModelBuilder.
  #
  some_property = GObject.Property (type = str)
  
  #
  # Implementing the adiabatic spectrum function
  #
  def do_lnSA_powspec_lnk (self, lnk):
    lnk0 = self.get_lnk_pivot ()
    lnka = lnk - lnk0
    As = self.props.As
    ns = self.props.ns
    
    return (ns - 1.0) * lnka + math.log (1.0e10 * As) - 10.0 * math.log (10.0)

#
# Register our new Python class PyHIPrimExample
#
GObject.type_register(PyHIPrimExample)

#
# Creating a new instance of PyHIPrimExample
#
prim = PyHIPrimExample ()

print "# As     = ", prim.props.As
print "# P(k=1) = ", prim.SA_powspec_k (1.0)


#
#  New CLASS backend precision object
#
cbe_prec = Nc.CBEPrecision.new ()

#
#  Set precision parameters
#
#ser.set_property_from_key_file (cbe_prec, "cl_ref.pre")

#
#  New CLASS backend object
#
cbe = Nc.HIPertBoltzmannCBE.prec_new (cbe_prec)
cbe.set_TT_lmax (1024)
cbe.set_target_Cls (Nc.DataCMBDataType.TT)
cbe.props.use_lensed_Cls = True

#
# Makes sure that cbe constain the default values
#
#cbe.props.precision.assert_default ()

#
#  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm
#
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm")
cosmo.omega_x2omega_k ()
cosmo.param_set_by_name ("Omegak", 0.0)

#
# Preparing the Class backend object
#
cbe.prepare (prim, cosmo)
Cls = Ncm.Vector.new (1024)
cbe.get_TT_Cls (Cls)

l = 0
for Cl in Cls.dup_array ():
  print l, Cl
  l = l + 1


