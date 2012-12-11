#!/usr/bin/python2

from math import *
from gi.repository import Numcosmo as Nc
import matplotlib.pyplot as plt

#
#  Initialize the library objects, this must be called before 
#  any other library function.
#
Nc.cfg_init ()

#
#  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm 
#
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm")

#
#  Setting values for the cosmological model, those not set stay in the
#  default values. Remeber to use the _orig_ version to set the original
#  parameters in case when a reparametrization is used.
#

#
# OO-like
#
cosmo.props.H0      = 70.0
cosmo.props.Omegab  = 0.05
cosmo.props.Omegac  = 0.25
cosmo.props.Omegax  = 0.70
cosmo.props.Tgamma0 = 2.72
cosmo.props.ns      = 1.0
cosmo.props.sigma8  = 0.9
cosmo.props.w       = -1.0

#
#  Create a new Modelset and set cosmo as the HICosmo model to be used.
#
mset = Nc.MSet ()
mset.set (cosmo)

#
#  Set parameters Omega_c and w to be fitted and change parameter
#  Omega_x -> Omega_k.
#
cosmo.de_omega_x2omega_k ()
cosmo.props.Omegac_fit = True
cosmo.props.w_fit = True

#
#  Create a new Distance object optimized to redshift 2.
#
dist = Nc.Distance (zf = 2.0)

#
#  Create a new Data object from distance modulus catalogs.
#
snia = Nc.DataDistMu.new (dist, Nc.DataSNIAId.SIMPLE_UNION2_1)

#
#  Create a new Dataset and add snia to it.
#
dset = Nc.Dataset ()
dset.append_data (snia)

#
#  Create a Likelihood from the Dataset.
#
lh = Nc.Likelihood (dataset = dset)

#
#  Create a Fit object of type NLOPT using the fitting algorithm ln-neldermead to
#  fit the Modelset mset using the Likelihood lh and using a numerical differentiation
#  algorithm (NUMDIFF_FORWARD) to obtain the gradient (if needed).
#
fit = Nc.Fit.new (Nc.FitType.NLOPT, "ln-neldermead", lh, mset, Nc.FitGradType.NUMDIFF_FORWARD)

#
#  Run the fitter printing messages.
#
fit.run (Nc.FitRunMsgs.SIMPLE)

#
#  Print fitting informations.
#
fit.log_info ()

#
#  Calculates the parameters covariance using numerical differentiation.
#
fit.numdiff_m2lnL_covar ()

#
#  Print the covariance matrix.
# 
fit.log_covar ()
