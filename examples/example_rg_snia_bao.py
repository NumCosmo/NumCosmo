#!/usr/bin/python2

from math import *
from gi.repository import Numcosmo as Nc
import matplotlib.pyplot as plt

#
#  Initializing the library objects, this must be called before 
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
#  A new Modelset with cosmo as the HICosmo model to be used.
#
mset = Nc.MSet ()
mset.set (cosmo)

#
#  Setting parameters Omega_c, Omega_x and w to be fitted (and change 
#  parameter Omega_x -> Omega_k).
#
#cosmo.de_omega_x2omega_k ()
cosmo.props.Omegac_fit = True
cosmo.props.Omegax_fit = True
cosmo.props.w_fit = True

#
#  A new Distance object optimized to redshift 2.
#
dist = Nc.Distance (zf = 2.0)

#
#  A new Data object from distance modulus catalogs.
#  A new Data object from BAO catalogs.
#
snia = Nc.DataDistMu.new (dist, Nc.DataSNIAId.SIMPLE_UNION2_1)
bao = Nc.Data.bao_new (dist, Nc.DataBaoId.A_EISENSTEIN2005)

#
#  A new Dataset with snia and bao set.
#
dset = Nc.Dataset ()
dset.append_data (snia)
dset.append_data (bao)

#
#  Creating a Likelihood from the Dataset.
#
lh = Nc.Likelihood (dataset = dset)

#
#  Creating a Fit object of type NLOPT using the fitting algorithm ln-neldermead to
#  fit the Modelset mset using the Likelihood lh and using a numerical differentiation
#  algorithm (NUMDIFF_FORWARD) to obtain the gradient (if needed).
#
fit = Nc.Fit.new (Nc.FitType.NLOPT, "ln-neldermead", lh, mset, Nc.FitGradType.NUMDIFF_FORWARD)

#
#  Running the fitter printing messages.
#
fit.run (Nc.FitRunMsgs.SIMPLE)

#
#  Printing fitting informations.
#
fit.log_info ()

#
#  Calculating the parameters covariance using numerical differentiation.
#
fit.numdiff_m2lnL_covar ()

#
#  Printing the covariance matrix.
# 
fit.log_covar ()

#
#  Creating a new Likelihood ratio test object.
#  First we create two PIndex indicating which parameter
#    we are going to study.
# 
p1 = Nc.MSetPIndex.new (cosmo.id (), Nc.HICosmoDEParams.OMEGA_C)
p2 = Nc.MSetPIndex.new (cosmo.id (), Nc.HICosmoDEXCDMParams.W)

lhr2d = Nc.LHRatio2d.new (fit, p1, p2)

#
#  Calculating the confidence region using the Likelihood ratio test.
#  Also calculate using the Fisher matrix approach.
#
cr_rg = lhr2d.conf_region (0.6826, 300.0, Nc.FitRunMsgs.SIMPLE)
fisher_rg = lhr2d.fisher_border (0.6826, 300.0, Nc.FitRunMsgs.SIMPLE)

cr_p1array = cr_rg.p1.dup_array ()
cr_p2array = cr_rg.p2.dup_array ()

fisher_p1array = fisher_rg.p1.dup_array ()
fisher_p2array = fisher_rg.p2.dup_array ()

#
#  Ploting the confidence regions obtained from both methods.
#

plt.title ("Confidence regions (%f)" % (cr_rg.clevel * 100.0))
plt.plot (cr_p1array, cr_p2array, 'r', label="Likelihood Ratio")
plt.plot (fisher_p1array, fisher_p2array, 'b-', label="Fisher Matrix")

plt.xlabel(r'$\Omega_c$')
plt.ylabel(r'$w$')

plt.legend(loc=4)

plt.savefig ("snia_bao_rg_omegac_w.png")
