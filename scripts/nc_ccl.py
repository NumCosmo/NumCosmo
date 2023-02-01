try:
    import gi
    gi.require_version('NumCosmo', '1.0')
    gi.require_version('NumCosmoMath', '1.0')
except:
    pass

from gi.repository import GObject
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

import pyccl
import math
import numpy as np

def create_nc_obj(ccl_cosmo, prec = 1.0e-7, dist_z_max = 15.0, ps_nln_z_max = 10.0, k_min = 1.0e-6, k_max = 1.0e3):

    # Creates the base HI object    
    cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDECpl{'massnu-length':<0>}")
    cosmo.omega_x2omega_k ()
    cosmo.param_set_by_name ("H0",        ccl_cosmo['h']*100)
    cosmo.param_set_by_name ("Omegak",    ccl_cosmo['Omega_k'])
    cosmo.param_set_by_name ("w0",        ccl_cosmo['w0'])
    cosmo.param_set_by_name ("w1",        ccl_cosmo['wa'])
    cosmo.param_set_by_name ("Omegab",    ccl_cosmo['Omega_b'])
    cosmo.param_set_by_name ("Omegac",    ccl_cosmo['Omega_c'])
    cosmo.param_set_by_name ("ENnu",      ccl_cosmo['Neff'])
    cosmo.param_set_by_name ("Tgamma0",   ccl_cosmo['T_CMB'])

    # Creates the HI Primordial object
    hiprim = Nc.HIPrimPowerLaw.new ()
    hiprim.param_set_by_name ("n_SA", ccl_cosmo['n_s'])

    cosmo.add_submodel (hiprim)

    # Creates the distance object optimized up to z = 15.0
    dist = Nc.Distance.new (15.0)
    dist.prepare (cosmo)

    # Checking if neutrinos are compatible
    if isinstance (ccl_cosmo['m_nu'],(list,np.ndarray)):
        for m_nu_i in ccl_cosmo['m_nu']:
            if m_nu_i != 0.0:
                raise ValueError ("Massive neutrinos are not supported")
    else:
        if ccl_cosmo['m_nu'] != 0:
            raise ValueError ("Massive neutrinos are not supported")

    # Creating the transfer/linear power spectrum
    tf = None
    ps_lin = None

    if ccl_cosmo._config_init_kwargs['transfer_function'] == "eisenstein_hu":
        tf = Nc.TransferFuncEH.new ()
        tf.props.CCL_comp = True
        
        ps_lin = Nc.PowspecMLTransfer.new (tf)
        ps_lin.prepare (cosmo)
    else:
        raise ValueError ("Transfer function type `' not supported" % (ccl_cosmo._config_init_kwargs['transfer_function']))    
    
    if not math.isnan (ccl_cosmo['A_s']):
        hiprim.param_set_by_name ("ln10e10ASA", math.log (1.0e10 * ccl_cosmo['A_s']))
    else:
        A_s = math.exp (hiprim.param_get_by_name ("ln10e10ASA")) * 1.0e-10
        fact = (ccl_cosmo['sigma8'] / ps_lin.sigma_tophat_R (cosmo, prec, 0.0, 8.0 / cosmo.h ()))**2    
        hiprim.param_set_by_name ("ln10e10ASA", math.log (1.0e10 * A_s * fact))
        
    ps_nln = None
    if ccl_cosmo._config_init_kwargs['matter_power_spectrum'] == 'halofit':
        ps_nln = Nc.PowspecMNLHaloFit.new (ps_lin, ps_nln_z_max, prec)
        ps_nln.prepare (cosmo)

    hmfunc = None

    if ps_lin:            
        psf = Ncm.PowspecFilter.new (ps_lin, Ncm.PowspecFilterType.TOPHAT)
        psf.set_best_lnr0 ()
    
        if ccl_cosmo._config_init_kwargs['mass_function'] == 'tinker10':
            hmf_T10 = Nc.MultiplicityFuncTinkerMeanNormalized.new ()
            hmfunc  = Nc.HaloMassFunction.new (dist, psf, hmf_T10)
            hmfunc.prepare (cosmo)

    return cosmo, dist, ps_lin, ps_nln, hmfunc


def ccl_cosmo_set_high_prec(ccl_cosmo):
    pyccl.gsl_params.INTEGRATION_EPSREL        = 1.0e-13
    pyccl.gsl_params.ODE_GROWTH_EPSREL         = 1.0e-13
    pyccl.gsl_params.N_ITERATION               = 10000
    pyccl.gsl_params.INTEGRATION_SIGMAR_EPSREL = 1.0e-9
    pyccl.spline_params.A_SPLINE_NLOG          = 1000
    pyccl.spline_params.A_SPLINE_NA            = 1000
    pyccl.spline_params.A_SPLINE_NA_PK         = 1000
    pyccl.spline_params.A_SPLINE_NLOG_PK       = 1000
    pyccl.spline_params.N_K                    = 1000
    pyccl.spline_params.K_MIN                  = 1.0e-6
    pyccl.spline_params.K_MAX                  = 1.0e3
    
    pyccl.spline_params.A_SPLINE_NLOG_SM       = 100
    pyccl.spline_params.A_SPLINE_NA_SM         = 100
    pyccl.spline_params.LOGM_SPLINE_NM         = 300

if __name__ == "__main__":
    
    Omega_c = 0.25
    Omega_b = 0.05
    Omega_k = 0.0
    h       = 0.7
    A_s     = 2.1e-9
    n_s     = 0.96
    Neff    = 3.046
    w0      = -1.0
    wa      = 0.0
     
    ccl_cosmo = pyccl.Cosmology(
        Omega_c=Omega_c, Omega_b=Omega_b, Neff=Neff,
        h=h, A_s = A_s, n_s=n_s, Omega_k=Omega_k,
        w0=w0, wa=wa, transfer_function='eisenstein_hu')
    
    Ncm.cfg_init ()
    
    create_nc_obj (ccl_cosmo)


# Missing function in CCL

def dsigmaM_dlnM(cosmo, M, a):

    cosmo.compute_sigma()

    logM = np.log10(np.atleast_1d(M))
    status = 0
    dsigMdlnM, status = pyccl.lib.dlnsigM_dlogM_vec(cosmo.cosmo, a, logM,
                                               len(logM), status)
    if np.ndim(M) == 0:
        dsigMdlnM = dsigMdlnM[0]
    return dsigMdlnM
