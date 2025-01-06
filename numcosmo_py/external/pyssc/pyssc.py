# Modified version of the original PySSC code, available at
# https://github.com/fabienlacasa/PySSC, LICENSE and README files are included in
# the same directory.

#
# IMPORT NECESSARY MODULES
#

"""Module to compute the Sij matrix for a given set of redshifts and kernels."""


import math
import copy

import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.special import spherical_jn as jn
from scipy.special import jv as Jn

try:
    import healpy
except ImportError:
    healpy = None

from . import numcosmo_class

##################################################

# Default values for redshift bin, cosmo parameters etc
pi = math.pi


def Sij(
    z_arr,
    kernels,
    cosmo,
    order=2,
    sky="full",
    method="classic",
    convention=0,
    precision=10,
    clmask=None,
    mask=None,
    mask2=None,
    var_tol=0.05,
    verbose=False,
    debug=False,
):
    """Wrapper routine to compute the Sij matrix.
    It calls different routines depending on the inputs : fullsky or partial sky,
    computation method.

    """

    test_zw(z_arr, kernels)

    # Full sky
    if sky.casefold() in ["full", "fullsky", "full sky", "full-sky"]:
        if method.casefold() in ["classic", "standard", "default", "std"]:
            Sij_array = Sij_fullsky(
                z_arr,
                kernels,
                cosmo,
                order=order,
                convention=convention,
                precision=precision,
            )
        elif method.casefold() in ["alternative", "alt"]:
            Sij_array = Sij_alt_fullsky(
                z_arr,
                kernels,
                cosmo,
                order=order,
                convention=convention,
            )
        else:
            raise ValueError(
                "Invalid string given for method parameter. Main possibilities: "
                "classic, alternative or AngPow (or variants, see code for details)."
            )
    # Partial sky
    elif sky.casefold() in ["psky", "partial sky", "partial-sky", "partial", "masked"]:
        test_mask(mask, clmask, mask2=mask2)
        if method.casefold() in ["classic", "standard", "default"]:
            Sij_array = Sij_psky(
                z_arr,
                kernels,
                cosmo,
                order=order,
                clmask=clmask,
                mask=mask,
                mask2=mask2,
                convention=convention,
                precision=precision,
                var_tol=var_tol,
                verbose=verbose,
                debug=debug,
            )
        elif method.casefold() in ["alt", "alternative"]:
            raise ValueError(
                "No implementation of the alternative method for partial sky. "
                "Use classic instead."
            )
        else:
            raise ValueError(
                "Invalid string given for method parameter. Main possibilities: "
                "classic, alternative or AngPow (or variants, see code for details)."
            )
    # Wrong geometry
    else:
        raise ValueError(
            "Invalid string given for sky geometry parameter. Main possibilities : "
            "full sky or partial sky (or abbreviations, see code for details)."
        )

    return Sij_array


# FULL SKY ROUTINES


def Sij_fullsky(
    z_arr,
    kernels,
    cosmo,
    order=2,
    convention=0,
    precision=10,
):
    """Routine to compute the Sij matrix in full sky. Standard computation method."""

    # Find number of redshifts and bins
    nbins = kernels.shape[0]

    # Get cosmology, comoving distances etc from dedicated auxiliary routine
    h, comov_dist, dcomov_dist, growth, psml = numcosmo_class.get_numcosmo(z_arr, cosmo)
    psml.prepare(cosmo)
    # Get element of z integration, depending on kernel convention
    dX_dz = get_dX_dz(comov_dist, dcomov_dist, convention=convention)
    # Compute normalisations
    Inorm = np.zeros(nbins)
    for i1 in range(nbins):
        integrand = dX_dz * kernels[i1, :] ** order
        Inorm[i1] = integrate.simpson(y=integrand, x=z_arr)

    # Compute U(i,k), numerator of Sij (integral of kernels**2 * matter )
    keq = 0.02 / h  # Equality matter radiation in 1/Mpc (more or less)
    klogwidth = 10  # Factor of width of the integration range. 10 seems ok
    kmin = min(keq, 1.0 / comov_dist.max()) / klogwidth  # 1e-4
    kmax = max(keq, 1.0 / comov_dist.min()) * klogwidth  # 1
    nk = 2**precision
    logkmin = np.log(kmin)
    logkmax = np.log(kmax)
    logk = np.linspace(logkmin, logkmax, num=nk)
    kk = np.exp(logk)  # logarithmic grid on k
    Pk = np.zeros(nk)
    for ik in range(nk):
        Pk[ik] = psml.eval(cosmo, 0.0, kk[ik])  # In Mpc^3
    Uarr = np.zeros((nbins, nk))
    for ibin in range(nbins):
        for ik in range(nk):
            kr = kk[ik] * comov_dist
            integrand = dX_dz * kernels[ibin, :] ** order * growth * np.sin(kr) / kr
            Uarr[ibin, ik] = integrate.simpson(y=integrand, x=z_arr)
    # Compute Sij finally
    Cl_zero = np.zeros((nbins, nbins))
    # For i<=j
    for ibin in range(nbins):
        U1 = Uarr[ibin, :] / Inorm[ibin]
        for jbin in range(ibin, nbins):
            U2 = Uarr[jbin, :] / Inorm[jbin]
            integrand = kk**2 * Pk * U1 * U2
            # Cl_zero[ibin,jbin] = 2/pi * integrate.simpson(integrand,kk)
            # linear integration
            Cl_zero[ibin, jbin] = (2 / pi) * integrate.simpson(
                y=integrand * kk, x=logk
            )  # log integration
    # Fill by symmetry
    for ibin in range(nbins):
        for jbin in range(nbins):
            Cl_zero[ibin, jbin] = Cl_zero[min(ibin, jbin), max(ibin, jbin)]
    Sij_array = Cl_zero / (4 * pi)

    return Sij_array


def Sij_alt_fullsky(
    z_arr,
    kernels,
    cosmo,
    order=2,
    convention=0,
):
    """Alternative routine to compute the Sij matrix in full sky."""

    # Find number of redshifts and bins
    nz = z_arr.size
    nbins = kernels.shape[0]

    # Get cosmology, comoving distances etc from dedicated auxiliary routine
    h, comov_dist, dcomov_dist, growth, _ = numcosmo_class.get_numcosmo(z_arr, cosmo)

    # Get element of z integration, depending on kernel convention
    dX_dz = get_dX_dz(comov_dist, dcomov_dist, convention=convention)

    keq = 0.02 / h  # Equality matter radiation in 1/Mpc (more or less)
    klogwidth = 10  # Factor of width of the integration range.
    # 10 seems ok ; going higher needs to increase nk_fft to reach convergence
    # (fine cancellation issue noted in Lacasa & Grain)
    kmin = min(keq, 1.0 / comov_dist.max()) / klogwidth
    kmax = max(keq, 1.0 / comov_dist.min()) * klogwidth
    nk_fft = (
        2**11
    )  # seems to be enough. Increase to test precision, reduce to speed up.
    k_4fft = np.linspace(
        kmin, kmax, nk_fft
    )  # linear grid on k, as we need to use an FFT
    Deltak = kmax - kmin
    Dk = Deltak / nk_fft
    Pk_4fft = np.zeros(nk_fft)
    for ik in range(nk_fft):
        Pk_4fft[ik] = cosmo.pk(k_4fft[ik], 0.0)  # In Mpc^3
    dr_fft = np.linspace(0, nk_fft // 2, nk_fft // 2 + 1) * 2 * pi / Deltak

    # Compute necessary FFTs and make interpolation functions
    fft0 = np.fft.rfft(Pk_4fft) * Dk
    dct0 = fft0.real
    Pk_dct = interp1d(dr_fft, dct0, kind="cubic")

    # Compute sigma^2(z1,z2)
    sigma2_nog = np.zeros((nz, nz))
    # First with P(k,z=0) and z1<=z2
    for iz in range(nz):
        r1 = comov_dist[iz]
        for jz in range(iz, nz):
            r2 = comov_dist[jz]
            rsum = r1 + r2
            rdiff = abs(r1 - r2)
            Icp0 = Pk_dct(rsum)
            Icm0 = Pk_dct(rdiff)
            sigma2_nog[iz, jz] = (Icm0 - Icp0) / (4 * pi**2 * r1 * r2)
    # Now fill by symmetry and put back growth functions
    sigma2 = np.zeros((nz, nz))
    for iz in range(nz):
        growth1 = growth[iz]
        for jz in range(nz):
            growth2 = growth[jz]
            sigma2[iz, jz] = sigma2_nog[min(iz, jz), max(iz, jz)] * growth1 * growth2

    # Compute normalisations
    Inorm = np.zeros(nbins)
    for i1 in range(nbins):
        integrand = dX_dz * kernels[i1, :] ** order
        Inorm[i1] = integrate.simpson(y=integrand, x=z_arr)

    # Compute Sij finally
    prefactor = sigma2 * (dX_dz * dX_dz[:, None])
    Sij_array = np.zeros((nbins, nbins))
    # For i<=j
    for i1 in range(nbins):
        for i2 in range(i1, nbins):
            integrand = prefactor * (
                kernels[i1, :] ** order * kernels[i2, :, None] ** order
            )
            Sij_array[i1, i2] = integrate.simpson(
                y=integrate.simpson(y=integrand, x=z_arr), x=z_arr
            ) / (Inorm[i1] * Inorm[i2])
    # Fill by symmetry
    for i1 in range(nbins):
        for i2 in range(nbins):
            Sij_array[i1, i2] = Sij_array[min(i1, i2), max(i1, i2)]

    return Sij_array


# PARTIAL SKY ROUTINES


def Sij_psky(
    z_arr,
    kernels,
    cosmo,
    order=2,
    clmask=None,
    mask=None,
    mask2=None,
    convention=0,
    precision=10,
    var_tol=0.05,
    verbose=False,
    debug=False,
):
    """Routine to compute the Sij matrix in partial sky."""

    # Find number of redshifts and bins
    nbins = kernels.shape[0]

    # compute Cl(mask) and fsky computed from user input (mask(s) or clmask)
    ell, cl_mask, fsky = get_mask_quantities(
        clmask=clmask, mask=mask, mask2=mask2, verbose=verbose
    )

    # Search of the best lmax for all later sums on ell
    lmax = find_lmax(ell, cl_mask, var_tol, debug=debug)
    if verbose:
        print(f"lmax = {lmax}")

    # Cut ell and Cl_mask to lmax, for all later computations
    cl_mask = cl_mask[: (lmax + 1)]
    nell = lmax + 1
    ell = np.arange(nell)  # 0..lmax

    # Get cosmology, comoving distances etc from dedicated auxiliary routine
    h, comov_dist, dcomov_dist, growth, psml = numcosmo_class.get_numcosmo(z_arr, cosmo)
    psml.prepare(cosmo)
    # Get element of z integration, depending on kernel convention
    dX_dz = get_dX_dz(comov_dist, dcomov_dist, convention=convention)

    # Compute normalisations
    Inorm = np.zeros(nbins)
    for i1 in range(nbins):
        integrand = dX_dz * kernels[i1, :] ** order
        Inorm[i1] = integrate.simpson(y=integrand, x=z_arr)

    # Full sky computation for debugging
    if debug:
        keq = 0.02 / h  # Equality matter radiation in 1/Mpc (more or less)
        klogwidth = 10  # Factor of width of the integration range. 10 seems ok
        kmin = min(keq, 1.0 / comov_dist.max()) / klogwidth
        kmax = max(keq, 1.0 / comov_dist.min()) * klogwidth
        nk = (
            2**precision
        )  # 10 seems to be enough. Increase to test precision, reduce to speed up.
        # kk          = np.linspace(kmin,kmax,num=nk)
        # linear grid on k
        logkmin = np.log(kmin)
        logkmax = np.log(kmax)
        logk = np.linspace(logkmin, logkmax, num=nk)
        kk = np.exp(logk)  # logarithmic grid on k
        Pk = np.zeros(nk)
        for ik in range(nk):
            # Pk[ik] = cosmo.pk(kk[ik],0.)                              #In Mpc^3
            Pk[ik] = psml.eval(cosmo, 0.0, kk[ik])  # In Mpc^3
        Uarr = np.zeros((nbins, nk))
        for ibin in range(nbins):
            for ik in range(nk):
                kr = kk[ik] * comov_dist
                integrand = dX_dz * kernels[ibin, :] ** order * growth * np.sin(kr) / kr
                Uarr[ibin, ik] = integrate.simpson(y=integrand, x=z_arr)
        Cl_zero = np.zeros((nbins, nbins))
        # For i<=j
        for ibin in range(nbins):
            U1 = Uarr[ibin, :] / Inorm[ibin]
            for jbin in range(ibin, nbins):
                U2 = Uarr[jbin, :] / Inorm[jbin]
                integrand = kk**2 * Pk * U1 * U2
                # Cl_zero[ibin,jbin] = 2/pi * integrate.simpson(integrand,kk)
                # linear integration
                Cl_zero[ibin, jbin] = (
                    2 / pi * integrate.simpson(y=integrand * kk, x=logk)
                )  # log integration
        # Fill by symmetry
        for ibin in range(nbins):
            for jbin in range(nbins):
                Cl_zero[ibin, jbin] = Cl_zero[min(ibin, jbin), max(ibin, jbin)]

    # Compute U(i;k,ell) = int dX kernels(i,z)^order growth(z) j_ell(kk*r)
    keq = 0.02 / h  # Equality matter radiation in 1/Mpc (more or less)
    klogwidth = 10  # Factor of width of the integration range. 10 seems ok
    kmin = min(keq, 1.0 / comov_dist.max()) / klogwidth
    kmax = max(keq, 1.0 / comov_dist.min()) * klogwidth
    # print(kmax)
    # kmax = 0.005
    nk = (
        2**precision
    )  # 10 seems to be enough. Increase to test precision, reduce to speed up.
    # kk          = np.linspace(kmin,kmax,num=nk)                   #linear grid on k
    logkmin = np.log(kmin)
    logkmax = np.log(kmax)
    logk = np.linspace(logkmin, logkmax, num=nk)
    kk = np.exp(logk)  # logarithmic grid on k
    Pk = np.zeros(nk)
    for ik in range(nk):
        # Pk[ik] = cosmo.pk(kk[ik],0.)                              #In Mpc^3
        Pk[ik] = psml.eval(cosmo, 0.0, kk[ik])  # In Mpc^3
    Uarr = np.zeros((nbins, nk, nell))
    for ik in range(nk):
        kr = kk[ik] * comov_dist
        for ll in ell:
            bessel_jl = jn(ll, kr)
            for ibin in range(nbins):
                integrand = dX_dz * kernels[ibin, :] ** order * growth * bessel_jl
                Uarr[ibin, ik, ll] = integrate.simpson(y=integrand, x=z_arr)

    # Compute Cl(X,Y) = 2/pi \int kk^2 dkk P(kk) U(i;kk,ell)/
    # I_\mathrm{norm}(i) U(j;kk,ell)/I_\mathrm{norm}(j)
    Cl_XY = np.zeros((nbins, nbins, nell))
    for ll in ell:
        # For i<=j
        for ibin in range(nbins):
            U1 = Uarr[ibin, :, ll] / Inorm[ibin]
            for jbin in range(ibin, nbins):
                U2 = Uarr[jbin, :, ll] / Inorm[jbin]
                integrand = kk**2 * Pk * U1 * U2
                # Cl_XY[ibin,jbin,ll] = 2/pi * integrate.simpson(integrand,kk)
                # linear integration
                Cl_XY[ibin, jbin, ll] = (
                    2 / pi * integrate.simpson(y=integrand * kk, x=logk)
                )  # log integration
        # Fill by symmetry
        for ibin in range(nbins):
            for jbin in range(nbins):
                Cl_XY[ibin, jbin, ll] = Cl_XY[min(ibin, jbin), max(ibin, jbin), ll]

    if debug:
        truc = (Cl_zero - Cl_XY[:, :, 0]) / Cl_zero
        print(
            "Debug: minmax of relative difference Cl_zero vs Cl_XY(ell=0)",
            truc.min(),
            truc.max(),
        )

    # Finally sum over ell to get
    # Sij = sum_ell (2ell+1) C(ell,i,j) C(ell,mask) /(4pi*fsky)^2
    Sij_array = np.zeros((nbins, nbins))
    # For i<=j
    for ibin in range(nbins):
        for jbin in range(ibin, nbins):
            Sij_array[ibin, jbin] = (
                np.sum((2 * ell + 1) * cl_mask * Cl_XY[ibin, jbin, :])
                / (4 * pi * fsky) ** 2
            )
            if debug and ibin == 0 and jbin == 0:
                print("Debug: fsky,ell,cl_mask", fsky, ell, cl_mask)
                print(
                    "Debug: Sij computation",
                    Cl_XY[ibin, jbin, 0] / (4 * pi),
                    np.sum((2 * ell + 1) * cl_mask * Cl_XY[ibin, jbin, :])
                    / (4 * pi * fsky) ** 2,
                )

    # Fill by symmetry
    for ibin in range(nbins):
        for jbin in range(nbins):
            Sij_array[ibin, jbin] = Sij_array[min(ibin, jbin), max(ibin, jbin)]

    return Sij_array


def Sij_flatsky(
    z_arr,
    kernels,
    cosmo,
    bin_centres,
    theta,
    verbose=False,
):
    """Routine to compute Sij according to the flat-sky approximation
    See Eq. 16 of arXiv:1612.05958

    """

    # Check inputs
    test_zw(z_arr, kernels)
    # Find number of redshifts and bins
    nbins = kernels.shape[0]

    theta = theta * np.pi / 180.0  # converts in radians

    # Get cosmology, comoving distances etc from dedicated auxiliary routine

    h, comov_dist, _, growth, _ = numcosmo_class.get_numcosmo(z_arr, cosmo)

    keq = 0.02 / h  # Equality matter radiation in 1/Mpc (more or less)
    klogwidth = 10  # Factor of width of the integration range. 10 seems ok
    kmin = min(keq, 1.0 / comov_dist.max()) / klogwidth
    kmax = max(keq, 1.0 / comov_dist.min()) * klogwidth

    # kperp array
    kmin_perp = kmin
    kmax_perp = kmax
    nk_perp = 50
    lnkperp_arr = np.linspace(np.log10(kmin_perp), np.log10(kmax_perp), nk_perp)
    kperp_arr = 10 ** (lnkperp_arr)
    # kpar array
    kmin_par = kmin
    kmax_par = kmax
    nk_par = 50
    lnkpar_arr = np.linspace(np.log10(kmin_par), np.log10(kmax_par), nk_par)
    kpar_arr = 10 ** (lnkpar_arr)
    # k2 = kperp2 + kpar2
    k_arr = np.sqrt(kperp_arr[:, None] ** 2 + kpar_arr[None, :] ** 2)
    # growth      = np.zeros(nz)                              #Growth factor
    # for iz in range(nbins):
    #     growth[iz] = cosmo.scale_independent_growth_factor(z_arr[iz])

    if verbose:
        print("Computing flat-sky approximation")
    Sij_array = np.zeros((nbins, nbins))
    for ibin in range(nbins):
        z1 = bin_centres[ibin]
        r1 = cosmo.z_of_r([z1])[0][0]
        dr1 = (
            comov_dist[kernels[ibin, :] != 0].max()
            - comov_dist[kernels[ibin, :] != 0].min()
        )  # width of function function

        for jbin in range(nbins):
            z2 = bin_centres[jbin]
            r2 = cosmo.z_of_r([z2])[0][0]
            dr2 = (
                comov_dist[kernels[jbin, :] != 0].max()
                - comov_dist[kernels[jbin, :] != 0].min()
            )  # width of kernel

            z12 = np.mean([bin_centres[ibin], bin_centres[jbin]])

            growth = cosmo.scale_independent_growth_factor(z12)
            Pk = np.array(
                [
                    cosmo.pk(k_arr[i, j], 0.0)
                    for i in range(nk_perp)
                    for j in range(nk_par)
                ]
            )
            Pk = Pk.reshape(k_arr.shape)

            integ_kperp = (
                kperp_arr
                * 4.0
                * (Jn(1, kperp_arr * theta * r1) / kperp_arr / theta / r1)
                * (Jn(1, kperp_arr * theta * r2) / kperp_arr / theta / r2)
            )
            integ_kpar = jn(0, kpar_arr * dr1 / 2) * jn(0, kpar_arr * dr2 / 2)
            dSij = (
                integ_kperp[:, None]
                * integ_kpar[None, :]
                * growth
                * Pk
                * np.cos(kpar_arr[None, :] * abs(r1 - r2))
            )

            Sij_array[ibin, jbin] = (
                integrate.simpson(y=integrate.simpson(y=dSij, x=kpar_arr), x=kperp_arr)
                / 2.0
                / np.pi**2
            )

    return Sij_array


def get_dX_dz(comov_dist, dcomov_dist, convention=0):
    """Auxiliary routine to compute the element of integration for z integrals in
    Sij routines."""
    if (
        convention == 0
    ):  # Default: comoving volume dV = r^2 dr. Convention of Lacasa & Grain 2019
        dX_dz = comov_dist**2 * dcomov_dist
    elif (
        convention == 1
    ):  # Convention of Cosmosis, Euclid Forecasts (Blanchard et al 2020)
        dX_dz = dcomov_dist / comov_dist**2
    else:
        raise ValueError("convention must be either 0 or 1")

    return dX_dz


def get_mask_quantities(clmask=None, mask=None, mask2=None, verbose=False):
    """Auxiliary routine to compute different mask quantities (ell,Cl,fsky) for
    partial sky Sij routines.

    """
    if healpy is None:
        raise ImportError("Healpy is needed to compute Sij for partial sky.")

    if mask is None:  # User gives Cl(mask)
        if verbose:
            print("Using given Cls")
        if isinstance(clmask, str):
            cl_mask = healpy.fitsfunc.read_cl(str(clmask))
        elif isinstance(clmask, np.ndarray):
            cl_mask = clmask
        ell = np.arange(len(cl_mask))
        lmaxofcl = ell.max()
    else:  # User gives mask as a map
        if verbose:
            print("Using given mask map")
        if isinstance(mask, str):
            map_mask = healpy.read_map(mask, dtype=np.float64)
        elif isinstance(mask, np.ndarray):
            map_mask = mask
        nside = healpy.pixelfunc.get_nside(map_mask)
        lmaxofcl = 2 * nside
        if mask2 is None:
            map_mask2 = copy.copy(map_mask)
        else:
            if isinstance(mask2, str):
                map_mask2 = healpy.read_map(mask2, dtype=np.float64)
            elif isinstance(mask2, np.ndarray):
                map_mask2 = mask2
        cl_mask = healpy.anafast(map_mask, map2=map_mask2, lmax=lmaxofcl)
        ell = np.arange(lmaxofcl + 1)

    # Compute fsky from the mask
    fsky = np.sqrt(cl_mask[0] / (4 * np.pi))
    if verbose:
        print(f"f_sky = {fsky:.4f}")

    return ell, cl_mask, fsky


def find_lmax(ell, cl_mask, var_tol, debug=False):
    """Auxiliary routine to search the best lmax for all later sums on multipoles."""

    assert ell.ndim == 1, "ell must be a 1-dimensional array"
    assert cl_mask.ndim == 1, "cl_mask must be a 1-dimensional array"
    assert len(ell) == len(cl_mask), "ell and cl_mask must have the same size"
    lmaxofcl = ell.max()
    summand = (2 * ell + 1) / (4 * pi) * cl_mask
    var_target = np.sum(summand)
    # Initialisation
    lmax = 0
    var_est = np.sum(summand[: (lmax + 1)])
    while abs(var_est - var_target) / var_target > var_tol and lmax < lmaxofcl:
        lmax = lmax + 1
        var_est = np.sum(summand[: (lmax + 1)])
        if debug:
            print(
                "In lmax search",
                lmax,
                abs(var_est - var_target) / var_target,
                var_target,
                var_est,
            )
    lmax = min(lmax, lmaxofcl)  # make sure we didnt overshoot at the last iteration
    return lmax


def test_z(z_arr):
    """
    Assert redshift array has the good type, shape and values.
    """
    assert isinstance(z_arr, np.ndarray), "z_arr must be a numpy array."
    assert z_arr.ndim == 1, "z_arr must be a 1-dimensional array."
    assert z_arr.min() > 0, "z_arr must have values > 0."


def test_w(kernels, nz):
    """
    Assert kernels array has the good type and shape.
    """
    assert isinstance(kernels, np.ndarray), "kernels must be a numpy array."
    assert kernels.ndim == 2, "kernels must be a 2-dimensional array."
    assert kernels.shape[1] == nz, "kernels must have shape (nbins,nz)."


def test_zw(z_arr, kernels):
    """
    Assert redshift and kernels arrays have the good type, shape and values.
    """
    test_z(z_arr)
    nz = len(z_arr)
    test_w(kernels, nz)


# test_mask
def test_mask(mask, clmask, mask2=None):
    """Assert that either the mask or its Cl has been provided. If two masks are
    provided, check that they have the same resolution.
    """
    if healpy is None:
        raise ImportError("Healpy is needed to compute Sij for partial sky.")

    assert (mask is not None) or (
        clmask is not None
    ), "You need to provide either the mask or its angular power spectrum Cl."
    if mask is not None:
        assert isinstance(mask, str) or isinstance(
            mask, np.ndarray
        ), "mask needs to be either a filename or a numpy array"
    if clmask is not None:
        assert isinstance(clmask, str) or isinstance(
            clmask, np.ndarray
        ), "Clmask needs to be either a filename or a numpy array"
    if mask2 is not None:
        if isinstance(mask, str):
            map_mask = healpy.read_map(mask, dtype=np.float64)
        elif isinstance(mask, np.ndarray):
            map_mask = mask
        nside = healpy.pixelfunc.get_nside(map_mask)
        if isinstance(mask2, str):
            map_mask2 = healpy.read_map(mask2, dtype=np.float64)
        elif isinstance(mask2, np.ndarray):
            map_mask2 = mask2
        nside2 = healpy.pixelfunc.get_nside(map_mask2)
        assert (
            nside == nside2
        ), "The resolutions (nside) of both masks need to be the same."


# test_multimask
def test_multimask(multimask):
    """
    Assert that multimask has the good structure: a list of dictionnaries,
    all of the form {'mask':'mask.fits', 'kernels':kernels_array}.
    """
    assert isinstance(
        multimask, list
    ), "multimask needs to be a list (of dictionnaries)."
    for dico in multimask:
        assert isinstance(
            dico, dict
        ), "The elements of multimask must be dictionnaries."
        assert "mask" in dico.keys(), 'The dictionnaries must contain the key "mask".'
        assert (
            "kernels" in dico.keys()
        ), 'The dictionnaries must contain the key "kernels".'
        assert isinstance(
            dico.mask, str
        ), 'The key "mask" must contain a string (pointing to a healpix fits file).'
        assert isinstance(
            dico.kernels, np.ndarray
        ), 'The key "kernels" must contain a numpy array.'


# End of PySSC.py
