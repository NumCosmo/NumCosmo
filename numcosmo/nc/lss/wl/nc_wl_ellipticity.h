/***************************************************************************
 *            nc_wl_ellipticity.h
 *
 *  Tue Jun 24 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_wl_ellipticity.h
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _NC_WL_ELLIPTICITY_H_
#define _NC_WL_ELLIPTICITY_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_util.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#include <math.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

/**
 * NcWLEllipticityFrame:
 * @NC_WL_ELLIPTICITY_FRAME_CELESTIAL: celestial frame (reference frame).
 * @NC_WL_ELLIPTICITY_FRAME_CARTESIAN: Earth-view ("Euclidean"/image) frame.
 *
 * Handedness of the orthonormal frame whose basis vectors are used to resolve a
 * complex spin-2 ellipticity (or reduced shear) into its two components.
 *
 * Every sky position in NumCosmo is given as right ascension and declination, so
 * the position angle is always *measured* in the celestial convention ($\phi$
 * increases eastward, i.e. with increasing RA). This enum does not change how the
 * angle is measured; it only selects the basis onto which the spin-2 ellipticity
 * components are projected:
 *
 * - #NC_WL_ELLIPTICITY_FRAME_CELESTIAL: the holonomic frame of the celestial
 *   convention. The position angle $\phi_C$ increases eastward. Taken as the
 *   reference frame (identity).
 *
 * - #NC_WL_ELLIPTICITY_FRAME_CARTESIAN: the frame an observer on the ground sees
 *   when facing the object with their head pointing to the North celestial pole.
 *   The right-hand direction then points East, so the image $x$-axis increases
 *   towards the West (decreasing RA) -- the opposite handedness. This is the
 *   "Euclidean" image-plane view (the legacy name was EUCLIDEAN).
 *
 * The two frames differ by a parity (handedness) flip. The position angle
 * reverses sense,
 * $$ \phi_E = \pi - \phi_C \quad (\equiv -\phi_C \bmod \pi), $$
 * and the spin-2 ellipticity transforms by complex conjugation, i.e. its
 * second/cross component changes sign,
 * $$ \epsilon_E = \epsilon_C^{*} . $$
 * Since $e^{2 i \phi_E} = e^{-2 i \phi_C}$, flipping $\phi$ and conjugating
 * $\epsilon$ are the same operation; the parity is its own inverse, so a single
 * map converts CELESTIAL -> CARTESIAN and CARTESIAN -> CELESTIAL alike.
 *
 * The integer values match the legacy #NcGalaxyWLObsCoord (CELESTIAL = 0,
 * EUCLIDEAN = CARTESIAN = 1), so previously serialized data is read back
 * unchanged.
 *
 */
typedef enum _NcWLEllipticityFrame
{
  NC_WL_ELLIPTICITY_FRAME_CELESTIAL = 0,
  NC_WL_ELLIPTICITY_FRAME_CARTESIAN,
} NcWLEllipticityFrame;

/*
 * Reduced-shear transformations of the complex ellipticity, one set per
 * convention. The convention suffix matches #NcGalaxyWLObsEllipConv:
 *
 *   _trace      -> NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE     (distortion chi)
 *   _trace_det  -> NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET (ellipticity epsilon)
 *
 * In every function @g is the reduced shear, with the |g| > 1 (strong-lensing)
 * branch handled explicitly. apply_shear maps intrinsic -> observed,
 * apply_shear_inv maps observed -> intrinsic, and lndet_jac returns the natural
 * log of |det J| of the intrinsic -> observed map evaluated at the observed
 * ellipticity.
 *
 * Two flavours of each kernel are provided, following the ncm_complex_set /
 * ncm_complex_set_c convention:
 *
 *   - the plain-named functions take #NcmComplex and are GObject-introspectable
 *     (usable from Python);
 *   - the _c-suffixed functions take native C99 complex double, are inlined,
 *     and are meant for the C hot loops.
 */

/* Introspectable NcmComplex API (TRACE convention: distortion chi). */
void nc_wl_ellipticity_apply_shear_trace (const NcmComplex *g, const NcmComplex *chi, NcmComplex *chi_obs);
void nc_wl_ellipticity_apply_shear_inv_trace (const NcmComplex *g, const NcmComplex *chi_obs, NcmComplex *chi);
gdouble nc_wl_ellipticity_lndet_jac_trace (const NcmComplex *g, const NcmComplex *chi_obs);

/* Introspectable NcmComplex API (TRACE_DET convention: ellipticity epsilon). */
void nc_wl_ellipticity_apply_shear_trace_det (const NcmComplex *g, const NcmComplex *e, NcmComplex *e_obs);
void nc_wl_ellipticity_apply_shear_inv_trace_det (const NcmComplex *g, const NcmComplex *e_obs, NcmComplex *e);
gdouble nc_wl_ellipticity_lndet_jac_trace_det (const NcmComplex *g, const NcmComplex *e_obs);

#ifndef NUMCOSMO_GIR_SCAN

/* Express a celestial-frame complex spin-2 ellipticity @e in @frame:
 * epsilon_E = conj (epsilon_C) for CARTESIAN, identity for CELESTIAL. The map is
 * its own inverse (a parity flip), so it equally takes a @frame ellipticity back
 * to celestial; a separate frame_to_celestial entry point can be added if a call
 * site needs to read in the opposite direction. */
NCM_INLINE complex double nc_wl_ellipticity_celestial_to_frame_c (NcWLEllipticityFrame frame, complex double e);

/* Re-express the celestial position angle phi_C (as returned by
 * nc_halo_position_polar_angles) in @frame: CARTESIAN reverses the handedness,
 * phi_E = pi - phi_C; CELESTIAL leaves it unchanged. Using this angle in the
 * spin-2 factor e^{2 i phi} is equivalent to conjugating the ellipticity, so it
 * stays consistent with nc_wl_ellipticity_celestial_to_frame_c(). */
NCM_INLINE gdouble nc_wl_ellipticity_celestial_to_frame_angle (NcWLEllipticityFrame frame, gdouble phi);

/* Inline complex double kernels (TRACE convention: distortion chi). */
NCM_INLINE complex double nc_wl_ellipticity_apply_shear_trace_c (complex double g, complex double chi);
NCM_INLINE complex double nc_wl_ellipticity_apply_shear_inv_trace_c (complex double g, complex double chi_obs);
NCM_INLINE gdouble nc_wl_ellipticity_lndet_jac_trace_c (complex double g, complex double chi_obs);

/* Inline complex double kernels (TRACE_DET convention: ellipticity epsilon). */
NCM_INLINE complex double nc_wl_ellipticity_apply_shear_trace_det_c (complex double g, complex double e);
NCM_INLINE complex double nc_wl_ellipticity_apply_shear_inv_trace_det_c (complex double g, complex double e_obs);
NCM_INLINE gdouble nc_wl_ellipticity_lndet_jac_trace_det_c (complex double g, complex double e_obs);

#endif /* NUMCOSMO_GIR_SCAN */

G_END_DECLS

#endif /* _NC_WL_ELLIPTICITY_H_ */

#ifndef _NC_WL_ELLIPTICITY_INLINE_H_
#define _NC_WL_ELLIPTICITY_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__
#ifndef NUMCOSMO_GIR_SCAN

G_BEGIN_DECLS

/* Coordinate-frame handedness (parity). */

NCM_INLINE complex double
nc_wl_ellipticity_celestial_to_frame_c (NcWLEllipticityFrame frame, complex double e)
{
  return (frame == NC_WL_ELLIPTICITY_FRAME_CARTESIAN) ? conj (e) : e;
}

NCM_INLINE gdouble
nc_wl_ellipticity_celestial_to_frame_angle (NcWLEllipticityFrame frame, gdouble phi)
{
  return (frame == NC_WL_ELLIPTICITY_FRAME_CARTESIAN) ? (M_PI - phi) : phi;
}

/* TRACE convention (distortion chi). */

NCM_INLINE complex double
nc_wl_ellipticity_apply_shear_trace_c (complex double g, complex double chi)
{
  return (chi + g * (g * conj (chi) + 2.0)) /
         (1.0 + g * conj (g) + 2.0 * creal (g * conj (chi)));
}

NCM_INLINE complex double
nc_wl_ellipticity_apply_shear_inv_trace_c (complex double g, complex double chi_obs)
{
  return (chi_obs + g * (g * conj (chi_obs) - 2.0)) /
         (1.0 + g * conj (g) - 2.0 * creal (g * conj (chi_obs)));
}

NCM_INLINE gdouble
nc_wl_ellipticity_lndet_jac_trace_c (complex double g, complex double chi_obs)
{
  const complex double g_conj = conj (g);
  const gdouble abs_g2        = creal (g * g_conj);
  const gdouble lndet_jac_den = 3.0 * log (fabs (1.0 - 2.0 * creal (g_conj * chi_obs) + abs_g2));

  if (abs_g2 <= 1.0)
    return 3.0 * log1p (-abs_g2) - lndet_jac_den;
  else
    return 3.0 * log (abs_g2 - 1.0) - lndet_jac_den;
}

/* TRACE_DET convention (ellipticity epsilon). */

NCM_INLINE complex double
nc_wl_ellipticity_apply_shear_trace_det_c (complex double g, complex double e)
{
  if (cabs (g) <= 1.0)
    return (e + g) / (1.0 + conj (g) * e);
  else
    return (1.0 + g * conj (e)) / (conj (e) + conj (g));
}

NCM_INLINE complex double
nc_wl_ellipticity_apply_shear_inv_trace_det_c (complex double g, complex double e_obs)
{
  if (cabs (g) <= 1.0)
    return (e_obs - g) / (1.0 - conj (g) * e_obs);
  else
    return (1.0 - g * conj (e_obs)) / (conj (e_obs) - conj (g));
}

NCM_INLINE gdouble
nc_wl_ellipticity_lndet_jac_trace_det_c (complex double g, complex double e_obs)
{
  const complex double g_conj     = conj (g);
  const complex double e_obs_conj = conj (e_obs);
  const gdouble abs_g2            = creal (g * g_conj);

  if (abs_g2 <= 1.0)
  {
    const gdouble abs_e_obs2 = creal (e_obs * e_obs_conj);
    const gdouble ln_jac_num = 2.0 * log1p (-abs_g2);
    const gdouble ln_jac_den = 2.0 * log1p (-2.0 * creal (g_conj * e_obs) + abs_g2 * abs_e_obs2);

    return ln_jac_num - ln_jac_den;
  }
  else
  {
    const complex double e_obs_m_g = e_obs - g;
    const gdouble abs_e_obs_m_g2   = creal (e_obs_m_g * conj (e_obs_m_g));
    const gdouble ln_jac_num       = 2.0 * log (abs_g2 - 1.0);
    const gdouble ln_jac_den       = 2.0 * log (abs_e_obs_m_g2);

    return ln_jac_num - ln_jac_den;
  }
}

G_END_DECLS

#endif /* NUMCOSMO_GIR_SCAN */
#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_WL_ELLIPTICITY_INLINE_H_ */

