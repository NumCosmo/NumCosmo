/***************************************************************************
 *            nc_wl_ellipticity.c
 *
 *  Tue Jun 24 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_wl_ellipticity.c
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

/**
 * NcWLEllipticity:
 *
 * Reduced-shear transformations of the complex ellipticity.
 *
 * Collects the pure, stateless transformations that map an intrinsic complex
 * ellipticity to the sheared (observed) one and back, together with the
 * log-determinant of the corresponding Jacobian. A separate set of functions is
 * provided for each ellipticity convention (#NcGalaxyWLObsEllipConv) instead of
 * selecting it through object state, so the math can be inlined directly into
 * the weak-lensing hot loops.
 *
 * Each transformation comes in two flavours, following the
 * ncm_complex_set() / ncm_complex_set_c() convention: the plain-named functions
 * documented here take #NcmComplex and are introspectable, while the
 * _c-suffixed counterparts take a native C99 complex double and are inlined.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/wl/nc_wl_ellipticity.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

/**
 * nc_wl_ellipticity_apply_shear_trace:
 * @g: reduced shear as a #NcmComplex
 * @chi: intrinsic distortion as a #NcmComplex
 * @chi_obs: output observed distortion as a #NcmComplex
 *
 * Applies the reduced shear @g to the intrinsic distortion @chi in the trace
 * (distortion) convention, storing the result in @chi_obs.
 */
void
nc_wl_ellipticity_apply_shear_trace (const NcmComplex *g, const NcmComplex *chi, NcmComplex *chi_obs)
{
  ncm_complex_set_c (chi_obs, nc_wl_ellipticity_apply_shear_trace_c (ncm_complex_c (g), ncm_complex_c (chi)));
}

/**
 * nc_wl_ellipticity_apply_shear_inv_trace:
 * @g: reduced shear as a #NcmComplex
 * @chi_obs: observed distortion as a #NcmComplex
 * @chi: output intrinsic distortion as a #NcmComplex
 *
 * Recovers the intrinsic distortion @chi from the observed distortion @chi_obs
 * under the reduced shear @g, in the trace (distortion) convention.
 */
void
nc_wl_ellipticity_apply_shear_inv_trace (const NcmComplex *g, const NcmComplex *chi_obs, NcmComplex *chi)
{
  ncm_complex_set_c (chi, nc_wl_ellipticity_apply_shear_inv_trace_c (ncm_complex_c (g), ncm_complex_c (chi_obs)));
}

/**
 * nc_wl_ellipticity_lndet_jac_trace:
 * @g: reduced shear as a #NcmComplex
 * @chi_obs: observed distortion as a #NcmComplex
 *
 * Computes the natural logarithm of the absolute value of the Jacobian
 * determinant of the intrinsic -> observed map at @chi_obs, in the trace
 * (distortion) convention.
 *
 * Returns: the log-determinant of the shear Jacobian.
 */
gdouble
nc_wl_ellipticity_lndet_jac_trace (const NcmComplex *g, const NcmComplex *chi_obs)
{
  return nc_wl_ellipticity_lndet_jac_trace_c (ncm_complex_c (g), ncm_complex_c (chi_obs));
}

/**
 * nc_wl_ellipticity_apply_shear_trace_det:
 * @g: reduced shear as a #NcmComplex
 * @e: intrinsic ellipticity as a #NcmComplex
 * @e_obs: output observed ellipticity as a #NcmComplex
 *
 * Applies the reduced shear @g to the intrinsic ellipticity @e in the
 * trace-determinant (ellipticity) convention, storing the result in @e_obs.
 */
void
nc_wl_ellipticity_apply_shear_trace_det (const NcmComplex *g, const NcmComplex *e, NcmComplex *e_obs)
{
  ncm_complex_set_c (e_obs, nc_wl_ellipticity_apply_shear_trace_det_c (ncm_complex_c (g), ncm_complex_c (e)));
}

/**
 * nc_wl_ellipticity_apply_shear_inv_trace_det:
 * @g: reduced shear as a #NcmComplex
 * @e_obs: observed ellipticity as a #NcmComplex
 * @e: output intrinsic ellipticity as a #NcmComplex
 *
 * Recovers the intrinsic ellipticity @e from the observed ellipticity @e_obs
 * under the reduced shear @g, in the trace-determinant (ellipticity)
 * convention.
 */
void
nc_wl_ellipticity_apply_shear_inv_trace_det (const NcmComplex *g, const NcmComplex *e_obs, NcmComplex *e)
{
  ncm_complex_set_c (e, nc_wl_ellipticity_apply_shear_inv_trace_det_c (ncm_complex_c (g), ncm_complex_c (e_obs)));
}

/**
 * nc_wl_ellipticity_lndet_jac_trace_det:
 * @g: reduced shear as a #NcmComplex
 * @e_obs: observed ellipticity as a #NcmComplex
 *
 * Computes the natural logarithm of the absolute value of the Jacobian
 * determinant of the intrinsic -> observed map at @e_obs, in the
 * trace-determinant (ellipticity) convention.
 *
 * Returns: the log-determinant of the shear Jacobian.
 */
gdouble
nc_wl_ellipticity_lndet_jac_trace_det (const NcmComplex *g, const NcmComplex *e_obs)
{
  return nc_wl_ellipticity_lndet_jac_trace_det_c (ncm_complex_c (g), ncm_complex_c (e_obs));
}
