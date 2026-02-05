/***************************************************************************
 *            ncm_spectral.h
 *
 *  Tue Feb 04 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NCM_SPECTRAL_H_
#define _NCM_SPECTRAL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>

G_BEGIN_DECLS

/**
 * NcmSpectralF:
 * @user_data: user data
 * @x: point to evaluate
 *
 * Function to be used in spectral computations.
 *
 * Returns: the value of the function at @x
 */
typedef gdouble (*NcmSpectralF) (gpointer user_data, gdouble x);

#define NCM_TYPE_SPECTRAL (ncm_spectral_get_type ())

G_DECLARE_FINAL_TYPE (NcmSpectral, ncm_spectral, NCM, SPECTRAL, GObject)

NcmSpectral *ncm_spectral_new (void);
NcmSpectral *ncm_spectral_ref (NcmSpectral *spectral);

void ncm_spectral_free (NcmSpectral *spectral);
void ncm_spectral_clear (NcmSpectral **spectral);

void ncm_spectral_compute_chebyshev_coeffs (NcmSpectral *spectral, NcmSpectralF F, gdouble a, gdouble b, NcmVector *coeffs, gpointer user_data);

void ncm_spectral_chebT_to_gegenbauer_alpha1 (NcmVector *c, NcmVector *g);
void ncm_spectral_chebT_to_gegenbauer_alpha2 (NcmVector *c, NcmVector *g);

gdouble ncm_spectral_gegenbauer_alpha1_eval (NcmVector *c, gdouble x);
gdouble ncm_spectral_gegenbauer_alpha2_eval (NcmVector *c, gdouble x);

gdouble ncm_spectral_chebyshev_eval (NcmVector *a, gdouble t);
gdouble ncm_spectral_chebyshev_deriv (NcmVector *a, gdouble t);

G_END_DECLS

#endif /* _NCM_SPECTRAL_H_ */

