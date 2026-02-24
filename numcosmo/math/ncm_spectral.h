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
#include <numcosmo/math/ncm_matrix.h>

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
NcmSpectral *ncm_spectral_new_with_max_order (guint max_order);
NcmSpectral *ncm_spectral_ref (NcmSpectral *spectral);

void ncm_spectral_free (NcmSpectral *spectral);
void ncm_spectral_clear (NcmSpectral **spectral);

void ncm_spectral_compute_chebyshev_coeffs (NcmSpectral *spectral, NcmSpectralF F, gdouble a, gdouble b, NcmVector *coeffs, gpointer user_data);
NcmVector *ncm_spectral_compute_chebyshev_coeffs_adaptive (NcmSpectral *spectral, NcmSpectralF F, gdouble a, gdouble b, guint k_min, gdouble tol, gpointer user_data);

void ncm_spectral_chebT_to_gegenbauer_alpha1 (NcmVector *c, NcmVector *g);
void ncm_spectral_chebT_to_gegenbauer_alpha2 (NcmVector *c, NcmVector *g);

gdouble ncm_spectral_gegenbauer_alpha1_eval (NcmVector *c, gdouble x);
gdouble ncm_spectral_gegenbauer_alpha2_eval (NcmVector *c, gdouble x);

gdouble ncm_spectral_chebyshev_eval (NcmVector *a, gdouble t);
gdouble ncm_spectral_chebyshev_deriv (NcmVector *a, gdouble t);

NcmMatrix *ncm_spectral_get_proj_matrix (guint N);
NcmMatrix *ncm_spectral_get_x_matrix (guint N);
NcmMatrix *ncm_spectral_get_x2_matrix (guint N);
NcmMatrix *ncm_spectral_get_d_matrix (guint N);
NcmMatrix *ncm_spectral_get_x_d_matrix (guint N);
NcmMatrix *ncm_spectral_get_d2_matrix (guint N);
NcmMatrix *ncm_spectral_get_x_d2_matrix (guint N);
NcmMatrix *ncm_spectral_get_x2_d2_matrix (guint N);

NCM_INLINE void ncm_spectral_compute_proj_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff);
NCM_INLINE void ncm_spectral_compute_x_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff);
NCM_INLINE void ncm_spectral_compute_x2_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff);
NCM_INLINE void ncm_spectral_compute_d_row (gdouble * restrict row_data, glong offset, gdouble coeff);
NCM_INLINE void ncm_spectral_compute_x_d_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff);
NCM_INLINE void ncm_spectral_compute_d2_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff);
NCM_INLINE void ncm_spectral_compute_x_d2_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff);
NCM_INLINE void ncm_spectral_compute_x2_d2_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff);

G_END_DECLS

#endif /* _NCM_SPECTRAL_H_ */


#ifndef _NCM_SPECTRAL_INLINE_H_
#define _NCM_SPECTRAL_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__
#ifndef NUMCOSMO_GIR_SCAN

G_BEGIN_DECLS

NCM_INLINE void
ncm_spectral_compute_proj_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* General formula: three entries with rational function coefficients */
  const gdouble value_0 = coeff / (2.0 * (kd + 1.0));
  const gdouble value_1 = -coeff * (kd + 2.0) / ((kd + 1.0) * (kd + 3.0));
  const gdouble value_2 = coeff / (2.0 * (kd + 3.0));

  row_data[offset]     += value_0;
  row_data[offset + 2] += value_1;
  row_data[offset + 4] += value_2;

  /* Special case: k=0 has an additional contribution */
  if (k == 0)
    row_data[offset] += coeff * 0.5;
}

NCM_INLINE void
ncm_spectral_compute_x_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* General formula: up to four entries */
  if (k >= 1)
    row_data[offset - 1] += coeff / (4.0 * (kd + 1.0));

  row_data[offset + 1] += -coeff / (4.0 * (kd + 3.0));
  row_data[offset + 3] += -coeff / (4.0 * (kd + 1.0));
  row_data[offset + 5] += coeff / (4.0 * (kd + 3.0));

  /* Special cases */
  if (k == 0)
    row_data[offset + 1] += coeff * 0.25;  /* Additional 1/4 at column 1 = k+1 */
  else if (k == 1)
    row_data[offset - 1] += coeff * 0.125;  /* Additional 1/8 at column 0 = k-1 */
}

NCM_INLINE void
ncm_spectral_compute_x2_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff)
{
  const gdouble kd  = (gdouble) k;
  const gdouble kp2 = kd + 2.0;

  /* General formula: up to five entries */
  if (k >= 2)
    row_data[offset - 2] += coeff / (8.0 * (kd + 1.0));

  row_data[offset]     += coeff / (4.0 * (kd + 1.0) * (kd + 3.0));
  row_data[offset + 2] += -coeff * kp2 / (4.0 * (kd + 1.0) * (kd + 3.0));
  row_data[offset + 4] += -coeff / (4.0 * (kd + 1.0) * (kd + 3.0));
  row_data[offset + 6] += coeff / (8.0 * (kd + 3.0));

  /* Special cases */
  if (k == 0)
  {
    row_data[offset]     += coeff / 12.0; /* Additional 1/12 at column 0 = k */
    row_data[offset + 2] += coeff / 8.0;  /* Additional 1/8 at column 2 = k+2 */
  }
  else if (k == 1)
  {
    row_data[offset] += coeff / 16.0; /* Additional 1/16 at column 1 = k */
  }
  else if (k == 2)
  {
    row_data[offset - 2] += coeff / 24.0; /* Additional 1/24 at column 0 = k-2 */
  }
}

NCM_INLINE void
ncm_spectral_compute_d_row (gdouble * restrict row_data, glong offset, gdouble coeff)
{
  /* Two non-zero entries with opposite signs */
  row_data[offset + 1] += coeff;
  row_data[offset + 3] += -coeff;
}

NCM_INLINE void
ncm_spectral_compute_x_d_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* Three non-zero entries with rational function coefficients */
  const gdouble value_0 = coeff * kd / (2.0 * (kd + 1.0));
  const gdouble value_1 = coeff * (kd + 2.0) / ((kd + 1.0) * (kd + 3.0));
  const gdouble value_2 = -coeff * (kd + 4.0) / (2.0 * (kd + 3.0));

  row_data[offset]     += value_0;
  row_data[offset + 2] += value_1;
  row_data[offset + 4] += value_2;
}

NCM_INLINE void
ncm_spectral_compute_d2_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  row_data[offset + 2] += coeff * 2.0 * (kd + 2.0);
}

NCM_INLINE void
ncm_spectral_compute_x_d2_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff)
{
  const gdouble kd = (gdouble) k;

  /* Two non-zero entries in this row */
  row_data[offset + 1] += coeff * kd;
  row_data[offset + 3] += coeff * (kd + 4.0);
}

NCM_INLINE void
ncm_spectral_compute_x2_d2_row (gdouble * restrict row_data, glong k, glong offset, gdouble coeff)
{
  const gdouble kd  = (gdouble) k;
  const gdouble kp2 = kd + 2.0;

  /* Three non-zero entries with rational function coefficients */
  const gdouble value_0 = coeff * kd * (kd - 1.0) / (2.0 * (kd + 1.0));
  const gdouble value_1 = coeff * kp2 * (kp2 * kp2 - 3.0) / ((kd + 1.0) * (kd + 3.0));
  const gdouble value_2 = coeff * (kd + 4.0) * (kd + 5.0) / (2.0 * (kd + 3.0));

  row_data[offset]     += value_0;
  row_data[offset + 2] += value_1;
  row_data[offset + 4] += value_2;
}

G_END_DECLS

#endif /* NUMCOSMO_GIR_SCAN */
#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_SPECTRAL_INLINE_H_ */

