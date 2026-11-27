/***************************************************************************
 *            ncm_spectral_panels.h
 *
 *  Sat Oct 10 2026
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

#ifndef _NCM_SPECTRAL_PANELS_H_
#define _NCM_SPECTRAL_PANELS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/algebra/ncm_vector.h>
#include <numcosmo/ncm/algebra/ncm_matrix.h>
#include <numcosmo/ncm/algebra/ncm_spectral.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPECTRAL_PANELS (ncm_spectral_panels_get_type ())

/* Largest N of a uniform set evaluated in the power basis */
#define NCM_SPECTRAL_PANELS_POWER_N_MAX (9)

G_DECLARE_FINAL_TYPE (NcmSpectralPanels, ncm_spectral_panels, NCM, SPECTRAL_PANELS, GObject)

struct _NcmSpectralPanels
{
  /*< private >*/
  GObject parent_instance;
  guint fdim;
  guint n_panels;
  GArray *edges;          /* gdouble, n_panels + 1 entries once a panel is added */
  GPtrArray *pool;        /* GArray of gdouble, the first n_panels hold the panels' fdim rows of N_i */
  gdouble **panel_data;   /* pool->len entries, the first n_panels the data of those arrays */
  GVariant *init_panels;  /* the construct-only property, consumed by constructed() */
  gboolean init_uniform;  /* the construct-only property, applied by constructed() */
  gboolean uniform;       /* set by ncm_spectral_panels_set_uniform(), with the fields below */
  gdouble uniform_x0;     /* the first edge */
  gdouble uniform_inv_h;  /* the number of panels over the range */
  gdouble uniform_last;   /* the index of the last panel */
  guint uniform_N;        /* the number of coefficients in every row */
  gdouble *uniform_power; /* N <= NCM_SPECTRAL_PANELS_POWER_N_MAX: power-basis rows, row r of panel i at (i fdim + r) N */
};

NcmSpectralPanels *ncm_spectral_panels_new (guint fdim);
NcmSpectralPanels *ncm_spectral_panels_ref (NcmSpectralPanels *sp);

void ncm_spectral_panels_free (NcmSpectralPanels *sp);
void ncm_spectral_panels_clear (NcmSpectralPanels **sp);

void ncm_spectral_panels_append (NcmSpectralPanels *sp, gdouble a, gdouble b, NcmMatrix *coeffs);
void ncm_spectral_panels_append_data (NcmSpectralPanels *sp, gdouble a, gdouble b, const gdouble *c, guint N, guint tda);
void ncm_spectral_panels_reset (NcmSpectralPanels *sp);
void ncm_spectral_panels_set_uniform (NcmSpectralPanels *sp);
gboolean ncm_spectral_panels_is_uniform (NcmSpectralPanels *sp);

guint ncm_spectral_panels_get_fdim (NcmSpectralPanels *sp);
guint ncm_spectral_panels_get_n_panels (NcmSpectralPanels *sp);
void ncm_spectral_panels_get_range (NcmSpectralPanels *sp, gdouble *a, gdouble *b);
const gdouble *ncm_spectral_panels_peek_panel (NcmSpectralPanels *sp, guint i, gdouble *a, gdouble *b, guint *N);
NcmMatrix *ncm_spectral_panels_get_panel (NcmSpectralPanels *sp, guint i, gdouble *a, gdouble *b);
guint ncm_spectral_panels_find_panel (NcmSpectralPanels *sp, gdouble x);

void ncm_spectral_panels_eval (NcmSpectralPanels *sp, gdouble x, NcmVector *fval);

NCM_INLINE gdouble ncm_spectral_panels_eval_row (NcmSpectralPanels *sp, guint r, gdouble x);

gdouble ncm_spectral_panels_eval_row_on (NcmSpectralPanels *sp, guint i, guint r, gdouble x);
void ncm_spectral_panels_square_cheb_moments (NcmSpectralPanels *sp, guint r, gdouble a, gdouble b, guint n, GArray **moments);
gboolean ncm_spectral_panels_restrict (NcmSpectralPanels *sp, NcmSpectral *spectral, gdouble a, gdouble b, NcmMatrix **coeffs);
gboolean ncm_spectral_panels_restrict_rows (NcmSpectralPanels *sp, NcmSpectral *spectral, gdouble a, gdouble b, guint r0, guint n, NcmMatrix **coeffs);
NcmSpectralPanels *ncm_spectral_panels_new_sum_rows (NcmSpectral *spectral, guint n, NcmSpectralPanels * const *sps, const guint *rows, const gdouble *scales, const gdouble *x_min, const gdouble *x_max, gdouble abstol);

G_END_DECLS

#endif /* _NCM_SPECTRAL_PANELS_H_ */

#ifndef _NCM_SPECTRAL_PANELS_CLENSHAW_H_
#define _NCM_SPECTRAL_PANELS_CLENSHAW_H_
#ifndef __GTK_DOC_IGNORE__
#ifndef NUMCOSMO_GIR_SCAN

/* The series of N coefficients c at s by the Clenshaw recurrence; always inlined, so that
 * each constant N gets its own unrolled copy */
static inline __attribute__ ((always_inline)) gdouble

_ncm_spectral_panels_clenshaw_n (const gdouble *c, const guint N, const gdouble s)
{
  const gdouble two_s = 2.0 * s;
  gdouble b_1         = 0.0;
  gdouble b_2         = 0.0;
  gint n;

  for (n = (gint) N - 1; n >= 1; n--)
  {
    const gdouble b_0 = two_s * b_1 - b_2 + c[n];

    b_2 = b_1;
    b_1 = b_0;
  }

  return s * b_1 - b_2 + c[0];
}

/* The polynomial with the N power-basis coefficients a at s by Horner's rule; always inlined */
static inline __attribute__ ((always_inline)) gdouble

_ncm_spectral_panels_horner_n (const gdouble *a, const guint N, const gdouble s)
{
  gdouble p = a[N - 1];
  gint k;

  for (k = (gint) N - 2; k >= 0; k--)
    p = p * s + a[k];

  return p;
}

#endif /* NUMCOSMO_GIR_SCAN */
#endif /* __GTK_DOC_IGNORE__ */
#endif /* _NCM_SPECTRAL_PANELS_CLENSHAW_H_ */

#ifndef _NCM_SPECTRAL_PANELS_INLINE_H_
#define _NCM_SPECTRAL_PANELS_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble
ncm_spectral_panels_eval_row (NcmSpectralPanels *sp, guint r, gdouble x)
{
  if (sp->uniform)
  {
    const gdouble u  = (x - sp->uniform_x0) * sp->uniform_inv_h;
    const gdouble u0 = (u > 0.0) ? u : 0.0;
    const gint i     = (gint) ((u0 < sp->uniform_last) ? u0 : sp->uniform_last);
    const guint N    = sp->uniform_N;
    const gsize p    = ((gsize) i * sp->fdim + r) * N;
    const gdouble s  = 2.0 * (u - i) - 1.0;

    /* The power-basis copy for N <= NCM_SPECTRAL_PANELS_POWER_N_MAX, else the Chebyshev series */
    switch (N)
    {
      case 1:
        return sp->uniform_power[p];

      case 2:
        return _ncm_spectral_panels_horner_n (sp->uniform_power + p, 2, s);

      case 3:
        return _ncm_spectral_panels_horner_n (sp->uniform_power + p, 3, s);

      case 4:
        return _ncm_spectral_panels_horner_n (sp->uniform_power + p, 4, s);

      case 5:
        return _ncm_spectral_panels_horner_n (sp->uniform_power + p, 5, s);

      case 6:
        return _ncm_spectral_panels_horner_n (sp->uniform_power + p, 6, s);

      case 7:
        return _ncm_spectral_panels_horner_n (sp->uniform_power + p, 7, s);

      case 8:
        return _ncm_spectral_panels_horner_n (sp->uniform_power + p, 8, s);

      case 9:
        return _ncm_spectral_panels_horner_n (sp->uniform_power + p, 9, s);

      default:
        return _ncm_spectral_panels_clenshaw_n (sp->panel_data[i] + (gsize) r * N, N, s);
    }
  }

  return ncm_spectral_panels_eval_row_on (sp, ncm_spectral_panels_find_panel (sp, x), r, x);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_SPECTRAL_PANELS_INLINE_H_ */

