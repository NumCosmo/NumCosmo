/***************************************************************************
 *            ncm_laurent_series.h
 *
 *  Tue Jul 8 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_laurent_series.h
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
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _NCM_LAURENT_SERIES_H_
#define _NCM_LAURENT_SERIES_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/ncm/core/ncm_util.h>
#include <numcosmo/ncm/algebra/ncm_complex.h>

G_BEGIN_DECLS

#define NCM_TYPE_LAURENT_SERIES (ncm_laurent_series_get_type ())

GType ncm_laurent_series_get_type (void) G_GNUC_CONST;

typedef struct _NcmLaurentSeries NcmLaurentSeries;

/**
 * NcmLaurentSeries:
 * @hmin: lowest power $h_\mathrm{min}$
 * @hmax: highest power $h_\mathrm{max}$
 * @c_cap: allocated length of @c, at least $h_\mathrm{max} - h_\mathrm{min} + 1$
 * @c: the coefficients $c_{h_\mathrm{min}}, \dots, c_{h_\mathrm{max}}$
 * @ref_count: the reference count
 *
 * A Laurent polynomial, see the class documentation in ncm_laurent_series.c.
 */
struct _NcmLaurentSeries
{
  gint hmin;
  gint hmax;
  gint c_cap;
  NcmComplex *c;
  gatomicrefcount ref_count;
};

NcmLaurentSeries *ncm_laurent_series_new_single (gint h, NcmComplex val);
NcmLaurentSeries *ncm_laurent_series_new (gint hmin, gint hmax);
NcmLaurentSeries *ncm_laurent_series_copy (const NcmLaurentSeries *a);
NcmLaurentSeries *ncm_laurent_series_ref (NcmLaurentSeries *a);
void ncm_laurent_series_free (NcmLaurentSeries *a);
void ncm_laurent_series_clear (NcmLaurentSeries **a);

void ncm_laurent_series_reset (NcmLaurentSeries *a, gint hmin, gint hmax);
gint ncm_laurent_series_get_hmin (const NcmLaurentSeries *a);
gint ncm_laurent_series_get_hmax (const NcmLaurentSeries *a);
void ncm_laurent_series_get_ptr (const NcmLaurentSeries *a, gint h, NcmComplex *out);
void ncm_laurent_series_set_ptr (NcmLaurentSeries *a, gint h, const NcmComplex *val);

NcmLaurentSeries *ncm_laurent_series_add (const NcmLaurentSeries *a, const NcmLaurentSeries *b, gdouble sb);
NcmLaurentSeries *ncm_laurent_series_scale_ptr (const NcmLaurentSeries *a, const NcmComplex *s);
NcmLaurentSeries *ncm_laurent_series_conv (const NcmLaurentSeries *a, const NcmLaurentSeries *b);
NcmLaurentSeries *ncm_laurent_series_conj (const NcmLaurentSeries *a);

void ncm_laurent_series_eval_ptr (const NcmLaurentSeries *a, const NcmComplex *w, NcmComplex *out);
gdouble ncm_laurent_series_jacobi_anger_reduce (const NcmLaurentSeries *cm, gdouble phi, const gdouble *Ik, gint n_Ik);
void ncm_laurent_series_jacobi_anger_accumulate (const NcmLaurentSeries *cm, const gdouble *Ik, gint n_Ik, gdouble scale, NcmComplex *H);
gdouble ncm_laurent_series_jacobi_anger_eval (const NcmComplex *H, gint n_H, gdouble phi);
NcmComplex ncm_laurent_series_get (const NcmLaurentSeries *a, gint h);

void ncm_laurent_series_set (NcmLaurentSeries *a, gint h, NcmComplex val);
NcmLaurentSeries *ncm_laurent_series_scale (const NcmLaurentSeries *a, NcmComplex s);
NcmComplex ncm_laurent_series_eval (const NcmLaurentSeries *a, NcmComplex w);

void ncm_laurent_series_set_single_into (NcmLaurentSeries *out, gint h, NcmComplex val);
void ncm_laurent_series_add_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, const NcmLaurentSeries *b, gdouble sb);
void ncm_laurent_series_conv_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, const NcmLaurentSeries *b);
void ncm_laurent_series_scale_into (NcmLaurentSeries *out, const NcmLaurentSeries *a, NcmComplex s);
void ncm_laurent_series_conj_into (NcmLaurentSeries *out, const NcmLaurentSeries *a);

typedef struct _NcmLaurentSeriesTPS NcmLaurentSeriesTPS;

#define NCM_TYPE_LAURENT_SERIES_TPS (ncm_laurent_series_tps_get_type ())

GType ncm_laurent_series_tps_get_type (void) G_GNUC_CONST;

NcmLaurentSeriesTPS *ncm_laurent_series_tps_new (guint order);
NcmLaurentSeriesTPS *ncm_laurent_series_tps_ref (NcmLaurentSeriesTPS *tps);
void ncm_laurent_series_tps_unref (NcmLaurentSeriesTPS *tps);
void ncm_laurent_series_tps_clear (NcmLaurentSeriesTPS **tps);

guint ncm_laurent_series_tps_order (const NcmLaurentSeriesTPS *tps);
NcmLaurentSeries *ncm_laurent_series_tps_get (const NcmLaurentSeriesTPS *tps, guint n);

void ncm_laurent_series_tps_eval_ptr (const NcmLaurentSeriesTPS *tps, const NcmComplex *w, const NcmComplex *g, NcmComplex *out);

void ncm_laurent_series_tps_pow (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, gdouble p);

void ncm_laurent_series_tps_conv (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, const NcmLaurentSeriesTPS *b);
void ncm_laurent_series_tps_conj (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a);
void ncm_laurent_series_tps_add (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, const NcmLaurentSeriesTPS *b, gdouble sb);

void ncm_laurent_series_tps_scale (NcmLaurentSeriesTPS *out, const NcmLaurentSeriesTPS *a, NcmComplex s);
NcmComplex ncm_laurent_series_tps_eval (const NcmLaurentSeriesTPS *tps, NcmComplex w, NcmComplex g);

G_END_DECLS

#endif /* _NCM_LAURENT_SERIES_H_ */

