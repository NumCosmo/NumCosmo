/***************************************************************************
 *            ncm_fftlog.h
 *
 *  Fri May 18 16:44:28 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#ifndef _NCM_FFTLOG_H_
#define _NCM_FFTLOG_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_spline.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NCM_TYPE_FFTLOG (ncm_fftlog_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmFftlog, ncm_fftlog, NCM, FFTLOG, GObject)

struct _NcmFftlogClass
{
  /*< private >*/
  GObjectClass parent_class;
  const gchar *name;

  void (*compute_Ym) (NcmFftlog *fftlog, gpointer Ym_0);
  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[16];
};

typedef gdouble (*NcmFftlogFunc) (const gdouble x, gpointer user_data);

NcmFftlog *ncm_fftlog_ref (NcmFftlog *fftlog);

void ncm_fftlog_free (NcmFftlog *fftlog);
void ncm_fftlog_clear (NcmFftlog **fftlog);

const gchar *ncm_fftlog_peek_name (NcmFftlog *fftlog);

void ncm_fftlog_reset (NcmFftlog *fftlog);

void ncm_fftlog_set_nderivs (NcmFftlog *fftlog, guint nderivs);
guint ncm_fftlog_get_nderivs (NcmFftlog *fftlog);

void ncm_fftlog_set_lnr0 (NcmFftlog *fftlog, const gdouble lnr0);
gdouble ncm_fftlog_get_lnr0 (NcmFftlog *fftlog);

void ncm_fftlog_set_lnk0 (NcmFftlog *fftlog, const gdouble lnk0);
gdouble ncm_fftlog_get_lnk0 (NcmFftlog *fftlog);

void ncm_fftlog_set_size (NcmFftlog *fftlog, guint n);

void ncm_fftlog_set_padding (NcmFftlog *fftlog, gdouble pad_p);
gdouble ncm_fftlog_get_padding (NcmFftlog *fftlog);

void ncm_fftlog_set_noring (NcmFftlog *fftlog, gboolean active);
gboolean ncm_fftlog_get_noring (NcmFftlog *fftlog);

void ncm_fftlog_set_length (NcmFftlog *fftlog, gdouble Lk);

void ncm_fftlog_use_eval_interval (NcmFftlog *fftlog, gboolean use_eval_interal);
void ncm_fftlog_use_smooth_padding (NcmFftlog *fftlog, gboolean use_smooth_padding);
void ncm_fftlog_set_smooth_padding_scale (NcmFftlog *fftlog, gdouble log10sc);
gdouble ncm_fftlog_get_smooth_padding_scale (NcmFftlog *fftlog);

void ncm_fftlog_set_eval_r_min (NcmFftlog *fftlog, const gdouble eval_r_min);
void ncm_fftlog_set_eval_r_max (NcmFftlog *fftlog, const gdouble eval_r_max);
gdouble ncm_fftlog_get_eval_r_min (NcmFftlog *fftlog);
gdouble ncm_fftlog_get_eval_r_max (NcmFftlog *fftlog);

gdouble *ncm_fftlog_get_Ym (NcmFftlog *fftlog, guint *size);

void ncm_fftlog_get_lnk_vector (NcmFftlog *fftlog, NcmVector *lnk);
void ncm_fftlog_eval_by_vector (NcmFftlog *fftlog, NcmVector *Fk);
void ncm_fftlog_eval_by_function (NcmFftlog *fftlog, NcmFftlogFunc Fk, gpointer user_data);
void ncm_fftlog_eval_by_gsl_function (NcmFftlog *fftlog, gsl_function *Fk);

void ncm_fftlog_prepare_splines (NcmFftlog *fftlog);

NcmVector *ncm_fftlog_get_vector_lnr (NcmFftlog *fftlog);
NcmVector *ncm_fftlog_get_vector_Gr (NcmFftlog *fftlog, guint nderiv);

NcmSpline *ncm_fftlog_peek_spline_Gr (NcmFftlog *fftlog, guint nderiv);

gdouble ncm_fftlog_eval_output (NcmFftlog *fftlog, guint nderiv, const gdouble lnr);

void ncm_fftlog_calibrate_size (NcmFftlog *fftlog, NcmFftlogFunc Fk, gpointer user_data, const gdouble reltol);
void ncm_fftlog_calibrate_size_gsl (NcmFftlog *fftlog, gsl_function *Fk, const gdouble reltol);

guint ncm_fftlog_get_size (NcmFftlog *fftlog);
gint ncm_fftlog_get_full_size (NcmFftlog *fftlog);
gdouble ncm_fftlog_get_norma (NcmFftlog *fftlog);
gdouble ncm_fftlog_get_length (NcmFftlog *fftlog);
gdouble ncm_fftlog_get_full_length (NcmFftlog *fftlog);

gint ncm_fftlog_get_mode_index (NcmFftlog *fftlog, gint i);
gint ncm_fftlog_get_array_index (NcmFftlog *fftlog, gint phys_i);

NcmVector *ncm_fftlog_peek_output_vector (NcmFftlog *fftlog, guint nderiv);

G_END_DECLS

#endif /* _NCM_FFTLOG_H_ */

