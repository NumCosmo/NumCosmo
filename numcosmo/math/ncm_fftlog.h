/***************************************************************************
 *            ncm_fftlog.h
 *
 *  Fri May 18 16:44:28 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
#include <gsl/gsl_math.h>
#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */
#ifndef NUMCOSMO_GIR_SCAN
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#endif

G_BEGIN_DECLS

#define NCM_TYPE_FFTLOG             (ncm_fftlog_get_type ())
#define NCM_FFTLOG(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_FFTLOG, NcmFftlog))
#define NCM_FFTLOG_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_FFTLOG, NcmFftlogClass))
#define NCM_IS_FFTLOG(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_FFTLOG))
#define NCM_IS_FFTLOG_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_FFTLOG))
#define NCM_FFTLOG_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_FFTLOG, NcmFftlogClass))

typedef struct _NcmFftlogClass NcmFftlogClass;
typedef struct _NcmFftlog NcmFftlog;

struct _NcmFftlogClass
{
  /*< private >*/
  GObjectClass parent_class;
  const gchar *name;
  void (*get_Ym) (NcmFftlog *fftlog, gpointer Ym_0);
};

struct _NcmFftlog
{
  /*< private >*/
  GObject parent_instance;
  gint Nr;
  gint N;
  gint N_2;
  gint Nf;
  gint Nf_2;
  guint nderivs;
  guint pad;
  gdouble lnk0;
  gdouble lnr0;
  gdouble Lk;
  gdouble Lk_N;
  gdouble pad_p;
  gboolean prepared;
  gboolean evaluated;
  NcmVector *lnr_vec;
  GPtrArray *Gr_vec;
  GPtrArray *Gr_s;
#ifdef NUMCOSMO_HAVE_FFTW3
  fftw_complex *Fk;
  fftw_complex *Cm;
  fftw_complex *Gr;
  fftw_complex *CmYm;
  GPtrArray *Ym;
  fftw_plan p_Fk2Cm;
  fftw_plan p_CmYm2Gr;
#endif /* NUMCOSMO_HAVE_FFTW3 */
};

GType ncm_fftlog_get_type (void) G_GNUC_CONST;

NcmFftlog *ncm_fftlog_ref (NcmFftlog *fftlog);
void ncm_fftlog_free (NcmFftlog *fftlog);
void ncm_fftlog_clear (NcmFftlog **fftlog);

void ncm_fftlog_set_name (NcmFftlog *fftlog, const gchar *name);
const gchar *ncm_fftlog_peek_name (NcmFftlog *fftlog);

void ncm_fftlog_set_nderivs (NcmFftlog *fftlog, guint nderivs);
guint ncm_fftlog_get_nderivs (NcmFftlog *fftlog);

void ncm_fftlog_set_lnr0 (NcmFftlog *fftlog, const gdouble lnr0);
gdouble ncm_fftlog_get_lnr0 (NcmFftlog *fftlog);

void ncm_fftlog_set_lnk0 (NcmFftlog *fftlog, const gdouble lnk0);
gdouble ncm_fftlog_get_lnk0 (NcmFftlog *fftlog);

void ncm_fftlog_set_size (NcmFftlog *fftlog, guint n);

void ncm_fftlog_set_padding (NcmFftlog *fftlog, gdouble pad_p);
gdouble ncm_fftlog_get_padding (NcmFftlog *fftlog);

void ncm_fftlog_set_length (NcmFftlog *fftlog, gdouble Lk);

void ncm_fftlog_eval_by_vector (NcmFftlog *fftlog, NcmVector *Fk);
void ncm_fftlog_eval_by_function (NcmFftlog *fftlog, gsl_function *Fk);

void ncm_fftlog_prepare_splines (NcmFftlog *fftlog);

NcmVector *ncm_fftlog_get_vector_lnr (NcmFftlog *fftlog);
NcmVector *ncm_fftlog_get_vector_Gr (NcmFftlog *fftlog, guint nderiv);

NcmSpline *ncm_fftlog_peek_spline_Gr (NcmFftlog *fftlog, guint nderiv);

gdouble ncm_fftlog_eval_output (NcmFftlog *fftlog, guint nderiv, const gdouble lnr);

void ncm_fftlog_calibrate_size (NcmFftlog *fftlog, gsl_function *Fk, gdouble reltol);

G_INLINE_FUNC guint ncm_fftlog_get_size (NcmFftlog *fftlog);
G_INLINE_FUNC gint ncm_fftlog_get_full_size (NcmFftlog *fftlog);
G_INLINE_FUNC gdouble ncm_fftlog_get_norma (NcmFftlog *fftlog);
G_INLINE_FUNC gdouble ncm_fftlog_get_length (NcmFftlog *fftlog);
G_INLINE_FUNC gdouble ncm_fftlog_get_full_length (NcmFftlog *fftlog);

G_INLINE_FUNC gint ncm_fftlog_get_mode_index (NcmFftlog *fftlog, gint i);

G_INLINE_FUNC NcmVector *ncm_fftlog_peek_output_vector (NcmFftlog *fftlog, guint nderiv);

G_END_DECLS

#endif /* _NCM_FFTLOG_H_ */

#ifndef _NC_FFTLOG_INLINE_H_
#define _NC_FFTLOG_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC guint
ncm_fftlog_get_size (NcmFftlog *fftlog)
{
  return fftlog->N;
}

G_INLINE_FUNC gint
ncm_fftlog_get_full_size (NcmFftlog *fftlog)
{
  return fftlog->Nf;
}

G_INLINE_FUNC gdouble
ncm_fftlog_get_norma (NcmFftlog *fftlog)
{
  return ncm_fftlog_get_full_size (fftlog);
}

G_INLINE_FUNC gdouble
ncm_fftlog_get_length (NcmFftlog *fftlog)
{
  return fftlog->Lk;
}

G_INLINE_FUNC gdouble 
ncm_fftlog_get_full_length (NcmFftlog *fftlog)
{
  return fftlog->Lk + 2.0 * fftlog->Lk_N * fftlog->pad;
}

G_INLINE_FUNC gint 
ncm_fftlog_get_mode_index (NcmFftlog *fftlog, gint i)
{
  return (i > fftlog->Nf_2) ? i - fftlog->Nf : i;
}

G_INLINE_FUNC NcmVector *
ncm_fftlog_peek_output_vector (NcmFftlog *fftlog, guint nderiv)
{
  return g_ptr_array_index (fftlog->Gr_vec, nderiv);
}

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NC_FFTLOG_INLINE_H_ */
