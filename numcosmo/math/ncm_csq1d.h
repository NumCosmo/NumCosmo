/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_csq1d.h
 *
 *  Mon September 09 13:56:11 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_csq1d.h
 * Copyright (C) 2019 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_CSQ1D_H_
#define _NCM_CSQ1D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>

G_BEGIN_DECLS

#define NCM_TYPE_CSQ1D             (ncm_csq1d_get_type ())
#define NCM_CSQ1D(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_CSQ1D, NcmCSQ1D))
#define NCM_CSQ1D_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_CSQ1D, NcmCSQ1DClass))
#define NCM_IS_CSQ1D(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_CSQ1D))
#define NCM_IS_CSQ1D_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_CSQ1D))
#define NCM_CSQ1D_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_CSQ1D, NcmCSQ1DClass))

typedef struct _NcmCSQ1DClass NcmCSQ1DClass;
typedef struct _NcmCSQ1D NcmCSQ1D;
typedef struct _NcmCSQ1DPrivate NcmCSQ1DPrivate;
typedef struct _NcmCSQ1DSingFitUp NcmCSQ1DSingFitUp;
typedef struct _NcmCSQ1DSingFitUm NcmCSQ1DSingFitUm;

struct _NcmCSQ1DClass
{
  /*< private >*/
  GObjectClass parent_class;
  gdouble (*eval_xi)  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_dxi) (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_nu)  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_nu2) (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_m)   (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_dm)  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_F1)  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_F2)  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_FN)  (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t, const gdouble k);
	gdouble (*eval_powspec_factor) (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k);
  void (*prepare) (NcmCSQ1D *csq1d, NcmModel *model);
};

struct _NcmCSQ1D
{
  /*< private >*/
  GObject parent_instance;
  NcmCSQ1DPrivate *priv;
};

/**
 * NcmCSQ1DSingFitUp:
 * 
 * Struct containig the $\Upsilon_+$ fitting model.
 * 
 */
struct _NcmCSQ1DSingFitUp
{
  /*< private >*/  
  gint chi_dim;
  gint Up_dim;
  NcmVector *chi_c;
  NcmVector *Up_c;
};

GType ncm_csq1d_sing_fit_up_get_type (void) G_GNUC_CONST;
GType ncm_csq1d_get_type (void) G_GNUC_CONST;

NcmCSQ1DSingFitUp *ncm_csq1d_sing_fit_up_new (const gint chi_dim, const gint Up_dim);
NcmCSQ1DSingFitUp *ncm_csq1d_sing_fit_up_dup (NcmCSQ1DSingFitUp *sing_up);
void ncm_csq1d_sing_fit_up_free (NcmCSQ1DSingFitUp *sing_up);

void ncm_csq1d_sing_fit_up_fit (NcmCSQ1DSingFitUp *sing_up, NcmCSQ1D *csq1d, NcmModel *model, NcmVector *t, NcmVector *chim_t, NcmVector *exp_Up);

gdouble ncm_csq1d_sing_fit_up_eval_chi (NcmCSQ1DSingFitUp *sing_up, NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
gdouble ncm_csq1d_sing_fit_up_eval_exp_Up (NcmCSQ1DSingFitUp *sing_up, NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);

gdouble ncm_csq1d_sing_fit_up_eval_dchi (NcmCSQ1DSingFitUp *sing_up, NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);
gdouble ncm_csq1d_sing_fit_up_eval_dexp_Up (NcmCSQ1DSingFitUp *sing_up, NcmCSQ1D *csq1d, NcmModel *model, const gdouble t);

NcmCSQ1D *ncm_csq1d_ref (NcmCSQ1D *csq1d);
void ncm_csq1d_free (NcmCSQ1D *csq1d);
void ncm_csq1d_clear (NcmCSQ1D **csq1d);

void ncm_csq1d_set_reltol (NcmCSQ1D *csq1d, const gdouble reltol);
void ncm_csq1d_set_abstol (NcmCSQ1D *csq1d, const gdouble abstol);
void ncm_csq1d_set_k (NcmCSQ1D *csq1d, const gdouble k);
void ncm_csq1d_set_ti (NcmCSQ1D *csq1d, const gdouble ti);
void ncm_csq1d_set_tf (NcmCSQ1D *csq1d, const gdouble tf);
void ncm_csq1d_set_adiab_threshold (NcmCSQ1D *csq1d, const gdouble adiab_threshold);
void ncm_csq1d_set_save_evol (NcmCSQ1D *csq1d, const gboolean save);
void ncm_csq1d_set_sing_detect (NcmCSQ1D *csq1d, const gboolean enable);

gdouble ncm_csq1d_get_reltol (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_abstol (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_k (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_ti (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_tf (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_adiab_threshold (NcmCSQ1D *csq1d);
gboolean ncm_csq1d_get_save_evol (NcmCSQ1D *csq1d);
gboolean ncm_csq1d_get_sing_detect (NcmCSQ1D *csq1d);

NCM_INLINE gdouble ncm_csq1d_eval_xi  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_dxi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_nu  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_nu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_m   (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_dm  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_F1  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_F2  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_FN  (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k);

void ncm_csq1d_prepare (NcmCSQ1D *csq1d, NcmModel *model);

GArray *ncm_csq1d_get_time_array (NcmCSQ1D *csq1d, gdouble *smallest_t);

gboolean ncm_csq1d_find_adiab_time_limit (NcmCSQ1D *csq1d, NcmModel *model, gdouble t0, gdouble t1, const gdouble reltol, gdouble *ti);

void ncm_csq1d_eval_adiab_at (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *alpha, gdouble *dgamma, gdouble *alpha_reltol, gdouble *dgamma_reltol);
void ncm_csq1d_eval_at (NcmCSQ1D *csq1d, const gdouble t, gdouble *alpha, gdouble *dgamma);

void ncm_csq1d_alpha_dgamma_to_phi_Pphi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble alpha, const gdouble dgamma, gdouble *phi, gdouble *Pphi);
void ncm_csq1d_get_J_at (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *J11, gdouble *J12, gdouble *J22);

G_END_DECLS

#endif /* _NCM_CSQ1D_H_ */
#ifndef _NCM_CSQ1D_INLINE_H_
#define _NCM_CSQ1D_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble 
ncm_csq1d_eval_xi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_xi (csq1d, model, t, k);
}

NCM_INLINE gdouble 
ncm_csq1d_eval_dxi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_dxi (csq1d, model, t, k);
}

NCM_INLINE gdouble 
ncm_csq1d_eval_nu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_nu (csq1d, model, t, k);
}

NCM_INLINE gdouble 
ncm_csq1d_eval_nu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_nu2 (csq1d, model, t, k);
}

NCM_INLINE gdouble 
ncm_csq1d_eval_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_m (csq1d, model, t, k);
}

NCM_INLINE gdouble 
ncm_csq1d_eval_dm (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_dm (csq1d, model, t, k);
}

NCM_INLINE gdouble 
ncm_csq1d_eval_F1 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_F1 (csq1d, model, t, k);
}

NCM_INLINE gdouble 
ncm_csq1d_eval_F2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_F2 (csq1d, model, t, k);
}

NCM_INLINE gdouble 
ncm_csq1d_eval_FN (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_FN (csq1d, model, n, t, k);
}

NCM_INLINE gdouble 
ncm_csq1d_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_powspec_factor (csq1d, model, k);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_CSQ1D_INLINE_H_ */
