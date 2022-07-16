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
  
  gdouble (*eval_xi)         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_dxi)        (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_nu)         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_nu2)        (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_m)          (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_mnu)        (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_dlnmnu)     (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_V)          (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_system)     (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_int_1_m)    (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_int_mnu2)   (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_int_qmnu2)  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_int_q2mnu2) (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_dm)         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_F1)         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_F2)         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_FN)         (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t, const gdouble k);
  gdouble (*eval_powspec_factor) (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k);
  void (*prepare) (NcmCSQ1D *csq1d, NcmModel *model);
};

/**
 * NcmCSQ1DEvolState:
 * @NCM_CSQ1D_EVOL_STATE_INVALID: Invalid state
 * @NCM_CSQ1D_EVOL_STATE_ADIABATIC: Adiabatic state variables $(\alpha,\, \delta\gamma)$
 * @NCM_CSQ1D_EVOL_STATE_UP: $(\chi,\, U_+)$ state variables
 * @NCM_CSQ1D_EVOL_STATE_UM: $(\chi,\, U_-)$ state variables
 *
 * Variables describing the system state.
 */
typedef enum _NcmCSQ1DEvolState
{
  NCM_CSQ1D_EVOL_STATE_INVALID = 0,
  NCM_CSQ1D_EVOL_STATE_ADIABATIC,
  NCM_CSQ1D_EVOL_STATE_UP,
  NCM_CSQ1D_EVOL_STATE_UM,
} NcmCSQ1DEvolState;

struct _NcmCSQ1D
{
  /*< private >*/
  GObject parent_instance;
  NcmCSQ1DPrivate *priv;
};

GType ncm_csq1d_get_type (void) G_GNUC_CONST;

NcmCSQ1D *ncm_csq1d_ref (NcmCSQ1D *csq1d);
void ncm_csq1d_free (NcmCSQ1D *csq1d);
void ncm_csq1d_clear (NcmCSQ1D **csq1d);

void ncm_csq1d_set_reltol (NcmCSQ1D *csq1d, const gdouble reltol);
void ncm_csq1d_set_abstol (NcmCSQ1D *csq1d, const gdouble abstol);
void ncm_csq1d_set_k (NcmCSQ1D *csq1d, const gdouble k);
void ncm_csq1d_set_ti (NcmCSQ1D *csq1d, const gdouble ti);
void ncm_csq1d_set_tf (NcmCSQ1D *csq1d, const gdouble tf);
void ncm_csq1d_set_adiab_threshold (NcmCSQ1D *csq1d, const gdouble adiab_threshold);
void ncm_csq1d_set_prop_threshold (NcmCSQ1D *csq1d, const gdouble prop_threshold);
void ncm_csq1d_set_save_evol (NcmCSQ1D *csq1d, const gboolean save);
void ncm_csq1d_set_init_cond (NcmCSQ1D *csq1d, NcmCSQ1DEvolState state, const gdouble ti, const gdouble x, const gdouble y);
void ncm_csq1d_set_init_cond_adiab (NcmCSQ1D *csq1d, NcmModel *model, const gdouble ti);

gdouble ncm_csq1d_get_reltol (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_abstol (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_k (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_ti (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_tf (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_adiab_threshold (NcmCSQ1D *csq1d);
gdouble ncm_csq1d_get_prop_threshold (NcmCSQ1D *csq1d);
gboolean ncm_csq1d_get_save_evol (NcmCSQ1D *csq1d);

NCM_INLINE gdouble ncm_csq1d_eval_xi         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_dxi        (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_nu         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_nu2        (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_m          (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_mnu        (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_dlnmnu     (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_V          (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_system     (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_int_1_m    (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_int_mnu2   (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_int_qmnu2  (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_int_q2mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_dm         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_F1         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_F2         (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_FN         (NcmCSQ1D *csq1d, NcmModel *model, const gint n, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_csq1d_eval_powspec_factor (NcmCSQ1D *csq1d, NcmModel *model, const gdouble k);

void ncm_csq1d_prepare (NcmCSQ1D *csq1d, NcmModel *model);

GArray *ncm_csq1d_get_time_array (NcmCSQ1D *csq1d, gdouble *smallest_t);

gboolean ncm_csq1d_find_adiab_time_limit (NcmCSQ1D *csq1d, NcmModel *model, gdouble t0, gdouble t1, const gdouble reltol, gdouble *ti);
gdouble ncm_csq1d_find_adiab_max (NcmCSQ1D *csq1d, NcmModel *model, gdouble t0, gdouble t1, const gdouble border_eps, gdouble *F1_min, gdouble *t_Bl, gdouble *t_Bu);

void ncm_csq1d_eval_adiab_at (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *alpha, gdouble *dgamma, gdouble *alpha_reltol, gdouble *dgamma_reltol);
void ncm_csq1d_eval_nonadiab_at (NcmCSQ1D *csq1d, NcmModel *model, guint nonadiab_frame, const gdouble t, gdouble *chi, gdouble *Up);
void ncm_csq1d_eval_at (NcmCSQ1D *csq1d, const gdouble t, gdouble *alpha, gdouble *dgamma);

void ncm_csq1d_alpha_dgamma_to_phi_Pphi (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble alpha, const gdouble dgamma, gdouble *phi, gdouble *Pphi);
void ncm_csq1d_get_J_at (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *J11, gdouble *J12, gdouble *J22);

void ncm_csq1d_get_H_poincare_hp (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny);
void ncm_csq1d_get_Hadiab_poincare_hp (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny);
void ncm_csq1d_get_poincare_hp (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny);
void ncm_csq1d_get_poincare_hp_frame (NcmCSQ1D *csq1d, NcmModel *model, guint adiab_frame, const gdouble t, gdouble *x, gdouble *lny);

void ncm_csq1d_alpha_dgamma_to_minkowski_frame (NcmCSQ1D *csq1d, NcmModel *model, guint adiab_frame, const gdouble t, const gdouble alpha, const gdouble dgamma, gdouble *x1, gdouble *x2);
void ncm_csq1d_get_minkowski_frame (NcmCSQ1D *csq1d, NcmModel *model, guint adiab_frame, const gdouble t, gdouble *x1, gdouble *x2);

void ncm_csq1d_get_H_poincare_disc (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny);
void ncm_csq1d_get_Hadiab_poincare_disc (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny);
void ncm_csq1d_get_poincare_disc (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *x, gdouble *lny);

void ncm_csq1d_prepare_prop (NcmCSQ1D *csq1d, NcmModel *model, const gdouble ti, const gdouble tii, const gdouble tf);

gdouble ncm_csq1d_get_tf_prop (NcmCSQ1D *csq1d);
void ncm_csq1d_get_prop_vector_chi_Up (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, gdouble *chi, gdouble *Up);
void ncm_csq1d_evolve_prop_vector_chi_Up (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const guint nonadiab_frame, gdouble chi_i, gdouble Up_i, gdouble *chi, gdouble *Up);

void ncm_csq1d_alpha_gamma_circle (NcmCSQ1D *csq1d, NcmModel *model, const gdouble alpha, const gdouble gamma, const gdouble r, const gdouble theta, gdouble *alphap, gdouble *gammap);

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
ncm_csq1d_eval_mnu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_mnu (csq1d, model, t, k);
}

NCM_INLINE gdouble
ncm_csq1d_eval_dlnmnu (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_dlnmnu (csq1d, model, t, k);
}

NCM_INLINE gdouble
ncm_csq1d_eval_V (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_V (csq1d, model, t, k);
}

NCM_INLINE gdouble
ncm_csq1d_eval_system (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_system (csq1d, model, t, k);
}

NCM_INLINE gdouble
ncm_csq1d_eval_int_1_m (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_int_1_m (csq1d, model, t, k);
}

NCM_INLINE gdouble
ncm_csq1d_eval_int_mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_int_mnu2 (csq1d, model, t, k);
}

NCM_INLINE gdouble
ncm_csq1d_eval_int_qmnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_int_qmnu2 (csq1d, model, t, k);
}

NCM_INLINE gdouble
ncm_csq1d_eval_int_q2mnu2 (NcmCSQ1D *csq1d, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_CSQ1D_GET_CLASS (csq1d)->eval_int_q2mnu2 (csq1d, model, t, k);
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

