/***************************************************************************
 *            ncm_hoaa.h
 *
 *  Fri November 04 13:27:40 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_hoaa.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_HOAA_H_
#define _NCM_HOAA_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>

G_BEGIN_DECLS

#define NCM_TYPE_HOAA             (ncm_hoaa_get_type ())
#define NCM_HOAA(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_HOAA, NcmHOAA))
#define NCM_HOAA_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_HOAA, NcmHOAAClass))
#define NCM_IS_HOAA(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_HOAA))
#define NCM_IS_HOAA_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_HOAA))
#define NCM_HOAA_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_HOAA, NcmHOAAClass))

typedef struct _NcmHOAAClass NcmHOAAClass;
typedef struct _NcmHOAA NcmHOAA;
typedef struct _NcmHOAAPrivate NcmHOAAPrivate;

/**
 * NcmHOAAOpt:
 * @NCM_HOAA_OPT_FULL: FIXME
 * @NCM_HOAA_OPT_V_ONLY: FIXME
 * @NCM_HOAA_OPT_DLNMNU_ONLY: FIXME
 * @NCM_HOAA_OPT_INVALID: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcmHOAAOpt
{
  NCM_HOAA_OPT_FULL = 0,
  NCM_HOAA_OPT_V_ONLY,
  NCM_HOAA_OPT_DLNMNU_ONLY,
  NCM_HOAA_OPT_INVALID,
} NcmHOAAOpt;

/**
 * NcmHOAASingType:
 * @NCM_HOAA_SING_TYPE_ZERO: FIXME
 * @NCM_HOAA_SING_TYPE_INF: FIXME
 * @NCM_HOAA_SING_TYPE_INVALID: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcmHOAASingType
{
  NCM_HOAA_SING_TYPE_ZERO = 0,
  NCM_HOAA_SING_TYPE_INF,
  NCM_HOAA_SING_TYPE_INVALID,
} NcmHOAASingType;

struct _NcmHOAAClass
{
  /*< private >*/
  GObjectClass parent_class;
  gdouble (*eval_mnu) (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_nu) (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_dlnmnu) (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
  gdouble (*eval_V) (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
  void (*eval_system) (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu);
  guint (*nsing) (NcmHOAA *hoaa, NcmModel *model, const gdouble k);
  void (*get_sing_info) (NcmHOAA *hoaa, NcmModel *model, const gdouble k, const guint sing, gdouble *ts, gdouble *dts_i, gdouble *dts_f, NcmHOAASingType *st);
  gdouble (*eval_sing_mnu) (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing);
  gdouble (*eval_sing_dlnmnu) (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing);
  gdouble (*eval_sing_V) (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing);
  void (*eval_sing_system) (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu);
	gdouble (*eval_powspec_factor) (NcmHOAA *hoaa, NcmModel *model);
  void (*prepare) (NcmHOAA *hoaa, NcmModel *model);
};

struct _NcmHOAA
{
  /*< private >*/
  GObject parent_instance;
  NcmHOAAPrivate *priv;
  gdouble k;
};

/**
 * NcmHOAAVar:
 * @NCM_HOAA_VAR_QBAR: FIXME
 * @NCM_HOAA_VAR_PBAR: FIXME
 * @NCM_HOAA_VAR_UPSILON: FIXME
 * @NCM_HOAA_VAR_GAMMA: FIXME
 * @NCM_HOAA_VAR_SYS_SIZE: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcmHOAAVar
{
  NCM_HOAA_VAR_QBAR = 0,
  NCM_HOAA_VAR_PBAR,
  NCM_HOAA_VAR_UPSILON,
  NCM_HOAA_VAR_GAMMA,
  NCM_HOAA_VAR_SYS_SIZE,  
} NcmHOAAVar;

GType ncm_hoaa_get_type (void) G_GNUC_CONST;

NcmHOAA *ncm_hoaa_ref (NcmHOAA *hoaa);
void ncm_hoaa_free (NcmHOAA *hoaa);
void ncm_hoaa_clear (NcmHOAA **hoaa);

void ncm_hoaa_set_reltol (NcmHOAA *hoaa, const gdouble reltol);
void ncm_hoaa_set_abstol (NcmHOAA *hoaa, const gdouble abstol);
void ncm_hoaa_set_k (NcmHOAA *hoaa, const gdouble k);
void ncm_hoaa_set_ti (NcmHOAA *hoaa, const gdouble ti);
void ncm_hoaa_set_tf (NcmHOAA *hoaa, const gdouble tf);

void ncm_hoaa_save_evol (NcmHOAA *hoaa, gboolean save_evol);
void ncm_hoaa_prepare (NcmHOAA *hoaa, NcmModel *model);

void ncm_hoaa_get_t0_t1 (NcmHOAA *hoaa, NcmModel *model, gdouble *t0, gdouble *t1);

void ncm_hoaa_eval_adiabatic_approx (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *thetab, gdouble *upsilon, gdouble *gamma);
void ncm_hoaa_eval_adiabatic_LnI_approx (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble theta, const gdouble psi, gdouble *LnI, gdouble *LnJ);
void ncm_hoaa_eval_AA (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *upsilon, gdouble *gamma, gdouble *qbar, gdouble *pbar);
void ncm_hoaa_eval_AA2QV (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble upsilon, const gdouble gamma, const gdouble qbar, const gdouble pbar, gdouble *q, gdouble *v, gdouble *Pq, gdouble *Pv);
void ncm_hoaa_eval_QV2AA (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble q, const gdouble v, const gdouble Pq, const gdouble Pv, gdouble *upsilon, gdouble *gamma, gdouble *qbar, gdouble *pbar);
void ncm_hoaa_eval_QV (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *q, gdouble *v, gdouble *Pq, gdouble *Pv);
void ncm_hoaa_eval_Delta (NcmHOAA *hoaa, NcmModel *model, const gdouble t, gdouble *Delta_phi, gdouble *Delta_Pphi);
void ncm_hoaa_eval_solution (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble S, const gdouble PS, gdouble *Aq, gdouble *Av);

NCM_INLINE gdouble ncm_hoaa_eval_nu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_hoaa_eval_mnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_hoaa_eval_dlnmnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE gdouble ncm_hoaa_eval_V (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
NCM_INLINE void ncm_hoaa_eval_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu);

NCM_INLINE guint ncm_hoaa_nsing (NcmHOAA *hoaa, NcmModel *model, const gdouble k);
NCM_INLINE void ncm_hoaa_get_sing_info (NcmHOAA *hoaa, NcmModel *model, const gdouble k, const guint sing, gdouble *ts, gdouble *dts_i, gdouble *dts_f, NcmHOAASingType *st);
NCM_INLINE gdouble ncm_hoaa_eval_sing_mnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing);
NCM_INLINE gdouble ncm_hoaa_eval_sing_dlnmnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing);
NCM_INLINE gdouble ncm_hoaa_eval_sing_V (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing);
NCM_INLINE void ncm_hoaa_eval_sing_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu);

NCM_INLINE gdouble ncm_hoaa_eval_powspec_factor (NcmHOAA *hoaa, NcmModel *model);

#define NCM_HOAA_TIME_FRAC (1.0e-13)
#define NCM_HOAA_DEBUG_EVOL (FALSE)
#define NCM_HOAA_DEBUG_SING (FALSE)
#define NCM_HOAA_DEBUG_EVOL_SING (FALSE)
#define NCM_HOAA_PARABOLIC_MIN_POINTS (4)
#define NCM_HOAA_PARABOLIC_TRIG_ONE (0.999)

G_END_DECLS

#endif /* _NCM_HOAA_H_ */

#ifndef _NCM_HOAA_INLINE_H_
#define _NCM_HOAA_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble 
ncm_hoaa_eval_nu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_HOAA_GET_CLASS (hoaa)->eval_nu (hoaa, model, t, k);
}

NCM_INLINE gdouble 
ncm_hoaa_eval_mnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_HOAA_GET_CLASS (hoaa)->eval_mnu (hoaa, model, t, k);
}

NCM_INLINE gdouble 
ncm_hoaa_eval_dlnmnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_HOAA_GET_CLASS (hoaa)->eval_dlnmnu (hoaa, model, t, k);
}

NCM_INLINE gdouble 
ncm_hoaa_eval_V (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k)
{
  return NCM_HOAA_GET_CLASS (hoaa)->eval_V (hoaa, model, t, k);
}

NCM_INLINE void 
ncm_hoaa_eval_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu)
{
  NCM_HOAA_GET_CLASS (hoaa)->eval_system (hoaa, model, t, k, nu, dlnmnu, Vnu);
}

NCM_INLINE guint 
ncm_hoaa_nsing (NcmHOAA *hoaa, NcmModel *model, const gdouble k)
{
  return NCM_HOAA_GET_CLASS (hoaa)->nsing (hoaa, model, k);
}

NCM_INLINE void 
ncm_hoaa_get_sing_info (NcmHOAA *hoaa, NcmModel *model, const gdouble k, const guint sing, gdouble *ts, gdouble *dts_i, gdouble *dts_f, NcmHOAASingType *st)
{
  return NCM_HOAA_GET_CLASS (hoaa)->get_sing_info (hoaa, model, k, sing, ts, dts_i, dts_f, st);
}

NCM_INLINE gdouble 
ncm_hoaa_eval_sing_mnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing)
{
  return NCM_HOAA_GET_CLASS (hoaa)->eval_sing_mnu (hoaa, model, t_m_ts, k, sing);
}

NCM_INLINE gdouble 
ncm_hoaa_eval_sing_dlnmnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing)
{
  return NCM_HOAA_GET_CLASS (hoaa)->eval_sing_dlnmnu (hoaa, model, t_m_ts, k, sing);
}

NCM_INLINE gdouble 
ncm_hoaa_eval_sing_V (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing)
{
  return NCM_HOAA_GET_CLASS (hoaa)->eval_sing_V (hoaa, model, t_m_ts, k, sing);
}

NCM_INLINE void 
ncm_hoaa_eval_sing_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu)
{
  NCM_HOAA_GET_CLASS (hoaa)->eval_sing_system (hoaa, model, t_m_ts, k, sing, nu, dlnmnu, Vnu);
}

NCM_INLINE gdouble 
ncm_hoaa_eval_powspec_factor (NcmHOAA *hoaa, NcmModel *model)
{
  return NCM_HOAA_GET_CLASS (hoaa)->eval_powspec_factor (hoaa, model);
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_HOAA_INLINE_H_ */
