/***************************************************************************
 *            ncm_hoaa.h
 *
 *  Fri November 04 13:27:40 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_hoaa.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#define NCM_TYPE_HOAA (ncm_hoaa_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmHOAA, ncm_hoaa, NCM, HOAA, GObject)

/**
 * NcmHOAAOpt:
 * @NCM_HOAA_OPT_FULL: System with all variables non-zero
 * @NCM_HOAA_OPT_V_ONLY: System with only V non-zero
 * @NCM_HOAA_OPT_DLNMNU_ONLY: System with only dlnmnu non-zero
 * @NCM_HOAA_OPT_INVALID: Invalid option
 *
 * When solving the system of equations, the user can choose to solve the
 * system with all variables non-zero, only V non-zero or only dlnmnu non-zero.
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
 * @NCM_HOAA_SING_TYPE_ZERO: mass pass through zero
 * @NCM_HOAA_SING_TYPE_INF: mass pass through infinity
 * @NCM_HOAA_SING_TYPE_INVALID: Invalid option
 *
 * The singularity type.
 *
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

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[5];
};

/**
 * NcmHOAAVar:
 * @NCM_HOAA_VAR_QBAR: Field variable
 * @NCM_HOAA_VAR_PBAR: Momentum variable
 * @NCM_HOAA_VAR_UPSILON: adiabatic parameter
 * @NCM_HOAA_VAR_GAMMA: adiabatic parameter
 * @NCM_HOAA_VAR_SYS_SIZE: Number of variables in the system
 *
 * The variables in the system.
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

gdouble ncm_hoaa_eval_nu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
gdouble ncm_hoaa_eval_mnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
gdouble ncm_hoaa_eval_dlnmnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
gdouble ncm_hoaa_eval_V (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k);
void ncm_hoaa_eval_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t, const gdouble k, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu);

guint ncm_hoaa_nsing (NcmHOAA *hoaa, NcmModel *model, const gdouble k);
void ncm_hoaa_get_sing_info (NcmHOAA *hoaa, NcmModel *model, const gdouble k, const guint sing, gdouble *ts, gdouble *dts_i, gdouble *dts_f, NcmHOAASingType *st);
gdouble ncm_hoaa_eval_sing_mnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing);
gdouble ncm_hoaa_eval_sing_dlnmnu (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing);
gdouble ncm_hoaa_eval_sing_V (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing);
void ncm_hoaa_eval_sing_system (NcmHOAA *hoaa, NcmModel *model, const gdouble t_m_ts, const gdouble k, const guint sing, gdouble *nu, gdouble *dlnmnu, gdouble *Vnu);

gdouble ncm_hoaa_eval_powspec_factor (NcmHOAA *hoaa, NcmModel *model);

#define NCM_HOAA_TIME_FRAC (1.0e-13)
#define NCM_HOAA_DEBUG_EVOL (FALSE)
#define NCM_HOAA_DEBUG_SING (FALSE)
#define NCM_HOAA_DEBUG_EVOL_SING (FALSE)
#define NCM_HOAA_PARABOLIC_MIN_POINTS (4)
#define NCM_HOAA_PARABOLIC_TRIG_ONE (0.999)

G_END_DECLS

#endif /* _NCM_HOAA_H_ */

