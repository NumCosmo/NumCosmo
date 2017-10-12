/***************************************************************************
 *            quantum-gravity.h
 *
 *  Wed Jan 28 11:23:51 2009
 *  Copyright  2009  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NC_QUANTUM_GRAVITY_H
#define _NC_QUANTUM_GRAVITY_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/math/ncm_spline.h>
#include <numcosmo/nc_hicosmo.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#endif

G_BEGIN_DECLS

NcHICosmo *nc_hicosmo_qg_new (void);
void nc_hicosmo_qg_max_z (NcmModel *model, gdouble *max, gdouble *trans);
gdouble nc_hicosmo_qg_get_eta_b (NcmModel *model, gpointer userdata);

gdouble nc_hicosmo_qg_gbar2 (NcmModel *model, gdouble x);
gdouble nc_hicosmo_qg_gbarbar (NcmModel *model, gdouble x);
gdouble nc_hicosmo_qg_xbar (NcmModel *model, gdouble x, gpointer userdata);
gdouble nc_hicosmo_qg_xxbarzeta2 (NcmModel *model, gdouble x, gpointer userdata);
gdouble nc_hicosmo_qg_dxxbarzeta2_xxbarzeta2 (NcmModel *model, gdouble x, gpointer userdata);

gdouble nc_hicosmo_qg_d2sqrtxxbarzeta_sqrtxxbarzeta (NcmModel *model, gdouble x, gpointer userdata);
gdouble nc_hicosmo_qg_cs2_xxbar2 (NcmModel *model, gdouble x, gpointer userdata);

gdouble nc_hicosmo_qg_cs2 (NcmModel *model, gdouble x, gpointer userdata);
gdouble nc_hicosmo_qg_dcs2 (NcmModel *model, gdouble x, gpointer userdata);
gdouble nc_hicosmo_qg_beta (NcmModel *model, gdouble x, gpointer userdata);
gdouble nc_hicosmo_qg_zeta (NcmModel *model, gdouble x, gpointer userdata);
gdouble nc_hicosmo_qg_dzeta_zeta (NcmModel *model, gdouble x, gpointer userdata);
gdouble nc_hicosmo_qg_ddzeta_zeta (NcmModel *model, gdouble x, gpointer userdata);
gdouble nc_hicosmo_qg_xddzeta_zeta_mxdzeta_zeta2_dzeta_zeta (NcmModel *model, gdouble x, gpointer userdata);
gdouble nc_hicosmo_qg_V (NcmModel *model, gdouble x, gpointer userdata);

gdouble nc_hicosmo_qg_alphaprime2 (NcmModel *model, gdouble alpha, gpointer data);
gdouble nc_hicosmo_qg_dalphaprime2_dalpha (NcmModel *model, gdouble alpha, gpointer data);

gdouble nc_hicosmo_qg_lambda_x (NcmModel *model, gdouble x, gpointer userdata);
gboolean nc_hicosmo_qg_past_sol (NcmModel *model, gdouble k, gdouble lambda, gsl_matrix *sol);
void nc_hicosmo_qg_h_to_R_matrix (NcmModel *model, gdouble x, gsl_matrix *T);
void nc_hicosmo_qg_R_to_h_matrix (NcmModel *model, gdouble x, gsl_matrix *T);

gdouble nc_hicosmo_qg_get_lambda_f (NcmModel *model, gpointer userdata);
gdouble nc_hicosmo_qg_get_lambda_i (NcmModel *model, gpointer userdata);
gdouble nc_hicosmo_qg_get_lambda_d (NcmModel *model, gpointer userdata);

gdouble nc_hicosmo_qg_eta_lambda (NcmModel *model, gdouble lambda, gboolean deriv);
gdouble nc_hicosmo_qg_x_lambda (NcmModel *model, gdouble lambda, gboolean deriv);
gdouble nc_hicosmo_qg_gbar_lambda (NcmModel *model, gdouble lambda, gboolean deriv);
gdouble nc_hicosmo_qg_gbarbar_lambda (NcmModel *model, gdouble lambda, gboolean deriv);
gdouble nc_hicosmo_qg_int_1_zeta2_lambda (NcmModel *model, gdouble lambda, gboolean deriv);
gdouble nc_hicosmo_qg_cs2zeta2_int_1_zeta2_lambda (NcmModel *model, gdouble lambda, gboolean deriv);

gdouble nc_hicosmo_qg_lambda_k_cross (NcmModel *model, gdouble lambda, gboolean deriv);

gdouble nc_hicosmo_qg_cs2_lambda (NcmModel *model, gdouble lambda, gboolean deriv);
gdouble nc_hicosmo_qg_V_lambda (NcmModel *model, gdouble lambda, gboolean deriv);

/**
 * NcHICosmoQGPertType:
 * @NC_HICOSMO_QG_PERT_CURVATURE: FIXME
 * @NC_HICOSMO_QG_PERT_H: FIXME
 * 
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_QG_PERT_TYPE >*/
{
  NC_HICOSMO_QG_PERT_CURVATURE = 0,
  NC_HICOSMO_QG_PERT_H,  
} NcHICosmoQGPertType;

typedef struct _NcHICosmoQGMode NcHICosmoQGMode;

/**
 * NcHICosmoQGMode:
 * 
 * FIXME
 */
struct _NcHICosmoQGMode
{
  /*< private >*/
  NcmModel *model;
  NcHICosmoQGPertType type;
  gpointer cvode;
  gpointer cvode_R;
  gpointer cvode_h;
  gboolean init_R;
  gboolean init_h;
  gboolean *initialized;
  CVDlsDenseJacFn jac;
  CVRhsFn f;
  N_Vector y;
  N_Vector yQ;
  gdouble reltol;
  gdouble abstol;
  gdouble ax_i;
  gdouble x_i;
  gdouble t_i;
  gdouble t;
  long double k;
  long double alpha0;
  long double alphai;
  long double alphaf;
  NcmSpline *pw_spline;
};

NcHICosmoQGMode *nc_hicosmo_qg_pert_new (NcmModel *model, gdouble ax_i, gdouble x_i);
gboolean nc_hicosmo_qg_pert_set_opts (NcHICosmoQGMode *qgmode);
gboolean nc_hicosmo_qg_pert_init (NcHICosmoQGMode *qgmode, gdouble k);
gboolean nc_hicosmo_qg_pert_evolve (NcHICosmoQGMode *qgmode);
gboolean nc_hicosmo_qg_pert_prepare_pw_spline (NcHICosmoQGMode *qgmode, gboolean verbose);
gdouble nc_hicosmo_qg_pert_powerspectrum (NcHICosmoQGMode *qgmode, gdouble x, gdouble *R);

NcHICosmoQGMode *nc_hicosmo_qg_modefunc (NcmModel *model, long double k, long double x0, long double xf);
gboolean nc_hicosmo_qg_modefunc_set_opts (NcHICosmoQGMode *qgmode);
gboolean nc_hicosmo_qg_modefunc_init (NcHICosmoQGMode *qgmode);
gboolean nc_hicosmo_qg_modefuncm_cvode_init (NcHICosmoQGMode *qgmode);
gboolean nc_hicosmo_qg_modefunc_evolve (NcHICosmoQGMode *qgmode);
void nc_hicosmo_qg_evolfunc (NcmModel *model, long double x, long double *x2d2sqrtxxbarzeta_sqrtxxbarzeta, long double *x2cs2_xxbar2);
void nc_hicosmo_qg_modefunc_sol (NcHICosmoQGMode *qgmode, long double x, long double x0, long double *Re_u, long double *Im_u, long double *Re_up, long double *Im_up);
void nc_hicosmo_qg_pert_R_to_h (NcHICosmoQGMode *qgmode, gdouble x, gdouble *R);
void nc_hicosmo_qg_pert_h_to_R (NcHICosmoQGMode *qgmode, gdouble x, gdouble *h);

G_END_DECLS

#endif /* NC_QUANTUM_GRAVITY_H */
