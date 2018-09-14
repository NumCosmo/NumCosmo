/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hiqg_1d.h
 *
 *  Thu February 15 14:45:15 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hiqg_1d.h
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_HIQG_1D_H_
#define _NC_HIQG_1D_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>

G_BEGIN_DECLS

#define NC_TYPE_HIQG_1D             (nc_hiqg_1d_get_type ())
#define NC_HIQG_1D(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HIQG_1D, NcHIQG1D))
#define NC_HIQG_1D_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HIQG_1D, NcHIQG1DClass))
#define NC_IS_HIQG_1D(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HIQG_1D))
#define NC_IS_HIQG_1D_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HIQG_1D))
#define NC_HIQG_1D_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HIQG_1D, NcHIQG1DClass))

typedef struct _NcHIQG1DClass NcHIQG1DClass;
typedef struct _NcHIQG1D NcHIQG1D;
typedef struct _NcHIQG1DPrivate NcHIQG1DPrivate;
typedef struct _NcHIQG1DGauss NcHIQG1DGauss;
typedef struct _NcHIQG1DExp NcHIQG1DExp;

struct _NcHIQG1DClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcHIQG1D
{
  /*< private >*/
  GObject parent_instance;
  NcHIQG1DPrivate *priv;
};

/**
 * NcHIQG1DGauss:
 * 
 * Gaussian wave-function.
 * 
 */
struct _NcHIQG1DGauss
{
  /*< private >*/
  gdouble mean;
  gdouble alpha; 
  gdouble sigma;
  gdouble Hi;
  gdouble lnNorm;
};

/**
 * NcHIQG1DExp:
 * 
 * Exponential wave-function.
 * 
 */
struct _NcHIQG1DExp
{
  /*< private >*/
  gdouble n;
  gdouble V; 
  gdouble pV;
  gdouble lnNorm;
};

/**
 * NcHIQG1DPsi:
 * @psi_data: object pointer
 * @x: eval point $x$
 * @psi: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $\psi$
 * 
 * Wave-function
 * 
 */
typedef void (*NcHIQG1DPsi) (gpointer psi_data, const gdouble x, gdouble *psi);

GType nc_hiqg_1d_get_type (void) G_GNUC_CONST;
GType nc_hiqg_1d_gauss_get_type (void) G_GNUC_CONST;
GType nc_hiqg_1d_exp_get_type (void) G_GNUC_CONST;

NcHIQG1DGauss *nc_hiqg_1d_gauss_new (const gdouble mean, const gdouble alpha, const gdouble sigma, const gdouble Hi);
NcHIQG1DGauss *nc_hiqg_1d_gauss_dup (NcHIQG1DGauss *qm_gauss);
void nc_hiqg_1d_gauss_free (NcHIQG1DGauss *qm_gauss);

void nc_hiqg_1d_gauss_eval (NcHIQG1DGauss *qm_gauss, const gdouble x, gdouble *psi);
void nc_hiqg_1d_gauss_eval_hermit (NcHIQG1DGauss *qm_gauss, const gdouble x, gdouble *psi);
void nc_hiqg_1d_gauss_eval_lnRS (NcHIQG1DGauss *qm_gauss, const gdouble x, gdouble *lnRS);

NcHIQG1DExp *nc_hiqg_1d_exp_new (const gdouble n, const gdouble V, const gdouble pV);
NcHIQG1DExp *nc_hiqg_1d_exp_dup (NcHIQG1DExp *qm_exp);
void nc_hiqg_1d_exp_free (NcHIQG1DExp *qm_exp);

void nc_hiqg_1d_exp_eval (NcHIQG1DExp *qm_exp, const gdouble x, gdouble *psi);
void nc_hiqg_1d_exp_eval_lnRS (NcHIQG1DExp *qm_exp, const gdouble x, gdouble *lnRS);

NcHIQG1D *nc_hiqg_1d_new (void);
NcHIQG1D *nc_hiqg_1d_new_full (guint nknots, gdouble lambda);
NcHIQG1D *nc_hiqg_1d_ref (NcHIQG1D *qg1d);

void nc_hiqg_1d_free (NcHIQG1D *qg1d);
void nc_hiqg_1d_clear (NcHIQG1D **qg1d);

void nc_hiqg_1d_set_nknots (NcHIQG1D *qg1d, const guint nknots);
guint nc_hiqg_1d_get_nknots (NcHIQG1D *qg1d);

void nc_hiqg_1d_set_init_cond (NcHIQG1D *qg1d, NcHIQG1DPsi psi0_lnRS, gpointer psi_data, const gdouble xi, const gdouble xf);
void nc_hiqg_1d_set_init_cond_gauss (NcHIQG1D *qg1d, NcHIQG1DGauss *qm_gauss, const gdouble xi, const gdouble xf);
void nc_hiqg_1d_set_init_cond_exp (NcHIQG1D *qg1d, NcHIQG1DExp *qm_exp, const gdouble xi, const gdouble xf);

gdouble nc_hiqg_1d_basis (NcHIQG1D *qg1d, const gdouble x, const gdouble y, const gdouble h, const gdouble a);
gdouble nc_hiqg_1d_Hbasis (NcHIQG1D *qg1d, const gdouble x, const gdouble y, const gdouble h, const gdouble a);
gdouble nc_hiqg_1d_Sbasis_x3 (NcHIQG1D *qg1d, const gdouble x, const gdouble y1, const gdouble y2, const gdouble h, const gdouble a);

gdouble nc_hiqg_1d_get_lambda (NcHIQG1D *qg1d);
gdouble nc_hiqg_1d_get_basis_a (NcHIQG1D *qg1d);
gdouble nc_hiqg_1d_get_acs_a (NcHIQG1D *qg1d);
gdouble nc_hiqg_1d_get_nu (NcHIQG1D *qg1d);
gdouble nc_hiqg_1d_get_mu (NcHIQG1D *qg1d);

void nc_hiqg_1d_prepare (NcHIQG1D *qg1d);

NcmVector *nc_hiqg_1d_peek_knots (NcHIQG1D *qg1d);
gdouble nc_hiqg_1d_eval_ev (NcHIQG1D *qg1d, const gint i, const gdouble x);
void nc_hiqg_1d_eval_psi0 (NcHIQG1D *qg1d, const gdouble x, gdouble *psi0);
void nc_hiqg_1d_evol (NcHIQG1D *qg1d, const gdouble t);
void nc_hiqg_1d_eval_psi (NcHIQG1D *qg1d, const gdouble x, gdouble *psi);
gdouble nc_hiqg_1d_eval_dS (NcHIQG1D *qg1d, const gdouble x);
gdouble nc_hiqg_1d_int_rho_0_inf (NcHIQG1D *qg1d);
gdouble nc_hiqg_1d_int_xrho_0_inf (NcHIQG1D *qg1d);

gint nc_hiqg_1d_nBohm (NcHIQG1D *qg1d);
gdouble nc_hiqg_1d_Bohm (NcHIQG1D *qg1d, gint i);
gdouble nc_hiqg_1d_Bohm_p (NcHIQG1D *qg1d, gint i);

G_END_DECLS

#endif /* _NC_HIQG_1D_H_ */
