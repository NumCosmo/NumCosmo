/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_qm_prop.h
 *
 *  Thu February 15 14:45:15 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_qm_prop.h
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

#ifndef _NCM_QM_PROP_H_
#define _NCM_QM_PROP_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_spline.h>

G_BEGIN_DECLS

#define NCM_TYPE_QM_PROP             (ncm_qm_prop_get_type ())
#define NCM_QM_PROP(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_QM_PROP, NcmQMProp))
#define NCM_QM_PROP_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_QM_PROP, NcmQMPropClass))
#define NCM_IS_QM_PROP(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_QM_PROP))
#define NCM_IS_QM_PROP_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_QM_PROP))
#define NCM_QM_PROP_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_QM_PROP, NcmQMPropClass))

typedef struct _NcmQMPropClass NcmQMPropClass;
typedef struct _NcmQMProp NcmQMProp;
typedef struct _NcmQMPropPrivate NcmQMPropPrivate;
typedef struct _NcmQMPropGauss NcmQMPropGauss;
typedef struct _NcmQMPropExp NcmQMPropExp;

struct _NcmQMPropClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmQMProp
{
  /*< private >*/
  GObject parent_instance;
  NcmQMPropPrivate *priv;
};

/**
 * NcmQMPropGauss:
 * 
 * Gaussian wave-function.
 * 
 */
struct _NcmQMPropGauss
{
  /*< private >*/
  gdouble mean;
  gdouble alpha; 
  gdouble sigma;
  gdouble Hi;
  gdouble lnNorm;
};

/**
 * NcmQMPropExp:
 * 
 * Exponential wave-function.
 * 
 */
struct _NcmQMPropExp
{
  /*< private >*/
  gdouble n;
  gdouble V; 
  gdouble pV;
  gdouble lnNorm;
};

/**
 * NcmQMPropPsi:
 * @psi_data: object pointer
 * @x: eval point $x$
 * @psi: (out caller-allocates) (array fixed-size=2) (element-type gdouble): $\psi$
 * 
 * Wave-function
 * 
 */
typedef void (*NcmQMPropPsi) (gpointer psi_data, const gdouble x, gdouble *psi);

GType ncm_qm_prop_get_type (void) G_GNUC_CONST;
GType ncm_qm_prop_gauss_get_type (void) G_GNUC_CONST;
GType ncm_qm_prop_exp_get_type (void) G_GNUC_CONST;

NcmQMPropGauss *ncm_qm_prop_gauss_new (const gdouble mean, const gdouble alpha, const gdouble sigma, const gdouble Hi);
NcmQMPropGauss *ncm_qm_prop_gauss_dup (NcmQMPropGauss *qm_gauss);
void ncm_qm_prop_gauss_free (NcmQMPropGauss *qm_gauss);

void ncm_qm_prop_gauss_eval (NcmQMPropGauss *qm_gauss, const gdouble x, gdouble *psi);
void ncm_qm_prop_gauss_eval_hermit (NcmQMPropGauss *qm_gauss, const gdouble x, gdouble *psi);
void ncm_qm_prop_gauss_eval_lnRS (NcmQMPropGauss *qm_gauss, const gdouble x, gdouble *lnRS);

NcmQMPropExp *ncm_qm_prop_exp_new (const gdouble n, const gdouble V, const gdouble pV);
NcmQMPropExp *ncm_qm_prop_exp_dup (NcmQMPropExp *qm_exp);
void ncm_qm_prop_exp_free (NcmQMPropExp *qm_exp);

void ncm_qm_prop_exp_eval (NcmQMPropExp *qm_exp, const gdouble x, gdouble *psi);
void ncm_qm_prop_exp_eval_lnRS (NcmQMPropExp *qm_exp, const gdouble x, gdouble *lnRS);

NcmQMProp *ncm_qm_prop_new (void);
NcmQMProp *ncm_qm_prop_ref (NcmQMProp *qm_prop);

void ncm_qm_prop_free (NcmQMProp *qm_prop);
void ncm_qm_prop_clear (NcmQMProp **qm_prop);

void ncm_qm_prop_set_nknots (NcmQMProp *qm_prop, const guint nknots);
guint ncm_qm_prop_get_nknots (NcmQMProp *qm_prop);

void ncm_qm_prop_eval (NcmQMProp *qm_prop, const gdouble x, const gdouble y, const gdouble t, gdouble *G);
void ncm_qm_prop_eval_array (NcmQMProp *qm_prop, const gdouble x, const gdouble *ya, gsize n, const gdouble t, gdouble *G);

void ncm_qm_prop_gauss_ini (NcmQMProp *qm_prop, const gdouble mean, const gdouble alpha, const gdouble sigma, const gdouble Hi);
void ncm_qm_prop_propto (NcmQMProp *qm_prop, const gdouble x, const gdouble t, gdouble *psi);

gdouble ncm_qm_prop_propto_norm (NcmQMProp *qm_prop, const gdouble t);

void ncm_qm_prop_set_init_cond (NcmQMProp *qm_prop, NcmQMPropPsi psi0_lnRS, gpointer psi_data, const gdouble xi, const gdouble xf);
void ncm_qm_prop_set_init_cond_gauss (NcmQMProp *qm_prop, NcmQMPropGauss *qm_gauss, const gdouble xi, const gdouble xf);
void ncm_qm_prop_set_init_cond_exp (NcmQMProp *qm_prop, NcmQMPropExp *qm_exp, const gdouble xi, const gdouble xf);

void ncm_qm_prop_evolve (NcmQMProp *qm_prop, const gdouble tf);
void ncm_qm_prop_evolve_spec (NcmQMProp *qm_prop, const gdouble t);

GArray *ncm_qm_prop_eval_psi (NcmQMProp *qm_prop, const gdouble *x, const guint len);
GArray *ncm_qm_prop_eval_rho (NcmQMProp *qm_prop, const gdouble *x, const guint len);
GArray *ncm_qm_prop_eval_dS (NcmQMProp *qm_prop, const gdouble *x, const guint len);

NcmSpline *ncm_qm_prop_peek_rho_s (NcmQMProp *qm_prop);

G_END_DECLS

#endif /* _NCM_QM_PROP_H_ */
