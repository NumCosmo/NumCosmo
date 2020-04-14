/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hicosmo_qrbf.h
 *
 *  Fri November 01 14:17:49 2019
 *  Copyright  2019  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hicosmo_qrbf.h
 *
 * Copyright (C) 2019 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef _NC_HICOSMO_QRBF_H_
#define _NC_HICOSMO_QRBF_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_fit.h>
#include <numcosmo/math/ncm_prior.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_multifit.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_QRBF             (nc_hicosmo_qrbf_get_type ())
#define NC_HICOSMO_QRBF(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_QRBF, NcHICosmoQRBF))
#define NC_HICOSMO_QRBF_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_QRBF, NcHICosmoQRBFClass))
#define NC_IS_HICOSMO_QRBF(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_QRBF))
#define NC_IS_HICOSMO_QRBF_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_QRBF))
#define NC_HICOSMO_QRBF_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_QRBF, NcHICosmoQRBFClass))

typedef struct _NcHICosmoQRBFClass NcHICosmoQRBFClass;
typedef struct _NcHICosmoQRBF NcHICosmoQRBF;
typedef struct _NcHICosmoQRBFPrivate NcHICosmoQRBFPrivate;

/**
 * NcHICosmoQRBFSParams:
 * @NC_HICOSMO_QRBF_H0: Hubble constant
 * @NC_HICOSMO_QRBF_OMEGA_T: Total energy density of the universe
 * @NC_HICOSMO_QRBF_AS_DRAG: FIXME
 * @NC_HICOSMO_QRBF_RBF_H: the rbf wave-length 
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_QRBF_SPARAMS >*/
{
  NC_HICOSMO_QRBF_H0 = 0,
  NC_HICOSMO_QRBF_OMEGA_T,    
  NC_HICOSMO_QRBF_AS_DRAG,
  NC_HICOSMO_QRBF_RBF_H,
  /* < private > */
  NC_HICOSMO_QRBF_SPARAM_LEN, /*< skip >*/
} NcHICosmoQRBFSParams;

/**
 * NcHICosmoQRBFVParams:
 * @NC_HICOSMO_QRBF_RBF_CENTERS: RBF centers
 * @NC_HICOSMO_QRBF_RBF_COEFFS: RBF coefficients
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_QRBF_VPARAMS >*/
{
  NC_HICOSMO_QRBF_RBF_CENTERS,
  NC_HICOSMO_QRBF_RBF_COEFFS,
  /* < private > */
  NC_HICOSMO_QRBF_VPARAM_LEN, /*< skip >*/
} NcHICosmoQRBFVParams;

#define NC_HICOSMO_QRBF_DEFAULT_H0              (73.0)
#define NC_HICOSMO_QRBF_DEFAULT_OMEGA_T         ( 1.0)
#define NC_HICOSMO_QRBF_DEFAULT_AS_DRAG         ( 0.035)
#define NC_HICOSMO_QRBF_DEFAULT_RBF_H           ( 0.5)
#define NC_HICOSMO_QRBF_DEFAULT_RBF_CENTERS     ( 1.0)
#define NC_HICOSMO_QRBF_DEFAULT_RBF_CENTERS_LEN ( 3  )
#define NC_HICOSMO_QRBF_DEFAULT_RBF_COEFFS      (+0.5)
#define NC_HICOSMO_QRBF_DEFAULT_RBF_COEFFS_LEN  ( 3  )

struct _NcHICosmoQRBFClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

struct _NcHICosmoQRBF
{
  /*< private >*/
  NcHICosmo parent_instance;
  NcHICosmoQRBFPrivate *priv;
};

GType nc_hicosmo_qrbf_get_type (void) G_GNUC_CONST;

NcHICosmoQRBF *nc_hicosmo_qrbf_new (gsize np, gdouble z_f);

void nc_hicosmo_qrbf_set_z_f (NcHICosmoQRBF *qrbf, const gdouble z_f);
gdouble nc_hicosmo_qrbf_q_roughness (NcHICosmoQRBF *qrbf);

/*----------------------------------------------------------------------------*/

#define NC_TYPE_HICOSMO_QRBF_RPRIOR             (nc_hicosmo_qrbf_rprior_get_type ())
#define NC_HICOSMO_QRBF_RPRIOR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_QRBF_RPRIOR, NcHICosmoQRBFRprior))
#define NC_HICOSMO_QRBF_RPRIOR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_QRBF_RPRIOR, NcHICosmoQRBFRpriorClass))
#define NC_IS_HICOSMO_QRBF_RPRIOR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_QRBF_RPRIOR))
#define NC_IS_HICOSMO_QRBF_RPRIOR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_QRBF_RPRIOR))
#define NC_HICOSMO_QRBF_RPRIOR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_QRBF_RPRIOR, NcHICosmoQRBFRpriorClass))

typedef struct _NcHICosmoQRBFRpriorClass NcHICosmoQRBFRpriorClass;
typedef struct _NcHICosmoQRBFRprior NcHICosmoQRBFRprior;
typedef struct _NcHICosmoQRBFRpriorPrivate NcHICosmoQRBFRpriorPrivate;

struct _NcHICosmoQRBFRpriorClass
{
  /*< private >*/
  NcmPriorClass parent_class;
};

struct _NcHICosmoQRBFRprior
{
  /*< private >*/
  NcmPrior parent_instance;
  NcHICosmoQRBFRpriorPrivate *priv;
};

GType nc_hicosmo_qrbf_rprior_get_type (void) G_GNUC_CONST;

NcHICosmoQRBFRprior *nc_hicosmo_qrbf_rprior_new (gdouble lambda);

G_END_DECLS

#endif /* _NC_HICOSMO_QRBF_H_ */
