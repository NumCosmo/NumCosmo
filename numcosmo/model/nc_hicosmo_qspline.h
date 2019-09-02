/***************************************************************************
 *            nc_hicosmo_qspline.h
 *
 *  Wed February 15 11:31:28 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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

#ifndef _NC_HICOSMO_QSPLINE_H_
#define _NC_HICOSMO_QSPLINE_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/nc_hicosmo.h>
#include <numcosmo/math/ncm_fit.h>
#include <numcosmo/math/ncm_ode_spline.h>

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_multifit.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

#define NC_TYPE_HICOSMO_QSPLINE             (nc_hicosmo_qspline_get_type ())
#define NC_HICOSMO_QSPLINE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_QSPLINE, NcHICosmoQSpline))
#define NC_HICOSMO_QSPLINE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_QSPLINE, NcHICosmoQSplineClass))
#define NC_IS_HICOSMO_QSPLINE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_QSPLINE))
#define NC_IS_HICOSMO_QSPLINE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_QSPLINE))
#define NC_HICOSMO_QSPLINE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_QSPLINE, NcHICosmoQSplineClass))

typedef struct _NcHICosmoQSplineClass NcHICosmoQSplineClass;
typedef struct _NcHICosmoQSpline NcHICosmoQSpline;
typedef struct _NcHICosmoQSplineContPriorClass NcHICosmoQSplineContPriorClass;
typedef struct _NcHICosmoQSplineContPrior NcHICosmoQSplineContPrior;

/**
 * NcHICosmoQSplineSParams:
 * @NC_HICOSMO_QSPLINE_H0: Hubble constant
 * @NC_HICOSMO_QSPLINE_OMEGA_T: Total energy density of the universe
 * @NC_HICOSMO_QSPLINE_AS_DRAG: FIXME
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_QSPLINE_SPARAMS >*/
{
  NC_HICOSMO_QSPLINE_H0 = 0,
  NC_HICOSMO_QSPLINE_OMEGA_T,    
  NC_HICOSMO_QSPLINE_AS_DRAG,    
  /* < private > */
  NC_HICOSMO_QSPLINE_SPARAM_LEN, /*< skip >*/
} NcHICosmoQSplineSParams;

/**
 * NcHICosmoQSplineVParams:
 * @NC_HICOSMO_QSPLINE_Q: FIXME
 *
 */
typedef enum /*< enum,underscore_name=NC_HICOSMO_QSPLINE_VPARAMS >*/
{
  NC_HICOSMO_QSPLINE_Q,          
  /* < private > */
  NC_HICOSMO_QSPLINE_VPARAM_LEN, /*< skip >*/
} NcHICosmoQSplineVParams;

#define NC_HICOSMO_QSPLINE_DEFAULT_H0      ncm_c_hubble_cte_wmap ()
#define NC_HICOSMO_QSPLINE_DEFAULT_OMEGA_T ( 1.0)
#define NC_HICOSMO_QSPLINE_DEFAULT_AS_DRAG ( 0.035)
#define NC_HICOSMO_QSPLINE_DEFAULT_Q       (-0.5)
#define NC_HICOSMO_QSPLINE_DEFAULT_Q_LEN   (3)

struct _NcHICosmoQSplineClass
{
  /*< private >*/
  NcHICosmoClass parent_class;
};

struct _NcHICosmoQSpline
{
  /*< private >*/
  NcHICosmo parent_instance;
  guint nknots;
  guint size;
  gdouble z_f;
  NcmSpline *q_z;
  NcmOdeSpline *E2_z;
};

GType nc_hicosmo_qspline_get_type (void) G_GNUC_CONST;

NcHICosmoQSpline *nc_hicosmo_qspline_new (NcmSpline *s, gsize np, gdouble z_f);

NcHICosmoQSplineContPrior *nc_hicosmo_qspline_add_continuity_priors (NcHICosmoQSpline *qspline, NcmLikelihood *lh, gdouble sigma, gdouble abstol);

/****************************** Continuity Prior Model ******************************/

#define NC_TYPE_HICOSMO_QSPLINE_CONT_PRIOR             (nc_hicosmo_qspline_cont_prior_get_type ())
#define NC_HICOSMO_QSPLINE_CONT_PRIOR(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HICOSMO_QSPLINE_CONT_PRIOR, NcHICosmoQSplineContPrior))
#define NC_HICOSMO_QSPLINE_CONT_PRIOR_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HICOSMO_QSPLINE_CONT_PRIOR, NcHICosmoQSplineContPriorClass))
#define NC_IS_HICOSMO_QSPLINE_CONT_PRIOR(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HICOSMO_QSPLINE_CONT_PRIOR))
#define NC_IS_HICOSMO_QSPLINE_CONT_PRIOR_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HICOSMO_QSPLINE_CONT_PRIOR))
#define NC_HICOSMO_QSPLINE_CONT_PRIOR_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HICOSMO_QSPLINE_CONT_PRIOR, NcHICosmoQSplineContPriorClass))

struct _NcHICosmoQSplineContPrior
{
  /*< private >*/
  NcmModel parent_instance;
};

struct _NcHICosmoQSplineContPriorClass
{
  /*< private >*/
  NcmModelClass parent_class;
};

GType nc_hicosmo_qspline_cont_prior_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_hicosmo_qspline_cont_prior);

#define NC_HICOSMO_QSPLINE_CONT_PRIOR_ABSTOL (0)
#define NC_HICOSMO_QSPLINE_CONT_PRIOR_LNSIGMA (0)

NcHICosmoQSplineContPrior *nc_hicosmo_qspline_cont_prior_new (guint npriors);
NcHICosmoQSplineContPrior *nc_hicosmo_qspline_cont_prior_ref (NcHICosmoQSplineContPrior *qspline_cp);
void nc_hicosmo_qspline_cont_prior_free (NcHICosmoQSplineContPrior *qspline_cp);
void nc_hicosmo_qspline_cont_prior_set_lnsigma (NcHICosmoQSplineContPrior *qspline_cp, guint i, gdouble ln_sigma);
void nc_hicosmo_qspline_cont_prior_set_all_lnsigma (NcHICosmoQSplineContPrior *qspline_cp, gdouble ln_sigma);
gdouble nc_hicosmo_qspline_cont_prior_get_lnsigma (NcHICosmoQSplineContPrior *qspline_cp, guint i);

void nc_hicosmo_qspline_cont_prior_set_abstol (NcHICosmoQSplineContPrior *qspline_cp, gdouble abstol);
gdouble nc_hicosmo_qspline_cont_prior_get_abstol (NcHICosmoQSplineContPrior *qspline_cp);

/****************************** Continuity Prior ******************************/

#define NC_TYPE_PRIOR_QSPLINE_CONT             (nc_prior_qspline_cont_get_type ())
#define NC_PRIOR_QSPLINE_CONT(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_PRIOR_QSPLINE_CONT, NcPriorQSplineCont))
#define NC_PRIOR_QSPLINE_CONT_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_PRIOR_QSPLINE_CONT, NcPriorQSplineContClass))
#define NC_IS_PRIOR_QSPLINE_CONT(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_PRIOR_QSPLINE_CONT))
#define NC_IS_PRIOR_QSPLINE_CONT_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_PRIOR_QSPLINE_CONT))
#define NC_PRIOR_QSPLINE_CONT_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_PRIOR_QSPLINE_CONT, NcPriorQSplineContClass))

typedef struct _NcPriorQSplineContClass NcPriorQSplineContClass;
typedef struct _NcPriorQSplineCont NcPriorQSplineCont;

struct _NcPriorQSplineContClass
{
  /*< private >*/
  NcmPriorClass parent_class;
};

struct _NcPriorQSplineCont
{
  /*< private >*/
  NcmPrior parent_instance;
};

GType nc_prior_qspline_cont_get_type (void) G_GNUC_CONST;

NcPriorQSplineCont *nc_prior_qspline_cont_new (void);

G_END_DECLS

#endif /* _NC_HICOSMO_QSPLINE_H_ */
