/***************************************************************************
 *            nc_powspec_mnl_halofit.c
 *
 *  Thu March 17 14:57:40 2016
 *  Copyright  2016  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * nc_powspec_mnl_halofit.c
 * Copyright (C) 2016 Cyrille Doux <cdoux@apc.in2p3.fr>
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

/**
 * SECTION:nc_powspec_mnl_halofit
 * @title: NcPowspecMNLHaloFit
 * @short_description: nonlinear matter power spectrum from Halofit model.
 *
 * Provides the nonlinear matter power spectrum using Halofit model [Smith et al (2003)][XSmith2003]
 * and [Takahashi et al. (2012)][XTakahashi2012] FIXME.
 * 
 * For PKEqual see [Casarini et al. (2009)][XCasarini2009] and [Casarini et al. (2016)][XCasarini2016].
 * 
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hicosmo.h"
#include "nc_powspec_mnl_halofit.h"

#include "model/nc_hicosmo_de.h"
#include "model/nc_hicosmo_de_cpl.h"
#include "nc_distance.h"

#include "math/integral.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_spline_cubic_notaknot.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_exp.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcPowspecMNLHaloFitPrivate
{
  NcPowspecML *psml;
  gdouble zmaxnl;
  gdouble znl;
  gdouble reltol;
  NcmSpline *Rsigma;
  NcmSpline *neff;
  NcmSpline *Cur;
  NcmPowspecFilter *psml_gauss;
	gdouble z;
	gdouble ksigma;
	gdouble an;
	gdouble bn;
	gdouble cn;
	gdouble gamman;
	gdouble alphan;
	gdouble betan;
	gdouble nun;
	gdouble f1;
	gdouble f2;
	gdouble f3;
	gdouble mnu_corr_halo;
	gdouble fnu;
	gsl_root_fdfsolver* linear_scale_solver;
	gsl_root_fsolver* znl_solver;
  gboolean pkequal;
  NcHICosmo *cpl;
  NcDistance *cpl_dist;
};

enum
{
	PROP_0,
	PROP_PSML,
	PROP_ZMAXNL,
	PROP_RELTOL,
  PROP_PKEQUAL,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcPowspecMNLHaloFit, nc_powspec_mnl_halofit, NC_TYPE_POWSPEC_MNL);

static void
nc_powspec_mnl_halofit_init (NcPowspecMNLHaloFit *pshf)
{
	NcPowspecMNLHaloFitPrivate * const self = pshf->priv = nc_powspec_mnl_halofit_get_instance_private (pshf);

  self->psml       = NULL;
	self->zmaxnl     = 0.0;
	self->znl        = 0.0;
	self->reltol     = 0.0;
	self->Rsigma     = NULL;
	self->neff       = NULL;
	self->Cur        = NULL;
	self->psml_gauss = NULL;

	self->linear_scale_solver = gsl_root_fdfsolver_alloc (gsl_root_fdfsolver_steffenson);
	self->znl_solver          = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);

	self->z        = HUGE_VAL;
  self->pkequal  = FALSE;
  self->cpl      = NULL;
  self->cpl_dist = NULL;
}

static void
_nc_powspec_mnl_halofit_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
	NcPowspecMNLHaloFit *pshf = NC_POWSPEC_MNL_HALOFIT (object);
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;
  g_return_if_fail (NC_IS_POWSPEC_MNL_HALOFIT (object));

  switch (prop_id)
  {
    case PROP_PSML:
      self->psml = g_value_dup_object (value);
      self->psml_gauss = ncm_powspec_filter_new (NCM_POWSPEC (self->psml), NCM_POWSPEC_FILTER_TYPE_GAUSS);
      break;
    case PROP_ZMAXNL:
      self->zmaxnl = g_value_get_double (value);
      break;
    case PROP_RELTOL:
      self->reltol = g_value_get_double (value);
      break;
    case PROP_PKEQUAL:
      nc_powspec_mnl_halofit_pkequal (pshf, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_powspec_mnl_halofit_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
	NcPowspecMNLHaloFit *pshf = NC_POWSPEC_MNL_HALOFIT (object);
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;
	g_return_if_fail (NC_IS_POWSPEC_MNL_HALOFIT (object));

	switch (prop_id)
	{
    case PROP_PSML:
      g_value_set_object (value, self->psml);
      break;
    case PROP_ZMAXNL:
      g_value_set_double (value, self->zmaxnl);
      break;
    case PROP_RELTOL:
      g_value_set_double (value, self->reltol);
      break;
    case PROP_PKEQUAL:
      g_value_set_boolean (value, self->pkequal);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_powspec_mnl_halofit_constructed (GObject* object)
{
	/* Chain up : start */
	G_OBJECT_CLASS (nc_powspec_mnl_halofit_parent_class)->constructed (object);
	{
		NcPowspecMNLHaloFit *pshf = NC_POWSPEC_MNL_HALOFIT (object);
    NcPowspecMNLHaloFitPrivate * const self = pshf->priv;

		ncm_powspec_require_zf (NCM_POWSPEC (pshf), self->zmaxnl);

		ncm_powspec_require_zi (NCM_POWSPEC (self->psml), ncm_powspec_get_zi (NCM_POWSPEC (pshf)));
		ncm_powspec_require_zf (NCM_POWSPEC (self->psml), self->zmaxnl);

		self->Rsigma = ncm_spline_cubic_notaknot_new ();
		self->neff   = ncm_spline_cubic_notaknot_new ();
		self->Cur    = ncm_spline_cubic_notaknot_new ();
	}
}

static void
_nc_powspec_mnl_halofit_dispose (GObject* object)
{
	NcPowspecMNLHaloFit *pshf = NC_POWSPEC_MNL_HALOFIT (object);
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;

	nc_powspec_ml_clear (&self->psml);
	ncm_spline_clear (&self->Rsigma);
	ncm_spline_clear (&self->neff);
	ncm_spline_clear (&self->Cur);

	ncm_powspec_filter_clear (&self->psml_gauss);

  nc_hicosmo_clear (&self->cpl);
  nc_distance_clear (&self->cpl_dist);

	/* Chain up : end */
	G_OBJECT_CLASS (nc_powspec_mnl_halofit_parent_class)->dispose (object);
}

static void
_nc_powspec_mnl_halofit_finalize (GObject* object)
{
	NcPowspecMNLHaloFit *pshf = NC_POWSPEC_MNL_HALOFIT (object);
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;

	gsl_root_fdfsolver_free (self->linear_scale_solver);
	gsl_root_fsolver_free (self->znl_solver);

	/* Chain up : end */
	G_OBJECT_CLASS (nc_powspec_mnl_halofit_parent_class)->finalize (object);
}

static void _nc_powspec_mnl_halofit_prepare (NcmPowspec* powspec, NcmModel *model);
static gdouble _nc_powspec_mnl_halofit_eval (NcmPowspec* powspec, NcmModel *model, const gdouble z, const gdouble k);
static void _nc_powspec_mnl_halofit_eval_vec (NcmPowspec* powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk);
static void _nc_powspec_mnl_halofit_get_nknots (NcmPowspec* powspec, guint *Nz, guint *Nk);

static void
nc_powspec_mnl_halofit_class_init (NcPowspecMNLHaloFitClass* klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcmPowspecClass* powspec_class = NCM_POWSPEC_CLASS (klass);

	object_class->set_property = &_nc_powspec_mnl_halofit_set_property;
	object_class->get_property = &_nc_powspec_mnl_halofit_get_property;

	object_class->constructed = &_nc_powspec_mnl_halofit_constructed;
	object_class->dispose = &_nc_powspec_mnl_halofit_dispose;
	object_class->finalize = &_nc_powspec_mnl_halofit_finalize;

	g_object_class_install_property (object_class,
	                                 PROP_PSML,
	                                 g_param_spec_object ("power-spec",
	                                                      NULL,
	                                                      "Linear power spectrum.",
	                                                      NC_TYPE_POWSPEC_ML,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_ZMAXNL,
	                                 g_param_spec_double ("zmaxnl",
	                                                      NULL,
	                                                      "Max redshift for halofit correction",
	                                                      0.0, 10000.0, 10.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_RELTOL,
	                                 g_param_spec_double ("reltol",
	                                                      NULL,
	                                                      "Relative tolerance (precision) for halofit computations",
	                                                      GSL_DBL_EPSILON, 1.0, 1e-3,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_PKEQUAL,
                                   g_param_spec_boolean ("use-pkequal",
                                                         NULL,
                                                         "Whether to use PKEqual",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
	powspec_class->prepare    = &_nc_powspec_mnl_halofit_prepare;
	powspec_class->eval       = &_nc_powspec_mnl_halofit_eval;
	powspec_class->eval_vec   = &_nc_powspec_mnl_halofit_eval_vec;
	powspec_class->get_nknots = &_nc_powspec_mnl_halofit_get_nknots;
}

typedef struct _int_var_moment_params
{
	const gint n;
	const gdouble R;
	const gdouble z;
	NcPowspecML *ps;
	NcHICosmo *cosmo;
} int_var_moment_params;

///////////////////// INTEGRATION OVER K*R ////////////////////////////
static gdouble
_nc_powspec_mnl_halofit_var_moment_integrand (gdouble kR, gpointer params)
{
	int_var_moment_params* ts = (int_var_moment_params *)params;
	const gdouble k        = kR / ts->R;
	const gdouble matter_P = ncm_powspec_eval (NCM_POWSPEC (ts->ps), NCM_MODEL (ts->cosmo), ts->z, k);
	const gdouble kR2      = kR * kR;
	const gdouble W2       = exp (-kR2);
/*	printf ("NC: % 22.15g % 22.15g % 22.15g % 22.15g % 22.15g %i % 22.15g % 22.15g\n", k, ts->R, ts->z, matter_P, kR2, ts->n, W2, 
          matter_P * gsl_pow_int (kR2, ts->n + 1) * W2 / (gsl_pow_3 (ts->R) * (gdouble)ncm_c_2_pi_2 ()));*/
	return matter_P * gsl_pow_int (kR2, ts->n + 1) * W2;
}

static gdouble
_nc_powspec_mnl_halofit_var_moment (NcPowspecML *ps, NcHICosmo *cosmo, const gdouble R, const gdouble z, const gint n)
{
	gdouble result, error;
	gsl_function F;
	int_var_moment_params ivmps = { n, R, z, ps, cosmo };

	gsl_integration_workspace *w = gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);

	F.function = &_nc_powspec_mnl_halofit_var_moment_integrand;
	F.params = &ivmps;

	//gsl_integration_qagiu (&F, NCM_DEFAULT_PRECISION, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, &result, &error);

  //gsl_integration_qagiu (&F, NCM_DEFAULT_PRECISION, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, &result, &error);

  gsl_integration_qag (&F, 7.1972528136862e-06 * R, 57.2224908217619 * R, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, w, &result, &error);

 
  
	gsl_integration_workspace_free (w);

	return result / (gsl_pow_3 (R) * ncm_c_2_pi_2 ());
}

typedef struct _var_params
{
	NcPowspecMNLHaloFit *pshf;
	NcHICosmo *cosmo;
	const gdouble z;
	const gdouble R_min;
} var_params;

static gdouble
_nc_powspec_mnl_halofit_varm1 (gdouble lnR, gpointer params)
{
	var_params *vps = (var_params *) params;
  NcPowspecMNLHaloFitPrivate * const self = vps->pshf->priv;
	/*
  const gdouble R = gsl_sf_exp (lnR);
  const gdouble sigma2_0 = _nc_powspec_mnl_halofit_var_moment (vps->self->psml, vps->cosmo, R, vps->z, 0);
  return sigma2_0 - 1.0;
*/
	return ncm_powspec_filter_eval_lnvar_lnr (self->psml_gauss, vps->z, lnR);
}

static gdouble
_nc_powspec_mnl_halofit_varm1_deriv (gdouble lnR, gpointer params)
{
	var_params *vps = (var_params *) params;
  NcPowspecMNLHaloFitPrivate * const self = vps->pshf->priv;
	/*
  const gdouble R = gsl_sf_exp (lnR);
  const gdouble sigma2_1 = _nc_powspec_mnl_halofit_var_moment (vps->pshf->psml, vps->cosmo, R, vps->z, 1);
  return -2.0 * sigma2_1;
*/

	return ncm_powspec_filter_eval_dlnvar_dlnr (self->psml_gauss, vps->z, lnR);
}

static void
_nc_powspec_mnl_halofit_varm1_fdf (gdouble lnR, gpointer params, gdouble* varm1, gdouble* dvarm1)
{
	var_params *vps = (var_params *) params;
  NcPowspecMNLHaloFitPrivate * const self = vps->pshf->priv;

	/*
  const gdouble R = gsl_sf_exp (lnR);
  const gdouble sigma2_0 = _nc_powspec_mnl_halofit_var_moment (vps->pshf->psml, vps->cosmo, R, vps->z, 0);
  const gdouble sigma2_1 = _nc_powspec_mnl_halofit_var_moment (vps->pshf->psml, vps->cosmo, R, vps->z, 1);
*/
	*varm1 = ncm_powspec_filter_eval_lnvar_lnr (self->psml_gauss, vps->z, lnR);
	*dvarm1 = ncm_powspec_filter_eval_dlnvar_dlnr (self->psml_gauss, vps->z, lnR);
}

typedef struct _root_params
{
	NcPowspecML *ps;
	NcHICosmo *cosmo;
} root_params;

static gdouble
_nc_powspec_mnl_halofit_linear_scale (NcPowspecMNLHaloFit *pshf, NcHICosmo *cosmo, const gdouble z)
{
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;
  
	gdouble lnR0         = 0.0;
	gdouble lnR          = (-z / 2.0 < NC_POWSPEC_MNL_HALOFIT_LOGRMIN) ? NC_POWSPEC_MNL_HALOFIT_LOGRMIN : -z / 2.0 + NCM_DEFAULT_PRECISION;
	const gdouble reltol = self->reltol / 10.0;
	gdouble res          = 0.0;
  gint max_iter        = 20000;
	gint iter            = 0;
	gint status;

	gsl_function_fdf FDF;

	var_params vps = { pshf, cosmo, z, 0.0 };

	FDF.f      = &_nc_powspec_mnl_halofit_varm1;
	FDF.df     = &_nc_powspec_mnl_halofit_varm1_deriv;
	FDF.fdf    = &_nc_powspec_mnl_halofit_varm1_fdf;
	FDF.params = &vps;

	gsl_root_fdfsolver_set (self->linear_scale_solver, &FDF, lnR);

	do
	{
		iter++;
		status = gsl_root_fdfsolver_iterate (self->linear_scale_solver);

		lnR0 = lnR;
		lnR = gsl_root_fdfsolver_root (self->linear_scale_solver);

		res = gsl_expm1 (lnR0 - lnR);

		status = gsl_root_test_residual (res, reltol); //Compares R vs R0 !

	} while (status == GSL_CONTINUE && iter < max_iter);

	res = exp (lnR); // Now res is the result

	if (iter >= max_iter)
		g_warning ("_nc_powspec_mnl_halofit_linear_scale: maximum number of iteration reached (%u), non-linear scale found R(z=%.3f).", max_iter, z);

	if (!(gsl_finite (res)))
		g_warning ("_nc_powspec_mnl_halofit_linear_scale: non-linear scale found R(z=%.3f) not finite.", z);

	if ((iter >= max_iter) || !(gsl_finite (res)))
	{
		g_message ("_nc_powspec_mnl_halofit_linear_scale: running the bracketing solver...");

		iter = 0;

		gsl_root_fsolver* s = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);

		gsl_function F;
		F.function = &_nc_powspec_mnl_halofit_varm1;
		F.params = &vps;

		gdouble lnRlo = -z / 2.0 - 10.0;
		gdouble lnRup = -z / 2.0 + 10.0;

		gsl_root_fsolver_set (s, &F, lnRlo, lnRup);

		do
		{
			iter++;
			status = gsl_root_fsolver_iterate (s);

			lnR = gsl_root_fsolver_root (s);
			lnRlo = gsl_root_fsolver_x_lower (s);
			lnRup = gsl_root_fsolver_x_upper (s);

			status = gsl_root_test_interval (lnRlo, lnRup, reltol, 0.0); // dlnR = dR/R, so test with abstol

		} while (status == GSL_CONTINUE && iter < max_iter);

		gsl_root_fsolver_free (s);

		res = exp (lnR);
	}

	return res;
}

static gdouble
_nc_powspec_mnl_halofit_linear_scale_z (gdouble z, gpointer params)
{
	var_params *vps = (var_params*)params;

	return _nc_powspec_mnl_halofit_linear_scale (vps->pshf, vps->cosmo, z);
}

static gdouble
_nc_powspec_mnl_halofit_linear_scale_z_znl (gdouble z, gpointer params)
{
	var_params *vps = (var_params*)params;

	return _nc_powspec_mnl_halofit_linear_scale (vps->pshf, vps->cosmo, z) - vps->R_min;
}

static void
_nc_powspec_mnl_halofit_prepare_nl (NcPowspecMNLHaloFit *pshf, NcmModel *model)
{
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;
	NcHICosmo *cosmo = NC_HICOSMO (model);
	guint i;

	self->z = HUGE_VAL;
  
	ncm_powspec_filter_set_zi (self->psml_gauss, ncm_powspec_get_zi (NCM_POWSPEC (pshf)));
	ncm_powspec_filter_set_zf (self->psml_gauss, self->zmaxnl);

	ncm_powspec_filter_set_best_lnr0 (self->psml_gauss);
	ncm_powspec_filter_prepare_if_needed (self->psml_gauss, model);

	{
		const gdouble R_min = ncm_powspec_filter_get_r_min (self->psml_gauss);
		var_params vps      = {pshf, cosmo, 0.0, R_min};
		gsl_function Fznl;

    Fznl.function = &_nc_powspec_mnl_halofit_linear_scale_z_znl;
		Fznl.params   = &vps;

		if (_nc_powspec_mnl_halofit_linear_scale_z (self->zmaxnl, Fznl.params) > R_min)
			self->znl = self->zmaxnl;
		else
		{
			const gdouble step = 1.0;
			gdouble z0         = 0.0;
			gdouble z1         = z0 + step;
			gdouble R0         = _nc_powspec_mnl_halofit_linear_scale_z (z0, Fznl.params);
			gdouble R1         = _nc_powspec_mnl_halofit_linear_scale_z (z1, Fznl.params);
			gboolean found     = FALSE;

			if (R0 < R_min)
				g_error ("_nc_powspec_mnl_halofit_prepare_nl: linear universe or too large R_min, in the latter case increase k_max (R0 == % 21.15g, R_min == % 21.15g).", R0, R_min);

			while (R1 > R_min)
			{
				z0 = z1;
				z1 += step;
				R1 = _nc_powspec_mnl_halofit_linear_scale_z (z1, Fznl.params);
				if (z1 > self->zmaxnl)
				{
					self->znl = self->zmaxnl;
					found = TRUE;
					break;
				}
			}

			if (!found)
			{
				gint status;
				gint iter = 0, max_iter = 20000;

				gsl_root_fsolver_set (self->znl_solver, &Fznl, z0, z1);
				do
				{
					iter++;
					status = gsl_root_fsolver_iterate (self->znl_solver);

					self->znl = gsl_root_fsolver_root (self->znl_solver);
					z0 = gsl_root_fsolver_x_lower (self->znl_solver);
					z1 = gsl_root_fsolver_x_upper (self->znl_solver);
					status = gsl_root_test_interval (z0, z1, 0.0, 1.0e-3);
				} while (status == GSL_CONTINUE && iter < max_iter);

				if (iter >= max_iter)
					g_warning ("_nc_powspec_mnl_halofit_prepare_nl: maximum number of iteration reached (%u), giving up.", max_iter);
			}
		}

		{
			gsl_function F;
			F.function = &_nc_powspec_mnl_halofit_linear_scale_z;
			F.params = &vps;

			ncm_spline_set_func (self->Rsigma, NCM_SPLINE_FUNCTION_SPLINE, &F, 0.0, self->znl, 0, self->reltol);
		}
	}

	{
		const guint len  = ncm_vector_len (self->Rsigma->xv);
		NcmVector* neffv = ncm_vector_new (len);
		NcmVector* Curv  = ncm_vector_new (len);

		for (i = 0; i < len; i++)
		{
			if (FALSE)
			{
				const gdouble z = ncm_vector_get (self->Rsigma->xv, i);
				const gdouble R = ncm_vector_get (self->Rsigma->yv, i);
				const gdouble sigma2_0 = _nc_powspec_mnl_halofit_var_moment (self->psml, cosmo, R, z, 0);
				const gdouble sigma2_1 = _nc_powspec_mnl_halofit_var_moment (self->psml, cosmo, R, z, 1);
				const gdouble sigma2_2 = _nc_powspec_mnl_halofit_var_moment (self->psml, cosmo, R, z, 2);
				const gdouble d1 = -2.0 * sigma2_1 / sigma2_0;
				const gdouble d2 = -d1 * d1 + 4.0 * (-sigma2_1 + sigma2_2) / sigma2_0;

				ncm_vector_set (neffv, i, -3.0 - d1);
				ncm_vector_set (Curv, i, -d2);
/*
				printf ("# z = % 20.15g, R = % 20.15g | % 20.15g % 20.15g % 20.15g | %e %e %e\n",
				        z, R,
				        sigma2_0, d1, d2,
				        fabs ((ncm_powspec_filter_eval_var (self->psml_gauss, z, R) - sigma2_0) / sigma2_0),
				        fabs ((ncm_powspec_filter_eval_dlnvar_dlnr (self->psml_gauss, z, log (R)) - d1) / d1),
				        fabs ((ncm_powspec_filter_eval_dnlnvar_dlnrn (self->psml_gauss, z, log (R), 2) - d2) / d2));
*/
      }
			else
			{
				const gdouble z = ncm_vector_get (self->Rsigma->xv, i);
				const gdouble R = ncm_vector_get (self->Rsigma->yv, i);
				const gdouble lnR = log (R);
				const gdouble d1 = ncm_powspec_filter_eval_dlnvar_dlnr (self->psml_gauss, z, lnR);
				const gdouble d2 = ncm_powspec_filter_eval_dnlnvar_dlnrn (self->psml_gauss, z, lnR, 2);

				ncm_vector_set (neffv, i, -3.0 - d1);
				ncm_vector_set (Curv, i, -d2);
			}
		}

		ncm_spline_set (self->neff, self->Rsigma->xv, neffv, TRUE);
		ncm_spline_set (self->Cur, self->Rsigma->xv, Curv, TRUE);

		ncm_vector_free (neffv);
		ncm_vector_free (Curv);
	}
}

static void
_nc_powspec_mnl_halofit_prepare (NcmPowspec* powspec, NcmModel *model)
{
	NcPowspecMNLHaloFit *pshf = NC_POWSPEC_MNL_HALOFIT (powspec);
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;

	g_assert (NC_IS_HICOSMO (model));

	ncm_powspec_require_zi (NCM_POWSPEC (self->psml), ncm_powspec_get_zi (NCM_POWSPEC (pshf)));
	ncm_powspec_require_zf (NCM_POWSPEC (self->psml), ncm_powspec_get_zf (NCM_POWSPEC (pshf)));

	ncm_powspec_require_kmin (NCM_POWSPEC (self->psml), ncm_powspec_get_kmin (NCM_POWSPEC (pshf)));
	ncm_powspec_require_kmax (NCM_POWSPEC (self->psml), ncm_powspec_get_kmax (NCM_POWSPEC (pshf)));

	ncm_powspec_prepare (NCM_POWSPEC (self->psml), model);

	_nc_powspec_mnl_halofit_prepare_nl (pshf, model);
}

static gdouble
_nc_powspec_mnl_halofit_xi_cmp (gdouble w, gpointer params)
{
	var_params *vps = (var_params *)params;
  NcPowspecMNLHaloFitPrivate * const self = vps->pshf->priv;

  ncm_model_orig_param_set (NCM_MODEL (vps->cosmo), NC_HICOSMO_DE_CPL_W0, w);
  /*ncm_model_params_log_all (NCM_MODEL (vps->cosmo));*/

  {
    gdouble zdrag = nc_distance_drag_redshift (self->cpl_dist, vps->cosmo);
    gdouble xi    = nc_distance_comoving (self->cpl_dist, vps->cosmo, zdrag) - nc_distance_comoving (self->cpl_dist, vps->cosmo, vps->z);

    /*printf ("CMP % 22.15g % 22.15g | % 22.15g % 22.15g\n", vps->z, zdrag, xi, vps->R_min);*/

    return xi / vps->R_min - 1.0;
  }
}

static void
_nc_powspec_mnl_halofit_preeval (NcPowspecMNLHaloFit *pshf, NcHICosmo *cosmo, const gdouble z)
{
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;
  
	const gdouble E2       = nc_hicosmo_E2 (cosmo, z);
	const gdouble Rsigma   = ncm_spline_eval (self->Rsigma, z);

	const gdouble neff     = ncm_spline_eval (self->neff, z);
	const gdouble Cur      = ncm_spline_eval (self->Cur, z);

	const gdouble neff2    = neff * neff;
	const gdouble neff3    = neff2 * neff;
	const gdouble neff4    = neff2 * neff2;

	const gdouble Omega_m  = nc_hicosmo_E2Omega_m (cosmo, z) / E2;
	const gdouble fnu      = nc_hicosmo_E2Omega_mnu (cosmo, z) / nc_hicosmo_E2Omega_m (cosmo, z);
	const gdouble frac     = nc_hicosmo_de_E2Omega_de (NC_HICOSMO_DE (cosmo), z) / (E2 - nc_hicosmo_E2Omega_m (cosmo, z));

	gdouble Omega_de_onepw = 0.0;
/*
  printf ("#NC:  z % 22.15g E2 % 22.15g Rsigma % 22.15g neff % 22.15g Cur % 22.15g Omega_m % 22.15g fnu % 22.15g frac % 22.15g | % 22.15g\n",
          z, E2, Rsigma, neff, Cur, Omega_m, fnu, frac, 
          _nc_powspec_mnl_halofit_var_moment (self->psml, cosmo, Rsigma, z, 0)
          );
*/
  if (NC_IS_HICOSMO_DE (cosmo))
  {
    gdouble wa;
    if (self->pkequal && NC_IS_HICOSMO_DE_CPL (cosmo) && ((wa = ncm_model_orig_param_get (NCM_MODEL (cosmo), NC_HICOSMO_DE_CPL_W1)) != 0.0) )
    {
      /*printf ("Using pkequal!\n");*/
      if (self->cpl == NULL)
      {
        self->cpl      = NC_HICOSMO (nc_hicosmo_de_cpl_new ());
        self->cpl_dist = nc_distance_new (1.0);
      }

      ncm_vector_memcpy (ncm_model_orig_params_peek_vector (NCM_MODEL (self->cpl)), ncm_model_orig_params_peek_vector (NCM_MODEL (cosmo)));
      ncm_model_orig_param_set (NCM_MODEL (self->cpl), NC_HICOSMO_DE_CPL_W1, 0.0);
      
      {
        gdouble zdrag_cosmo, xi_cosmo; 
        nc_distance_prepare (self->cpl_dist, cosmo);

        zdrag_cosmo = nc_distance_drag_redshift (self->cpl_dist, cosmo);
        xi_cosmo    = nc_distance_comoving (self->cpl_dist, cosmo, zdrag_cosmo) - nc_distance_comoving (self->cpl_dist, cosmo, z);

        {
          var_params vps = { pshf, self->cpl, z, xi_cosmo };

          gsl_function Fxi;

          gint status;
          gint iter = 0, max_iter = 20000;
          gdouble w0 = -2.0, w1 = 0.0, wfinal = 0.0;

          Fxi.function = &_nc_powspec_mnl_halofit_xi_cmp;
          Fxi.params = &vps;
          
          gsl_root_fsolver_set (self->znl_solver, &Fxi, w0, w1);
          do
          {
            iter++;
            status = gsl_root_fsolver_iterate (self->znl_solver);

            wfinal = gsl_root_fsolver_root (self->znl_solver);
            w0     = gsl_root_fsolver_x_lower (self->znl_solver);
            w1     = gsl_root_fsolver_x_upper (self->znl_solver);
            status = gsl_root_test_interval (w0, w1, 0.0, 1.0e-6);
          } while (status == GSL_CONTINUE && iter < max_iter);

          ncm_model_orig_param_set (NCM_MODEL (self->cpl), NC_HICOSMO_DE_CPL_W0, wfinal);
          if (iter >= max_iter)
            g_warning ("_nc_powspec_mnl_halofit_preeval: maximum number of iteration reached (%u), giving up.", max_iter);
/*
          printf ("# for z = % 22.15g, using wequiv = % 22.15g vs (w0 = % 22.15g, w1 = % 22.15g)\n", 
                  z, wfinal, 
                  ncm_model_orig_param_get (NCM_MODEL (cosmo), NC_HICOSMO_DE_CPL_W0),
                  ncm_model_orig_param_get (NCM_MODEL (cosmo), NC_HICOSMO_DE_CPL_W1)
                  );
*/          
          Omega_de_onepw = nc_hicosmo_de_E2Omega_de_onepw (NC_HICOSMO_DE (self->cpl), z) / nc_hicosmo_E2 (self->cpl, z);
        }        
        
      }
    }
    else
    {
      /*printf ("Not using pkequal!\n");*/
      Omega_de_onepw = nc_hicosmo_de_E2Omega_de_onepw (NC_HICOSMO_DE (cosmo), z) / E2;
    }
  }
  
	self->z             = z;
	self->ksigma        = 1.0 / Rsigma;

	self->an            = ncm_util_exp10 (+1.5222 + 2.8553 * neff + 2.3706 * neff2 + 0.9903 * neff3 + 0.2250 * neff4 - 0.6038 * Cur + 0.1749 * Omega_de_onepw);
	self->bn            = ncm_util_exp10 (-0.5642 + 0.5864 * neff + 0.5716 * neff2 - 1.5474 * Cur + 0.2279 * Omega_de_onepw);
	self->cn            = ncm_util_exp10 (+0.3698 + 2.0404 * neff + 0.8161 * neff2 + 0.5869 * Cur);
	self->gamman        = 0.1971 - 0.0843 * neff + 0.8460 * Cur;
	self->alphan        = fabs (6.0835 + 1.3373 * neff - 0.1959 * neff2 - 5.5274 * Cur);
	self->betan         = 2.0379 - 0.7354 * neff + 0.3157 * neff2 + 1.2490 * neff3 + 0.3980 * neff4 - 0.1682 * Cur + fnu * (1.081 + 0.395 * neff2);
	self->nun           = ncm_util_exp10 (5.2105 + 3.6902 * neff);

	self->f1            = frac * pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F1bPOW) + (1.0 - frac) * pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F1aPOW);
	self->f2            = frac * pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F2bPOW) + (1.0 - frac) * pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F2aPOW);
	self->f3            = frac * pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F3bPOW) + (1.0 - frac) * pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F3aPOW);

	self->mnu_corr_halo = 1.0 + fnu * (0.977 - 18.015 * (nc_hicosmo_Omega_m0 (cosmo) - 0.3));
	self->fnu           = fnu;
}

static gdouble _nc_powspec_mnl_halofit_Pklin2Pknln (NcPowspecMNLHaloFit *pshf, NcHICosmo *cosmo, const gdouble k, const gdouble Pklin)
{
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;
  
	const gdouble kh2          = gsl_pow_2 (k / nc_hicosmo_h (cosmo));
	const gdouble k3           = gsl_pow_3 (k);
	const gdouble k3o2pi2      = k3 / ncm_c_2_pi_2 ();
	const gdouble Delta_lin    = k3o2pi2 * Pklin;

	const gdouble y            = k / self->ksigma;

	const gdouble Delta_lin_nu = Delta_lin * (1.0 + self->fnu * 47.48 * kh2 / (1.0 + 1.5 * kh2));
	const gdouble P_Q          = Pklin * (pow (1.0 + Delta_lin_nu, self->betan) / (1.0 + self->alphan * Delta_lin_nu)) * exp (-y / 4.0 - y * y / 8.0);

	const gdouble Delta_Hprime = self->an * pow (y, 3.0 * self->f1) / (1.0 + self->bn * pow (y, self->f2) + pow (self->cn * self->f3 * y, 3.0 - self->gamman));
	const gdouble Delta_H      = Delta_Hprime / (1.0 + self->nun / (y * y)) * self->mnu_corr_halo;

	const gdouble P_H          = Delta_H / k3o2pi2;

	return P_Q + P_H;
}

static gdouble
_nc_powspec_mnl_halofit_eval (NcmPowspec* powspec, NcmModel *model, const gdouble z, const gdouble k)
{
	NcHICosmo *cosmo           = NC_HICOSMO (model);
	NcPowspecMNLHaloFit *pshf  = NC_POWSPEC_MNL_HALOFIT (powspec);
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;
	const gboolean linscale    = (z > self->znl);
	const gboolean applysmooth = (z + 1.0 > self->znl);
	const gdouble zhf          = linscale ? self->znl : z;
	const gdouble Pklin        = ncm_powspec_eval (NCM_POWSPEC (self->psml), model, z, k);
	gdouble Pknln;

	if (zhf != self->z)
	{
		_nc_powspec_mnl_halofit_preeval (pshf, cosmo, zhf);
	}

	Pknln = _nc_powspec_mnl_halofit_Pklin2Pknln (pshf, cosmo, k, Pklin);

	if (applysmooth)
		Pknln = ncm_util_smooth_trans (Pknln, Pklin, self->znl, 1.0, z);

	return Pknln;
}

static void
_nc_powspec_mnl_halofit_eval_vec (NcmPowspec* powspec, NcmModel *model, const gdouble z, NcmVector *k, NcmVector *Pk)
{
	NcHICosmo *cosmo           = NC_HICOSMO (model);
	NcPowspecMNLHaloFit *pshf  = NC_POWSPEC_MNL_HALOFIT (powspec);
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;
	const gboolean linscale    = (z > self->znl);
	const gboolean applysmooth = (z + 1.0 > self->znl);
	const gdouble zhf          = linscale ? self->znl : z;
	gdouble theta0, theta1;

	ncm_powspec_eval_vec (NCM_POWSPEC (self->psml), model, z, k, Pk);

	if (applysmooth)
		ncm_util_smooth_trans_get_theta (self->znl, 1.0, z, &theta0, &theta1);

	if (zhf != self->z)
	{
		_nc_powspec_mnl_halofit_preeval (pshf, cosmo, zhf);
	}

	{
		const guint len = ncm_vector_len (k);
		guint i;

		for (i = 0; i < len; i++)
		{
			const gdouble ki = ncm_vector_get (k, i);
			const gdouble Pklin = ncm_vector_get (Pk, i);

			const gdouble Pknln = _nc_powspec_mnl_halofit_Pklin2Pknln (pshf, cosmo, ki, Pklin);

			if (applysmooth)
				ncm_vector_set (Pk, i, theta0 * Pknln + theta1 * Pklin);
			else
				ncm_vector_set (Pk, i, Pknln);
		}
	}
}

static void
_nc_powspec_mnl_halofit_get_nknots (NcmPowspec* powspec, guint *Nz, guint *Nk)
{
	NcPowspecMNLHaloFit *pshf = NC_POWSPEC_MNL_HALOFIT (powspec);
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;
  
	ncm_powspec_get_nknots (NCM_POWSPEC (self->psml), Nz, Nk);
	*Nz = ncm_vector_len (self->Rsigma->xv);
}

/**
 * nc_powspec_mnl_halofit_new:
 * @psml: a #NcPowspecML
 * @zmaxnl: a gdouble
 * @reltol: a gdouble
 *
 * Creates a new #NcPowspecMNLHaloFit from the transfer
 * function @tf.
 *
 * Returns: (transfer full): the newly created #NcPowspecMNLHaloFit.
 */
NcPowspecMNLHaloFit*
nc_powspec_mnl_halofit_new (NcPowspecML *psml, gdouble zmaxnl, gdouble reltol)
{
	NcPowspecMNLHaloFit *pshf = g_object_new (NC_TYPE_POWSPEC_MNL_HALOFIT,
	                                          "power-spec", psml,
	                                          "zmaxnl", zmaxnl,
	                                          "reltol", reltol,
	                                          NULL);

	return pshf;
}

/**
 * nc_powspec_mnl_halofit_set_kbounds_from_ml:
 * @pshf: a #NcPowspecMNLHaloFit
 *
 * Sets mode $k$ boundaries from the linear matter power spectrum.
 *
 */
void 
nc_powspec_mnl_halofit_set_kbounds_from_ml (NcPowspecMNLHaloFit *pshf)
{
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;
  
	ncm_powspec_set_kmin (NCM_POWSPEC (pshf), ncm_powspec_get_kmin (NCM_POWSPEC (self->psml)));
	ncm_powspec_set_kmax (NCM_POWSPEC (pshf), ncm_powspec_get_kmax (NCM_POWSPEC (self->psml)));
}

/**
 * nc_powspec_mnl_halofit_pkequal:
 * @pshf: a #NcPowspecMNLHaloFit
 * @on: a boolean
 *
 * Whether to use PKEqual to adjust the HaloFit formula when using a #NcHICosmoDECpl
 * model, see [Casarini et al. (2009)][XCasarini2009] and [Casarini et al. (2016)][XCasarini2016].
 *
 */
void 
nc_powspec_mnl_halofit_pkequal (NcPowspecMNLHaloFit *pshf, gboolean on)
{
  NcPowspecMNLHaloFitPrivate * const self = pshf->priv;

  if (on)
  {
    if (!self->pkequal)
    {
      self->pkequal = on;
      ncm_model_ctrl_force_update (NCM_POWSPEC (pshf)->ctrl);
    }
  }
  else
  {
    if (self->pkequal)
    {
      self->pkequal = on;
      ncm_model_ctrl_force_update (NCM_POWSPEC (pshf)->ctrl);
    }
  }
}
