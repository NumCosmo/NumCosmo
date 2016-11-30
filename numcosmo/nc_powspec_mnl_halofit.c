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
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_hicosmo.h"
#include "nc_powspec_mnl_halofit.h"

#include "model/nc_hicosmo_de.h"

#include "math/integral.h"
#include "math/memory_pool.h"
#include "math/ncm_spline_cubic_notaknot.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_exp.h>

struct _NcPowspecMNLHaloFitPrivate
{
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
	gsl_root_fdfsolver* linear_scale_solver;
	gsl_root_fsolver* znl_solver;
};

enum
{
	PROP_0,
	PROP_PSML,
	PROP_ZMAXNL,
	PROP_RELTOL,
};

G_DEFINE_TYPE (NcPowspecMNLHaloFit, nc_powspec_mnl_halofit, NC_TYPE_POWSPEC_MNL);

static void
nc_powspec_mnl_halofit_init (NcPowspecMNLHaloFit* pshf)
{
	pshf->psml = NULL;
	pshf->zmaxnl = 0.0;
	pshf->znl = 0.0;
	pshf->reltol = 0.0;
	pshf->Rsigma = NULL;
	pshf->neff = NULL;
	pshf->Cur = NULL;
	pshf->psml_gauss = NULL;
	pshf->priv = G_TYPE_INSTANCE_GET_PRIVATE (pshf, NC_TYPE_POWSPEC_MNL_HALOFIT, NcPowspecMNLHaloFitPrivate);

	pshf->priv->linear_scale_solver = gsl_root_fdfsolver_alloc (gsl_root_fdfsolver_steffenson);
	pshf->priv->znl_solver = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);

	pshf->priv->z = HUGE_VAL;
}

static void
_nc_powspec_mnl_halofit_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
	NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (object);
	g_return_if_fail (NC_IS_POWSPEC_MNL_HALOFIT (object));

	switch (prop_id)
	{
	case PROP_PSML:
		pshf->psml = g_value_dup_object (value);
		pshf->psml_gauss = ncm_powspec_filter_new (NCM_POWSPEC (pshf->psml), NCM_POWSPEC_FILTER_TYPE_GAUSS);
		break;
	case PROP_ZMAXNL:
		pshf->zmaxnl = g_value_get_double (value);
		break;
	case PROP_RELTOL:
		pshf->reltol = g_value_get_double (value);
		break;
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
_nc_powspec_mnl_halofit_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
	NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (object);
	g_return_if_fail (NC_IS_POWSPEC_MNL_HALOFIT (object));

	switch (prop_id)
	{
	case PROP_PSML:
		g_value_set_object (value, pshf->psml);
		break;
	case PROP_ZMAXNL:
		g_value_set_double (value, pshf->zmaxnl);
		break;
	case PROP_RELTOL:
		g_value_set_double (value, pshf->reltol);
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
	G_OBJECT_CLASS (nc_powspec_mnl_halofit_parent_class)
	->constructed (object);
	{
		NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (object);

		ncm_powspec_require_zf (NCM_POWSPEC (pshf), pshf->zmaxnl);

		ncm_powspec_require_zi (NCM_POWSPEC (pshf->psml), ncm_powspec_get_zi (NCM_POWSPEC (pshf)));
		ncm_powspec_require_zf (NCM_POWSPEC (pshf->psml), pshf->zmaxnl);

		pshf->Rsigma = ncm_spline_cubic_notaknot_new ();
		pshf->neff = ncm_spline_cubic_notaknot_new ();
		pshf->Cur = ncm_spline_cubic_notaknot_new ();
	}
}

static void
_nc_powspec_mnl_halofit_dispose (GObject* object)
{
	NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (object);

	nc_powspec_ml_clear (&pshf->psml);
	ncm_spline_clear (&pshf->Rsigma);
	ncm_spline_clear (&pshf->neff);
	ncm_spline_clear (&pshf->Cur);

	ncm_powspec_filter_clear (&pshf->psml_gauss);

	/* Chain up : end */
	G_OBJECT_CLASS (nc_powspec_mnl_halofit_parent_class)
	->dispose (object);
}

static void
_nc_powspec_mnl_halofit_finalize (GObject* object)
{
	NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (object);

	gsl_root_fdfsolver_free (pshf->priv->linear_scale_solver);
	gsl_root_fsolver_free (pshf->priv->znl_solver);

	/* Chain up : end */
	G_OBJECT_CLASS (nc_powspec_mnl_halofit_parent_class)
	->finalize (object);
}

static void _nc_powspec_mnl_halofit_prepare (NcmPowspec* powspec, NcmModel* model);
static gdouble _nc_powspec_mnl_halofit_eval (NcmPowspec* powspec, NcmModel* model, const gdouble z, const gdouble k);
static void _nc_powspec_mnl_halofit_eval_vec (NcmPowspec* powspec, NcmModel* model, const gdouble z, NcmVector* k, NcmVector* Pk);
static void _nc_powspec_mnl_halofit_get_nknots (NcmPowspec* powspec, guint* Nz, guint* Nk);

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

	g_type_class_add_private (klass, sizeof (NcPowspecMNLHaloFitPrivate));

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


	powspec_class->prepare = &_nc_powspec_mnl_halofit_prepare;
	powspec_class->eval = &_nc_powspec_mnl_halofit_eval;
	powspec_class->eval_vec = &_nc_powspec_mnl_halofit_eval_vec;
	powspec_class->get_nknots = &_nc_powspec_mnl_halofit_get_nknots;
}

typedef struct _int_var_moment_params
{
	const gint n;
	const gdouble R;
	const gdouble z;
	NcPowspecML* ps;
	NcHICosmo* cosmo;
} int_var_moment_params;

///////////////////// INTEGRATION OVER K*R ////////////////////////////
static gdouble
_nc_powspec_mnl_halofit_var_moment_integrand (gdouble kR, gpointer params)
{
	int_var_moment_params* ts = (int_var_moment_params*)params;
	const gdouble k = kR / ts->R;
	const gdouble matter_P = ncm_powspec_eval (NCM_POWSPEC (ts->ps), NCM_MODEL (ts->cosmo), ts->z, k);
	const gdouble kR2 = kR * kR;
	const gdouble W2 = exp (-kR2);
	// printf ("%g %g %g %g %g %i %g\n", k, ts->R, ts->z, matter_P, kR2, ts->n, W2);
	return matter_P * gsl_pow_int (kR2, ts->n + 1) * W2;
}

static gdouble
_nc_powspec_mnl_halofit_var_moment (NcPowspecML* ps, NcHICosmo* cosmo, const gdouble R, const gdouble z, const gint n)
{
	gdouble result, error;
	gsl_function F;
	int_var_moment_params ivmps = { n, R, z, ps, cosmo };

	gsl_integration_workspace* w = gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);

	F.function = &_nc_powspec_mnl_halofit_var_moment_integrand;
	F.params = &ivmps;

	gsl_integration_qagiu (&F, NCM_DEFAULT_PRECISION, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, w, &result, &error);

	gsl_integration_workspace_free (w);

	return result / (gsl_pow_3 (R) * ncm_c_2_pi_2 ());
}

typedef struct _var_params
{
	NcPowspecMNLHaloFit* pshf;
	NcHICosmo* cosmo;
	const gdouble z;
	const gdouble R_min;
} var_params;

static gdouble
_nc_powspec_mnl_halofit_varm1 (gdouble lnR, gpointer params)
{
	var_params* vps = (var_params*)params;
	/*
  const gdouble R = gsl_sf_exp (lnR);
  const gdouble sigma2_0 = _nc_powspec_mnl_halofit_var_moment (vps->pshf->psml, vps->cosmo, R, vps->z, 0);
  return sigma2_0 - 1.0;
*/
	return ncm_powspec_filter_eval_lnvar_lnr (vps->pshf->psml_gauss, vps->z, lnR);
}

static gdouble
_nc_powspec_mnl_halofit_varm1_deriv (gdouble lnR, gpointer params)
{
	var_params* vps = (var_params*)params;
	/*
  const gdouble R = gsl_sf_exp (lnR);
  const gdouble sigma2_1 = _nc_powspec_mnl_halofit_var_moment (vps->pshf->psml, vps->cosmo, R, vps->z, 1);
  return -2.0 * sigma2_1;
*/

	return ncm_powspec_filter_eval_dlnvar_dlnr (vps->pshf->psml_gauss, vps->z, lnR);
}

static void
_nc_powspec_mnl_halofit_varm1_fdf (gdouble lnR, gpointer params, gdouble* varm1, gdouble* dvarm1)
{
	var_params* vps = (var_params*)params;
	/*
  const gdouble R = gsl_sf_exp (lnR);
  const gdouble sigma2_0 = _nc_powspec_mnl_halofit_var_moment (vps->pshf->psml, vps->cosmo, R, vps->z, 0);
  const gdouble sigma2_1 = _nc_powspec_mnl_halofit_var_moment (vps->pshf->psml, vps->cosmo, R, vps->z, 1);
*/
	*varm1 = ncm_powspec_filter_eval_lnvar_lnr (vps->pshf->psml_gauss, vps->z, lnR);
	*dvarm1 = ncm_powspec_filter_eval_dlnvar_dlnr (vps->pshf->psml_gauss, vps->z, lnR);
}

typedef struct _root_params
{
	NcPowspecML* ps;
	NcHICosmo* cosmo;
} root_params;

static gdouble
_nc_powspec_mnl_halofit_linear_scale (NcPowspecMNLHaloFit* pshf, NcHICosmo* cosmo, const gdouble z)
{
	gint status;
	gint iter = 0, max_iter = 20000;

	gdouble lnR0 = 0.0;
	gdouble lnR = (-z / 2.0 < NC_POWSPEC_MNL_HALOFIT_LOGRMIN) ? NC_POWSPEC_MNL_HALOFIT_LOGRMIN : -z / 2.0 + NCM_DEFAULT_PRECISION;
	const gdouble reltol = pshf->reltol / 10.0;
	gdouble res = 0.0;

	gsl_function_fdf FDF;

	var_params vps = { pshf, cosmo, z, 0.0 };

	FDF.f = &_nc_powspec_mnl_halofit_varm1;
	FDF.df = &_nc_powspec_mnl_halofit_varm1_deriv;
	FDF.fdf = &_nc_powspec_mnl_halofit_varm1_fdf;
	FDF.params = &vps;

	gsl_root_fdfsolver_set (pshf->priv->linear_scale_solver, &FDF, lnR);

	do
	{
		iter++;
		status = gsl_root_fdfsolver_iterate (pshf->priv->linear_scale_solver);

		lnR0 = lnR;
		lnR = gsl_root_fdfsolver_root (pshf->priv->linear_scale_solver);

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

		gdouble lnRlo = -z / 2. - 10.;
		gdouble lnRup = -z / 2. + 10.;

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
	var_params* vps = (var_params*)params;

	return _nc_powspec_mnl_halofit_linear_scale (vps->pshf, vps->cosmo, z);
}

static gdouble
_nc_powspec_mnl_halofit_linear_scale_z_znl (gdouble z, gpointer params)
{
	var_params* vps = (var_params*)params;

	return _nc_powspec_mnl_halofit_linear_scale (vps->pshf, vps->cosmo, z) - vps->R_min;
}

static void
_nc_powspec_mnl_halofit_prepare_nl (NcPowspecMNLHaloFit* pshf, NcmModel* model)
{
	NcHICosmo* cosmo = NC_HICOSMO (model);
	guint i;

	pshf->priv->z = HUGE_VAL;

	ncm_powspec_require_zi (NCM_POWSPEC (pshf->psml), ncm_powspec_get_zi (NCM_POWSPEC (pshf)));
	ncm_powspec_require_zf (NCM_POWSPEC (pshf->psml), ncm_powspec_get_zf (NCM_POWSPEC (pshf)));

	ncm_powspec_require_kmin (NCM_POWSPEC (pshf->psml), ncm_powspec_get_kmin (NCM_POWSPEC (pshf)));
	ncm_powspec_require_kmax (NCM_POWSPEC (pshf->psml), ncm_powspec_get_kmax (NCM_POWSPEC (pshf)));

	ncm_powspec_filter_set_zi (pshf->psml_gauss, ncm_powspec_get_zi (NCM_POWSPEC (pshf)));
	ncm_powspec_filter_set_zf (pshf->psml_gauss, pshf->zmaxnl);

	ncm_powspec_filter_set_best_lnr0 (pshf->psml_gauss);
	ncm_powspec_filter_prepare_if_needed (pshf->psml_gauss, model);

	{
		const gdouble R_min = ncm_powspec_filter_get_r_min (pshf->psml_gauss);
		var_params vps = { pshf, cosmo, 0.0, R_min };

		gsl_function Fznl;

		Fznl.function = &_nc_powspec_mnl_halofit_linear_scale_z_znl;
		Fznl.params = &vps;

		if (_nc_powspec_mnl_halofit_linear_scale_z (pshf->zmaxnl, Fznl.params) > R_min)
			pshf->znl = pshf->zmaxnl;
		else
		{
			const gdouble step = 1.0;
			gdouble z0 = 0.0;
			gdouble z1 = z0 + step;
			gdouble R0 = _nc_powspec_mnl_halofit_linear_scale_z (z0, Fznl.params);
			gdouble R1 = _nc_powspec_mnl_halofit_linear_scale_z (z1, Fznl.params);
			gboolean found = FALSE;

			if (R0 < R_min)
				g_error ("_nc_powspec_mnl_halofit_prepare_nl: linear universe or too large R_min, in the latter case increase k_max (R0 == % 21.15g, R_min == % 21.15g).", R0, R_min);

			while (R1 > R_min)
			{
				z0 = z1;
				z1 += step;
				R1 = _nc_powspec_mnl_halofit_linear_scale_z (z1, Fznl.params);
				if (z1 > pshf->zmaxnl)
				{
					pshf->znl = pshf->zmaxnl;
					found = TRUE;
					break;
				}
			}

			if (!found)
			{
				gint status;
				gint iter = 0, max_iter = 20000;

				gsl_root_fsolver_set (pshf->priv->znl_solver, &Fznl, z0, z1);
				do
				{
					iter++;
					status = gsl_root_fsolver_iterate (pshf->priv->znl_solver);

					pshf->znl = gsl_root_fsolver_root (pshf->priv->znl_solver);
					z0 = gsl_root_fsolver_x_lower (pshf->priv->znl_solver);
					z1 = gsl_root_fsolver_x_upper (pshf->priv->znl_solver);
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

			ncm_spline_set_func (pshf->Rsigma, NCM_SPLINE_FUNCTION_SPLINE, &F, 0.0, pshf->znl, 0, pshf->reltol);
		}
	}

	{
		const guint len = ncm_vector_len (pshf->Rsigma->xv);
		NcmVector* neffv = ncm_vector_new (len);
		NcmVector* Curv = ncm_vector_new (len);

		for (i = 0; i < len; i++)
		{
			if (FALSE)
			{
				const gdouble z = ncm_vector_get (pshf->Rsigma->xv, i);
				const gdouble R = ncm_vector_get (pshf->Rsigma->yv, i);
				const gdouble sigma2_0 = _nc_powspec_mnl_halofit_var_moment (pshf->psml, cosmo, R, z, 0);
				const gdouble sigma2_1 = _nc_powspec_mnl_halofit_var_moment (pshf->psml, cosmo, R, z, 1);
				const gdouble sigma2_2 = _nc_powspec_mnl_halofit_var_moment (pshf->psml, cosmo, R, z, 2);
				const gdouble d1 = -2.0 * sigma2_1 / sigma2_0;
				const gdouble d2 = -d1 * d1 + 4.0 * (-sigma2_1 + sigma2_2) / sigma2_0;

				ncm_vector_set (neffv, i, -3.0 - d1);
				ncm_vector_set (Curv, i, -d2);

				printf ("# z = % 20.15g, R = % 20.15g | % 20.15g % 20.15g % 20.15g | %e %e %e\n",
				        z, R,
				        sigma2_0, d1, d2,
				        fabs ((ncm_powspec_filter_eval_var (pshf->psml_gauss, z, R) - sigma2_0) / sigma2_0),
				        fabs ((ncm_powspec_filter_eval_dlnvar_dlnr (pshf->psml_gauss, z, log (R)) - d1) / d1),
				        fabs ((ncm_powspec_filter_eval_dnlnvar_dlnrn (pshf->psml_gauss, z, log (R), 2) - d2) / d2));
			}
			else
			{
				const gdouble z = ncm_vector_get (pshf->Rsigma->xv, i);
				const gdouble R = ncm_vector_get (pshf->Rsigma->yv, i);
				const gdouble lnR = log (R);
				const gdouble d1 = ncm_powspec_filter_eval_dlnvar_dlnr (pshf->psml_gauss, z, lnR);
				const gdouble d2 = ncm_powspec_filter_eval_dnlnvar_dlnrn (pshf->psml_gauss, z, lnR, 2);

				ncm_vector_set (neffv, i, -3.0 - d1);
				ncm_vector_set (Curv, i, -d2);
			}
		}

		ncm_spline_set (pshf->neff, pshf->Rsigma->xv, neffv, TRUE);
		ncm_spline_set (pshf->Cur, pshf->Rsigma->xv, Curv, TRUE);

		ncm_vector_free (neffv);
		ncm_vector_free (Curv);
	}
}

static void
_nc_powspec_mnl_halofit_prepare (NcmPowspec* powspec, NcmModel* model)
{
	NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (powspec);

	g_assert (NC_IS_HICOSMO (model));

	ncm_powspec_prepare (NCM_POWSPEC (pshf->psml), model);

	_nc_powspec_mnl_halofit_prepare_nl (pshf, model);
}

static void
_nc_powspec_mnl_halofit_preeval (NcPowspecMNLHaloFit* pshf, NcHICosmo* cosmo, const gdouble z)
{
	const gdouble E2 = nc_hicosmo_E2 (cosmo, z);
	const gdouble Rsigma = ncm_spline_eval (pshf->Rsigma, z);

	const gdouble neff = ncm_spline_eval (pshf->neff, z);
	const gdouble Cur = ncm_spline_eval (pshf->Cur, z);

	const gdouble neff2 = neff * neff;
	const gdouble neff3 = neff2 * neff;
	const gdouble neff4 = neff2 * neff2;

	const gdouble Omega_de_onepp = NC_IS_HICOSMO_DE (cosmo) ? nc_hicosmo_de_E2Omega_de_onepw (NC_HICOSMO_DE (cosmo), z) / E2 : 0.0;
	const gdouble Omega_m = nc_hicosmo_Omega_m0 (cosmo) * gsl_pow_3 (1.0 + z) / E2;

	pshf->priv->z = z;

	pshf->priv->ksigma = 1.0 / Rsigma;

	pshf->priv->an = ncm_util_exp10 (1.5222 + 2.8553 * neff + 2.3706 * neff2 + 0.9903 * neff3 + 0.2250 * neff4 - 0.6038 * Cur + 0.1749 * Omega_de_onepp);
	pshf->priv->bn = ncm_util_exp10 (-0.5642 + 0.5864 * neff + 0.5716 * neff2 - 1.5474 * Cur + 0.2279 * Omega_de_onepp);
	pshf->priv->cn = ncm_util_exp10 (0.3698 + 2.0404 * neff + 0.8161 * neff2 + 0.5869 * Cur);
	pshf->priv->gamman = 0.1971 - 0.0843 * neff + 0.8460 * Cur;
	pshf->priv->alphan = fabs (6.0835 + 1.3373 * neff - 0.1959 * neff2 - 5.5274 * Cur);
	pshf->priv->betan = 2.0379 - 0.7354 * neff + 0.3157 * neff2 + 1.2490 * neff3 + 0.3980 * neff4 - 0.1682 * Cur; // + fnu*(1.081 + 0.395*pow(rneff,2)
	pshf->priv->nun = ncm_util_exp10 (5.2105 + 3.6902 * neff);

	pshf->priv->f1 = pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F1POW);
	pshf->priv->f2 = pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F2POW);
	pshf->priv->f3 = pow (Omega_m, NC_POWSPEC_MNL_HALOFIT_F3POW);
}

static gdouble
_nc_powspec_mnl_halofit_eval (NcmPowspec* powspec, NcmModel* model, const gdouble z, const gdouble k)
{
	NcHICosmo* cosmo = NC_HICOSMO (model);
	NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (powspec);
	const gboolean linscale = (z > pshf->znl);
	const gboolean applysmooth = (z + 1.0 > pshf->znl);
	const gdouble zhf = linscale ? pshf->znl : z;
	const gdouble Pklin = ncm_powspec_eval (NCM_POWSPEC (pshf->psml), model, z, k);
	gdouble Pknln;

	if (zhf != pshf->priv->z)
	{
		_nc_powspec_mnl_halofit_preeval (pshf, cosmo, zhf);
	}
	{
		const gdouble k3 = gsl_pow_3 (k);
		const gdouble k3o2pi2 = k3 / ncm_c_2_pi_2 ();
		const gdouble Delta_lin = k3o2pi2 * Pklin;
		const gdouble y = k / pshf->priv->ksigma;
		const gdouble P_Q = Pklin * (pow (1.0 + Delta_lin, pshf->priv->betan) / (1.0 + pshf->priv->alphan * Delta_lin)) * exp (-y / 4.0 - y * y / 8.0);

		const gdouble Delta_Hprime = pshf->priv->an * pow (y, 3.0 * pshf->priv->f1) / (1.0 + pshf->priv->bn * pow (y, pshf->priv->f2) + pow (pshf->priv->cn * pshf->priv->f3 * y, 3.0 - pshf->priv->gamman));
		const gdouble Delta_H = Delta_Hprime / (1.0 + pshf->priv->nun / (y * y));

		const gdouble P_H = Delta_H / k3o2pi2;

		Pknln = P_Q + P_H;
	}

	if (applysmooth)
		Pknln = ncm_util_smooth_trans (Pknln, Pklin, pshf->znl, 1.0, z);

	return Pknln;
}

static void
_nc_powspec_mnl_halofit_eval_vec (NcmPowspec* powspec, NcmModel* model, const gdouble z, NcmVector* k, NcmVector* Pk)
{
	NcHICosmo* cosmo = NC_HICOSMO (model);
	NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (powspec);
	const gboolean linscale = (z > pshf->znl);
	const gboolean applysmooth = (z + 1.0 > pshf->znl);
	const gdouble zhf = linscale ? pshf->znl : z;
	gdouble theta0, theta1;

	ncm_powspec_eval_vec (NCM_POWSPEC (pshf->psml), model, z, k, Pk);

	if (applysmooth)
		ncm_util_smooth_trans_get_theta (pshf->znl, 1.0, z, &theta0, &theta1);

	if (zhf != pshf->priv->z)
	{
		_nc_powspec_mnl_halofit_preeval (pshf, cosmo, zhf);
	}

	{
		const guint len = ncm_vector_len (k);
		guint i;

		for (i = 0; i < len; i++)
		{
			const gdouble ki = ncm_vector_get (k, i);
			const gdouble ki3 = gsl_pow_3 (ki);
			const gdouble ki3o2pi2 = ki3 / ncm_c_2_pi_2 ();
			const gdouble Pklin = ncm_vector_get (Pk, i);
			const gdouble Delta_lin = ki3o2pi2 * Pklin;
			const gdouble y = ki / pshf->priv->ksigma;

			const gdouble P_Q = Pklin * (pow (1.0 + Delta_lin, pshf->priv->betan) / (1.0 + pshf->priv->alphan * Delta_lin)) * exp (-y / 4.0 - y * y / 8.0);
			const gdouble Delta_Hprime = pshf->priv->an * pow (y, 3.0 * pshf->priv->f1) / (1.0 + pshf->priv->bn * pow (y, pshf->priv->f2) + pow (pshf->priv->cn * pshf->priv->f3 * y, 3.0 - pshf->priv->gamman));
			const gdouble Delta_H = Delta_Hprime / (1.0 + pshf->priv->nun / (y * y));
			const gdouble P_H = Delta_H / ki3o2pi2;

			const gdouble Pknln = P_Q + P_H;

			if (applysmooth)
				ncm_vector_set (Pk, i, theta0 * Pknln + theta1 * Pklin);
			else
				ncm_vector_set (Pk, i, Pknln);
		}
	}
}

static void
_nc_powspec_mnl_halofit_get_nknots (NcmPowspec* powspec, guint* Nz, guint* Nk)
{
	NcPowspecMNLHaloFit* pshf = NC_POWSPEC_MNL_HALOFIT (powspec);
	ncm_powspec_get_nknots (NCM_POWSPEC (pshf->psml), Nz, Nk);
	*Nz = ncm_vector_len (pshf->Rsigma->xv);
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
nc_powspec_mnl_halofit_new (NcPowspecML* psml, gdouble zmaxnl, gdouble reltol)
{
	NcPowspecMNLHaloFit* pshf = g_object_new (NC_TYPE_POWSPEC_MNL_HALOFIT,
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
void nc_powspec_mnl_halofit_set_kbounds_from_ml (NcPowspecMNLHaloFit* pshf)
{
	ncm_powspec_set_kmin (NCM_POWSPEC (pshf), ncm_powspec_get_kmin (NCM_POWSPEC (pshf->psml)));
	ncm_powspec_set_kmax (NCM_POWSPEC (pshf), ncm_powspec_get_kmax (NCM_POWSPEC (pshf->psml)));
}
