/***************************************************************************
 *            nc_xcor.h
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
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
 * SECTION:nc_xcor
 * @title: Cross-correlations
 * @short_description: Cross-spectra using the Limber approximation
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/integral.h"
#include "math/memory_pool.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "lss/nc_window_tophat.h"
#include "lss/nc_matter_var.h"
#include "xcor/nc_xcor.h"

enum
{
	PROP_0,
	PROP_DISTANCE,
	PROP_TRANSFER_FUNC,
	PROP_GROWTH_FUNC,
	PROP_ZL,
	PROP_ZU,
};

G_DEFINE_TYPE (NcXcor, nc_xcor, G_TYPE_OBJECT);

/**
 * nc_xcor_new:
 * @dist: a #NcDistance
 * @tf: a #NcTransferFunc
 * @gf: a #NcGrowthFunc
 * @zl: a #gdouble
 * @zu: a #gdouble
 *
 * FIXME
 *
 * Returns: FIXME
 *
*/
NcXcor* nc_xcor_new (NcDistance* dist, NcTransferFunc* tf, NcGrowthFunc* gf, gdouble zl, gdouble zu)
{
	return g_object_new (NC_TYPE_XCOR, "distance", dist, "transfer-func", tf, "growth-func", gf, "zu", zu, "zl", zl, NULL);
}

/**
 * nc_xcor_ref:
 * @xcl: a #NcXcor.
 *
 * FIXME
 *
 * Returns: (transfer full): @xcl.
 */
NcXcor* nc_xcor_ref (NcXcor* xcl)
{
	return g_object_ref (xcl);
}

/**
 * nc_xcor_free:
 * @xcl: a #NcXcor.
 *
 * FIXME
 *
 */
void nc_xcor_free (NcXcor* xcl)
{
	g_object_unref (xcl);
}

/**
 * nc_xcor_clear:
 * @xcl: a #NcXcor.
 *
 * FIXME
 *
 */
void nc_xcor_clear (NcXcor** xcl)
{
	g_clear_object (xcl);
}

void nc_xcor_prepare (NcXcor *xc, NcHIReion *reion, NcHICosmo *cosmo)
{
	nc_transfer_func_prepare (xc->tf, reion, cosmo);
	nc_growth_func_prepare (xc->gf, cosmo);
	nc_distance_prepare_if_needed (xc->dist, cosmo);

	NcWindow* wp = nc_window_tophat_new ();
	NcMatterVar* var = nc_matter_var_new (NC_MATTER_VAR_FFT, wp, xc->tf);
	gdouble ln8_hMpc = log (8.0);
	gdouble sigma82_no_norm = nc_matter_var_var0 (var, cosmo, ln8_hMpc);
	gdouble sigma82 = pow (nc_hicosmo_sigma_8 (cosmo), 2.0);
	xc->normPS = sigma82 / sigma82_no_norm;
}


typedef struct _xcor_limber_cross_cl_int
{
	NcHICosmo* cosmo;
	NcTransferFunc* tf;
	NcDistance* dist;
	NcGrowthFunc* gf;

	NcXcorLimber* xcl1;
	NcXcorLimber* xcl2;
	gint l;

} xcor_limber_cross_cl_int;

static gdouble _xcor_limber_cross_cl_int_z (gdouble z, gpointer ptr)
{
	xcor_limber_cross_cl_int* xcli = (xcor_limber_cross_cl_int*)ptr;

	gdouble k1z = nc_xcor_limber_eval_kernel (xcli->xcl1, xcli->cosmo, z, xcli->l);
	gdouble k2z = nc_xcor_limber_eval_kernel (xcli->xcl2, xcli->cosmo, z, xcli->l);

	gdouble xi_z = nc_distance_comoving (xcli->dist, xcli->cosmo, z); // dimensionless, ie it's in unit of hubble radius
	gdouble kh = xcli->l / (xi_z * ncm_c_hubble_radius ()); // in h Mpc-1, hubble
	gdouble power_spec_init = nc_hicosmo_powspec (xcli->cosmo, kh); // k^ns
	gdouble t_k = nc_transfer_func_eval (xcli->tf, xcli->cosmo, kh);
	gdouble D_z = nc_growth_func_eval (xcli->gf, xcli->cosmo, z);
	gdouble power_spec = power_spec_init * t_k * t_k * D_z * D_z;

	gdouble E_z = nc_hicosmo_E (xcli->cosmo, z);

	return E_z * k1z * k2z * power_spec / (xi_z * xi_z);
}

/**
 * nc_xcor_limber_cross_cl:
 * @xc: a #NcXcor
 * @xcl1: a #NcXcorLimber
 * @xcl2: a #NcXcorLimber
 * @reion: a #NcHIReion
 * @cosmo: a #NcHICosmo
 * @ell: a #NcmVector
 * @vp: a #NcmVector
 * @lmin_idx: a #guint
 *
 * FIXME
 *
 * Returns: FIXME
 *
*/
void 
nc_xcor_limber_cross_cl (NcXcor* xc, NcXcorLimber* xcl1, NcXcorLimber* xcl2, NcHIReion *reion, NcHICosmo* cosmo, NcmVector* ell, NcmVector* vp, guint lmin_idx)
{
	guint nell = ncm_vector_len (ell);
	guint nvp = ncm_vector_len (vp);
	if (nvp < nell + lmin_idx)
	{
		g_error ("vector vp too short");
	}

	gboolean xcl1_up = ncm_model_ctrl_update (xc->ctrl1, NCM_MODEL (xcl1));
	gboolean xcl2_up = ncm_model_ctrl_update (xc->ctrl2, NCM_MODEL (xcl2));
	gboolean cosmo_up = ncm_model_ctrl_update (xc->ctrlcosmo, NCM_MODEL (cosmo));

	if (xcl1_up || xcl2_up || cosmo_up)
	{
		nc_xcor_limber_prepare (xcl1, cosmo);
		nc_xcor_limber_prepare (xcl2, cosmo);
		nc_xcor_prepare (xc, reion, cosmo);
	}

	gdouble cl, cons_factor, err;
	xcor_limber_cross_cl_int xcli;
	gsl_function F;
	gsl_integration_workspace** w = ncm_integral_get_workspace ();

	xcli.xcl1 = xcl1;
	xcli.xcl2 = xcl2;
	xcli.cosmo = cosmo;
	xcli.dist = xc->dist;
	// xcli.l     = l;
	xcli.tf = xc->tf;
	xcli.gf = xc->gf;


	F.function = &_xcor_limber_cross_cl_int_z;
	F.params = &xcli;

	// H0/c in h Mpc-1:
	gdouble H0_c = 1e5 / ncm_c_c ();
	cons_factor = xcl1->cons_factor * xcl2->cons_factor * pow (H0_c, 3.0); // power spectrum is in (h^-1 Mpc)^3 integrand.xi_lss = nc_distance_comoving_lss (dist, model);

	gdouble r = 0.0;
	guint i;

	for (i = 0; i < nell; i++)
	{
		xcli.l = ncm_vector_get (ell, i);
		gsl_integration_qag (&F, xc->zl, xc->zu, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &cl, &err);
		r = cl * cons_factor * xc->normPS;
		ncm_vector_set (vp, i + lmin_idx, r);
	}

	ncm_memory_pool_return (w);
}

typedef struct _xcor_limber_auto_cl_int
{
	NcHICosmo* cosmo;
	NcTransferFunc* tf;
	NcDistance* dist;
	NcGrowthFunc* gf;

	NcXcorLimber* xcl;
	gint l;

} xcor_limber_auto_cl_int;

static gdouble _xcor_limber_auto_cl_int_z (gdouble z, gpointer ptr)
{
	xcor_limber_auto_cl_int* xcli = (xcor_limber_auto_cl_int*)ptr;

	gdouble k1z = nc_xcor_limber_eval_kernel (xcli->xcl, xcli->cosmo, z, xcli->l);

	gdouble xi_z = nc_distance_comoving (xcli->dist, xcli->cosmo, z); // dimensionless, ie it's in unit of hubble radius
	gdouble kh = (xcli->l + 0.5) / (xi_z * ncm_c_hubble_radius ()); // in h Mpc-1, hubble
	// radius = c[m/s] / (the +0.5 is not to be forgotten, cf. LoVerde, M., & Afshordi, N. (2008). Extended Limber approximation. Physical Review D, 78(1), 123506. http://doi.org/10.1103/PhysRevD.78.123506)
	gdouble power_spec_init = nc_hicosmo_powspec (xcli->cosmo, kh); // k^ns
	gdouble t_k = nc_transfer_func_eval (xcli->tf, xcli->cosmo, kh);
	gdouble D_z = nc_growth_func_eval (xcli->gf, xcli->cosmo, z);
	gdouble power_spec = power_spec_init * t_k * t_k * D_z * D_z;

	gdouble E_z = nc_hicosmo_E (xcli->cosmo, z);

	return E_z * k1z * k1z * power_spec / (xi_z * xi_z);
}

/**
 * nc_xcor_limber_auto_cl:
 * @xc: a #NcXcor
 * @xcl: a #NcXcorLimber
 * @reion: a #NcHIReion
 * @cosmo: a #NcHICosmo
 * @ell: a #NcmVector
 * @vp: a #NcmVector
 * @lmin_idx: a #guint
 * @withnoise: a #gboolean
 *
 * FIXME
 *
 * Returns: FIXME
 *
*/
void 
nc_xcor_limber_auto_cl (NcXcor* xc, NcXcorLimber* xcl, NcHIReion *reion, NcHICosmo* cosmo, NcmVector* ell, NcmVector* vp, guint lmin_idx, gboolean withnoise)
{
	guint nell = ncm_vector_len (ell);
	guint nvp = ncm_vector_len (vp);
	if (nvp < nell + lmin_idx)
	{
		g_error ("vector vp too short");
	}

	gboolean xcl_up = ncm_model_ctrl_update (xc->ctrl1, NCM_MODEL (xcl));
	gboolean cosmo_up = ncm_model_ctrl_update (xc->ctrlcosmo, NCM_MODEL (cosmo));

	if (xcl_up || cosmo_up)
	{
		nc_xcor_limber_prepare (xcl, cosmo);
		nc_xcor_prepare (xc, reion, cosmo);
	}

	gdouble cl, cons_factor, err;
	xcor_limber_auto_cl_int xcli;
	gsl_function F;
	gsl_integration_workspace** w = ncm_integral_get_workspace ();

	xcli.xcl = xcl;
	xcli.cosmo = cosmo;
	xcli.dist = xc->dist;
	xcli.tf = xc->tf;
	xcli.gf = xc->gf;


	F.function = &_xcor_limber_auto_cl_int_z;
	F.params = &xcli;

	// H0/c in h Mpc-1:
	gdouble H0_c = 1e5 / ncm_c_c ();
	cons_factor = pow (xcl->cons_factor, 2.0) * pow (H0_c, 3.0); // power spectrum is in (h^-1 Mpc)^3 integrand.xi_lss = nc_distance_comoving_lss (dist, model);

	gdouble r = 0.0;
	guint i, l;

	for (i = 0; i < nell; i++)
	{
		l = ncm_vector_get (ell, i);
		xcli.l = l;
		gsl_integration_qag (&F, xc->zl, xc->zu, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &cl, &err);
		r = cl * cons_factor * xc->normPS;
		if (withnoise)
		{
			r += nc_xcor_limber_noise_spec (xcl, l);
		}
		ncm_vector_set (vp, i + lmin_idx, r);
	}

	ncm_memory_pool_return (w);
}

static void nc_xcor_init (NcXcor* xc)
{
	xc->ctrlcosmo = ncm_model_ctrl_new (NULL);
	xc->ctrl1 = ncm_model_ctrl_new (NULL);
	xc->ctrl2 = ncm_model_ctrl_new (NULL);
	xc->tf = NULL;
	xc->gf = NULL;
	xc->dist = NULL;
	xc->zl = 0.0;
	xc->zu = 0.0;
	xc->normPS = 0.0;
}

static void _nc_xcor_dispose (GObject* object)
{
	NcXcor* xc = NC_XCOR (object);

	ncm_model_ctrl_clear (&xc->ctrlcosmo);
	ncm_model_ctrl_clear (&xc->ctrl1);
	ncm_model_ctrl_clear (&xc->ctrl2);


	/* Chain up : end */
	G_OBJECT_CLASS (nc_xcor_parent_class)
	->dispose (object);
}

static void _nc_xcor_finalize (GObject* object)
{

	/* Chain up : end */
	G_OBJECT_CLASS (nc_xcor_parent_class)
	->finalize (object);
}

static void
_nc_xcor_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
	NcXcor* xc = NC_XCOR (object);
	g_return_if_fail (NC_IS_XCOR (object));

	switch (prop_id)
	{
	case PROP_DISTANCE:
		xc->dist = g_value_get_object (value);
		break;
	case PROP_TRANSFER_FUNC:
		xc->tf = g_value_get_object (value);
		break;
	case PROP_GROWTH_FUNC:
		xc->gf = g_value_get_object (value);
		break;
	case PROP_ZL:
		xc->zl = g_value_get_double (value);
		break;
	case PROP_ZU:
		xc->zu = g_value_get_double (value);
		break;
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
_nc_xcor_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
	NcXcor* xc = NC_XCOR (object);
	g_return_if_fail (NC_IS_XCOR (object));

	switch (prop_id)
	{
	case PROP_DISTANCE:
		g_value_set_object (value, xc->dist);
		break;
	case PROP_TRANSFER_FUNC:
		g_value_set_object (value, xc->tf);
		break;
	case PROP_GROWTH_FUNC:
		g_value_set_object (value, xc->gf);
		break;
	case PROP_ZL:
		g_value_set_double (value, xc->zl);
		break;
	case PROP_ZU:
		g_value_set_double (value, xc->zu);
		break;
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
nc_xcor_class_init (NcXcorClass* klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	//GObjectClass* parent_class = G_OBJECT_CLASS (klass);

	object_class->dispose = _nc_xcor_dispose;
	object_class->finalize = _nc_xcor_finalize;
	object_class->set_property = _nc_xcor_set_property;
	object_class->get_property = _nc_xcor_get_property;

	/**
   * NcXcor:distance:
   *
   * This property keeps the distance object.
   */
	g_object_class_install_property (object_class,
	                                 PROP_DISTANCE,
	                                 g_param_spec_object ("distance",
	                                                      NULL,
	                                                      "Distance.",
	                                                      NC_TYPE_DISTANCE,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	/**
   * NcXcor:transfer_func:
   *
   * This property keeps the transfer function object.
   */
	g_object_class_install_property (object_class,
	                                 PROP_TRANSFER_FUNC,
	                                 g_param_spec_object ("transfer-func",
	                                                      NULL,
	                                                      "Transfer Function.",
	                                                      NC_TYPE_TRANSFER_FUNC,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	/**
   * NcXcor:growth:
   *
   * This property keeps the growth function object.
   */
	g_object_class_install_property (object_class,
	                                 PROP_GROWTH_FUNC,
	                                 g_param_spec_object ("growth-func",
	                                                      NULL,
	                                                      "Growth function.",
	                                                      NC_TYPE_GROWTH_FUNC,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	/**
	* NcXcor:zl:
	*
	* This property sets lower redshift bound.
	*/
	g_object_class_install_property (object_class,
	                                 PROP_ZL,
	                                 g_param_spec_double ("zl",
	                                                      NULL,
	                                                      "Lower redshift integration bound",
	                                                      0.0, G_MAXDOUBLE, 0.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	/**
	* NcXcor:zu:
	*
	* This property sets the upper redshift bound.
	*/
	g_object_class_install_property (object_class,
	                                 PROP_ZU,
	                                 g_param_spec_double ("zu",
	                                                      NULL,
	                                                      "Upper",
	                                                      0.0, G_MAXDOUBLE, 0.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}
