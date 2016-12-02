/***************************************************************************
 *            nc_xcor_limber_kernel_weak_lensing.c
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
 * SECTION:nc_xcor_limber_kernel_weak_lensing
 * @title: Cross-correlations Galaxy Kernel
 * @short_description: Galaxy implementation of NcNcXcorLimberKernel
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_cfg.h"
#include "xcor/nc_xcor.h"
#include "xcor/nc_xcor_limber_kernel_weak_lensing.h"

#include "math/integral.h"
#include "math/ncm_spline_gsl.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "math/ncm_spline_func.h"

G_DEFINE_TYPE (NcXcorLimberKernelWeakLensing, nc_xcor_limber_kernel_weak_lensing, NC_TYPE_XCOR_LIMBER_KERNEL);

#define VECTOR (NCM_MODEL (xclkg)->params)
#define MAG_BIAS (ncm_vector_get (VECTOR, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_MAG_BIAS))
#define NOISE_BIAS (ncm_vector_get (VECTOR, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_NOISE_BIAS))

enum
{
	PROP_0,
	PROP_DN_DZ,
	PROP_NBAR,
	PROP_INTR_SHEAR,
	PROP_DIST,
	PROP_SIZE,
};

static void
nc_xcor_limber_kernel_weak_lensing_init (NcXcorLimberKernelWeakLensing* xclkg)
{
	xclkg->dn_dz = NULL;

	xclkg->dist = NULL;

	xclkg->src_int = NULL;

	xclkg->nbar = 0.0;
	xclkg->intr_shear = 0.0;

	xclkg->noise = 0.0;
}

static void
_nc_xcor_limber_kernel_weak_lensing_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
	NcXcorLimberKernelWeakLensing* xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (object);
	g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_WEAK_LENSING (object));

	switch (prop_id)
	{
	case PROP_DN_DZ:
		xclkg->dn_dz = g_value_dup_object (value);
		break;
	case PROP_NBAR:
		xclkg->nbar = g_value_get_double (value);
		break;
	case PROP_INTR_SHEAR:
		xclkg->intr_shear = g_value_get_double (value);
		break;
	case PROP_DIST:
		xclkg->dist = g_value_dup_object (value);
		break;
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
_nc_xcor_limber_kernel_weak_lensing_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
	NcXcorLimberKernelWeakLensing* xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (object);
	g_return_if_fail (NC_IS_XCOR_LIMBER_KERNEL_WEAK_LENSING (object));

	switch (prop_id)
	{
	case PROP_DN_DZ:
		g_value_set_object (value, xclkg->dn_dz);
		break;
	case PROP_NBAR:
		g_value_set_double (value, xclkg->nbar);
		break;
	case PROP_INTR_SHEAR:
		g_value_set_double (value, xclkg->intr_shear);
		break;
	case PROP_DIST:
		g_value_set_object (value, xclkg->dist);
		break;

	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
_nc_xcor_limber_kernel_weak_lensing_constructed (GObject* object)
{
	/* Chain up : start */
	G_OBJECT_CLASS (nc_xcor_limber_kernel_weak_lensing_parent_class)->constructed (object);
	{
		NcXcorLimberKernelWeakLensing* xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (object);
		NcXcorLimberKernel* xclk = NC_XCOR_LIMBER_KERNEL (xclkg);
		// NcmModel* model = NCM_MODEL (xclkg);
		// guint i;

		// Initialize g function spline for magnification bias
		// if (xclkg->domagbias)
		// {
		// 	NcmVector* gzv;
		// 	ncm_spline_clear (&xclkg->g_func);
		//
		// 	xclkg->g_func = ncm_spline_cubic_notaknot_new ();
		// 	ncm_spline_set_len (xclkg->g_func, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_G_FUNC_LEN);
		// 	gzv = ncm_spline_get_xv (xclkg->g_func);
		//
		// 	for (i = 0; i < NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_G_FUNC_LEN; i++)
		// 	{
		// 		ncm_vector_set (gzv, i, xclk->zmin + (xclk->zmax - xclk->zmin) / (NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_G_FUNC_LEN - 1) * i);
		// 	}
		//
		// 	ncm_vector_free (gzv);
		// 	/*ncm_spline_acc (xclkg->g_func, TRUE);*/
		// }

		// Initialize the source integral
		xclkg->src_int = ncm_spline_cubic_notaknot_new ();

		// Normalize the redshift distribution
		ncm_spline_prepare (xclkg->dn_dz);
		gdouble ngal = ncm_spline_eval_integ (xclkg->dn_dz, xclk->zmin, xclk->zmax);
		NcmVector* yv = ncm_spline_get_yv (xclkg->dn_dz);

		ncm_vector_scale (yv, 1.0 / ngal);
		ncm_spline_prepare (xclkg->dn_dz);
		ncm_vector_free (yv);

		// Prepare the bias spline and link it to the bias parameter vector
		// NcmVector *zv, *bv;
		// guint bvi;
		// guint bz_size = ncm_model_vparam_len (model, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_BIAS);
		// xclkg->nknots = bz_size;
		//
		// bvi = ncm_model_vparam_index (model, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_BIAS, 0);
		// zv = ncm_vector_new (bz_size);
		// bv = ncm_vector_get_subvector (model->params, bvi, bz_size);
		//
		// switch (bz_size)
		// {
		// case 1:
		// 	xclkg->bias_spline = NULL;
		// 	xclkg->bias = ncm_vector_ptr (bv, 0);
		// 	break;
		// case 2:
		// 	ncm_vector_set (zv, 0, xclk->zmin);
		// 	ncm_vector_set (zv, 1, xclk->zmax);
		// 	xclkg->bias_spline = ncm_spline_gsl_new_full (gsl_interp_linear, zv, bv, FALSE);
		// 	break;
		// default:
		// {
		// 	for (i = 0; i < bz_size; i++)
		// 	{
		// 		gdouble zi = xclk->zmin + (xclk->zmax - xclk->zmin) / (bz_size - 1) * i;
		// 		ncm_vector_set (zv, i, zi);
		// 	}
		// 	xclkg->bias_spline = ncm_spline_gsl_new_full (gsl_interp_polynomial, zv, bv, FALSE);
		// 	break;
		// }
		// }
		// ncm_vector_free (bv);
		// ncm_vector_free (zv);

		// If bias is constant, it's just a multiplicative factor, so possible to accelerate recomputation of C_l's
		// xclkg->fast_update = (bz_size == 1) && !(xclkg->domagbias);

		// Noise level
		xclkg->noise = gsl_pow_2 (xclkg->intr_shear) / xclkg->nbar;
	}
}

static void
_nc_xcor_limber_kernel_weak_lensing_dispose (GObject* object)
{
	NcXcorLimberKernelWeakLensing* xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (object);

	ncm_spline_clear (&xclkg->dn_dz);
	ncm_spline_clear (&xclkg->src_int);

	nc_distance_clear (&xclkg->dist);

	/* Chain up : end */
	G_OBJECT_CLASS (nc_xcor_limber_kernel_weak_lensing_parent_class)->dispose (object);
}

static void
_nc_xcor_limber_kernel_weak_lensing_finalize (GObject* object)
{

	/* Chain up : end */
	G_OBJECT_CLASS (nc_xcor_limber_kernel_weak_lensing_parent_class)->finalize (object);
}

static gdouble _nc_xcor_limber_kernel_weak_lensing_eval (NcXcorLimberKernel* xclk, NcHICosmo* cosmo, gdouble z, const NcXcorKinetic* xck, gint l);
static void _nc_xcor_limber_kernel_weak_lensing_prepare (NcXcorLimberKernel* xclk, NcHICosmo* cosmo);
/*static gdouble _nc_xcor_limber_kernel_weak_lensing_noise_spec (NcXcorLimberKernel* xclk, guint l);*/
static void _nc_xcor_limber_kernel_weak_lensing_add_noise (NcXcorLimberKernel* xclk, NcmVector* vp1, NcmVector* vp2, guint lmin);
guint _nc_xcor_limber_kernel_weak_lensing_obs_len (NcXcorLimberKernel* xclk);
guint _nc_xcor_limber_kernel_weak_lensing_obs_params_len (NcXcorLimberKernel* xclk);

static void
nc_xcor_limber_kernel_weak_lensing_class_init (NcXcorLimberKernelWeakLensingClass* klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcXcorLimberKernelClass* parent_class = NC_XCOR_LIMBER_KERNEL_CLASS (klass);
	NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

	model_class->set_property = &_nc_xcor_limber_kernel_weak_lensing_set_property;
	model_class->get_property = &_nc_xcor_limber_kernel_weak_lensing_get_property;
	object_class->constructed = &_nc_xcor_limber_kernel_weak_lensing_constructed;
	object_class->dispose = &_nc_xcor_limber_kernel_weak_lensing_dispose;
	object_class->finalize = &_nc_xcor_limber_kernel_weak_lensing_finalize;

	ncm_model_class_set_name_nick (model_class, "Xcor limber weak lensing", "Xcor-WL");
	ncm_model_class_add_params (model_class, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_SPARAM_LEN, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_VPARAM_LEN, PROP_SIZE);

	// g_object_class_install_property (object_class,
	//                                  PROP_Z_MIN,
	//                                  g_param_spec_double ("zmin",
	//                                                       NULL,
	//                                                       "Minimum redshift",
	//                                                       0.0, 20.0, 0.0,
	//                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	// g_object_class_install_property (object_class,
	//                                  PROP_Z_MAX,
	//                                  g_param_spec_double ("zmax",
	//                                                       NULL,
	//                                                       "Maximum redshift",
	//                                                       0.0, 20.0, 10.0,
	//                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_DN_DZ,
	                                 g_param_spec_object ("dndz",
	                                                      NULL,
	                                                      "Source redshift distribution",
	                                                      NCM_TYPE_SPLINE,
	                                                      G_PARAM_CONSTRUCT | G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); //G_PARAM_CONSTRUCT_ONLY |

	g_object_class_install_property (object_class,
	                                 PROP_DIST,
	                                 g_param_spec_object ("dist",
	                                                      NULL,
	                                                      "Distance object",
	                                                      NC_TYPE_DISTANCE,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


	g_object_class_install_property (object_class,
	                                 PROP_NBAR,
	                                 g_param_spec_double ("nbar",
	                                                      NULL,
	                                                      "nbar (galaxy angular density)",
	                                                      0.0, GSL_POSINF, 0.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_INTR_SHEAR,
	                                 g_param_spec_double ("intr-shear",
	                                                      NULL,
	                                                      "Intrinsic galaxy shear",
	                                                      0.0, 20.0, 0.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


	/* Check for errors in parameters initialization */
	ncm_model_class_check_params_info (model_class);

	parent_class->eval = &_nc_xcor_limber_kernel_weak_lensing_eval;
	parent_class->prepare = &_nc_xcor_limber_kernel_weak_lensing_prepare;
	//parent_class->noise_spec     = &_nc_xcor_limber_kernel_weak_lensing_noise_spec;
	parent_class->add_noise = &_nc_xcor_limber_kernel_weak_lensing_add_noise;

	parent_class->obs_len = &_nc_xcor_limber_kernel_weak_lensing_obs_len;
	parent_class->obs_params_len = &_nc_xcor_limber_kernel_weak_lensing_obs_params_len;

	parent_class->impl = NC_XCOR_LIMBER_KERNEL_IMPL_ALL;
}

/**
 * nc_xcor_limber_kernel_weak_lensing_new:
 * @zmin: a gdouble
 * @zmax: a gdouble
 * @dn_dz: a #NcmSpline
 * @nbar: a gdouble, gal density
 * @intr_shear: a gdouble, intrinsic galaxy shear
 * @dist: a #NcDistance
 *
 * Returns: a #NcXcorLimberKernelWeakLensing
 */
NcXcorLimberKernelWeakLensing*
nc_xcor_limber_kernel_weak_lensing_new (gdouble zmin, gdouble zmax, NcmSpline* dn_dz, gdouble nbar, gdouble intr_shear, NcDistance* dist)
{
	NcXcorLimberKernelWeakLensing* xclkg = g_object_new (NC_TYPE_XCOR_LIMBER_KERNEL_WEAK_LENSING,
	                                                     "zmin", zmin,
	                                                     "zmax", zmax,
	                                                     "dndz", dn_dz,
	                                                     "nbar", nbar,
	                                                     "intr-shear", intr_shear,
	                                                     "dist", dist,
	                                                     NULL);
	return xclkg;
}


typedef struct _int_src_int_params
{
	gdouble chiz;
	NcDistance* dist;
	NcHICosmo* cosmo;
	NcmSpline* dn_dz;
} int_src_int_params;

static gdouble
_nc_xcor_limber_kernel_weak_lensing_src_int_integrand (gdouble zz, gpointer params)
{
	int_src_int_params* ts = (int_src_int_params*)params;
	const gdouble a = 1.0 - ts->chiz / nc_distance_comoving (ts->dist, ts->cosmo, zz);
	const gdouble dn_dz_zz = ncm_spline_eval (ts->dn_dz, zz);

	return a * dn_dz_zz;
}

typedef struct _src_int_params
{
	NcXcorLimberKernelWeakLensing* xclkg;
	NcHICosmo* cosmo;
} src_int_params;

static gdouble _nc_xcor_limber_kernel_weak_lensing_src_int (gdouble z, gpointer params)
{
	src_int_params* ts = (src_int_params*)params;
	NcXcorLimberKernelWeakLensing* xclkg = ts->xclkg;
	NcXcorLimberKernel* xclk = NC_XCOR_LIMBER_KERNEL (xclkg);
	NcHICosmo* cosmo = ts->cosmo;

	if (z > xclk->zmax)
		return 0.0;

	gdouble result, error;
	// static gsl_integration_workspace* w = NULL; //gsl_integration_workspace_alloc (INT_PARTITION);
	gsl_integration_workspace* w = gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);

	gsl_function F;
	int_src_int_params int_ts;

	const gdouble chiz = nc_distance_comoving (xclkg->dist, cosmo, z);
	int_ts.chiz = chiz;
	int_ts.dist = xclkg->dist;
	int_ts.cosmo = cosmo;
	int_ts.dn_dz = xclkg->dn_dz;

	// w = gsl_integration_workspace_alloc (NCM_INTEGRAL_PARTITION);

	F.function = &_nc_xcor_limber_kernel_weak_lensing_src_int_integrand;
	F.params = &int_ts;

	gsl_integration_qag (&F, GSL_MAX (z, xclk->zmin), xclk->zmax, 0., NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, w, &result, &error);

	// printf ("integration result = %g with integration error = %g \n", result, error);

	gsl_integration_workspace_free (w);
	return result;
}

static void
_nc_xcor_limber_kernel_weak_lensing_prepare (NcXcorLimberKernel* xclk, NcHICosmo* cosmo)
{
	NcXcorLimberKernelWeakLensing* xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (xclk);
	// NcmModel* model = NCM_MODEL (xclk);

	xclk->cons_factor = 1.5 * nc_hicosmo_Omega_m0(cosmo);

	xclk->zmin = 0.0;

	nc_distance_prepare_if_needed (xclkg->dist, cosmo);


	{
		// guint i;
		//
		// if (xclkg->g_func == NULL)
		// {
		// 	NcmVector* gzv;
		//
		// 	xclkg->g_func = ncm_spline_cubic_notaknot_new ();
		// 	ncm_spline_set_len (xclkg->g_func, NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_G_FUNC_LEN);
		//
		// 	gzv = ncm_spline_get_xv (xclkg->g_func);
		//
		// 	for (i = 0; i < NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_G_FUNC_LEN; i++)
		// 	{
		// 		ncm_vector_set (gzv, i, xclk->zmin + (xclk->zmax - xclk->zmin) / (NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_G_FUNC_LEN - 1) * i);
		// 	}
		//
		// 	ncm_vector_free (gzv);
		// 	/*ncm_spline_acc (xclkg->g_func, TRUE);*/
		// }

		gsl_function F;

		src_int_params ts;

		ts.xclkg = xclkg;
		ts.cosmo = cosmo;

		F.function = &_nc_xcor_limber_kernel_weak_lensing_src_int;
		F.params = &ts;

		ncm_spline_set_func (xclkg->src_int, NCM_SPLINE_FUNCTION_SPLINE, &F, xclk->zmin, xclk->zmax, 1000, 1e-5);

		// {
		// 	g_func_params ts;
		// 	NcmVector* xv = ncm_spline_get_xv (xclkg->g_func);
		// 	NcmVector* yv = ncm_spline_get_yv (xclkg->g_func);
		//
		// 	ts.xclkg = xclkg;
		// 	ts.cosmo = cosmo;
		//
		// 	for (i = 0; i < NC_XCOR_LIMBER_KERNEL_WEAK_LENSING_G_FUNC_LEN; i++)
		// 	{
		// 		ncm_vector_set (yv, i, _nc_xcor_limber_kernel_weak_lensing_g_func (ncm_vector_get (xv, i), &ts));
		// 	}
		//
		// 	ncm_vector_free (xv);
		// 	ncm_vector_free (yv);
		// }

		ncm_spline_prepare (xclkg->src_int);
	}
}

// /**
//  * nc_xcor_limber_kernel_weak_lensing_set_dndz:
//  * @xclkg: a #NcXcorLimberKernelWeakLensing
//  * @z: (element-type double): a #GArray
//  * @dn_dz_array: (element-type double): a #GArray
//  *
//  * FIXME
//  *
//  * Returns: FIXME
//  *
// */
// void
// nc_xcor_limber_kernel_weak_lensing_set_dndz (NcXcorLimberKernelWeakLensing* xclkg, GArray* z, GArray* dn_dz_array)
// {
// 	ncm_spline_set_array (xclkg->dn_dz, z, dn_dz_array, TRUE);
//
// 	// NcXcorLimberKernel* xclk = NC_XCOR_LIMBER_KERNEL (xclkg);
// 	// xclk->zmin = g_array_index (z, gdouble, 0);
// 	// xclk->zmax = g_array_index (z, gdouble, z->len - 1);
// }


static gdouble
_nc_xcor_limber_kernel_weak_lensing_eval (NcXcorLimberKernel* xclk, NcHICosmo* cosmo, gdouble z, const NcXcorKinetic* xck, gint l)
{
	NcXcorLimberKernelWeakLensing* xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (xclk);

	NCM_UNUSED (l);
	NCM_UNUSED (cosmo);

	// printf("%g\n", xck->xi_z);
	// printf("%g\n", xck->E_z);
	// printf("%g\n", ncm_spline_eval( xclkg->src_int, z));

	return xck->xi_z / xck->E_z * (1. + z) * ncm_spline_eval( xclkg->src_int, z);
}


static void
_nc_xcor_limber_kernel_weak_lensing_add_noise (NcXcorLimberKernel* xclk, NcmVector* vp1, NcmVector* vp2, guint lmin)
{
	NcXcorLimberKernelWeakLensing* xclkg = NC_XCOR_LIMBER_KERNEL_WEAK_LENSING (xclk);

	// ncm_vector_add_constant (vp1, NOISE_BIAS); /* take noise_bias into account */
	ncm_vector_memcpy (vp2, vp1);
	ncm_vector_add_constant (vp2, xclkg->noise);
}

guint _nc_xcor_limber_kernel_weak_lensing_obs_len (NcXcorLimberKernel* xclk)
{
	NCM_UNUSED (xclk);
	return 2;
}

guint _nc_xcor_limber_kernel_weak_lensing_obs_params_len (NcXcorLimberKernel* xclk)
{
	NCM_UNUSED (xclk);
	return 1;
}
