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
 * 
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "xcor/nc_xcor.h"
#include "math/integral.h"
#include "math/memory_pool.h"
#include "math/ncm_cfg.h"
#include "math/ncm_serialize.h"

#include <cuba.h>
#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>

enum
{
	PROP_0,
	PROP_DISTANCE,
	PROP_MATTER_POWER_SPECTRUM,
};

G_DEFINE_TYPE (NcXcor, nc_xcor, G_TYPE_OBJECT);

static void 
nc_xcor_init (NcXcor *xc)
{
	xc->ps   = NULL;
	xc->dist = NULL;
	xc->RH   = 0.0;
}

static void
_nc_xcor_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
	NcXcor* xc = NC_XCOR (object);
	g_return_if_fail (NC_IS_XCOR (object));

	switch (prop_id)
	{
    case PROP_DISTANCE:
      xc->dist = g_value_dup_object (value);
      break;
    case PROP_MATTER_POWER_SPECTRUM:
      xc->ps   = g_value_dup_object (value);
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
    case PROP_MATTER_POWER_SPECTRUM:
      g_value_set_object (value, xc->ps);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void 
_nc_xcor_dispose (GObject* object)
{
	NcXcor *xc = NC_XCOR (object);

	nc_distance_clear (&xc->dist);
	ncm_powspec_clear (&xc->ps);

	/* Chain up : end */
	G_OBJECT_CLASS (nc_xcor_parent_class)->dispose (object);
}

static void 
_nc_xcor_finalize (GObject* object)
{
  
	/* Chain up : end */
	G_OBJECT_CLASS (nc_xcor_parent_class)->finalize (object);
}

static void
nc_xcor_class_init (NcXcorClass *klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	//GObjectClass* parent_class = G_OBJECT_CLASS (klass);

	object_class->set_property = &_nc_xcor_set_property;
	object_class->get_property = &_nc_xcor_get_property;
	object_class->dispose      = &_nc_xcor_dispose;
	object_class->finalize     = &_nc_xcor_finalize;

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
   * NcXcor:power-spec:
   *
   * This property keeps the matter power spectrum object.
   */
	g_object_class_install_property (object_class,
	                                 PROP_MATTER_POWER_SPECTRUM,
	                                 g_param_spec_object ("power-spec",
	                                                      NULL,
	                                                      "Matter power spectrum.",
	                                                      NCM_TYPE_POWSPEC,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_xcor_new:
 * @dist: a #NcDistance
 * @ps: a #NcmPowspec
 *
 * FIXME
 *
 * Returns: FIXME
 *
*/
NcXcor *
nc_xcor_new (NcDistance* dist, NcmPowspec* ps)
{
	return g_object_new (NC_TYPE_XCOR, 
                       "distance", dist, 
                       "power-spec", ps, 
                       NULL);
}

/**
 * nc_xcor_ref:
 * @xc: a #NcXcor
 *
 * Returns: (transfer full): @xc
 */
NcXcor *
nc_xcor_ref (NcXcor* xc)
{
	return g_object_ref (xc);
}

/**
 * nc_xcor_free:
 * @xc: a #NcXcor
 *
 * FIXME
 *
 */
void 
nc_xcor_free (NcXcor *xc)
{
	g_object_unref (xc);
}

/**
 * nc_xcor_clear:
 * @xc: a #NcXcor
 *
 * FIXME
 *
 */
void 
nc_xcor_clear (NcXcor **xc)
{
	g_clear_object (xc);
}

/**
 * nc_xcor_prepare:
 * @xc: a #NcXcor
 * @cosmo: a #NcHICosmo
 * 
 * FIXME
 * 
 */
void 
nc_xcor_prepare (NcXcor *xc, NcHICosmo *cosmo)
{
	nc_distance_prepare_if_needed (xc->dist, cosmo);
	ncm_powspec_prepare_if_needed (xc->ps, NCM_MODEL (cosmo));

	xc->RH = nc_hicosmo_RH_Mpc (cosmo);
}

#ifdef HAVE_SUNDIALS
typedef struct _xcor_limber_cvode
{
	gboolean isauto;
	guint lmin, lmax;
	guint nell;
	NcXcor *xc;
	NcXcorLimberKernel *xclk1;
	NcXcorLimberKernel *xclk2;
	NcHICosmo *cosmo;
} xcor_limber_cvode;

static gint
_xcor_limber_cvode_int (realtype z, N_Vector y, N_Vector ydot, gpointer params)
{
	xcor_limber_cvode* xclc = (xcor_limber_cvode*)params;
	const gdouble xi_z      = nc_distance_comoving (xclc->xc->dist, xclc->cosmo, z); // in units of Hubble radius
	const gdouble xi_z_phys = xi_z * xclc->xc->RH; // in Mpc-1
	const gdouble E_z       = nc_hicosmo_E (xclc->cosmo, z);
	const gdouble k1z       = nc_xcor_limber_kernel_eval (xclc->xclk1, xclc->cosmo, z, 0);
	const guint nell        = xclc->nell;
	NcmVector* Pk           = ncm_vector_new (nell);
	NcmVector* k            = ncm_vector_new (nell);

  gdouble geoW1W2;
  guint i, l;

	if (xclc->isauto)
	{
		geoW1W2 = E_z * gsl_pow_2 (k1z / xi_z);
	}
	else
	{
		const gdouble k2z = nc_xcor_limber_kernel_eval (xclc->xclk2, xclc->cosmo, z, 0);
		geoW1W2 = E_z * k1z * k2z / (xi_z * xi_z);
	}

	for (i = 0; i < nell; i++)
	{
		l = i + xclc->lmin;
		ncm_vector_set (k, i, (l + 0.5) / xi_z_phys);
	}

	ncm_powspec_eval_vec (xclc->xc->ps, NCM_MODEL (xclc->cosmo), z, k, Pk);

	for (i = 0; i < nell; i++)
	{
		NV_Ith_S (ydot, i) = ncm_vector_get (Pk, i) * geoW1W2;
	}

	ncm_vector_free (Pk);
	ncm_vector_free (k);

	return 0;
}

static void 
_nc_xcor_limber_cvode (NcXcor* xc, NcXcorLimberKernel* xclk1, NcXcorLimberKernel* xclk2, NcHICosmo* cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector* vp)
{
	const guint nell   = lmax - lmin + 1;
	const gdouble init = NCM_DEFAULT_PRECISION;
	const gdouble zmid = exp ((log1p (zmin) + log1p (zmax)) / 2.0) - 1.0;
	N_Vector yv        = N_VNew_Serial (nell);
	gpointer cvode     = CVodeCreate (CV_ADAMS, CV_FUNCTIONAL);
	gpointer cvodefunc = &_xcor_limber_cvode_int; //isauto ? &_xcor_limber_cvode_auto_int : &_xcor_limber_cvode_cross_int;
	gdouble z;
	gint flag;
	guint i;

  for (i = 0; i < nell; i++)
	{
		NV_Ith_S (yv, i) = init;
	}
  
	xcor_limber_cvode xclc = { isauto, lmin, lmax, nell, xc, xclk1, xclk2, cosmo }; //, cons_factor };

	/* First integrate from zmid to zmin*/
	CVodeInit (cvode, cvodefunc, zmid, yv);

	CVodeSStolerances (cvode, NCM_DEFAULT_PRECISION, 0.0);
	CVodeSetMaxNumSteps (cvode, 5000);
	CVodeSetUserData (cvode, &xclc);

	flag = CVode (cvode, zmin, yv, &z, CV_NORMAL);
	NCM_CVODE_CHECK (&flag, "CVode", 1, );

	for (i = 0; i < nell; i++)
	{
		NV_Ith_S (yv, i) = -NV_Ith_S (yv, i) + init; //integration done backwards, hence the minus sign
	}

	/* Then integrate from zmid to zmax*/
	CVodeInit (cvode, cvodefunc, zmid, yv);

	CVodeSStolerances (cvode, NCM_DEFAULT_PRECISION, 0.0);
	CVodeSetMaxNumSteps (cvode, 5000);
	CVodeSetUserData (cvode, &xclc);

	flag = CVode (cvode, zmax, yv, &z, CV_NORMAL);
	NCM_CVODE_CHECK (&flag, "CVode", 1, );

	for (i = 0; i < nell; i++)
	{
		ncm_vector_set (vp, i, NV_Ith_S (yv, i));
	}

	CVodeFree (&cvode);
	N_VDestroy (yv);
}

#endif /* HAVE_SUNDIALS_2_5_0 */

typedef struct _xcor_limber_gsl
{
	NcHICosmo* cosmo;
	NcDistance* dist;
	NcmPowspec* ps;

	NcXcorLimberKernel *xclk1;
	NcXcorLimberKernel *xclk2;
	guint l;

	gdouble RH;

} xcor_limber_gsl;

static gdouble 
_xcor_limber_gsl_cross_int (gdouble z, gpointer ptr)
{
	xcor_limber_gsl *xclki   = (xcor_limber_gsl*)ptr;
	const gdouble xi_z       = nc_distance_comoving (xclki->dist, xclki->cosmo, z); // in units of Hubble radius
	const gdouble xi_z_phys  = xi_z * xclki->RH; // in Mpc-1
	const gdouble E_z        = nc_hicosmo_E (xclki->cosmo, z);
	const gdouble k          = (xclki->l + 0.5) / (xi_z_phys); // in Mpc-1
	const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xclki->ps), NCM_MODEL (xclki->cosmo), k, z);

	const gdouble k1z = nc_xcor_limber_kernel_eval (xclki->xclk1, xclki->cosmo, z, xclki->l);
	const gdouble k2z = nc_xcor_limber_kernel_eval (xclki->xclk2, xclki->cosmo, z, xclki->l);

	return E_z * k1z * k2z * power_spec / (xi_z * xi_z);
}

static gdouble 
_xcor_limber_gsl_auto_int (gdouble z, gpointer ptr)
{
	xcor_limber_gsl *xclki   = (xcor_limber_gsl*)ptr;
	const gdouble xi_z       = nc_distance_comoving (xclki->dist, xclki->cosmo, z); // in units of Hubble radius
	const gdouble xi_z_phys  = xi_z * xclki->RH; // in Mpc-1
	const gdouble E_z        = nc_hicosmo_E (xclki->cosmo, z);
	const gdouble k          = (xclki->l + 0.5) / (xi_z_phys); // in Mpc-1
	const gdouble power_spec = ncm_powspec_eval (NCM_POWSPEC (xclki->ps), NCM_MODEL (xclki->cosmo), k, z);
	const gdouble k1z        = nc_xcor_limber_kernel_eval (xclki->xclk1, xclki->cosmo, z, xclki->l);

	return E_z * gsl_pow_2 (k1z / xi_z) * power_spec;
}

static void 
_nc_xcor_limber_gsl (NcXcor* xc, NcXcorLimberKernel* xclk1, NcXcorLimberKernel* xclk2, NcHICosmo* cosmo, guint lmin, guint lmax, gdouble zmin, gdouble zmax, gboolean isauto, NcmVector* vp)
{
	xcor_limber_gsl xclki;
	gdouble r, err;
	gsl_function F;
	guint i;

	xclki.xclk1 = xclk1;
	xclki.xclk2 = xclk2;
	xclki.cosmo = cosmo;
	xclki.dist  = xc->dist;
	xclki.ps    = xc->ps;
	xclki.RH    = xc->RH;

	if (isauto)
	{
		F.function = &_xcor_limber_gsl_auto_int;
	}
	else
	{
		F.function = &_xcor_limber_gsl_cross_int;
	}

  F.params = &xclki;

	gsl_integration_workspace **w = ncm_integral_get_workspace ();

	for (i = 0; i < lmax - lmin + 1; i++)
	{
		xclki.l = lmin + i;
		gsl_integration_qag (&F, zmin, zmax, 0.0, NCM_DEFAULT_PRECISION, NCM_INTEGRAL_PARTITION, 6, *w, &r, &err);
		ncm_vector_set (vp, i, r);
	}

	ncm_memory_pool_return (w);
}

/**
 * nc_xcor_limber:
 * @xc: a #NcXcor
 * @xclk1: a #NcXcorLimberKernel
 * @xclk2: a #NcXcorLimberKernel
 * @cosmo: a #NcHICosmo
 * @lmin: a #guint
 * @lmax: a #guint
 * @vp: a #NcmVector
 * @meth: a #NcXcorLimberMethod
 *
 * FIXME
 *
 */
void 
nc_xcor_limber (NcXcor* xc, NcXcorLimberKernel* xclk1, NcXcorLimberKernel* xclk2, NcHICosmo* cosmo, guint lmin, guint lmax, NcmVector* vp, NcXcorLimberMethod meth)
{
	const guint nell          = ncm_vector_len (vp);
	const gboolean isauto     = (xclk2 == NULL);
	const gdouble cons_factor = ((isauto) ? gsl_pow_2 (xclk1->cons_factor) : xclk1->cons_factor * xclk2->cons_factor) / gsl_pow_3 (xc->RH);
	gdouble zmin, zmax;

  if (nell != lmax - lmin + 1)
    g_error ("nc_xcor_limber: vector size does not match multipole limits");

	if (isauto)
	{
		zmin = xclk1->zmin;
		zmax = xclk1->zmax;
	}
	else
	{
		zmin = GSL_MAX (xclk1->zmin, xclk2->zmin);
		zmax = GSL_MIN (xclk1->zmax, xclk2->zmax);
	}

	if (zmin < zmax)
	{
		switch (meth)
		{
      case NC_XCOR_LIMBER_METHOD_CVODE:
        _nc_xcor_limber_cvode (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
      case NC_XCOR_LIMBER_METHOD_GSL:
        _nc_xcor_limber_gsl (xc, xclk1, xclk2, cosmo, lmin, lmax, zmin, zmax, isauto, vp);
        break;
    }
		ncm_vector_scale (vp, cons_factor);
	}
	else
	{
		ncm_vector_set_zero (vp);
	}
}
