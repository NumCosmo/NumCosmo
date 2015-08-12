/***************************************************************************
 *            nc_xcor_limber_gal.c
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "xcor/nc_xcor_limber_gal.h"
#include "math/ncm_cfg.h"

#include <gsl/gsl_randist.h>

G_DEFINE_TYPE (NcXcorLimberGal, nc_xcor_limber_gal, NC_TYPE_XCOR_LIMBER);

#define VECTOR (NCM_MODEL (xclg)->params)
#define BIAS (ncm_vector_get (VECTOR, NC_XCOR_LIMBER_GAL_BIAS))

enum
{
	PROP_0,
	PROP_DN_DZ,
	PROP_Z_MIN,
	PROP_Z_MAX,
	PROP_NBARM1,
	PROP_SIZE,
};

static gdouble _nc_xcor_limber_gal_eval_kernel (NcXcorLimber* xcl, NcHICosmo* cosmo, gdouble z, gint l)
{
	NcXcorLimberGal* xclg = NC_XCOR_LIMBER_GAL (xcl);

	NCM_UNUSED (l);
	NCM_UNUSED (cosmo);

	if (z < xclg->z_min || z > xclg->z_max)
	{
		return 0.0;
	}
	else
	{
		const gdouble dN_dz_z = ncm_spline_eval (xclg->dN_dz, z);
		return BIAS * dN_dz_z;
	}
}

static void _nc_xcor_limber_gal_prepare (NcXcorLimber* xcl, NcHICosmo* cosmo)
{
	xcl->cons_factor = 1.0;
	NCM_UNUSED (cosmo);
}

static gdouble _nc_xcor_limber_gal_noise_spec (NcXcorLimber* xcl, guint l)
{
	NcXcorLimberGal* xclg = NC_XCOR_LIMBER_GAL (xcl);

	return xclg->nbarm1;
}

/**
 * nc_xcor_limber_gal_set_dNdz:
 * @xclg: a #NcXcorLimberGal
 * @z: (element-type double): a #GArray
 * @dN_dz_array: (element-type double): a #GArray
 *
 * FIXME
 *
 * Returns: FIXME
 *
*/
void nc_xcor_limber_gal_set_dNdz (NcXcorLimberGal* xclg, GArray* z, GArray* dN_dz_array)
{
	ncm_spline_set_array (xclg->dN_dz, z, dN_dz_array, TRUE);
}

guint _nc_xcor_limber_gal_obs_len (NcXcorLimber* xcl)
{
	NCM_UNUSED (xcl);
	return 1;
}
guint _nc_xcor_limber_gal_obs_params_len (NcXcorLimber* xcl)
{
	NCM_UNUSED (xcl);
	return 0;
}

static void
_nc_xcor_limber_gal_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
	NcXcorLimberGal* xclg = NC_XCOR_LIMBER_GAL (object);
	g_return_if_fail (NC_IS_XCOR_LIMBER_GAL (object));

	switch (prop_id)
	{
	case PROP_DN_DZ:
		xclg->dN_dz = g_value_get_object (value);
		break;
	case PROP_Z_MAX:
		xclg->z_max = g_value_get_double (value);
		break;
	case PROP_Z_MIN:
		xclg->z_min = g_value_get_double (value);
		break;
	case PROP_NBARM1:
		xclg->nbarm1 = g_value_get_double (value);
		break;

	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
_nc_xcor_limber_gal_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
	NcXcorLimberGal* xclg = NC_XCOR_LIMBER_GAL (object);
	g_return_if_fail (NC_IS_XCOR_LIMBER_GAL (object));

	switch (prop_id)
	{
	case PROP_DN_DZ:
		g_value_set_object (value, xclg->dN_dz);
		break;
	case PROP_Z_MIN:
		g_value_set_double (value, xclg->z_min);
		break;
	case PROP_Z_MAX:
		g_value_set_double (value, xclg->z_max);
		break;
	case PROP_NBARM1:
		g_value_set_double (value, xclg->nbarm1);
		break;

	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
nc_xcor_limber_gal_init (NcXcorLimberGal* xclg)
{
	xclg->z_min = 0.0;
	xclg->z_max = 20.0;
	xclg->nbarm1 = 0.0;
	xclg->dN_dz = ncm_spline_cubic_notaknot_new ();
	(xclg->parent_instance).cons_factor = 0.0;
}

static void
_nc_xcor_limber_gal_finalize (GObject* object)
{
	NcXcorLimberGal* xclg = NC_XCOR_LIMBER_GAL (object);
	ncm_spline_clear (&xclg->dN_dz);

	/* Chain up : end */
	G_OBJECT_CLASS (nc_xcor_limber_gal_parent_class)
	->finalize (object);
}

static void
nc_xcor_limber_gal_class_init (NcXcorLimberGalClass* klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcXcorLimberClass* parent_class = NC_XCOR_LIMBER_CLASS (klass);
	NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

	parent_class->eval_kernel = &_nc_xcor_limber_gal_eval_kernel;
	parent_class->prepare = &_nc_xcor_limber_gal_prepare;
	parent_class->noise_spec = &_nc_xcor_limber_gal_noise_spec;

	parent_class->obs_len = &_nc_xcor_limber_gal_obs_len;
	parent_class->obs_params_len = &_nc_xcor_limber_gal_obs_params_len;

	parent_class->impl = NC_XCOR_LIMBER_IMPL_ALL;

	object_class->finalize = &_nc_xcor_limber_gal_finalize;

	model_class->set_property = &_nc_xcor_limber_gal_set_property;
	model_class->get_property = &_nc_xcor_limber_gal_get_property;

	ncm_model_class_add_params (model_class, 1, 0, PROP_SIZE);
	ncm_model_class_set_name_nick (model_class, "Xcor quasar distribution", "Xcor-gal");

	/**
  * NcXcorLimberGal:dN_dz:
  *
  * FIXME Set correct values (limits)
  */
	g_object_class_install_property (object_class,
	                                 PROP_DN_DZ,
	                                 g_param_spec_object ("dNdz",
	                                                      NULL,
	                                                      "Quasar distribution",
	                                                      NCM_TYPE_SPLINE,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB)); //G_PARAM_CONSTRUCT_ONLY |

	/**
   * NcClusterMassLnnormal:z_max:
   *
   * FIXME Set correct values (limits)
   */
	g_object_class_install_property (object_class,
	                                 PROP_Z_MAX,
	                                 g_param_spec_double ("zmax",
	                                                      NULL,
	                                                      "Maximum redshift",
	                                                      0.0, 20.0, 10.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	/**
    * NcClusterMassLnnormal:z_min:
    *
    * FIXME Set correct values (limits)
    */
	g_object_class_install_property (object_class,
	                                 PROP_Z_MIN,
	                                 g_param_spec_double ("zmin",
	                                                      NULL,
	                                                      "Minimum redshift",
	                                                      0.0, 20.0, 0.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	/**
	* NcClusterMassLnnormal:nbarm1:
	*
	* FIXME Set correct values (limits)
	*/
	g_object_class_install_property (object_class,
	                                 PROP_NBARM1,
	                                 g_param_spec_double ("nbarm1",
	                                                      NULL,
	                                                      "One over nbar (galaxy angular density)",
	                                                      0.0, 20.0, 0.0,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	/*
    * Distribution's bias: bias.
    * FIXME Set correct values (limits)
    */
	ncm_model_class_set_sparam (model_class, NC_XCOR_LIMBER_GAL_BIAS, "bias", "bias",
	                            0.0, 10.0, 1.0e-2,
	                            NC_XCOR_LIMBER_GAL_DEFAULT_PARAMS_ABSTOL, NC_XCOR_LIMBER_GAL_DEFAULT_BIAS,
	                            NCM_PARAM_TYPE_FIXED);

	/* Check for errors in parameters initialization */
	ncm_model_class_check_params_info (model_class);
}
