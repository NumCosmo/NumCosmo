/***************************************************************************
 *            nc_xcor_limber_lensing.c
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

#include "math/ncm_cfg.h"
#include "xcor/nc_xcor_limber_lensing.h"
#include <gsl/gsl_randist.h>

G_DEFINE_TYPE (NcXcorLimberLensing, nc_xcor_limber_lensing, NC_TYPE_XCOR_LIMBER);

#define VECTOR (NCM_MODEL (xcll)->params)

enum
{
	PROP_0,
	PROP_DIST,
	PROP_NL,
	PROP_SIZE,
};

static gdouble
_nc_xcor_limber_lensing_eval_kernel (NcXcorLimber* xcl, NcHICosmo* cosmo, gdouble z, gint l)
{
	NcXcorLimberLensing* xcll = NC_XCOR_LIMBER_LENSING (xcl);

	NCM_UNUSED (l);

	gdouble xi_z = nc_distance_comoving (xcll->dist, cosmo, z);
	gdouble xi_lss = nc_distance_comoving_lss (xcll->dist, cosmo);
	gdouble E_z = nc_hicosmo_E (cosmo, z);

	return ((1.0 + z) * xi_z * (xi_lss - xi_z)) / (E_z * xi_lss);
}

static void _nc_xcor_limber_lensing_prepare (NcXcorLimber* xcl, NcHICosmo* cosmo)
{
	NcXcorLimberLensing* xcll = NC_XCOR_LIMBER_LENSING (xcl);

	nc_distance_prepare_if_needed (xcll->dist, cosmo);

	xcl->cons_factor = (3.0 * nc_hicosmo_Omega_m (cosmo)) / 2.0;
}

static gdouble _nc_xcor_limber_lensing_noise_spec (NcXcorLimber* xcl, guint l)
{
	NcXcorLimberLensing* xcll = NC_XCOR_LIMBER_LENSING (xcl);

	if (xcll->Nl == NULL) g_error ("nc_xcor_limber_lensing_noise_spec : noise spectrum empty");

	if (l > xcll->Nlmax) g_error ("nc_xcor_limber_lensing_noise_spec : too high multipole");

	return ncm_vector_get (xcll->Nl, l);
}

guint _nc_xcor_limber_lensing_obs_len (NcXcorLimber* xcl)
{
	NCM_UNUSED (xcl);
	return 1;
}
guint _nc_xcor_limber_lensing_obs_params_len (NcXcorLimber* xcl)
{
	NCM_UNUSED (xcl);
	return 0;
}

static void
_nc_xcor_limber_lensing_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
	NcXcorLimberLensing* xcll = NC_XCOR_LIMBER_LENSING (object);
	g_return_if_fail (NC_IS_XCOR_LIMBER_LENSING (object));

	switch (prop_id)
	{
	case PROP_DIST:
		xcll->dist = g_value_get_object (value);
		break;
	case PROP_NL:
		xcll->Nl = g_value_get_object (value);
		xcll->Nlmax = ncm_vector_len (xcll->Nl);
		break;
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
_nc_xcor_limber_lensing_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
	NcXcorLimberLensing* xcll = NC_XCOR_LIMBER_LENSING (object);
	g_return_if_fail (NC_IS_XCOR_LIMBER_LENSING (object));

	switch (prop_id)
	{
	case PROP_DIST:
		g_value_set_object (value, xcll->dist);
		break;
	case PROP_NL:
		g_value_set_object (value, xcll->Nl);
		break;
	default:
		G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
		break;
	}
}

static void
nc_xcor_limber_lensing_init (NcXcorLimberLensing* xcll)
{
	xcll->dist = NULL;
	xcll->Nl = NULL;
	xcll->Nlmax = 0;
	(xcll->parent_instance).cons_factor = 0.0;
}

static void
_nc_xcor_limber_lensing_finalize (GObject* object)
{
	/* Chain up : end */
	G_OBJECT_CLASS (nc_xcor_limber_lensing_parent_class)
	->finalize (object);
}

static void
nc_xcor_limber_lensing_class_init (NcXcorLimberLensingClass* klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	NcXcorLimberClass* parent_class = NC_XCOR_LIMBER_CLASS (klass);
	NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

	parent_class->eval_kernel = &_nc_xcor_limber_lensing_eval_kernel;
	parent_class->prepare = &_nc_xcor_limber_lensing_prepare;
	parent_class->noise_spec = &_nc_xcor_limber_lensing_noise_spec;

	parent_class->obs_len = &_nc_xcor_limber_lensing_obs_len;
	parent_class->obs_params_len = &_nc_xcor_limber_lensing_obs_params_len;

	parent_class->impl = NC_XCOR_LIMBER_IMPL_ALL;

	object_class->finalize = &_nc_xcor_limber_lensing_finalize;

	model_class->set_property = &_nc_xcor_limber_lensing_set_property;
	model_class->get_property = &_nc_xcor_limber_lensing_get_property;

	ncm_model_class_add_params (model_class, 0, 0, PROP_SIZE);
	ncm_model_class_set_name_nick (model_class, "Xcor lensing distribution", "Xcor-lensing");

	/**
     * NcXcorLimberLensing:dist:
     *
     * FIXME Set correct values (limits)
     */
	g_object_class_install_property (object_class,
	                                 PROP_DIST,
	                                 g_param_spec_object ("dist",
	                                                      NULL,
	                                                      "Distance object",
	                                                      NC_TYPE_DISTANCE,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	/**
     * NcXcorLimberLensing:Nl:
     *
     * FIXME Set correct values (limits)
     */
	g_object_class_install_property (object_class,
	                                 PROP_NL,
	                                 g_param_spec_object ("Nl",
	                                                      NULL,
	                                                      "Noise spectrum",
	                                                      NCM_TYPE_VECTOR,
	                                                      G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	/* Check for errors in parameters initialization */
	ncm_model_class_check_params_info (model_class);
}
