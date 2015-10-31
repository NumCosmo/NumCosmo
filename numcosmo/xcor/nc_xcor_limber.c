/***************************************************************************
 *            nc_xcor_limber.c
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

#include "math/integral.h"
#include "math/memory_pool.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "lss/nc_window_tophat.h"
#include "lss/nc_matter_var.h"
#include "xcor/nc_xcor_limber.h"

G_DEFINE_ABSTRACT_TYPE (NcXcorLimber, nc_xcor_limber, NCM_TYPE_MODEL);

/**
 * nc_xcor_limber_new_from_name:
 * @xcor_name: string which specifies the type of the observable.
 *
 * This function returns a new #NcXcorLimber whose type is defined by
 *@xcor_name.
 *
 * Returns: A new #NcXcorLimber.
 */
NcXcorLimber* nc_xcor_limber_new_from_name (gchar* xcor_name)
{
	GObject* obj = ncm_serialize_global_from_string (xcor_name);
	GType xcor_type = G_OBJECT_TYPE (obj);
	if (!g_type_is_a (xcor_type, NC_TYPE_XCOR_LIMBER))
		g_error ("nc_xcor_limber_new_from_name: NcXcorLimber %s do not "
		         "descend from %s.",
		         xcor_name, g_type_name (NC_TYPE_XCOR_LIMBER));
	return NC_XCOR_LIMBER (obj);
}

/**
 * nc_xcor_limber_ref:
 * @xcl: a #NcXcorLimber.
 *
 * FIXME
 *
 * Returns: (transfer full): @xcl.
 */
NcXcorLimber* nc_xcor_limber_ref (NcXcorLimber* xcl)
{
	return g_object_ref (xcl);
}

/**
 * nc_xcor_limber_free:
 * @xcl: a #NcXcorLimber.
 *
 * FIXME
 *
 */
void nc_xcor_limber_free (NcXcorLimber* xcl)
{
	g_object_unref (xcl);
}

/**
 * nc_xcor_limber_clear:
 * @xcl: a #NcXcorLimber.
 *
 * FIXME
 *
 */
void nc_xcor_limber_clear (NcXcorLimber** xcl)
{
	g_clear_object (xcl);
}

/**
 * nc_xcor_limber_impl:
 * @xcl: a #NcXcorLimber.
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcXcorLimberImpl nc_xcor_limber_impl (NcXcorLimber* xcl)
{
	return NC_XCOR_LIMBER_GET_CLASS (xcl)->impl;
}

/**
 * nc_xcor_limber_obs_len:
 * @xcl: a #NcXcorLimber.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint nc_xcor_limber_obs_len (NcXcorLimber* xcl)
{
	return NC_XCOR_LIMBER_GET_CLASS (xcl)->obs_len (xcl);
}

/**
 * nc_xcor_limber_obs_params_len:
 * @xcl: a #NcXcorLimber.
 *
 * FIXME
 *
 * Returns: FIXME
 */
guint nc_xcor_limber_obs_params_len (NcXcorLimber* xcl)
{
	return NC_XCOR_LIMBER_GET_CLASS (xcl)->obs_params_len (xcl);
}

/**
 * nc_xcor_limber_eval_kernel:
 * @xcl: a #NcXcorLimber
 * @cosmo: a #NcHICosmo
 * @z: a #gdouble
 * @l: a #gint
 *
 * FIXME
 *
 * Returns: FIXME
 *
*/
gdouble nc_xcor_limber_eval_kernel (NcXcorLimber* xcl, NcHICosmo* cosmo, gdouble z, gint l)
{
	return NC_XCOR_LIMBER_GET_CLASS (xcl)->eval_kernel (xcl, cosmo, z, l);
}

/**
 * nc_xcor_limber_noise_spec:
 * @xcl: a #NcXcorLimber
 * @l: a #guint
 *
 * FIXME
 *
 * Returns: FIXME
 *
*/
gdouble nc_xcor_limber_noise_spec (NcXcorLimber* xcl, guint l)
{
	return NC_XCOR_LIMBER_GET_CLASS (xcl)->noise_spec (xcl, l);
}


/**
 * nc_xcor_limber_prepare:
 * @xcl: a NcXcorLimber
 * @cosmo: a NcHICosmo
 *
 * FIXME
 *
 * Returns: FIXME
*/
void nc_xcor_limber_prepare (NcXcorLimber* xcl, NcHICosmo* cosmo)
{
	if (ncm_model_ctrl_update (xcl->cosmo_ctrl, NCM_MODEL (cosmo)))
	{
		return NC_XCOR_LIMBER_GET_CLASS (xcl)->prepare (xcl, cosmo);
	}
}

static void _nc_xcor_limber_log_all_models_go (GType model_type, guint n)
{
	guint nc, i, j;
	GType* models = g_type_children (model_type, &nc);
	for (i = 0; i < nc; i++)
	{
		guint ncc;
		GType* modelsc = g_type_children (models[i], &ncc);

		g_message ("#  ");
		for (j = 0; j < n; j++)
			g_message (" ");
		g_message ("%s\n", g_type_name (models[i]));
		if (ncc)
			_nc_xcor_limber_log_all_models_go (models[i], n + 2);

		g_free (modelsc);
	}
	g_free (models);
}

/**
 * nc_xcor_limber_log_all_models:
 *
 * FIXME
 *
 */
void nc_xcor_limber_log_all_models (void)
{
	g_message ("# Registred NcXcorLimber:%s are:\n",
	           g_type_name (NC_TYPE_XCOR_LIMBER));
	_nc_xcor_limber_log_all_models_go (NC_TYPE_XCOR_LIMBER, 0);
}

static void
nc_xcor_limber_init (NcXcorLimber* xcl)
{
	xcl->cosmo_ctrl = ncm_model_ctrl_new(NULL);
	xcl->cons_factor = 0.0;
}

static void nc_xcor_limber_finalize (GObject* object)
{

	/* Chain up : end */
	G_OBJECT_CLASS (nc_xcor_limber_parent_class)->finalize (object);
}

NCM_MSET_MODEL_REGISTER_ID (nc_xcor_limber, NC_TYPE_XCOR_LIMBER);

static void nc_xcor_limber_class_init (NcXcorLimberClass* klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

	object_class->finalize = nc_xcor_limber_finalize;

	ncm_model_class_set_name_nick (model_class, "Xcor limber", "Xcor-limber");
	ncm_model_class_add_params (model_class, 0, 0, 1);

	ncm_mset_model_register_id (model_class, "NcXcorLimber", "Cluster mass observable relation models (this at "
                                                           "line 297 of nc_cor_limber.c, to be modified "
	                                                         "accordingly).",
	                            NULL);
	ncm_model_class_check_params_info (NCM_MODEL_CLASS (klass));
}
