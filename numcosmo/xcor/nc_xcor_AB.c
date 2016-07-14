/***************************************************************************
 *            nc_xcor_AB.h
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
 * SECTION:nc_xcor_AB
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
#include "math/ncm_cfg.h"
#include "math/ncm_serialize.h"
#include "xcor/nc_xcor_AB.h"

enum
{
	PROP_0,
	PROP_A,
	PROP_B,
	PROP_ELL_TH_CUT_OFF,
	PROP_ELL_LIK_MIN,
	PROP_ELL_LIK_MAX,
	PROP_MIXING,
	PROP_CL_TH,
	PROP_CL_OBS,
	PROP_SIZE,
};

G_DEFINE_TYPE (NcXcorAB, nc_xcor_AB, G_TYPE_OBJECT);

static void 
nc_xcor_AB_init (NcXcorAB* xcab)
{
	xcab->a = 99;
	xcab->b = 99;

	xcab->ell_th_cut_off = 0;
	xcab->ell_lik_min    = 0;
	xcab->ell_lik_max    = 0;
	xcab->nell_lik       = 0;

	xcab->mixing = NULL;
	xcab->cl_th  = NULL; //column 0 : C_l^th, 1 : C_l^th+N_l, 2 : mixed C_l
	xcab->cl_obs = NULL;
}

static void
_nc_xcor_AB_set_property (GObject* object, guint prop_id, const GValue* value, GParamSpec* pspec)
{
	NcXcorAB* xcab = NC_XCOR_AB (object);
	g_return_if_fail (NC_IS_XCOR_AB (object));

	switch (prop_id)
  {
    case PROP_A:
      xcab->a = g_value_get_uint (value);
      break;
    case PROP_B:
      xcab->b = g_value_get_uint (value);
      break;
    case PROP_ELL_TH_CUT_OFF:
      xcab->ell_th_cut_off = g_value_get_uint (value);
      break;
    case PROP_ELL_LIK_MIN:
      xcab->ell_lik_min = g_value_get_uint (value);
      xcab->nell_lik = xcab->ell_lik_max - xcab->ell_lik_min + 1;
      break;
    case PROP_ELL_LIK_MAX:
      xcab->ell_lik_max = g_value_get_uint (value);
      xcab->nell_lik = xcab->ell_lik_max - xcab->ell_lik_min + 1;
      break;
    case PROP_MIXING:
      xcab->mixing = g_value_dup_object (value);
      break;
    case PROP_CL_TH:
      xcab->cl_th = g_value_dup_object (value);
      break;
    case PROP_CL_OBS:
      xcab->cl_obs = g_value_dup_object (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_xcor_AB_get_property (GObject* object, guint prop_id, GValue* value, GParamSpec* pspec)
{
  NcXcorAB* xcab = NC_XCOR_AB (object);
  g_return_if_fail (NC_IS_XCOR_AB (object));

  switch (prop_id)
  {
    case PROP_A:
      g_value_set_uint (value, xcab->a);
      break;
    case PROP_B:
      g_value_set_uint (value, xcab->b);
      break;
    case PROP_ELL_TH_CUT_OFF:
      g_value_set_uint (value, xcab->ell_th_cut_off);
      break;
    case PROP_ELL_LIK_MIN:
      g_value_set_uint (value, xcab->ell_lik_min);
      break;
    case PROP_ELL_LIK_MAX:
      g_value_set_uint (value, xcab->ell_lik_max);
      break;
    case PROP_MIXING:
      g_value_set_object (value, xcab->mixing);
      break;
    case PROP_CL_TH:
      g_value_set_object (value, xcab->cl_th);
      break;
    case PROP_CL_OBS:
      g_value_set_object (value, xcab->cl_obs);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void 
_nc_xcor_AB_dispose (GObject* object)
{
	NcXcorAB* xcab = NC_XCOR_AB (object);

	ncm_matrix_clear (&xcab->mixing);
	ncm_matrix_clear (&xcab->cl_th);
	ncm_vector_clear (&xcab->cl_obs);

	/* Chain up : end */
	G_OBJECT_CLASS (nc_xcor_AB_parent_class)->dispose (object);
}

static void 
_nc_xcor_AB_finalize (GObject* object)
{
  
	/* Chain up : end */
	G_OBJECT_CLASS (nc_xcor_AB_parent_class)->finalize (object);
}

static void
nc_xcor_AB_class_init (NcXcorABClass* klass)
{
	GObjectClass* object_class = G_OBJECT_CLASS (klass);
	//GObjectClass* parent_class = G_OBJECT_CLASS (klass);

	object_class->dispose = _nc_xcor_AB_dispose;
	object_class->finalize = _nc_xcor_AB_finalize;
	object_class->set_property = _nc_xcor_AB_set_property;
	object_class->get_property = _nc_xcor_AB_get_property;

	g_object_class_install_property (object_class,
	                                 PROP_A,
	                                 g_param_spec_uint ("a",
	                                                    NULL,
	                                                    "a",
	                                                    0, 99, 99,
	                                                    G_PARAM_CONSTRUCT_ONLY | G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_B,
	                                 g_param_spec_uint ("b",
	                                                    NULL,
	                                                    "b",
	                                                    0, 99, 99,
	                                                    G_PARAM_CONSTRUCT_ONLY | G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_ELL_TH_CUT_OFF,
	                                 g_param_spec_uint ("ell-th-cut-off",
	                                                    NULL,
	                                                    "ell_th_cut_off",
	                                                    0, 10000, 0,
	                                                    G_PARAM_CONSTRUCT | G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_ELL_LIK_MIN,
	                                 g_param_spec_uint ("ell-lik-min",
	                                                    NULL,
	                                                    "ell_lik_min",
	                                                    0, 10000, 0,
	                                                    G_PARAM_CONSTRUCT | G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

	g_object_class_install_property (object_class,
	                                 PROP_ELL_LIK_MAX,
	                                 g_param_spec_uint ("ell-lik-max",
	                                                    NULL,
	                                                    "ell_lik_max",
	                                                    0, 10000, 0,
	                                                    G_PARAM_CONSTRUCT | G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_MIXING,
	                                 g_param_spec_object ("mixing",
	                                                      NULL,
	                                                      "mixing",
	                                                      NCM_TYPE_MATRIX,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_CL_TH,
	                                 g_param_spec_object ("cl-th",
	                                                      NULL,
	                                                      "cl_th",
	                                                      NCM_TYPE_MATRIX,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
	g_object_class_install_property (object_class,
	                                 PROP_CL_OBS,
	                                 g_param_spec_object ("cl-obs",
	                                                      NULL,
	                                                      "cl_obs",
	                                                      NCM_TYPE_VECTOR,
	                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * nc_xcor_AB_new:
 * @a: a #guint
 * @b: a #guint
 * @ell_th_cut_off: a #guint
 * @ell_lik_min: a #guint
 * @ell_lik_max: a #guint
 * @clobs_filename: a #gchar
 * @mixing_filename: a #gchar
 * @mixing_filelength: a #guint
 *
 * FIXME
 *
 * Returns: FIXME
 *
 */
NcXcorAB *
nc_xcor_AB_new (guint a, guint b, guint ell_th_cut_off, guint ell_lik_min, guint ell_lik_max, const gchar *clobs_filename, const gchar *mixing_filename, const guint mixing_filelength)
{
	NcmMatrix *mixing = ncm_matrix_new (ell_th_cut_off + 1, ell_th_cut_off + 1);
	NcmMatrix *cl_th  = ncm_matrix_new (ell_th_cut_off + 1, 2);
	NcmVector* cl_obs = ncm_vector_new (ell_lik_max + 1);

  ncm_matrix_set_zero (cl_th);

  if (b < a)
		g_error ("nc_xcor_AB_new: b must be greater or equal to a");

	if (clobs_filename != NULL)
	{
    FILE *f = fopen (clobs_filename, "r");
    gsl_vector_fscanf (f, ncm_vector_gsl (cl_obs));
    fclose (f);
  }
	else
    cl_obs = NULL;

	if (mixing_filename != NULL)
	{
		NcmMatrix *mixing_full = ncm_matrix_new (mixing_filelength, mixing_filelength);
		FILE* g = fopen (mixing_filename, "r");
    
		gsl_matrix_fscanf (g, ncm_matrix_gsl (mixing_full));

    fclose (g);

		ncm_matrix_memcpy (mixing, ncm_matrix_get_submatrix (mixing_full, 0, 0, ell_th_cut_off + 1, ell_th_cut_off + 1));
    ncm_matrix_free (mixing_full);
	}
	else
    mixing = NULL;

  {
    NcXcorAB* xcab = g_object_new (NC_TYPE_XCOR_AB,
                                   "a", a,
                                   "b", b,
                                   "ell-th-cut-off", ell_th_cut_off,
                                   "ell-lik-min", ell_lik_min,
                                   "ell-lik-max", ell_lik_max,
                                   "mixing", mixing,
                                   "cl-th", cl_th,
                                   "cl-obs", cl_obs,
                                   NULL);

    xcab->nell_lik = ell_lik_max - ell_lik_min + 1;

    return xcab;
  }
}

/**
 * nc_xcor_AB_ref:
 * @xcab: a #NcXcorAB
 * *
 * Returns: (transfer full): @xcab
 */
NcXcorAB* nc_xcor_AB_ref (NcXcorAB* xcab)
{
	return g_object_ref (xcab);
}

/**
 * nc_xcor_AB_free:
 * @xcab: a #NcXcorAB
 *
 * FIXME
 *
 */
void nc_xcor_AB_free (NcXcorAB* xcab)
{
	g_object_unref (xcab);
}

/**
 * nc_xcor_AB_clear:
 * @xcab: a #NcXcorAB
 *
 * FIXME
 *
 */
void nc_xcor_AB_clear (NcXcorAB** xcab)
{
	g_clear_object (xcab);
}
