/***************************************************************************
 *            nc_halo_bias.c
 *
 *  Tue June 28 15:41:57 2011
 *  Copyright  2011  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
 * SECTION:nc_halo_bias
 * @title: NcHaloBias
 * @short_description: Abstract class for halo bias function type.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_bias.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_ABSTRACT_TYPE (NcHaloBias, nc_halo_bias, G_TYPE_OBJECT);


static void
nc_halo_bias_init (NcHaloBias *mbias)
{
  mbias->mfp = NULL;
  mbias->bias = NULL;
}

static void
_nc_halo_bias_dispose (GObject *object)
{
  NcHaloBias *mbias = nc_halo_bias (object);
  
  nc_halo_mass_function_clear (&mbias->mfp);
  nc_halo_bias_clear (&mbias->bias);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_bias_parent_class)->dispose (object);
}

static void
_nc_halo_bias_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_bias_parent_class)->finalize (object);
}

static void
_nc_halo_bias_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHaloBias *mbias = NC_HALO_BIAS (object);
  g_return_if_fail (NC_IS_HALO_BIAS (object));

  switch (prop_id)
  {
	case PROP_MASS_FUNCTION:
	  mbias->mfp = g_value_dup_object (value);
	  break;
	default:
	  G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
	  break;
  }
}

static void
_nc_halo_bias_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHaloBias *mbias = NC_HALO_BIAS (object);
  g_return_if_fail (NC_IS_HALO_BIAS (object));

  switch (prop_id)
  {
	case PROP_MASS_FUNCTION:
	  g_value_set_object (value, mbias->mfp);
	  break;
	default:
	  G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
	  break;
  }
}

static void
nc_halo_bias_class_init (NcHaloBiasClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  //GObjectClass* parent_class = G_OBJECT_CLASS (klass);

  object_class->dispose = _nc_halo_bias_func_dispose;
  object_class->finalize = _nc_halo_bias_func_finalize;
  object_class->set_property = _nc_halo_bias_func_set_property;
  object_class->get_property = _nc_halo_bias_func_get_property;

  /**
   * NcHaloBias:mass-function:
   *
   * This property keeps the mass function object.
   */
  g_object_class_install_property (object_class,
                                   PROP_MASS_FUNCTION,
                                   g_param_spec_object ("mass-function",
                                                        NULL,
                                                        "Mass Function.",
                                                        NC_TYPE_HALO_MASS_FUNCTION,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


}


/**
 * nc_halo_bias_new_from_name:
 * @bias_name: string which specifies the multiplicity function type.
 *
 * This function returns a new #NcMultiplicityFunc whose type is defined by @multiplicity_name.
 *
 * Returns: A new #NcHaloBIas.
 */
NcHaloBias *
nc_halo_bias_new_from_name (gchar *bias_name)
{
  GObject *obj = ncm_serialize_global_from_string (bias_name);
  GType bias_type = G_OBJECT_TYPE (obj);

  if (!g_type_is_a (bias_type, NC_TYPE_HALO_BIAS))
	  g_error ("nc_halo_bias_new_from_name: NcHaloBias %s do not descend from %s.", bias_name, g_type_name (NC_TYPE_HALO_BIAS));
  return NC_HALO_BIAS (obj);
}

/**
 * nc_halo_bias_eval:
 * @bias: a #NcHaloBias.
 * @sigma: FIXME
 * @z: redshift.
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_halo_bias_eval (NcHaloBias *bias, gdouble sigma, gdouble z)
{
  return NC_HALO_BIAS_GET_CLASS (bias)->eval (bias, sigma, z);
}

/**
 * nc_halo_bias_free:
 * @bias: a #NcHaloBias.
 *
 * Atomically decrements the reference count of @bias by one. If the reference count drops to 0,
 * all memory allocated by @bias is released.
 *
 */
void
nc_halo_bias_free (NcHaloBias *bias)
{
  g_object_unref (bias);
}

/**
 * nc_halo_bias_clear:
 * @bias: a #NcHaloBias.
 *
 * Atomically decrements the reference count of @bias by one. If the reference count drops to 0,
 * all memory allocated by @bias is released. Set pointer to NULL.
 *
 */
void
nc_halo_bias_clear (NcHaloBias **bias)
{
  g_clear_object (bias);
}
