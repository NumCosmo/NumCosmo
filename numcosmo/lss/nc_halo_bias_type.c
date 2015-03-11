/***************************************************************************
 *            nc_halo_bias_type.c
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
 * SECTION:nc_halo_bias_type
 * @title: NcHaloBiasType
 * @short_description: Abstract class for halo bias function type.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_halo_bias_type.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"
#include "math/ncm_util.h"

G_DEFINE_ABSTRACT_TYPE (NcHaloBiasType, nc_halo_bias_type, G_TYPE_OBJECT);

/**
 * nc_halo_bias_type_new_from_name:
 * @bias_name: string which specifies the multiplicity function type.
 *
 * This function returns a new #NcMultiplicityFunc whose type is defined by @multiplicity_name.
 *
 * Returns: A new #NcHaloBiasType.
 */
NcHaloBiasType *
nc_halo_bias_type_new_from_name (gchar *bias_name)
{
  GObject *obj = ncm_serialize_global_from_string (bias_name);
  GType bias_type = G_OBJECT_TYPE (obj);

  if (!g_type_is_a (bias_type, NC_TYPE_HALO_BIAS_TYPE))
	g_error ("nc_halo_bias_type_new_from_name: NcHaloBiasType %s do not descend from %s.", bias_name, g_type_name (NC_TYPE_HALO_BIAS_TYPE));
  return NC_HALO_BIAS_TYPE (obj);
}

/**
 * nc_halo_bias_type_eval:
 * @biasf: a #NcHaloBiasType.
 * @sigma: FIXME
 * @z: redshift.
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_halo_bias_type_eval (NcHaloBiasType *biasf, gdouble sigma, gdouble z)
{
  return NC_HALO_BIAS_TYPE_GET_CLASS (biasf)->eval (biasf, sigma, z);
}

/**
 * nc_halo_bias_type_free:
 * @biasf: a #NcHaloBiasType.
 *
 * Atomically decrements the reference count of @biasf by one. If the reference count drops to 0,
 * all memory allocated by @biasf is released.
 *
 */
void
nc_halo_bias_type_free (NcHaloBiasType *biasf)
{
  g_object_unref (biasf);
}

/**
 * nc_halo_bias_type_clear:
 * @biasf: a #NcHaloBiasType.
 *
 * Atomically decrements the reference count of @biasf by one. If the reference count drops to 0,
 * all memory allocated by @biasf is released. Set pointer to NULL.
 *
 */
void
nc_halo_bias_type_clear (NcHaloBiasType **biasf)
{
  g_clear_object (biasf);
}

static void
nc_halo_bias_type_init (NcHaloBiasType *nc_halo_bias_type)
{
  NCM_UNUSED (nc_halo_bias_type);
}

static void
_nc_halo_bias_type_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_halo_bias_type_parent_class)->finalize (object);
}

static void
nc_halo_bias_type_class_init (NcHaloBiasTypeClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  //GObjectClass* parent_class = G_OBJECT_CLASS (klass);

  object_class->finalize = _nc_halo_bias_type_finalize;
}

