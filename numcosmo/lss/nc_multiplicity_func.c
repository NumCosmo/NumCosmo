/***************************************************************************
 *            nc_multiplicity_func.c
 *
 *  Mon Jun 28 15:09:13 2010
 *  Copyright  2010  Mariana Penna Lima
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
 * SECTION:nc_multiplicity_func
 * @title: NcMultiplicityFunc
 * @short_description: Dark matter halo multiplicity function.
 *
 * FIXME
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "lss/nc_multiplicity_func.h"
#include "math/ncm_serialize.h"
#include "math/ncm_cfg.h"

G_DEFINE_ABSTRACT_TYPE (NcMultiplicityFunc, nc_multiplicity_func, G_TYPE_OBJECT);

/**
 * nc_multiplicity_func_new_from_name:
 * @multiplicity_name: string which specifies the multiplicity function type.
 *
 * This function returns a new #NcMultiplicityFunc whose type is defined by @multiplicity_name.
 *
 * Returns: A new #NcMultiplicityFunc.
 */
NcMultiplicityFunc *
nc_multiplicity_func_new_from_name (gchar *multiplicity_name)
{
  GObject *obj = ncm_serialize_global_from_string (multiplicity_name);
  GType multiplicity_type = G_OBJECT_TYPE (obj);
  if (!g_type_is_a (multiplicity_type, NC_TYPE_MULTIPLICITY_FUNC))
	g_error ("nc_multiplicity_func_new_from_name: NcMultiplicityFunc %s do not descend from %s.", multiplicity_name, g_type_name (NC_TYPE_MULTIPLICITY_FUNC));
  return NC_MULTIPLICITY_FUNC (obj);
}

/**
 * nc_multiplicity_func_eval:
 * @mulf: a #NcMultiplicityFunc.
 * @cosmo: a #NcHICosmo.
 * @sigma: FIXME
 * @z: redshift.
 *
 * FIXME
 *
 * Returns: FIXME
*/
gdouble
nc_multiplicity_func_eval (NcMultiplicityFunc *mulf, NcHICosmo *cosmo, gdouble sigma, gdouble z)
{
  return NC_MULTIPLICITY_FUNC_GET_CLASS (mulf)->eval (mulf, cosmo, sigma, z);
}

/**
 * nc_multiplicity_func_free:
 * @mulf: a #NcMultiplicityFunc.
 *
 * Atomically decrements the reference count of @mulf by one. If the reference count drops to 0,
 * all memory allocated by @mulf is released.
 *
*/
void
nc_multiplicity_func_free (NcMultiplicityFunc *mulf)
{
  g_object_unref (mulf);
}

/**
 * nc_multiplicity_func_clear:
 * @mulf: a #NcMultiplicityFunc.
 *
 * Atomically decrements the reference count of @mulf by one. If the reference count drops to 0,
 * all memory allocated by @mulf is released. Set pointer to NULL.
 *
*/
void
nc_multiplicity_func_clear (NcMultiplicityFunc **mulf)
{
  g_clear_object (mulf);
}

static void
nc_multiplicity_func_init (NcMultiplicityFunc *mulf)
{
  NCM_UNUSED (mulf);
}

static void
_nc_multiplicity_func_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_multiplicity_func_parent_class)->finalize (object);
}

static void
nc_multiplicity_func_class_init (NcMultiplicityFuncClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  //GObjectClass* parent_class = G_OBJECT_CLASS (klass);

  object_class->finalize = _nc_multiplicity_func_finalize;
}

