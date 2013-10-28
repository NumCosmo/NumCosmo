/***************************************************************************
 *            ncm_obj_array.c
 *
 *  Wed October 16 11:04:01 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_obj_array.c
 *
 * Copyright (C) 2013 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:ncm_obj_array
 * @title: GObjects array.
 * @short_description: GObjects array with serialization support. 
 *
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_obj_array.h"

G_DEFINE_BOXED_TYPE (NcmObjArray, ncm_obj_array, ncm_obj_array_ref, ncm_obj_array_unref);

/**
 * ncm_obj_array_new:
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmObjArray *
ncm_obj_array_new ()
{
  GPtrArray *oa = g_ptr_array_new ();
  g_ptr_array_set_free_func (oa, g_object_unref);
  return (NcmObjArray *) oa;
}

/**
 * ncm_obj_array_new_from_variant:
 * @ser: a #NcmSerialize.
 * @var: a #GVariant containing an array of objects.
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmObjArray *
ncm_obj_array_new_from_variant (NcmSerialize *ser, GVariant *var)
{
  g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE (NCM_OBJ_ARRAY_TYPE)));
  {
    guint i, n = g_variant_n_children (var);
    NcmObjArray *oa = ncm_obj_array_sized_new (n);
    for (i = 0; i < n; i++)
    {
      GVariant *cvar = g_variant_get_child_value (var, i);
      GObject *cobj = ncm_serialize_from_variant (ser, cvar);
      ncm_obj_array_add (oa, cobj);
      g_object_unref (cobj);
      g_variant_unref (cvar);
    }
    return oa;
  }
}

/**
 * ncm_obj_array_sized_new:
 * @n: initial allocation size.
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmObjArray *
ncm_obj_array_sized_new (guint n)
{
  GPtrArray *oa = g_ptr_array_sized_new (n);
  g_ptr_array_set_free_func (oa, g_object_unref);
  return (NcmObjArray *) oa;
}

/**
 * ncm_obj_array_ref:
 * @oa: a #NcmObjArray.
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
NcmObjArray *
ncm_obj_array_ref (NcmObjArray *oa)
{
  return (NcmObjArray *)g_ptr_array_ref ((GPtrArray *)oa);
}

/**
 * ncm_obj_array_unref:
 * @oa: a #NcmObjArray.
 * 
 * FIXME
 *
 */
void 
ncm_obj_array_unref (NcmObjArray *oa)
{
  g_ptr_array_unref ((GPtrArray *)oa);
}

/**
 * ncm_obj_array_clear:
 * @oa: a pointer to a #NcmObjArray.
 * 
 * FIXME
 *
 */
void 
ncm_obj_array_clear (NcmObjArray **oa)
{
  g_clear_pointer (oa, &ncm_obj_array_unref);
}

GVariant *_ncm_serialize_to_variant (NcmSerialize *ser, GObject *obj);

GVariant *
_ncm_obj_array_ser (NcmObjArray *oa, NcmSerialize *ser)
{
  GVariantBuilder *builder;
  guint i;

  builder = g_variant_builder_new (G_VARIANT_TYPE (NCM_OBJ_ARRAY_TYPE));

  for (i = 0; i < oa->len; i++)
  {
    GVariant *cvar = _ncm_serialize_to_variant (ser, ncm_obj_array_peek (oa, i));
    g_variant_builder_add_value (builder, cvar);
  }

  return g_variant_ref_sink (g_variant_builder_end (builder));
}


/**
 * ncm_obj_array_ser:
 * @oa: a #NcmObjArray.
 * @ser: a #NcmSerialize.
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
GVariant *
ncm_obj_array_ser (NcmObjArray *oa, NcmSerialize *ser)
{
  GVariantBuilder *builder;
  guint i;

  builder = g_variant_builder_new (G_VARIANT_TYPE (NCM_OBJ_ARRAY_TYPE));

  for (i = 0; i < oa->len; i++)
  {
    GVariant *cvar = ncm_serialize_to_variant (ser, ncm_obj_array_peek (oa, i));
    g_variant_builder_add_value (builder, cvar);
  }

  return g_variant_ref_sink (g_variant_builder_end (builder));
}

/**
 * ncm_obj_array_add:
 * @oa: a #NcmObjArray.
 * @obj: a #GObject.
 * 
 * FIXME
 *
 */
void 
ncm_obj_array_add (NcmObjArray *oa, GObject *obj)
{
  g_ptr_array_add ((GPtrArray *)oa, g_object_ref (obj));
}

/**
 * ncm_obj_array_set:
 * @oa: a #NcmObjArray.
 * @i: object index.
 * @obj: a #GObject.
 * 
 * FIXME
 *
 */
void 
ncm_obj_array_set (NcmObjArray *oa, guint i, GObject *obj)
{
  g_assert_cmpuint (i, <, oa->len);
  g_object_unref (g_ptr_array_index ((GPtrArray *)oa, i));
  g_ptr_array_index ((GPtrArray *)oa, i) = g_object_ref (obj);
}

/**
 * ncm_obj_array_get:
 * @oa: a #NcmObjArray.
 * @i: object index.
 * 
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
GObject *
ncm_obj_array_get (NcmObjArray *oa, guint i)
{
  return g_object_ref (ncm_obj_array_peek (oa, i));
}

/**
 * ncm_obj_array_peek:
 * @oa: a #NcmObjArray.
 * @i: object index.
 * 
 * FIXME
 *
 * Returns: (transfer none): FIXME
 */
GObject *
ncm_obj_array_peek (NcmObjArray *oa, guint i)
{
  g_assert_cmpuint (i, <, oa->len);
  return g_ptr_array_index ((GPtrArray *)oa, i);
}

