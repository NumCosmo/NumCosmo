/***************************************************************************
 *            ncm_obj_array.h
 *
 *  Wed October 16 11:04:36 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_obj_array.h
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

#ifndef NCM_OBJ_ARRAY_H
#define NCM_OBJ_ARRAY_H

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_serialize.h>

G_BEGIN_DECLS

#define NCM_TYPE_OBJ_ARRAY (ncm_obj_array_get_type ())

typedef struct _GPtrArray NcmObjArray;

GType ncm_obj_array_get_type (void) G_GNUC_CONST;

NcmObjArray *ncm_obj_array_new (void);
NcmObjArray *ncm_obj_array_new_from_variant (NcmSerialize *ser, GVariant *var);
NcmObjArray *ncm_obj_array_sized_new (guint n);
NcmObjArray *ncm_obj_array_ref (NcmObjArray *oa);
void ncm_obj_array_unref (NcmObjArray *oa);
void ncm_obj_array_clear (NcmObjArray **oa);
GVariant *ncm_obj_array_ser (NcmObjArray *oa, NcmSerialize *ser);

void ncm_obj_array_add (NcmObjArray *oa, GObject *obj);
void ncm_obj_array_set (NcmObjArray *oa, guint i, GObject *obj);
GObject *ncm_obj_array_get (NcmObjArray *oa, guint i);
GObject *ncm_obj_array_peek (NcmObjArray *oa, guint i);

#define NCM_OBJ_ARRAY_TYPE "a"NCM_SERIALIZE_OBJECT_TYPE

G_END_DECLS

#endif /* NCM_OBJ_ARRAY_H */
