/***************************************************************************
 *            ncm_obj_array.h
 *
 *  Wed October 16 11:04:36 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
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

G_BEGIN_DECLS

#define NCM_TYPE_OBJ_ARRAY (ncm_obj_array_get_type ())
#define NCM_TYPE_OBJ_DICT_STR (ncm_obj_dict_str_get_type ())
#define NCM_TYPE_OBJ_DICT_INT (ncm_obj_dict_int_get_type ())

#define NCM_OBJ_ARRAY(obj) ((NcmObjArray *) (obj))
#define NCM_OBJ_DICT_STR(obj) ((NcmObjDictStr *) (obj))
#define NCM_OBJ_DICT_INT(obj) ((NcmObjDictInt *) (obj))

GType ncm_obj_array_get_type (void) G_GNUC_CONST;
GType ncm_obj_dict_str_get_type (void) G_GNUC_CONST;
GType ncm_obj_dict_int_get_type (void) G_GNUC_CONST;

/*
 * NcmObjArray
 */

typedef struct _GPtrArray NcmObjArray;

NcmObjArray *ncm_obj_array_new (void);
NcmObjArray *ncm_obj_array_sized_new (guint n);
NcmObjArray *ncm_obj_array_ref (NcmObjArray *oa);
void ncm_obj_array_unref (NcmObjArray *oa);
void ncm_obj_array_clear (NcmObjArray **oa);

void ncm_obj_array_add (NcmObjArray *oa, GObject *obj);
void ncm_obj_array_set (NcmObjArray *oa, guint i, GObject *obj);
GObject *ncm_obj_array_get (NcmObjArray *oa, guint i);
GObject *ncm_obj_array_peek (NcmObjArray *oa, guint i);

guint ncm_obj_array_len (NcmObjArray *oa);

/*
 * NcmObjDictStr
 */

typedef struct _GHashTable NcmObjDictStr;

NcmObjDictStr *ncm_obj_dict_str_new (void);
NcmObjDictStr *ncm_obj_dict_str_ref (NcmObjDictStr *od);
void ncm_obj_dict_str_unref (NcmObjDictStr *od);
void ncm_obj_dict_str_clear (NcmObjDictStr **od);

void ncm_obj_dict_str_add (NcmObjDictStr *od, const gchar *key, GObject *obj);
void ncm_obj_dict_str_set (NcmObjDictStr *od, const gchar *key, GObject *obj);
GObject *ncm_obj_dict_str_get (NcmObjDictStr *od, const gchar *key);
GObject *ncm_obj_dict_str_peek (NcmObjDictStr *od, const gchar *key);

guint ncm_obj_dict_str_len (NcmObjDictStr *od);
GStrv ncm_obj_dict_str_keys (NcmObjDictStr *od);

/*
 * NcmObjDictInt
 */

typedef struct _GHashTable NcmObjDictInt;

NcmObjDictInt *ncm_obj_dict_int_new (void);
NcmObjDictInt *ncm_obj_dict_int_ref (NcmObjDictInt *od);
void ncm_obj_dict_int_unref (NcmObjDictInt *od);
void ncm_obj_dict_int_clear (NcmObjDictInt **od);

void ncm_obj_dict_int_add (NcmObjDictInt *od, gint key, GObject *obj);
void ncm_obj_dict_int_set (NcmObjDictInt *od, gint key, GObject *obj);
GObject *ncm_obj_dict_int_get (NcmObjDictInt *od, gint key);
GObject *ncm_obj_dict_int_peek (NcmObjDictInt *od, gint key);

guint ncm_obj_dict_int_len (NcmObjDictInt *od);
GArray *ncm_obj_dict_int_keys (NcmObjDictInt *od);

G_END_DECLS

#endif /* NCM_OBJ_ARRAY_H */

