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
#define NCM_TYPE_VAR_DICT (ncm_var_dict_get_type ())

#define NCM_OBJ_ARRAY(obj) ((NcmObjArray *) (obj))
#define NCM_OBJ_DICT_STR(obj) ((NcmObjDictStr *) (obj))
#define NCM_OBJ_DICT_INT(obj) ((NcmObjDictInt *) (obj))
#define NCM_VAR_DICT(obj) ((NcmVarDict *) (obj))

GType ncm_obj_array_get_type (void) G_GNUC_CONST;
GType ncm_obj_dict_str_get_type (void) G_GNUC_CONST;
GType ncm_obj_dict_int_get_type (void) G_GNUC_CONST;
GType ncm_var_dict_get_type (void) G_GNUC_CONST;

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
NcmObjDictStr *ncm_obj_dict_str_ref (NcmObjDictStr *ods);
void ncm_obj_dict_str_unref (NcmObjDictStr *ods);
void ncm_obj_dict_str_clear (NcmObjDictStr **ods);

void ncm_obj_dict_str_add (NcmObjDictStr *ods, const gchar *key, GObject *obj);
void ncm_obj_dict_str_set (NcmObjDictStr *ods, const gchar *key, GObject *obj);
GObject *ncm_obj_dict_str_get (NcmObjDictStr *ods, const gchar *key);
GObject *ncm_obj_dict_str_peek (NcmObjDictStr *ods, const gchar *key);

guint ncm_obj_dict_str_len (NcmObjDictStr *ods);
GStrv ncm_obj_dict_str_keys (NcmObjDictStr *ods);

/*
 * NcmObjDictInt
 */

typedef struct _GHashTable NcmObjDictInt;

NcmObjDictInt *ncm_obj_dict_int_new (void);
NcmObjDictInt *ncm_obj_dict_int_ref (NcmObjDictInt *odi);
void ncm_obj_dict_int_unref (NcmObjDictInt *odi);
void ncm_obj_dict_int_clear (NcmObjDictInt **odi);

void ncm_obj_dict_int_add (NcmObjDictInt *odi, gint key, GObject *obj);
void ncm_obj_dict_int_set (NcmObjDictInt *odi, gint key, GObject *obj);
GObject *ncm_obj_dict_int_get (NcmObjDictInt *odi, gint key);
GObject *ncm_obj_dict_int_peek (NcmObjDictInt *odi, gint key);

guint ncm_obj_dict_int_len (NcmObjDictInt *odi);
GArray *ncm_obj_dict_int_keys (NcmObjDictInt *odi);

/*
 * NcmVarDict
 */
typedef struct _GHashTable NcmVarDict;

NcmVarDict *ncm_var_dict_new (void);
NcmVarDict *ncm_var_dict_ref (NcmVarDict *vd);
void ncm_var_dict_unref (NcmVarDict *vd);
void ncm_var_dict_clear (NcmVarDict **vd);

void ncm_var_dict_set_string (NcmVarDict *vd, const gchar *key, const gchar *value);
void ncm_var_dict_set_int (NcmVarDict *vd, const gchar *key, gint value);
void ncm_var_dict_set_double (NcmVarDict *vd, const gchar *key, gdouble value);
void ncm_var_dict_set_boolean (NcmVarDict *vd, const gchar *key, gboolean value);
void ncm_var_dict_set_int_array (NcmVarDict *vd, const gchar *key, GArray *value);
void ncm_var_dict_set_double_array (NcmVarDict *vd, const gchar *key, GArray *value);
void ncm_var_dict_set_boolean_array (NcmVarDict *vd, const gchar *key, GArray *value);
void ncm_var_dict_set_variant (NcmVarDict *vd, const gchar *key, GVariant *value);

gboolean ncm_var_dict_has_key (NcmVarDict *vd, const gchar *key);
gboolean ncm_var_dict_get_string (NcmVarDict *vd, const gchar *key, gchar **value);
gboolean ncm_var_dict_get_int (NcmVarDict *vd, const gchar *key, gint *value);
gboolean ncm_var_dict_get_double (NcmVarDict *vd, const gchar *key, gdouble *value);
gboolean ncm_var_dict_get_boolean (NcmVarDict *vd, const gchar *key, gboolean *value);
gboolean ncm_var_dict_get_int_array (NcmVarDict *vd, const gchar *key, GArray **value);
gboolean ncm_var_dict_get_double_array (NcmVarDict *vd, const gchar *key, GArray **value);
gboolean ncm_var_dict_get_boolean_array (NcmVarDict *vd, const gchar *key, GArray **value);
gboolean ncm_var_dict_get_variant (NcmVarDict *vd, const gchar *key, GVariant **value);

guint ncm_var_dict_len (NcmVarDict *vd);
GStrv ncm_var_dict_keys (NcmVarDict *vd);


G_END_DECLS

#endif /* NCM_OBJ_ARRAY_H */

