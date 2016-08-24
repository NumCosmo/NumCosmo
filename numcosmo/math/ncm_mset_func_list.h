/***************************************************************************
 *            ncm_mset_func_list.h
 *
 *  Mon August 08 17:29:24 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_mset_func_list.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_MSET_FUNC_LIST_H_
#define _NCM_MSET_FUNC_LIST_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset_func.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_FUNC_LIST             (ncm_mset_func_list_get_type ())
#define NCM_MSET_FUNC_LIST(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MSET_FUNC_LIST, NcmMSetFuncList))
#define NCM_MSET_FUNC_LIST_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MSET_FUNC_LIST, NcmMSetFuncListClass))
#define NCM_IS_MSET_FUNC_LIST(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MSET_FUNC_LIST))
#define NCM_IS_MSET_FUNC_LIST_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MSET_FUNC_LIST))
#define NCM_MSET_FUNC_LIST_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MSET_FUNC_LIST, NcmMSetFuncListClass))

typedef struct _NcmMSetFuncListClass NcmMSetFuncListClass;
typedef struct _NcmMSetFuncList NcmMSetFuncList;

typedef void (*NcmMSetFuncListN) (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *res);

struct _NcmMSetFuncListClass
{
  /*< private >*/
  NcmMSetFuncClass parent_class;
  GArray *func_array;
  GHashTable *ns_hash;
};

typedef struct _NcmMSetFuncListStruct
{
  gchar *name;
  gchar *symbol;
  gchar *ns;
  gchar *desc;
  GType obj_type;
  NcmMSetFuncListN func;
  guint nvar;
  guint dim;
  guint pos;
} NcmMSetFuncListStruct;

struct _NcmMSetFuncList
{
  /*< private >*/
  NcmMSetFunc parent_instance;
  GType obj_type;
  NcmMSetFuncListN func;
  GObject *obj;
};

GType ncm_mset_func_list_get_type (void) G_GNUC_CONST;

void ncm_mset_func_list_register (const gchar *name, const gchar *symbol, const gchar *ns, const gchar *desc, GType obj_type, NcmMSetFuncListN func, guint nvar, guint dim);
GArray *ncm_mset_func_list_select (const gchar *ns, gint nvar, gint dim);

NcmMSetFuncList *ncm_mset_func_list_new (const gchar *full_name, gpointer obj);
NcmMSetFuncList *ncm_mset_func_list_new_ns_name (const gchar *ns, const gchar *name, gpointer obj);

gboolean ncm_mset_func_list_has_ns_name (const gchar *ns, const gchar *name);
gboolean ncm_mset_func_list_has_full_name (const gchar *full_name);

G_END_DECLS

#endif /* _NCM_MSET_FUNC_LIST_H_ */

