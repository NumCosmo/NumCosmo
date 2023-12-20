/***************************************************************************
 *            ncm_mset_func_list.h
 *
 *  Mon August 08 17:29:24 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_func_list.h
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#define NCM_TYPE_MSET_FUNC_LIST (ncm_mset_func_list_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmMSetFuncList, ncm_mset_func_list, NCM, MSET_FUNC_LIST, NcmMSetFunc)

typedef void (*NcmMSetFuncListN) (NcmMSetFuncList *flist, NcmMSet *mset, const gdouble *x, gdouble *res);

struct _NcmMSetFuncListClass
{
  /*< private >*/
  NcmMSetFuncClass parent_class;
  GArray *func_array;
  GHashTable *ns_hash;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[16];
};

/**
 * NcmMSetFuncListStruct:
 * @name: FIXME
 * @symbol: FIXME
 * @ns: FIXME
 * @desc: FIXME
 *
 * FIXME
 *
 */
typedef struct _NcmMSetFuncListStruct
{
  gchar *name;
  gchar *symbol;
  gchar *ns;
  gchar *desc;
  /*< private >*/
  GType obj_type;
  NcmMSetFuncListN func;
  guint nvar;
  guint dim;
  guint pos;
} NcmMSetFuncListStruct;

void ncm_mset_func_list_register (const gchar *name, const gchar *symbol, const gchar *ns, const gchar *desc, GType obj_type, NcmMSetFuncListN func, guint nvar, guint dim);
GArray *ncm_mset_func_list_select (const gchar *ns, gint nvar, gint dim);

NcmMSetFuncList *ncm_mset_func_list_new (const gchar *full_name, GObject *obj);
NcmMSetFuncList *ncm_mset_func_list_new_ns_name (const gchar *ns, const gchar *name, GObject *obj);

gboolean ncm_mset_func_list_has_ns_name (const gchar *ns, const gchar *name);
gboolean ncm_mset_func_list_has_full_name (const gchar *full_name);
GObject *ncm_mset_func_list_peek_obj (NcmMSetFuncList *flist);

G_END_DECLS

#endif /* _NCM_MSET_FUNC_LIST_H_ */

