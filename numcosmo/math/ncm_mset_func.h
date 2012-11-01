/***************************************************************************
 *            ncm_mset_func.h
 *
 *  Wed June 06 15:32:36 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#ifndef _NCM_MSET_FUNC_H_
#define _NCM_MSET_FUNC_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_mset.h>

G_BEGIN_DECLS

#define NCM_TYPE_GMSET_FUNC             (ncm_mset_func_get_type ())
#define NCM_MSET_FUNC(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_GMSET_FUNC, NcmMSetFunc))
#define NCM_MSET_FUNC_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_GMSET_FUNC, NcmMSetFuncClass))
#define NCM_IS_GMSET_FUNC(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_GMSET_FUNC))
#define NCM_IS_GMSET_FUNC_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_GMSET_FUNC))
#define NCM_MSET_FUNC_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_GMSET_FUNC, NcmMSetFuncClass))

typedef struct _NcmMSetFuncClass NcmMSetFuncClass;
typedef struct _NcmMSetFunc NcmMSetFunc;
typedef void (*NcmMSetFuncN) (NcmMSet *mset, gpointer obj, const gdouble *x, gdouble *f);

struct _NcmMSetFuncClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmMSetFunc
{
  /*< private >*/
  GObject parent_instance;
  NcmMSetFuncN func;
  gpointer obj;
  GDestroyNotify free;
  guint np;
  guint dim;
};

GType ncm_mset_func_get_type (void) G_GNUC_CONST;

NcmMSetFunc *ncm_mset_func_new (NcmMSetFuncN func, guint np, guint dim, gpointer obj, GDestroyNotify free);
NcmMSetFunc *ncm_mset_func_ref (NcmMSetFunc *func);

void ncm_mset_func_free (NcmMSetFunc *func);

GPtrArray *ncm_mset_func_array_new (void);

gdouble ncm_mset_func_eval0 (NcmMSetFunc *func, NcmMSet *mset);
gdouble ncm_mset_func_eval1 (NcmMSetFunc *func, NcmMSet *mset, const gdouble x);

G_END_DECLS

#endif /* _NCM_MSET_FUNC_H_ */
