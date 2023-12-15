/***************************************************************************
 *            ncm_mset_func.h
 *
 *  Wed June 06 15:32:36 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_mset.h>
#include <numcosmo/math/ncm_diff.h>

G_BEGIN_DECLS

#define NCM_TYPE_MSET_FUNC (ncm_mset_func_get_type ())

G_DECLARE_DERIVABLE_TYPE (NcmMSetFunc, ncm_mset_func, NCM, MSET_FUNC, GObject)

typedef void (*NcmMSetFuncN) (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res);

struct _NcmMSetFuncClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcmMSetFuncN eval;

  /* Padding to allow 18 virtual functions without breaking ABI. */
  gpointer padding[17];
};

NcmMSetFunc *ncm_mset_func_ref (NcmMSetFunc *func);

void ncm_mset_func_free (NcmMSetFunc *func);
void ncm_mset_func_clear (NcmMSetFunc **func);

GPtrArray *ncm_mset_func_array_new (void);

void ncm_mset_func_eval (NcmMSetFunc *func, NcmMSet *mset, gdouble *x, gdouble *res);
gdouble ncm_mset_func_eval_nvar (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x);
gdouble ncm_mset_func_eval0 (NcmMSetFunc *func, NcmMSet *mset);
gdouble ncm_mset_func_eval1 (NcmMSetFunc *func, NcmMSet *mset, const gdouble x);
void ncm_mset_func_eval_vector (NcmMSetFunc *func, NcmMSet *mset, NcmVector *x_v, NcmVector *res_v);

void ncm_mset_func_set_eval_x (NcmMSetFunc *func, gdouble *x, guint len);

gboolean ncm_mset_func_is_scalar (NcmMSetFunc *func);
gboolean ncm_mset_func_is_vector (NcmMSetFunc *func, guint dim);
gboolean ncm_mset_func_is_const (NcmMSetFunc *func);
gboolean ncm_mset_func_has_nvar (NcmMSetFunc *func, guint nvar);

guint ncm_mset_func_get_nvar (NcmMSetFunc *func);
guint ncm_mset_func_get_dim (NcmMSetFunc *func);

const gchar *ncm_mset_func_peek_name (NcmMSetFunc *func);
const gchar *ncm_mset_func_peek_symbol (NcmMSetFunc *func);
const gchar *ncm_mset_func_peek_ns (NcmMSetFunc *func);
const gchar *ncm_mset_func_peek_desc (NcmMSetFunc *func);

const gchar *ncm_mset_func_peek_uname (NcmMSetFunc *func);
const gchar *ncm_mset_func_peek_usymbol (NcmMSetFunc *func);

void ncm_mset_func_numdiff_fparams (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, NcmVector **out);

G_END_DECLS

#endif /* _NCM_MSET_FUNC_H_ */
