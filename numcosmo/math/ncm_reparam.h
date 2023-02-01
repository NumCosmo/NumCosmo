/***************************************************************************
 *            ncm_reparam.h
 *
 *  Thu March 08 00:36:24 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

#ifndef _NCM_REPARAM_H_
#define _NCM_REPARAM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_sparam.h>
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>
#include <numcosmo/math/ncm_serialize.h>

G_BEGIN_DECLS

#define NCM_TYPE_REPARAM             (ncm_reparam_get_type ())
#define NCM_REPARAM(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_REPARAM, NcmReparam))
#define NCM_REPARAM_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_REPARAM, NcmReparamClass))
#define NCM_IS_REPARAM(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_REPARAM))
#define NCM_IS_REPARAM_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_REPARAM))
#define NCM_REPARAM_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_REPARAM, NcmReparamClass))

typedef struct _NcmReparamClass NcmReparamClass;
typedef struct _NcmReparam NcmReparam;

struct _NcmModel;

/**
 * NcmReparamV:
 * @reparam: FIXME
 * @model: FIXME
 *
 * FIXME
 */
typedef gboolean (*NcmReparamV) (NcmReparam *reparam, struct _NcmModel *model);

/**
 * NcmReparamJ:
 * @reparam: FIXME
 * @model: FIXME
 * @jac: FIXME
 *
 * FIXME
 */
typedef gboolean (*NcmReparamJ) (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac);

struct _NcmReparamClass
{
  /*< private >*/
  GObjectClass parent_class;
  NcmReparamV old2new;
  NcmReparamV new2old;
  NcmReparamJ jac;
};

struct _NcmReparam
{
  /*< private >*/
  GObject parent_instance;
  guint length;
  NcmVector *new_params;
  GPtrArray *sparams;
  GHashTable *sparams_name_id;
  GType compat_type;
};

GType ncm_reparam_get_type (void) G_GNUC_CONST;

NcmReparam *ncm_reparam_ref (NcmReparam *reparam);
void ncm_reparam_free (NcmReparam *reparam);
void ncm_reparam_clear (NcmReparam **reparam);
GType ncm_reparam_get_compat_type (NcmReparam *reparam);

void ncm_reparam_old2new (NcmReparam *reparam, struct _NcmModel *model);
void ncm_reparam_new2old (NcmReparam *reparam, struct _NcmModel *model);
void ncm_reparam_jac (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac);
void ncm_reparam_grad_old2new (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac, NcmVector *old_grad, NcmVector *new_grad);
void ncm_reparam_M_old2new (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac, NcmMatrix *old_M, NcmMatrix *new_M);

GVariant *ncm_reparam_get_params_desc_dict (NcmReparam *reparam);
void ncm_reparam_set_params_desc_dict (NcmReparam *reparam, GVariant *pdesc_dict);
void ncm_reparam_set_param_desc (NcmReparam *reparam, guint i, NcmSParam *sp);
NcmSParam *ncm_reparam_peek_param_desc (NcmReparam *reparam, guint i);
NcmSParam *ncm_reparam_get_param_desc (NcmReparam *reparam, guint i);
void ncm_reparam_set_param_desc_full (NcmReparam *reparam, guint i, const gchar *name, const gchar *symbol, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_val, NcmParamType ftype);
gboolean ncm_reparam_index_from_name (NcmReparam *reparam, const gchar *param_name, guint *i);

#define NCM_REPARAM_PARAMS_DESC_DICT_TYPE "a{u"NCM_SERIALIZE_OBJECT_TYPE"}"

G_END_DECLS

#endif /* _NCM_REPARAM_H_ */
