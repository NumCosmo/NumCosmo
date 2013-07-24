/***************************************************************************
 *            ncm_reparam.h
 *
 *  Thu March 08 00:36:24 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
#include <numcosmo/math/ncm_vector.h>
#include <numcosmo/math/ncm_matrix.h>

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
 * @src: FIXME
 * @dest: FIXME
 *
 * FIXME
 */
typedef gboolean (*NcmReparamV) (NcmReparam *reparam, struct _NcmModel *model, NcmVector *src, NcmVector *dest);

/**
 * NcmReparamJ:
 * @reparam: FIXME
 * @model: FIXME
 * @jac: FIXME
 *
 * FIXME
 */
typedef gboolean (*NcmReparamJ) (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac);

struct _NcmReparam
{
  /*< private >*/
  GObject parent_instance;
  guint length;
  NcmVector *new_params;
  GPtrArray *sparams;
};

struct _NcmReparamClass
{
  /*< private >*/
  GObjectClass parent_class;
  void (*copyto) (NcmReparam *reparam, NcmReparam *reparam_dest);
  NcmReparam *(*copy) (NcmReparam *reparam);
  NcmReparamV old2new;
  NcmReparamV new2old;
  NcmReparamJ jac;
};

GType ncm_reparam_get_type (void) G_GNUC_CONST;

NcmReparam *ncm_reparam_copy (NcmReparam *reparam);
NcmReparam *ncm_reparam_ref (NcmReparam *reparam);
void ncm_reparam_copyto (NcmReparam *reparam, NcmReparam *reparam_dest);
void ncm_reparam_old2new (NcmReparam *reparam, struct _NcmModel *model, NcmVector *src, NcmVector *dest);
void ncm_reparam_new2old (NcmReparam *reparam, struct _NcmModel *model, NcmVector *src, NcmVector *dest);
void ncm_reparam_jac (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac);
void ncm_reparam_grad_old2new (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac, NcmVector *old_grad, NcmVector *new_grad);
void ncm_reparam_M_old2new (NcmReparam *reparam, struct _NcmModel *model, NcmMatrix *jac, NcmMatrix *old_M, NcmMatrix *new_M);

void ncm_reparam_free (NcmReparam *reparam);
void ncm_reparam_clear (NcmReparam **reparam);

G_END_DECLS

#endif /* _NCM_REPARAM_H_ */
