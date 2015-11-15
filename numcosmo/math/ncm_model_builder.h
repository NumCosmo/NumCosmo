/***************************************************************************
 *            ncm_model_builder.h
 *
 *  Fri November 06 12:19:13 2015
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_model_builder.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NCM_MODEL_BUILDER_H_
#define _NCM_MODEL_BUILDER_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_cfg.h>
#include <numcosmo/math/ncm_sparam.h>
#include <numcosmo/math/ncm_vparam.h>

G_BEGIN_DECLS

#define NCM_TYPE_MODEL_BUILDER             (ncm_model_builder_get_type ())
#define NCM_MODEL_BUILDER(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_MODEL_BUILDER, NcmModelBuilder))
#define NCM_MODEL_BUILDER_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_MODEL_BUILDER, NcmModelBuilderClass))
#define NCM_IS_MODEL_BUILDER(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_MODEL_BUILDER))
#define NCM_IS_MODEL_BUILDER_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_MODEL_BUILDER))
#define NCM_MODEL_BUILDER_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_MODEL_BUILDER, NcmModelBuilderClass))

typedef struct _NcmModelBuilderClass NcmModelBuilderClass;
typedef struct _NcmModelBuilder NcmModelBuilder;

struct _NcmModelBuilderClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmModelBuilder
{
  /*< private >*/
  GObject parent_instance;
  gchar *name;
  gchar *desc;
  GType ptype;
  GType type;
  GPtrArray *sparams;
  GPtrArray *vparams;
  gboolean created;
};

GType ncm_model_builder_get_type (void) G_GNUC_CONST;

NcmModelBuilder *ncm_model_builder_new (GType ptype, const gchar *name, const gchar *desc);
NcmModelBuilder *ncm_model_builder_ref (NcmModelBuilder *mb);

void ncm_model_builder_add_sparam_obj (NcmModelBuilder *mb, NcmSParam *sparam);
void ncm_model_builder_add_vparam_obj (NcmModelBuilder *mb, NcmVParam *vparam);

void ncm_model_builder_add_sparam (NcmModelBuilder *mb, const gchar *symbol, const gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt);
void ncm_model_builder_add_vparam (NcmModelBuilder *mb, guint default_length, const gchar *symbol, const gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt);

GType ncm_model_builder_create (NcmModelBuilder *mb);

G_END_DECLS

#endif /* _NCM_MODEL_BUILDER_H_ */
