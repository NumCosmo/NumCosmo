/***************************************************************************
 *            ncm_model_builder.h
 *
 *  Fri November 06 12:19:13 2015
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_model_builder.h
 * Copyright (C) 2015 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
#include <numcosmo/math/ncm_obj_array.h>

G_BEGIN_DECLS

#define NCM_TYPE_MODEL_BUILDER (ncm_model_builder_get_type ())

G_DECLARE_FINAL_TYPE (NcmModelBuilder, ncm_model_builder, NCM, MODEL_BUILDER, GObject)

NcmModelBuilder *ncm_model_builder_new (GType ptype, const gchar *name, const gchar *desc);
NcmModelBuilder *ncm_model_builder_ref (NcmModelBuilder *mb);

void ncm_model_builder_add_sparam_obj (NcmModelBuilder *mb, NcmSParam *sparam);
void ncm_model_builder_add_vparam_obj (NcmModelBuilder *mb, NcmVParam *vparam);

void ncm_model_builder_add_sparam (NcmModelBuilder *mb, const gchar *symbol, const gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt);
void ncm_model_builder_add_vparam (NcmModelBuilder *mb, guint default_length, const gchar *symbol, const gchar *name, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_value, NcmParamType ppt);

GType ncm_model_builder_create (NcmModelBuilder *mb);

void ncm_model_builder_add_sparams (NcmModelBuilder *mb, NcmObjArray *sparams);
NcmObjArray *ncm_model_builder_get_sparams (NcmModelBuilder *mb);


G_END_DECLS

#endif /* _NCM_MODEL_BUILDER_H_ */

