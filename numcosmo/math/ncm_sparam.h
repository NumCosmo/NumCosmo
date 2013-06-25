/***************************************************************************
 *            ncm_sparam.h
 *
 *  Fri February 24 20:14:00 2012
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

#ifndef _NCM_SPARAM_H_
#define _NCM_SPARAM_H_

#include <glib-object.h>
#include <numcosmo/build_cfg.h>

G_BEGIN_DECLS

#define NCM_TYPE_SPARAM             (ncm_sparam_get_type ())
#define NCM_SPARAM(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_SPARAM, NcmSParam))
#define NCM_SPARAM_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_SPARAM, NcmSParamClass))
#define NCM_IS_SPARAM(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_SPARAM))
#define NCM_IS_SPARAM_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_SPARAM))
#define NCM_SPARAM_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_SPARAM, NcmSParamClass))

typedef struct _NcmSParamClass NcmSParamClass;
typedef struct _NcmSParam NcmSParam;

struct _NcmSParamClass
{
  /*< private >*/
  GObjectClass parent_class;
};

/**
 * NcmParamType:
 * @NCM_PARAM_TYPE_FREE: FIXME
 * @NCM_PARAM_TYPE_FIXED: FIXME
 *
 * FIXME
 */
typedef enum _NcmParamType
{
  NCM_PARAM_TYPE_FREE = 0,
  NCM_PARAM_TYPE_FIXED,
} NcmParamType;

struct _NcmSParam
{
  /*< private >*/
  GObject parent_instance;
  gchar *name;
  gchar *symbol;
  gdouble lower_bound;
  gdouble upper_bound;
  gdouble scale;
  gdouble abstol;
  gdouble default_val;
  NcmParamType ftype;
};

GType ncm_sparam_get_type (void) G_GNUC_CONST;

NcmSParam *ncm_sparam_new (gchar *name, gchar *symbol, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_val, NcmParamType ftype);
NcmSParam *ncm_sparam_copy (NcmSParam *sparam);
NcmSParam *ncm_sparam_ref (NcmSParam *sparam);
void ncm_sparam_free (NcmSParam *sparam);
void ncm_sparam_clear (NcmSParam **sparam);

void ncm_sparam_set_lower_bound (NcmSParam *sparam, const gdouble lb);
void ncm_sparam_set_upper_bound (NcmSParam *sparam, const gdouble ub);
void ncm_sparam_set_scale (NcmSParam *sparam, const gdouble scale);
void ncm_sparam_set_absolute_tolerance (NcmSParam *sparam, const gdouble abstol);
void ncm_sparam_set_default_value (NcmSParam *sparam, const gdouble default_val);
void ncm_sparam_set_fit_type (NcmSParam *sparam, const NcmParamType ftype);

void ncm_sparam_take_name (NcmSParam *sparam, gchar *name);
void ncm_sparam_take_symbol (NcmSParam *sparam, gchar *symbol);
const gchar *ncm_sparam_name (const NcmSParam *sparam);
const gchar *ncm_sparam_symbol (const NcmSParam *sparam);

gdouble ncm_sparam_get_lower_bound (const NcmSParam *sparam);
gdouble ncm_sparam_get_upper_bound (const NcmSParam *sparam);
gdouble ncm_sparam_get_scale (const NcmSParam *sparam);
gdouble ncm_sparam_get_absolute_tolerance (const NcmSParam *sparam);
gdouble ncm_sparam_get_default_value (const NcmSParam *sparam);
NcmParamType ncm_sparam_get_fit_type (const NcmSParam *sparam);

G_END_DECLS

#endif /* _NCM_SPARAM_H_ */
