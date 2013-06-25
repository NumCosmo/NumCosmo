/***************************************************************************
 *            ncm_vparam.h
 *
 *  Thu May 10 15:50:55 2012
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

#ifndef _NCM_VPARAM_H_
#define _NCM_VPARAM_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_sparam.h>

G_BEGIN_DECLS

#define NCM_TYPE_VPARAM             (ncm_vparam_get_type ())
#define NCM_VPARAM(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NCM_TYPE_VPARAM, NcmVParam))
#define NCM_VPARAM_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NCM_TYPE_VPARAM, NcmVParamClass))
#define NCM_IS_VPARAM(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NCM_TYPE_VPARAM))
#define NCM_IS_VPARAM_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NCM_TYPE_VPARAM))
#define NCM_VPARAM_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NCM_TYPE_VPARAM, NcmVParamClass))

typedef struct _NcmVParamClass NcmVParamClass;
typedef struct _NcmVParam NcmVParam;

struct _NcmVParamClass
{
  /*< private >*/
  GObjectClass parent_class;
};

struct _NcmVParam
{
  /*< private >*/
  GObject parent_instance;
  guint len;
  NcmSParam *default_sparam;
  GPtrArray *sparam;
};

GType ncm_vparam_get_type (void) G_GNUC_CONST;

NcmVParam *ncm_vparam_new (guint len, NcmSParam *default_param);
NcmVParam *ncm_vparam_full_new (guint len, gchar *name, gchar *symbol, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_val, NcmParamType ftype);
NcmVParam *ncm_vparam_copy (NcmVParam *vparam);
void ncm_vparam_free (NcmVParam *vparam);
void ncm_vparam_clear (NcmVParam **vparam);

void ncm_vparam_set_len (NcmVParam *vparam, guint len);
guint ncm_vparam_get_len (NcmVParam *vparam);
void ncm_vparam_set_sparam (NcmVParam *vparam, guint n, NcmSParam *spn);
void ncm_vparam_set_sparam_full (NcmVParam *vparam, guint n, gchar *name, gchar *symbol, gdouble lower_bound, gdouble upper_bound, gdouble scale, gdouble abstol, gdouble default_val, NcmParamType ftype);
NcmSParam *ncm_vparam_peek_sparam (const NcmVParam *vparam, guint n);
NcmSParam *ncm_vparam_get_sparam (NcmVParam *vparam, guint n);

void ncm_vparam_set_lower_bound (NcmVParam *vparam, guint n, const gdouble lb);
void ncm_vparam_set_upper_bound (NcmVParam *vparam, guint n, const gdouble ub);
void ncm_vparam_set_scale (NcmVParam *vparam, guint n, const gdouble scale);
void ncm_vparam_set_absolute_tolerance (NcmVParam *vparam, guint n, const gdouble abstol);
void ncm_vparam_set_default_value (NcmVParam *vparam, guint n, const gdouble default_val);
void ncm_vparam_set_fit_type (NcmVParam *vparam, guint n, const NcmParamType ftype);

gdouble ncm_vparam_get_lower_bound (const NcmVParam *vparam, guint n);
gdouble ncm_vparam_get_upper_bound (const NcmVParam *vparam, guint n);
gdouble ncm_vparam_get_scale (const NcmVParam *vparam, guint n);
gdouble ncm_vparam_get_absolute_tolerance (const NcmVParam *vparam, guint n);
gdouble ncm_vparam_get_default_value (const NcmVParam *vparam, guint n);
NcmParamType ncm_vparam_get_fit_type (const NcmVParam *vparam, guint n);

G_END_DECLS

#endif /* _NCM_VPARAM_H_ */
