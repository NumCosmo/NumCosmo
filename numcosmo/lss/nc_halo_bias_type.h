/***************************************************************************
 *            nc_halo_bias_type.h
 *
 *  Tue June 28 15:41:57 2011
 *  Copyright  2011  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_HALO_BIAS_TYPE_H_
#define _NC_HALO_BIAS_TYPE_H_

#include <glib-object.h>

G_BEGIN_DECLS

#define NC_TYPE_HALO_BIAS_TYPE             (nc_halo_bias_type_get_type ())
#define NC_HALO_BIAS_TYPE(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_HALO_BIAS_TYPE, NcHaloBiasType))
#define NC_HALO_BIAS_TYPE_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_HALO_BIAS_TYPE, NcHaloBiasTypeClass))
#define NC_IS_HALO_BIAS_TYPE(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_HALO_BIAS_TYPE))
#define NC_IS_HALO_BIAS_TYPE_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_HALO_BIAS_TYPE))
#define NC_HALO_BIAS_TYPE_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_HALO_BIAS_TYPE, NcHaloBiasTypeClass))

typedef struct _NcHaloBiasTypeClass NcHaloBiasTypeClass;
typedef struct _NcHaloBiasType NcHaloBiasType;

struct _NcHaloBiasTypeClass
{
  /*< private >*/
  GObjectClass parent_class;
  gdouble (*eval) (NcHaloBiasType *biasf, gdouble sigma, gdouble z); 
};

struct _NcHaloBiasType
{
  /*< private >*/
  GObject parent_instance;
};

GType nc_halo_bias_type_get_type (void) G_GNUC_CONST;

NcHaloBiasType *nc_halo_bias_type_new_from_name (gchar *bias_name);
gdouble nc_halo_bias_type_eval (NcHaloBiasType *biasf, gdouble sigma, gdouble z);
void nc_halo_bias_type_free (NcHaloBiasType *biasf);
void nc_halo_bias_type_clear (NcHaloBiasType **biasf);

G_END_DECLS

#endif /* _NC_HALO_BIAS_TYPE_H_ */
