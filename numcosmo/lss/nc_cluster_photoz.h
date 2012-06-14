/***************************************************************************
 *            nc_cluster_photoz.h
 *
 *  Tue Apr 20 10:59:01 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com> & <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
/*
 * numcosmo
 * Copyright (C) Mariana Penna lima 2012 <pennalima@gmail.com>
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

#ifndef _NC_CLUSTER_PHOTOZ_H_
#define _NC_CLUSTER_PHOTOZ_H_

#include <glib-object.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_PHOTOZ             (nc_cluster_photoz_get_type ())
#define NC_CLUSTER_PHOTOZ(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_PHOTOZ, NcClusterPhotoz))
#define NC_CLUSTER_PHOTOZ_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_PHOTOZ, NcClusterPhotozClass))
#define NC_IS_CLUSTER_PHOTOZ(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_PHOTOZ))
#define NC_IS_CLUSTER_PHOTOZ_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_PHOTOZ))
#define NC_CLUSTER_PHOTOZ_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_PHOTOZ, NcClusterPhotozClass))

typedef struct _NcClusterPhotozClass NcClusterPhotozClass;
typedef struct _NcClusterPhotoz NcClusterPhotoz;

struct _NcClusterPhotozClass
{
  /*< private >*/
  GObjectClass parent_class;
  gdouble (*dist_eval) (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_photo, gdouble z_real);
  void (*resample) (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_real, gdouble* result);
  void (*integ_limits) (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_photo, gdouble *z_lower, gdouble *z_upper);
};

struct _NcClusterPhotoz
{
  /*< private >*/
  GObject parent_instance;
};

GType nc_cluster_photoz_get_type (void) G_GNUC_CONST;

NcClusterPhotoz *nc_cluster_photoz_new_from_name (gchar *photoz_name);
gdouble nc_cluster_photoz_dist_eval (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_photo, gdouble z_real);
void nc_cluster_photoz_resample (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_real, gdouble *result);
void nc_cluster_photoz_integ_limits (NcClusterPhotoz *photo, gdouble *pz_data, gdouble z_photo, gdouble *z_lower, gdouble *z_upper);
void nc_cluster_photoz_free (NcClusterPhotoz *photo);

G_END_DECLS

#endif /* _NC_CLUSTER_PHOTOZ_H_ */
