/***************************************************************************
 *            nc_cluster_photoz_gauss.h
 *
 *  Tue Apr 20 10:59:01 2010
 *  Copyright  2010  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com> & <vitenti@uel.br>
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

#ifndef _NC_CLUSTER_PHOTOZ_GAUSS_H_
#define _NC_CLUSTER_PHOTOZ_GAUSS_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_redshift.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_PHOTOZ_GAUSS             (nc_cluster_photoz_gauss_get_type ())
#define NC_CLUSTER_PHOTOZ_GAUSS(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_PHOTOZ_GAUSS, NcClusterPhotozGauss))
#define NC_CLUSTER_PHOTOZ_GAUSS_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_PHOTOZ_GAUSS, NcClusterPhotozGaussClass))
#define NC_IS_CLUSTER_PHOTOZ_GAUSS(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_PHOTOZ_GAUSS))
#define NC_IS_CLUSTER_PHOTOZ_GAUSS_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_PHOTOZ_GAUSS))
#define NC_CLUSTER_PHOTOZ_GAUSS_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_PHOTOZ_GAUSS, NcClusterPhotozGaussClass))

typedef struct _NcClusterPhotozGaussClass NcClusterPhotozGaussClass;
typedef struct _NcClusterPhotozGauss NcClusterPhotozGauss;

struct _NcClusterPhotozGaussClass
{
  /*< private >*/
  NcClusterRedshiftClass parent_class;
};

struct _NcClusterPhotozGauss
{
  /*< private >*/
  NcClusterRedshift parent_instance;
  gdouble pz_max;
  gdouble pz_min;
};

GType nc_cluster_photoz_gauss_get_type (void) G_GNUC_CONST;

NcClusterRedshift *nc_cluster_photoz_gauss_new (void);

#define NC_CLUSTER_PHOTOZ_GAUSS_BIAS  (0)
#define NC_CLUSTER_PHOTOZ_GAUSS_SIGMA (1)

G_END_DECLS

#endif /* _NC_CLUSTER_PHOTOZ_GAUSS_H_ */

