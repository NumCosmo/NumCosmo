/***************************************************************************
 *            nc_cluster_photoz_gauss_global.h
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

#ifndef _NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_H_
#define _NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/lss/nc_cluster_redshift.h>

G_BEGIN_DECLS

#define NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL             (nc_cluster_photoz_gauss_global_get_type ())
#define NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL, NcClusterPhotozGaussGlobal))
#define NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL, NcClusterPhotozGaussGlobalClass))
#define NC_IS_CLUSTER_PHOTOZ_GAUSS_GLOBAL(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL))
#define NC_IS_CLUSTER_PHOTOZ_GAUSS_GLOBAL_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL))
#define NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL, NcClusterPhotozGaussGlobalClass))

typedef struct _NcClusterPhotozGaussGlobalClass NcClusterPhotozGaussGlobalClass;
typedef struct _NcClusterPhotozGaussGlobal NcClusterPhotozGaussGlobal;

/**
 * NcClusterPhotozGaussGlobalSParams:
 * @NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_Z_BIAS: FIXME
 * @NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_SIGMA0: standard deviation of the gaussian distribution
 *
 * FIXME
 */
typedef enum /*< enum,underscore_name=NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_SPARAMS >*/
{
  NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_Z_BIAS = 0,
  NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_SIGMA0,  
  /* < private > */
  NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_SPARAM_LEN, /*< skip >*/
} NcClusterPhotozGaussGlobalSParams;

struct _NcClusterPhotozGaussGlobalClass
{
  /*< private >*/
  NcClusterRedshiftClass parent_class;
};

struct _NcClusterPhotozGaussGlobal
{
  /*< private >*/
  NcClusterRedshift parent_instance;
  gdouble pz_min;
  gdouble pz_max;
};

GType nc_cluster_photoz_gauss_global_get_type (void) G_GNUC_CONST;

NcClusterRedshift *nc_cluster_photoz_gauss_global_new (gdouble pz_min, gdouble pz_max, gdouble z_bias, gdouble sigma0);
void nc_cluster_photoz_gauss_global_set_z_bias (NcClusterPhotozGaussGlobal *pzg_global, gdouble z_bias);
gdouble nc_cluster_photoz_gauss_global_get_z_bias (const NcClusterPhotozGaussGlobal *pzg_global);
void nc_cluster_photoz_gauss_global_set_sigma0 (NcClusterPhotozGaussGlobal *pzg_global, gdouble sigma0);
gdouble nc_cluster_photoz_gauss_global_get_sigma0 (const NcClusterPhotozGaussGlobal *pzg_global);

#define NC_CLUSTER_REDSHIFT_PHOTOZ_GAUSS_GLOBAL_DEFAULT_BIAS (0.0)
#define NC_CLUSTER_REDSHIFT_PHOTOZ_GAUSS_GLOBAL_DEFAULT_SIGMA0 (0.03)
#define NC_CLUSTER_REDSHIFT_PHOTOZ_GAUSS_GLOBAL_DEFAULT_PARAMS_ABSTOL (0.0)

G_END_DECLS

#endif /* _NC_CLUSTER_PHOTOZ_GAUSS_GLOBAL_H_ */
