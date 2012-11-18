/***************************************************************************
 *            nc_data_distance_mu.h
 *
 *  Thu Apr 22 10:37:39 2010
 *  Copyright  2010  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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

#ifndef _NC_DATA_DIST_MU_H_
#define _NC_DATA_DIST_MU_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/math/ncm_data_gauss_diag.h>
#include <numcosmo/nc_distance.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_DIST_MU             (nc_data_dist_mu_get_type ())
#define NC_DATA_DIST_MU(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_DIST_MU, NcDataDistMu))
#define NC_DATA_DIST_MU_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_DIST_MU, NcDataDistMuClass))
#define NC_IS_DATA_DIST_MU(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_DIST_MU))
#define NC_IS_DATA_DIST_MU_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_DIST_MU))
#define NC_DATA_DIST_MU_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_DIST_MU, NcDataDistMuClass))

typedef struct _NcDataDistMuClass NcDataDistMuClass;
typedef struct _NcDataDistMu NcDataDistMu;

/**
 * NcDataDistMuId:
 * @NC_DATA_DIST_MU_SNIA_GOLD_157: FIXME
 * @NC_DATA_DIST_MU_SNIA_GOLD_182: FIXME
 * @NC_DATA_DIST_MU_SNIA_GOLD_182_FULL: FIXME
 * @NC_DATA_DIST_MU_SNIA_ESSENCE: FIXME
 * @NC_DATA_DIST_MU_SNIA_LEGACY: FIXME
 * @NC_DATA_DIST_MU_SNIA_UNION: FIXME
 * @NC_DATA_DIST_MU_SNIA_CfA3: FIXME
 * @NC_DATA_DIST_MU_SNIA_UNION2: FIXME
 * @NC_DATA_DIST_MU_SNIA_UNION2_1: FIXME
 */
typedef enum _NcDataDistMuId
{
  NC_DATA_DIST_MU_SNIA_GOLD_157 = 0,
  NC_DATA_DIST_MU_SNIA_GOLD_182,
  NC_DATA_DIST_MU_SNIA_GOLD_182_FULL,
  NC_DATA_DIST_MU_SNIA_ESSENCE,
  NC_DATA_DIST_MU_SNIA_LEGACY,
  NC_DATA_DIST_MU_SNIA_UNION,
  NC_DATA_DIST_MU_SNIA_CfA3,
  NC_DATA_DIST_MU_SNIA_UNION2,
  NC_DATA_DIST_MU_SNIA_UNION2_1, /*< private >*/
  NC_DATA_DIST_MU_NSAMPLES,      /*< skip >*/
} NcDataDistMuId;

struct _NcDataDistMuClass
{
  /*< private >*/
  NcmDataGaussDiagClass parent_class;
};

struct _NcDataDistMu
{
  /*< private >*/
  NcmDataGaussDiag parent_instance;
  NcDistance *dist;
  NcmVector *x;
  NcDataDistMuId id;
};

GType nc_data_dist_mu_get_type (void) G_GNUC_CONST;

NcmData *nc_data_dist_mu_new (NcDistance *dist, NcDataDistMuId id);

void nc_data_dist_mu_set_size (NcDataDistMu *dist_mu, guint np);
guint nc_data_dist_mu_get_size (NcDataDistMu *dist_mu);

void nc_data_dist_mu_set_sample (NcDataDistMu *dist_mu, NcDataDistMuId id);
NcDataDistMuId nc_data_dist_mu_get_sample (NcDataDistMu *dist_mu);

G_END_DECLS

#endif /* _NC_DATA_DIST_MU_H_ */

