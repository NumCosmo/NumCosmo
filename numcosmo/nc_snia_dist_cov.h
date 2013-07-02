/***************************************************************************
 *            nc_snia_dist_cov.h
 *
 *  Mon December 03 19:34:20 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
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

#ifndef _NC_SNIA_DIST_COV_H_
#define _NC_SNIA_DIST_COV_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_model.h>
#include <numcosmo/nc_data_snia_cov.h>
#include <numcosmo/nc_distance.h>

G_BEGIN_DECLS

#define NC_TYPE_SNIA_DIST_COV             (nc_snia_dist_cov_get_type ())
#define NC_SNIA_DIST_COV(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_SNIA_DIST_COV, NcSNIADistCov))
#define NC_SNIA_DIST_COV_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_SNIA_DIST_COV, NcSNIADistCovClass))
#define NC_IS_SNIA_DIST_COV(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_SNIA_DIST_COV))
#define NC_IS_SNIA_DIST_COV_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_SNIA_DIST_COV))
#define NC_SNIA_DIST_COV_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_SNIA_DIST_COV, NcSNIADistCovClass))

typedef struct _NcSNIADistCovClass NcSNIADistCovClass;
typedef struct _NcSNIADistCov NcSNIADistCov;

/**
 * NcSNIADistCovParams:
 * @NC_SNIA_DIST_COV_ALPHA: FIXME
 * @NC_SNIA_DIST_COV_BETA: FIXME
 * @NC_SNIA_DIST_COV_M1: FIXME
 * @NC_SNIA_DIST_COV_M2: FIXME
 *
 * FIXME
 * 
 */
typedef enum _NcSNIADistCovParams
{
  NC_SNIA_DIST_COV_ALPHA = 0,
  NC_SNIA_DIST_COV_BETA,
  NC_SNIA_DIST_COV_M1,
  NC_SNIA_DIST_COV_M2,         /*< private >*/
  NC_SNIA_DIST_COV_SPARAM_LEN, /*< skip >*/
} NcSNIADistCovParams;

struct _NcSNIADistCovClass
{
  /*< private >*/
  NcmModelClass parent_class;
};

struct _NcSNIADistCov
{
  /*< private >*/
  NcmModel parent_instance;
  NcDistance *dist;
};

GType nc_snia_dist_cov_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_snia_dist_cov);

NcSNIADistCov *nc_snia_dist_cov_new (NcDistance *dist);
NcSNIADistCov *nc_snia_dist_cov_ref (NcSNIADistCov *dcov);
void nc_snia_dist_cov_free (NcSNIADistCov *dcov);
void nc_snia_dist_cov_clear (NcSNIADistCov **dcov);

void nc_snia_dist_cov_prepare (NcSNIADistCov *dcov, NcmMSet *mset);
void nc_snia_dist_cov_prepare_if_needed (NcSNIADistCov *dcov, NcmMSet *mset);

void nc_snia_dist_cov_calc (NcSNIADistCov *dcov, NcDataSNIACov *snia_cov, NcmMatrix *cov);
void nc_snia_dist_cov_mean (NcSNIADistCov *dcov, NcHICosmo *cosmo, NcDataSNIACov *snia_cov, NcmVector *y);

#define NC_SNIA_DIST_COV_DEFAULT_ALPHA (1.45)
#define NC_SNIA_DIST_COV_DEFAULT_BETA (3.16)
#define NC_SNIA_DIST_COV_DEFAULT_M1 (-19.1686133146)
#define NC_SNIA_DIST_COV_DEFAULT_M2 (-19.1856133146)
#define NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL (0.0)

#define NC_SNIA_DIST_COV_SIGMA_INT (0)
#define NC_SNIA_DIST_COV_SIGMA_INT_DEFAULT_LEN (4)
#define NC_SNIA_DIST_COV_DEFAULT_SIGMA_INT (0.0989)
#define NC_SNIA_DIST_COV_VPARAM_LEN (1)

G_END_DECLS

#endif /* _NC_SNIA_DIST_COV_H_ */

