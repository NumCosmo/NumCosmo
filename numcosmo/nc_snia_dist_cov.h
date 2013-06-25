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

/**
 * NcSNIADistCovData:
 * @NC_SNIA_DIST_COV_ZCMB: FIXME
 * @NC_SNIA_DIST_COV_ZHE: FIXME
 * @NC_SNIA_DIST_COV_SIGMA_Z: FIXME
 * @NC_SNIA_DIST_COV_MAG: FIXME
 * @NC_SNIA_DIST_COV_SIGMA_MAG: FIXME
 * @NC_SNIA_DIST_COV_WIDTH: FIXME
 * @NC_SNIA_DIST_COV_SIGMA_WIDTH: FIXME
 * @NC_SNIA_DIST_COV_COLOUR: FIXME
 * @NC_SNIA_DIST_COV_SIGMA_COLOUR: FIXME
 * @NC_SNIA_DIST_COV_THIRDPAR: FIXME
 * @NC_SNIA_DIST_COV_SIGMA_THIRDPAR: FIXME
 * @NC_SNIA_DIST_COV_DIAG_MAG_WIDTH: FIXME
 * @NC_SNIA_DIST_COV_DIAG_MAG_COLOUR: FIXME
 * @NC_SNIA_DIST_COV_DIAG_WIDTH_COLOUR: FIXME
 * @NC_SNIA_DIST_COV_ABSMAG_SET: FIXME
 * @NC_SNIA_DIST_COV_VAR_MAG: FIXME
 * @NC_SNIA_DIST_COV_VAR_WIDTH: FIXME
 * @NC_SNIA_DIST_COV_VAR_COLOUR: FIXME
 * @NC_SNIA_DIST_COV_VAR_MAG_WIDTH: FIXME
 * @NC_SNIA_DIST_COV_VAR_MAG_COLOUR: FIXME
 * @NC_SNIA_DIST_COV_VAR_WIDTH_COLOUR: FIXME
 * 
 * FIXME
 * 
 */
typedef enum _NcSNIADistCovData
{
  NC_SNIA_DIST_COV_ZCMB = 0,
  NC_SNIA_DIST_COV_ZHE,
  NC_SNIA_DIST_COV_SIGMA_Z,
  NC_SNIA_DIST_COV_MAG,
  NC_SNIA_DIST_COV_SIGMA_MAG,
  NC_SNIA_DIST_COV_WIDTH,
  NC_SNIA_DIST_COV_SIGMA_WIDTH,
  NC_SNIA_DIST_COV_COLOUR,
  NC_SNIA_DIST_COV_SIGMA_COLOUR,
  NC_SNIA_DIST_COV_THIRDPAR,
  NC_SNIA_DIST_COV_SIGMA_THIRDPAR,
  NC_SNIA_DIST_COV_DIAG_MAG_WIDTH,
  NC_SNIA_DIST_COV_DIAG_MAG_COLOUR,
  NC_SNIA_DIST_COV_DIAG_WIDTH_COLOUR,
  NC_SNIA_DIST_COV_ABSMAG_SET, 
  NC_SNIA_DIST_COV_VAR_MAG,
  NC_SNIA_DIST_COV_VAR_WIDTH,
  NC_SNIA_DIST_COV_VAR_COLOUR,
  NC_SNIA_DIST_COV_VAR_MAG_WIDTH,
  NC_SNIA_DIST_COV_VAR_MAG_COLOUR,
  NC_SNIA_DIST_COV_VAR_WIDTH_COLOUR, /*< private >*/
  NC_SNIA_DIST_COV_TOTAL_LENGTH,     /*< skip >*/
} NcSNIADistCovData;

#define NC_SNIA_DIST_COV_LENGTH NC_SNIA_DIST_COV_ABSMAG_SET

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
  guint mu_len;
  NcmVector *z_cmb;
  NcmVector *z_he;
  NcmVector *mag;
  NcmVector *width;
  NcmVector *colour;
  NcmVector *thirdpar;
  NcmVector *sigma_z;
  NcmVector *sigma_mag;
  NcmVector *sigma_width;
  NcmVector *sigma_colour;
  NcmVector *sigma_thirdpar;
  NcmVector *diag_mag_width;
  NcmVector *diag_mag_colour;
  NcmVector *diag_width_colour;
  NcmMatrix *var_mag;
  NcmMatrix *var_width;
  NcmMatrix *var_colour;
  NcmMatrix *var_mag_width;
  NcmMatrix *var_mag_colour;
  NcmMatrix *var_width_colour;
  NcmVector *sigma_int;
  GArray *dataset;
  gdouble sigma_pecz;
};

GType nc_snia_dist_cov_get_type (void) G_GNUC_CONST;

NCM_MSET_MODEL_DECLARE_ID (nc_snia_dist_cov);

NcSNIADistCov *nc_snia_dist_cov_new (NcDistance *dist, guint mu_len);
NcSNIADistCov *nc_snia_dist_cov_ref (NcSNIADistCov *dcov);
void nc_snia_dist_cov_free (NcSNIADistCov *dcov);
void nc_snia_dist_cov_clear (NcSNIADistCov **dcov);

void nc_snia_dist_cov_prepare (NcSNIADistCov *dcov, NcmMSet *mset);
void nc_snia_dist_cov_prepare_if_needed (NcSNIADistCov *dcov, NcmMSet *mset);

void nc_snia_dist_cov_set_size (NcSNIADistCov *dcov, guint mu_len);
void nc_snia_dist_cov_calc (NcSNIADistCov *dcov, NcmMatrix *cov);

void nc_snia_dist_cov_mean (NcSNIADistCov *dcov, NcHICosmo *cosmo, NcmVector *y);

void nc_snia_dist_cov_load_txt (NcSNIADistCov *dcov, const gchar *filename);

#ifdef NUMCOSMO_HAVE_CFITSIO
void nc_snia_dist_cov_load (NcSNIADistCov *dcov, const gchar *filename);
void nc_snia_dist_cov_save (NcSNIADistCov *dcov, const gchar *filename, gboolean overwrite);
#endif /* NUMCOSMO_HAVE_CFITSIO */

#define NC_SNIA_DIST_COV_DEFAULT_ALPHA (1.45)
#define NC_SNIA_DIST_COV_DEFAULT_BETA (3.16)
#define NC_SNIA_DIST_COV_DEFAULT_M1 (-19.1686133146)
#define NC_SNIA_DIST_COV_DEFAULT_M2 (-19.1856133146)
#define NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL (0.0)

#define NC_SNIA_DIST_COV_SIGMA_INT (0)
#define NC_SNIA_DIST_COV_SIGMA_INT_DEFAULT_LEN (4)
#define NC_SNIA_DIST_COV_DEFAULT_SIGMA_INT (0.0989)
#define NC_SNIA_DIST_COV_VPARAM_LEN (1)

#define NC_SNIA_DIST_COV_DATA_GROUP "Supernovae Ia Data"
#define NC_SNIA_DIST_COV_DATA_LEN_KEY "data-length"
#define NC_SNIA_DIST_COV_DATA_KEY "snia-data"

#define NC_SNIA_DIST_COV_MAG_KEY "magnitude"
#define NC_SNIA_DIST_COV_WIDTH_KEY "width"
#define NC_SNIA_DIST_COV_COLOUR_KEY "colour"

#define NC_SNIA_DIST_COV_MAG_WIDTH_KEY "magnitude-width"
#define NC_SNIA_DIST_COV_MAG_COLOUR_KEY "magnitude-colour"
#define NC_SNIA_DIST_COV_WIDTH_COLOUR_KEY "width-colour"

G_END_DECLS

#endif /* _NC_SNIA_DIST_COV_H_ */

