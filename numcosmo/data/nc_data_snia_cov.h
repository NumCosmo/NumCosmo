/***************************************************************************
 *            nc_data_snia_cov.h
 *
 *  Sat December 08 15:58:20 2012
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

#ifndef _NC_DATA_SNIA_COV_H_
#define _NC_DATA_SNIA_COV_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>
#include <numcosmo/math/ncm_data_gauss_cov.h>
#include <numcosmo/math/ncm_model_ctrl.h>
#include <numcosmo/data/nc_data_snia.h>

G_BEGIN_DECLS

#define NC_TYPE_DATA_SNIA_COV             (nc_data_snia_cov_get_type ())
#define NC_DATA_SNIA_COV(obj)             (G_TYPE_CHECK_INSTANCE_CAST ((obj), NC_TYPE_DATA_SNIA_COV, NcDataSNIACov))
#define NC_DATA_SNIA_COV_CLASS(klass)     (G_TYPE_CHECK_CLASS_CAST ((klass), NC_TYPE_DATA_SNIA_COV, NcDataSNIACovClass))
#define NC_IS_DATA_SNIA_COV(obj)          (G_TYPE_CHECK_INSTANCE_TYPE ((obj), NC_TYPE_DATA_SNIA_COV))
#define NC_IS_DATA_SNIA_COV_CLASS(klass)  (G_TYPE_CHECK_CLASS_TYPE ((klass), NC_TYPE_DATA_SNIA_COV))
#define NC_DATA_SNIA_COV_GET_CLASS(obj)   (G_TYPE_INSTANCE_GET_CLASS ((obj), NC_TYPE_DATA_SNIA_COV, NcDataSNIACovClass))

typedef struct _NcDataSNIACovClass NcDataSNIACovClass;
typedef struct _NcDataSNIACov NcDataSNIACov;

/**
 * NcDataSNIACovData:
 * @NC_DATA_SNIA_COV_ZCMB: FIXME
 * @NC_DATA_SNIA_COV_ZHE: FIXME
 * @NC_DATA_SNIA_COV_SIGMA_Z: FIXME
 * @NC_DATA_SNIA_COV_MAG: FIXME
 * @NC_DATA_SNIA_COV_SIGMA_MAG: FIXME
 * @NC_DATA_SNIA_COV_WIDTH: FIXME
 * @NC_DATA_SNIA_COV_SIGMA_WIDTH: FIXME
 * @NC_DATA_SNIA_COV_COLOUR: FIXME
 * @NC_DATA_SNIA_COV_SIGMA_COLOUR: FIXME
 * @NC_DATA_SNIA_COV_THIRDPAR: FIXME
 * @NC_DATA_SNIA_COV_SIGMA_THIRDPAR: FIXME
 * @NC_DATA_SNIA_COV_DIAG_MAG_WIDTH: FIXME
 * @NC_DATA_SNIA_COV_DIAG_MAG_COLOUR: FIXME
 * @NC_DATA_SNIA_COV_DIAG_WIDTH_COLOUR: FIXME
 * @NC_DATA_SNIA_COV_ABSMAG_SET: FIXME
 * @NC_DATA_SNIA_COV_VAR_MAG: FIXME
 * @NC_DATA_SNIA_COV_VAR_WIDTH: FIXME
 * @NC_DATA_SNIA_COV_VAR_COLOUR: FIXME
 * @NC_DATA_SNIA_COV_VAR_MAG_WIDTH: FIXME
 * @NC_DATA_SNIA_COV_VAR_MAG_COLOUR: FIXME
 * @NC_DATA_SNIA_COV_VAR_WIDTH_COLOUR: FIXME
 * 
 * Data ordering of Version 0 (V0) data format.
 * 
 */
typedef enum _NcDataSNIACovData
{
  NC_DATA_SNIA_COV_ZCMB = 0,
  NC_DATA_SNIA_COV_ZHE,
  NC_DATA_SNIA_COV_SIGMA_Z,
  NC_DATA_SNIA_COV_MAG,
  NC_DATA_SNIA_COV_SIGMA_MAG,
  NC_DATA_SNIA_COV_WIDTH,
  NC_DATA_SNIA_COV_SIGMA_WIDTH,
  NC_DATA_SNIA_COV_COLOUR,
  NC_DATA_SNIA_COV_SIGMA_COLOUR,
  NC_DATA_SNIA_COV_THIRDPAR,
  NC_DATA_SNIA_COV_SIGMA_THIRDPAR,
  NC_DATA_SNIA_COV_DIAG_MAG_WIDTH,
  NC_DATA_SNIA_COV_DIAG_MAG_COLOUR,
  NC_DATA_SNIA_COV_DIAG_WIDTH_COLOUR,
  NC_DATA_SNIA_COV_ABSMAG_SET, 
  NC_DATA_SNIA_COV_VAR_MAG,
  NC_DATA_SNIA_COV_VAR_WIDTH,
  NC_DATA_SNIA_COV_VAR_COLOUR,
  NC_DATA_SNIA_COV_VAR_MAG_WIDTH,
  NC_DATA_SNIA_COV_VAR_MAG_COLOUR,
  NC_DATA_SNIA_COV_VAR_WIDTH_COLOUR,
  /* < private > */
  NC_DATA_SNIA_COV_TOTAL_LENGTH, /*< skip >*/
} NcDataSNIACovData;

#define NC_DATA_SNIA_COV_LENGTH NC_DATA_SNIA_COV_ABSMAG_SET

/**
 * NcDataSNIACovDataV1:
 * @NC_DATA_SNIA_COV_V1_ZCMB: Redshift in the CMB frame.
 * @NC_DATA_SNIA_COV_V1_ZHE: Redshift in sun's frame.
 * @NC_DATA_SNIA_COV_V1_SIGMA_Z: Redshift error.
 * @NC_DATA_SNIA_COV_V1_MAG: Magnitude.
 * @NC_DATA_SNIA_COV_V1_WIDTH: Width (strecth).
 * @NC_DATA_SNIA_COV_V1_COLOUR: Colour.
 * @NC_DATA_SNIA_COV_V1_THIRDPAR: Third parameter.
 * @NC_DATA_SNIA_COV_V1_SIGMA_THIRDPAR: Error on third parameter.
 * @NC_DATA_SNIA_COV_V1_ABSMAG_SET: Data set index.
 * @NC_DATA_SNIA_COV_V1_MAG_MAG: Covariance mag-mag.
 * @NC_DATA_SNIA_COV_V1_MAG_WIDTH: Covariance mag-width.
 * @NC_DATA_SNIA_COV_V1_MAG_COLOUR: Covariance mag-colour.
 * @NC_DATA_SNIA_COV_V1_WIDTH_WIDTH: Covariance width-width.
 * @NC_DATA_SNIA_COV_V1_WIDTH_COLOUR: Covariance width-colour.
 * @NC_DATA_SNIA_COV_V1_COLOUR_COLOUR: Covariance colour-colour.
 * 
 * Data ordering of Version 1 (V1) data format.
 * 
 */
typedef enum _NcDataSNIACovDataV1
{
  NC_DATA_SNIA_COV_V1_ZCMB = 0,
  NC_DATA_SNIA_COV_V1_ZHE,
  NC_DATA_SNIA_COV_V1_SIGMA_Z,
  NC_DATA_SNIA_COV_V1_MAG,
  NC_DATA_SNIA_COV_V1_WIDTH,
  NC_DATA_SNIA_COV_V1_COLOUR,
  NC_DATA_SNIA_COV_V1_THIRDPAR,
  NC_DATA_SNIA_COV_V1_SIGMA_THIRDPAR,
  NC_DATA_SNIA_COV_V1_ABSMAG_SET, 
  NC_DATA_SNIA_COV_V1_MAG_MAG,
  NC_DATA_SNIA_COV_V1_MAG_WIDTH,
  NC_DATA_SNIA_COV_V1_MAG_COLOUR,
  NC_DATA_SNIA_COV_V1_WIDTH_WIDTH,
  NC_DATA_SNIA_COV_V1_WIDTH_COLOUR,
  NC_DATA_SNIA_COV_V1_COLOUR_COLOUR,
  /* < private > */
  NC_DATA_SNIA_COV_V1_TOTAL_LENGTH, /*< skip >*/
} NcDataSNIACovDataV1;

#define NC_DATA_SNIA_COV_V1_LENGTH NC_DATA_SNIA_COV_V1_ABSMAG_SET

/**
 * NcDataSNIACovDataInit:
 * @NC_DATA_SNIA_COV_INIT_ZCMB: Redshift in the CMB frame.
 * @NC_DATA_SNIA_COV_INIT_ZHE: Redshift in sun's frame.
 * @NC_DATA_SNIA_COV_INIT_SIGMA_Z: Redshift error.
 * @NC_DATA_SNIA_COV_INIT_MAG: Magnitude.
 * @NC_DATA_SNIA_COV_INIT_WIDTH: Width (strecth).
 * @NC_DATA_SNIA_COV_INIT_COLOUR: Colour.
 * @NC_DATA_SNIA_COV_INIT_THIRDPAR: Third parameter.
 * @NC_DATA_SNIA_COV_INIT_ABSMAG_SET: Data set index.
 * @NC_DATA_SNIA_COV_INIT_COV_FULL: Full covariance matrix.
 * 
 * Bitwise control of data initialization.
 * 
 */
typedef enum _NcDataSNIACovDataInit
{
  NC_DATA_SNIA_COV_INIT_ZCMB       = 1 << 0,
  NC_DATA_SNIA_COV_INIT_ZHE        = 1 << 1,
  NC_DATA_SNIA_COV_INIT_SIGMA_Z    = 1 << 2,
  NC_DATA_SNIA_COV_INIT_MAG        = 1 << 3,
  NC_DATA_SNIA_COV_INIT_WIDTH      = 1 << 4,
  NC_DATA_SNIA_COV_INIT_COLOUR     = 1 << 5,
  NC_DATA_SNIA_COV_INIT_THIRDPAR   = 1 << 6,
  NC_DATA_SNIA_COV_INIT_ABSMAG_SET = 1 << 7, 
  NC_DATA_SNIA_COV_INIT_COV_FULL   = 1 << 8,
} NcDataSNIACovDataInit;

#define NC_DATA_SNIA_COV_INIT_ALL (NC_DATA_SNIA_COV_INIT_ZCMB | \
                                   NC_DATA_SNIA_COV_INIT_ZHE | \
                                   NC_DATA_SNIA_COV_INIT_SIGMA_Z | \
                                   NC_DATA_SNIA_COV_INIT_MAG | \
                                   NC_DATA_SNIA_COV_INIT_WIDTH | \
                                   NC_DATA_SNIA_COV_INIT_COLOUR | \
                                   NC_DATA_SNIA_COV_INIT_THIRDPAR | \
                                   NC_DATA_SNIA_COV_INIT_ABSMAG_SET | \
                                   NC_DATA_SNIA_COV_INIT_COV_FULL)

/**
 * NcDataSNIACovOrder:
 * @NC_DATA_SNIA_COV_ORDER_MAG_MAG: mag-mag.
 * @NC_DATA_SNIA_COV_ORDER_MAG_WIDTH: mag-width.
 * @NC_DATA_SNIA_COV_ORDER_MAG_COLOUR: mag-colour.
 * @NC_DATA_SNIA_COV_ORDER_WIDTH_WIDTH: width-width.
 * @NC_DATA_SNIA_COV_ORDER_WIDTH_COLOUR: width-colour.
 * @NC_DATA_SNIA_COV_ORDER_COLOUR_COLOUR: colour-colour.
 * 
 * Data ordering for covariance.
 * 
 */
typedef enum _NcDataSNIACovOrder
{
  NC_DATA_SNIA_COV_ORDER_MAG_MAG = 0,
  NC_DATA_SNIA_COV_ORDER_MAG_WIDTH,
  NC_DATA_SNIA_COV_ORDER_MAG_COLOUR,
  NC_DATA_SNIA_COV_ORDER_WIDTH_WIDTH,
  NC_DATA_SNIA_COV_ORDER_WIDTH_COLOUR,
  NC_DATA_SNIA_COV_ORDER_COLOUR_COLOUR,
  /* < private > */
  NC_DATA_SNIA_COV_ORDER_LENGTH, /*< skip >*/
} NcDataSNIACovOrder;

struct _NcDataSNIACovClass
{
  /*< private >*/
  NcmDataGaussCovClass parent_class;
};

struct _NcDataSNIACov
{
  /*< private >*/
  NcmDataGaussCov parent_instance;
  guint mu_len;
  guint uppertri_len;
  gdouble mag_cut;
  NcmVector *z_cmb;
  NcmVector *z_he;
  NcmVector *mag;
  NcmVector *width;
  NcmVector *colour;
  NcmVector *thirdpar;
  NcmVector *width_true;
  NcmVector *colour_true;
  NcmVector *mag_width_colour;
  NcmVector *sigma_z;
  NcmVector *sigma_thirdpar;
  NcmVector *cov_packed;
  NcmMatrix *cov_full;
  NcmVector *cov_full_diag;
  NcmMatrix *inv_cov_mm;
  NcmMatrix *inv_cov_mm_LU;
  gboolean has_complete_cov;
  guint cov_full_state;
  gboolean has_true_wc;
  GArray *dataset;
  GArray *dataset_size;
  guint dataset_len;
  guint data_init;
  NcmModelCtrl *cosmo_resample_ctrl;
  NcmModelCtrl *dcov_resample_ctrl;
  NcmModelCtrl *dcov_cov_full_ctrl;
};

GType nc_data_snia_cov_get_type (void) G_GNUC_CONST;

NcmData *nc_data_snia_cov_new (gboolean use_norma);
NcmData *nc_data_snia_cov_new_full (gchar *filename, gboolean use_norma);

guint nc_data_snia_cov_sigma_int_len (NcDataSNIACov *snia_cov);

NcmVector *nc_data_snia_cov_peek_z_cmb (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_z_he (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_sigma_z (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_mag (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_width (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_colour (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_thirdpar (NcDataSNIACov *snia_cov);
GArray *nc_data_snia_cov_peek_abs_mag_set (NcDataSNIACov *snia_cov);
NcmMatrix *nc_data_snia_cov_peek_cov_full (NcDataSNIACov *snia_cov);

void nc_data_snia_cov_set_mag_cut (NcDataSNIACov *snia_cov, const gdouble mag_cut);
gdouble nc_data_snia_cov_get_mag_cut (NcDataSNIACov *snia_cov);

void nc_data_snia_cov_set_z_cmb (NcDataSNIACov *snia_cov, NcmVector *z_cmb);
void nc_data_snia_cov_set_z_he (NcDataSNIACov *snia_cov, NcmVector *z_he);
void nc_data_snia_cov_set_sigma_z (NcDataSNIACov *snia_cov, NcmVector *sigma_z);
void nc_data_snia_cov_set_mag (NcDataSNIACov *snia_cov, NcmVector *mag);
void nc_data_snia_cov_set_width (NcDataSNIACov *snia_cov, NcmVector *width);
void nc_data_snia_cov_set_colour (NcDataSNIACov *snia_cov, NcmVector *colour);
void nc_data_snia_cov_set_thirdpar (NcDataSNIACov *snia_cov, NcmVector *thirdpar);
void nc_data_snia_cov_set_abs_mag_set (NcDataSNIACov *snia_cov, GArray *abs_mag_set);
void nc_data_snia_cov_set_cov_full (NcDataSNIACov *snia_cov, NcmMatrix *cov_full);

void nc_data_snia_cov_load_txt (NcDataSNIACov *snia_cov, const gchar *filename);
#ifdef NUMCOSMO_HAVE_CFITSIO
void nc_data_snia_cov_load (NcDataSNIACov *snia_cov, const gchar *filename);
void nc_data_snia_cov_save (NcDataSNIACov *snia_cov, const gchar *filename, gboolean overwrite);
#endif /* NUMCOSMO_HAVE_CFITSIO */

gdouble nc_data_snia_cov_estimate_width_colour (NcDataSNIACov *snia_cov, NcmMSet *mset);
NcmVector *nc_data_snia_cov_get_estimated_mag (NcDataSNIACov *snia_cov, NcmMSet *mset);
NcmVector *nc_data_snia_cov_get_estimated_width (NcDataSNIACov *snia_cov, NcmMSet *mset);
NcmVector *nc_data_snia_cov_get_estimated_colour (NcDataSNIACov *snia_cov, NcmMSet *mset);

void nc_data_snia_cov_load_cat (NcDataSNIACov *snia_cov, NcDataSNIAId id);
gchar *nc_data_snia_cov_get_fits (const gchar *filename, gboolean check_size);
gchar *nc_data_snia_cov_get_catalog (gchar *id);
gchar *nc_data_snia_cov_get_catalog_by_id (NcDataSNIAId id);

#define NC_DATA_SNIA_COV_SYMM_TOL (1.0e-13)

#define NC_DATA_SNIA_COV_CAT_LAST_VERSION 1

#define NC_DATA_SNIA_COV_CAT_DESC "DESC"
#define NC_DATA_SNIA_COV_DATA_DESC "Description"
#define NC_DATA_SNIA_COV_CAT_DESC_COMMENT "Catalog data description"

#define NC_DATA_SNIA_COV_CAT_VERSION "VERSION"
#define NC_DATA_SNIA_COV_CAT_VERSION_COMMENT "Version number"

#define NC_DATA_SNIA_COV_DATA_GROUP "Supernovae Ia Data"
#define NC_DATA_SNIA_COV_DATA_LEN_KEY "data-length"
#define NC_DATA_SNIA_COV_DATA_KEY "snia-data"

#define NC_DATA_SNIA_COV_DATA_HAS_COMPLETE_COV_KEY "has-complete-cov"
#define NC_DATA_SNIA_COV_CAT_HAS_COMPLETE_COV "CMPL_COV"
#define NC_DATA_SNIA_COV_CAT_HAS_COMPLETE_COV_COMMENT "Whether the covariance matrix is complete"

#define NC_DATA_SNIA_COV_MAG_CUT "MAG_CUT"
#define NC_DATA_SNIA_COV_MAG_CUT_DEFAULT (10.0)
#define NC_DATA_SNIA_COV_MAG_CUT_COMMENT "Absolute magnitude cut"

#define NC_DATA_SNIA_COV_MAG_KEY "magnitude"
#define NC_DATA_SNIA_COV_WIDTH_KEY "width"
#define NC_DATA_SNIA_COV_COLOUR_KEY "colour"

#define NC_DATA_SNIA_COV_MAG_WIDTH_KEY "magnitude-width"
#define NC_DATA_SNIA_COV_MAG_COLOUR_KEY "magnitude-colour"
#define NC_DATA_SNIA_COV_WIDTH_COLOUR_KEY "width-colour"

G_END_DECLS

#endif /* _NC_DATA_SNIA_COV_H_ */

