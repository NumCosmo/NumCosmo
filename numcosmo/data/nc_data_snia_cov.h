/***************************************************************************
 *            nc_data_snia_cov.h
 *
 *  Sat December 08 15:58:20 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#define NC_DATA_SNIA_COV_ERROR (nc_data_snia_cov_error_quark())

typedef struct _NcDataSNIACovClass NcDataSNIACovClass;
typedef struct _NcDataSNIACov NcDataSNIACov;
typedef struct _NcDataSNIACovPrivate NcDataSNIACovPrivate;

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
  NcDataSNIACovPrivate *priv;
};

/**
 * NcDataSNIACovError:
 * @NC_DATA_SNIA_COV_ERROR_ID_NOT_FOUND: id not found
 * @NC_DATA_SNIA_COV_ERROR_INVALID_ID: invalid id
 * @NC_DATA_SNIA_COV_ERROR_INVALID_SAMPLE: invalid sample
 *
 * #NcDataSNIACov error messages.
 *
 */
typedef enum /*< enum,underscore_name=NC_DATA_SNIA_COV_ERROR >*/
{
  NC_DATA_SNIA_COV_ERROR_ID_NOT_FOUND,
  NC_DATA_SNIA_COV_ERROR_INVALID_ID,
  NC_DATA_SNIA_COV_ERROR_INVALID_SAMPLE,
  /* < private > */
  NC_DATA_SNIA_COV_ERROR_LENGTH, /*< skip >*/
} NcDataSNIACovError;

GType nc_data_snia_cov_get_type (void) G_GNUC_CONST;
GQuark nc_data_snia_cov_error_quark (void) G_GNUC_CONST;

NcDataSNIACov *nc_data_snia_cov_new (gboolean use_norma, guint cat_version);
NcDataSNIACov *nc_data_snia_cov_new_full (const gchar *filename, gboolean use_norma);
NcDataSNIACov *nc_data_snia_cov_new_from_cat_id (NcDataSNIAId id, gboolean use_norma);

guint nc_data_snia_cov_sigma_int_len (NcDataSNIACov *snia_cov);
guint nc_data_snia_cov_snia_len (NcDataSNIACov *snia_cov);

NcmVector *nc_data_snia_cov_peek_z_hd (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_z_cmb (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_z_he (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_sigma_z (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_mag (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_width (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_colour (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_thirdpar (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_ceph_dist (NcDataSNIACov *snia_cov);
GArray *nc_data_snia_cov_peek_abs_mag_set (NcDataSNIACov *snia_cov);
NcmMatrix *nc_data_snia_cov_peek_cov_full (NcDataSNIACov *snia_cov);
NcmVector *nc_data_snia_cov_peek_cov_packed (NcDataSNIACov *snia_cov);
NcmMatrix *nc_data_snia_cov_peek_cov_mbc_mbc (NcDataSNIACov *snia_cov);
GArray *nc_data_snia_cov_peek_dataset (NcDataSNIACov *snia_cov);
GArray *nc_data_snia_cov_peek_is_calib (NcDataSNIACov *snia_cov);
GArray *nc_data_snia_cov_peek_used_in_sh0es (NcDataSNIACov *snia_cov);

void nc_data_snia_cov_set_mag_cut (NcDataSNIACov *snia_cov, const gdouble mag_cut);
gdouble nc_data_snia_cov_get_mag_cut (NcDataSNIACov *snia_cov);

void nc_data_snia_cov_set_z_hd (NcDataSNIACov *snia_cov, NcmVector *z_hd);
void nc_data_snia_cov_set_z_cmb (NcDataSNIACov *snia_cov, NcmVector *z_cmb);
void nc_data_snia_cov_set_z_he (NcDataSNIACov *snia_cov, NcmVector *z_he);
void nc_data_snia_cov_set_sigma_z (NcDataSNIACov *snia_cov, NcmVector *sigma_z);
void nc_data_snia_cov_set_mag (NcDataSNIACov *snia_cov, NcmVector *mag);
void nc_data_snia_cov_set_mag_b_corr (NcDataSNIACov *snia_cov, NcmVector *mag_b_corr);
void nc_data_snia_cov_set_ceph_dist (NcDataSNIACov *snia_cov, NcmVector *ceph_dist);
void nc_data_snia_cov_set_width (NcDataSNIACov *snia_cov, NcmVector *width);
void nc_data_snia_cov_set_colour (NcDataSNIACov *snia_cov, NcmVector *colour);
void nc_data_snia_cov_set_thirdpar (NcDataSNIACov *snia_cov, NcmVector *thirdpar);
void nc_data_snia_cov_set_abs_mag_set (NcDataSNIACov *snia_cov, GArray *abs_mag_set);
void nc_data_snia_cov_set_is_calib (NcDataSNIACov *snia_cov, GArray *is_calib);
void nc_data_snia_cov_set_used_in_sh0es (NcDataSNIACov *snia_cov, GArray *used_in_sh0es);
void nc_data_snia_cov_set_cov_full (NcDataSNIACov *snia_cov, NcmMatrix *cov_full);
void nc_data_snia_cov_set_cov_mbc_mbc (NcDataSNIACov *snia_cov, NcmMatrix *cov_mbc_mbc);

void nc_data_snia_cov_load_txt (NcDataSNIACov *snia_cov, const gchar *filename);
void nc_data_snia_cov_save (NcDataSNIACov *snia_cov, const gchar *filename, gboolean overwrite);

gdouble nc_data_snia_cov_estimate_width_colour (NcDataSNIACov *snia_cov, NcmMSet *mset);
NcmVector *nc_data_snia_cov_get_estimated_mag (NcDataSNIACov *snia_cov, NcmMSet *mset);
NcmVector *nc_data_snia_cov_get_estimated_width (NcDataSNIACov *snia_cov, NcmMSet *mset);
NcmVector *nc_data_snia_cov_get_estimated_colour (NcDataSNIACov *snia_cov, NcmMSet *mset);

gchar *nc_data_snia_cov_get_fits (const gchar *filename, gboolean check_size);
NcDataSNIAId nc_data_snia_cov_get_catalog_id (gchar *id, GError **err);
gchar *nc_data_snia_cov_get_catalog (gchar *id, GError **err);
gchar *nc_data_snia_cov_get_catalog_by_id (NcDataSNIAId id);

NcDataSNIACov *nc_data_snia_cov_apply_filter_sh0es_z (NcDataSNIACov *snia_cov, const gdouble z_min, const gboolean use_calib, GError **err);

G_END_DECLS

#endif /* _NC_DATA_SNIA_COV_H_ */

