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
 * FIXME
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
  NC_DATA_SNIA_COV_VAR_WIDTH_COLOUR, /*< private >*/
  NC_DATA_SNIA_COV_TOTAL_LENGTH,     /*< skip >*/
} NcDataSNIACovData;

#define NC_DATA_SNIA_COV_LENGTH NC_DATA_SNIA_COV_ABSMAG_SET

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
  guint dataset_len;
  gdouble sigma_pecz;
  gchar *filename;
};

GType nc_data_snia_cov_get_type (void) G_GNUC_CONST;

NcmData *nc_data_snia_cov_new (gboolean use_det);
NcmData *nc_data_snia_cov_new_full (gchar *filename, gboolean use_det);

void nc_data_snia_cov_set_size (NcDataSNIACov *snia_cov, guint mu_len);

void nc_data_snia_cov_load_txt (NcDataSNIACov *snia_cov, const gchar *filename);

#ifdef NUMCOSMO_HAVE_CFITSIO
void nc_data_snia_cov_load (NcDataSNIACov *snia_cov, const gchar *filename);
void nc_data_snia_cov_save (NcDataSNIACov *snia_cov, const gchar *filename, gboolean overwrite);
#endif /* NUMCOSMO_HAVE_CFITSIO */

#define NC_DATA_SNIA_COV_DATA_GROUP "Supernovae Ia Data"
#define NC_DATA_SNIA_COV_DATA_LEN_KEY "data-length"
#define NC_DATA_SNIA_COV_DATA_KEY "snia-data"

#define NC_DATA_SNIA_COV_MAG_KEY "magnitude"
#define NC_DATA_SNIA_COV_WIDTH_KEY "width"
#define NC_DATA_SNIA_COV_COLOUR_KEY "colour"

#define NC_DATA_SNIA_COV_MAG_WIDTH_KEY "magnitude-width"
#define NC_DATA_SNIA_COV_MAG_COLOUR_KEY "magnitude-colour"
#define NC_DATA_SNIA_COV_WIDTH_COLOUR_KEY "width-colour"

G_END_DECLS

#endif /* _NC_DATA_SNIA_COV_H_ */

