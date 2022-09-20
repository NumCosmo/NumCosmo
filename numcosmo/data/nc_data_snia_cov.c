/***************************************************************************
 *            nc_data_snia_cov.c
 *
 *  Sat December 08 15:58:15 2012
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

/**
 * SECTION:nc_data_snia_cov
 * @title: NcDataSNIACov
 * @short_description: Type Ia supernovae data with covariance error matrix
 *
 * See #NcSNIADistCov.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "data/nc_data_snia_cov.h"
#include "nc_snia_dist_cov.h"
#include "nc_enum_types.h"
#include "math/ncm_model_ctrl.h"
#include "math/ncm_lapack.h"
#include "math/ncm_iset.h"
#include "math/ncm_cfg.h"

#include <gio/gio.h>
#include <glib/gstdio.h>

#ifndef NUMCOSMO_GIR_SCAN
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#endif /* NUMCOSMO_GIR_SCAN */

/*
 * Data ordering of Version 0 (V0) data format.
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

/*
 * Data ordering of Version 1 (V1) data format.
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

#define NC_DATA_SNIA_COV_V1_LENGTH (NC_DATA_SNIA_COV_V1_ABSMAG_SET)
#define NC_DATA_SNIA_COV_V1_COV_FULL (NC_DATA_SNIA_COV_V1_TOTAL_LENGTH)

/**
 * NcDataSNIACovDataV2:
 * @NC_DATA_SNIA_COV_V2_ZHD: Hubble Diagram Redshift (with CMB and VPEC corrections)
 * @NC_DATA_SNIA_COV_V2_SIGMA_ZHD: Hubble Diagram Redshift Uncertainty
 * @NC_DATA_SNIA_COV_V2_ZCMB: CMB Corrected Redshift
 * @NC_DATA_SNIA_COV_V2_SIGMA_ZCMB: CMB Corrected Redshift Uncertainty
 * @NC_DATA_SNIA_COV_V2_ZHE: Heliocentric Redshift
 * @NC_DATA_SNIA_COV_V2_SIGMA_ZHE: Heliocentric Redshift Uncertainty
 * @NC_DATA_SNIA_COV_V2_MAG_B_CORR: Tripp1998 corrected/standardized m_b magnitude
 * @NC_DATA_SNIA_COV_V2_SIGMA_MAG_B_CORR: Tripp1998 corrected/standardized m_b magnitude uncertainty as determined from the diagonal of the covariance matrix. **WARNING, DO NOT FIT COSMOLOGICAL PARAMETERS WITH THESE UNCERTAINTIES. YOU MUST USE THE FULL COVARIANCE. THIS IS ONLY FOR PLOTTING/VISUAL PURPOSES**
 * @NC_DATA_SNIA_COV_V2_MU_CORR: Tripp1998 corrected/standardized distance modulus where fiducial SNIa magnitude (M) has been determined from SH0ES 2021 Cepheid host distances.
 * @NC_DATA_SNIA_COV_V2_SIGMA_MU_CORR: Uncertainty on MU_SH0ES as determined from the diagonal of the covariance matrix. **WARNING, DO NOT FIT COSMOLOGICAL PARAMETERS WITH THESE UNCERTAINTIES. YOU MUST USE THE FULL COVARIANCE. THIS IS ONLY FOR PLOTTING/VISUAL PURPOSES**
 * @NC_DATA_SNIA_COV_V2_CEPH_DIST: Cepheid calculated absolute distance to host (uncertainty is incorporated in covariance matrix .cov)
 * @NC_DATA_SNIA_COV_V2_IS_CALIB: binary to designate if this SN is in a host that has an associated cepheid distance
 * @NC_DATA_SNIA_COV_V2_USED_IN_SH0ES: 1 if used in SH0ES 2021 Hubble Flow dataset. 0 if not included.
 * @NC_DATA_SNIA_COV_V2_COLOUR: SALT2 color
 * @NC_DATA_SNIA_COV_V2_SIGMA_COLOUR: SALT2 color uncertainty
 * @NC_DATA_SNIA_COV_V2_WIDTH: SALT2 stretch
 * @NC_DATA_SNIA_COV_V2_SIGMA_WIDTH: SALT2 stretch uncertainty
 * @NC_DATA_SNIA_COV_V2_MAG: SALT2 uncorrected brightness
 * @NC_DATA_SNIA_COV_V2_SIGMA_MAG: SALT2 uncorrected brightness uncertainty
 * @NC_DATA_SNIA_COV_V2_AMPL: SALT2 light curve amplitude
 * @NC_DATA_SNIA_COV_V2_SIGMA_AMPL: SALT2 light curve amplitude uncertainty
 * @NC_DATA_SNIA_COV_V2_WIDTH_COLOUR: SALT2 fit covariance between x1 and c
 * @NC_DATA_SNIA_COV_V2_WIDTH_AMPL: SALT2 fit covariance between x1 and x0
 * @NC_DATA_SNIA_COV_V2_COLOUR_AMPL: SALT2 fit covariance between c and x0
 * @NC_DATA_SNIA_COV_V2_RA: Right Ascension
 * @NC_DATA_SNIA_COV_V2_DEC: Declination
 * @NC_DATA_SNIA_COV_V2_HOST_RA: Host Galaxy RA
 * @NC_DATA_SNIA_COV_V2_HOST_DEC: Host Galaxy DEC
 * @NC_DATA_SNIA_COV_V2_HOST_ANG_SEP: Angular separation between SN and host (arcsec)
 * @NC_DATA_SNIA_COV_V2_VPEC: Peculiar velocity (km/s)
 * @NC_DATA_SNIA_COV_V2_SIGMA_VPEC: Peculiar velocity uncertainty (km/s)
 * @NC_DATA_SNIA_COV_V2_MAG_B_CORR_MAG_B_CORR: Full covariance MAG_B_CORR
 *
 * Data ordering of Version 2 (V2) data format.
 *
 */
typedef enum _NcDataSNIACovDataV2
{
  NC_DATA_SNIA_COV_V2_ZHD = 0,
  NC_DATA_SNIA_COV_V2_SIGMA_ZHD,
  NC_DATA_SNIA_COV_V2_ZCMB,
  NC_DATA_SNIA_COV_V2_SIGMA_ZCMB,
  NC_DATA_SNIA_COV_V2_ZHE,
  NC_DATA_SNIA_COV_V2_SIGMA_ZHE,
  NC_DATA_SNIA_COV_V2_MAG_B_CORR,
  NC_DATA_SNIA_COV_V2_SIGMA_MAG_B_CORR,
  NC_DATA_SNIA_COV_V2_MU_CORR,
  NC_DATA_SNIA_COV_V2_SIGMA_MU_CORR,
  NC_DATA_SNIA_COV_V2_CEPH_DIST,
  NC_DATA_SNIA_COV_V2_IS_CALIB,
  NC_DATA_SNIA_COV_V2_USED_IN_SH0ES,
  NC_DATA_SNIA_COV_V2_COLOUR,
  NC_DATA_SNIA_COV_V2_SIGMA_COLOUR,
  NC_DATA_SNIA_COV_V2_WIDTH,
  NC_DATA_SNIA_COV_V2_SIGMA_WIDTH,
  NC_DATA_SNIA_COV_V2_MAG,
  NC_DATA_SNIA_COV_V2_SIGMA_MAG,
  NC_DATA_SNIA_COV_V2_AMPL,
  NC_DATA_SNIA_COV_V2_SIGMA_AMPL,
  NC_DATA_SNIA_COV_V2_WIDTH_COLOUR,
  NC_DATA_SNIA_COV_V2_WIDTH_AMPL,
  NC_DATA_SNIA_COV_V2_COLOUR_AMPL,
  NC_DATA_SNIA_COV_V2_RA,
  NC_DATA_SNIA_COV_V2_DEC,
  NC_DATA_SNIA_COV_V2_HOST_RA,
  NC_DATA_SNIA_COV_V2_HOST_DEC,
  NC_DATA_SNIA_COV_V2_HOST_ANG_SEP,
  NC_DATA_SNIA_COV_V2_VPEC,
  NC_DATA_SNIA_COV_V2_SIGMA_VPEC,
  NC_DATA_SNIA_COV_V2_MAG_B_CORR_MAG_B_CORR,
  /* < private > */
  NC_DATA_SNIA_COV_V2_TOTAL_LENGTH, /*< skip >*/
} NcDataSNIACovDataV2;

#define NC_DATA_SNIA_COV_V2_LENGTH (NC_DATA_SNIA_COV_V2_SIGMA_VPEC + 1)

#define NC_DATA_SNIA_COV_V1_INIT(D) (1UL << NC_DATA_SNIA_COV_V1_##D)
#define NC_DATA_SNIA_COV_V2_INIT(D) (1UL << NC_DATA_SNIA_COV_V2_##D)

#define NC_DATA_SNIA_COV_V1_ALL \
  (NC_DATA_SNIA_COV_V1_INIT (ZCMB) | \
   NC_DATA_SNIA_COV_V1_INIT (ZHE) | \
   NC_DATA_SNIA_COV_V1_INIT (SIGMA_Z) | \
   NC_DATA_SNIA_COV_V1_INIT (MAG) | \
   NC_DATA_SNIA_COV_V1_INIT (WIDTH) | \
   NC_DATA_SNIA_COV_V1_INIT (COLOUR) | \
   NC_DATA_SNIA_COV_V1_INIT (THIRDPAR) | \
   NC_DATA_SNIA_COV_V1_INIT (ABSMAG_SET) | \
   NC_DATA_SNIA_COV_V1_INIT (COV_FULL))

#define NC_DATA_SNIA_COV_V2_ALL ((1UL << NC_DATA_SNIA_COV_V2_TOTAL_LENGTH) - 1UL)

#define NC_DATA_SNIA_COV_SYMM_TOL (1.0e-13)

#define NC_DATA_SNIA_COV_CAT_LAST_VERSION 2

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

#define NC_DATA_SNIA_COV_MBC_MBC_KEY "mbc-mbc"

struct _NcDataSNIACovPrivate
{
  guint mu_len;
  guint uppertri_len;
  gdouble mag_cut;
  NcmVector *z_hd;
  NcmVector *sigma_z_hd;
  NcmVector *z_cmb;
  NcmVector *sigma_z_cmb;
  NcmVector *z_he;
  NcmVector *sigma_z_he;
  NcmVector *mag_b_corr;
  NcmVector *sigma_mag_b_corr;
  NcmVector *mu_corr;
  NcmVector *sigma_mu_corr;
  NcmVector *ceph_dist;
  GArray *is_calib;
  GArray *used_in_sh0es;
  NcmVector *mag;
  NcmVector *sigma_mag;
  NcmVector *width;
  NcmVector *sigma_width;
  NcmVector *colour;
  NcmVector *sigma_colour;
  NcmVector *thirdpar;
  NcmVector *ampl;
  NcmVector *sigma_ampl;
  NcmVector *width_true;
  NcmVector *colour_true;
  NcmVector *mag_width_colour;
  NcmVector *sigma_z;
  NcmVector *sigma_thirdpar;
  NcmVector *cov_width_colour;
  NcmVector *cov_width_ampl;
  NcmVector *cov_colour_ampl;
  NcmVector *ra;
  NcmVector *dec;
  NcmVector *host_ra;
  NcmVector *host_dec;
  NcmVector *host_ang_sep;
  NcmVector *vpec;
  NcmVector *sigma_vpec;
  NcmMatrix *cov_mbc_mbc;
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
  guint64 data_init;
  guint cat_version;
  NcmModelCtrl *cosmo_resample_ctrl;
  NcmModelCtrl *dcov_resample_ctrl;
  NcmModelCtrl *dcov_cov_full_ctrl;
};

enum
{
  PROP_0,
  PROP_CAT_VERSION,
  PROP_MAG_CUT,
  PROP_ZHD,
  PROP_ZCMB,
  PROP_ZHE,
  PROP_SIGMA_Z,
  PROP_MAG,
  PROP_MAG_B_CORR,
  PROP_CEPH_DIST,
  PROP_WIDTH,
  PROP_COLOUR,
  PROP_THIRDPAR,
  PROP_SIGMA_THIRDPAR,
  PROP_ABSMAG_SET,
  PROP_IS_CALIB,
  PROP_USED_IN_SH0ES,
  PROP_COV_FULL,
  PROP_HAS_COMPLETE_COV,
  PROP_COV_MBC_MBC,
  PROP_SIZE,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcDataSNIACov, nc_data_snia_cov, NCM_TYPE_DATA_GAUSS_COV);
G_DEFINE_QUARK(nc-data-snia-cov-error-quark, nc_data_snia_cov_error);

enum
{
  NC_DATA_SNIA_COV_PREP_TO_NOTHING = 0,
  NC_DATA_SNIA_COV_PREP_TO_RESAMPLE,
  NC_DATA_SNIA_COV_PREP_TO_ESTIMATE,
};

static void
nc_data_snia_cov_init (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv = nc_data_snia_cov_get_instance_private (snia_cov);

  self->mu_len       = 0;
  self->uppertri_len = 0;

  self->mag_cut = 0.0;

  self->z_hd  = NULL;
  self->z_cmb = NULL;
  self->z_he  = NULL;

  self->mag_b_corr = NULL;
  self->mu_corr    = NULL;
  self->ceph_dist  = NULL;
  self->mag        = NULL;
  self->width      = NULL;
  self->colour     = NULL;
  self->ampl       = NULL;
  self->vpec       = NULL;
  self->thirdpar   = NULL;

  self->ra           = NULL;
  self->dec          = NULL;
  self->host_ra      = NULL;
  self->host_dec     = NULL;
  self->host_ang_sep = NULL;

  self->width_true       = NULL;
  self->colour_true      = NULL;
  self->mag_width_colour = NULL;

  self->sigma_z_hd  = NULL;
  self->sigma_z_cmb = NULL;
  self->sigma_z_he  = NULL;

  self->sigma_z        = NULL;
  self->sigma_thirdpar = NULL;

  self->sigma_mag_b_corr = NULL;
  self->sigma_mu_corr    = NULL;
  self->sigma_colour     = NULL;
  self->sigma_width      = NULL;
  self->sigma_mag        = NULL;
  self->sigma_ampl       = NULL;
  self->sigma_vpec       = NULL;

  self->cov_width_colour = NULL;
  self->cov_width_ampl   = NULL;
  self->cov_colour_ampl  = NULL;

  self->cov_mbc_mbc   = NULL;
  self->cov_full      = NULL;
  self->cov_full_diag = NULL;
  self->cov_packed    = NULL;

  self->inv_cov_mm    = NULL;
  self->inv_cov_mm_LU = NULL;

  self->has_complete_cov = FALSE;
  self->cov_full_state   = NC_DATA_SNIA_COV_PREP_TO_NOTHING;
  self->has_true_wc      = FALSE;

  self->dataset      = NULL;
  self->dataset_size = NULL;
  self->dataset_len  = 0;

  self->is_calib      = NULL;
  self->used_in_sh0es = NULL;

  self->data_init = 0;

  self->cosmo_resample_ctrl = ncm_model_ctrl_new (NULL);
  self->dcov_resample_ctrl  = ncm_model_ctrl_new (NULL);
  self->dcov_cov_full_ctrl  = ncm_model_ctrl_new (NULL);
}

static void
nc_data_snia_cov_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (object);
  NcDataSNIACovPrivate * const self = snia_cov->priv;

  g_return_if_fail (NC_IS_DATA_SNIA_COV (object));

  switch (prop_id)
  {
    case PROP_CAT_VERSION:
      self->cat_version = g_value_get_uint (value);
      break;
    case PROP_MAG_CUT:
      nc_data_snia_cov_set_mag_cut (snia_cov, g_value_get_double (value));
      break;
    case PROP_ZHD:
      nc_data_snia_cov_set_z_hd (snia_cov, g_value_get_object (value));
      break;
    case PROP_ZCMB:
      nc_data_snia_cov_set_z_cmb (snia_cov, g_value_get_object (value));
      break;
    case PROP_ZHE:
      nc_data_snia_cov_set_z_he (snia_cov, g_value_get_object (value));
      break;
    case PROP_SIGMA_Z:
      nc_data_snia_cov_set_sigma_z (snia_cov, g_value_get_object (value));
      break;
    case PROP_MAG:
      nc_data_snia_cov_set_mag (snia_cov, g_value_get_object (value));
      break;
    case PROP_MAG_B_CORR:
      nc_data_snia_cov_set_mag_b_corr (snia_cov, g_value_get_object (value));
      break;
    case PROP_CEPH_DIST:
      nc_data_snia_cov_set_ceph_dist (snia_cov, g_value_get_object (value));
      break;
    case PROP_WIDTH:
      nc_data_snia_cov_set_width (snia_cov, g_value_get_object (value));
      break;
    case PROP_COLOUR:
      nc_data_snia_cov_set_colour (snia_cov, g_value_get_object (value));
      break;
    case PROP_THIRDPAR:
      nc_data_snia_cov_set_thirdpar (snia_cov, g_value_get_object (value));
      break;
    case PROP_SIGMA_THIRDPAR:
      ncm_vector_substitute (&self->sigma_thirdpar, g_value_get_object (value), TRUE);
      break;
    case PROP_ABSMAG_SET:
    {
      GVariant *var = g_value_get_variant (value);

      if (var != NULL)
      {
        ncm_cfg_array_set_variant (self->dataset, var);
        nc_data_snia_cov_set_abs_mag_set (snia_cov, self->dataset);
      }

      break;
    }
    case PROP_IS_CALIB:
    {
      GVariant *var = g_value_get_variant (value);

      if (var != NULL)
      {
        ncm_cfg_array_set_variant (self->is_calib, var);
        nc_data_snia_cov_set_is_calib (snia_cov, self->is_calib);
      }

      break;
    }
    case PROP_USED_IN_SH0ES:
    {
      GVariant *var = g_value_get_variant (value);

      if (var != NULL)
      {
        ncm_cfg_array_set_variant (self->used_in_sh0es, var);
        nc_data_snia_cov_set_used_in_sh0es (snia_cov, self->used_in_sh0es);
      }

      break;
    }
    case PROP_COV_FULL:
      nc_data_snia_cov_set_cov_full (snia_cov, g_value_get_object (value));
      break;
    case PROP_HAS_COMPLETE_COV:
      self->has_complete_cov = g_value_get_boolean (value);
      break;
    case PROP_COV_MBC_MBC:
      nc_data_snia_cov_set_cov_mbc_mbc (snia_cov, g_value_get_object (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_snia_cov_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (object);
  NcDataSNIACovPrivate * const self = snia_cov->priv;

  g_return_if_fail (NC_IS_DATA_SNIA_COV (object));

  switch (prop_id)
  {
    case PROP_CAT_VERSION:
      g_value_set_uint (value, self->cat_version);
      break;
    case PROP_MAG_CUT:
      g_value_set_double (value, nc_data_snia_cov_get_mag_cut (snia_cov));
      break;
    case PROP_ZHD:
      g_value_set_object (value, self->z_hd);
      break;
    case PROP_ZCMB:
      g_value_set_object (value, self->z_cmb);
      break;
    case PROP_ZHE:
      g_value_set_object (value, self->z_he);
      break;
    case PROP_SIGMA_Z:
      g_value_set_object (value, self->sigma_z);
      break;
    case PROP_MAG:
      g_value_set_object (value, self->mag);
      break;
    case PROP_MAG_B_CORR:
      g_value_set_object (value, self->mag_b_corr);
      break;
    case PROP_CEPH_DIST:
      g_value_set_object (value, self->ceph_dist);
      break;
    case PROP_WIDTH:
      g_value_set_object (value, self->width);
      break;
    case PROP_COLOUR:
      g_value_set_object (value, self->colour);
      break;
    case PROP_THIRDPAR:
      g_value_set_object (value, self->thirdpar);
      break;
    case PROP_SIGMA_THIRDPAR:
      g_value_set_object (value, self->sigma_thirdpar);
      break;
    case PROP_ABSMAG_SET:
      g_value_take_variant (value, ncm_cfg_array_to_variant (self->dataset, G_VARIANT_TYPE ("u")));
      break;
    case PROP_IS_CALIB:
      g_value_take_variant (value, ncm_cfg_array_to_variant (self->is_calib, G_VARIANT_TYPE ("u")));
      break;
    case PROP_USED_IN_SH0ES:
      g_value_take_variant (value, ncm_cfg_array_to_variant (self->used_in_sh0es, G_VARIANT_TYPE ("u")));
      break;
    case PROP_COV_FULL:
      g_value_set_object (value, nc_data_snia_cov_peek_cov_full (snia_cov));
      break;
    case PROP_HAS_COMPLETE_COV:
      g_value_set_boolean (value, self->has_complete_cov);
      break;
    case PROP_COV_MBC_MBC:
      g_value_set_object (value, nc_data_snia_cov_peek_cov_mbc_mbc (snia_cov));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void _nc_data_snia_cov_prep_to_resample (NcDataSNIACov *snia_cov, NcSNIADistCov *dcov);
static void _nc_data_snia_cov_prep_to_estimate (NcDataSNIACov *snia_cov, NcSNIADistCov *dcov);
static void _nc_data_snia_cov_restore_full_cov (NcDataSNIACov *snia_cov);

static void
nc_data_snia_cov_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_snia_cov_parent_class)->constructed (object);
  {
  }
}

static void
nc_data_snia_cov_dispose (GObject *object)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (object);
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  NcmDataGaussCov *gauss  = NCM_DATA_GAUSS_COV (object);

  ncm_data_gauss_cov_set_size (gauss, 0);

  ncm_model_ctrl_clear (&self->cosmo_resample_ctrl);
  ncm_model_ctrl_clear (&self->dcov_resample_ctrl);
  ncm_model_ctrl_clear (&self->dcov_cov_full_ctrl);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_snia_cov_parent_class)->dispose (object);
}

static void
nc_data_snia_cov_finalize (GObject *object)
{
  /* NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (object); */

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_snia_cov_parent_class)->finalize (object);
}

static void _nc_data_snia_cov_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_snia_cov_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng);
static void _nc_data_snia_cov_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);
static gboolean _nc_data_snia_cov_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *cov);
static void _nc_data_snia_cov_set_size (NcmDataGaussCov *gauss, guint np);

static void
nc_data_snia_cov_class_init (NcDataSNIACovClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class          = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass *gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &nc_data_snia_cov_set_property;
  object_class->get_property = &nc_data_snia_cov_get_property;
  object_class->constructed  = &nc_data_snia_cov_constructed;
  object_class->dispose      = &nc_data_snia_cov_dispose;
  object_class->finalize     = &nc_data_snia_cov_finalize;

  g_object_class_install_property (object_class,
                                   PROP_CAT_VERSION,
                                   g_param_spec_uint ("cat-version",
                                                      NULL,
                                                      "Catalog version",
                                                      0, NC_DATA_SNIA_COV_CAT_LAST_VERSION, NC_DATA_SNIA_COV_CAT_LAST_VERSION,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MAG_CUT,
                                   g_param_spec_double ("magnitude-cut",
                                                        NULL,
                                                        "Threshold where to change absolute magnitude",
                                                        0.0, G_MAXDOUBLE, NC_DATA_SNIA_COV_MAG_CUT_DEFAULT,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ZHD,
                                   g_param_spec_object ("z-hd",
                                                        NULL,
                                                        "Data CMB redshifts (peculiar velocity corrected)",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ZCMB,
                                   g_param_spec_object ("z-cmb",
                                                        NULL,
                                                        "Data cmb redshifts",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ZHE,
                                   g_param_spec_object ("z-He",
                                                        NULL,
                                                        "Data He redshifts",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SIGMA_Z,
                                   g_param_spec_object ("sigma-z",
                                                        NULL,
                                                        "Redshifts standard deviation",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MAG,
                                   g_param_spec_object ("magnitudes",
                                                        NULL,
                                                        "Magnitudes",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MAG_B_CORR,
                                   g_param_spec_object ("magnitude-b-corrected",
                                                        NULL,
                                                        "Magnitude B corrected",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_CEPH_DIST,
                                   g_param_spec_object ("ceph-dist",
                                                        NULL,
                                                        "Cepheid distance",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_WIDTH,
                                   g_param_spec_object ("width",
                                                        NULL,
                                                        "Width",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_COLOUR,
                                   g_param_spec_object ("colour",
                                                        NULL,
                                                        "Colour",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_THIRDPAR,
                                   g_param_spec_object ("thirdpar",
                                                        NULL,
                                                        "Thirdpar",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SIGMA_THIRDPAR,
                                   g_param_spec_object ("sigma-thirdpar",
                                                        NULL,
                                                        "Thirdpar standard deviation",
                                                        NCM_TYPE_VECTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ABSMAG_SET,
                                   g_param_spec_variant ("absmag-set",
                                                         NULL,
                                                         "Absolute magnitude set",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_IS_CALIB,
                                   g_param_spec_variant ("is-calib",
                                                         NULL,
                                                         "Whether the SNIa is a calibrator",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_USED_IN_SH0ES,
                                   g_param_spec_variant ("used-in-sh0es",
                                                         NULL,
                                                         "Whether the SNIa was used in SH0ES",
                                                         G_VARIANT_TYPE ("au"), NULL,
                                                         G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_COV_FULL,
                                   g_param_spec_object ("cov-full",
                                                        NULL,
                                                        "Full covariance matrix",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_HAS_COMPLETE_COV,
                                   g_param_spec_boolean ("has-complete-cov",
                                                         NULL,
                                                         "Whether the covariance matrix is complete",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_COV_MBC_MBC,
                                   g_param_spec_object ("cov-mbc-mbc",
                                                        NULL,
                                                        "Covariance matrix for mag b corr",
                                                        NCM_TYPE_MATRIX,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->resample = &_nc_data_snia_cov_resample;
  data_class->prepare  = &_nc_data_snia_cov_prepare;

  gauss_class->mean_func = &_nc_data_snia_cov_mean_func;
  gauss_class->cov_func  = &_nc_data_snia_cov_func;
  gauss_class->set_size  = &_nc_data_snia_cov_set_size;

}

static void
_nc_data_snia_cov_prepare (NcmData *data, NcmMSet *mset)
{
  NcSNIADistCov *dcov = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));

  g_assert (dcov != NULL);

  nc_snia_dist_cov_prepare_if_needed (dcov, mset);
}

static void
_nc_data_snia_cov_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (gauss);
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  NcHICosmo *cosmo        = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcSNIADistCov *dcov     = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));

  switch (self->cat_version)
  {
    case 0:
    case 1:
      nc_snia_dist_cov_mean (dcov, cosmo, snia_cov, vp);
      break;
    case 2:
      nc_snia_dist_cov_mean_V2 (dcov, cosmo, snia_cov, vp);
      break;
    default:
      g_assert_not_reached ();
      break;
  }
}

static gboolean
_nc_data_snia_cov_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *cov)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (gauss);
  NcDataSNIACovPrivate * const self = snia_cov->priv;

  switch (self->cat_version)
  {
    case 0:
    case 1:
    {
      NcSNIADistCov *dcov = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
      return nc_snia_dist_cov_calc (dcov, snia_cov, cov);
    }
    case 2:
      ncm_matrix_memcpy (cov, self->cov_mbc_mbc);
      return FALSE;
      break;
    default:
      g_assert_not_reached ();
      return FALSE;
      break;
  }
}


/**
 * nc_data_snia_cov_new:
 * @use_norma: Whether to use the correct Likelihood normalization
 * @cat_version: Catalog version
 *
 * Creates a new empty #NcDataSNIACov object. If @use_norma is
 * true the object will use the correct Likelihood normalization
 * when calculating $-2\ln(L)$.
 *
 * @cat_version indicates the data format used by the likelihood.
 *
 * Returns: (transfer full): the newly created instance of #NcDataSNIACov.
 */
NcDataSNIACov *
nc_data_snia_cov_new (gboolean use_norma, guint cat_version)
{
  return g_object_new (NC_TYPE_DATA_SNIA_COV,
                       "use-norma", use_norma,
                       "cat-version", cat_version,
                       NULL);
}

#ifdef NUMCOSMO_HAVE_CFITSIO
static void nc_data_snia_cov_load_V0 (NcDataSNIACov *snia_cov, fitsfile *fptr);
#endif /* NUMCOSMO_HAVE_CFITSIO */

static void _nc_data_snia_cov_load_snia_data_V01 (NcDataSNIACov *snia_cov, const gchar *filename);
static void _nc_data_snia_cov_load_snia_data_V2 (NcDataSNIACov *snia_cov, const gchar *filename);
static void _nc_data_snia_cov_load_matrix (const gchar *filename, NcmMatrix *data);
static void _nc_data_snia_cov_matrix_to_cov_full (NcDataSNIACov *snia_cov, NcmMatrix *cov, guint i, guint j);
static void _nc_data_snia_cov_diag_to_full_cov (NcDataSNIACov *snia_cov,
                                                NcmVector     *sigma_mag,
                                                NcmVector     *sigma_width,
                                                NcmVector     *sigma_colour,
                                                NcmVector     *diag_mag_width,
                                                NcmVector     *diag_mag_colour,
                                                NcmVector     *diag_width_colour);
static void _nc_data_snia_cov_save_cov_lowertri (NcDataSNIACov *snia_cov);

#define _NC_DATA_SNIA_COV_SET_DATA_INIT(D) \
G_STMT_START { \
  guint64 data_bw = 0; \
  switch (self->cat_version) \
  { \
    case 0: \
    case 1: \
      data_bw = NC_DATA_SNIA_COV_V1_INIT (D); \
      self->data_init = self->data_init | data_bw; \
      if ((self->data_init & NC_DATA_SNIA_COV_V1_ALL) == self->data_init) \
        ncm_data_set_init (NCM_DATA (snia_cov), TRUE); \
      else \
        ncm_data_set_init (NCM_DATA (snia_cov), FALSE); \
      break; \
    case 2: \
      data_bw = NC_DATA_SNIA_COV_V2_INIT (D); \
      self->data_init = self->data_init | data_bw; \
      if ((self->data_init & NC_DATA_SNIA_COV_V2_ALL) == self->data_init) \
        ncm_data_set_init (NCM_DATA (snia_cov), TRUE); \
      else \
        ncm_data_set_init (NCM_DATA (snia_cov), FALSE); \
      break; \
    default: \
      g_assert_not_reached (); \
      break; \
  } \
} G_STMT_END

#define _NC_DATA_SNIA_COV_SET_DATA_INIT_V01(D) \
G_STMT_START { \
  guint64 data_bw = 0; \
  switch (self->cat_version) \
  { \
    case 0: \
    case 1: \
      data_bw = NC_DATA_SNIA_COV_V1_INIT (D); \
      self->data_init = self->data_init | data_bw; \
      if ((self->data_init & NC_DATA_SNIA_COV_V1_ALL) == self->data_init) \
        ncm_data_set_init (NCM_DATA (snia_cov), TRUE); \
      else \
        ncm_data_set_init (NCM_DATA (snia_cov), FALSE); \
      break; \
    case 2: \
      /*g_warning ("_NC_DATA_SNIA_COV_SET_DATA_INIT_V01: element `" #D "' does not exists in version 2");*/ \
      break; \
    default: \
      g_assert_not_reached (); \
      break; \
  } \
} G_STMT_END

#define _NC_DATA_SNIA_COV_SET_DATA_INIT_V2(D) \
G_STMT_START { \
  guint64 data_bw = 0; \
  switch (self->cat_version) \
  { \
    case 0: \
    case 1: \
      /*g_warning ("_NC_DATA_SNIA_COV_SET_DATA_INIT_V2: element `" #D "' does not exists in version 1");*/ \
      break; \
    case 2: \
      data_bw = NC_DATA_SNIA_COV_V2_INIT (D); \
      self->data_init = self->data_init | data_bw; \
      if ((self->data_init & NC_DATA_SNIA_COV_V2_ALL) == self->data_init) \
        ncm_data_set_init (NCM_DATA (snia_cov), TRUE); \
      else \
        ncm_data_set_init (NCM_DATA (snia_cov), FALSE); \
      break; \
    default: \
      g_assert_not_reached (); \
      break; \
  } \
} G_STMT_END

#define _NC_DATA_SNIA_COV_SET_DATA_INIT_ALL \
G_STMT_START { \
  guint64 data_bw = 0; \
  switch (self->cat_version) \
  { \
    case 0: \
    case 1: \
      data_bw = NC_DATA_SNIA_COV_V1_ALL; \
      self->data_init = self->data_init | data_bw; \
      if ((self->data_init & NC_DATA_SNIA_COV_V1_ALL) == self->data_init) \
        ncm_data_set_init (NCM_DATA (snia_cov), TRUE); \
      else \
        ncm_data_set_init (NCM_DATA (snia_cov), FALSE); \
      break; \
    case 2: \
      data_bw = NC_DATA_SNIA_COV_V2_ALL; \
      self->data_init = self->data_init | data_bw; \
      if ((self->data_init & NC_DATA_SNIA_COV_V2_ALL) == self->data_init) \
        ncm_data_set_init (NCM_DATA (snia_cov), TRUE); \
      else \
        ncm_data_set_init (NCM_DATA (snia_cov), FALSE); \
      break; \
    default: \
      g_assert_not_reached (); \
      break; \
  } \
} G_STMT_END


/**
 * nc_data_snia_cov_new_full: (constructor)
 * @filename: catalog file name
 * @use_norma: Whether to use the correct Likelihood normalization
 *
 * Creates a new #NcDataSNIACov object and load with the catalog
 * in @filename. If @use_norma is true the object will use the
 * correct Likelihood normalization when calculating $-2\ln(L)$
 *
 * Returns: (transfer full): the newly created instance of #NcDataSNIACov.
 */
NcDataSNIACov *
nc_data_snia_cov_new_full (const gchar *filename, gboolean use_norma)
{
  NcDataSNIACov *snia_cov = NULL;
  NcDataSNIACovPrivate * self = NULL;
  fitsfile *fptr;
  gchar *desc = NULL;
  gchar comment[FLEN_COMMENT];
  glong nrows, cat_version;
  gint hdutype;
  gdouble mag_cut = NC_DATA_SNIA_COV_MAG_CUT_DEFAULT;
  gint status     = 0;

  if (filename == NULL)
    g_error ("nc_data_snia_cov_load: null filename");

  fits_open_file (&fptr, filename, READONLY, &status);
  NCM_FITS_ERROR (status);

  fits_movabs_hdu (fptr, 2, &hdutype, &status);
  NCM_FITS_ERROR (status);

  if (hdutype != BINARY_TBL)
    g_error ("nc_data_snia_cov_load: NcSNIADistCov catalog is not binary.");

  fits_read_key_lng (fptr, NC_DATA_SNIA_COV_CAT_VERSION, &cat_version, NULL, &status);

  if (status)
  {
    cat_version = 0;
    status      = 0;
  }

  snia_cov = nc_data_snia_cov_new (use_norma, cat_version);
  self = snia_cov->priv;

  if (cat_version != self->cat_version)
    g_error ("nc_data_snia_cov_load: cannot load catalog `%s', object was created for version %d but catalog is version %ld",
        filename,
        self->cat_version,
        cat_version);

  fits_read_key_log (fptr, NC_DATA_SNIA_COV_CAT_HAS_COMPLETE_COV, &self->has_complete_cov, NULL, &status);

  if (status)
  {
    self->has_complete_cov = FALSE;
    status                     = 0;
  }

  fits_read_key_longstr (fptr, NC_DATA_SNIA_COV_CAT_DESC, &desc, comment, &status);

  if (!status)
  {
    NCM_FITS_ERROR (status);

    ncm_data_set_desc (NCM_DATA (snia_cov), desc);
    fits_free_memory (desc, &status);
    NCM_FITS_ERROR (status);
  }
  else
  {
    status = 0;
  }

  fits_read_key_dbl (fptr, NC_DATA_SNIA_COV_MAG_CUT, &mag_cut, comment, &status);

  if (!status)
  {
    NCM_FITS_ERROR (status);
    nc_data_snia_cov_set_mag_cut (snia_cov, mag_cut);
  }
  else
  {
    nc_data_snia_cov_set_mag_cut (snia_cov, NC_DATA_SNIA_COV_MAG_CUT_DEFAULT);
    status = 0;
  }

  fits_get_num_rows (fptr, &nrows, &status);
  NCM_FITS_ERROR (status);
  ncm_data_gauss_cov_set_size (NCM_DATA_GAUSS_COV (snia_cov), nrows);

  /* Setting everything to zero */
  ncm_matrix_set_zero (self->cov_mbc_mbc);
  ncm_matrix_set_zero (self->cov_full);
  ncm_vector_set_zero (self->cov_full_diag);
  ncm_vector_set_zero (self->cov_packed);

  if (self->cat_version == 0)
  {
    printf ("AQUI\n");fflush (stdout);
    nc_data_snia_cov_load_V0 (snia_cov, fptr);
  }
  else if (self->cat_version == 1)
  {
    NcmMatrix *cov = ncm_matrix_new (self->mu_len, self->mu_len);
    NcmVector *data[NC_DATA_SNIA_COV_V1_LENGTH];
    guint i;

    data[NC_DATA_SNIA_COV_V1_ZCMB]           = self->z_cmb;
    data[NC_DATA_SNIA_COV_V1_ZHE]            = self->z_he;
    data[NC_DATA_SNIA_COV_V1_SIGMA_Z]        = self->sigma_z;
    data[NC_DATA_SNIA_COV_V1_MAG]            = self->mag;
    data[NC_DATA_SNIA_COV_V1_WIDTH]          = self->width;
    data[NC_DATA_SNIA_COV_V1_COLOUR]         = self->colour;
    data[NC_DATA_SNIA_COV_V1_THIRDPAR]       = self->thirdpar;
    data[NC_DATA_SNIA_COV_V1_SIGMA_THIRDPAR] = self->sigma_thirdpar;

    for (i = 0; i < NC_DATA_SNIA_COV_V1_LENGTH; i++)
    {
      fits_read_col_dbl (fptr, i + 1, 1, 1, self->mu_len, GSL_NAN,
                         ncm_vector_ptr (data[i], 0), NULL,
                         &status);
      NCM_FITS_ERROR (status);
    }

    fits_read_col_uint (fptr, NC_DATA_SNIA_COV_V1_ABSMAG_SET + 1, 1, 1,
                        self->mu_len,
                        0, &g_array_index (self->dataset, guint32, 0),
                        NULL, &status);
    NCM_FITS_ERROR (status);
    nc_data_snia_cov_set_abs_mag_set (snia_cov, self->dataset);

    for (i = NC_DATA_SNIA_COV_V1_MAG_MAG; i < NC_DATA_SNIA_COV_V1_TOTAL_LENGTH; i++)
    {
      fits_read_col_dbl (fptr, i + 1, 1, 1, self->mu_len * self->mu_len,
                         GSL_NAN, ncm_matrix_ptr (cov, 0, 0),
                         NULL, &status);
      NCM_FITS_ERROR (status);

      switch (i)
      {
        case NC_DATA_SNIA_COV_V1_MAG_MAG:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 0, 0);
          break;
        case NC_DATA_SNIA_COV_V1_MAG_WIDTH:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 0, 1);
          break;
        case NC_DATA_SNIA_COV_V1_MAG_COLOUR:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 0, 2);
          break;
        case NC_DATA_SNIA_COV_V1_WIDTH_WIDTH:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 1, 1);
          break;
        case NC_DATA_SNIA_COV_V1_WIDTH_COLOUR:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 1, 2);
          break;
        case NC_DATA_SNIA_COV_V1_COLOUR_COLOUR:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 2, 2);
          break;
        default:
          g_assert_not_reached ();
          break;
      }
    }
  }
  else if (self->cat_version == 2)
  {
    NcmVector *data[NC_DATA_SNIA_COV_V2_LENGTH];
    guint i;

    data[NC_DATA_SNIA_COV_V2_ZHD]              = self->z_hd;
    data[NC_DATA_SNIA_COV_V2_SIGMA_ZHD]        = self->sigma_z_hd;
    data[NC_DATA_SNIA_COV_V2_ZCMB]             = self->z_cmb;
    data[NC_DATA_SNIA_COV_V2_SIGMA_ZCMB]       = self->sigma_z_cmb;
    data[NC_DATA_SNIA_COV_V2_ZHE]              = self->z_he;
    data[NC_DATA_SNIA_COV_V2_SIGMA_ZHE]        = self->sigma_z_he;
    data[NC_DATA_SNIA_COV_V2_MAG_B_CORR]       = self->mag_b_corr;
    data[NC_DATA_SNIA_COV_V2_SIGMA_MAG_B_CORR] = self->sigma_mag_b_corr;
    data[NC_DATA_SNIA_COV_V2_MU_CORR]          = self->mu_corr;
    data[NC_DATA_SNIA_COV_V2_SIGMA_MU_CORR]    = self->sigma_mu_corr;
    data[NC_DATA_SNIA_COV_V2_CEPH_DIST]        = self->ceph_dist;
    data[NC_DATA_SNIA_COV_V2_COLOUR]           = self->colour;
    data[NC_DATA_SNIA_COV_V2_SIGMA_COLOUR]     = self->sigma_colour;
    data[NC_DATA_SNIA_COV_V2_WIDTH]            = self->width;
    data[NC_DATA_SNIA_COV_V2_SIGMA_WIDTH]      = self->sigma_width;
    data[NC_DATA_SNIA_COV_V2_MAG]              = self->mag;
    data[NC_DATA_SNIA_COV_V2_SIGMA_MAG]        = self->sigma_mag;
    data[NC_DATA_SNIA_COV_V2_AMPL]             = self->ampl;
    data[NC_DATA_SNIA_COV_V2_SIGMA_AMPL]       = self->sigma_ampl;
    data[NC_DATA_SNIA_COV_V2_WIDTH_COLOUR]     = self->cov_width_colour;
    data[NC_DATA_SNIA_COV_V2_WIDTH_AMPL]       = self->cov_width_ampl;
    data[NC_DATA_SNIA_COV_V2_COLOUR_AMPL]      = self->cov_colour_ampl;
    data[NC_DATA_SNIA_COV_V2_RA]               = self->ra;
    data[NC_DATA_SNIA_COV_V2_DEC]              = self->dec;
    data[NC_DATA_SNIA_COV_V2_HOST_RA]          = self->host_ra;
    data[NC_DATA_SNIA_COV_V2_HOST_DEC]         = self->host_dec;
    data[NC_DATA_SNIA_COV_V2_HOST_ANG_SEP]     = self->host_ang_sep;
    data[NC_DATA_SNIA_COV_V2_VPEC]             = self->vpec;
    data[NC_DATA_SNIA_COV_V2_SIGMA_VPEC]       = self->sigma_vpec;

    for (i = 0; i < NC_DATA_SNIA_COV_V2_LENGTH; i++)
    {
      if ((i == NC_DATA_SNIA_COV_V2_IS_CALIB) || (i == NC_DATA_SNIA_COV_V2_USED_IN_SH0ES))
        continue;

      g_assert (data[i] != NULL);
      g_assert (NCM_IS_VECTOR (data[i]));

      fits_read_col_dbl (fptr, i + 1, 1, 1, self->mu_len, GSL_NAN,
                         ncm_vector_ptr (data[i], 0), NULL,
                         &status);
      NCM_FITS_ERROR (status);
    }

    fits_read_col_uint (fptr, NC_DATA_SNIA_COV_V2_IS_CALIB + 1, 1, 1,
                        self->mu_len,
                        0, &g_array_index (self->is_calib, guint32, 0),
                        NULL, &status);
    NCM_FITS_ERROR (status);

    fits_read_col_uint (fptr, NC_DATA_SNIA_COV_V2_USED_IN_SH0ES + 1, 1, 1,
                        self->mu_len,
                        0, &g_array_index (self->used_in_sh0es, guint32, 0),
                        NULL, &status);
    NCM_FITS_ERROR (status);

    fits_read_col_dbl (fptr, NC_DATA_SNIA_COV_V2_MAG_B_CORR_MAG_B_CORR + 1, 1, 1, self->mu_len * self->mu_len,
        GSL_NAN, ncm_matrix_ptr (self->cov_mbc_mbc, 0, 0),
        NULL, &status);
    NCM_FITS_ERROR (status);
  }
  else
  {
    g_assert_not_reached ();
  }

  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);

  _NC_DATA_SNIA_COV_SET_DATA_INIT_ALL;
  nc_data_snia_cov_set_cov_full (snia_cov, self->cov_full);

  return snia_cov;
}

/**
 * nc_data_snia_cov_new_from_cat_id: (constructor)
 * @id: Catalog id
 * @use_norma: Whether to use the correct Likelihood normalization
 *
 * Creates a new NcDataSNIACov using the catalog from @id.
 *
 * Returns: (transfer full): the newly created instance of #NcDataSNIACov.
 */
NcDataSNIACov *
nc_data_snia_cov_new_from_cat_id (NcDataSNIAId id, gboolean use_norma)
{
  NcDataSNIACov *snia_cov = NULL;
  gchar *full_filename = nc_data_snia_cov_get_catalog_by_id (id);

#ifdef NUMCOSMO_HAVE_CFITSIO
  snia_cov = nc_data_snia_cov_new_full (full_filename, FALSE);
#else
  g_error ("nc_data_snia_cov_new_from_cat_id: catalog load not available, recompile "PACKAGE_NAME " with cfitsio support.");
#endif /* NUMCOSMO_HAVE_CFITSIO */
  g_free (full_filename);

  return snia_cov;
}


/**
 * nc_data_snia_cov_sigma_int_len:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the number of different intrinsic sigma parameters in the
 * catalog.
 *
 * Returns: The number of different sigma_int.
 */
guint
nc_data_snia_cov_sigma_int_len (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->dataset_len;
}

/**
 * nc_data_snia_cov_snia_len:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the number of supernovae in the
 * catalog.
 *
 * Returns: The total number of supernovae in the catalog.
 */
guint
nc_data_snia_cov_snia_len (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->mu_len;
}

static void
_nc_data_snia_cov_set_size (NcmDataGaussCov *gauss, guint mu_len)
{
  /* Chain up : start */
  NCM_DATA_GAUSS_COV_CLASS (nc_data_snia_cov_parent_class)->set_size (gauss, mu_len);
  {
    NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (gauss);
    NcDataSNIACovPrivate * const self = snia_cov->priv;

    if ((mu_len == 0) || (mu_len != self->mu_len))
    {
      ncm_vector_clear (&self->z_hd);
      ncm_vector_clear (&self->z_cmb);
      ncm_vector_clear (&self->z_he);

      ncm_vector_clear (&self->mag_b_corr);
      ncm_vector_clear (&self->mu_corr);
      ncm_vector_clear (&self->ceph_dist);
      ncm_vector_clear (&self->ampl);
      ncm_vector_clear (&self->vpec);
      ncm_vector_clear (&self->mag);
      ncm_vector_clear (&self->width);
      ncm_vector_clear (&self->colour);
      ncm_vector_clear (&self->thirdpar);

      ncm_vector_clear (&self->ra);
      ncm_vector_clear (&self->dec);
      ncm_vector_clear (&self->host_ra);
      ncm_vector_clear (&self->host_dec);
      ncm_vector_clear (&self->host_ang_sep);

      ncm_vector_clear (&self->width_true);
      ncm_vector_clear (&self->colour_true);
      ncm_vector_clear (&self->mag_width_colour);

      ncm_vector_clear (&self->sigma_z_hd);
      ncm_vector_clear (&self->sigma_z_cmb);
      ncm_vector_clear (&self->sigma_z_he);
      ncm_vector_clear (&self->sigma_z);

      ncm_vector_clear (&self->sigma_mag_b_corr);
      ncm_vector_clear (&self->sigma_mu_corr);
      ncm_vector_clear (&self->sigma_colour);
      ncm_vector_clear (&self->sigma_width);
      ncm_vector_clear (&self->sigma_mag);
      ncm_vector_clear (&self->sigma_ampl);
      ncm_vector_clear (&self->sigma_vpec);
      ncm_vector_clear (&self->sigma_thirdpar);

      ncm_vector_clear (&self->cov_width_colour);
      ncm_vector_clear (&self->cov_width_ampl);
      ncm_vector_clear (&self->cov_colour_ampl);

      ncm_vector_clear (&self->cov_packed);
      ncm_matrix_clear (&self->cov_mbc_mbc);
      ncm_vector_clear (&self->cov_full_diag);
      ncm_matrix_clear (&self->cov_full);

      ncm_matrix_clear (&self->inv_cov_mm_LU);
      ncm_matrix_clear (&self->inv_cov_mm);

      g_clear_pointer (&self->is_calib, g_array_unref);
      g_clear_pointer (&self->used_in_sh0es, g_array_unref);

      if (self->dataset != NULL)
      {
        g_array_unref (self->dataset);
        g_array_unref (self->dataset_size);
        self->dataset      = NULL;
        self->dataset_size = NULL;
        self->dataset_len  = 0;
      }

      ncm_data_set_init (NCM_DATA (snia_cov), FALSE);
    }

    if ((mu_len > 0) && (mu_len != self->mu_len))
    {
      self->mu_len       = mu_len;
      self->uppertri_len = (mu_len * (mu_len + 1)) / 2;

      g_assert_cmpuint (ncm_vector_len (gauss->y), ==, mu_len);

      switch (self->cat_version)
      {
        case 0:
        case 1:
          self->mag        = ncm_vector_ref (gauss->y);
          self->mag_b_corr = ncm_vector_new (mu_len);
          break;
        case 2:
          self->mag        = ncm_vector_new (mu_len);
          self->mag_b_corr = ncm_vector_ref (gauss->y);
          break;
        default:
          g_assert_not_reached ();
          break;
      }

      self->z_hd  = ncm_vector_new (mu_len);
      self->z_cmb = ncm_vector_new (mu_len);
      self->z_he  = ncm_vector_new (mu_len);

      self->mu_corr   = ncm_vector_new (mu_len);
      self->ceph_dist = ncm_vector_new (mu_len);
      self->ampl      = ncm_vector_new (mu_len);
      self->vpec      = ncm_vector_new (mu_len);
      self->width     = ncm_vector_new (mu_len);
      self->colour    = ncm_vector_new (mu_len);
      self->thirdpar  = ncm_vector_new (mu_len);

      self->ra           = ncm_vector_new (mu_len);
      self->dec          = ncm_vector_new (mu_len);
      self->host_ra      = ncm_vector_new (mu_len);
      self->host_dec     = ncm_vector_new (mu_len);
      self->host_ang_sep = ncm_vector_new (mu_len);

      self->width_true       = ncm_vector_new (mu_len);
      self->colour_true      = ncm_vector_new (mu_len);
      self->mag_width_colour = ncm_vector_new (mu_len * 3);

      self->sigma_z_hd  = ncm_vector_new (mu_len);
      self->sigma_z_cmb = ncm_vector_new (mu_len);
      self->sigma_z_he  = ncm_vector_new (mu_len);
      self->sigma_z     = ncm_vector_new (mu_len);

      self->sigma_mag_b_corr = ncm_vector_new (mu_len);
      self->sigma_mu_corr    = ncm_vector_new (mu_len);
      self->sigma_colour     = ncm_vector_new (mu_len);
      self->sigma_width      = ncm_vector_new (mu_len);
      self->sigma_mag        = ncm_vector_new (mu_len);
      self->sigma_ampl       = ncm_vector_new (mu_len);
      self->sigma_vpec       = ncm_vector_new (mu_len);
      self->sigma_thirdpar   = ncm_vector_new (mu_len);

      self->cov_width_colour = ncm_vector_new (mu_len);
      self->cov_width_ampl   = ncm_vector_new (mu_len);
      self->cov_colour_ampl  = ncm_vector_new (mu_len);

      self->cov_full_diag = ncm_vector_new (3 * mu_len);
      self->cov_mbc_mbc   = ncm_matrix_new (mu_len, mu_len);
      self->cov_full      = ncm_matrix_new (3 * mu_len, 3 * mu_len);
      self->cov_packed    = ncm_vector_new (self->uppertri_len * NC_DATA_SNIA_COV_ORDER_LENGTH);

      self->inv_cov_mm    = ncm_matrix_new (mu_len, mu_len);
      self->inv_cov_mm_LU = ncm_matrix_new (mu_len, mu_len);

      self->is_calib      = g_array_sized_new (FALSE, FALSE, sizeof (guint32), mu_len);
      self->used_in_sh0es = g_array_sized_new (FALSE, FALSE, sizeof (guint32), mu_len);

      self->dataset      = g_array_sized_new (FALSE, FALSE, sizeof (guint32), mu_len);
      self->dataset_size = g_array_new (FALSE, TRUE, sizeof (guint32));
      self->dataset_len  = 0;

      g_array_set_size (self->dataset, mu_len);
      g_array_set_size (self->is_calib, mu_len);
      g_array_set_size (self->used_in_sh0es, mu_len);

      ncm_data_set_init (NCM_DATA (snia_cov), FALSE);
    }
  }
}

/**
 * nc_data_snia_cov_peek_z_hd:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the redshift in the CMB frame corrected for peculiar velocities
 * $z_\mathrm{hd}$ #NcmVector.
 *
 * Returns: (transfer none): the $z_\mathrm{hd}$ #NcmVector
 */
NcmVector *
nc_data_snia_cov_peek_z_hd (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->z_hd;
}

/**
 * nc_data_snia_cov_peek_z_cmb:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the $z_\mathrm{cmb}$ #NcmVector.
 *
 * Returns: (transfer none): the $z_\mathrm{cmb}$ #NcmVector
 */
NcmVector *
nc_data_snia_cov_peek_z_cmb (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->z_cmb;
}

/**
 * nc_data_snia_cov_peek_z_he:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the $z_\mathrm{he}$ #NcmVector.
 *
 * Returns: (transfer none): the $z_\mathrm{he}$ #NcmVector
 */
NcmVector *
nc_data_snia_cov_peek_z_he (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->z_he;
}

/**
 * nc_data_snia_cov_peek_sigma_z:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the $\sigma_z$ #NcmVector.
 *
 * Returns: (transfer none): the $\sigma_z$ #NcmVector
 */
NcmVector *
nc_data_snia_cov_peek_sigma_z (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->sigma_z;
}

/**
 * nc_data_snia_cov_peek_mag:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the magnitude #NcmVector.
 *
 * Returns: (transfer none): the magnitude #NcmVector
 */
NcmVector *
nc_data_snia_cov_peek_mag (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->mag;
}

/**
 * nc_data_snia_cov_peek_width:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the width #NcmVector.
 *
 * Returns: (transfer none): the width #NcmVector
 */
NcmVector *
nc_data_snia_cov_peek_width (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->width;
}

/**
 * nc_data_snia_cov_peek_colour:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the colour #NcmVector.
 *
 * Returns: (transfer none): the colour #NcmVector
 */
NcmVector *
nc_data_snia_cov_peek_colour (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->colour;
}

/**
 * nc_data_snia_cov_peek_thirdpar:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the thirdpar #NcmVector.
 *
 * Returns: (transfer none): the thirdpar #NcmVector
 */
NcmVector *
nc_data_snia_cov_peek_thirdpar (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->thirdpar;
}

/**
 * nc_data_snia_cov_peek_ceph_dist:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the Cepheid distances #NcmVector.
 *
 * Returns: (transfer none): the Cepheid distances #NcmVector.
 */
NcmVector *
nc_data_snia_cov_peek_ceph_dist (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->ceph_dist;
}

/**
 * nc_data_snia_cov_peek_abs_mag_set:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the dataset array.
 *
 * Returns: (transfer none) (array) (element-type guint32): the dataset array
 */
GArray *
nc_data_snia_cov_peek_abs_mag_set (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->dataset;
}

/**
 * nc_data_snia_cov_set_mag_cut:
 * @snia_cov: a #NcDataSNIACov
 * @mag_cut: a double
 *
 * Sets the absolute magnitude cut value.
 *
 */
void
nc_data_snia_cov_set_mag_cut (NcDataSNIACov *snia_cov, const gdouble mag_cut)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  self->mag_cut = mag_cut;
}

/**
 * nc_data_snia_cov_get_mag_cut:
 * @snia_cov: a #NcDataSNIACov
 *
 * Returns: the current absolute magnitude cut value.
 */
gdouble
nc_data_snia_cov_get_mag_cut (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  return self->mag_cut;
}

/**
 * nc_data_snia_cov_peek_cov_full:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the cov_full #NcmMatrix.
 *
 * Returns: (transfer none): the cov_full #NcmMatrix
 */
NcmMatrix *
nc_data_snia_cov_peek_cov_full (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (self->mu_len > 0)
    _nc_data_snia_cov_restore_full_cov (snia_cov);

  return self->cov_full;
}

/**
 * nc_data_snia_cov_peek_cov_packed:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the packed version of the full covariance.
 *
 * Returns: (transfer none): an #NcmVector containing the packed covariance matrix.
 */
NcmVector *
nc_data_snia_cov_peek_cov_packed (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;

  return self->cov_packed;
}

/**
 * nc_data_snia_cov_peek_cov_mbc_mbc:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets the cov_full #NcmMatrix.
 *
 * Returns: (transfer none): the covariance of mag b corr #NcmMatrix
 */
NcmMatrix *
nc_data_snia_cov_peek_cov_mbc_mbc (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;

  return self->cov_mbc_mbc;
}

/**
 * nc_data_snia_cov_peek_dataset:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets an array containing the dataset id for each SNIa.
 *
 * Returns: (transfer none) (element-type guint32): an #GArray containing the dataset id for each SNIa.
 */
GArray *
nc_data_snia_cov_peek_dataset (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;

  return self->dataset;
}

/**
 * nc_data_snia_cov_peek_is_calib:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets an array containing whether the SNIa is a calibrator.
 *
 * Returns: (transfer none) (element-type guint32): an #GArray containing whether the SNIa is a calibrator.
 */
GArray *
nc_data_snia_cov_peek_is_calib (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;

  return self->is_calib;
}

/**
 * nc_data_snia_cov_peek_used_in_sh0es:
 * @snia_cov: a #NcDataSNIACov
 *
 * Gets an array containing whether the SNIa was used in SH0ES.
 *
 * Returns: (transfer none) (element-type guint32): an #GArray containing containing whether the SNIa was used in SH0ES.
 */
GArray *
nc_data_snia_cov_peek_used_in_sh0es (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;

  return self->used_in_sh0es;
}

/**
 * nc_data_snia_cov_set_z_hd:
 * @snia_cov: a #NcDataSNIACov
 * @z_hd: the $z_\mathrm{hd}$ #NcmVector
 *
 * Sets the $z_\mathrm{hd}$ vector to @z_cmb.
 *
 */
void
nc_data_snia_cov_set_z_hd (NcDataSNIACov *snia_cov, NcmVector *z_hd)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (z_hd != self->z_hd)
  {
    g_assert_cmpuint (self->mu_len, ==, ncm_vector_len (z_hd));
    ncm_vector_free (self->z_hd);
    self->z_hd = ncm_vector_ref (z_hd);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT_V2 (ZHD);
}

/**
 * nc_data_snia_cov_set_z_cmb:
 * @snia_cov: a #NcDataSNIACov
 * @z_cmb: the $z_\mathrm{cmb}$ #NcmVector
 *
 * Sets the $z_\mathrm{cmb}$ vector to @z_cmb.
 *
 */
void
nc_data_snia_cov_set_z_cmb (NcDataSNIACov *snia_cov, NcmVector *z_cmb)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (z_cmb != self->z_cmb)
  {
    g_assert_cmpuint (self->mu_len, ==, ncm_vector_len (z_cmb));
    ncm_vector_free (self->z_cmb);
    self->z_cmb = ncm_vector_ref (z_cmb);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT (ZCMB);
}

/**
 * nc_data_snia_cov_set_z_he:
 * @snia_cov: a #NcDataSNIACov
 * @z_he: the $z_\mathrm{he}$ #NcmVector
 *
 * Sets the $z_\mathrm{he}$ vector to @z_he.
 *
 */
void
nc_data_snia_cov_set_z_he (NcDataSNIACov *snia_cov, NcmVector *z_he)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (z_he != self->z_he)
  {
    g_assert_cmpuint (self->mu_len, ==, ncm_vector_len (z_he));
    ncm_vector_free (self->z_he);
    self->z_he = ncm_vector_ref (z_he);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT (ZHE);
}

/**
 * nc_data_snia_cov_set_sigma_z:
 * @snia_cov: a #NcDataSNIACov
 * @sigma_z: the $\sigma_z$ #NcmVector
 *
 * Sets the $\sigma_z$ vector to @sigma_z.
 *
 */
void
nc_data_snia_cov_set_sigma_z (NcDataSNIACov *snia_cov, NcmVector *sigma_z)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (sigma_z != self->sigma_z)
  {
    g_assert_cmpuint (self->mu_len, ==, ncm_vector_len (sigma_z));
    ncm_vector_free (self->sigma_z);
    self->sigma_z = ncm_vector_ref (sigma_z);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT_V01 (SIGMA_Z);
}

/**
 * nc_data_snia_cov_set_mag:
 * @snia_cov: a #NcDataSNIACov
 * @mag: the magnitude #NcmVector
 *
 * Sets the magnitude vector to @mag.
 *
 */
void
nc_data_snia_cov_set_mag (NcDataSNIACov *snia_cov, NcmVector *mag)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (mag != self->mag)
  {
    g_assert_cmpuint (self->mu_len, ==, ncm_vector_len (mag));
    ncm_vector_free (self->mag);
    self->mag = ncm_vector_ref (mag);

    if (self->cat_version != 2)
      ncm_vector_substitute (&NCM_DATA_GAUSS_COV (snia_cov)->y, mag, TRUE);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT (MAG);
}

/**
 * nc_data_snia_cov_set_mag_b_corr:
 * @snia_cov: a #NcDataSNIACov
 * @mag_b_corr: the magnitude b corrected #NcmVector
 *
 * Sets the magnitude b corrected vector to @mag_b_corr.
 *
 */
void
nc_data_snia_cov_set_mag_b_corr (NcDataSNIACov *snia_cov, NcmVector *mag_b_corr)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (mag_b_corr != self->mag_b_corr)
  {
    g_assert_cmpuint (self->mu_len, ==, ncm_vector_len (mag_b_corr));
    ncm_vector_free (self->mag_b_corr);
    self->mag_b_corr = ncm_vector_ref (mag_b_corr);

    if (self->cat_version == 2)
      ncm_vector_substitute (&NCM_DATA_GAUSS_COV (snia_cov)->y, mag_b_corr, TRUE);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT_V2 (MAG_B_CORR);
}

/**
 * nc_data_snia_cov_set_ceph_dist:
 * @snia_cov: a #NcDataSNIACov
 * @ceph_dist: the Cepheid distances #NcmVector
 *
 * Sets the Cepheid distances vector to @ceph_dist.
 *
 */
void
nc_data_snia_cov_set_ceph_dist (NcDataSNIACov *snia_cov, NcmVector *ceph_dist)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (ceph_dist != self->ceph_dist)
  {
    g_assert_cmpuint (self->mu_len, ==, ncm_vector_len (ceph_dist));
    ncm_vector_free (self->ceph_dist);
    self->ceph_dist = ncm_vector_ref (ceph_dist);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT_V2 (CEPH_DIST);
}


/**
 * nc_data_snia_cov_set_width:
 * @snia_cov: a #NcDataSNIACov
 * @width: the width #NcmVector
 *
 * Sets the width vector to @width.
 *
 */
void
nc_data_snia_cov_set_width (NcDataSNIACov *snia_cov, NcmVector *width)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (width != self->width)
  {
    g_assert_cmpuint (self->mu_len, ==, ncm_vector_len (width));
    ncm_vector_free (self->width);
    self->width = ncm_vector_ref (width);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT (WIDTH);
}

/**
 * nc_data_snia_cov_set_colour:
 * @snia_cov: a #NcDataSNIACov
 * @colour: the colour #NcmVector
 *
 * Sets the colour vector to @colour.
 *
 */
void
nc_data_snia_cov_set_colour (NcDataSNIACov *snia_cov, NcmVector *colour)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (colour != self->colour)
  {
    g_assert_cmpuint (self->mu_len, ==, ncm_vector_len (colour));
    ncm_vector_free (self->colour);
    self->colour = ncm_vector_ref (colour);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT (COLOUR);
}

/**
 * nc_data_snia_cov_set_thirdpar:
 * @snia_cov: a #NcDataSNIACov
 * @thirdpar: the thirdpar #NcmVector
 *
 * Sets the thirdpar vector to @thirdpar.
 *
 */
void
nc_data_snia_cov_set_thirdpar (NcDataSNIACov *snia_cov, NcmVector *thirdpar)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (thirdpar != self->thirdpar)
  {
    g_assert_cmpuint (self->mu_len, ==, ncm_vector_len (thirdpar));
    ncm_vector_free (self->thirdpar);
    self->thirdpar = ncm_vector_ref (thirdpar);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT_V01 (THIRDPAR);
}

/**
 * nc_data_snia_cov_set_abs_mag_set:
 * @snia_cov: a #NcDataSNIACov
 * @abs_mag_set: (in) (array) (element-type guint32): the absolute magnitude set
 *
 * Sets the array containing the indexes labeling to which set each SNIa
 * belongs.
 *
 */
void
nc_data_snia_cov_set_abs_mag_set (NcDataSNIACov *snia_cov, GArray *abs_mag_set)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  guint total;
  guint max_id = 0;
  guint i;

  g_assert_cmpuint (self->mu_len, ==, abs_mag_set->len);

  if (abs_mag_set != self->dataset)
  {
    g_assert_cmpuint (g_array_get_element_size (abs_mag_set), ==, sizeof (guint32));
    g_array_unref (self->dataset);
    self->dataset = g_array_ref (abs_mag_set);
  }

  for (i = 0; i < abs_mag_set->len; i++)
  {
    guint32 vi = g_array_index (abs_mag_set, guint32, i);

    max_id = GSL_MAX (vi, max_id);
  }

  self->dataset_len = max_id + 1;
  g_array_set_size (self->dataset_size, self->dataset_len);

  for (i = 0; i < abs_mag_set->len; i++)
  {
    guint32 vi = g_array_index (abs_mag_set, guint32, i);

    g_array_index (self->dataset_size, guint32, vi)++;
  }

  total = 0;

  for (i = 0; i < self->dataset_size->len; i++)
  {
    total += g_array_index (self->dataset_size, guint32, i);
  }

  g_assert_cmpuint (total, ==, self->mu_len);
  _NC_DATA_SNIA_COV_SET_DATA_INIT_V01 (ABSMAG_SET);
}

/**
 * nc_data_snia_cov_set_is_calib:
 * @snia_cov: a #NcDataSNIACov
 * @is_calib: (in) (array) (element-type guint32): whether the SNIa is a calibrator
 *
 * Sets the array containing whether the SNIa is a calibrator.
 *
 */
void
nc_data_snia_cov_set_is_calib (NcDataSNIACov *snia_cov, GArray *is_calib)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;

  g_assert_cmpuint (self->mu_len, ==, is_calib->len);

  if (is_calib != self->is_calib)
  {
    g_assert_cmpuint (g_array_get_element_size (is_calib), ==, sizeof (guint32));
    g_array_unref (self->is_calib);
    self->is_calib = g_array_ref (is_calib);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT_V2 (IS_CALIB);
}

/**
 * nc_data_snia_cov_set_used_in_sh0es:
 * @snia_cov: a #NcDataSNIACov
 * @used_in_sh0es: (in) (array) (element-type guint32): whether the SNIa was used in SH0ES
 *
 * Sets the array containing whether the SNIa was used in SH0ES.
 *
 */
void
nc_data_snia_cov_set_used_in_sh0es (NcDataSNIACov *snia_cov, GArray *used_in_sh0es)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;

  g_assert_cmpuint (self->mu_len, ==, used_in_sh0es->len);

  if (used_in_sh0es != self->used_in_sh0es)
  {
    g_assert_cmpuint (g_array_get_element_size (used_in_sh0es), ==, sizeof (guint32));
    g_array_unref (self->used_in_sh0es);
    self->used_in_sh0es = g_array_ref (used_in_sh0es);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT_V2 (USED_IN_SH0ES);
}

/**
 * nc_data_snia_cov_set_cov_full:
 * @snia_cov: a #NcDataSNIACov
 * @cov_full: the full covariance #NcmMatrix
 *
 * Sets the full covariance for the system, the size of @cov_full,
 * must match the system size.
 *
 */
void
nc_data_snia_cov_set_cov_full (NcDataSNIACov *snia_cov, NcmMatrix *cov_full)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  const guint mu_len = self->mu_len;
  const guint tmu_len = 3 * mu_len;
  guint i, j, ij;

  if (self->cov_full != cov_full)
  {
    g_assert_cmpuint (ncm_matrix_nrows (cov_full), ==, tmu_len);
    g_assert_cmpuint (ncm_matrix_ncols (cov_full), ==, tmu_len);

    ncm_matrix_memcpy (self->cov_full, cov_full);
  }

  /* Filling the cov_packed from cov_full. */
  {
    ij = 0;

    for (i = 0; i < self->mu_len; i++)
    {
      for (j = i; j < self->mu_len; j++)
      {
        const gdouble mag_mag_i = 0.5 * (ncm_matrix_get (self->cov_full, 0 * mu_len + i, 0 * mu_len + j) +
                                         ncm_matrix_get (self->cov_full, 0 * mu_len + j, 0 * mu_len + i));
        const gdouble mag_width_i = 0.5 * (ncm_matrix_get (self->cov_full, 0 * mu_len + i, 1 * mu_len + j) +
                                           ncm_matrix_get (self->cov_full, 0 * mu_len + j, 1 * mu_len + i));
        const gdouble mag_colour_i = 0.5 * (ncm_matrix_get (self->cov_full, 0 * mu_len + i, 2 * mu_len + j) +
                                            ncm_matrix_get (self->cov_full, 0 * mu_len + j, 2 * mu_len + i));
        const gdouble width_width_i = 0.5 * (ncm_matrix_get (self->cov_full, 1 * mu_len + i, 1 * mu_len + j) +
                                             ncm_matrix_get (self->cov_full, 1 * mu_len + j, 1 * mu_len + i));
        const gdouble width_colour_i = 0.5 * (ncm_matrix_get (self->cov_full, 1 * mu_len + i, 2 * mu_len + j) +
                                              ncm_matrix_get (self->cov_full, 1 * mu_len + j, 2 * mu_len + i));
        const gdouble colour_colour_i = 0.5 * (ncm_matrix_get (self->cov_full, 2 * mu_len + i, 2 * mu_len + j) +
                                               ncm_matrix_get (self->cov_full, 2 * mu_len + j, 2 * mu_len + i));

        ncm_vector_set (self->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_MAG_MAG,       mag_mag_i);
        ncm_vector_set (self->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_MAG_WIDTH,     mag_width_i);
        ncm_vector_set (self->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_MAG_COLOUR,    mag_colour_i);
        ncm_vector_set (self->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_WIDTH_WIDTH,   width_width_i);
        ncm_vector_set (self->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_WIDTH_COLOUR,  width_colour_i);
        ncm_vector_set (self->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_COLOUR_COLOUR, colour_colour_i);

        ij++;
      }
    }
  }

  _nc_data_snia_cov_save_cov_lowertri (snia_cov);

  if (self->has_complete_cov)
  {
    gint ret;

    ret = ncm_matrix_cholesky_decomp (self->cov_full, 'U');

    if (ret != 0)
      g_error ("nc_data_snia_cov_set_cov_full[ncm_matrix_cholesky_decomp]: %d.", ret);

    ret = ncm_matrix_cholesky_inverse (self->cov_full, 'U');

    if (ret != 0)
      g_error ("nc_data_snia_cov_set_cov_full[ncm_matrix_cholesky_inverse]: %d.", ret);

    ncm_matrix_set_zero (self->inv_cov_mm);

    for (i = 0; i < mu_len; i++)
    {
      for (j = i; j < mu_len; j++)
      {
        const gdouble inv_cov_full_ij = ncm_matrix_get (self->cov_full, i, j);

        ncm_matrix_set (self->inv_cov_mm, i, j, inv_cov_full_ij);
      }
    }
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT_V01 (COV_FULL);
}

/**
 * nc_data_snia_cov_set_cov_mbc_mbc:
 * @snia_cov: a #NcDataSNIACov
 * @cov_mbc_mbc: the mag b corr covariance #NcmMatrix
 *
 * Sets the the mag b corr covariance to @cov_mbc_mbc.
 *
 */
void
nc_data_snia_cov_set_cov_mbc_mbc (NcDataSNIACov *snia_cov, NcmMatrix *cov_mbc_mbc)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  if (cov_mbc_mbc != self->cov_mbc_mbc)
  {
    g_assert_cmpuint (self->mu_len, ==, ncm_matrix_nrows (cov_mbc_mbc));
    g_assert_cmpuint (self->mu_len, ==, ncm_matrix_ncols (cov_mbc_mbc));
    ncm_matrix_clear (&self->cov_mbc_mbc);
    self->cov_mbc_mbc = ncm_matrix_ref (cov_mbc_mbc);
  }

  _NC_DATA_SNIA_COV_SET_DATA_INIT_V2 (MAG_B_CORR_MAG_B_CORR);
}

/**
 * nc_data_snia_cov_load_txt:
 * @snia_cov: a #NcDataSNIACov
 * @filename: text file name
 *
 * Loads SNIa data from the key-file @filename.
 *
 */
void
nc_data_snia_cov_load_txt (NcDataSNIACov *snia_cov, const gchar *filename)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  GKeyFile *snia_keyfile = g_key_file_new ();
  GError *error          = NULL;
  NcmMatrix *cov         = NULL;

  if (!g_key_file_load_from_file (snia_keyfile, filename, G_KEY_FILE_NONE, &error))
    g_error ("nc_data_snia_cov_load: invalid configuration: %s %s",
             filename,
             error->message);

  if (!g_key_file_has_key (snia_keyfile,
                           NC_DATA_SNIA_COV_DATA_GROUP,
                           NC_DATA_SNIA_COV_DATA_LEN_KEY,
                           &error))
  {
    g_error ("nc_data_snia_cov_load: invalid %s key file, it must define at least the data length file key ["NC_DATA_SNIA_COV_DATA_LEN_KEY "]", filename);
  }
  else
  {
    guint64 mu_len = g_key_file_get_uint64 (snia_keyfile,
                                            NC_DATA_SNIA_COV_DATA_GROUP,
                                            NC_DATA_SNIA_COV_DATA_LEN_KEY,
                                            &error);

    ncm_data_gauss_cov_set_size (NCM_DATA_GAUSS_COV (snia_cov), mu_len);
    cov = ncm_matrix_new (mu_len, mu_len);
  }

  /* Setting everything to zero */
  ncm_matrix_set_zero (self->cov_mbc_mbc);
  ncm_matrix_set_zero (self->cov_full);
  ncm_vector_set_zero (self->cov_packed);

  if (!g_key_file_has_key (snia_keyfile,
                           NC_DATA_SNIA_COV_DATA_GROUP,
                           NC_DATA_SNIA_COV_DATA_HAS_COMPLETE_COV_KEY,
                           &error))
    self->has_complete_cov = FALSE;
  else
    self->has_complete_cov = g_key_file_get_boolean (snia_keyfile,
                                                         NC_DATA_SNIA_COV_DATA_GROUP,
                                                         NC_DATA_SNIA_COV_DATA_HAS_COMPLETE_COV_KEY,
                                                         &error);

  if (!g_key_file_has_key (snia_keyfile,
                           NC_DATA_SNIA_COV_DATA_GROUP,
                           NC_DATA_SNIA_COV_DATA_KEY,
                           &error))
  {
    g_error ("nc_data_snia_cov_load: invalid %s key file, it must define at least the data file key ["NC_DATA_SNIA_COV_DATA_KEY "]", filename);
  }
  else
  {
    gchar *datafile = g_key_file_get_string (snia_keyfile,
                                             NC_DATA_SNIA_COV_DATA_GROUP,
                                             NC_DATA_SNIA_COV_DATA_KEY,
                                             &error);

    switch (self->cat_version)
    {
      case 0:
      case 1:
        _nc_data_snia_cov_load_snia_data_V01 (snia_cov, datafile);
        break;
      case 2:
        _nc_data_snia_cov_load_snia_data_V2 (snia_cov, datafile);
        break;
      default:
        g_assert_not_reached ();
        break;
    }
    g_free (datafile);
  }

  if (g_key_file_has_key (snia_keyfile,
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_DATA_DESC,
                          &error))
  {
    gchar *desc = g_key_file_get_string (snia_keyfile,
                                         NC_DATA_SNIA_COV_DATA_GROUP,
                                         NC_DATA_SNIA_COV_DATA_DESC,
                                         &error);

    ncm_data_set_desc (NCM_DATA (snia_cov), desc);
    g_free (desc);
  }

  /* Get magnitude-cut */
  if (g_key_file_has_key (snia_keyfile,
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_MAG_CUT,
                          &error))
  {
    const gdouble mag_cut = g_key_file_get_double (snia_keyfile,
                                                   NC_DATA_SNIA_COV_DATA_GROUP,
                                                   NC_DATA_SNIA_COV_MAG_CUT,
                                                   &error);

    nc_data_snia_cov_set_mag_cut (snia_cov, mag_cut);
  }
  else
  {
    nc_data_snia_cov_set_mag_cut (snia_cov, NC_DATA_SNIA_COV_MAG_CUT_DEFAULT); /* Sets the default value */
  }

  /* Get magnitude cov matrix */

  if (g_key_file_has_key (snia_keyfile,
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_MAG_KEY,
                          &error))
  {
    gchar *datafile = g_key_file_get_string (snia_keyfile,
                                             NC_DATA_SNIA_COV_DATA_GROUP,
                                             NC_DATA_SNIA_COV_MAG_KEY,
                                             &error);

    _nc_data_snia_cov_load_matrix (datafile, cov);
    ncm_matrix_add_constant (cov, 0.0);
    _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 0, 0);
    g_free (datafile);
  }

  /* Get width cov matrix */

  if (g_key_file_has_key (snia_keyfile,
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_WIDTH_KEY,
                          &error))
  {
    gchar *datafile = g_key_file_get_string (snia_keyfile,
                                             NC_DATA_SNIA_COV_DATA_GROUP,
                                             NC_DATA_SNIA_COV_WIDTH_KEY,
                                             &error);

    _nc_data_snia_cov_load_matrix (datafile, cov);
    ncm_matrix_add_constant (cov, 0.0);
    _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 1, 1);
    g_free (datafile);
  }

  /* Get colour cov matrix */

  if (g_key_file_has_key (snia_keyfile,
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_COLOUR_KEY,
                          &error))
  {
    gchar *datafile = g_key_file_get_string (snia_keyfile,
                                             NC_DATA_SNIA_COV_DATA_GROUP,
                                             NC_DATA_SNIA_COV_COLOUR_KEY,
                                             &error);

    _nc_data_snia_cov_load_matrix (datafile, cov);
    ncm_matrix_add_constant (cov, 0.0);
    _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 2, 2);

    g_free (datafile);
  }

  /* Get magnitude-width cov matrix */

  if (g_key_file_has_key (snia_keyfile,
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_MAG_WIDTH_KEY,
                          &error))
  {
    gchar *datafile = g_key_file_get_string (snia_keyfile,
                                             NC_DATA_SNIA_COV_DATA_GROUP,
                                             NC_DATA_SNIA_COV_MAG_WIDTH_KEY,
                                             &error);

    _nc_data_snia_cov_load_matrix (datafile, cov);
    ncm_matrix_add_constant (cov, 0.0);

    _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 0, 1);

    g_free (datafile);
  }

  /* Get magnitude-colour cov matrix */

  if (g_key_file_has_key (snia_keyfile,
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_MAG_COLOUR_KEY,
                          &error))
  {
    gchar *datafile = g_key_file_get_string (snia_keyfile,
                                             NC_DATA_SNIA_COV_DATA_GROUP,
                                             NC_DATA_SNIA_COV_MAG_COLOUR_KEY,
                                             &error);

    _nc_data_snia_cov_load_matrix (datafile, cov);
    ncm_matrix_add_constant (cov, 0.0);

    _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 0, 2);
    g_free (datafile);
  }

  /* Get width-colour cov matrix */

  if (g_key_file_has_key (snia_keyfile,
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_WIDTH_COLOUR_KEY,
                          &error))
  {
    gchar *datafile = g_key_file_get_string (snia_keyfile,
                                             NC_DATA_SNIA_COV_DATA_GROUP,
                                             NC_DATA_SNIA_COV_WIDTH_COLOUR_KEY,
                                             &error);

    _nc_data_snia_cov_load_matrix (datafile, cov);
    ncm_matrix_add_constant (cov, 0.0);

    _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 1, 2);
    g_free (datafile);
  }

  /* Get mbc-mbc cov matrix */

  if (g_key_file_has_key (snia_keyfile,
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_MBC_MBC_KEY,
                          &error))
  {
    gchar *datafile = g_key_file_get_string (snia_keyfile,
                                             NC_DATA_SNIA_COV_DATA_GROUP,
                                             NC_DATA_SNIA_COV_MBC_MBC_KEY,
                                             &error);

    _nc_data_snia_cov_load_matrix (datafile, self->cov_mbc_mbc);
    ncm_matrix_add_constant (self->cov_mbc_mbc, 0.0);

    g_free (datafile);
  }

  {
    const guint mu_len = self->mu_len;
    const guint tmu_len = 3 * mu_len;
    glong i, j;

    for (i = 0; i < tmu_len; i++)
    {
      for (j = i + 1; j < tmu_len; j++)
      {
        const gdouble cov_ij = ncm_matrix_get (self->cov_full, i, j);
        const gdouble cov_ji = ncm_matrix_get (self->cov_full, j, i);

        if (((i / mu_len) == (j / mu_len)) && (ncm_cmp (cov_ij, cov_ji, NC_DATA_SNIA_COV_SYMM_TOL, 0.0) != 0))
          g_error ("nc_data_snia_cov_load_txt: full covariance matrix is not symmetric in %ld %ld [% 22.15g != % 22.15g]\n",
                   i, j, cov_ij, cov_ji);

        ncm_matrix_set (self->cov_full, j, i, cov_ij);
      }
    }
  }

  nc_data_snia_cov_set_cov_full (snia_cov, self->cov_full);
  _NC_DATA_SNIA_COV_SET_DATA_INIT_ALL;
  ncm_matrix_clear (&cov);
  g_key_file_free (snia_keyfile);
}

static void
_nc_data_snia_cov_load_snia_data_V01 (NcDataSNIACov *snia_cov, const gchar *filename)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  GArray *dataset = self->dataset;
  gchar *line = NULL;
  gsize len = 0, tpos = 0;
  GError *error                = NULL;
  GIOChannel *file             = g_io_channel_new_file (filename, "r", &error);
  NcmVector *sigma_mag         = ncm_vector_new (self->mu_len);
  NcmVector *sigma_width       = ncm_vector_new (self->mu_len);
  NcmVector *sigma_colour      = ncm_vector_new (self->mu_len);
  NcmVector *diag_mag_width    = ncm_vector_new (self->mu_len);
  NcmVector *diag_mag_colour   = ncm_vector_new (self->mu_len);
  NcmVector *diag_width_colour = ncm_vector_new (self->mu_len);
  guint i;

  if (file == NULL)
    g_error ("_nc_data_snia_cov_load_snia_data_V01: cannot open file %s: %s",
             filename, error->message);

  {
    GRegex *comment_line = g_regex_new ("^\\s*#", 0, 0, &error);
    guint n              = 0;
    guint nrow           = 0;
    guint max_dset_id    = 0;
    guint min_dset_id    = 1000;
    gboolean has_dataset = FALSE;

    while (g_io_channel_read_line (file, &line, &len, &tpos, &error) != G_IO_STATUS_EOF)
    {
      if (g_regex_match (comment_line, line, 0, NULL))
      {
        g_free (line);
        continue;
      }
      else
      {
        gchar **itens = g_regex_split_simple ("\\s+", line, 0, 0);

        if (n == 0)
        {
          guint pad_end = 0;

          n = g_strv_length (itens);

          if (*itens[n - 1] == '\0')
            pad_end++;

          if (n - 1 - pad_end == NC_DATA_SNIA_COV_V1_LENGTH)
          {
            has_dataset = FALSE;
          }
          else if (n - 2 - pad_end >= NC_DATA_SNIA_COV_V1_LENGTH)
          {
            has_dataset = TRUE;

            if (n - 2 - pad_end > NC_DATA_SNIA_COV_V1_LENGTH)
              g_warning ("_nc_data_snia_cov_load_snia_data_V01: data file [%s] must have %d or %d columns, it has %u, ignoring the extra columns!",
                         filename, NC_DATA_SNIA_COV_V1_LENGTH + 1, NC_DATA_SNIA_COV_V1_LENGTH + 2, n - pad_end);
          }
          else
          {
            g_error ("_nc_data_snia_cov_load_snia_data_V01: data file [%s] must have %d or %d columns, it has %u!",
                     filename, NC_DATA_SNIA_COV_V1_LENGTH + 1, NC_DATA_SNIA_COV_V1_LENGTH + 2, n - pad_end);
          }
        }
        else if (n != g_strv_length (itens))
        {
          g_error ("_nc_data_snia_cov_load_snia_data_V01: data file [%s] has different number of columns in different rows [%u]!", filename, nrow);
        }

        if (nrow >= self->mu_len)
          g_error ("_nc_data_snia_cov_load_snia_data_V01: cannot load data file [%s] expected nrows %u obtained >%u!",
                   filename, self->mu_len, nrow);

        {
          NcmVector *data[NC_DATA_SNIA_COV_V1_LENGTH];

          data[NC_DATA_SNIA_COV_ZCMB]              = self->z_cmb;
          data[NC_DATA_SNIA_COV_ZHE]               = self->z_he;
          data[NC_DATA_SNIA_COV_SIGMA_Z]           = self->sigma_z;
          data[NC_DATA_SNIA_COV_MAG]               = self->mag;
          data[NC_DATA_SNIA_COV_SIGMA_MAG]         = sigma_mag;
          data[NC_DATA_SNIA_COV_WIDTH]             = self->width;
          data[NC_DATA_SNIA_COV_SIGMA_WIDTH]       = sigma_width;
          data[NC_DATA_SNIA_COV_COLOUR]            = self->colour;
          data[NC_DATA_SNIA_COV_SIGMA_COLOUR]      = sigma_colour;
          data[NC_DATA_SNIA_COV_THIRDPAR]          = self->thirdpar;
          data[NC_DATA_SNIA_COV_SIGMA_THIRDPAR]    = self->sigma_thirdpar;
          data[NC_DATA_SNIA_COV_DIAG_MAG_WIDTH]    = diag_mag_width;
          data[NC_DATA_SNIA_COV_DIAG_MAG_COLOUR]   = diag_mag_colour;
          data[NC_DATA_SNIA_COV_DIAG_WIDTH_COLOUR] = diag_width_colour;

          for (i = 0; i < NC_DATA_SNIA_COV_V1_LENGTH; i++)
          {
            const gdouble val = g_ascii_strtod (itens[i + 1], NULL);

            ncm_vector_set (data[i], nrow, val);
          }
        }

        if (has_dataset)
        {
          gint64 dset_id = g_ascii_strtoll (itens[i + 1], NULL, 10);

          g_array_index (dataset, guint32, nrow) = dset_id;
          max_dset_id                            = GSL_MAX (max_dset_id, dset_id);
          min_dset_id                            = GSL_MIN (min_dset_id, dset_id);
        }
        else
        {
          g_array_index (dataset, guint32, nrow) = 0;
        }

        nrow++;
        g_strfreev (itens);
        g_free (line);
      }
    }

    if (nrow != self->mu_len)
      g_error ("_nc_data_snia_cov_load_snia_data_V01: cannot load data file [%s] expected nrows %u obtained %u",
               filename, self->mu_len, nrow);

    if (min_dset_id > 1)
    {
      g_error ("_nc_data_snia_cov_load_snia_data_V01: wrongly aligned dataset ids. First id = %u\n", min_dset_id);
    }
    else if (min_dset_id == 1)
    {
      for (i = 0; i < dataset->len; i++)
      {
        g_array_index (dataset, guint32, i)--;
      }

      max_dset_id--;
    }

    nc_data_snia_cov_set_abs_mag_set (snia_cov, self->dataset);

    g_regex_unref (comment_line);
  }

  _nc_data_snia_cov_diag_to_full_cov (snia_cov, sigma_mag, sigma_width, sigma_colour, diag_mag_width, diag_mag_colour, diag_width_colour);

  ncm_vector_free (sigma_mag);
  ncm_vector_free (sigma_width);
  ncm_vector_free (sigma_colour);
  ncm_vector_free (diag_mag_width);
  ncm_vector_free (diag_mag_colour);
  ncm_vector_free (diag_width_colour);

  g_io_channel_unref (file);
}

static void
_nc_data_snia_cov_load_snia_data_V2 (NcDataSNIACov *snia_cov, const gchar *filename)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  gchar *line = NULL;
  gsize len = 0, tpos = 0;
  GError *error                = NULL;
  GIOChannel *file             = g_io_channel_new_file (filename, "r", &error);
  guint i;

  if (file == NULL)
    g_error ("_nc_data_snia_cov_load_snia_data_V2: cannot open file %s: %s",
             filename, error->message);

  {
    GRegex *comment_line = g_regex_new ("^\\s*#", 0, 0, &error);
    guint n              = 0;
    guint nrow           = 0;

    while (g_io_channel_read_line (file, &line, &len, &tpos, &error) != G_IO_STATUS_EOF)
    {
      if (g_regex_match (comment_line, line, 0, NULL))
      {
        g_free (line);
        continue;
      }
      else
      {
        gchar **itens = g_regex_split_simple ("\\s+", line, 0, 0);

        if (n == 0)
        {
          guint pad_end = 0;

          n = g_strv_length (itens);

          if (*itens[n - 1] == '\0')
            pad_end++;

          if (n - pad_end != NC_DATA_SNIA_COV_V2_LENGTH)
          {
            g_error ("_nc_data_snia_cov_load_snia_data_V2: data file [%s] must have %d columns, it has %u!",
                     filename, NC_DATA_SNIA_COV_V2_LENGTH, n - pad_end);
          }
        }
        else if (n != g_strv_length (itens))
        {
          g_error ("_nc_data_snia_cov_load_snia_data_V2: data file [%s] has different number of columns in different rows [%u]!", filename, nrow);
        }

        if (nrow >= self->mu_len)
          g_error ("_nc_data_snia_cov_load_snia_data_V2: cannot load data file [%s] expected %u rows but obtained >%u!",
                   filename, self->mu_len, nrow);

        {
          gpointer data[NC_DATA_SNIA_COV_V2_LENGTH];

          data[NC_DATA_SNIA_COV_V2_ZHD]              = self->z_hd;
          data[NC_DATA_SNIA_COV_V2_SIGMA_ZHD]        = self->sigma_z_hd;
          data[NC_DATA_SNIA_COV_V2_ZCMB]             = self->z_cmb;
          data[NC_DATA_SNIA_COV_V2_SIGMA_ZCMB]       = self->sigma_z_cmb;
          data[NC_DATA_SNIA_COV_V2_ZHE]              = self->z_he;
          data[NC_DATA_SNIA_COV_V2_SIGMA_ZHE]        = self->sigma_z_he;
          data[NC_DATA_SNIA_COV_V2_MAG_B_CORR]       = self->mag_b_corr;
          data[NC_DATA_SNIA_COV_V2_SIGMA_MAG_B_CORR] = self->sigma_mag_b_corr;
          data[NC_DATA_SNIA_COV_V2_MU_CORR]          = self->mu_corr;
          data[NC_DATA_SNIA_COV_V2_SIGMA_MU_CORR]    = self->sigma_mu_corr;
          data[NC_DATA_SNIA_COV_V2_CEPH_DIST]        = self->ceph_dist;
          data[NC_DATA_SNIA_COV_V2_IS_CALIB]         = self->is_calib;
          data[NC_DATA_SNIA_COV_V2_USED_IN_SH0ES]    = self->used_in_sh0es;
          data[NC_DATA_SNIA_COV_V2_COLOUR]           = self->colour;
          data[NC_DATA_SNIA_COV_V2_SIGMA_COLOUR]     = self->sigma_colour;
          data[NC_DATA_SNIA_COV_V2_WIDTH]            = self->width;
          data[NC_DATA_SNIA_COV_V2_SIGMA_WIDTH]      = self->sigma_width;
          data[NC_DATA_SNIA_COV_V2_MAG]              = self->mag;
          data[NC_DATA_SNIA_COV_V2_SIGMA_MAG]        = self->sigma_mag;
          data[NC_DATA_SNIA_COV_V2_AMPL]             = self->ampl;
          data[NC_DATA_SNIA_COV_V2_SIGMA_AMPL]       = self->sigma_ampl;
          data[NC_DATA_SNIA_COV_V2_WIDTH_COLOUR]     = self->cov_width_colour;
          data[NC_DATA_SNIA_COV_V2_WIDTH_AMPL]       = self->cov_width_ampl;
          data[NC_DATA_SNIA_COV_V2_COLOUR_AMPL]      = self->cov_colour_ampl;
          data[NC_DATA_SNIA_COV_V2_RA]               = self->ra;
          data[NC_DATA_SNIA_COV_V2_DEC]              = self->dec;
          data[NC_DATA_SNIA_COV_V2_HOST_RA]          = self->host_ra;
          data[NC_DATA_SNIA_COV_V2_HOST_DEC]         = self->host_dec;
          data[NC_DATA_SNIA_COV_V2_HOST_ANG_SEP]     = self->host_ang_sep;
          data[NC_DATA_SNIA_COV_V2_VPEC]             = self->vpec;
          data[NC_DATA_SNIA_COV_V2_SIGMA_VPEC]       = self->sigma_vpec;

          for (i = 0; i < NC_DATA_SNIA_COV_V2_LENGTH; i++)
          {
            const gdouble val = g_ascii_strtod (itens[i], NULL);

            if ((i == NC_DATA_SNIA_COV_V2_IS_CALIB) || (i == NC_DATA_SNIA_COV_V2_USED_IN_SH0ES))
            {
              gint64 tmp1 = g_ascii_strtoll (itens[i], NULL, 10);
              GArray *a = data[i];

              g_array_index (a, guint32, nrow) = tmp1;
            }
            else
            {
              ncm_vector_set (data[i], nrow, val);
              ncm_vector_add_constant (data[i], 0.0);
            }
          }
        }

        nrow++;
        g_strfreev (itens);
        g_free (line);
      }
    }

    if (nrow != self->mu_len)
      g_error ("_nc_data_snia_cov_load_snia_data_V2: cannot load data file [%s] expected %u rows obtained %u",
               filename, self->mu_len, nrow);

    g_regex_unref (comment_line);
  }

  g_io_channel_unref (file);
}


static void
_nc_data_snia_cov_load_matrix (const gchar *filename, NcmMatrix *data)
{
  GError *error = NULL;
  gchar *file   = NULL;
  gsize length  = 0;
  gchar **itens;
  guint itens_len = 0;
  guint pad_start = 0;
  guint pad_end = 0;
  guint64 matrix_len = 0;
  guint i, j;

  if (!g_file_get_contents (filename, &file, &length, &error))
    g_error ("_nc_data_snia_cov_load_matrix: cannot open file %s: %s",
             filename, error->message);

  itens     = g_regex_split_simple ("\\s+", file, 0, 0);
  itens_len = g_strv_length (itens);

  if (*itens[0] == '\0')
    pad_start = 1;

  if (*itens[itens_len - 1] == '\0')
    pad_end = 1;

  if (pad_start)
    itens_len--;

  if (pad_end)
    itens_len--;

  matrix_len = g_ascii_strtoll (itens[pad_start], NULL, 10);
  pad_start++;
  itens_len--;

  if (matrix_len != ncm_matrix_nrows (data))
    g_error ("_nc_data_snia_cov_load_matrix: expected a %ux%u matrix but got %"G_GINT64_MODIFIER "dx%"G_GINT64_MODIFIER "d",
             ncm_matrix_nrows (data), ncm_matrix_nrows (data), matrix_len, matrix_len);

  if (matrix_len * matrix_len != itens_len)
    g_error ("_nc_data_snia_cov_load_matrix: matrix header say %"G_GINT64_MODIFIER "d length but got %u",
             matrix_len * matrix_len, itens_len);

  for (i = 0; i < matrix_len; i++)
  {
    for (j = 0; j < matrix_len; j++)
    {
      const gdouble val = g_ascii_strtod (itens[pad_start + i * matrix_len + j], NULL);

      ncm_matrix_set (data, i, j, val);
    }
  }

  g_strfreev (itens);
  g_free (file);
}

static void
_nc_data_snia_cov_diag_to_full_cov (NcDataSNIACov *snia_cov,
                                    NcmVector     *sigma_mag,
                                    NcmVector     *sigma_width,
                                    NcmVector     *sigma_colour,
                                    NcmVector     *diag_mag_width,
                                    NcmVector     *diag_mag_colour,
                                    NcmVector     *diag_width_colour)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  const guint mu_len = self->mu_len;
  guint i;

  for (i = 0; i < mu_len; i++)
  {
    const gdouble sigma_mag_i    = ncm_vector_get (sigma_mag, i);
    const gdouble sigma_width_i  = ncm_vector_get (sigma_width, i);
    const gdouble sigma_colour_i = ncm_vector_get (sigma_colour, i);

    const gdouble mag_mag_i       = sigma_mag_i * sigma_mag_i;
    const gdouble mag_width_i     = ncm_vector_get (diag_mag_width, i);
    const gdouble mag_colour_i    = ncm_vector_get (diag_mag_colour, i);
    const gdouble width_width_i   = sigma_width_i * sigma_width_i;
    const gdouble width_colour_i  = ncm_vector_get (diag_width_colour, i);
    const gdouble colour_colour_i = sigma_colour_i * sigma_colour_i;

    ncm_matrix_set (self->cov_full, 0 * mu_len + i, 0 * mu_len + i, mag_mag_i);
    ncm_matrix_set (self->cov_full, 0 * mu_len + i, 1 * mu_len + i, mag_width_i);
    ncm_matrix_set (self->cov_full, 0 * mu_len + i, 2 * mu_len + i, mag_colour_i);
    ncm_matrix_set (self->cov_full, 1 * mu_len + i, 1 * mu_len + i, width_width_i);
    ncm_matrix_set (self->cov_full, 1 * mu_len + i, 2 * mu_len + i, width_colour_i);
    ncm_matrix_set (self->cov_full, 2 * mu_len + i, 2 * mu_len + i, colour_colour_i);
  }
}

static void
_nc_data_snia_cov_matrix_to_cov_full (NcDataSNIACov *snia_cov, NcmMatrix *cov, guint i, guint j)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  const guint mu_len = self->mu_len;
  NcmMatrix *subcov  = ncm_matrix_get_submatrix (self->cov_full, i * mu_len, j * mu_len, mu_len, mu_len);

  ncm_matrix_add_mul (subcov, 1.0, cov);
  ncm_matrix_clear (&subcov);
}

#ifdef NUMCOSMO_HAVE_CFITSIO

static void
nc_data_snia_cov_load_V0 (NcDataSNIACov *snia_cov, fitsfile *fptr)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  NcmVector *sigma_mag         = ncm_vector_new (self->mu_len);
  NcmVector *sigma_width       = ncm_vector_new (self->mu_len);
  NcmVector *sigma_colour      = ncm_vector_new (self->mu_len);
  NcmVector *diag_mag_width    = ncm_vector_new (self->mu_len);
  NcmVector *diag_mag_colour   = ncm_vector_new (self->mu_len);
  NcmVector *diag_width_colour = ncm_vector_new (self->mu_len);
  NcmVector *data[NC_DATA_SNIA_COV_ABSMAG_SET];
  gint status = 0;
  gint i;

  data[NC_DATA_SNIA_COV_ZCMB]              = self->z_cmb;
  data[NC_DATA_SNIA_COV_ZHE]               = self->z_he;
  data[NC_DATA_SNIA_COV_SIGMA_Z]           = self->sigma_z;
  data[NC_DATA_SNIA_COV_MAG]               = self->mag;
  data[NC_DATA_SNIA_COV_SIGMA_MAG]         = sigma_mag;
  data[NC_DATA_SNIA_COV_WIDTH]             = self->width;
  data[NC_DATA_SNIA_COV_SIGMA_WIDTH]       = sigma_width;
  data[NC_DATA_SNIA_COV_COLOUR]            = self->colour;
  data[NC_DATA_SNIA_COV_SIGMA_COLOUR]      = sigma_colour;
  data[NC_DATA_SNIA_COV_THIRDPAR]          = self->thirdpar;
  data[NC_DATA_SNIA_COV_SIGMA_THIRDPAR]    = self->sigma_thirdpar;
  data[NC_DATA_SNIA_COV_DIAG_MAG_WIDTH]    = diag_mag_width;
  data[NC_DATA_SNIA_COV_DIAG_MAG_COLOUR]   = diag_mag_colour;
  data[NC_DATA_SNIA_COV_DIAG_WIDTH_COLOUR] = diag_width_colour;

  for (i = 0; i < NC_DATA_SNIA_COV_ABSMAG_SET; i++)
  {
    fits_read_col_dbl (fptr, i + 1, 1, 1, self->mu_len, GSL_NAN,
                       ncm_vector_ptr (data[i], 0), NULL,
                       &status);
    NCM_FITS_ERROR (status);
  }

  fits_read_col_uint (fptr, NC_DATA_SNIA_COV_ABSMAG_SET + 1, 1, 1,
                      self->mu_len,
                      0, &g_array_index (self->dataset, guint32, 0),
                      NULL, &status);
  NCM_FITS_ERROR (status);

  nc_data_snia_cov_set_abs_mag_set (snia_cov, self->dataset);

  _nc_data_snia_cov_diag_to_full_cov (snia_cov, sigma_mag, sigma_width, sigma_colour, diag_mag_width, diag_mag_colour, diag_width_colour);

  {
    NcmMatrix *cov = ncm_matrix_new (self->mu_len, self->mu_len);

    for (i = NC_DATA_SNIA_COV_VAR_MAG; i < NC_DATA_SNIA_COV_TOTAL_LENGTH; i++)
    {
      fits_read_col_dbl (fptr, i + 1, 1, 1, self->mu_len * self->mu_len,
                         GSL_NAN, ncm_matrix_ptr (cov, 0, 0),
                         NULL, &status);
      NCM_FITS_ERROR (status);

      switch (i)
      {
        case NC_DATA_SNIA_COV_VAR_MAG:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 0, 0);
          break;
        case NC_DATA_SNIA_COV_VAR_MAG_WIDTH:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 0, 1);
          break;
        case NC_DATA_SNIA_COV_VAR_MAG_COLOUR:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 0, 2);
          break;
        case NC_DATA_SNIA_COV_VAR_WIDTH:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 1, 1);
          break;
        case NC_DATA_SNIA_COV_VAR_WIDTH_COLOUR:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 1, 2);
          break;
        case NC_DATA_SNIA_COV_VAR_COLOUR:
          _nc_data_snia_cov_matrix_to_cov_full (snia_cov, cov, 2, 2);
          break;
        default:
          g_assert_not_reached ();
          break;
      }
    }

    ncm_matrix_free (cov);
  }

  ncm_vector_free (sigma_mag);
  ncm_vector_free (sigma_width);
  ncm_vector_free (sigma_colour);
  ncm_vector_free (diag_mag_width);
  ncm_vector_free (diag_mag_colour);
  ncm_vector_free (diag_width_colour);
}

/**
 * nc_data_snia_cov_save:
 * @snia_cov: a #NcDataSNIACov
 * @filename: file name of the catalog
 * @overwrite: whether to overwrite an already existing catalog
 *
 * Saves the catalog in fits (cfitsio) format using @filename.
 *
 */
void
nc_data_snia_cov_save (NcDataSNIACov *snia_cov, const gchar *filename, gboolean overwrite)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  /*******************************************************************/
  /* Create a binary table extension                                 */
  /*******************************************************************/
  fitsfile *fptr; /* pointer to the FITS file, defined in fitsio.h */
  gint status;
  gint tfields;

  gchar extname[]        = "NcDataSNIACov"; /* extension name */
  GPtrArray *ttype_array = g_ptr_array_sized_new (NC_DATA_SNIA_COV_V2_TOTAL_LENGTH);
  GPtrArray *tform_array = g_ptr_array_sized_new (NC_DATA_SNIA_COV_V2_TOTAL_LENGTH);
  GPtrArray *tunit_array = g_ptr_array_sized_new (NC_DATA_SNIA_COV_V2_TOTAL_LENGTH);

  if (ncm_data_gauss_cov_get_size (NCM_DATA_GAUSS_COV (snia_cov)) == 0)
    g_error ("nc_data_snia_cov_save: cannot save an empty catalog.");

  g_ptr_array_set_size (ttype_array, NC_DATA_SNIA_COV_V2_TOTAL_LENGTH);
  g_ptr_array_set_size (tform_array, NC_DATA_SNIA_COV_V2_TOTAL_LENGTH);
  g_ptr_array_set_size (tunit_array, NC_DATA_SNIA_COV_V2_TOTAL_LENGTH);

  g_ptr_array_set_free_func (tform_array, g_free);

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_ZHD) = "Z_HD";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_ZHD) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_ZHD) = "Hubble Diagram redshift";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_SIGMA_ZHD) = "SIGMA_Z_HD";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_SIGMA_ZHD) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_SIGMA_ZHD) = "Hubble Diagram redshift - sigma";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_ZCMB) = "Z_CMB";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_ZCMB) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_ZCMB) = "CMB corrected redshift";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_SIGMA_ZCMB) = "SIGMA_Z_CMB";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_SIGMA_ZCMB) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_SIGMA_ZCMB) = "CMB corrected redshift - sigma";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_ZHE) = "Z_HE";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_ZHE) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_ZHE) = "Heliocentric redshift";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_SIGMA_ZHE) = "SIGMA_ZHE";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_SIGMA_ZHE) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_SIGMA_ZHE) = "Heliocentric redshift - sigma";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_MAG_B_CORR) = "MAG_B_CORR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_MAG_B_CORR) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_MAG_B_CORR) = "Magnitude B corr";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_SIGMA_MAG_B_CORR) = "SIGMA_MAG_B_CORR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_SIGMA_MAG_B_CORR) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_SIGMA_MAG_B_CORR) = "Magnitude B corr - sigma";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_MU_CORR) = "MU_CORR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_MU_CORR) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_MU_CORR) = "Distance modulus corr";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_SIGMA_MU_CORR) = "SIGMA_MU_CORR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_SIGMA_MU_CORR) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_SIGMA_MU_CORR) = "Distance modulus corr - sigma";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_CEPH_DIST) = "CEPH_DIST";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_CEPH_DIST) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_CEPH_DIST) = "Cepheid absolute distance to host";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_IS_CALIB) = "IS_CALIB";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_IS_CALIB) = g_strdup ("1J");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_IS_CALIB) = "Whether to use as distance calibrator";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_USED_IN_SH0ES) = "USED_IN_SH0ES";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_USED_IN_SH0ES) = g_strdup ("1J");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_USED_IN_SH0ES) = "Whether it was used in SH0ES analysis";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_COLOUR) = "COLOUR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_COLOUR) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_COLOUR) = "COLOUR";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_SIGMA_COLOUR) = "SIGMA_COLOUR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_SIGMA_COLOUR) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_SIGMA_COLOUR) = "SIGMA_COLOUR";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_WIDTH) = "WIDTH";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_WIDTH) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_WIDTH) = "WIDTH";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_SIGMA_WIDTH) = "SIGMA_WIDTH";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_SIGMA_WIDTH) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_SIGMA_WIDTH) = "SIGMA_WIDTH";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_MAG) = "MAG_B";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_MAG) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_MAG) = "MAG_B";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_SIGMA_MAG) = "SIGMA_MAG_B";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_SIGMA_MAG) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_SIGMA_MAG) = "MAG_B SIGMA";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_AMPL) = "AMPL";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_AMPL) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_AMPL) = "AMPL";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_SIGMA_AMPL) = "SIGMA_AMPL";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_SIGMA_AMPL) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_SIGMA_AMPL) = "AMPL SIGMA";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_WIDTH_COLOUR) = "WIDTH_COLOUR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_WIDTH_COLOUR) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_WIDTH_COLOUR) = "WIDTH_COLOUR";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_WIDTH_AMPL) = "WIDTH_AMPL";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_WIDTH_AMPL) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_WIDTH_AMPL) = "WIDTH_AMPL";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_COLOUR_AMPL) = "COLOUR_AMPL";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_COLOUR_AMPL) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_COLOUR_AMPL) = "COLOUR_AMPL";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_RA) = "RA";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_RA) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_RA) = "RA";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_DEC) = "DEC";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_DEC) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_DEC) = "DEC";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_HOST_RA) = "HOST_RA";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_HOST_RA) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_HOST_RA) = "HOST_RA";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_HOST_DEC) = "HOST_DEC";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_HOST_DEC) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_HOST_DEC) = "HOST_DEC";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_HOST_ANG_SEP) = "ANG_SEP";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_HOST_ANG_SEP) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_HOST_ANG_SEP) = "ANG_SEP";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_VPEC) = "VPEC";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_VPEC) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_VPEC) = "VPEC";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_SIGMA_VPEC) = "SIGMA_VPEC";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_SIGMA_VPEC) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_SIGMA_VPEC) = "SIGMA_VPEC";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_V2_MAG_B_CORR_MAG_B_CORR) = "COV_MAG_B_CORR_MAG_B_CORR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_V2_MAG_B_CORR_MAG_B_CORR) = g_strdup_printf ("%uD", self->mu_len);
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_V2_MAG_B_CORR_MAG_B_CORR) = "MAG_B_CORR-MAG_B_CORR";

  tfields = ttype_array->len;

  /* initialize status before calling fitsio routines */
  status = 0;

  if (overwrite && g_file_test (filename, G_FILE_TEST_EXISTS))
    g_unlink (filename);

  /* create new FITS file */
  fits_create_file (&fptr, filename, &status);
  NCM_FITS_ERROR (status);

  /* append a new empty binary table onto the FITS file */
  fits_create_tbl (fptr, BINARY_TBL, self->mu_len, tfields, (gchar **) ttype_array->pdata, (gchar **) tform_array->pdata,
                   (gchar **) tunit_array->pdata, extname, &status);
  NCM_FITS_ERROR (status);

  {
    const gchar *desc = ncm_data_peek_desc (NCM_DATA (snia_cov));

    fits_write_key_longstr (fptr, NC_DATA_SNIA_COV_CAT_DESC, desc, NC_DATA_SNIA_COV_CAT_DESC_COMMENT, &status);
    NCM_FITS_ERROR (status);
  }

  {
    /* Always create catalogs of latest version */
    fits_write_key_lng (fptr, NC_DATA_SNIA_COV_CAT_VERSION, NC_DATA_SNIA_COV_CAT_LAST_VERSION, NC_DATA_SNIA_COV_CAT_VERSION_COMMENT, &status);
    NCM_FITS_ERROR (status);
  }

  {
    /* Set to true if there is a complete covariance matrix */
    fits_write_key_log (fptr, NC_DATA_SNIA_COV_CAT_HAS_COMPLETE_COV, self->has_complete_cov, NC_DATA_SNIA_COV_CAT_HAS_COMPLETE_COV_COMMENT, &status);
    NCM_FITS_ERROR (status);
  }

  {
    NcmVector *data[NC_DATA_SNIA_COV_V2_LENGTH] = {NULL, };
    guint i;

    data[NC_DATA_SNIA_COV_V2_ZHD]              = self->z_hd;
    data[NC_DATA_SNIA_COV_V2_SIGMA_ZHD]        = self->sigma_z_hd;
    data[NC_DATA_SNIA_COV_V2_ZCMB]             = self->z_cmb;
    data[NC_DATA_SNIA_COV_V2_SIGMA_ZCMB]       = self->sigma_z_cmb;
    data[NC_DATA_SNIA_COV_V2_ZHE]              = self->z_he;
    data[NC_DATA_SNIA_COV_V2_SIGMA_ZHE]        = self->sigma_z_he;
    data[NC_DATA_SNIA_COV_V2_MAG_B_CORR]       = self->mag_b_corr;
    data[NC_DATA_SNIA_COV_V2_SIGMA_MAG_B_CORR] = self->sigma_mag_b_corr;
    data[NC_DATA_SNIA_COV_V2_MU_CORR]          = self->mu_corr;
    data[NC_DATA_SNIA_COV_V2_SIGMA_MU_CORR]    = self->sigma_mu_corr;
    data[NC_DATA_SNIA_COV_V2_CEPH_DIST]        = self->ceph_dist;
    data[NC_DATA_SNIA_COV_V2_COLOUR]           = self->colour;
    data[NC_DATA_SNIA_COV_V2_SIGMA_COLOUR]     = self->sigma_colour;
    data[NC_DATA_SNIA_COV_V2_WIDTH]            = self->width;
    data[NC_DATA_SNIA_COV_V2_SIGMA_WIDTH]      = self->sigma_width;
    data[NC_DATA_SNIA_COV_V2_MAG]            = self->mag;
    data[NC_DATA_SNIA_COV_V2_SIGMA_MAG]      = self->sigma_mag;
    data[NC_DATA_SNIA_COV_V2_AMPL]             = self->ampl;
    data[NC_DATA_SNIA_COV_V2_SIGMA_AMPL]       = self->sigma_ampl;
    data[NC_DATA_SNIA_COV_V2_WIDTH_COLOUR]     = self->cov_width_colour;
    data[NC_DATA_SNIA_COV_V2_WIDTH_AMPL]       = self->cov_width_ampl;
    data[NC_DATA_SNIA_COV_V2_COLOUR_AMPL]      = self->cov_colour_ampl;
    data[NC_DATA_SNIA_COV_V2_RA]               = self->ra;
    data[NC_DATA_SNIA_COV_V2_DEC]              = self->dec;
    data[NC_DATA_SNIA_COV_V2_HOST_RA]          = self->host_ra;
    data[NC_DATA_SNIA_COV_V2_HOST_DEC]         = self->host_dec;
    data[NC_DATA_SNIA_COV_V2_HOST_ANG_SEP]     = self->host_ang_sep;
    data[NC_DATA_SNIA_COV_V2_VPEC]             = self->vpec;
    data[NC_DATA_SNIA_COV_V2_SIGMA_VPEC]       = self->sigma_vpec;

    for (i = 0; i < NC_DATA_SNIA_COV_V2_LENGTH; i++)
    {
      if ((i == NC_DATA_SNIA_COV_V2_IS_CALIB) || (i == NC_DATA_SNIA_COV_V2_USED_IN_SH0ES))
        continue;

      g_assert (data[i] != NULL);
      g_assert (NCM_IS_VECTOR (data[i]));

      fits_write_col_dbl (fptr, i + 1, 1, 1, self->mu_len, ncm_vector_ptr (data[i], 0), &status);
      NCM_FITS_ERROR (status);
    }

    fits_write_col_uint (fptr, NC_DATA_SNIA_COV_V2_IS_CALIB + 1, 1, 1,
                         self->mu_len, &g_array_index (self->is_calib, guint32, 0), &status);
    NCM_FITS_ERROR (status);

    fits_write_col_uint (fptr, NC_DATA_SNIA_COV_V2_USED_IN_SH0ES + 1, 1, 1,
                         self->mu_len, &g_array_index (self->used_in_sh0es, guint32, 0), &status);
    NCM_FITS_ERROR (status);

    fits_write_col_dbl (fptr, NC_DATA_SNIA_COV_V2_MAG_B_CORR_MAG_B_CORR + 1, 1, 1, self->mu_len * self->mu_len,
                        ncm_matrix_ptr (self->cov_mbc_mbc, 0, 0), &status);
    NCM_FITS_ERROR (status);
  }

  fits_write_chksum (fptr, &status);
  NCM_FITS_ERROR (status);

  fits_close_file (fptr, &status);
  NCM_FITS_ERROR (status);

  g_ptr_array_unref (ttype_array);
  g_ptr_array_unref (tform_array);
  g_ptr_array_unref (tunit_array);

  return;
}

#endif /* NUMCOSMO_HAVE_CFITSIO */

static void
_nc_data_snia_cov_restore_full_cov (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  const guint mu_len = self->mu_len;
  const guint tmu_len = 3 * mu_len;
  guint i, j;

  for (i = 0; i < tmu_len; i++)
  {
    const gdouble diag = ncm_vector_get (self->cov_full_diag, i);

    ncm_matrix_set (self->cov_full, i, i, diag);

    for (j = i + 1; j < tmu_len; j++)
    {
      ncm_matrix_set (self->cov_full, i, j, ncm_matrix_get (self->cov_full, j, i));
    }
  }
}

static void
_nc_data_snia_cov_prep_to_resample (NcDataSNIACov *snia_cov, NcSNIADistCov *dcov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  const guint mu_len = self->mu_len;
  const guint tmu_len = 3 * self->mu_len;
  guint i, j;
  gint ret;

  if ((mu_len == 0) || !self->has_complete_cov)
    g_error ("_nc_data_snia_cov_prep_to_resample: cannot prepare to resample, empty catalog %d or it hasn't a complete covariance %d.\n",
             mu_len == 0, !self->has_complete_cov);

  for (i = 0; i < tmu_len; i++)
  {
    const gdouble diag = ncm_vector_get (self->cov_full_diag, i);

    ncm_matrix_set (self->cov_full, i, i, diag);

    for (j = i + 1; j < tmu_len; j++)
    {
      ncm_matrix_set (self->cov_full, i, j, ncm_matrix_get (self->cov_full, j, i));
    }
  }

  ret = ncm_matrix_cholesky_decomp (self->cov_full, 'U');

  if (ret != 0)
    g_error ("_nc_data_snia_cov_prep_to_resample[ncm_matrix_cholesky_decomp]: %d.", ret);

  self->cov_full_state = NC_DATA_SNIA_COV_PREP_TO_RESAMPLE;
}

static void
_nc_data_snia_cov_prep_to_estimate (NcDataSNIACov *snia_cov, NcSNIADistCov *dcov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  const guint mu_len = self->mu_len;
  const guint tmu_len = 3 * self->mu_len;
  guint i, j;
  gint ret;

  if ((mu_len == 0) || !self->has_complete_cov)
    g_error ("_nc_data_snia_cov_prep_to_estimate: cannot prepare to estimate, empty catalog %d or it hasn't a complete covariance %d.\n",
             mu_len == 0, !self->has_complete_cov);

  for (i = 0; i < mu_len; i++)
  {
    const gdouble diag = ncm_vector_get (self->cov_full_diag, i) +
                         nc_snia_dist_cov_extra_var (dcov, snia_cov, i);

    ncm_matrix_set (self->cov_full, i, i, diag);

    for (j = i + 1; j < tmu_len; j++)
    {
      ncm_matrix_set (self->cov_full, i, j, ncm_matrix_get (self->cov_full, j, i));
    }
  }

  /* Continue for the last 2 * mu_len. */
  for (i = mu_len; i < tmu_len; i++)
  {
    const gdouble diag = ncm_vector_get (self->cov_full_diag, i);

    ncm_matrix_set (self->cov_full, i, i, diag);

    for (j = i + 1; j < tmu_len; j++)
    {
      ncm_matrix_set (self->cov_full, i, j, ncm_matrix_get (self->cov_full, j, i));
    }
  }

  /* Make the Cholesky decomposition substituting the upper triagle of cov_full. */
  ret = ncm_matrix_cholesky_decomp (self->cov_full, 'U');

  if (ret != 0)
    g_error ("_nc_data_snia_cov_prep_to_estimate[ncm_matrix_cholesky_decomp]: %d.", ret);

  self->cov_full_state = NC_DATA_SNIA_COV_PREP_TO_ESTIMATE;
}

static void
_nc_data_snia_cov_save_cov_lowertri (NcDataSNIACov *snia_cov)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  const guint mu_len = self->mu_len;
  const guint tmu_len = 3 * self->mu_len;
  guint i, j;

  if (self->cov_full_state != NC_DATA_SNIA_COV_PREP_TO_NOTHING)
    g_error ("_nc_data_snia_cov_save_cov_lowertri: cannot save since Cholesky decomposition was already done.");

  for (i = 0; i < mu_len; i++)
  {
    /* Save the original diagonal terms. */
    ncm_vector_set (self->cov_full_diag, i, ncm_matrix_get (self->cov_full, i, i));

    /* Save the original triangular superior terms in the lower triangle. */
    for (j = i + 1; j < tmu_len; j++)
    {
      ncm_matrix_set (self->cov_full, j, i, ncm_matrix_get (self->cov_full, i, j));
    }
  }

  /* Continue for the last 2 * mu_len. */
  for (i = mu_len; i < tmu_len; i++)
  {
    /* Save the original diagonal terms. */
    ncm_vector_set (self->cov_full_diag, i, ncm_matrix_get (self->cov_full, i, i));

    /* Save the original triangular superior terms in the lower triangle. */
    for (j = i + 1; j < tmu_len; j++)
    {
      ncm_matrix_set (self->cov_full, j, i, ncm_matrix_get (self->cov_full, i, j));
    }
  }
}

static void
_nc_data_snia_cov_resample (NcmData *data, NcmMSet *mset, NcmRNG *rng)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (data);
  NcDataSNIACovPrivate * const self = snia_cov->priv;

  if (self->has_complete_cov)
  {
    NcHICosmo *cosmo      = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
    NcSNIADistCov *dcov   = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
    const guint total_len = 3 * self->mu_len;

    gboolean dcov_resample_up  = ncm_model_ctrl_update (self->dcov_resample_ctrl, NCM_MODEL (dcov));
    gboolean cosmo_resample_up = ncm_model_ctrl_update (self->cosmo_resample_ctrl, NCM_MODEL (cosmo));
    gboolean dcov_cov_full_up  = ncm_model_ctrl_update (self->dcov_cov_full_ctrl, NCM_MODEL (dcov));

    gint ret;
    guint i;

    if (dcov_resample_up || cosmo_resample_up || !self->has_true_wc)
      nc_data_snia_cov_estimate_width_colour (snia_cov, mset);

    if (dcov_cov_full_up || (self->cov_full_state != NC_DATA_SNIA_COV_PREP_TO_RESAMPLE))
      _nc_data_snia_cov_prep_to_resample (snia_cov, dcov);

    ncm_rng_lock (rng);

    for (i = 0; i < total_len; i++)
    {
      const gdouble u_i = gsl_ran_ugaussian (rng->r);

      ncm_vector_set (self->mag_width_colour, i, u_i);
    }

    ret = gsl_blas_dtrmv (CblasUpper, CblasTrans, CblasNonUnit,
                          ncm_matrix_gsl (self->cov_full), ncm_vector_gsl (self->mag_width_colour));
    NCM_TEST_GSL_RESULT ("_nc_data_snia_cov_resample", ret);

    for (i = 0; i < self->mu_len; i++)
    {
      const gdouble width_th  = ncm_vector_get (self->width_true, i);
      const gdouble colour_th = ncm_vector_get (self->colour_true, i);
      const gdouble mag_th    = nc_snia_dist_cov_mag (dcov, cosmo, snia_cov, i, width_th, colour_th);
      const gdouble var_tot   = nc_snia_dist_cov_extra_var (dcov, snia_cov, i);

      const gdouble delta_mag    = ncm_vector_get (self->mag_width_colour, i + 0 * self->mu_len) + gsl_ran_ugaussian (rng->r) * sqrt (var_tot);
      const gdouble delta_width  = ncm_vector_get (self->mag_width_colour, i + 1 * self->mu_len);
      const gdouble delta_colour = ncm_vector_get (self->mag_width_colour, i + 2 * self->mu_len);

      ncm_vector_set (self->mag,    i, mag_th    - delta_mag);
      ncm_vector_set (self->width,  i, width_th  - delta_width);
      ncm_vector_set (self->colour, i, colour_th - delta_colour);
    }

    ncm_rng_unlock (rng);
  }
  else
  {
    /* Fallback to the default method. */
    g_warning ("_nc_data_snia_cov_resample: data set does not support full covariance resampling. Resampling from the reduced covariance.");
    NCM_DATA_CLASS (nc_data_snia_cov_parent_class)->resample (data, mset, rng);
  }
}

/**
 * nc_data_snia_cov_estimate_width_colour:
 * @snia_cov: a #NcDataSNIACov
 * @mset: a #NcmMSet
 *
 * Estimate the values of width and colour from the catalog using the models in @mset
 * and fitting the width and colour as free parameters.
 *
 * Returns: the value of the chisq for the fit.
 */
gdouble
nc_data_snia_cov_estimate_width_colour (NcDataSNIACov *snia_cov, NcmMSet *mset)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  const guint mu_len = self->mu_len;

  if (mu_len == 0)
  {
    g_error ("nc_data_snia_estimate_width_colour: empty catalog.");

    return 0.0;
  }
  else
  {
    NcHICosmo *cosmo                       = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
    NcSNIADistCov *dcov                    = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
    guint nobs                             = mu_len * 3;
    guint nparams                          = mu_len * 2;
    gsl_multifit_linear_workspace *mfit_ws = gsl_multifit_linear_alloc (nobs, nparams);
    NcmMatrix *cov                         = ncm_matrix_new (nparams, nparams);
    NcmMatrix *X                           = ncm_matrix_new (nobs, nparams);
    NcmVector *obs                         = ncm_vector_new (nobs);
    NcmVector *y                           = ncm_vector_new (nobs);
    NcmVector *params                      = ncm_vector_new (nparams);
    gdouble chisq                          = 0.0;

    _nc_data_snia_cov_prep_to_estimate (snia_cov, dcov);

#ifdef HAVE_LAPACK
    {
      NcmMatrix *L = ncm_matrix_new0 (nobs, nobs);
      GArray *ws   = ncm_lapack_dggglm_alloc (L, X, params, obs, y);
      nc_snia_dist_cov_mag_to_width_colour (dcov, cosmo, snia_cov, obs, X, TRUE);
      gint i, j, info;

      for (i = 0; i < nobs; i++)
      {
        for (j = i; j < nobs; j++)
        {
          ncm_matrix_set_colmajor (L, j, i, ncm_matrix_get (self->cov_full, i, j));
        }
      }

      info = ncm_lapack_dggglm_run (ws, L, X, params, obs, y);
      NCM_LAPACK_CHECK_INFO ("nc_data_snia_cov_estimate_width_colour", info);
      chisq = gsl_pow_2 (ncm_vector_dnrm2 (y));

      ncm_matrix_clear (&L);
      g_array_unref (ws);
    }

#else /* Fallback to GSL (much slower) */
    {
      gint ret;
      nc_snia_dist_cov_mag_to_width_colour (dcov, cosmo, snia_cov, obs, X, FALSE);
      ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit,
                            ncm_matrix_gsl (self->cov_full), ncm_vector_gsl (obs));
      NCM_TEST_GSL_RESULT ("nc_data_snia_estimate_width_colour", ret);

      ret = gsl_blas_dtrsm (CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, 1.0, ncm_matrix_gsl (self->cov_full), ncm_matrix_gsl (X));
      NCM_TEST_GSL_RESULT ("nc_data_snia_estimate_width_colour", ret);

      ret = gsl_multifit_linear (ncm_matrix_gsl (X),
                                 ncm_vector_gsl (obs),
                                 ncm_vector_gsl (params),
                                 ncm_matrix_gsl (cov),
                                 &chisq,
                                 mfit_ws);
      NCM_TEST_GSL_RESULT ("nc_data_snia_estimate_width_colour", ret);
    }
#endif

    {
      gint i;

      for (i = 0; i < mu_len; i++)
      {
        const gdouble width_true  = ncm_vector_get (params, i);
        const gdouble colour_true = ncm_vector_get (params, i + mu_len);

        ncm_vector_set (self->width_true, i, width_true);
        ncm_vector_set (self->colour_true, i, colour_true);
      }
    }
    self->has_true_wc = TRUE;

    if (chisq / (mu_len * 1.0) > 2.0)
      g_warning ("nc_data_snia_cov_estimate_width_colour: estimate with a very poor fit chisq = % 20.15g/%-20.15g = %20.15g", chisq, mu_len * 1.0, chisq / (mu_len * 1.0));

    ncm_vector_clear (&obs);
    ncm_vector_clear (&y);
    ncm_vector_clear (&params);
    ncm_matrix_clear (&X);
    ncm_matrix_clear (&cov);
    gsl_multifit_linear_free (mfit_ws);

    return chisq;
  }
}

/**
 * nc_data_snia_cov_get_estimated_mag:
 * @snia_cov: a #NcDataSNIACov
 * @mset: a #NcmMSet
 *
 * Estimate the values of width and colour from the catalog using the models in @mset
 * and fitting the width and colour as free parameters.
 *
 * Returns: (transfer full): the magnitude vector.
 */
NcmVector *
nc_data_snia_cov_get_estimated_mag (NcDataSNIACov *snia_cov, NcmMSet *mset)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  NcHICosmo *cosmo           = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcSNIADistCov *dcov        = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
  gboolean dcov_resample_up  = ncm_model_ctrl_update (self->dcov_resample_ctrl, NCM_MODEL (dcov));
  gboolean cosmo_resample_up = ncm_model_ctrl_update (self->cosmo_resample_ctrl, NCM_MODEL (cosmo));
  NcmVector *mag             = ncm_vector_new (self->mu_len);
  guint i;

  if (dcov_resample_up || cosmo_resample_up || !self->has_true_wc)
    nc_data_snia_cov_estimate_width_colour (snia_cov, mset);

  for (i = 0; i < self->mu_len; i++)
  {
    const gdouble width_th  = ncm_vector_get (self->width_true, i);
    const gdouble colour_th = ncm_vector_get (self->colour_true, i);
    const gdouble mag_th    = nc_snia_dist_cov_mag (dcov, cosmo, snia_cov, i, width_th, colour_th);

    ncm_vector_set (mag, i, mag_th);
  }

  return mag;
}

/**
 * nc_data_snia_cov_get_estimated_width:
 * @snia_cov: a #NcDataSNIACov
 * @mset: a #NcmMSet
 *
 * Estimate the values of width and colour from the catalog using the models in @mset
 * and fitting the width and colour as free parameters.
 *
 * Returns: (transfer full): the width vector.
 */
NcmVector *
nc_data_snia_cov_get_estimated_width (NcDataSNIACov *snia_cov, NcmMSet *mset)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  NcHICosmo *cosmo           = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcSNIADistCov *dcov        = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
  gboolean dcov_resample_up  = ncm_model_ctrl_update (self->dcov_resample_ctrl, NCM_MODEL (dcov));
  gboolean cosmo_resample_up = ncm_model_ctrl_update (self->cosmo_resample_ctrl, NCM_MODEL (cosmo));

  if (dcov_resample_up || cosmo_resample_up || !self->has_true_wc)
    nc_data_snia_cov_estimate_width_colour (snia_cov, mset);

  return ncm_vector_dup (self->width_true);
}

/**
 * nc_data_snia_cov_get_estimated_colour:
 * @snia_cov: a #NcDataSNIACov
 * @mset: a #NcmMSet
 *
 * Estimate the values of width and colour from the catalog using the models in @mset
 * and fitting the width and colour as free parameters.
 *
 * Returns: (transfer full): the colour vector.
 */
NcmVector *
nc_data_snia_cov_get_estimated_colour (NcDataSNIACov *snia_cov, NcmMSet *mset)
{
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  NcHICosmo *cosmo           = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcSNIADistCov *dcov        = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
  gboolean dcov_resample_up  = ncm_model_ctrl_update (self->dcov_resample_ctrl, NCM_MODEL (dcov));
  gboolean cosmo_resample_up = ncm_model_ctrl_update (self->cosmo_resample_ctrl, NCM_MODEL (cosmo));

  if (dcov_resample_up || cosmo_resample_up || !self->has_true_wc)
    nc_data_snia_cov_estimate_width_colour (snia_cov, mset);

  return ncm_vector_dup (self->colour_true);
}

static const gchar *_nc_data_snia_cats[] = {
  "snls3_conley_2011_sys_stat.fits",
  "snls3_conley_2011_stat_only.fits",
  "jla_snls3_sdss_sys_stat.fits",
  "jla_snls3_sdss_sys_stat_cmpl.fits",
  "snia_pantheon.fits",
  "snia_pantheon_plus_sh0es.fits",
  "snia_pantheon_plus_sh0es_statonly.fits",
};

/**
 * nc_data_snia_cov_get_catalog_id:
 * @id: string representation of the catalog id #NcDataSNIAId
 * @err: return location for a #GError, or NULL
 *
 * Gets catalog id from the @id string.
 *
 * Returns: catalog id
 */
NcDataSNIAId
nc_data_snia_cov_get_catalog_id (gchar *id, GError **err)
{
  const GEnumValue *snia_id =
    ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_SNIA_ID, id);

  if (snia_id == NULL)
  {
    g_set_error (err,
        NC_DATA_SNIA_COV_ERROR,
        NC_DATA_SNIA_COV_ERROR_ID_NOT_FOUND,
        "nc_data_snia_cov_get_catalog: Cannot find id ``%s'' in catalogs.", id);
    return -1;
  }

  return snia_id->value;
}

/**
 * nc_data_snia_cov_get_catalog:
 * @id: string representation of the catalog id #NcDataSNIAId
 * @err: return location for a #GError, or NULL
 *
 * Gets catalog filename from the @id.
 *
 * Returns: (transfer full): Catalog filename.
 */
gchar *
nc_data_snia_cov_get_catalog (gchar *id, GError **err)
{
  NcDataSNIAId idval = nc_data_snia_cov_get_catalog_id (id, err);
  g_return_val_if_fail (err == NULL, NULL);

  return nc_data_snia_cov_get_catalog_by_id (idval);
}

/**
 * nc_data_snia_cov_get_catalog_by_id:
 * @id: an id #NcDataSNIAId
 *
 * Gets catalog filename from the @id.
 *
 * Returns: (transfer full): Catalog filename.
 */
gchar *
nc_data_snia_cov_get_catalog_by_id (NcDataSNIAId id)
{
  g_assert ((id <= NC_DATA_SNIA_COV_END) && (id >= NC_DATA_SNIA_COV_START));
  {
    gint i               = id - NC_DATA_SNIA_COV_START;
    gchar *full_filename = nc_data_snia_cov_get_fits (_nc_data_snia_cats[i], FALSE);

    return full_filename;
  }
}

void
_nc_data_snia_copy_prog (goffset current_num_bytes, goffset total_num_bytes, gpointer user_data)
{
  gint *old_prog = (gint *) user_data;
  gint prog      = (100 * current_num_bytes) / total_num_bytes;

  if (prog > *old_prog)
  {
    ncm_message ("# % 3d%%\r", prog);
    *old_prog = prog;
  }
}

/**
 * nc_data_snia_cov_get_fits:
 * @filename: Catalog fits filename
 * @check_size: Whether to check if the file size matches
 *
 * Downloads catalog from the repository.
 *
 * Returns: (transfer full): Filename of the downloaded catalog.
 */
gchar *
nc_data_snia_cov_get_fits (const gchar *filename, gboolean check_size)
{
  gchar *full_filename = ncm_cfg_get_fullpath (filename);
  gchar *url_str       = g_strdup_printf ("https://github.com/NumCosmo/NumCosmo/releases/download/v"PACKAGE_VERSION"/%s", filename);
  GFile *local         = g_file_new_for_path (full_filename);
  GFile *remote        = g_file_new_for_uri (url_str);
  GError *error        = NULL;
  gint prog            = 0;
  gboolean download    = FALSE;

  if (g_file_test (full_filename, G_FILE_TEST_EXISTS))
  {
    if (check_size)
    {
      GFileInfo *local_info, *remote_info;

      local_info = g_file_query_info (local, G_FILE_ATTRIBUTE_STANDARD_SIZE,
                                      G_FILE_QUERY_INFO_NONE,
                                      NULL, &error);

      if (local_info == NULL)
        g_error ("nc_data_snia_cov_get_fits: cannot get info for %s: %s.", full_filename, error->message);

      remote_info = g_file_query_info (remote, G_FILE_ATTRIBUTE_STANDARD_SIZE,
                                       G_FILE_QUERY_INFO_NONE,
                                       NULL, &error);

      if (remote_info == NULL)
        g_error ("nc_data_snia_cov_get_fits: cannot get info for %s: %s."
                 " To use this catalog, download the file from the url and copy "
                 "to ~/.numcosmo directory.", url_str, error->message);

      if (g_file_info_get_size (local_info) != g_file_info_get_size (remote_info))
        download = TRUE;
    }
  }
  else
  {
    download = TRUE;
  }

  if (download)
  {
    ncm_message ("# Downloading file [%s]...\n", url_str);

    if (!g_file_copy (remote, local, G_FILE_COPY_OVERWRITE, NULL,
                      &_nc_data_snia_copy_prog, &prog, &error))
      g_error ("nc_data_snia_cov_get_fits: cannot get fits file from %s: %s."
               " To use this catalog, download the file from the url and copy "
               "to ~/.numcosmo directory.",
               url_str, error->message);
  }

  g_free (url_str);

  g_object_unref (local);
  g_object_unref (remote);

  return full_filename;
}

/**
 * nc_data_snia_cov_apply_filter_sh0es_z:
 * @snia_cov: a #NcDataSNIACov
 * @z_min: Minimum redshift
 * @use_calib: Whether to use calibrators
 * @err: return location for a #GError, or NULL
 *
 * Apply SH0ES+z_min filter to data.
 *
 * Returns: (transfer full): the new #NcDataSNIACov with filtered data.
 */
NcDataSNIACov *
nc_data_snia_cov_apply_filter_sh0es_z (NcDataSNIACov *snia_cov, const gdouble z_min, const gboolean use_calib, GError **err)
{
  NcDataSNIACov *snia_cov_filter = NULL;
  NcDataSNIACovPrivate * const self = snia_cov->priv;
  NcmISet *is = ncm_iset_new (self->mu_len);
  gint i;

  if (self->cat_version < 2)
  {
    g_set_error (err,
        NC_DATA_SNIA_COV_ERROR,
        NC_DATA_SNIA_COV_ERROR_ID_NOT_FOUND,
        "Filter not compatible with catalog version <2");
    ncm_iset_free (is);
    return NULL;
  }

  for (i = 0; i < self->mu_len; i++)
  {
    const gdouble z_hd = ncm_vector_get (self->z_hd, i);
    if ((z_hd > z_min) || (use_calib && g_array_index (self->is_calib, guint32, i)))
    {
      ncm_iset_add (is, i);
    }
  }

  {
    snia_cov_filter = nc_data_snia_cov_new (FALSE, self->cat_version);
    const guint new_len = ncm_iset_get_len (is);

    ncm_data_gauss_cov_set_size (NCM_DATA_GAUSS_COV (snia_cov_filter), new_len);
    ncm_data_set_desc (NCM_DATA (snia_cov_filter), ncm_data_peek_desc (NCM_DATA (snia_cov)));

    {
      NcDataSNIACovPrivate * const selff = snia_cov_filter->priv;
      NcmVector *vec_from[34] = {self->z_hd, self->z_cmb, self->z_he, self->mag_b_corr, self->mu_corr,
          self->ceph_dist, self->ampl, self->vpec, self->mag, self->width, self->colour, self->thirdpar,
          self->ra, self->dec, self->host_ra, self->host_dec, self->host_ang_sep, self->width_true,
          self->colour_true, self->sigma_z_hd, self->sigma_z_cmb, self->sigma_z_he,
          self->sigma_z, self->sigma_mag_b_corr, self->sigma_mu_corr, self->sigma_colour, self->sigma_width,
          self->sigma_mag, self->sigma_ampl, self->sigma_vpec, self->sigma_thirdpar, self->cov_width_colour,
          self->cov_width_ampl, self->cov_colour_ampl};
      NcmVector *vec_to[34] = {selff->z_hd, selff->z_cmb, selff->z_he, selff->mag_b_corr, selff->mu_corr,
          selff->ceph_dist, selff->ampl, selff->vpec, selff->mag, selff->width, selff->colour, selff->thirdpar,
          selff->ra, selff->dec, selff->host_ra, selff->host_dec, selff->host_ang_sep, selff->width_true,
          selff->colour_true, selff->sigma_z_hd, selff->sigma_z_cmb, selff->sigma_z_he,
          selff->sigma_z, selff->sigma_mag_b_corr, selff->sigma_mu_corr, selff->sigma_colour, selff->sigma_width,
          selff->sigma_mag, selff->sigma_ampl, selff->sigma_vpec, selff->sigma_thirdpar, selff->cov_width_colour,
          selff->cov_width_ampl, selff->cov_colour_ampl};

      for (i = 0; i < 34; i++)
      {
        ncm_iset_get_subvector (is, vec_from[i], vec_to[i]);
      }

      ncm_iset_get_submatrix (is, self->cov_mbc_mbc, selff->cov_mbc_mbc);
      ncm_iset_get_subarray (is, self->is_calib, selff->is_calib);
      ncm_iset_get_subarray (is, self->used_in_sh0es, selff->used_in_sh0es);
      ncm_iset_get_subarray (is, self->dataset, selff->dataset);

      /*printf ("New catalog has %u\n", selff->mu_len);*/
    }
  }

  ncm_iset_clear (&is);

  snia_cov = snia_cov_filter;
  _NC_DATA_SNIA_COV_SET_DATA_INIT_ALL;

  return snia_cov_filter;
}

