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
 * @title: Supernovae Ia Data -- Covariance
 * @short_description: SNIa data with covariance error matrix
 * 
 * See #NcSNIADistCov.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_snia_cov.h"
#include "nc_snia_dist_cov.h"
#include "math/ncm_model_ctrl.h"

#include <glib/gstdio.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

enum
{
  PROP_0,
  PROP_SIGMA_PECZ,
  PROP_DATA_FILENAME,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcDataSNIACov, nc_data_snia_cov, NCM_TYPE_DATA_GAUSS_COV);

static void
nc_data_snia_cov_init (NcDataSNIACov *snia_cov)
{
  snia_cov->mu_len            = 0;
  
  snia_cov->z_cmb             = NULL;
  snia_cov->z_he              = NULL;

  snia_cov->mag               = NULL;
  snia_cov->width             = NULL;
  snia_cov->colour            = NULL;
  snia_cov->thirdpar          = NULL;

  snia_cov->sigma_z           = NULL;
  snia_cov->sigma_mag         = NULL;
  snia_cov->sigma_width       = NULL;
  snia_cov->sigma_colour      = NULL;
  snia_cov->sigma_thirdpar    = NULL;

  snia_cov->diag_mag_width    = NULL;
  snia_cov->diag_mag_colour   = NULL;
  snia_cov->diag_width_colour = NULL;

  snia_cov->var_mag           = NULL;
  snia_cov->var_width         = NULL;
  snia_cov->var_colour        = NULL;
  snia_cov->var_mag_width     = NULL;
  snia_cov->var_mag_colour    = NULL;
  snia_cov->var_width_colour  = NULL;

  snia_cov->sigma_pecz        = 0.0;
  snia_cov->dataset           = NULL;
  snia_cov->dataset_len       = 0;

  snia_cov->filename          = NULL;
}

static void
nc_data_snia_cov_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (object);
  g_return_if_fail (NC_IS_DATA_SNIA_COV (object));

  switch (prop_id)
  {
    case PROP_SIGMA_PECZ:
      snia_cov->sigma_pecz = g_value_get_double (value);
      break;
    case PROP_DATA_FILENAME:
    {
      const gchar *filename = g_value_get_string (value);  
      if (filename != NULL)
        nc_data_snia_cov_load (snia_cov, filename);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_snia_cov_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (object);
  g_return_if_fail (NC_IS_DATA_SNIA_COV (object));

  switch (prop_id)
  {
    case PROP_SIGMA_PECZ:
      g_value_set_double (value, snia_cov->sigma_pecz);
      break;
    case PROP_DATA_FILENAME:
      g_value_set_string (value, snia_cov->filename);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_data_snia_cov_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_data_snia_cov_parent_class)->constructed (object);
  {


    ncm_data_set_init (NCM_DATA (object));
  }
}

static void
nc_data_snia_cov_dispose (GObject *object)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (object);
  nc_data_snia_cov_set_size (snia_cov, 0);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_snia_cov_parent_class)->dispose (object);
}

static void
nc_data_snia_cov_finalize (GObject *object)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (object);
  g_free (snia_cov->filename);
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_data_snia_cov_parent_class)->finalize (object);
}

static NcmData *_nc_data_snia_cov_dup (NcmData *data);
static void _nc_data_snia_cov_copyto (NcmData *data, NcmData *data_dest);
static void _nc_data_snia_cov_prepare (NcmData *data, NcmMSet *mset);
static void _nc_data_snia_cov_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);
static gboolean _nc_data_snia_cov_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *cov);

static void
nc_data_snia_cov_class_init (NcDataSNIACovClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmDataClass *data_class   = NCM_DATA_CLASS (klass);
  NcmDataGaussCovClass* gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->set_property = &nc_data_snia_cov_set_property;
  object_class->get_property = &nc_data_snia_cov_get_property;
  object_class->constructed  = &nc_data_snia_cov_constructed;
  object_class->dispose      = &nc_data_snia_cov_dispose;
  object_class->finalize     = &nc_data_snia_cov_finalize;

  g_object_class_install_property (object_class,
                                   PROP_SIGMA_PECZ,
                                   g_param_spec_double ("sigma-pecz",
                                                        NULL,
                                                        "Error from SN Ia peculiar velocity",
                                                        0.0, 1.0e1, 5.0e-4,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_DATA_FILENAME,
                                   g_param_spec_string ("data-filename",
                                                        NULL,
                                                        "Data filename",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  data_class->dup        = &_nc_data_snia_cov_dup;
  data_class->copyto     = &_nc_data_snia_cov_copyto;
  data_class->prepare    = &_nc_data_snia_cov_prepare;
  gauss_class->mean_func = &_nc_data_snia_cov_mean_func;
  gauss_class->cov_func  = &_nc_data_snia_cov_func;
  
}

static NcmData *
_nc_data_snia_cov_dup (NcmData *data)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (data);
  NcDataSNIACov *snia_cov_dup = g_object_new (NC_TYPE_DATA_SNIA_COV,
                                              NULL);
  NcmData *data_dest = NCM_DATA (snia_cov_dup);
  nc_data_snia_cov_set_size (snia_cov_dup, snia_cov->mu_len);

  ncm_data_copyto (data, data_dest);
  
  return data_dest;
}

static void 
_nc_data_snia_cov_copyto (NcmData *data, NcmData *data_dest)
{
  /* Chain up : start */
  NCM_DATA_CLASS (nc_data_snia_cov_parent_class)->copyto (data, data_dest);
  {
    NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (data);
    NcDataSNIACov *snia_cov_dest = NC_DATA_SNIA_COV (data_dest);

    g_assert_cmpuint (snia_cov->mu_len, ==, snia_cov_dest->mu_len);
#define _COPY_VEC(vec) ncm_vector_memcpy (snia_cov_dest->vec, snia_cov->vec)
#define _COPY_MAT(mat) ncm_matrix_memcpy (snia_cov_dest->mat, snia_cov->mat)
    _COPY_VEC (z_cmb);
    _COPY_VEC (z_he);
    _COPY_VEC (mag);
    _COPY_VEC (width);
    _COPY_VEC (colour);
    _COPY_VEC (thirdpar);
    _COPY_VEC (sigma_z);
    _COPY_VEC (sigma_mag);
    _COPY_VEC (sigma_mag);
    _COPY_VEC (sigma_width);
    _COPY_VEC (sigma_colour);
    _COPY_VEC (sigma_thirdpar);
    _COPY_VEC (diag_mag_width);
    _COPY_VEC (diag_mag_colour);
    _COPY_VEC (diag_width_colour);
    _COPY_MAT (var_mag);
    _COPY_MAT (var_width);
    _COPY_MAT (var_colour);
    _COPY_MAT (var_mag_width);
    _COPY_MAT (var_mag_colour);
    _COPY_MAT (var_width_colour);
    g_assert (snia_cov->dataset != NULL && snia_cov_dest->dataset != NULL);
    NCM_GARRAY_MEMCPY (snia_cov_dest->dataset, snia_cov->dataset);
    snia_cov_dest->dataset_len = snia_cov->dataset_len;
    snia_cov_dest->sigma_pecz  = snia_cov->sigma_pecz;
    if (snia_cov->filename != NULL)
      snia_cov_dest->filename = g_strdup (snia_cov->filename);
    else
      snia_cov_dest->filename = NULL;
#undef _COPY_VEC
#undef _COPY_MAT
  }
}

static void
_nc_data_snia_cov_prepare (NcmData *data, NcmMSet *mset)
{
  NcSNIADistCov *dcov = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
  nc_snia_dist_cov_prepare_if_needed (dcov, mset);
}

static void 
_nc_data_snia_cov_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (gauss);
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  NcSNIADistCov *dcov = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
  nc_snia_dist_cov_mean (dcov, cosmo, snia_cov, vp);
}

static gboolean 
_nc_data_snia_cov_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmMatrix *cov)
{
  NcDataSNIACov *snia_cov = NC_DATA_SNIA_COV (gauss);
  NcSNIADistCov *dcov = NC_SNIA_DIST_COV (ncm_mset_peek (mset, nc_snia_dist_cov_id ()));
  nc_snia_dist_cov_calc (dcov, snia_cov, cov);
  
  return TRUE;
}

/**
 * nc_data_snia_cov_new:
 * @use_det: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcmData *
nc_data_snia_cov_new (gboolean use_det)
{
  return g_object_new (NC_TYPE_DATA_SNIA_COV,
                       "use-det", use_det,
                       NULL);
}

/**
 * nc_data_snia_cov_new_full:
 * @filename: FIXME
 * @use_det: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
NcmData *
nc_data_snia_cov_new_full (gchar *filename, gboolean use_det)
{
  return g_object_new (NC_TYPE_DATA_SNIA_COV,
                       "use-det", use_det,
                       "data-filename", filename,
                       NULL);
}

/**
 * nc_data_snia_cov_set_size:
 * @snia_cov: FIXME
 * @mu_len: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_data_snia_cov_set_size (NcDataSNIACov *snia_cov, guint mu_len)
{
  NcmDataGaussCov *gauss = NCM_DATA_GAUSS_COV (snia_cov);
  
  if (mu_len == 0 || mu_len != snia_cov->mu_len)
  {
    ncm_data_gauss_cov_set_size (gauss, 0);
    ncm_vector_clear (&snia_cov->z_cmb);
    ncm_vector_clear (&snia_cov->z_he);

    ncm_vector_clear (&snia_cov->mag);
    ncm_vector_clear (&snia_cov->width);
    ncm_vector_clear (&snia_cov->colour);
    ncm_vector_clear (&snia_cov->thirdpar);

    ncm_vector_clear (&snia_cov->sigma_z);
    ncm_vector_clear (&snia_cov->sigma_mag);
    ncm_vector_clear (&snia_cov->sigma_width);
    ncm_vector_clear (&snia_cov->sigma_colour);
    ncm_vector_clear (&snia_cov->sigma_thirdpar);

    ncm_vector_clear (&snia_cov->diag_mag_width);
    ncm_vector_clear (&snia_cov->diag_mag_colour);
    ncm_vector_clear (&snia_cov->diag_width_colour);
    
    ncm_matrix_clear (&snia_cov->var_mag);
    ncm_matrix_clear (&snia_cov->var_width);
    ncm_matrix_clear (&snia_cov->var_colour);

    ncm_matrix_clear (&snia_cov->var_mag_width);
    ncm_matrix_clear (&snia_cov->var_mag_colour);
    ncm_matrix_clear (&snia_cov->var_width_colour);

    if (snia_cov->dataset != NULL)
    {
      g_array_unref (snia_cov->dataset);
      snia_cov->dataset     = NULL;
      snia_cov->dataset_len = 0;
    }
    NCM_DATA (snia_cov)->init = FALSE;
  }

  if (mu_len > 0 && mu_len != snia_cov->mu_len)
  {
    ncm_data_gauss_cov_set_size (gauss, mu_len);
    
    snia_cov->mu_len           = mu_len;

    snia_cov->z_cmb            = ncm_vector_new (mu_len);
    snia_cov->z_he             = ncm_vector_new (mu_len);

    g_assert_cmpuint (ncm_vector_len (gauss->y), ==, mu_len);
    
    snia_cov->mag              = ncm_vector_ref (gauss->y);
    snia_cov->width            = ncm_vector_new (mu_len);
    snia_cov->colour           = ncm_vector_new (mu_len);
    snia_cov->thirdpar         = ncm_vector_new (mu_len);

    snia_cov->sigma_z          = ncm_vector_new (mu_len);
    snia_cov->sigma_mag        = ncm_vector_new (mu_len);
    snia_cov->sigma_width      = ncm_vector_new (mu_len);
    snia_cov->sigma_colour     = ncm_vector_new (mu_len);
    snia_cov->sigma_thirdpar   = ncm_vector_new (mu_len);

    snia_cov->diag_mag_width    = ncm_vector_new (mu_len);
    snia_cov->diag_mag_colour   = ncm_vector_new (mu_len);
    snia_cov->diag_width_colour = ncm_vector_new (mu_len);

    snia_cov->var_mag          = ncm_matrix_new (mu_len, mu_len);
    snia_cov->var_width        = ncm_matrix_new (mu_len, mu_len);
    snia_cov->var_colour       = ncm_matrix_new (mu_len, mu_len);

    snia_cov->var_mag_width    = ncm_matrix_new (mu_len, mu_len);
    snia_cov->var_mag_colour   = ncm_matrix_new (mu_len, mu_len);
    snia_cov->var_width_colour = ncm_matrix_new (mu_len, mu_len);

    snia_cov->dataset          = g_array_sized_new (FALSE, FALSE, sizeof (guint32), mu_len);
    snia_cov->dataset_len      = 0;

    g_array_set_size (snia_cov->dataset, mu_len);

    NCM_DATA (snia_cov)->init = FALSE;
  }
}

static void _nc_data_snia_cov_load_snia_data (NcDataSNIACov *snia_cov, const gchar *filename);
static void _nc_data_snia_cov_load_matrix (const gchar *filename, NcmMatrix *data);

/**
 * nc_data_snia_cov_load_txt:
 * @snia_cov: FIXME
 * @filename: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_data_snia_cov_load_txt (NcDataSNIACov *snia_cov, const gchar *filename)
{
  GKeyFile *snia_keyfile = g_key_file_new ();
  GError *error   = NULL;
  gchar *datafile = NULL;
  guint64 mu_len;

  if (filename == NULL)
    g_error ("nc_data_snia_cov_load: null filename");
  else if (snia_cov->filename == NULL)
  {
    snia_cov->filename = g_strdup (filename);
  }
  else if (snia_cov->filename != filename && g_strcmp0 (filename, snia_cov->filename) != 0)
  {
    g_free (snia_cov->filename);
    snia_cov->filename = g_strdup (filename);
  }

  if (!g_key_file_load_from_file (snia_keyfile, filename, G_KEY_FILE_NONE, &error))
    g_error ("nc_data_snia_cov_load: invalid configuration: %s\n  %s\n", 
             filename, 
             error->message);

  if (!g_key_file_has_key (snia_keyfile, 
                           NC_DATA_SNIA_COV_DATA_GROUP,
                           NC_DATA_SNIA_COV_DATA_KEY,
                           &error))
    g_error ("nc_data_snia_cov_load: invalid %s key file, it must define at least the data file key ["NC_DATA_SNIA_COV_DATA_KEY"]", filename);
  if (!g_key_file_has_key (snia_keyfile, 
                           NC_DATA_SNIA_COV_DATA_GROUP,
                           NC_DATA_SNIA_COV_DATA_LEN_KEY,
                           &error))
    g_error ("nc_data_snia_cov_load: invalid %s key file, it must define at least the data length file key ["NC_DATA_SNIA_COV_DATA_LEN_KEY"]", filename);

  datafile = g_key_file_get_string (snia_keyfile, 
                                    NC_DATA_SNIA_COV_DATA_GROUP,
                                    NC_DATA_SNIA_COV_DATA_KEY,
                                    &error);

  mu_len = g_key_file_get_uint64 (snia_keyfile, 
                                  NC_DATA_SNIA_COV_DATA_GROUP,
                                  NC_DATA_SNIA_COV_DATA_LEN_KEY,
                                  &error);

  nc_data_snia_cov_set_size (snia_cov, mu_len);
  
  _nc_data_snia_cov_load_snia_data (snia_cov, datafile);

  /* Get magnitude cov matrix */
  
  if (g_key_file_has_key (snia_keyfile, 
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_MAG_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_DATA_SNIA_COV_DATA_GROUP,
                                      NC_DATA_SNIA_COV_MAG_KEY,
                                      &error);
    
    _nc_data_snia_cov_load_matrix (datafile, snia_cov->var_mag);
  }

  /* Get width cov matrix */

  if (g_key_file_has_key (snia_keyfile, 
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_WIDTH_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_DATA_SNIA_COV_DATA_GROUP,
                                      NC_DATA_SNIA_COV_WIDTH_KEY,
                                      &error);
    
    _nc_data_snia_cov_load_matrix (datafile, snia_cov->var_width);
  }

  /* Get colour cov matrix */

  if (g_key_file_has_key (snia_keyfile, 
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_COLOUR_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_DATA_SNIA_COV_DATA_GROUP,
                                      NC_DATA_SNIA_COV_COLOUR_KEY,
                                      &error);
    
    _nc_data_snia_cov_load_matrix (datafile, snia_cov->var_colour);
  }

  /* Get magnitude-width cov matrix */

  if (g_key_file_has_key (snia_keyfile, 
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_MAG_WIDTH_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_DATA_SNIA_COV_DATA_GROUP,
                                      NC_DATA_SNIA_COV_MAG_WIDTH_KEY,
                                      &error);
    
    _nc_data_snia_cov_load_matrix (datafile, snia_cov->var_mag_width);
  }

  /* Get magnitude-colour cov matrix */

  if (g_key_file_has_key (snia_keyfile, 
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_MAG_COLOUR_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_DATA_SNIA_COV_DATA_GROUP,
                                      NC_DATA_SNIA_COV_MAG_COLOUR_KEY,
                                      &error);
    
    _nc_data_snia_cov_load_matrix (datafile, snia_cov->var_mag_colour);
  }

  /* Get width-colour cov matrix */

  if (g_key_file_has_key (snia_keyfile, 
                          NC_DATA_SNIA_COV_DATA_GROUP,
                          NC_DATA_SNIA_COV_WIDTH_COLOUR_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_DATA_SNIA_COV_DATA_GROUP,
                                      NC_DATA_SNIA_COV_WIDTH_COLOUR_KEY,
                                      &error);
    
    _nc_data_snia_cov_load_matrix (datafile, snia_cov->var_width_colour);
  }

  ncm_data_set_init (NCM_DATA (snia_cov));
 
  g_key_file_free (snia_keyfile);
}

static void 
_nc_data_snia_cov_load_snia_data (NcDataSNIACov *snia_cov, const gchar *filename)
{
  GArray *dataset = snia_cov->dataset;
  gchar *line = NULL;
  gsize len = 0, tpos = 0;
  GError *error = NULL;
  GIOChannel *file = g_io_channel_new_file (filename, "r", &error);
  if (file == NULL)
    g_error ("_nc_data_snia_cov_load_snia_data: cannot open file %s: %s", 
             filename, error->message);

  {
    GRegex *comment_line = g_regex_new ("^\\s*#", 0, 0, &error);
    guint n = 0;
    guint nrow = 0;
    guint i;
    guint max_dset_id = 0;
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
          if (*itens[n-1] == '\0')
            pad_end++;
          if (n - 1 - pad_end == NC_DATA_SNIA_COV_LENGTH)
            has_dataset = FALSE;
          else if (n - 2 - pad_end == NC_DATA_SNIA_COV_LENGTH)
            has_dataset = TRUE;
          else
            g_error ("_nc_data_snia_cov_load_snia_data: data file [%s] must have %d or %d columns, it has %u", 
                     filename, NC_DATA_SNIA_COV_LENGTH + 1, NC_DATA_SNIA_COV_LENGTH + 2, n - pad_end);
        }
        else if (n != g_strv_length (itens))
          g_error ("_nc_data_snia_cov_load_snia_data: data file [%s] has different number of columns in different rows [%u]", filename, nrow);
        
        if (nrow >= snia_cov->mu_len)
          g_error ("_nc_data_snia_cov_load_snia_data: cannot load data file [%s] expected nrows %u obtained >%u\n", 
                   filename, snia_cov->mu_len, nrow);

        {
          NcmVector *data[NC_DATA_SNIA_COV_LENGTH];
          data[NC_DATA_SNIA_COV_ZCMB]              = snia_cov->z_cmb;
          data[NC_DATA_SNIA_COV_ZHE]               = snia_cov->z_he;
          data[NC_DATA_SNIA_COV_SIGMA_Z]           = snia_cov->sigma_z;
          data[NC_DATA_SNIA_COV_MAG]               = snia_cov->mag;
          data[NC_DATA_SNIA_COV_SIGMA_MAG]         = snia_cov->sigma_mag;
          data[NC_DATA_SNIA_COV_WIDTH]             = snia_cov->width;
          data[NC_DATA_SNIA_COV_SIGMA_WIDTH]       = snia_cov->sigma_width;
          data[NC_DATA_SNIA_COV_COLOUR]            = snia_cov->colour;
          data[NC_DATA_SNIA_COV_SIGMA_COLOUR]      = snia_cov->sigma_colour;
          data[NC_DATA_SNIA_COV_THIRDPAR]          = snia_cov->thirdpar;
          data[NC_DATA_SNIA_COV_SIGMA_THIRDPAR]    = snia_cov->sigma_thirdpar;
          data[NC_DATA_SNIA_COV_DIAG_MAG_WIDTH]    = snia_cov->diag_mag_width;
          data[NC_DATA_SNIA_COV_DIAG_MAG_COLOUR]   = snia_cov->diag_mag_colour;
          data[NC_DATA_SNIA_COV_DIAG_WIDTH_COLOUR] = snia_cov->diag_width_colour;

          for (i = 0; i < NC_DATA_SNIA_COV_LENGTH; i++)
          {
            const gdouble val = g_ascii_strtod (itens[i + 1], NULL);
            ncm_vector_set (data[i], nrow, val);
          }
        }
        
        if (has_dataset)
        {
          gint64 dset_id = g_ascii_strtoll (itens[i + 1], NULL, 10);
          g_array_index (dataset, guint32, nrow) = dset_id;
          max_dset_id = GSL_MAX (max_dset_id, dset_id);
        }
        else
          g_array_index (dataset, guint32, nrow) = 0;

        nrow++;
        g_strfreev (itens);
        g_free (line);
      }
    }

    if (nrow != snia_cov->mu_len)
      g_error ("_nc_data_snia_cov_load_snia_data: cannot load data file [%s] expected nrows %u obtained %u\n", 
               filename, snia_cov->mu_len, nrow);

    snia_cov->dataset_len = max_dset_id + 1;

    g_regex_unref (comment_line);
  }
  
  g_io_channel_unref (file);
}

static void 
_nc_data_snia_cov_load_matrix (const gchar *filename, NcmMatrix *data)
{
  GError *error = NULL;
  gchar *file = NULL;
  gsize length = 0;
  gchar **itens;
  guint itens_len = 0;
  guint pad_start = 0;
  guint pad_end = 0;
  gint64 matrix_len = 0;
  guint i, j;
  if (!g_file_get_contents (filename, &file, &length, &error))
    g_error ("_nc_data_snia_cov_load_matrix: cannot open file %s: %s", 
             filename, error->message);

  itens = g_regex_split_simple ("\\s+", file, 0, 0);
  itens_len = g_strv_length (itens);

  if (*itens[0] == '\0')
    pad_start = 1;
  if (*itens[itens_len - 1] == '\0')
    pad_end = 1;
  
  if (pad_start) itens_len--;
  if (pad_end) itens_len--;

  matrix_len = g_ascii_strtoll (itens[pad_start], NULL, 10);
  pad_start++;
  itens_len--;
  
  if (matrix_len != NCM_MATRIX_NROWS (data))
    g_error ("_nc_data_snia_cov_load_matrix: expected a %zux%zu matrix but got %"G_GINT64_MODIFIER"dx%"G_GINT64_MODIFIER"d", 
             NCM_MATRIX_NROWS (data), NCM_MATRIX_NROWS (data), matrix_len, matrix_len);
  
  if (matrix_len * matrix_len != itens_len)
    g_error ("_nc_data_snia_cov_load_matrix: matrix header say %"G_GINT64_MODIFIER"d length but got %u\n", 
             matrix_len * matrix_len , itens_len);

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

#ifdef NUMCOSMO_HAVE_CFITSIO

#define NC_FITS_ERROR(status) \
do { \
if (status) \
{ \
gchar errormsg[30]; \
fits_get_errstatus (status, errormsg); \
g_error ("FITS: %s", errormsg); \
} \
} while (FALSE);

/**
 * nc_data_snia_cov_load:
 * @snia_cov: FIXME
 * @filename: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_data_snia_cov_load (NcDataSNIACov *snia_cov, const gchar *filename)
{
  gint status, hdutype;
  fitsfile *fptr;

  status = 0;

  if (filename == NULL)
    g_error ("nc_data_snia_cov_load: null filename");
  else if (snia_cov->filename == NULL)
  {
    snia_cov->filename = g_strdup (filename);
  }
  else if (snia_cov->filename != filename && g_strcmp0 (filename, snia_cov->filename) != 0)
  {
    g_free (snia_cov->filename);
    snia_cov->filename = g_strdup (filename);
  }

  if (fits_open_file (&fptr, filename, READONLY, &status))
    NC_FITS_ERROR (status);

  if (fits_movabs_hdu (fptr, 2, &hdutype, &status))
    NC_FITS_ERROR (status);

  if (hdutype != BINARY_TBL)
    g_error ("nc_data_snia_cov_load: NcSNIADistCov catalog is not binary.");

  {
    glong nrows;
    if (fits_get_num_rows (fptr, &nrows, &status))
      NC_FITS_ERROR (status);
    nc_data_snia_cov_set_size (snia_cov, nrows);
  }

  {
    NcmVector *data[NC_DATA_SNIA_COV_LENGTH];
    gint i;
    
    data[NC_DATA_SNIA_COV_ZCMB]              = snia_cov->z_cmb;
    data[NC_DATA_SNIA_COV_ZHE]               = snia_cov->z_he;
    data[NC_DATA_SNIA_COV_SIGMA_Z]           = snia_cov->sigma_z;
    data[NC_DATA_SNIA_COV_MAG]               = snia_cov->mag;
    data[NC_DATA_SNIA_COV_SIGMA_MAG]         = snia_cov->sigma_mag;
    data[NC_DATA_SNIA_COV_WIDTH]             = snia_cov->width;
    data[NC_DATA_SNIA_COV_SIGMA_WIDTH]       = snia_cov->sigma_width;
    data[NC_DATA_SNIA_COV_COLOUR]            = snia_cov->colour;
    data[NC_DATA_SNIA_COV_SIGMA_COLOUR]      = snia_cov->sigma_colour;
    data[NC_DATA_SNIA_COV_THIRDPAR]          = snia_cov->thirdpar;
    data[NC_DATA_SNIA_COV_SIGMA_THIRDPAR]    = snia_cov->sigma_thirdpar;
    data[NC_DATA_SNIA_COV_DIAG_MAG_WIDTH]    = snia_cov->diag_mag_width;
    data[NC_DATA_SNIA_COV_DIAG_MAG_COLOUR]   = snia_cov->diag_mag_colour;
    data[NC_DATA_SNIA_COV_DIAG_WIDTH_COLOUR] = snia_cov->diag_width_colour;

    for (i = 0; i < NC_DATA_SNIA_COV_LENGTH; i++)
    {
      if (fits_read_col_dbl (fptr, i + 1, 1, 1, snia_cov->mu_len, GSL_NAN,
                             ncm_vector_ptr (data[i], 0), NULL, 
                             &status))
        NC_FITS_ERROR(status);
    }
  }

  if (fits_read_col_uint (fptr, NC_DATA_SNIA_COV_ABSMAG_SET + 1, 1, 1, 
                          snia_cov->mu_len, 
                          0, &g_array_index (snia_cov->dataset, guint32, 0), 
                          NULL, &status))
    NC_FITS_ERROR (status);
  {
    gint i;
    for (i = 0; i < snia_cov->mu_len; i++)
    {
      snia_cov->dataset_len = GSL_MAX (g_array_index (snia_cov->dataset, guint32, i) + 1, snia_cov->dataset_len);
    }
  }

  {
    NcmMatrix *data[NC_DATA_SNIA_COV_TOTAL_LENGTH];
    data[NC_DATA_SNIA_COV_VAR_MAG]          = snia_cov->var_mag;
    data[NC_DATA_SNIA_COV_VAR_WIDTH]        = snia_cov->var_width;
    data[NC_DATA_SNIA_COV_VAR_COLOUR]       = snia_cov->var_colour;
    data[NC_DATA_SNIA_COV_VAR_MAG_WIDTH]    = snia_cov->var_mag_width;
    data[NC_DATA_SNIA_COV_VAR_MAG_COLOUR]   = snia_cov->var_mag_colour;
    data[NC_DATA_SNIA_COV_VAR_WIDTH_COLOUR] = snia_cov->var_width_colour;
    
    gint i;
    for (i = NC_DATA_SNIA_COV_VAR_MAG; i < NC_DATA_SNIA_COV_TOTAL_LENGTH; i++)
    {
      if (fits_read_col_dbl (fptr, i + 1, 1, 1, snia_cov->mu_len * snia_cov->mu_len, 
                             GSL_NAN, ncm_matrix_ptr (data[i], 0, 0), 
                             NULL, &status))
        NC_FITS_ERROR(status);
    }
  }

  ncm_data_set_init (NCM_DATA (snia_cov));
}

/**
 * nc_data_snia_cov_save:
 * @snia_cov: FIXME
 * @filename: FIXME
 * @overwrite: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_data_snia_cov_save (NcDataSNIACov *snia_cov, const gchar *filename, gboolean overwrite)
{
  /*******************************************************************/
  /* Create a binary table extension                                 */
  /*******************************************************************/
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  gint status;
  gint tfields;

  gchar extname[] = "NcSNIADistCov";   /* extension name */
  GPtrArray *ttype_array = g_ptr_array_sized_new (NC_DATA_SNIA_COV_TOTAL_LENGTH);
  GPtrArray *tform_array = g_ptr_array_sized_new (NC_DATA_SNIA_COV_TOTAL_LENGTH);
  GPtrArray *tunit_array = g_ptr_array_sized_new (NC_DATA_SNIA_COV_TOTAL_LENGTH);

  g_ptr_array_set_size (ttype_array, NC_DATA_SNIA_COV_TOTAL_LENGTH);
  g_ptr_array_set_size (tform_array, NC_DATA_SNIA_COV_TOTAL_LENGTH);
  g_ptr_array_set_size (tunit_array, NC_DATA_SNIA_COV_TOTAL_LENGTH);
  
  g_ptr_array_set_free_func (tform_array, g_free);

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_ZCMB)              = "Z_CMB";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_ZCMB)              = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_ZCMB)              = "CMB FRAME REDSHIFT";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_ZHE)               = "Z_HE";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_ZHE)               = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_ZHE)               = "SUN FRAME REDSHIFT";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_SIGMA_Z)           = "SIGMA_Z";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_SIGMA_Z)           = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_SIGMA_Z)           = "Z STANDARD DEVIATION";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_MAG)               = "MAG";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_MAG)               = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_MAG)               = "MAGNITUDE";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_SIGMA_MAG)         = "SIGMA_MAG";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_SIGMA_MAG)         = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_SIGMA_MAG)         = "MAGNITUDE STANDARD DEVIATION";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_WIDTH)             = "WIDTH";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_WIDTH)             = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_WIDTH)             = "WIDTH";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_SIGMA_WIDTH)       = "SIGMA_WIDTH";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_SIGMA_WIDTH)       = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_SIGMA_WIDTH)       = "WIDTH STANDARD DEVIATION";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_COLOUR)            = "COLOUR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_COLOUR)            = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_COLOUR)            = "COLOUR";
  
  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_SIGMA_COLOUR)      = "SIGMA_COLOUR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_SIGMA_COLOUR)      = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_SIGMA_COLOUR)      = "COLOUR STANDARD DEVIATION";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_THIRDPAR)          = "THIRDPAR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_THIRDPAR)          = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_THIRDPAR)          = "THIRDPAR";
  
  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_SIGMA_THIRDPAR)    = "SIGMA_THIRDPAR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_SIGMA_THIRDPAR)    = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_SIGMA_THIRDPAR)    = "THIRDPAR STANDARD DEVIATION";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_DIAG_MAG_WIDTH)    = "DIAG_MAG_WIDTH";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_DIAG_MAG_WIDTH)    = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_DIAG_MAG_WIDTH)    = "DIAGONAL MAG WIDTH VARIANCE";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_DIAG_MAG_COLOUR)   = "DIAG_MAG_COLOUR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_DIAG_MAG_COLOUR)   = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_DIAG_MAG_COLOUR)   = "DIAGONAL MAG COLOUR VARIANCE";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_DIAG_WIDTH_COLOUR) = "DIAG_WIDTH_COLOUR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_DIAG_WIDTH_COLOUR) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_DIAG_WIDTH_COLOUR) = "DIAGONAL WIDTH COLOUR VARIANCE";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_ABSMAG_SET)        = "ABSMAG_SET";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_ABSMAG_SET)        = g_strdup ("1J");
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_ABSMAG_SET)        = "ABSOLUTE MAG SET";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_VAR_MAG)           = "VAR_MAG";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_VAR_MAG)           = g_strdup_printf ("%uD", snia_cov->mu_len);
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_VAR_MAG)           = "MAG VARIANCE";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_VAR_WIDTH)         = "VAR_WIDTH";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_VAR_WIDTH)         = g_strdup_printf ("%uD", snia_cov->mu_len);
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_VAR_WIDTH)         = "WIDTH VARIANCE";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_VAR_COLOUR)        = "VAR_COLOUR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_VAR_COLOUR)        = g_strdup_printf ("%uD", snia_cov->mu_len);
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_VAR_COLOUR)        = "COLOUR VARIANCE";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_VAR_MAG_WIDTH)     = "VAR_MAG_WIDTH";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_VAR_MAG_WIDTH)     = g_strdup_printf ("%uD", snia_cov->mu_len);
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_VAR_MAG_WIDTH)     = "MAG WIDTH VARIANCE";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_VAR_MAG_COLOUR)    = "VAR_MAG_COLOUR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_VAR_MAG_COLOUR)    = g_strdup_printf ("%uD", snia_cov->mu_len);
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_VAR_MAG_COLOUR)    = "MAG COLOUR VARIANCE";

  g_ptr_array_index (ttype_array, NC_DATA_SNIA_COV_VAR_WIDTH_COLOUR)  = "VAR_WIDTH_COLOUR";
  g_ptr_array_index (tform_array, NC_DATA_SNIA_COV_VAR_WIDTH_COLOUR)  = g_strdup_printf ("%uD", snia_cov->mu_len);
  g_ptr_array_index (tunit_array, NC_DATA_SNIA_COV_VAR_WIDTH_COLOUR)  = "WIDTH COLOUR VARIANCE";
  
  tfields = ttype_array->len;

  /* initialize status before calling fitsio routines */
  status = 0;

  if (overwrite && g_file_test (filename, G_FILE_TEST_EXISTS))
    g_unlink (filename);

  /* create new FITS file */
  if (fits_create_file (&fptr, filename, &status))
    NC_FITS_ERROR (status);

  /* append a new empty binary table onto the FITS file */
  if (fits_create_tbl (fptr, BINARY_TBL, snia_cov->mu_len, tfields, (gchar **)ttype_array->pdata, (gchar **)tform_array->pdata,
                       (gchar **)tunit_array->pdata, extname, &status) )
    NC_FITS_ERROR (status);

  {
    NcmVector *data[NC_DATA_SNIA_COV_LENGTH];
    gint i;
    
    data[NC_DATA_SNIA_COV_ZCMB]              = snia_cov->z_cmb;
    data[NC_DATA_SNIA_COV_ZHE]               = snia_cov->z_he;
    data[NC_DATA_SNIA_COV_SIGMA_Z]           = snia_cov->sigma_z;
    data[NC_DATA_SNIA_COV_MAG]               = snia_cov->mag;
    data[NC_DATA_SNIA_COV_SIGMA_MAG]         = snia_cov->sigma_mag;
    data[NC_DATA_SNIA_COV_WIDTH]             = snia_cov->width;
    data[NC_DATA_SNIA_COV_SIGMA_WIDTH]       = snia_cov->sigma_width;
    data[NC_DATA_SNIA_COV_COLOUR]            = snia_cov->colour;
    data[NC_DATA_SNIA_COV_SIGMA_COLOUR]      = snia_cov->sigma_colour;
    data[NC_DATA_SNIA_COV_THIRDPAR]          = snia_cov->thirdpar;
    data[NC_DATA_SNIA_COV_SIGMA_THIRDPAR]    = snia_cov->sigma_thirdpar;
    data[NC_DATA_SNIA_COV_DIAG_MAG_WIDTH]    = snia_cov->diag_mag_width;
    data[NC_DATA_SNIA_COV_DIAG_MAG_COLOUR]   = snia_cov->diag_mag_colour;
    data[NC_DATA_SNIA_COV_DIAG_WIDTH_COLOUR] = snia_cov->diag_width_colour;

    for (i = 0; i < NC_DATA_SNIA_COV_LENGTH; i++)
    {
      if (fits_write_col_dbl (fptr, i + 1, 1, 1, snia_cov->mu_len, ncm_vector_ptr (data[i], 0), &status))
        NC_FITS_ERROR(status);
    }
  }

  if (fits_write_col_uint (fptr, NC_DATA_SNIA_COV_ABSMAG_SET + 1, 1, 1, 
                           snia_cov->mu_len, &g_array_index (snia_cov->dataset, guint32, 0), &status))
    NC_FITS_ERROR(status);

  {
    NcmMatrix *data[NC_DATA_SNIA_COV_TOTAL_LENGTH];
    data[NC_DATA_SNIA_COV_VAR_MAG]          = snia_cov->var_mag;
    data[NC_DATA_SNIA_COV_VAR_WIDTH]        = snia_cov->var_width;
    data[NC_DATA_SNIA_COV_VAR_COLOUR]       = snia_cov->var_colour;
    data[NC_DATA_SNIA_COV_VAR_MAG_WIDTH]    = snia_cov->var_mag_width;
    data[NC_DATA_SNIA_COV_VAR_MAG_COLOUR]   = snia_cov->var_mag_colour;
    data[NC_DATA_SNIA_COV_VAR_WIDTH_COLOUR] = snia_cov->var_width_colour;
    
    gint i;
    for (i = NC_DATA_SNIA_COV_VAR_MAG; i < NC_DATA_SNIA_COV_TOTAL_LENGTH; i++)
    {
      if (fits_write_col_dbl (fptr, i + 1, 1, 1, snia_cov->mu_len * snia_cov->mu_len, 
                          ncm_matrix_ptr (data[i], 0, 0), &status))
        NC_FITS_ERROR(status);
    }
  }

  if ( fits_close_file(fptr, &status) )
    NC_FITS_ERROR(status);

  g_ptr_array_unref (ttype_array);
  g_ptr_array_unref (tform_array);
  g_ptr_array_unref (tunit_array);

  return; 
}

#endif /* NUMCOSMO_HAVE_CFITSIO */
