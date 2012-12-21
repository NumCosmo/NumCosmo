/***************************************************************************
 *            nc_snia_dist_cov.c
 *
 *  Mon December 03 19:34:29 2012
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
 * SECTION:nc_snia_dist_cov
 * @title: Supernovae Distance Covariance
 * @short_description: Calculates the covariance between distance estimates
 *
 * This object implements the calculation necessary to make a statistical
 * analysis using data from <link linkend="XConley2011">Conley et al. (2011)</link>
 * and <link linkend="XSullivan2011">Sullivan et al. (2011)</link>.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_snia_dist_cov.h"
#include "math/ncm_cfg.h"

#include <glib/gstdio.h>
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

#define VECTOR  (NCM_MODEL (dcov)->params)
#define ALPHA   (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_ALPHA))
#define BETA    (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_BETA))
#define ABSMAG1 (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_M1))
#define ABSMAG2 (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_M2))

enum
{
  PROP_0,
  PROP_DIST,
  PROP_MU_LEN,
  PROP_SIGMA_PECZ,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcSNIADistCov, nc_snia_dist_cov, NCM_TYPE_MODEL);

static void
nc_snia_dist_cov_init (NcSNIADistCov *dcov)
{
  dcov->dist              = NULL;
  dcov->mu_len            = 0;
  
  dcov->z_cmb             = NULL;
  dcov->z_he              = NULL;

  dcov->mag               = NULL;
  dcov->width             = NULL;
  dcov->colour            = NULL;
  dcov->thirdpar          = NULL;

  dcov->sigma_z           = NULL;
  dcov->sigma_mag         = NULL;
  dcov->sigma_width       = NULL;
  dcov->sigma_colour      = NULL;
  dcov->sigma_thirdpar    = NULL;

  dcov->diag_mag_width    = NULL;
  dcov->diag_mag_colour   = NULL;
  dcov->diag_width_colour = NULL;

  dcov->var_mag           = NULL;
  dcov->var_width         = NULL;
  dcov->var_colour        = NULL;
  dcov->var_mag_width     = NULL;
  dcov->var_mag_colour    = NULL;
  dcov->var_width_colour  = NULL;

  dcov->sigma_pecz        = 0.0;
  dcov->dataset           = NULL;
  dcov->sigma_int         = NULL;  
}

static void
nc_snia_dist_cov_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_snia_dist_cov_parent_class)->constructed (object);
  {
    NcmModel *model = NCM_MODEL (object);
    guint sigma_int_len = ncm_model_vparam_len (model, 0);

    switch (sigma_int_len)
    {
      case 4:
        ncm_model_orig_vparam_set (model, 0, 0, 0.0675);
        ncm_model_orig_vparam_set (model, 0, 1, 0.1133);
        ncm_model_orig_vparam_set (model, 0, 2, 0.0815);
        ncm_model_orig_vparam_set (model, 0, 3, 0.0989);
        break;
    }
  }
}

static void
_nc_snia_dist_cov_set_property (GObject * object, guint prop_id, const GValue * value, GParamSpec * pspec)
{
  NcSNIADistCov *dcov = NC_SNIA_DIST_COV (object);
  g_return_if_fail (NC_IS_SNIA_DIST_COV (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_distance_clear (&dcov->dist);
      dcov->dist = g_value_dup_object (value);
      break;
    case PROP_MU_LEN:
      nc_snia_dist_cov_set_size (dcov, g_value_get_uint (value));
      break;
    case PROP_SIGMA_PECZ:
      dcov->sigma_pecz = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_snia_dist_cov_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcSNIADistCov *dcov = NC_SNIA_DIST_COV (object);
  g_return_if_fail (NC_IS_SNIA_DIST_COV (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_set_object (value, dcov->dist);
      break;
    case PROP_MU_LEN:
      g_value_set_uint (value, dcov->mu_len);
      break;
    case PROP_SIGMA_PECZ:
      g_value_set_double (value, dcov->sigma_pecz);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_snia_dist_cov_dispose (GObject *object)
{
  NcSNIADistCov *dcov = NC_SNIA_DIST_COV (object);

  nc_snia_dist_cov_set_size (dcov, 0);
  nc_distance_clear (&dcov->dist);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_snia_dist_cov_parent_class)->dispose (object);
}

static void
nc_snia_dist_cov_finalize (GObject *object)
{
  
  /* Chain up : end */
  G_OBJECT_CLASS (nc_snia_dist_cov_parent_class)->finalize (object);
}

gint32 NC_SNIA_DIST_COV_ID = -1;

static void
nc_snia_dist_cov_class_init (NcSNIADistCovClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  object_class->constructed  = &nc_snia_dist_cov_constructed;
  object_class->set_property = &ncm_model_class_set_property;
  object_class->get_property = &ncm_model_class_get_property;  
  object_class->dispose      = &nc_snia_dist_cov_dispose;
  object_class->finalize     = &nc_snia_dist_cov_finalize;

  model_class->set_property = &_nc_snia_dist_cov_set_property;
  model_class->get_property = &_nc_snia_dist_cov_get_property;

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_MU_LEN,
                                   g_param_spec_uint ("mu-len",
                                                      NULL,
                                                      "Distance modulus length",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SIGMA_PECZ,
                                   g_param_spec_double ("sigma-pecz",
                                                      NULL,
                                                      "Error from SN Ia peculiar velocity",
                                                      0.0, 1.0e1, 5.0e-4,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT |G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_register_id (model_class);
  NC_SNIA_DIST_COV_ID = model_class->model_id;

  ncm_model_class_add_params (model_class, NC_SNIA_DIST_COV_SPARAM_LEN, NC_SNIA_DIST_COV_VPARAM_LEN, PROP_SIZE);
  ncm_model_class_set_name_nick (model_class, "Supernovae Ia Distance Covariance", "SNIaDistCov");

  ncm_model_class_set_sparam (model_class, NC_SNIA_DIST_COV_ALPHA, "alpha", "alpha",
                              -10.0, 10.0, 1.0e-1,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_ALPHA,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_SNIA_DIST_COV_BETA, "beta", "beta",
                              -10.0, 10.0, 1.0e-1,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_BETA,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_SNIA_DIST_COV_M1, "Absolute Magnitude 1", "M1",
                              -50.0, 10.0, 1.0,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_M1,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_SNIA_DIST_COV_M2, "Absolute Magnitude 2", "M2",
                              -50.0, 10.0, 1.0,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_M2,
                              NCM_PARAM_TYPE_FIXED);
  
  ncm_model_class_set_vparam (model_class, NC_SNIA_DIST_COV_SIGMA_INT, NC_SNIA_DIST_COV_SIGMA_INT_DEFAULT_LEN, 
                              "Sigma intrisic", "sigma_int",
                              0.0, 1.0e1, 1.0e-3, 
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_SIGMA_INT,
                              NCM_PARAM_TYPE_FIXED);

}

/**
 * nc_snia_dist_cov_new:
 * @dist: FIXME
 * @mu_len: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
NcSNIADistCov *
nc_snia_dist_cov_new (NcDistance *dist, guint mu_len)
{
  return g_object_new (NC_TYPE_SNIA_DIST_COV,
                       "mu-len", mu_len,
                       "dist", dist,
                       NULL);
}

/**
 * nc_snia_dist_cov_ref:
 * @dcov: FIXME
 * 
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
NcSNIADistCov *
nc_snia_dist_cov_ref (NcSNIADistCov *dcov)
{
  return g_object_ref (dcov);
}

/**
 * nc_snia_dist_cov_free:
 * @dcov: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_snia_dist_cov_free (NcSNIADistCov *dcov)
{
  g_object_unref (dcov);
}

/**
 * nc_snia_dist_cov_clear:
 * @dcov: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_snia_dist_cov_clear (NcSNIADistCov **dcov)
{
  g_clear_object (dcov);
}

void 
nc_snia_dist_cov_prepare (NcSNIADistCov *dcov, NcmMSet *mset)
{
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  nc_distance_prepare (dcov->dist, cosmo);
}

void 
nc_snia_dist_cov_prepare_if_needed (NcSNIADistCov *dcov, NcmMSet *mset)
{
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, NC_HICOSMO_ID));
  nc_distance_prepare_if_needed (dcov->dist, cosmo);
}

/**
 * nc_snia_dist_cov_set_size:
 * @dcov: FIXME
 * @mu_len: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_snia_dist_cov_set_size (NcSNIADistCov *dcov, guint mu_len)
{
  if (mu_len == 0 || mu_len != dcov->mu_len)
  {
    ncm_vector_clear (&dcov->z_cmb);
    ncm_vector_clear (&dcov->z_he);

    ncm_vector_clear (&dcov->mag);
    ncm_vector_clear (&dcov->width);
    ncm_vector_clear (&dcov->colour);
    ncm_vector_clear (&dcov->thirdpar);

    ncm_vector_clear (&dcov->sigma_z);
    ncm_vector_clear (&dcov->sigma_mag);
    ncm_vector_clear (&dcov->sigma_width);
    ncm_vector_clear (&dcov->sigma_colour);
    ncm_vector_clear (&dcov->sigma_thirdpar);

    ncm_vector_clear (&dcov->diag_mag_width);
    ncm_vector_clear (&dcov->diag_mag_colour);
    ncm_vector_clear (&dcov->diag_width_colour);
    
    ncm_matrix_clear (&dcov->var_mag);
    ncm_matrix_clear (&dcov->var_width);
    ncm_matrix_clear (&dcov->var_colour);

    ncm_matrix_clear (&dcov->var_mag_width);
    ncm_matrix_clear (&dcov->var_mag_colour);
    ncm_matrix_clear (&dcov->var_width_colour);

    if (dcov->dataset != NULL)
    {
      g_array_unref (dcov->dataset);
      dcov->dataset = NULL;
    }
  }

  if (mu_len > 0)
  {
    dcov->mu_len           = mu_len;

    dcov->z_cmb            = ncm_vector_new (mu_len);
    dcov->z_he             = ncm_vector_new (mu_len);

    dcov->mag              = ncm_vector_new (mu_len);
    dcov->width            = ncm_vector_new (mu_len);
    dcov->colour           = ncm_vector_new (mu_len);
    dcov->thirdpar         = ncm_vector_new (mu_len);

    dcov->sigma_z          = ncm_vector_new (mu_len);
    dcov->sigma_mag        = ncm_vector_new (mu_len);
    dcov->sigma_width      = ncm_vector_new (mu_len);
    dcov->sigma_colour     = ncm_vector_new (mu_len);
    dcov->sigma_thirdpar   = ncm_vector_new (mu_len);

    dcov->diag_mag_width    = ncm_vector_new (mu_len);
    dcov->diag_mag_colour   = ncm_vector_new (mu_len);
    dcov->diag_width_colour = ncm_vector_new (mu_len);

    dcov->var_mag          = ncm_matrix_new (mu_len, mu_len);
    dcov->var_width        = ncm_matrix_new (mu_len, mu_len);
    dcov->var_colour       = ncm_matrix_new (mu_len, mu_len);

    dcov->var_mag_width    = ncm_matrix_new (mu_len, mu_len);
    dcov->var_mag_colour   = ncm_matrix_new (mu_len, mu_len);
    dcov->var_width_colour = ncm_matrix_new (mu_len, mu_len);

    dcov->dataset          = g_array_sized_new (FALSE, FALSE, sizeof (guint32), mu_len);

    g_array_set_size (dcov->dataset, mu_len);
  }
}

/**
 * nc_snia_dist_cov_calc:
 * @dcov: FIXME
 * @cov: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_snia_dist_cov_calc (NcSNIADistCov *dcov, NcmMatrix *cov)
{
  NcmModel *model = NCM_MODEL (dcov);
  const gdouble alpha          = ALPHA;
  const gdouble beta           = BETA;
  const gdouble alpha2         = alpha * alpha;
  const gdouble beta2          = beta * beta;
  const gdouble two_alpha_beta = 2.0 * alpha * beta;
  const gdouble two_alpha      = 2.0 * alpha;
  const gdouble two_beta       = 2.0 * beta;
  const guint ij_len = dcov->mu_len * dcov->mu_len;
  register gint i, ij;

  for (ij = 0; ij < ij_len; ij++)
  {
    const gdouble var_mag          = ncm_matrix_fast_get (dcov->var_mag, ij);
    const gdouble var_width        = ncm_matrix_fast_get (dcov->var_width, ij);
    const gdouble var_colour       = ncm_matrix_fast_get (dcov->var_colour, ij);
    const gdouble var_mag_width    = ncm_matrix_fast_get (dcov->var_mag_width, ij);
    const gdouble var_mag_colour   = ncm_matrix_fast_get (dcov->var_mag_colour, ij);
    const gdouble var_width_colour = ncm_matrix_fast_get (dcov->var_width_colour, ij);
    
    ncm_matrix_fast_set (cov, ij, 
                         var_mag 
                         + alpha2 * var_width
                         + beta2 * var_colour
                         + two_alpha * var_mag_width
                         - two_beta * var_mag_colour
                         - two_alpha_beta * var_width_colour
                         );
  }
  for (i = 0; i < dcov->mu_len; i++)
  {
    const gdouble zfacsq    = (5.0 / M_LN10) * (5.0 / M_LN10);
    const guint dset_id     = g_array_index (dcov->dataset, guint32, i);
    const gdouble sigma_int = ncm_model_orig_vparam_get (model, 0, dset_id);
    const gdouble var_int   = sigma_int * sigma_int; 
    const gdouble z_cmb     = ncm_vector_get (dcov->z_cmb, i);
    const gdouble sigma_z   = ncm_vector_get (dcov->sigma_z, i);
    const gdouble emptyfac  = (1.0 + z_cmb) / (z_cmb * (1.0 + 0.5 * z_cmb));
    const gdouble var_pecz  = dcov->sigma_pecz * dcov->sigma_pecz;
    const gdouble var_z_tot = (var_pecz + sigma_z * sigma_z) * zfacsq * emptyfac * emptyfac; 

    const gdouble sigma_mag        =  ncm_vector_get (dcov->sigma_mag, i);
    const gdouble sigma_width      =  ncm_vector_get (dcov->sigma_width, i);
    const gdouble sigma_colour     =  ncm_vector_get (dcov->sigma_colour, i);    
    const gdouble var_mag          =  sigma_mag * sigma_mag;
    const gdouble var_width        =  alpha2 * sigma_width * sigma_width;
    const gdouble var_colour       =  beta2 * sigma_colour * sigma_colour;
    const gdouble var_mag_width    =  two_alpha * ncm_vector_get (dcov->diag_mag_width, i);
    const gdouble var_mag_colour   = -two_beta * ncm_vector_get (dcov->diag_mag_colour, i);
    const gdouble var_width_colour = -two_alpha_beta * ncm_vector_get (dcov->diag_width_colour, i);

    const gdouble var_tot = var_mag + var_z_tot + var_int + var_width + 
      var_colour + var_mag_width + var_mag_colour + var_width_colour;
    ncm_matrix_set (cov, i, i, ncm_matrix_get (cov, i, i) + var_tot);
  }
}

/**
 * nc_snia_dist_cov_mean:
 * @dcov: FIXME
 * @cosmo: FIXME
 * @y: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_snia_dist_cov_mean (NcSNIADistCov *dcov, NcHICosmo *cosmo, NcmVector *y)
{
  const gdouble alpha = ALPHA;
  const gdouble beta  = BETA;
  const gdouble DH    = nc_distance_hubble (dcov->dist, cosmo);
  const gdouble Mcal1 = ABSMAG1 + 5.0 * log10 (DH);
  const gdouble Mcal2 = ABSMAG2 + 5.0 * log10 (DH);
  gint i;

  for (i = 0; i < dcov->mu_len; i++)
  {
    const gdouble z_cmb    = ncm_vector_get (dcov->z_cmb, i);
    const gdouble z_he     = ncm_vector_get (dcov->z_he, i);
    const gdouble width    = ncm_vector_get (dcov->width, i);
    const gdouble colour   = ncm_vector_get (dcov->colour, i);
    const gdouble thirdpar = ncm_vector_get (dcov->thirdpar, i);
    const gdouble mu       = nc_distance_modulus_hef (dcov->dist, cosmo, z_he, z_cmb);
    const gdouble mag_th   = mu - alpha * (width - 1.0) + beta * colour + ((thirdpar < 10.0) ? Mcal1 : Mcal2);
    const gdouble y_i      = mag_th;

    ncm_vector_set (y, i, y_i);
  }
  
}

static void _nc_snia_dist_cov_load_snia_data (NcSNIADistCov *dcov, const gchar *filename);
static void _nc_snia_dist_cov_load_matrix (const gchar *filename, NcmMatrix *data);

/**
 * nc_snia_dist_cov_load_txt:
 * @dcov: FIXME
 * @filename: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_snia_dist_cov_load_txt (NcSNIADistCov *dcov, const gchar *filename)
{
  GKeyFile *snia_keyfile = g_key_file_new ();
  GError *error   = NULL;
  gchar *datafile = NULL;
  guint64 mu_len;

  if (!g_key_file_load_from_file (snia_keyfile, filename, G_KEY_FILE_NONE, &error))
    g_error ("nc_snia_dist_cov_load: invalid configuration: %s\n  %s\n", 
             filename, 
             error->message);

  if (!g_key_file_has_key (snia_keyfile, 
                           NC_SNIA_DIST_COV_DATA_GROUP,
                           NC_SNIA_DIST_COV_DATA_KEY,
                           &error))
    g_error ("nc_snia_dist_cov_load: invalid %s key file, it must define at least the data file key ["NC_SNIA_DIST_COV_DATA_KEY"]", filename);
  if (!g_key_file_has_key (snia_keyfile, 
                           NC_SNIA_DIST_COV_DATA_GROUP,
                           NC_SNIA_DIST_COV_DATA_LEN_KEY,
                           &error))
    g_error ("nc_snia_dist_cov_load: invalid %s key file, it must define at least the data length file key ["NC_SNIA_DIST_COV_DATA_LEN_KEY"]", filename);

  datafile = g_key_file_get_string (snia_keyfile, 
                                    NC_SNIA_DIST_COV_DATA_GROUP,
                                    NC_SNIA_DIST_COV_DATA_KEY,
                                    &error);

  mu_len = g_key_file_get_uint64 (snia_keyfile, 
                                  NC_SNIA_DIST_COV_DATA_GROUP,
                                  NC_SNIA_DIST_COV_DATA_LEN_KEY,
                                  &error);

  nc_snia_dist_cov_set_size (dcov, mu_len);
  
  _nc_snia_dist_cov_load_snia_data (dcov, datafile);

  /* Get magnitude cov matrix */
  
  if (g_key_file_has_key (snia_keyfile, 
                          NC_SNIA_DIST_COV_DATA_GROUP,
                          NC_SNIA_DIST_COV_MAG_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_SNIA_DIST_COV_DATA_GROUP,
                                      NC_SNIA_DIST_COV_MAG_KEY,
                                      &error);
    
    _nc_snia_dist_cov_load_matrix (datafile, dcov->var_mag);
  }

  /* Get width cov matrix */

  if (g_key_file_has_key (snia_keyfile, 
                          NC_SNIA_DIST_COV_DATA_GROUP,
                          NC_SNIA_DIST_COV_WIDTH_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_SNIA_DIST_COV_DATA_GROUP,
                                      NC_SNIA_DIST_COV_WIDTH_KEY,
                                      &error);
    
    _nc_snia_dist_cov_load_matrix (datafile, dcov->var_width);
  }

  /* Get colour cov matrix */

  if (g_key_file_has_key (snia_keyfile, 
                          NC_SNIA_DIST_COV_DATA_GROUP,
                          NC_SNIA_DIST_COV_COLOUR_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_SNIA_DIST_COV_DATA_GROUP,
                                      NC_SNIA_DIST_COV_COLOUR_KEY,
                                      &error);
    
    _nc_snia_dist_cov_load_matrix (datafile, dcov->var_colour);
  }

  /* Get magnitude-width cov matrix */

  if (g_key_file_has_key (snia_keyfile, 
                          NC_SNIA_DIST_COV_DATA_GROUP,
                          NC_SNIA_DIST_COV_MAG_WIDTH_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_SNIA_DIST_COV_DATA_GROUP,
                                      NC_SNIA_DIST_COV_MAG_WIDTH_KEY,
                                      &error);
    
    _nc_snia_dist_cov_load_matrix (datafile, dcov->var_mag_width);
  }

  /* Get magnitude-colour cov matrix */

  if (g_key_file_has_key (snia_keyfile, 
                          NC_SNIA_DIST_COV_DATA_GROUP,
                          NC_SNIA_DIST_COV_MAG_COLOUR_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_SNIA_DIST_COV_DATA_GROUP,
                                      NC_SNIA_DIST_COV_MAG_COLOUR_KEY,
                                      &error);
    
    _nc_snia_dist_cov_load_matrix (datafile, dcov->var_mag_colour);
  }

  /* Get width-colour cov matrix */

  if (g_key_file_has_key (snia_keyfile, 
                          NC_SNIA_DIST_COV_DATA_GROUP,
                          NC_SNIA_DIST_COV_WIDTH_COLOUR_KEY,
                          &error))
  {
    datafile = g_key_file_get_string (snia_keyfile, 
                                      NC_SNIA_DIST_COV_DATA_GROUP,
                                      NC_SNIA_DIST_COV_WIDTH_COLOUR_KEY,
                                      &error);
    
    _nc_snia_dist_cov_load_matrix (datafile, dcov->var_width_colour);
  }

  g_key_file_free (snia_keyfile);
}

static void 
_nc_snia_dist_cov_load_snia_data (NcSNIADistCov *dcov, const gchar *filename)
{
  GArray *dataset = dcov->dataset;
  gchar *line = NULL;
  gsize len = 0, tpos = 0;
  GError *error = NULL;
  GIOChannel *file = g_io_channel_new_file (filename, "r", &error);
  if (file == NULL)
    g_error ("_nc_snia_dist_cov_load_snia_data: cannot open file %s: %s", 
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
          if (n - 1 - pad_end == NC_SNIA_DIST_COV_LENGTH)
            has_dataset = FALSE;
          else if (n - 2 - pad_end == NC_SNIA_DIST_COV_LENGTH)
            has_dataset = TRUE;
          else
            g_error ("_nc_snia_dist_cov_load_snia_data: data file [%s] must have %d or %d columns, it has %u", 
                     filename, NC_SNIA_DIST_COV_LENGTH + 1, NC_SNIA_DIST_COV_LENGTH + 2, n - pad_end);
        }
        else if (n != g_strv_length (itens))
          g_error ("_nc_snia_dist_cov_load_snia_data: data file [%s] has different number of columns in different rows [%u]", filename, nrow);
        
        if (nrow >= dcov->mu_len)
          g_error ("_nc_snia_dist_cov_load_snia_data: cannot load data file [%s] expected nrows %u obtained >%u\n", 
                   filename, dcov->mu_len, nrow);

        {
          NcmVector *data[NC_SNIA_DIST_COV_LENGTH];
          data[NC_SNIA_DIST_COV_ZCMB]              = dcov->z_cmb;
          data[NC_SNIA_DIST_COV_ZHE]               = dcov->z_he;
          data[NC_SNIA_DIST_COV_SIGMA_Z]           = dcov->sigma_z;
          data[NC_SNIA_DIST_COV_MAG]               = dcov->mag;
          data[NC_SNIA_DIST_COV_SIGMA_MAG]         = dcov->sigma_mag;
          data[NC_SNIA_DIST_COV_WIDTH]             = dcov->width;
          data[NC_SNIA_DIST_COV_SIGMA_WIDTH]       = dcov->sigma_width;
          data[NC_SNIA_DIST_COV_COLOUR]            = dcov->colour;
          data[NC_SNIA_DIST_COV_SIGMA_COLOUR]      = dcov->sigma_colour;
          data[NC_SNIA_DIST_COV_THIRDPAR]          = dcov->thirdpar;
          data[NC_SNIA_DIST_COV_SIGMA_THIRDPAR]    = dcov->sigma_thirdpar;
          data[NC_SNIA_DIST_COV_DIAG_MAG_WIDTH]    = dcov->diag_mag_width;
          data[NC_SNIA_DIST_COV_DIAG_MAG_COLOUR]   = dcov->diag_mag_colour;
          data[NC_SNIA_DIST_COV_DIAG_WIDTH_COLOUR] = dcov->diag_width_colour;

          for (i = 0; i < NC_SNIA_DIST_COV_LENGTH; i++)
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
    if (nrow != dcov->mu_len)
      g_error ("_nc_snia_dist_cov_load_snia_data: cannot load data file [%s] expected nrows %u obtained %u\n", 
               filename, dcov->mu_len, nrow);

    if (max_dset_id + 1 != ncm_model_vparam_len (NCM_MODEL (dcov), 0))
      g_error ("_nc_snia_dist_cov_load_snia_data: model has %u different dataset but datafile [%s] has %d", 
               ncm_model_vparam_len (NCM_MODEL (dcov), 0), filename, max_dset_id + 1);

    g_regex_unref (comment_line);
  }
  
  g_io_channel_unref (file);
}

static void 
_nc_snia_dist_cov_load_matrix (const gchar *filename, NcmMatrix *data)
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
    g_error ("_nc_snia_dist_cov_load_matrix: cannot open file %s: %s", 
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
    g_error ("_nc_snia_dist_cov_load_matrix: expected a %zux%zu matrix but got %"G_GINT64_MODIFIER"dx%"G_GINT64_MODIFIER"d", 
             NCM_MATRIX_NROWS (data), NCM_MATRIX_NROWS (data), matrix_len, matrix_len);
  
  if (matrix_len * matrix_len != itens_len)
    g_error ("_nc_snia_dist_cov_load_matrix: matrix header say %"G_GINT64_MODIFIER"d length but got %u\n", 
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
 * nc_snia_dist_cov_load:
 * @dcov: FIXME
 * @filename: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_snia_dist_cov_load (NcSNIADistCov *dcov, const gchar *filename)
{
  gint status, hdutype;
  fitsfile *fptr;

  status = 0;

  if (filename == NULL)
    g_error ("nc_snia_dist_cov_load: null filename");

  if (fits_open_file (&fptr, filename, READONLY, &status))
    NC_FITS_ERROR (status);

  if (fits_movabs_hdu (fptr, 2, &hdutype, &status))
    NC_FITS_ERROR (status);

  if (hdutype != BINARY_TBL)
    g_error ("nc_snia_dist_cov_load: NcSNIADistCov catalog is not binary.");

  {
    glong nrows;
    if (fits_get_num_rows (fptr, &nrows, &status))
      NC_FITS_ERROR (status);
    nc_snia_dist_cov_set_size (dcov, nrows);
  }

  {
    NcmVector *data[NC_SNIA_DIST_COV_LENGTH];
    gint i;
    
    data[NC_SNIA_DIST_COV_ZCMB]              = dcov->z_cmb;
    data[NC_SNIA_DIST_COV_ZHE]               = dcov->z_he;
    data[NC_SNIA_DIST_COV_SIGMA_Z]           = dcov->sigma_z;
    data[NC_SNIA_DIST_COV_MAG]               = dcov->mag;
    data[NC_SNIA_DIST_COV_SIGMA_MAG]         = dcov->sigma_mag;
    data[NC_SNIA_DIST_COV_WIDTH]             = dcov->width;
    data[NC_SNIA_DIST_COV_SIGMA_WIDTH]       = dcov->sigma_width;
    data[NC_SNIA_DIST_COV_COLOUR]            = dcov->colour;
    data[NC_SNIA_DIST_COV_SIGMA_COLOUR]      = dcov->sigma_colour;
    data[NC_SNIA_DIST_COV_THIRDPAR]          = dcov->thirdpar;
    data[NC_SNIA_DIST_COV_SIGMA_THIRDPAR]    = dcov->sigma_thirdpar;
    data[NC_SNIA_DIST_COV_DIAG_MAG_WIDTH]    = dcov->diag_mag_width;
    data[NC_SNIA_DIST_COV_DIAG_MAG_COLOUR]   = dcov->diag_mag_colour;
    data[NC_SNIA_DIST_COV_DIAG_WIDTH_COLOUR] = dcov->diag_width_colour;

    for (i = 0; i < NC_SNIA_DIST_COV_LENGTH; i++)
    {
      if (fits_read_col_dbl (fptr, i + 1, 1, 1, dcov->mu_len, GSL_NAN,
                             ncm_vector_ptr (data[i], 0), NULL, 
                             &status))
        NC_FITS_ERROR(status);
    }
  }

  if (fits_read_col_uint (fptr, NC_SNIA_DIST_COV_ABSMAG_SET + 1, 1, 1, 
                          dcov->mu_len, 
                          0, &g_array_index (dcov->dataset, guint32, 0), 
                          NULL, &status))
    NC_FITS_ERROR(status);

  {
    NcmMatrix *data[NC_SNIA_DIST_COV_TOTAL_LENGTH];
    data[NC_SNIA_DIST_COV_VAR_MAG]          = dcov->var_mag;
    data[NC_SNIA_DIST_COV_VAR_WIDTH]        = dcov->var_width;
    data[NC_SNIA_DIST_COV_VAR_COLOUR]       = dcov->var_colour;
    data[NC_SNIA_DIST_COV_VAR_MAG_WIDTH]    = dcov->var_mag_width;
    data[NC_SNIA_DIST_COV_VAR_MAG_COLOUR]   = dcov->var_mag_colour;
    data[NC_SNIA_DIST_COV_VAR_WIDTH_COLOUR] = dcov->var_width_colour;
    
    gint i;
    for (i = NC_SNIA_DIST_COV_VAR_MAG; i < NC_SNIA_DIST_COV_TOTAL_LENGTH; i++)
    {
      if (fits_read_col_dbl (fptr, i + 1, 1, 1, dcov->mu_len * dcov->mu_len, 
                             GSL_NAN, ncm_matrix_ptr (data[i], 0, 0), 
                             NULL, &status))
        NC_FITS_ERROR(status);
    }
  }
}

/**
 * nc_snia_dist_cov_save:
 * @dcov: FIXME
 * @filename: FIXME
 * @overwrite: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_snia_dist_cov_save (NcSNIADistCov *dcov, const gchar *filename, gboolean overwrite)
{
  /*******************************************************************/
  /* Create a binary table extension                                 */
  /*******************************************************************/
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  gint status;
  gint tfields;

  gchar extname[] = "NcSNIADistCov";   /* extension name */
  GPtrArray *ttype_array = g_ptr_array_sized_new (NC_SNIA_DIST_COV_TOTAL_LENGTH);
  GPtrArray *tform_array = g_ptr_array_sized_new (NC_SNIA_DIST_COV_TOTAL_LENGTH);
  GPtrArray *tunit_array = g_ptr_array_sized_new (NC_SNIA_DIST_COV_TOTAL_LENGTH);

  g_ptr_array_set_size (ttype_array, NC_SNIA_DIST_COV_TOTAL_LENGTH);
  g_ptr_array_set_size (tform_array, NC_SNIA_DIST_COV_TOTAL_LENGTH);
  g_ptr_array_set_size (tunit_array, NC_SNIA_DIST_COV_TOTAL_LENGTH);
  
  g_ptr_array_set_free_func (tform_array, g_free);

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_ZCMB)              = "Z_CMB";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_ZCMB)              = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_ZCMB)              = "CMB FRAME REDSHIFT";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_ZHE)               = "Z_HE";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_ZHE)               = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_ZHE)               = "SUN FRAME REDSHIFT";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_SIGMA_Z)           = "SIGMA_Z";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_SIGMA_Z)           = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_SIGMA_Z)           = "Z STANDARD DEVIATION";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_MAG)               = "MAG";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_MAG)               = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_MAG)               = "MAGNITUDE";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_SIGMA_MAG)         = "SIGMA_MAG";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_SIGMA_MAG)         = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_SIGMA_MAG)         = "MAGNITUDE STANDARD DEVIATION";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_WIDTH)             = "WIDTH";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_WIDTH)             = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_WIDTH)             = "WIDTH";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_SIGMA_WIDTH)       = "SIGMA_WIDTH";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_SIGMA_WIDTH)       = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_SIGMA_WIDTH)       = "WIDTH STANDARD DEVIATION";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_COLOUR)            = "COLOUR";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_COLOUR)            = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_COLOUR)            = "COLOUR";
  
  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_SIGMA_COLOUR)      = "SIGMA_COLOUR";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_SIGMA_COLOUR)      = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_SIGMA_COLOUR)      = "COLOUR STANDARD DEVIATION";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_THIRDPAR)          = "THIRDPAR";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_THIRDPAR)          = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_THIRDPAR)          = "THIRDPAR";
  
  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_SIGMA_THIRDPAR)    = "SIGMA_THIRDPAR";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_SIGMA_THIRDPAR)    = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_SIGMA_THIRDPAR)    = "THIRDPAR STANDARD DEVIATION";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_DIAG_MAG_WIDTH)    = "DIAG_MAG_WIDTH";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_DIAG_MAG_WIDTH)    = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_DIAG_MAG_WIDTH)    = "DIAGONAL MAG WIDTH VARIANCE";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_DIAG_MAG_COLOUR)   = "DIAG_MAG_COLOUR";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_DIAG_MAG_COLOUR)   = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_DIAG_MAG_COLOUR)   = "DIAGONAL MAG COLOUR VARIANCE";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_DIAG_WIDTH_COLOUR) = "DIAG_WIDTH_COLOUR";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_DIAG_WIDTH_COLOUR) = g_strdup ("1D");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_DIAG_WIDTH_COLOUR) = "DIAGONAL WIDTH COLOUR VARIANCE";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_ABSMAG_SET)        = "ABSMAG_SET";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_ABSMAG_SET)        = g_strdup ("1J");
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_ABSMAG_SET)        = "ABSOLUTE MAG SET";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_VAR_MAG)           = "VAR_MAG";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_VAR_MAG)           = g_strdup_printf ("%uD", dcov->mu_len);
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_VAR_MAG)           = "MAG VARIANCE";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_VAR_WIDTH)         = "VAR_WIDTH";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_VAR_WIDTH)         = g_strdup_printf ("%uD", dcov->mu_len);
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_VAR_WIDTH)         = "WIDTH VARIANCE";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_VAR_COLOUR)        = "VAR_COLOUR";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_VAR_COLOUR)        = g_strdup_printf ("%uD", dcov->mu_len);
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_VAR_COLOUR)        = "COLOUR VARIANCE";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_VAR_MAG_WIDTH)     = "VAR_MAG_WIDTH";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_VAR_MAG_WIDTH)     = g_strdup_printf ("%uD", dcov->mu_len);
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_VAR_MAG_WIDTH)     = "MAG WIDTH VARIANCE";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_VAR_MAG_COLOUR)    = "VAR_MAG_COLOUR";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_VAR_MAG_COLOUR)    = g_strdup_printf ("%uD", dcov->mu_len);
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_VAR_MAG_COLOUR)    = "MAG COLOUR VARIANCE";

  g_ptr_array_index (ttype_array, NC_SNIA_DIST_COV_VAR_WIDTH_COLOUR)  = "VAR_WIDTH_COLOUR";
  g_ptr_array_index (tform_array, NC_SNIA_DIST_COV_VAR_WIDTH_COLOUR)  = g_strdup_printf ("%uD", dcov->mu_len);
  g_ptr_array_index (tunit_array, NC_SNIA_DIST_COV_VAR_WIDTH_COLOUR)  = "WIDTH COLOUR VARIANCE";
  
  tfields = ttype_array->len;

  /* initialize status before calling fitsio routines */
  status = 0;

  if (overwrite && g_file_test (filename, G_FILE_TEST_EXISTS))
    g_unlink (filename);

  /* create new FITS file */
  if (fits_create_file (&fptr, filename, &status))
    NC_FITS_ERROR (status);

  /* append a new empty binary table onto the FITS file */
  if (fits_create_tbl (fptr, BINARY_TBL, dcov->mu_len, tfields, (gchar **)ttype_array->pdata, (gchar **)tform_array->pdata,
                       (gchar **)tunit_array->pdata, extname, &status) )
    NC_FITS_ERROR (status);

  if (FALSE) /* Not necessarily a good idea, but think more about it later. */
  {
    gchar *dcov_ser = ncm_cfg_serialize_to_string (G_OBJECT (dcov), FALSE);

    if (fits_write_key_longstr (fptr, "DCOV", dcov_ser, "Serialized distance covariance object", &status))
      NC_FITS_ERROR(status);

    g_free (dcov_ser);
  }

  {
    NcmVector *data[NC_SNIA_DIST_COV_LENGTH];
    gint i;
    
    data[NC_SNIA_DIST_COV_ZCMB]              = dcov->z_cmb;
    data[NC_SNIA_DIST_COV_ZHE]               = dcov->z_he;
    data[NC_SNIA_DIST_COV_SIGMA_Z]           = dcov->sigma_z;
    data[NC_SNIA_DIST_COV_MAG]               = dcov->mag;
    data[NC_SNIA_DIST_COV_SIGMA_MAG]         = dcov->sigma_mag;
    data[NC_SNIA_DIST_COV_WIDTH]             = dcov->width;
    data[NC_SNIA_DIST_COV_SIGMA_WIDTH]       = dcov->sigma_width;
    data[NC_SNIA_DIST_COV_COLOUR]            = dcov->colour;
    data[NC_SNIA_DIST_COV_SIGMA_COLOUR]      = dcov->sigma_colour;
    data[NC_SNIA_DIST_COV_THIRDPAR]          = dcov->thirdpar;
    data[NC_SNIA_DIST_COV_SIGMA_THIRDPAR]    = dcov->sigma_thirdpar;
    data[NC_SNIA_DIST_COV_DIAG_MAG_WIDTH]    = dcov->diag_mag_width;
    data[NC_SNIA_DIST_COV_DIAG_MAG_COLOUR]   = dcov->diag_mag_colour;
    data[NC_SNIA_DIST_COV_DIAG_WIDTH_COLOUR] = dcov->diag_width_colour;

    for (i = 0; i < NC_SNIA_DIST_COV_LENGTH; i++)
    {
      if (fits_write_col_dbl (fptr, i + 1, 1, 1, dcov->mu_len, ncm_vector_ptr (data[i], 0), &status))
        NC_FITS_ERROR(status);
    }
  }

  if (fits_write_col_uint (fptr, NC_SNIA_DIST_COV_ABSMAG_SET + 1, 1, 1, 
                           dcov->mu_len, &g_array_index (dcov->dataset, guint32, 0), &status))
    NC_FITS_ERROR(status);

  {
    NcmMatrix *data[NC_SNIA_DIST_COV_TOTAL_LENGTH];
    data[NC_SNIA_DIST_COV_VAR_MAG]          = dcov->var_mag;
    data[NC_SNIA_DIST_COV_VAR_WIDTH]        = dcov->var_width;
    data[NC_SNIA_DIST_COV_VAR_COLOUR]       = dcov->var_colour;
    data[NC_SNIA_DIST_COV_VAR_MAG_WIDTH]    = dcov->var_mag_width;
    data[NC_SNIA_DIST_COV_VAR_MAG_COLOUR]   = dcov->var_mag_colour;
    data[NC_SNIA_DIST_COV_VAR_WIDTH_COLOUR] = dcov->var_width_colour;
    
    gint i;
    for (i = NC_SNIA_DIST_COV_VAR_MAG; i < NC_SNIA_DIST_COV_TOTAL_LENGTH; i++)
    {
      if (fits_write_col_dbl (fptr, i + 1, 1, 1, dcov->mu_len * dcov->mu_len, 
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
