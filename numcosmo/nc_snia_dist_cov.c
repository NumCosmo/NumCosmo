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

#define VECTOR  (NCM_MODEL (dcov)->params)
#define ALPHA   (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_ALPHA))
#define BETA    (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_BETA))
#define ABSMAG1 (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_M1))
#define ABSMAG2 (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_M2))

enum
{
  PROP_0,
  PROP_DIST,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcSNIADistCov, nc_snia_dist_cov, NCM_TYPE_MODEL);

static void
nc_snia_dist_cov_init (NcSNIADistCov *dcov)
{
  dcov->dist = NULL;
}

static void
nc_snia_dist_cov_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_snia_dist_cov_parent_class)->constructed (object);
  {
    NcmModel *model = NCM_MODEL (object);
    guint sigma_int_len = ncm_model_vparam_len (model, NC_SNIA_DIST_COV_SIGMA_INT);

    switch (sigma_int_len)
    {
      case 4:
        ncm_model_orig_vparam_set (model, NC_SNIA_DIST_COV_SIGMA_INT, 0, 0.0675);
        ncm_model_orig_vparam_set (model, NC_SNIA_DIST_COV_SIGMA_INT, 1, 0.1133);
        ncm_model_orig_vparam_set (model, NC_SNIA_DIST_COV_SIGMA_INT, 2, 0.0815);
        ncm_model_orig_vparam_set (model, NC_SNIA_DIST_COV_SIGMA_INT, 3, 0.0989);
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
nc_snia_dist_cov_dispose (GObject *object)
{
  NcSNIADistCov *dcov = NC_SNIA_DIST_COV (object);
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

NCM_MSET_MODEL_REGISTER_ID (nc_snia_dist_cov, NC_TYPE_SNIA_DIST_COV);

static void
nc_snia_dist_cov_class_init (NcSNIADistCovClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  object_class->constructed  = &nc_snia_dist_cov_constructed;
  object_class->dispose      = &nc_snia_dist_cov_dispose;
  object_class->finalize     = &nc_snia_dist_cov_finalize;

  model_class->set_property = &_nc_snia_dist_cov_set_property;
  model_class->get_property = &_nc_snia_dist_cov_get_property;

  ncm_model_class_add_params (model_class, NC_SNIA_DIST_COV_SPARAM_LEN, NC_SNIA_DIST_COV_VPARAM_LEN, PROP_SIZE);
  ncm_model_class_set_name_nick (model_class, "Supernovae Ia Distance Covariance", "SNIaDistCov");
  
  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));  

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

  ncm_model_class_set_vparam (model_class, NC_SNIA_DIST_COV_MU, NC_SNIA_DIST_COV_MU_DEFAULT_LEN, 
                              "Distance modulus", "mu",
                              -50.0, 50.0, 1.0e-2, 
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_MU,
                              NCM_PARAM_TYPE_FIXED);
  
  ncm_model_class_check_params_info (model_class);

  ncm_mset_model_register_id (model_class, 
                              "NcSNIADistCov",
                              "Supernovae distance models with errors covariance.",
                              NULL);
  
}

/**
 * nc_snia_dist_cov_new:
 * @dist: FIXME
 * 
 * FIXME
 * 
 * Returns: FIXME
 */
NcSNIADistCov *
nc_snia_dist_cov_new (NcDistance *dist)
{
  return g_object_new (NC_TYPE_SNIA_DIST_COV,
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
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare (dcov->dist, cosmo);
}

void 
nc_snia_dist_cov_prepare_if_needed (NcSNIADistCov *dcov, NcmMSet *mset)
{
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (dcov->dist, cosmo);
}

/**
 * nc_snia_dist_cov_calc:
 * @dcov: FIXME
 * @snia_cov: FIXME
 * @cov: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_snia_dist_cov_calc (NcSNIADistCov *dcov, NcDataSNIACov *snia_cov, NcmMatrix *cov)
{
  NcmModel *model = NCM_MODEL (dcov);
  const gdouble alpha          = ALPHA;
  const gdouble beta           = BETA;
  const gdouble alpha2         = alpha * alpha;
  const gdouble beta2          = beta * beta;
  const gdouble two_alpha_beta = 2.0 * alpha * beta;
  const gdouble two_alpha      = 2.0 * alpha;
  const gdouble two_beta       = 2.0 * beta;
  const guint ij_len = snia_cov->mu_len * snia_cov->mu_len;
  register gint i, ij;

  g_assert (NCM_DATA (snia_cov)->init);

  if (ncm_model_vparam_len (model, NC_SNIA_DIST_COV_SIGMA_INT) > snia_cov->dataset_len)
    g_warning ("nc_snia_dist_cov_calc: model dataset is larger then the used by the data: %u > %u.", 
               ncm_model_vparam_len (model, NC_SNIA_DIST_COV_SIGMA_INT), snia_cov->dataset_len);
  else if (ncm_model_vparam_len (model, NC_SNIA_DIST_COV_SIGMA_INT) < snia_cov->dataset_len)    
    g_error ("nc_snia_dist_cov_calc: model dataset is smaller then the used by the data: %u < %u.", 
             ncm_model_vparam_len (model, NC_SNIA_DIST_COV_SIGMA_INT), snia_cov->dataset_len);
  
  for (ij = 0; ij < ij_len; ij++)
  {
    const gdouble var_mag          = ncm_matrix_fast_get (snia_cov->var_mag, ij);
    const gdouble var_width        = ncm_matrix_fast_get (snia_cov->var_width, ij);
    const gdouble var_colour       = ncm_matrix_fast_get (snia_cov->var_colour, ij);
    const gdouble var_mag_width    = ncm_matrix_fast_get (snia_cov->var_mag_width, ij);
    const gdouble var_mag_colour   = ncm_matrix_fast_get (snia_cov->var_mag_colour, ij);
    const gdouble var_width_colour = ncm_matrix_fast_get (snia_cov->var_width_colour, ij);
    
    ncm_matrix_fast_set (cov, ij, 
                         var_mag 
                         + alpha2 * var_width
                         + beta2 * var_colour
                         + two_alpha * var_mag_width
                         - two_beta * var_mag_colour
                         - two_alpha_beta * var_width_colour
                         );
  }
  for (i = 0; i < snia_cov->mu_len; i++)
  {
    const gdouble zfacsq    = (5.0 / M_LN10) * (5.0 / M_LN10);
    const guint dset_id     = g_array_index (snia_cov->dataset, guint32, i);
    const gdouble sigma_int = ncm_model_orig_vparam_get (model, NC_SNIA_DIST_COV_SIGMA_INT, dset_id);
    const gdouble var_int   = sigma_int * sigma_int; 
    const gdouble z_cmb     = ncm_vector_get (snia_cov->z_cmb, i);
    const gdouble sigma_z   = ncm_vector_get (snia_cov->sigma_z, i);
    const gdouble emptyfac  = (1.0 + z_cmb) / (z_cmb * (1.0 + 0.5 * z_cmb));
    const gdouble var_pecz  = snia_cov->sigma_pecz * snia_cov->sigma_pecz;
    const gdouble var_z_tot = (var_pecz + sigma_z * sigma_z) * zfacsq * emptyfac * emptyfac; 

    const gdouble sigma_mag        =  ncm_vector_get (snia_cov->sigma_mag, i);
    const gdouble sigma_width      =  ncm_vector_get (snia_cov->sigma_width, i);
    const gdouble sigma_colour     =  ncm_vector_get (snia_cov->sigma_colour, i);    
    const gdouble var_mag          =  sigma_mag * sigma_mag;
    const gdouble var_width        =  alpha2 * sigma_width * sigma_width;
    const gdouble var_colour       =  beta2 * sigma_colour * sigma_colour;
    const gdouble var_mag_width    =  two_alpha * ncm_vector_get (snia_cov->diag_mag_width, i);
    const gdouble var_mag_colour   = -two_beta * ncm_vector_get (snia_cov->diag_mag_colour, i);
    const gdouble var_width_colour = -two_alpha_beta * ncm_vector_get (snia_cov->diag_width_colour, i);

    const gdouble var_tot = var_mag + var_z_tot + var_int + var_width + 
      var_colour + var_mag_width + var_mag_colour + var_width_colour;
    ncm_matrix_set (cov, i, i, ncm_matrix_get (cov, i, i) + var_tot);
  }
}

/**
 * nc_snia_dist_cov_mean:
 * @dcov: FIXME
 * @cosmo: FIXME
 * @snia_cov: FIXME
 * @y: FIXME
 * 
 * FIXME
 * 
 */
void 
nc_snia_dist_cov_mean (NcSNIADistCov *dcov, NcHICosmo *cosmo, NcDataSNIACov *snia_cov, NcmVector *y)
{
  NcmModel *model = NCM_MODEL (dcov);
  guint mu_len = ncm_model_vparam_len (model, NC_SNIA_DIST_COV_MU);
  if (mu_len > 0)
  {
    if (mu_len != snia_cov->mu_len)
      g_error ("nc_snia_dist_cov_mean: number of distance modulus variables in NcSNIADistCov don't match NcDataSNIACov %u != %u", mu_len, snia_cov->mu_len);
    else
    {
      guint i;
      const gdouble alpha = ALPHA;
      const gdouble beta  = BETA;
      const gdouble DH    = nc_distance_hubble (dcov->dist, cosmo);
      const gdouble Mcal1 = ABSMAG1 + 5.0 * log10 (DH);
      const gdouble Mcal2 = ABSMAG2 + 5.0 * log10 (DH);
      
      for (i = 0; i < mu_len; i++)
      {
        const gdouble width    = ncm_vector_get (snia_cov->width, i);
        const gdouble colour   = ncm_vector_get (snia_cov->colour, i);
        const gdouble thirdpar = ncm_vector_get (snia_cov->thirdpar, i);
        const gdouble mu = ncm_model_orig_vparam_get (model, NC_SNIA_DIST_COV_MU, i);
        const gdouble mag_th = mu - alpha * (width - 1.0) + beta * colour + ((thirdpar < 10.0) ? Mcal1 : Mcal2);
        ncm_vector_set (y, i, mag_th);
      }
    }
  }
  else
  {
    const gdouble alpha = ALPHA;
    const gdouble beta  = BETA;
    const gdouble DH    = nc_distance_hubble (dcov->dist, cosmo);
    const gdouble Mcal1 = ABSMAG1 + 5.0 * log10 (DH);
    const gdouble Mcal2 = ABSMAG2 + 5.0 * log10 (DH);
    gint i;

    g_assert (NCM_DATA (snia_cov)->init);

    for (i = 0; i < snia_cov->mu_len; i++)
    {
      const gdouble z_cmb    = ncm_vector_get (snia_cov->z_cmb, i);
      const gdouble z_he     = ncm_vector_get (snia_cov->z_he, i);
      const gdouble width    = ncm_vector_get (snia_cov->width, i);
      const gdouble colour   = ncm_vector_get (snia_cov->colour, i);
      const gdouble thirdpar = ncm_vector_get (snia_cov->thirdpar, i);
      const gdouble mu       = nc_distance_modulus_hef (dcov->dist, cosmo, z_he, z_cmb);
      const gdouble mag_th   = mu - alpha * (width - 1.0) + beta * colour + ((thirdpar < 10.0) ? Mcal1 : Mcal2);
      const gdouble y_i      = mag_th;

      ncm_vector_set (y, i, y_i);
    }
  }  
}
