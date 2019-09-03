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
 * @title: NcSNIADistCov
 * @short_description: Supernovae distance covariance between distance estimates.
 *
 * This object implements the calculation necessary to make a statistical
 * analysis using data from [Conley et al. (2011)][XConley2011]
 * and [Sullivan et al. (2011)][XSullivan2011].
 *
 * Is also supports [Betoule et al. (2014)][XBetoule2014].
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_snia_dist_cov.h"
#include "math/ncm_cfg.h"

#define VECTOR       (NCM_MODEL (dcov)->params)
#define ALPHA        (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_ALPHA))
#define BETA         (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_BETA))
#define ABSMAG1      (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_M1))
#define ABSMAG2      (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_M2))
#define LNSIGMA_PECZ (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_LNSIGMA_PECZ))
#define LNSIGMA_LENS (ncm_vector_get (VECTOR, NC_SNIA_DIST_COV_LNSIGMA_LENS))

enum
{
  PROP_0,
  PROP_DIST,
  PROP_EMPTY_FAC,
  PROP_SIZE,
};

G_DEFINE_TYPE (NcSNIADistCov, nc_snia_dist_cov, NCM_TYPE_MODEL);

static void
nc_snia_dist_cov_init (NcSNIADistCov *dcov)
{
  dcov->dist             = nc_distance_new (2.0);
  dcov->var_int          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  dcov->empty_fac        = FALSE;
  dcov->cov_cpu          = NULL;
  dcov->alpha_cpu        = GSL_POSINF; 
  dcov->beta_cpu         = GSL_POSINF; 
  dcov->lnsigma_pecz_cpu = GSL_POSINF;
  dcov->lnsigma_lens_cpu = GSL_POSINF;
}

static void
nc_snia_dist_cov_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (nc_snia_dist_cov_parent_class)->constructed (object);
  {
    NcmModel *model     = NCM_MODEL (object);
    guint sigma_int_len = ncm_model_vparam_len (model, NC_SNIA_DIST_COV_LNSIGMA_INT);
    NcSNIADistCov *dcov = NC_SNIA_DIST_COV (object);
    g_array_set_size (dcov->var_int, sigma_int_len);

    switch (sigma_int_len)
    {
      case 4:
        ncm_model_orig_vparam_set (model, NC_SNIA_DIST_COV_LNSIGMA_INT, 0, log (0.0675));
        ncm_model_orig_vparam_set (model, NC_SNIA_DIST_COV_LNSIGMA_INT, 1, log (0.1133));
        ncm_model_orig_vparam_set (model, NC_SNIA_DIST_COV_LNSIGMA_INT, 2, log (0.0815));
        ncm_model_orig_vparam_set (model, NC_SNIA_DIST_COV_LNSIGMA_INT, 3, log (0.0989));
        break;
      default:
      {
        guint i;
        for (i = 0; i < sigma_int_len; i++)
        {
          ncm_model_orig_vparam_set (model, NC_SNIA_DIST_COV_LNSIGMA_INT, i, -10.0 * M_LN10);
        }
      }
      break;
    }
    dcov->cov_cpu = NULL;
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
    {
      NcDistance *dist = g_value_dup_object (value);
      if (dist != NULL)
      {
        nc_distance_clear (&dcov->dist);
        dcov->dist = dist;
      }
      break;
    }
    case PROP_EMPTY_FAC:
    {
      gboolean empty_fac = g_value_get_boolean (value);
      if ((empty_fac && !dcov->empty_fac) || (!empty_fac && dcov->empty_fac))
      {
        dcov->empty_fac = empty_fac;
        dcov->cov_cpu = NULL;
      }
      break;
    }
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
    case PROP_EMPTY_FAC:
      g_value_set_boolean (value, dcov->empty_fac);
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
  g_clear_pointer (&dcov->var_int, g_array_unref);

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

  ncm_model_class_set_name_nick (model_class, "Supernovae Ia Distance Covariance", "SNIaDistCov");
  ncm_model_class_add_params (model_class, NC_SNIA_DIST_COV_SPARAM_LEN, NC_SNIA_DIST_COV_VPARAM_LEN, PROP_SIZE);

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("dist",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_EMPTY_FAC,
                                   g_param_spec_boolean ("empty-fac",
                                                         NULL,
                                                         "Empty universe approximation factor",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_set_sparam (model_class, NC_SNIA_DIST_COV_ALPHA, "\\alpha", "alpha",
                              0.0, 5.0, 1.0e-1,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_ALPHA,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_SNIA_DIST_COV_BETA, "\\beta", "beta",
                              0.0, 5.0, 1.0e-1,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_BETA,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_SNIA_DIST_COV_M1, "\\mathcal{M}_1", "M1",
                              -30.0, -10.0, 1.0,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_M1,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_SNIA_DIST_COV_M2, "\\mathcal{M}_2", "M2",
                              -30.0, -10.0, 1.0,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_M2,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_SNIA_DIST_COV_LNSIGMA_PECZ, "\\ln(\\sigma_{\\mathrm{pecz}})", "lnsigma_pecz",
                              -20.0 * M_LN10, 5.0 * M_LN10, 1.0e-3,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_LNSIGMA_PECZ,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_sparam (model_class, NC_SNIA_DIST_COV_LNSIGMA_LENS, "\\ln(\\sigma_{\\mathrm{lens}})", "lnsigma_lens",
                              -20.0 * M_LN10, 5.0 * M_LN10, 1.0e-3,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_LNSIGMA_LENS,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_vparam (model_class, NC_SNIA_DIST_COV_LNSIGMA_INT, NC_SNIA_DIST_COV_LNSIGMA_INT_DEFAULT_LEN,
                              "\\ln(\\sigma_{\\mathrm{int}})", "lnsigma_int",
                              -20.0 * M_LN10, 5.0 * M_LN10, 1.0e-3,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_LNSIGMA_INT,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_set_vparam (model_class, NC_SNIA_DIST_COV_MU, NC_SNIA_DIST_COV_MU_DEFAULT_LEN,
                              "\\mu", "mu",
                              -50.0, 50.0, 1.0e-2,
                              NC_SNIA_DIST_COV_DEFAULT_PARAMS_ABSTOL, NC_SNIA_DIST_COV_DEFAULT_MU,
                              NCM_PARAM_TYPE_FIXED);

  ncm_model_class_check_params_info (model_class);

  ncm_mset_model_register_id (model_class,
                              "NcSNIADistCov",
                              "Supernovae distance models with errors covariance.",
                              NULL,
                              FALSE,
                              NCM_MSET_MODEL_MAIN);

}

/**
 * nc_snia_dist_cov_new:
 * @dist: a #NcDistance
 * @sigma_int_len: length of the sigma_int dataset
 *
 * FIXME
 *
 * Returns: FIXME
 */
NcSNIADistCov *
nc_snia_dist_cov_new (NcDistance *dist, guint sigma_int_len)
{
  return g_object_new (NC_TYPE_SNIA_DIST_COV,
                       "dist", dist,
                       "lnsigma_int-length", sigma_int_len,
                       NULL);
}

/**
 * nc_snia_dist_cov_ref:
 * @dcov: a #NcSNIADistCov
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
 * @dcov: a #NcSNIADistCov
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
 * @dcov: a #NcSNIADistCov
 *
 * FIXME
 *
 */
void
nc_snia_dist_cov_clear (NcSNIADistCov **dcov)
{
  g_clear_object (dcov);
}

/**
 * nc_snia_dist_cov_set_empty_fac:
 * @dcov: a #NcSNIADistCov
 * @enable: FIXME
 *
 * FIXME
 *
 */
void
nc_snia_dist_cov_set_empty_fac (NcSNIADistCov *dcov, gboolean enable)
{
  dcov->empty_fac = enable;
}

/**
 * nc_snia_dist_cov_set_dist:
 * @dcov: a #NcSNIADistCov
 * @dist: a #NcDistance
 *
 * Sets the #NcDistance object to @dist.
 *
 */
void 
nc_snia_dist_cov_set_dist (NcSNIADistCov *dcov, NcDistance *dist)
{
  g_assert (dist != NULL);
  nc_distance_clear (&dcov->dist);
  
  dcov->dist = nc_distance_ref (dist);
}

/**
 * nc_snia_dist_cov_prepare:
 * @dcov: a #NcSNIADistCov
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 */
void
nc_snia_dist_cov_prepare (NcSNIADistCov *dcov, NcmMSet *mset)
{
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare (dcov->dist, cosmo);
}

/**
 * nc_snia_dist_cov_prepare_if_needed:
 * @dcov: a #NcSNIADistCov
 * @mset: a #NcmMSet
 *
 * FIXME
 *
 */
void
nc_snia_dist_cov_prepare_if_needed (NcSNIADistCov *dcov, NcmMSet *mset)
{
  NcHICosmo *cosmo = NC_HICOSMO (ncm_mset_peek (mset, nc_hicosmo_id ()));
  nc_distance_prepare_if_needed (dcov->dist, cosmo);
}

static gdouble
_nc_snia_dist_cov_calc_empty_fac (NcSNIADistCov *dcov, gdouble z_cmb)
{
  static const gdouble zfac = (5.0 / M_LN10);
  if (dcov->empty_fac)
    return  ((1.0 + z_cmb) / (z_cmb * (1.0 + 0.5 * z_cmb))) * zfac;
  else
    return zfac / z_cmb;
}

/**
 * nc_snia_dist_cov_calc:
 * @dcov: a #NcSNIADistCov
 * @snia_cov: a #NcDataSNIACov
 * @cov: a #NcmMatrix
 *
 * FIXME
 *
 * Returns: whether the covariance was computed.
 */
gboolean
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
  const gdouble lnsigma_pecz   = LNSIGMA_PECZ;
  const gdouble lnsigma_lens   = LNSIGMA_LENS;
  const gdouble var_pecz       = exp (2.0 * lnsigma_pecz);
  const gdouble var_lens       = exp (2.0 * lnsigma_lens);
  const guint mu_len           = snia_cov->mu_len;
  gboolean needs_update        = FALSE;
  register guint i, j, ij;

  g_assert (NCM_DATA (snia_cov)->init);

  if (ncm_model_vparam_len (model, NC_SNIA_DIST_COV_LNSIGMA_INT) > snia_cov->dataset_len)
    g_warning ("nc_snia_dist_cov_calc: model dataset is larger then the used by the data: %u > %u.",
               ncm_model_vparam_len (model, NC_SNIA_DIST_COV_LNSIGMA_INT), snia_cov->dataset_len);
  else if (ncm_model_vparam_len (model, NC_SNIA_DIST_COV_LNSIGMA_INT) < snia_cov->dataset_len)
    g_error ("nc_snia_dist_cov_calc: model dataset is smaller then the used by the data: %u < %u.",
             ncm_model_vparam_len (model, NC_SNIA_DIST_COV_LNSIGMA_INT), snia_cov->dataset_len);
  ij = 0;

  for (i = 0; i < dcov->var_int->len; i++)
  {
    const gdouble var_int_i = exp (2.0 * ncm_model_orig_vparam_get (model, NC_SNIA_DIST_COV_LNSIGMA_INT, i));

    needs_update = needs_update || (g_array_index (dcov->var_int, gdouble, i) != var_int_i);
    
    g_array_index (dcov->var_int, gdouble, i) = var_int_i;
  }

  needs_update = needs_update || (cov          != dcov->cov_cpu);
  needs_update = needs_update || (alpha        != dcov->alpha_cpu);
  needs_update = needs_update || (beta         != dcov->beta_cpu);
  needs_update = needs_update || (lnsigma_pecz != dcov->lnsigma_pecz_cpu);
  needs_update = needs_update || (lnsigma_lens != dcov->lnsigma_lens_cpu);
  
  if (!needs_update)
    return FALSE;
  
  for (i = 0; i < mu_len; i++)
  {
    for (j = i; j < mu_len; j++)
    {
      const gdouble mag_mag       = ncm_vector_fast_get (snia_cov->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_MAG_MAG);
      const gdouble mag_width     = ncm_vector_fast_get (snia_cov->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_MAG_WIDTH);
      const gdouble mag_colour    = ncm_vector_fast_get (snia_cov->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_MAG_COLOUR);
      const gdouble width_width   = ncm_vector_fast_get (snia_cov->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_WIDTH_WIDTH);
      const gdouble width_colour  = ncm_vector_fast_get (snia_cov->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_WIDTH_COLOUR);
      const gdouble colour_colour = ncm_vector_fast_get (snia_cov->cov_packed, NC_DATA_SNIA_COV_ORDER_LENGTH * ij + NC_DATA_SNIA_COV_ORDER_COLOUR_COLOUR);
      ncm_matrix_set (cov, i, j,
                      mag_mag
                      + alpha2 * width_width
                      + beta2 * colour_colour
                      + two_alpha * mag_width
                      - two_beta * mag_colour
                      - two_alpha_beta * width_colour
                      );
      ij++;
    }

    {
      const guint dset_id      = g_array_index (snia_cov->dataset, guint32, i);
      const gdouble var_int    = g_array_index (dcov->var_int, gdouble, dset_id);
      const gdouble z_cmb      = ncm_vector_get (snia_cov->z_cmb, i);
      const gdouble sigma_z    = ncm_vector_get (snia_cov->sigma_z, i);
      const gdouble emptyfac   = _nc_snia_dist_cov_calc_empty_fac (dcov, z_cmb);
      const gdouble var_z_tot  = (var_pecz + sigma_z * sigma_z) * emptyfac * emptyfac;
      const gdouble var_lens_z = var_lens * z_cmb * z_cmb;
      const gdouble var_tot    = var_z_tot + var_int + var_lens_z;

      ncm_matrix_addto (cov, i, i, var_tot);
    }
  }

  dcov->cov_cpu          = cov;
  dcov->alpha_cpu        = alpha; 
  dcov->beta_cpu         = beta; 
  dcov->lnsigma_pecz_cpu = lnsigma_pecz;
  dcov->lnsigma_lens_cpu = lnsigma_lens;

  return TRUE;
}

/**
 * nc_snia_dist_cov_mean:
 * @dcov: a #NcSNIADistCov
 * @cosmo: a #NcHICosmo
 * @snia_cov: a #NcDataSNIACov
 * @y: a #NcmVector
 *
 * FIXME
 *
 */
void
nc_snia_dist_cov_mean (NcSNIADistCov *dcov, NcHICosmo *cosmo, NcDataSNIACov *snia_cov, NcmVector *y)
{
  NcmModel *model       = NCM_MODEL (dcov);
  const gdouble mag_cut = nc_data_snia_cov_get_mag_cut (snia_cov);
  const guint mu_len    = ncm_model_vparam_len (model, NC_SNIA_DIST_COV_MU);

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
        const gdouble mag_th = mu - alpha * (width - 1.0) + beta * colour + ((thirdpar < mag_cut) ? Mcal1 : Mcal2);
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
    guint i;

    g_assert (NCM_DATA (snia_cov)->init);

    for (i = 0; i < snia_cov->mu_len; i++)
    {
      const gdouble z_cmb    = ncm_vector_get (snia_cov->z_cmb, i);
      const gdouble z_he     = ncm_vector_get (snia_cov->z_he, i);
      const gdouble width    = ncm_vector_get (snia_cov->width, i);
      const gdouble colour   = ncm_vector_get (snia_cov->colour, i);
      const gdouble thirdpar = ncm_vector_get (snia_cov->thirdpar, i);
      const gdouble dmu      = nc_distance_dmodulus_hef (dcov->dist, cosmo, z_he, z_cmb);
      const gdouble mag_th   = dmu - alpha * (width - 1.0) + beta * colour + ((thirdpar < mag_cut) ? Mcal1 : Mcal2);
      const gdouble y_i      = mag_th;

      ncm_vector_set (y, i, y_i);
    }
  }
}

/**
 * nc_snia_dist_cov_mag:
 * @dcov: a #NcSNIADistCov
 * @cosmo: a #NcHICosmo
 * @snia_cov: a #NcDataSNIACov
 * @i: the distance index
 * @width_th: the true width
 * @colour_th: the true colour
 *
 * Computes the apparent magniture from model, width and colour.
 *
 * Returns: the apparent magniture.
 */
gdouble
nc_snia_dist_cov_mag (NcSNIADistCov *dcov, NcHICosmo *cosmo, NcDataSNIACov *snia_cov, guint i, gdouble width_th, gdouble colour_th)
{
  const gdouble alpha   = ALPHA;
  const gdouble beta    = BETA;
  const gdouble DH      = nc_distance_hubble (dcov->dist, cosmo);
  const gdouble Mcal1   = ABSMAG1 + 5.0 * log10 (DH);
  const gdouble Mcal2   = ABSMAG2 + 5.0 * log10 (DH);
  const gdouble mag_cut = nc_data_snia_cov_get_mag_cut (snia_cov);

  g_assert (NCM_DATA (snia_cov)->init);

  {
    const gdouble z_cmb    = ncm_vector_get (snia_cov->z_cmb, i);
    const gdouble z_he     = ncm_vector_get (snia_cov->z_he, i);
    const gdouble thirdpar = ncm_vector_get (snia_cov->thirdpar, i);
    const gdouble dmu      = nc_distance_dmodulus_hef (dcov->dist, cosmo, z_he, z_cmb);
    const gdouble mag_th   = dmu - alpha * (width_th - 1.0) + beta * colour_th + ((thirdpar < mag_cut) ? Mcal1 : Mcal2);

    return mag_th;
  }
}

/**
 * nc_snia_dist_cov_mag_to_width_colour:
 * @dcov: a #NcSNIADistCov
 * @cosmo: a #NcHICosmo
 * @snia_cov: a #NcDataSNIACov
 * @obs: a #NcmVector
 * @X: a #NcmMatrix
 * @colmajor: whether to fill the matrices in a col-major format
 *
 * Computes effective observed vector @obs, the first @snia_cov->mu_len params
 * are set to the width colour combination using the values of the distance
 * modulus from the model @cosmo and the SNIa model @dcov, i.e.,
 * $-\alpha{}w_i+\beta{}c_i = m_{\mathrm{B},i} - \mu_{\mathrm{th},i}-\alpha-\mathcal{M}_i$.
 * The next 2 * @snia_cov->mu_len are the observed widths and then the observed colours.
 *
 * The vector @obs must be of size 3 * @snia_cov->mu_len.
 *
 */
void
nc_snia_dist_cov_mag_to_width_colour (NcSNIADistCov *dcov, NcHICosmo *cosmo, NcDataSNIACov *snia_cov, NcmVector *obs, NcmMatrix *X, gboolean colmajor)
{
  const gdouble alpha   = ALPHA;
  const gdouble beta    = BETA;
  const gdouble DH      = nc_distance_hubble (dcov->dist, cosmo);
  const gdouble Mcal1   = ABSMAG1 + 5.0 * log10 (DH);
  const gdouble Mcal2   = ABSMAG2 + 5.0 * log10 (DH);
  const gdouble mag_cut = nc_data_snia_cov_get_mag_cut (snia_cov);
  const guint mu_len    = snia_cov->mu_len;
  guint i;

  g_assert (NCM_DATA (snia_cov)->init);

  ncm_matrix_set_zero (X);

  if (colmajor)
  {
    for (i = 0; i < mu_len; i++)
    {
      const gdouble z_cmb    = ncm_vector_get (snia_cov->z_cmb, i);
      const gdouble z_he     = ncm_vector_get (snia_cov->z_he, i);
      const gdouble thirdpar = ncm_vector_get (snia_cov->thirdpar, i);
      const gdouble dmu      = nc_distance_dmodulus_hef (dcov->dist, cosmo, z_he, z_cmb);
      const gdouble M        = ((thirdpar < mag_cut) ? Mcal1 : Mcal2);
      const gdouble m_obs_i  = ncm_vector_get (snia_cov->mag, i);
      const gdouble w_obs_i  = ncm_vector_get (snia_cov->width, i);
      const gdouble c_obs_i  = ncm_vector_get (snia_cov->colour, i);

      ncm_vector_set (obs, i + 0 * mu_len, m_obs_i - dmu - alpha - M);
      ncm_vector_set (obs, i + 1 * mu_len, w_obs_i);
      ncm_vector_set (obs, i + 2 * mu_len, c_obs_i);

      ncm_matrix_set_colmajor (X, i + 0 * mu_len, i + 0 * mu_len, -alpha);
      ncm_matrix_set_colmajor (X, i + 0 * mu_len, i + 1 * mu_len, beta);
      ncm_matrix_set_colmajor (X, i + 1 * mu_len, i + 0 * mu_len, 1.0);
      ncm_matrix_set_colmajor (X, i + 2 * mu_len, i + 1 * mu_len, 1.0);
    }
  }
  else
  {
    for (i = 0; i < mu_len; i++)
    {
      const gdouble z_cmb    = ncm_vector_get (snia_cov->z_cmb, i);
      const gdouble z_he     = ncm_vector_get (snia_cov->z_he, i);
      const gdouble thirdpar = ncm_vector_get (snia_cov->thirdpar, i);
      const gdouble dmu      = nc_distance_dmodulus_hef (dcov->dist, cosmo, z_he, z_cmb);
      const gdouble M        = ((thirdpar < mag_cut) ? Mcal1 : Mcal2);
      const gdouble m_obs_i  = ncm_vector_get (snia_cov->mag, i);
      const gdouble w_obs_i  = ncm_vector_get (snia_cov->width, i);
      const gdouble c_obs_i  = ncm_vector_get (snia_cov->colour, i);

      ncm_vector_set (obs, i + 0 * mu_len, m_obs_i - dmu - alpha - M);
      ncm_vector_set (obs, i + 1 * mu_len, w_obs_i);
      ncm_vector_set (obs, i + 2 * mu_len, c_obs_i);

      ncm_matrix_set (X, i + 0 * mu_len, i + 0 * mu_len, -alpha);
      ncm_matrix_set (X, i + 0 * mu_len, i + 1 * mu_len, beta);
      ncm_matrix_set (X, i + 1 * mu_len, i + 0 * mu_len, 1.0);
      ncm_matrix_set (X, i + 2 * mu_len, i + 1 * mu_len, 1.0);
    }
  }
}

/**
 * nc_snia_dist_cov_extra_var:
 * @dcov: a #NcSNIADistCov
 * @snia_cov: a #NcDataSNIACov
 * @i: the distance index
 *
 * Computes the total variance of the @i-th distance, not related to the
 * magnitute, width or colour errors.
 *
 * Returns: the variance
 */
gdouble
nc_snia_dist_cov_extra_var (NcSNIADistCov *dcov, NcDataSNIACov *snia_cov, guint i)
{
  NcmModel *model = NCM_MODEL (dcov);
  g_assert (NCM_DATA (snia_cov)->init);

  {
    const guint dset_id       = g_array_index (snia_cov->dataset, guint32, i);
    const gdouble lnsigma_int = ncm_model_orig_vparam_get (model, NC_SNIA_DIST_COV_LNSIGMA_INT, dset_id);
    const gdouble var_int     = exp (2.0 * lnsigma_int);
    const gdouble z_cmb       = ncm_vector_get (snia_cov->z_cmb, i);
    const gdouble sigma_z     = ncm_vector_get (snia_cov->sigma_z, i);
    const gdouble emptyfac    = _nc_snia_dist_cov_calc_empty_fac (dcov, z_cmb);
    const gdouble var_pecz    = exp (2.0 * LNSIGMA_PECZ);
    const gdouble var_lens    = exp (2.0 * LNSIGMA_LENS);
    const gdouble var_lens_z  = var_lens * z_cmb * z_cmb;
    const gdouble var_z_tot   = (var_pecz + sigma_z * sigma_z) * emptyfac * emptyfac;

    const gdouble var_tot = var_z_tot + var_int + var_lens_z;

    return var_tot;
  }
}

/**
 * nc_snia_dist_cov_alpha_beta:
 * @dcov: a #NcSNIADistCov
 * @alpha: (out caller-allocates): value of alpha
 * @beta: (out caller-allocates): value of beta
 *
 * FIXME
 *
 */
void
nc_snia_dist_cov_alpha_beta (NcSNIADistCov *dcov, gdouble *alpha, gdouble *beta)
{
  *alpha = ALPHA;
  *beta  = BETA;
}

/**
 * nc_snia_dist_cov_set_default_params_by_id:
 * @dcov: a #NcSNIADistCov
 * @snia_id: a #NcDataSNIAId
 *
 * Sets the default parameters appropriated to the sample defined by @snia_id.
 *
 */
void 
nc_snia_dist_cov_set_default_params_by_id (NcSNIADistCov *dcov, NcDataSNIAId snia_id)
{
  const gdouble lnsigma0 = -15.0 * M_LN10;
  
  switch (snia_id)
  {
    case NC_DATA_SNIA_COV_SNLS3_SYS_STAT:
      g_assert_not_reached ();
      break;
    case NC_DATA_SNIA_COV_SNLS3_STAT_ONLY:
      g_assert_not_reached ();
      break;
    case NC_DATA_SNIA_COV_JLA_SNLS3_SDSS_SYS_STAT:
    {
      gdouble lnsigma_int_data[4] = {lnsigma0, lnsigma0, lnsigma0, lnsigma0};
      NcmVector *lnsigma_int      = ncm_vector_new_data_static (lnsigma_int_data, 4, 1);
      g_object_set (dcov, 
                    "alpha",        0.141,
                    "beta",         2.60,
                    "M1",           -19.0497380934588,
                    "M2",           -19.1196618476607,
                    "lnsigma_pecz", lnsigma0,
                    "lnsigma_lens", lnsigma0,
                    "lnsigma_int",  lnsigma_int,
                    "alpha-fit",    TRUE,
                    "beta-fit",     TRUE,
                    "M1-fit",       TRUE,
                    "M2-fit",       TRUE,
                    NULL
                    );
      ncm_vector_free (lnsigma_int);
      nc_snia_dist_cov_set_empty_fac (dcov, FALSE);
      break;
    }
    case NC_DATA_SNIA_COV_JLA_SNLS3_SDSS_SYS_STAT_CMPL:
    {
      gdouble lnsigma_int_data[4] = {log (0.08), log (0.108), log (0.134), log (0.1)};
      NcmVector *lnsigma_int      = ncm_vector_new_data_static (lnsigma_int_data, 4, 1);
      g_object_set (dcov, 
                    "alpha",        0.141,
                    "beta",         2.60,
                    "M1",           -19.0497090404172,
                    "M2",           -19.1196499780362,
                    "lnsigma_pecz", log (150000.0 / ncm_c_c ()),
                    "lnsigma_lens", log (0.055),
                    "lnsigma_int",  lnsigma_int,
                    "alpha-fit",    TRUE,
                    "beta-fit",     TRUE,
                    "M1-fit",       TRUE,
                    "M2-fit",       TRUE,
                    NULL
                    );
      ncm_vector_free (lnsigma_int);
      nc_snia_dist_cov_set_empty_fac (dcov, FALSE);
      break;
    }
    case NC_DATA_SNIA_COV_PANTHEON:
    {
      gdouble lnsigma_int_data[1] = {lnsigma0};
      NcmVector *lnsigma_int      = ncm_vector_new_data_static (lnsigma_int_data, 1, 1);
      g_object_set (dcov, 
                    "alpha",        0.0,
                    "beta",         0.0,
                    "M1",           -19.0497380934588,
                    "M2",           -19.1196618476607,
                    "lnsigma_pecz", lnsigma0,
                    "lnsigma_lens", lnsigma0,
                    "lnsigma_int",  lnsigma_int,
                    "alpha-fit",    FALSE,
                    "beta-fit",     FALSE,
                    "M1-fit",       TRUE,
                    "M2-fit",       FALSE,
                    NULL
                    );
      ncm_vector_free (lnsigma_int);
      nc_snia_dist_cov_set_empty_fac (dcov, FALSE);
      break;
    }
    default:
      g_assert_not_reached ();
      break;      
  }
}
