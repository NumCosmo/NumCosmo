/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_data_gauss_cov_mvnd.c
 *
 *  Sun February 04 15:07:28 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_data_gauss_cov_mvnd.c
 * Copyright (C) 2018 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_data_gauss_cov_mvnd
 * @title: NcmDataGaussCov
 * @short_description: Multivariate Normal Distribution -- covariance provided.
 *
 * Multivariate Normal distribution which uses the covariance matrix as input.
 * It should be used with its companion object #NcmModelMVND.
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_data_gauss_cov_mvnd.h"
#include "math/ncm_model_mvnd.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmDataGaussCovMVND
{
  NcmDataGaussCov parent_instance;
};

G_DEFINE_TYPE (NcmDataGaussCovMVND, ncm_data_gauss_cov_mvnd, NCM_TYPE_DATA_GAUSS_COV);

static void
ncm_data_gauss_cov_mvnd_init (NcmDataGaussCovMVND *gauss_mvnd)
{
}

static void
ncm_data_gauss_cov_mvnd_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_data_gauss_cov_mvnd_parent_class)->finalize (object);
}

static void _ncm_data_gauss_cov_mvnd_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp);

static void
ncm_data_gauss_cov_mvnd_class_init (NcmDataGaussCovMVNDClass *klass)
{
  GObjectClass *object_class        = G_OBJECT_CLASS (klass);
  NcmDataGaussCovClass *gauss_class = NCM_DATA_GAUSS_COV_CLASS (klass);

  object_class->finalize = ncm_data_gauss_cov_mvnd_finalize;

  /*data_class->prepare    = &_ncm_data_gauss_cov_test_prepare;*/
  gauss_class->mean_func = &_ncm_data_gauss_cov_mvnd_mean_func;
  gauss_class->cov_func  = NULL;
}

static void
_ncm_data_gauss_cov_mvnd_mean_func (NcmDataGaussCov *gauss, NcmMSet *mset, NcmVector *vp)
{
  NcmModelMVND *model_mvnd = NCM_MODEL_MVND (ncm_mset_peek (mset, ncm_model_mvnd_id ()));

  ncm_model_mvnd_mean (model_mvnd, vp);
}

/**
 * ncm_data_gauss_cov_mvnd_new:
 * @dim: dimension of the MVND
 *
 * Creates a new @dim-dimensional MVND.
 *
 * Returns: the newly created object.
 */
NcmDataGaussCovMVND *
ncm_data_gauss_cov_mvnd_new (const guint dim)
{
  NcmDataGaussCovMVND *gauss_mvnd = g_object_new (NCM_TYPE_DATA_GAUSS_COV_MVND,
                                                  "n-points", dim,
                                                  "use-norma", TRUE,
                                                  NULL);

  return gauss_mvnd;
}

/**
 * ncm_data_gauss_cov_mvnd_new_full:
 * @dim: dimension of the MVND
 * @sigma_min: minimum value of $\sigma_i$
 * @sigma_max: maximum value of $\sigma_i$
 * @cor_level: correlation level
 * @mean_min: minimum mean $\mu_i$
 * @mean_max: maximum mean $\mu_i$
 * @rng: a #NcmRNG
 *
 * Creates a new @dim-dimensional MVND and generate using @rng a mean
 * and correlation matrix using the parameters above.
 *
 * Returns: the newly created object.
 */
NcmDataGaussCovMVND *
ncm_data_gauss_cov_mvnd_new_full (const guint dim, const gdouble sigma_min, const gdouble sigma_max, const gdouble cor_level, const gdouble mean_min, const gdouble mean_max, NcmRNG *rng)
{
  NcmDataGaussCovMVND *gauss_mvnd = ncm_data_gauss_cov_mvnd_new (dim);

  ncm_data_gauss_cov_mvnd_gen_cov_mean (gauss_mvnd, sigma_min, sigma_max, cor_level, mean_min, mean_max, rng);

  return gauss_mvnd;
}

/**
 * ncm_data_gauss_cov_mvnd_ref:
 * @data_mvnd: a #NcmDataGaussCovMVND
 *
 * Increases the reference count of @data_mvnd by one.
 *
 * Returns: (transfer full): @data_mvnd
 */
NcmDataGaussCovMVND *
ncm_data_gauss_cov_mvnd_ref (NcmDataGaussCovMVND *data_mvnd)
{
  return g_object_ref (data_mvnd);
}

/**
 * ncm_data_gauss_cov_mvnd_free:
 * @data_mvnd: a #NcmDataGaussCovMVND
 *
 * Decreases the reference count of @data_mvnd by one.
 *
 */
void
ncm_data_gauss_cov_mvnd_free (NcmDataGaussCovMVND *data_mvnd)
{
  g_object_unref (data_mvnd);
}

/**
 * ncm_data_gauss_cov_mvnd_clear:
 * @data_mvnd: a #NcmDataGaussCovMVND
 *
 * If @data_mvnd is different from NULL, decreases the reference count of
 * @data_mvnd by one and sets @data_mvnd to NULL.
 *
 */
void
ncm_data_gauss_cov_mvnd_clear (NcmDataGaussCovMVND **data_mvnd)
{
  g_clear_object (data_mvnd);
}

/**
 * ncm_data_gauss_cov_mvnd_gen_cov_mean:
 * @data_mvnd: a #NcmDataGaussCovMVND
 * @sigma_min: minimum value of $\sigma_i$
 * @sigma_max: maximum value of $\sigma_i$
 * @cor_level: correlation level
 * @mean_min: minimum mean $\mu_i$
 * @mean_max: maximum mean $\mu_i$
 * @rng: a #NcmRNG
 *
 * Generates using @rng the mean and correlation matrix using
 * the parameters above.
 *
 */
void
ncm_data_gauss_cov_mvnd_gen_cov_mean (NcmDataGaussCovMVND *data_mvnd, const gdouble sigma_min, const gdouble sigma_max, const gdouble cor_level, const gdouble mean_min, const gdouble mean_max, NcmRNG *rng)
{
  NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);
  NcmVector *y   = ncm_data_gauss_cov_peek_mean (gcov);
  NcmMatrix *cov = ncm_data_gauss_cov_peek_cov (gcov);
  const guint np = ncm_data_gauss_cov_get_size (gcov);
  gint i;

  g_assert_cmpfloat (mean_min, <=, mean_max);
  g_assert_cmpint (np, >, 0);
  g_assert (y != NULL);
  g_assert (cov != NULL);

  ncm_matrix_fill_rand_cov (cov, sigma_min, sigma_max, cor_level, rng);

  if (mean_min == mean_max)
  {
    ncm_vector_set_all (y, mean_min);
  }
  else
  {
    for (i = 0; i < np; i++)
    {
      const gdouble mean = ncm_rng_uniform_gen (rng, mean_min, mean_max);

      ncm_vector_set (y, i, mean);
    }
  }

  ncm_data_set_init (NCM_DATA (gcov), TRUE);
}

/**
 * ncm_data_gauss_cov_mvnd_set_cov_mean:
 * @data_mvnd: a #NcmDataGaussCovMVND
 * @mean: a #NcmVector
 * @cov: a #NcmMatrix
 *
 * Sets the mean and covariance of @data_mvnd.
 *
 */
void
ncm_data_gauss_cov_mvnd_set_cov_mean (NcmDataGaussCovMVND *data_mvnd, NcmVector *mean, NcmMatrix *cov)
{
  NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);
  NcmVector *cy         = ncm_data_gauss_cov_peek_mean (gcov);
  NcmMatrix *ccov       = ncm_data_gauss_cov_peek_cov (gcov);

  g_assert_cmpuint (ncm_vector_len (mean), ==, ncm_matrix_nrows (cov));
  g_assert_cmpuint (ncm_vector_len (mean), ==, ncm_matrix_ncols (cov));
  g_assert (cy != NULL);
  g_assert (cov != NULL);

  ncm_matrix_memcpy (ccov, cov);
  ncm_vector_memcpy (cy, mean);

  ncm_data_set_init (NCM_DATA (gcov), TRUE);
}

/**
 * ncm_data_gauss_cov_mvnd_peek_mean:
 * @data_mvnd: a #NcmDataGaussCovMVND
 *
 * Peeks current mean vector.
 *
 * Returns: (transfer none): the current mean vector.
 */
NcmVector *
ncm_data_gauss_cov_mvnd_peek_mean (NcmDataGaussCovMVND *data_mvnd)
{
  NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (gcov);

  return y;
}

/**
 * ncm_data_gauss_cov_mvnd_gen:
 * @data_mvnd: a #NcmDataGaussCovMVND
 * @mset: a #NcmMSet
 * @obj: (allow-none): a pointer to use in @bound
 * @bound: (scope call) (allow-none): a NcmDataGaussCovMVNDBound
 * @rng: a #NcmRNG
 * @N: (out): number of realizations necessary to generate a valid one
 *
 * Generates one realization of the MVND. If @bound is not NULL,
 * generates realizations untill @bound returns TRUE.
 *
 * Returns: (transfer none): a #NcmVector (should not be modified)
 */
NcmVector *
ncm_data_gauss_cov_mvnd_gen (NcmDataGaussCovMVND *data_mvnd, NcmMSet *mset, gpointer obj, NcmDataGaussCovMVNDBound bound, NcmRNG *rng, gulong *N)
{
  NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);
  NcmData *data         = NCM_DATA (data_mvnd);
  NcmVector *y          = ncm_data_gauss_cov_peek_mean (gcov);
  gulong maxiter        = 100000000;

  N[0] = 0;

  if (bound != NULL)
  {
    do {
      ncm_data_resample (data, mset, rng);
      N[0]++;

      if (N[0] > maxiter)
      {
        g_error ("ncm_data_gauss_cov_mvnd_gen: too many interations, cannot find a valid realization!");
        break;
      }
    } while (!bound (obj, y));
  }
  else
  {
    ncm_data_resample (data, mset, rng);
  }

  return y;
}

/**
 * ncm_data_gauss_cov_mvnd_est_ratio:
 * @data_mvnd: a #NcmDataGaussCovMVND
 * @mset: a #NcmMSet
 * @obj: a pointer to use in @bound
 * @bound: (scope call): a NcmDataGaussCovMVNDBound
 * @N: (inout): total number of realizations
 * @Nin: (inout): number of realizations accepted
 * @reltol: relative tolerance
 * @rng: a #NcmRNG
 *
 * Estimate the ratio between accepted realizations and total number
 * of realizations. The variable @reltol controls the relative tolerance
 * on the ratio estimate. The variables @N and @Nin can be used to inform
 * previous number of realizations.
 *
 * Returns: estimated ratio
 */
gdouble
ncm_data_gauss_cov_mvnd_est_ratio (NcmDataGaussCovMVND *data_mvnd, NcmMSet *mset, gpointer obj, NcmDataGaussCovMVNDBound bound, gulong *N, gulong *Nin, const gdouble reltol, NcmRNG *rng)
{
  gulong NN   = (N   != NULL) ? *N   : 0;
  gulong NNin = (Nin != NULL) ? *Nin : 0;
  gdouble ratio;
  gulong maxiter = 100000000;
  glong miniter  = MIN (10000, (glong) (1.0 / reltol));

  while (TRUE)
  {
    gdouble err_rel;
    gulong Ni;

    ncm_data_gauss_cov_mvnd_gen (data_mvnd, mset, obj, bound, rng, &Ni);

    NN += Ni;
    NNin++;

    ratio   = NNin * 1.0 / (1.0 * NN);
    err_rel = sqrt ((1.0 - ratio) / (NN * ratio));

    if ((miniter < NN) && (err_rel < reltol))
      break;

    if (NNin > maxiter)
    {
      g_warning ("ncm_data_gauss_cov_mvnd_est_ratio: too many iterations, result is not trustworthy!");
      break;
    }
  }

  if (N != NULL)
    *N = NN;

  if (Nin != NULL)
    *Nin = NNin;

  return ratio;
}

/**
 * ncm_data_gauss_cov_mvnd_log_info:
 * @data_mvnd: a #NcmDataGaussCovMVND
 *
 * Logs mean and covariance matrix.
 *
 */
void
ncm_data_gauss_cov_mvnd_log_info (NcmDataGaussCovMVND *data_mvnd)
{
  NcmDataGaussCov *gcov = NCM_DATA_GAUSS_COV (data_mvnd);
  NcmVector *y         = ncm_data_gauss_cov_peek_mean (gcov);
  NcmMatrix *cov       = ncm_data_gauss_cov_peek_cov (gcov);

  g_assert (y != NULL);
  g_assert (cov != NULL);

  ncm_vector_log_vals (y,   "# NcmDataGaussCovMVND data mean: ", "% 12.5g", TRUE);
  ncm_matrix_log_vals (cov, "# NcmDataGaussCovMVND data cov: ", "% 12.5g");
}

