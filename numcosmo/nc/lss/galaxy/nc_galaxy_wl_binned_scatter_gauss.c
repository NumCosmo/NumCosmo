/***************************************************************************
 *            nc_galaxy_wl_binned_scatter_gauss.c
 *
 *  Thu November 26 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * nc_galaxy_wl_binned_scatter_gauss.c
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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
 * NcGalaxyWLBinnedScatterGauss:
 *
 * Gaussian scatter of a binned weak-lensing profile.
 *
 * The measured profile $d$ is Gaussian around the model $f$ with the covariance $C$ of the
 * profile (#NcGalaxyWLObsBinned:cov),
 * $$-2\ln L = (d - f)^T C^{-1} (d - f),$$
 * as in Umetsu et al. 2020 (Eq. 21). The normalization
 * $\ln\det(2\pi C)$ does not depend on the parameters and is left out unless
 * #NcGalaxyWLBinnedScatterGauss:include-lnnorm is %TRUE. The model has no parameters.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_wl_binned_scatter_gauss.h"
#include "ncm/core/ncm_c.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <math.h>
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcGalaxyWLBinnedScatterGauss
{
  /*< private >*/
  NcGalaxyWLBinnedScatter parent_instance;
  gboolean include_lnnorm;
  NcmMatrix *chol;
  NcmVector *diff;
};

enum
{
  PROP_0,
  PROP_INCLUDE_LNNORM,
  PROP_LEN,
};

G_DEFINE_TYPE (NcGalaxyWLBinnedScatterGauss, nc_galaxy_wl_binned_scatter_gauss, NC_TYPE_GALAXY_WL_BINNED_SCATTER)

static void
nc_galaxy_wl_binned_scatter_gauss_init (NcGalaxyWLBinnedScatterGauss *sg)
{
  sg->include_lnnorm = FALSE;
  sg->chol           = NULL;
  sg->diff           = NULL;
}

static void
_nc_galaxy_wl_binned_scatter_gauss_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLBinnedScatterGauss *sg = NC_GALAXY_WL_BINNED_SCATTER_GAUSS (object);

  switch (prop_id)
  {
    case PROP_INCLUDE_LNNORM:
      nc_galaxy_wl_binned_scatter_gauss_set_include_lnnorm (sg, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_binned_scatter_gauss_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcGalaxyWLBinnedScatterGauss *sg = NC_GALAXY_WL_BINNED_SCATTER_GAUSS (object);

  switch (prop_id)
  {
    case PROP_INCLUDE_LNNORM:
      g_value_set_boolean (value, nc_galaxy_wl_binned_scatter_gauss_get_include_lnnorm (sg));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_nc_galaxy_wl_binned_scatter_gauss_dispose (GObject *object)
{
  NcGalaxyWLBinnedScatterGauss *sg = NC_GALAXY_WL_BINNED_SCATTER_GAUSS (object);

  ncm_matrix_clear (&sg->chol);
  ncm_vector_clear (&sg->diff);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_galaxy_wl_binned_scatter_gauss_parent_class)->dispose (object);
}

/*
 * Cholesky decomposition of the covariance of binned into sg->chol (lower triangle in the
 * row-major convention of NcmMatrix, 'U' for LAPACK). The covariance has a handful of bins,
 * so it is decomposed at every call: no cache to keep in step with the profile.
 */
static void
_nc_galaxy_wl_binned_scatter_gauss_decomp (NcGalaxyWLBinnedScatterGauss *sg, NcGalaxyWLObsBinned *binned)
{
  NcmMatrix *cov     = nc_galaxy_wl_obs_binned_peek_cov (binned);
  const guint n_bins = ncm_matrix_nrows (cov);
  gint ret;

  if ((sg->chol == NULL) || (ncm_matrix_nrows (sg->chol) != n_bins))
  {
    ncm_matrix_clear (&sg->chol);
    ncm_vector_clear (&sg->diff);
    sg->chol = ncm_matrix_new (n_bins, n_bins);
    sg->diff = ncm_vector_new (n_bins);
  }

  ncm_matrix_memcpy (sg->chol, cov);
  ret = ncm_matrix_cholesky_decomp (sg->chol, 'U');

  if (ret != 0)
    g_error ("nc_galaxy_wl_binned_scatter_gauss: the covariance of the profile is not positive definite (%d).", ret);
}

static gdouble
_nc_galaxy_wl_binned_scatter_gauss_m2lnL (NcGalaxyWLBinnedScatter *scatter, NcGalaxyWLObsBinned *binned, NcmVector *model)
{
  NcGalaxyWLBinnedScatterGauss *sg = NC_GALAXY_WL_BINNED_SCATTER_GAUSS (scatter);
  const guint n_bins               = nc_galaxy_wl_obs_binned_len (binned);
  gdouble m2lnL                    = 0.0;
  guint b;

  _nc_galaxy_wl_binned_scatter_gauss_decomp (sg, binned);

  for (b = 0; b < n_bins; b++)
  {
    const gdouble d = nc_galaxy_wl_obs_binned_get (binned, NC_GALAXY_WL_OBS_BINNED_COL_VALUE, b);

    ncm_vector_set (sg->diff, b, d - ncm_vector_get (model, b));
  }

  {
    NcmVector *sol = ncm_vector_dup (sg->diff);
    gint ret       = ncm_matrix_cholesky_solve2 (sg->chol, sol, 'U');

    if (ret != 0)
      g_error ("nc_galaxy_wl_binned_scatter_gauss: failed to solve the covariance system (%d).", ret);

    m2lnL = ncm_vector_dot (sg->diff, sol);
    ncm_vector_free (sol);
  }

  if (sg->include_lnnorm)
    m2lnL += n_bins * ncm_c_ln2pi () + ncm_matrix_cholesky_lndet (sg->chol);

  return m2lnL;
}

static void
_nc_galaxy_wl_binned_scatter_gauss_gen (NcGalaxyWLBinnedScatter *scatter, NcGalaxyWLObsBinned *binned, NcmVector *model, NcmRNG *rng, NcmVector *data)
{
  NcGalaxyWLBinnedScatterGauss *sg = NC_GALAXY_WL_BINNED_SCATTER_GAUSS (scatter);
  const guint n_bins               = nc_galaxy_wl_obs_binned_len (binned);
  guint i, j;

  _nc_galaxy_wl_binned_scatter_gauss_decomp (sg, binned);

  ncm_rng_lock (rng);

  for (i = 0; i < n_bins; i++)
    ncm_vector_set (sg->diff, i, ncm_rng_gaussian_gen (rng, 0.0, 1.0));

  ncm_rng_unlock (rng);

  /* data = model + L z, with C = L L^T and L = chol^T (chol holds U = L^T in row-major) */
  for (i = 0; i < n_bins; i++)
  {
    gdouble v = ncm_vector_get (model, i);

    for (j = 0; j <= i; j++)
      v += ncm_matrix_get (sg->chol, j, i) * ncm_vector_get (sg->diff, j);

    ncm_vector_set (data, i, v);
  }
}

static void
nc_galaxy_wl_binned_scatter_gauss_class_init (NcGalaxyWLBinnedScatterGaussClass *klass)
{
  GObjectClass *object_class                  = G_OBJECT_CLASS (klass);
  NcmModelClass *model_class                  = NCM_MODEL_CLASS (klass);
  NcGalaxyWLBinnedScatterClass *scatter_class = NC_GALAXY_WL_BINNED_SCATTER_CLASS (klass);

  model_class->set_property = &_nc_galaxy_wl_binned_scatter_gauss_set_property;
  model_class->get_property = &_nc_galaxy_wl_binned_scatter_gauss_get_property;
  object_class->dispose     = &_nc_galaxy_wl_binned_scatter_gauss_dispose;

  ncm_model_class_set_name_nick (model_class, "Gaussian binned weak-lensing profile scatter", "WLBinnedScatterGauss");
  ncm_model_class_add_params (model_class, 0, 0, PROP_LEN);

  /**
   * NcGalaxyWLBinnedScatterGauss:include-lnnorm:
   *
   * Whether $-2\ln L$ includes the normalization $\ln\det(2\pi C)$.
   */
  g_object_class_install_property (object_class,
                                   PROP_INCLUDE_LNNORM,
                                   g_param_spec_boolean ("include-lnnorm",
                                                         NULL,
                                                         "Whether to include the normalization",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  ncm_model_class_check_params_info (model_class);

  scatter_class->m2lnL = &_nc_galaxy_wl_binned_scatter_gauss_m2lnL;
  scatter_class->gen   = &_nc_galaxy_wl_binned_scatter_gauss_gen;
}

/**
 * nc_galaxy_wl_binned_scatter_gauss_new:
 *
 * Creates a new #NcGalaxyWLBinnedScatterGauss.
 *
 * Returns: (transfer full): a new #NcGalaxyWLBinnedScatterGauss.
 */
NcGalaxyWLBinnedScatterGauss *
nc_galaxy_wl_binned_scatter_gauss_new (void)
{
  return g_object_new (NC_TYPE_GALAXY_WL_BINNED_SCATTER_GAUSS, NULL);
}

/**
 * nc_galaxy_wl_binned_scatter_gauss_ref:
 * @sg: a #NcGalaxyWLBinnedScatterGauss
 *
 * Increases the reference count of @sg by one.
 *
 * Returns: (transfer full): @sg.
 */
NcGalaxyWLBinnedScatterGauss *
nc_galaxy_wl_binned_scatter_gauss_ref (NcGalaxyWLBinnedScatterGauss *sg)
{
  return g_object_ref (sg);
}

/**
 * nc_galaxy_wl_binned_scatter_gauss_free:
 * @sg: a #NcGalaxyWLBinnedScatterGauss
 *
 * Decreases the reference count of @sg by one.
 */
void
nc_galaxy_wl_binned_scatter_gauss_free (NcGalaxyWLBinnedScatterGauss *sg)
{
  g_object_unref (sg);
}

/**
 * nc_galaxy_wl_binned_scatter_gauss_clear:
 * @sg: a #NcGalaxyWLBinnedScatterGauss
 *
 * Decreases the reference count of *@sg by one and sets *@sg to %NULL.
 */
void
nc_galaxy_wl_binned_scatter_gauss_clear (NcGalaxyWLBinnedScatterGauss **sg)
{
  g_clear_object (sg);
}

/**
 * nc_galaxy_wl_binned_scatter_gauss_set_include_lnnorm:
 * @sg: a #NcGalaxyWLBinnedScatterGauss
 * @include_lnnorm: whether to include the normalization
 *
 * Sets #NcGalaxyWLBinnedScatterGauss:include-lnnorm.
 */
void
nc_galaxy_wl_binned_scatter_gauss_set_include_lnnorm (NcGalaxyWLBinnedScatterGauss *sg, gboolean include_lnnorm)
{
  sg->include_lnnorm = include_lnnorm;
}

/**
 * nc_galaxy_wl_binned_scatter_gauss_get_include_lnnorm:
 * @sg: a #NcGalaxyWLBinnedScatterGauss
 *
 * Returns: #NcGalaxyWLBinnedScatterGauss:include-lnnorm.
 */
gboolean
nc_galaxy_wl_binned_scatter_gauss_get_include_lnnorm (NcGalaxyWLBinnedScatterGauss *sg)
{
  return sg->include_lnnorm;
}

