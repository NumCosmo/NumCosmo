/***************************************************************************
 *            ncm_stats_dist.c
 *
 *  Wed November 07 16:02:36 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_stats_dist.c
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
 * NcmStatsDist:
 *
 * Base class for N-dimensional probability distributions reconstructed from samples.
 *
 * Represents a distribution as a mixture of radial-basis kernels placed on a set of
 * sample points, evaluates it and draws from it. The theory, and what the shrinkage
 * and kernel-selection options do, are on the <a
 * href="../../theory/stats_dist.html">Kernel Mixture Densities</a> page.
 *
 * This is an abstract class. #NcmStatsDistKDE gives all kernels a common bandwidth,
 * #NcmStatsDistVKDE lets it vary between sample points; both build the interpolation
 * matrix and the covariance decompositions. The kernel profile comes from a
 * #NcmStatsDistKernel, either #NcmStatsDistKernelGauss or #NcmStatsDistKernelST.
 *
 * Add the sample points with ncm_stats_dist_add_obs(), then call either
 * ncm_stats_dist_prepare(), which weights every kernel equally, or
 * ncm_stats_dist_prepare_interp() with the values $-2\ln g(x_i)$, which solves a
 * non-negative least-squares problem for the weights. Both leave the object ready for
 * ncm_stats_dist_eval() and ncm_stats_dist_sample().
 *
 * Nothing else has to be set. #NcmStatsDist:over-smooth sets the bandwidth,
 * #NcmStatsDist:CV-type the cross-validation that can fit it and
 * #NcmStatsDist:split-frac the fraction held out; #NcmStatsDist:center-shrink
 * contracts the mixture so its covariance matches the sample, and
 * #NcmStatsDist:auto-kernel fits the kernel tail along with the bandwidth, which
 * requires #NCM_STATS_DIST_CV_SPLIT_NOFIT and replaces the kernel the object was built
 * with.
 *
 * #NcmStatsDist:defensive-frac mixes a wide Student-t component into the proposal so
 * that its density has no holes; it is off by default.
 *
 * Center shrinkage requires a kernel with a finite covariance: a #NcmStatsDistKernelST
 * with $\nu \leq 2$, the Cauchy kernel included, aborts. It also raises the lower
 * bound on a fitted $\nu$ from $1$ to $2.5$.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/stats/ncm_stats_dist.h"
#include "ncm/stats/ncm_stats_dist_kde.h"
#include "ncm/stats/ncm_stats_dist_kernel_gauss.h"
#include "ncm/stats/ncm_stats_dist_kernel_st.h"
#include <nlopt.h>
#include "ncm/core/ncm_iset.h"
#include "ncm/core/ncm_rng.h"
#include "ncm/algebra/ncm_lapack.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sort.h>
#include "external/levmar/levmar.h"
#endif /* NUMCOSMO_GIR_SCAN */

#include "ncm/stats/ncm_stats_dist_private.h"

enum
{
  PROP_0,
  PROP_KERNEL,
  PROP_SAMPLE_SIZE,
  PROP_OVER_SMOOTH,
  PROP_CV_TYPE,
  PROP_USE_THREADS,
  PROP_SPLIT_FRAC,
  PROP_PRINT_FIT,
  PROP_CENTER_SHRINK,
  PROP_AUTO_KERNEL,
  PROP_DEFENSIVE_FRAC,
  PROP_DEFENSIVE_SCALE,
  PROP_DEFENSIVE_NU,
  PROP_UNIFORM_WEIGHTS,
};

G_DEFINE_ABSTRACT_TYPE_WITH_PRIVATE (NcmStatsDist, ncm_stats_dist, G_TYPE_OBJECT)

#define NCM_NNLS_SOLVE ncm_nnls_solve

static void _ncm_stats_dist_set_href (NcmStatsDist *sd, const gdouble href);

/*
 * Nominal bandwidth, computed by the subclass from over_smooth and the kernel rule of
 * thumb. The bandwidth actually applied to the kernels is set by
 * _ncm_stats_dist_set_href() and read back with ncm_stats_dist_get_href().
 */
static gdouble
_ncm_stats_dist_calc_href (NcmStatsDist *sd)
{
  return NCM_STATS_DIST_GET_CLASS (sd)->get_href (sd);
}

static void
ncm_stats_dist_init (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->kernel               = NULL;
  self->sample_array         = g_ptr_array_new ();
  self->weights              = NULL;
  self->wcum                 = NULL;
  self->wcum_ready           = FALSE;
  self->print_fit            = FALSE;
  self->over_smooth          = 0.0;
  self->cv_type              = NCM_STATS_DIST_CV_LEN;
  self->use_threads          = FALSE;
  self->split_frac           = 0.0;
  self->min_m2lnp            = 0.0;
  self->max_m2lnp            = 0.0;
  self->href                 = 0.0;
  self->rnorm                = 0.0;
  self->center_shrink        = FALSE;
  self->auto_kernel          = FALSE;
  self->center_array         = g_ptr_array_new ();
  self->center_mean          = NULL;
  self->center_C_decomp      = NULL;
  self->center_mean_cov      = NULL;
  self->center_A             = NULL;
  self->center_Ahat          = NULL;
  self->center_Ahat_identity = TRUE;
  self->center_a             = 1.0;
  self->defensive_frac       = 0.0;
  self->defensive_scale      = 4.0;
  self->defensive_nu         = 3.0;
  self->defensive_kernel     = NULL;
  self->defensive_decomp     = NULL;
  self->defensive_lnnorm     = 0.0;
  self->n_obs                = 0;
  self->n_kernels            = 0;
  self->alloc_n_obs          = 0;
  self->alloc_n_kernels      = 0;
  self->alloc_subs           = FALSE;
  self->d                    = 0;
  self->sampling             = g_array_new (FALSE, FALSE, sizeof (guint));
  self->nnls                 = NULL;
  self->IM                   = NULL;
  self->sub_IM               = NULL;
  self->sub_x                = NULL;
  self->f                    = NULL;
  self->f1                   = NULL;
  self->cv_m2lnp             = NULL;
  self->cv_m2lnL_sample      = NULL;
  self->cv_w                 = NULL;
  self->uniform_weights      = FALSE;
  self->levmar_workz         = NULL;
  self->levmar_n             = 0;
  self->fmin                 = gsl_multimin_fminimizer_alloc (gsl_multimin_fminimizer_nmsimplex2, 1);
  self->m2lnp_sort           = g_array_new (FALSE, FALSE, sizeof (size_t));
  self->m2lnp                = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->rng                  = ncm_rng_seeded_new (NULL, 0);

  g_ptr_array_set_free_func (self->sample_array, (GDestroyNotify) ncm_vector_free);
  g_ptr_array_set_free_func (self->center_array, (GDestroyNotify) ncm_vector_free);
}

static void
_ncm_stats_dist_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmStatsDist *sd = NCM_STATS_DIST (object);

  /*g_return_if_fail (NCM_IS_STATS_DIST (object));*/

  switch (prop_id)
  {
    case PROP_KERNEL:
      ncm_stats_dist_set_kernel (sd, g_value_get_object (value));
      break;
    case PROP_SAMPLE_SIZE: /* LCOV_EXCL_BR_LINE */
      g_assert_not_reached ();
      break;
    case PROP_OVER_SMOOTH:
      ncm_stats_dist_set_over_smooth (sd, g_value_get_double (value));
      break;
    case PROP_CV_TYPE:
      ncm_stats_dist_set_cv_type (sd, g_value_get_enum (value));
      break;
    case PROP_USE_THREADS:
      ncm_stats_dist_set_use_threads (sd, g_value_get_boolean (value));
      break;
    case PROP_SPLIT_FRAC:
      ncm_stats_dist_set_split_frac (sd, g_value_get_double (value));
      break;
    case PROP_PRINT_FIT:
      ncm_stats_dist_set_print_fit (sd, g_value_get_boolean (value));
      break;
    case PROP_CENTER_SHRINK:
      ncm_stats_dist_set_center_shrink (sd, g_value_get_boolean (value));
      break;
    case PROP_AUTO_KERNEL:
      ncm_stats_dist_set_auto_kernel (sd, g_value_get_boolean (value));
      break;
    case PROP_DEFENSIVE_FRAC:
      ncm_stats_dist_set_defensive_frac (sd, g_value_get_double (value));
      break;
    case PROP_DEFENSIVE_SCALE:
      ncm_stats_dist_set_defensive_scale (sd, g_value_get_double (value));
      break;
    case PROP_DEFENSIVE_NU:
      ncm_stats_dist_set_defensive_nu (sd, g_value_get_double (value));
      break;
    case PROP_UNIFORM_WEIGHTS:
      ncm_stats_dist_set_uniform_weights (sd, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_dist_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmStatsDist *sd                 = NCM_STATS_DIST (object);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  g_return_if_fail (NCM_IS_STATS_DIST (object));

  switch (prop_id)
  {
    case PROP_KERNEL:
      g_value_set_object (value, ncm_stats_dist_peek_kernel (sd));
      break;
    case PROP_SAMPLE_SIZE:                               /* LCOV_EXCL_BR_LINE */
      g_value_set_uint (value, self->sample_array->len); /* LCOV_EXCL_LINE */
      break;
    case PROP_OVER_SMOOTH:
      g_value_set_double (value, ncm_stats_dist_get_over_smooth (sd));
      break;
    case PROP_CV_TYPE:
      g_value_set_enum (value, ncm_stats_dist_get_cv_type (sd));
      break;
    case PROP_USE_THREADS:
      g_value_set_boolean (value, ncm_stats_dist_get_use_threads (sd));
      break;
    case PROP_SPLIT_FRAC:
      g_value_set_double (value, ncm_stats_dist_get_split_frac (sd));
      break;
    case PROP_PRINT_FIT:
      g_value_set_boolean (value, ncm_stats_dist_get_print_fit (sd));
      break;
    case PROP_CENTER_SHRINK:
      g_value_set_boolean (value, ncm_stats_dist_get_center_shrink (sd));
      break;
    case PROP_AUTO_KERNEL:
      g_value_set_boolean (value, ncm_stats_dist_get_auto_kernel (sd));
      break;
    case PROP_DEFENSIVE_FRAC:
      g_value_set_double (value, ncm_stats_dist_get_defensive_frac (sd));
      break;
    case PROP_DEFENSIVE_SCALE:
      g_value_set_double (value, ncm_stats_dist_get_defensive_scale (sd));
      break;
    case PROP_DEFENSIVE_NU:
      g_value_set_double (value, ncm_stats_dist_get_defensive_nu (sd));
      break;
    case PROP_UNIFORM_WEIGHTS:
      g_value_set_boolean (value, ncm_stats_dist_get_uniform_weights (sd));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_stats_dist_dispose (GObject *object)
{
  NcmStatsDist *sd                 = NCM_STATS_DIST (object);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  ncm_stats_dist_kernel_clear (&self->kernel);

  g_clear_pointer (&self->sample_array, g_ptr_array_unref); /* LCOV_EXCL_BR_LINE */
  g_clear_pointer (&self->center_array, g_ptr_array_unref); /* LCOV_EXCL_BR_LINE */
  ncm_vector_clear (&self->center_mean);
  ncm_matrix_clear (&self->center_C_decomp);
  ncm_matrix_clear (&self->center_mean_cov);
  ncm_matrix_clear (&self->center_A);
  ncm_matrix_clear (&self->center_Ahat);
  ncm_stats_dist_kernel_clear (&self->defensive_kernel);
  ncm_matrix_clear (&self->defensive_decomp);
  ncm_vector_clear (&self->weights);
  ncm_vector_clear (&self->wcum);

  g_clear_pointer (&self->sampling, g_array_unref); /* LCOV_EXCL_BR_LINE */

  ncm_nnls_clear (&self->nnls);

  ncm_matrix_clear (&self->IM);
  ncm_matrix_clear (&self->sub_IM);
  ncm_vector_clear (&self->sub_x);
  ncm_vector_clear (&self->f);
  ncm_vector_clear (&self->f1);
  ncm_vector_clear (&self->cv_m2lnp);
  ncm_vector_clear (&self->cv_m2lnL_sample);
  ncm_vector_clear (&self->cv_w);

  ncm_rng_clear (&self->rng);

  g_clear_pointer (&self->m2lnp_sort, g_array_unref); /* LCOV_EXCL_BR_LINE */
  g_clear_pointer (&self->m2lnp, g_array_unref);      /* LCOV_EXCL_BR_LINE */

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_finalize (GObject *object)
{
  NcmStatsDist *sd                 = NCM_STATS_DIST (object);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->levmar_n = 0;
  g_clear_pointer (&self->levmar_workz, g_free);
  g_clear_pointer (&self->fmin, gsl_multimin_fminimizer_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_parent_class)->finalize (object);
}

static void _ncm_stats_dist_set_dim (NcmStatsDist *sd, const guint dim);
static gdouble _ncm_stats_dist_get_href (NcmStatsDist *sd);

/* LCOV_EXCL_START these should be overwritten and never executed */

static void
_ncm_stats_dist_prepare_kernel (NcmStatsDist *sd, GPtrArray *sample_array)
{
  g_error ("method prepare_kernel not implemented by %s.", G_OBJECT_TYPE_NAME (sd));
}

static void _ncm_stats_dist_prepare (NcmStatsDist *sd);
static void _ncm_stats_dist_prepare_interp (NcmStatsDist *sd, NcmVector *m2lnp);

static void
_ncm_stats_dist_compute_IM (NcmStatsDist *sd, NcmMatrix *IM)
{
  g_error ("method compute_IM not implemented by %s.", G_OBJECT_TYPE_NAME (sd));
}

static NcmMatrix *
_ncm_stats_dist_peek_cov_decomp (NcmStatsDist *sd, guint i)
{
  g_error ("method peek_cov_decomp not implemented by %s.", G_OBJECT_TYPE_NAME (sd));

  return NULL;
}

static NcmMatrix *
_ncm_stats_dist_peek_full_cov_decomp (NcmStatsDist *sd)
{
  g_error ("method peek_full_cov_decomp not implemented by %s.", G_OBJECT_TYPE_NAME (sd));

  return NULL;
}

static NcmMatrix *
_ncm_stats_dist_peek_full_cov (NcmStatsDist *sd)
{
  g_error ("method peek_full_cov not implemented by %s.", G_OBJECT_TYPE_NAME (sd));

  return NULL;
}

static gdouble
_ncm_stats_dist_get_lnnorm (NcmStatsDist *sd, guint i)
{
  g_error ("method get_lnnorm not implemented by %s.", G_OBJECT_TYPE_NAME (sd));

  return 0.0;
}

static gdouble
_ncm_stats_dist_eval_weights (NcmStatsDist *sd, NcmVector *weights, NcmVector *x)
{
  g_error ("method eval_weights not implemented by %s.", G_OBJECT_TYPE_NAME (sd));

  return 0.0;
}

static gdouble
_ncm_stats_dist_eval_weights_m2lnp (NcmStatsDist *sd, NcmVector *weights, NcmVector *x)
{
  g_error ("method eval_weights_m2lnp not implemented by %s.", G_OBJECT_TYPE_NAME (sd));

  return 0.0;
}

static void _ncm_stats_dist_reset (NcmStatsDist *sd);
static void _ncm_stats_dist_update_centers (NcmStatsDist *sd);

/* LCOV_EXCL_STOP */

static void
ncm_stats_dist_class_init (NcmStatsDistClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_stats_dist_set_property;
  object_class->get_property = &_ncm_stats_dist_get_property;
  object_class->dispose      = &_ncm_stats_dist_dispose;
  object_class->finalize     = &_ncm_stats_dist_finalize;

  g_object_class_install_property (object_class,
                                   PROP_KERNEL,
                                   g_param_spec_object ("kernel",
                                                        NULL,
                                                        "Interpolating kernel",
                                                        NCM_TYPE_STATS_DIST_KERNEL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SAMPLE_SIZE,
                                   g_param_spec_uint ("N",
                                                      NULL,
                                                      "sample size",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READABLE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_OVER_SMOOTH,
                                   g_param_spec_double ("over-smooth",
                                                        NULL,
                                                        "Over-smooth distribution",
                                                        1.0e-5, G_MAXDOUBLE, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_CV_TYPE,
                                   g_param_spec_enum ("CV-type",
                                                      NULL,
                                                      "Cross-validation method",
                                                      NCM_TYPE_STATS_DIST_CV, NCM_STATS_DIST_CV_NONE,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_USE_THREADS,
                                   g_param_spec_boolean ("use-threads",
                                                         NULL,
                                                         "Whether to use OpenMP threads during computation",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_SPLIT_FRAC,
                                   g_param_spec_double ("split-frac",
                                                        NULL,
                                                        "Fraction to use in the split cross-validation",
                                                        0.10, 0.95, 0.5,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_PRINT_FIT,
                                   g_param_spec_boolean ("print-fit",
                                                         NULL,
                                                         "Whether to print the fitting process",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmStatsDist:center-shrink:
   *
   * Whether to shrink the kernel centers toward the sample mean so that the
   * covariance of the kernel mixture matches the sample covariance for any
   * bandwidth, see the class description. Default: FALSE.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_CENTER_SHRINK,
                                   g_param_spec_boolean ("center-shrink",
                                                         NULL,
                                                         "Whether to shrink the kernel centers toward the sample mean",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmStatsDist:auto-kernel:
   *
   * Whether the kernel is chosen together with the over-smooth factor, by the same
   * held-out objective. Requires a cross-validation that fits the bandwidth, currently
   * #NCM_STATS_DIST_CV_SPLIT_NOFIT; it is ignored otherwise. The kernel is a
   * #NcmStatsDistKernelST whose degrees of freedom $\nu$ are fitted jointly with the
   * over-smooth factor, over $\nu \in [\nu_\mathrm{min}, 10^4]$, with
   * $\nu_\mathrm{min} = 2.5$ when #NcmStatsDist:center-shrink is set and $1$
   * otherwise. Reaching the upper bound installs a #NcmStatsDistKernelGauss. See the
   * class description. Default: FALSE.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_AUTO_KERNEL,
                                   g_param_spec_boolean ("auto-kernel",
                                                         NULL,
                                                         "Whether to choose the kernel with the bandwidth",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmStatsDist:defensive-frac:
   *
   * Weight $\epsilon$ of a wide Student-t component mixed into the proposal,
   * $q = (1 - \epsilon)\, q_\mathrm{mixture} + \epsilon\, t_\nu(\mu, c\,C)$, with $\mu$
   * and $C$ the sample mean and covariance, $c$ #NcmStatsDist:defensive-scale and $\nu$
   * #NcmStatsDist:defensive-nu. It bounds the proposal density from below where the
   * kernels leave holes, so that a walker there can still move. Zero disables it and
   * leaves every evaluation and draw unchanged. Default: 0.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DEFENSIVE_FRAC,
                                   g_param_spec_double ("defensive-frac",
                                                        NULL,
                                                        "Weight of the wide Student-t component in the proposal",
                                                        0.0, 1.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmStatsDist:defensive-scale:
   *
   * Factor $c$ multiplying the sample covariance in the wide component of
   * #NcmStatsDist:defensive-frac. Default: 4.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DEFENSIVE_SCALE,
                                   g_param_spec_double ("defensive-scale",
                                                        NULL,
                                                        "Covariance factor of the wide component",
                                                        1.0e-2, 1.0e4, 4.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmStatsDist:defensive-nu:
   *
   * Degrees of freedom $\nu$ of the wide Student-t component of
   * #NcmStatsDist:defensive-frac. Default: 3.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_DEFENSIVE_NU,
                                   g_param_spec_double ("defensive-nu",
                                                        NULL,
                                                        "Degrees of freedom of the wide component",
                                                        1.0, 1.0e4, 3.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmStatsDist:uniform-weights:
   *
   * Whether ncm_stats_dist_prepare_interp() keeps uniform kernel weights instead of
   * fitting them by NNLS. The sample's $-2\ln L$ is then used only by the bandwidth
   * objectives that need it and by the outlier cut. Default: FALSE.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_UNIFORM_WEIGHTS,
                                   g_param_spec_boolean ("uniform-weights",
                                                         NULL,
                                                         "Keep uniform kernel weights in prepare_interp instead of the NNLS fit",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  klass->set_dim              = &_ncm_stats_dist_set_dim;
  klass->get_href             = &_ncm_stats_dist_get_href;
  klass->prepare_kernel       = &_ncm_stats_dist_prepare_kernel;
  klass->prepare              = &_ncm_stats_dist_prepare;
  klass->prepare_interp       = &_ncm_stats_dist_prepare_interp;
  klass->compute_IM           = &_ncm_stats_dist_compute_IM;
  klass->peek_cov_decomp      = &_ncm_stats_dist_peek_cov_decomp;
  klass->peek_full_cov_decomp = &_ncm_stats_dist_peek_full_cov_decomp;
  klass->peek_full_cov        = &_ncm_stats_dist_peek_full_cov;
  klass->get_lnnorm           = &_ncm_stats_dist_get_lnnorm;
  klass->eval_weights         = &_ncm_stats_dist_eval_weights;
  klass->eval_weights_m2lnp   = &_ncm_stats_dist_eval_weights_m2lnp;
  klass->reset                = &_ncm_stats_dist_reset;
  klass->update_centers       = &_ncm_stats_dist_update_centers;
}

static void
_ncm_stats_dist_set_dim (NcmStatsDist *sd, const guint dim)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->d = dim;
  g_ptr_array_set_size (self->center_array, 0);
  ncm_vector_clear (&self->center_mean);
}

static void
_ncm_stats_dist_update_centers (NcmStatsDist *sd)
{
  /* Nothing to do: subclasses that cache transformed centers override this. */
}

#define NCM_STATS_DIST_CENTER_NEARPD_MAXITER (200)

static void _ncm_stats_dist_update_defensive (NcmStatsDist *sd);
gdouble _ncm_stats_dist_accept (const gsl_vector *v, void *params);
gdouble _ncm_stats_dist_loo_m2lnp (const gsl_vector *v, void *params);

static inline gdouble
_ncm_stats_dist_logaddexp (const gdouble a, const gdouble b)
{
  if (a == GSL_NEGINF)
    return b;

  if (b == GSL_NEGINF)
    return a;

  return GSL_MAX (a, b) + log1p (exp (-fabs (a - b)));
}

static void
_ncm_stats_dist_cholesky_upper (NcmMatrix *U, NcmMatrix *M)
{
  ncm_matrix_memcpy (U, M);

  if (ncm_matrix_cholesky_decomp (U, 'U') != 0)
  {
    ncm_matrix_memcpy (U, M);

    if (ncm_matrix_nearPD (U, 'U', TRUE, NCM_STATS_DIST_CENTER_NEARPD_MAXITER) != 0)
      g_error ("_ncm_stats_dist_cholesky_upper: matrix is not positive definite.");
  }
}

/* Cholesky leaves the strict lower triangle untouched; BLAS products need it zero. */
void
_ncm_stats_dist_zero_strict_lower (NcmMatrix *U)
{
  const guint d = ncm_matrix_nrows (U);
  guint a, b;

  for (a = 1; a < d; a++)
    for (b = 0; b < a; b++)
      ncm_matrix_set (U, a, b, 0.0);
}

/*
 * The two d x d matrices a subclass hands over in prepare_kernel(): the upper
 * Cholesky factor of the sample covariance and the mean kernel scale matrix.
 */
void
_ncm_stats_dist_center_matrices (NcmStatsDist *sd, NcmMatrix **C_decomp, NcmMatrix **mean_cov)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  if ((self->center_C_decomp == NULL) || (ncm_matrix_nrows (self->center_C_decomp) != self->d))
  {
    ncm_matrix_clear (&self->center_C_decomp);
    ncm_matrix_clear (&self->center_mean_cov);
    self->center_C_decomp = ncm_matrix_new (self->d, self->d);
    self->center_mean_cov = ncm_matrix_new (self->d, self->d);
  }

  C_decomp[0] = self->center_C_decomp;
  mean_cov[0] = self->center_mean_cov;
}

gboolean
_ncm_stats_dist_center_transform_is_identity (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->center_Ahat_identity;
}

/*
 * U <- upper Cholesky factor of Ahat Sigma Ahat^T, with Sigma = U0^T U0 the untransformed
 * kernel scale matrix. Ahat has unit determinant, so the kernel normalization is unchanged
 * up to rounding. When the transform is the identity the factor is copied unchanged.
 */
void
_ncm_stats_dist_refactor_decomp (NcmStatsDist *sd, NcmMatrix *U0, NcmMatrix *U)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  if (self->center_Ahat_identity)
  {
    ncm_matrix_memcpy (U, U0);

    return;
  }

  {
    NcmMatrix *B = ncm_matrix_dup (U0);
    NcmMatrix *M = ncm_matrix_new (self->d, self->d);

    _ncm_stats_dist_zero_strict_lower (B);
    ncm_matrix_dgemm (M, 'N', 'T', 1.0, B, self->center_Ahat, 0.0); /* M = U0 Ahat^T      */
    ncm_matrix_dgemm (B, 'T', 'N', 1.0, M, M, 0.0);                 /* B = M^T M = Ahat Sigma Ahat^T */
    _ncm_stats_dist_cholesky_upper (U, B);

    ncm_matrix_free (B);
    ncm_matrix_free (M);
  }
}

/*
 * Sets the bandwidth and recomputes the kernel centers. With center shrinkage the
 * centers are c_i = mu + A (x_i - mu) and the kernel scale matrices Ahat Sigma_i Ahat^T,
 * where A = a Ahat, det Ahat = 1, and A solves A (C + kappa h^2 <Sigma>) A^T = C for the
 * sample covariance C = U_C^T U_C and the mean kernel scale matrix <Sigma>, both handed
 * over by the subclass in prepare_kernel(). The bandwidth stored and applied is a h.
 * Without center shrinkage A is the identity. Every place that changes href must go
 * through here so that centers, factors and bandwidth stay consistent.
 */
static void
_ncm_stats_dist_set_href (NcmStatsDist *sd, const gdouble href)
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  const guint d                    = self->d;
  gdouble a                        = 1.0;
  guint i;

  g_assert_cmpuint (self->sample_array->len, >=, self->n_kernels);

  if ((self->center_A == NULL) || (ncm_matrix_nrows (self->center_A) != d))
  {
    ncm_matrix_clear (&self->center_A);
    ncm_matrix_clear (&self->center_Ahat);
    self->center_A    = ncm_matrix_new (d, d);
    self->center_Ahat = ncm_matrix_new (d, d);
  }

  if (self->center_shrink)
  {
    const gdouble kappa = ncm_stats_dist_kernel_get_var_factor (self->kernel);
    NcmMatrix *UC       = NULL;
    NcmMatrix *M        = ncm_matrix_new (d, d);
    NcmMatrix *UM       = ncm_matrix_new (d, d);
    gint ret;

    if (!gsl_finite (kappa))
      g_error ("_ncm_stats_dist_set_href: center shrinkage requires a kernel with a finite "
               "covariance, but %s has none (a Student-t kernel needs nu > 2). Use a "
               "Gaussian kernel or a larger nu, or disable center-shrink.",
               G_OBJECT_TYPE_NAME (self->kernel));

    g_assert (self->center_C_decomp != NULL);
    g_assert (self->center_mean_cov != NULL);

    UC = ncm_matrix_dup (self->center_C_decomp);
    _ncm_stats_dist_zero_strict_lower (UC);

    /* M = C + kappa h^2 <Sigma> = U_M^T U_M */
    ncm_matrix_dgemm (M, 'T', 'N', 1.0, UC, UC, 0.0);
    ncm_matrix_add_mul (M, kappa * href * href, self->center_mean_cov);
    _ncm_stats_dist_cholesky_upper (UM, M);

    /* A = U_C^T U_M^{-T}: solve X U_M^T = U_C^T. */
    gsl_matrix_transpose_memcpy (ncm_matrix_gsl (self->center_A), ncm_matrix_gsl (UC));
    ret = gsl_blas_dtrsm (CblasRight, CblasUpper, CblasTrans, CblasNonUnit,
                          1.0, ncm_matrix_gsl (UM), ncm_matrix_gsl (self->center_A));
    NCM_TEST_GSL_RESULT ("_ncm_stats_dist_set_href", ret);

    /* a = det(A)^(1/d); Ahat = A / a, snapped to the identity when it is one to rounding. */
    a = exp ((ncm_matrix_cholesky_lndet (UC) - ncm_matrix_cholesky_lndet (UM)) / (2.0 * d));
    ncm_matrix_memcpy (self->center_Ahat, self->center_A);
    gsl_matrix_scale (ncm_matrix_gsl (self->center_Ahat), 1.0 / a);
    {
      gdouble dev = 0.0;
      guint p, q;

      for (p = 0; p < d; p++)
        for (q = 0; q < d; q++)
          dev = GSL_MAX (dev, fabs (ncm_matrix_get (self->center_Ahat, p, q) - ((p == q) ? 1.0 : 0.0)));

      self->center_Ahat_identity = (dev < 1.0e-12);

      if (self->center_Ahat_identity)
        gsl_matrix_set_identity (ncm_matrix_gsl (self->center_Ahat));
    }

    ncm_matrix_free (UC);
    ncm_matrix_free (M);
    ncm_matrix_free (UM);
  }
  else
  {
    gsl_matrix_set_identity (ncm_matrix_gsl (self->center_A));
    gsl_matrix_set_identity (ncm_matrix_gsl (self->center_Ahat));
    self->center_Ahat_identity = TRUE;
  }

  self->href     = a * href;
  self->center_a = a;

  if ((self->center_mean == NULL) || (ncm_vector_len (self->center_mean) != d))
  {
    ncm_vector_clear (&self->center_mean);
    self->center_mean = ncm_vector_new (d);
  }

  {
    const guint cur_len = self->center_array->len;

    g_ptr_array_set_size (self->center_array, self->n_kernels);

    for (i = cur_len; i < self->n_kernels; i++)
      g_ptr_array_index (self->center_array, i) = ncm_vector_new (d);
  }

  ncm_vector_set_zero (self->center_mean);

  for (i = 0; i < self->n_kernels; i++)
    ncm_vector_add (self->center_mean, g_ptr_array_index (self->sample_array, i));

  ncm_vector_scale (self->center_mean, 1.0 / (1.0 * self->n_kernels));

  if (self->center_shrink)
  {
    NcmVector *dx = ncm_vector_new (d);

    for (i = 0; i < self->n_kernels; i++)
    {
      NcmVector *x_i = g_ptr_array_index (self->sample_array, i);
      NcmVector *c_i = g_ptr_array_index (self->center_array, i);
      gint ret;

      ncm_vector_memcpy (dx, x_i);
      ncm_vector_sub (dx, self->center_mean);
      ret = gsl_blas_dgemv (CblasNoTrans, 1.0, ncm_matrix_gsl (self->center_A), ncm_vector_gsl (dx), 0.0, ncm_vector_gsl (c_i));
      NCM_TEST_GSL_RESULT ("_ncm_stats_dist_set_href", ret);
      ncm_vector_add (c_i, self->center_mean);
    }

    ncm_vector_free (dx);
  }
  else
  {
    for (i = 0; i < self->n_kernels; i++)
      ncm_vector_memcpy (g_ptr_array_index (self->center_array, i), g_ptr_array_index (self->sample_array, i));
  }

  sd_class->update_centers (sd);
  _ncm_stats_dist_update_defensive (sd);
}

/*
 * The wide component: a Student-t centered on the sample mean with c times the sample
 * covariance as scale matrix, weight epsilon in the proposal. Rebuilt with every bandwidth
 * since the sample mean and covariance factor are refreshed there.
 */
static void
_ncm_stats_dist_update_defensive (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  if (self->defensive_frac <= 0.0)
    return;

  g_assert (self->center_C_decomp != NULL);
  g_assert (self->center_mean != NULL);

  if ((self->defensive_kernel == NULL) || (ncm_stats_dist_kernel_get_dim (self->defensive_kernel) != self->d))
  {
    ncm_stats_dist_kernel_clear (&self->defensive_kernel);
    self->defensive_kernel = NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_st_new (self->d, self->defensive_nu));
  }
  else
  {
    ncm_stats_dist_kernel_st_set_nu (NCM_STATS_DIST_KERNEL_ST (self->defensive_kernel), self->defensive_nu);
  }

  if ((self->defensive_decomp == NULL) || (ncm_matrix_nrows (self->defensive_decomp) != self->d))
  {
    ncm_matrix_clear (&self->defensive_decomp);
    self->defensive_decomp = ncm_matrix_new (self->d, self->d);
  }

  ncm_matrix_memcpy (self->defensive_decomp, self->center_C_decomp);
  _ncm_stats_dist_zero_strict_lower (self->defensive_decomp);
  gsl_matrix_scale (ncm_matrix_gsl (self->defensive_decomp), sqrt (self->defensive_scale));
  self->defensive_lnnorm = ncm_stats_dist_kernel_get_lnnorm (self->defensive_kernel, self->defensive_decomp);
}

/* -2 ln of the wide component's density at x. */
static gdouble
_ncm_stats_dist_defensive_m2lnK (NcmStatsDist *sd, NcmVector *x)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  NcmVector *dx                    = ncm_vector_dup (x);
  gdouble chi2, m2lnK;
  gint ret;

  ncm_vector_sub (dx, self->center_mean);
  ret = gsl_blas_dtrsv (CblasUpper, CblasTrans, CblasNonUnit, ncm_matrix_gsl (self->defensive_decomp), ncm_vector_gsl (dx));
  NCM_TEST_GSL_RESULT ("_ncm_stats_dist_defensive_m2lnK", ret);
  chi2  = ncm_vector_dot (dx, dx);
  m2lnK = 2.0 * self->defensive_lnnorm - 2.0 * log (ncm_stats_dist_kernel_eval_unnorm (self->defensive_kernel, chi2));
  ncm_vector_free (dx);

  return m2lnK;
}

static gdouble
_ncm_stats_dist_get_href (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->over_smooth * ncm_stats_dist_kernel_get_rot_bandwidth (self->kernel, self->n_kernels);
}

gdouble
_ncm_stats_dist_m2lnp (const gsl_vector *v, void *params)
{
  NcmStatsDist *sd                 = NCM_STATS_DIST (params);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  const double lnos                = gsl_vector_get (v, 0);
  gdouble m2lnp                    = 0.0;
  gint i;

  self->over_smooth = exp (lnos);
  _ncm_stats_dist_set_href (sd, _ncm_stats_dist_calc_href (sd));

  /* Per-point values first, summed serially afterwards: a shared accumulator inside the
   * parallel loop was a data race, and even a reduction would round in thread order. */
  if ((self->cv_m2lnp == NULL) || (ncm_vector_len (self->cv_m2lnp) != self->n_obs))
  {
    ncm_vector_clear (&self->cv_m2lnp);
    self->cv_m2lnp = ncm_vector_new (self->n_obs);
  }

  #pragma omp parallel for if (self->use_threads)

  for (i = self->n_kernels; i < self->n_obs; i++)
  {
    NcmVector *x_i = g_ptr_array_index (self->sample_array, i);

    ncm_vector_set (self->cv_m2lnp, i, ncm_stats_dist_eval_m2lnp (sd, x_i));
  }

  for (i = self->n_kernels; i < (gint) self->n_obs; i++)
    m2lnp += ncm_vector_get (self->cv_m2lnp, i);

  if (self->print_fit)
    ncm_message ("# over-smooth: % 22.15g, m2lnp = % 22.15g\n",
                 self->over_smooth, m2lnp);

  return m2lnp;
}

/*
 * Held-out estimate of the independence-sampler acceptance (auto_tuning.md). With
 * r = ln q - ln pi at the kernel points (k) and at the held-out points (j), the mean over j
 * of E_{x' ~ q}[min (1, exp (r_j - r_{x'}))] is estimated by self-normalized importance
 * sampling with the kernel points as draws from pi and weights proportional to exp (r_k).
 * Returns -ln of the estimate. Needs the sample's -2ln(L) (cv_m2lnL_sample).
 */
gdouble
_ncm_stats_dist_accept (const gsl_vector *v, void *params)
{
  NcmStatsDist *sd                 = NCM_STATS_DIST (params);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  const double lnos                = gsl_vector_get (v, 0);
  const gint nk                    = self->n_kernels;
  const gint no                    = self->n_obs;
  gdouble lse_k                    = GSL_NEGINF;
  gdouble acc                      = 0.0;
  gint i, j;

  self->over_smooth = exp (lnos);
  _ncm_stats_dist_set_href (sd, _ncm_stats_dist_calc_href (sd));

  if ((self->cv_m2lnp == NULL) || (ncm_vector_len (self->cv_m2lnp) != self->n_obs))
  {
    ncm_vector_clear (&self->cv_m2lnp);
    self->cv_m2lnp = ncm_vector_new (self->n_obs);
  }

  /* r_i = ln q (x_i) - ln pi (x_i) = (m2lnL_i - m2lnq_i) / 2 at every sample point. */
  #pragma omp parallel for if (self->use_threads)

  for (i = 0; i < no; i++)
  {
    NcmVector *x_i = g_ptr_array_index (self->sample_array, i);

    ncm_vector_set (self->cv_m2lnp, i, 0.5 * (ncm_vector_get (self->cv_m2lnL_sample, i) - ncm_stats_dist_eval_m2lnp (sd, x_i)));
  }

  for (i = 0; i < nk; i++)
    lse_k = _ncm_stats_dist_logaddexp (lse_k, ncm_vector_get (self->cv_m2lnp, i));

  /* Effective sample size of the importance weights w_k = exp (r_k - lse_k): when a few
   * kernel points carry all the weight the estimate is meaningless and, left alone, the
   * optimizer runs to the largest bandwidth. Refuse such trials. */
  {
    gdouble lse2 = GSL_NEGINF;
    gdouble ess;

    for (i = 0; i < nk; i++)
      lse2 = _ncm_stats_dist_logaddexp (lse2, 2.0 * (ncm_vector_get (self->cv_m2lnp, i) - lse_k));

    ess = exp (-lse2);

    if (ess < GSL_MAX (10.0, 0.01 * nk))
    {
      if (self->print_fit)
        ncm_message ("# over-smooth: % 22.15g, importance ESS % 8.2f of %d: rejected\n", self->over_smooth, ess, nk);

      return GSL_POSINF;
    }
  }

  for (j = nk; j < no; j++)
  {
    const gdouble r_j = ncm_vector_get (self->cv_m2lnp, j);
    gdouble a_j       = GSL_NEGINF;

    for (i = 0; i < nk; i++)
    {
      const gdouble r_k = ncm_vector_get (self->cv_m2lnp, i);

      a_j = _ncm_stats_dist_logaddexp (a_j, (r_k - lse_k) + GSL_MIN (0.0, r_j - r_k));
    }

    acc += exp (a_j);
  }

  acc /= (no - nk);

  if (self->print_fit)
    ncm_message ("# over-smooth: % 22.15g, acceptance = % 22.15g\n", self->over_smooth, acc);

  return -log (acc);
}

/*
 * Leave-one-out likelihood cross-validation, -2 sum_i ln q_{-i} (x_i): every point is a
 * kernel and q_{-i} is the mixture with weight i zeroed and the rest rescaled.
 */
gdouble
_ncm_stats_dist_loo_m2lnp (const gsl_vector *v, void *params)
{
  NcmStatsDist *sd                 = NCM_STATS_DIST (params);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  const double lnos                = gsl_vector_get (v, 0);
  const gint n                     = self->n_kernels;
  gdouble m2lnp                    = 0.0;
  gint i;

  self->over_smooth = exp (lnos);
  _ncm_stats_dist_set_href (sd, _ncm_stats_dist_calc_href (sd));

  if ((self->cv_m2lnp == NULL) || (ncm_vector_len (self->cv_m2lnp) != self->n_obs))
  {
    ncm_vector_clear (&self->cv_m2lnp);
    self->cv_m2lnp = ncm_vector_new (self->n_obs);
  }

  if ((self->cv_w == NULL) || (ncm_vector_len (self->cv_w) != self->n_kernels))
  {
    ncm_vector_clear (&self->cv_w);
    self->cv_w = ncm_vector_new (self->n_kernels);
  }

  ncm_vector_memcpy (self->cv_w, self->weights);

  /* Serial: cv_w is modified in place, one weight at a time. */
  for (i = 0; i < n; i++)
  {
    NcmVector *x_i    = g_ptr_array_index (self->sample_array, i);
    const gdouble w_i = ncm_vector_get (self->weights, i);
    gdouble q_m;

    ncm_vector_set (self->cv_w, i, 0.0);
    q_m = sd_class->eval_weights_m2lnp (sd, self->cv_w, x_i) + 2.0 * log1p (-w_i);
    ncm_vector_set (self->cv_w, i, w_i);
    ncm_vector_set (self->cv_m2lnp, i, q_m);
  }

  for (i = 0; i < n; i++)
    m2lnp += ncm_vector_get (self->cv_m2lnp, i);

  if (self->print_fit)
    ncm_message ("# over-smooth: % 22.15g, loo m2lnp = % 22.15g\n", self->over_smooth, m2lnp);

  return m2lnp;
}

gdouble
_ncm_stats_dist_amise_kde_gauss (const gsl_vector *v, void *params)
{
  NcmStatsDist *sd                 = NCM_STATS_DIST (params);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  const double lnos                = gsl_vector_get (v, 0);
  gdouble amise                    = 0.0;
  guint i, j;

  self->over_smooth = exp (lnos);
  _ncm_stats_dist_set_href (sd, _ncm_stats_dist_calc_href (sd));
  self->href = sqrt (2.0) * self->href;

  sd_class->compute_IM (sd, self->IM);

  for (i = 0; i < self->n_kernels; i++)
  {
    for (j = 0; j < self->n_kernels; j++)
    {
      amise += ncm_matrix_get (self->IM, i, j) / gsl_pow_2 (self->n_kernels);
    }
  }

  self->over_smooth = exp (lnos);
  _ncm_stats_dist_set_href (sd, _ncm_stats_dist_calc_href (sd));

  sd_class->compute_IM (sd, self->IM);

  for (i = 0; i < self->n_kernels; i++)
  {
    for (j = 0; j < i; j++)
    {
      amise -= 2.0 * ncm_matrix_get (self->IM, i, j) / (self->n_kernels * (self->n_kernels - 1));
    }

    for (j = i + 1; j < self->n_kernels; j++)
    {
      amise -= 2.0 * ncm_matrix_get (self->IM, i, j) / (self->n_kernels * (self->n_kernels - 1));
    }
  }

  if (self->print_fit)
    ncm_message ("# over-smooth: % 22.15g, amise = % 22.15g\n",
                 self->over_smooth, amise);

  return amise;
}

void ncm_stats_dist_sample2 (NcmStatsDist *sd, NcmVector *x1, NcmVector *x2, NcmRNG *rng);

gdouble
_ncm_stats_dist_amise (const gsl_vector *v, void *params)
{
  NcmStatsDist *sd                 = NCM_STATS_DIST (params);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  const double lnos                = gsl_vector_get (v, 0);
  gdouble amise                    = 0.0;
  guint i, j;

  self->over_smooth = exp (lnos);
  _ncm_stats_dist_set_href (sd, _ncm_stats_dist_calc_href (sd));

  sd_class->compute_IM (sd, self->IM);

  g_array_set_size (self->m2lnp_sort, self->n_kernels);
  g_array_set_size (self->m2lnp, self->n_kernels);

  for (i = 0; i < self->n_kernels; i++)
  {
    register gdouble row_sum = 0.0;

    for (j = 0; j < i; j++)
    {
      row_sum += ncm_matrix_get (self->IM, i, j);
    }

    for (j = i + 1; j < self->n_kernels; j++)
    {
      row_sum += ncm_matrix_get (self->IM, i, j);
    }

    amise -= 2.0 * row_sum / (self->n_kernels * (self->n_kernels - 1));

    g_array_index (self->m2lnp, gdouble, i) = (row_sum + ncm_matrix_get (self->IM, i, i)) / self->n_kernels;
  }

  gsl_sort_index (&g_array_index (self->m2lnp_sort, size_t, 0),
                  (gdouble *) self->m2lnp->data, 1, self->m2lnp->len);

  {
    NcmRNG *rng        = ncm_rng_seeded_new (NULL, 0);
    NcmVector *x1      = ncm_vector_new (self->d);
    NcmVector *x2      = ncm_vector_new (self->d);
    NcmStatsVec *stats = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
    guint max_iter     = 100000000;
    gdouble mean;
    gdouble p1, p2;

    for (i = 0; i < 100; i++)
    {
      ncm_stats_dist_sample2 (sd, x1, x2, rng);
      p1 = ncm_stats_dist_eval (sd, x1);
      p2 = ncm_stats_dist_eval (sd, x2);

      ncm_stats_vec_set (stats, 0, p1);
      ncm_stats_vec_set (stats, 1, p2);
      ncm_stats_vec_update (stats);
    }

    for (i = 0; i < max_iter; i++)
    {
      ncm_stats_dist_sample2 (sd, x1, x2, rng);
      p1 = ncm_stats_dist_eval (sd, x1);
      p2 = ncm_stats_dist_eval (sd, x2);

      ncm_stats_vec_set (stats, 0, p1);
      ncm_stats_vec_set (stats, 1, p2);
      ncm_stats_vec_update (stats);

      mean = 0.5 * (ncm_stats_vec_get_mean (stats, 0) + ncm_stats_vec_get_mean (stats, 1));

      {
        const gdouble var = 0.25 * (ncm_stats_vec_get_var (stats, 0) + ncm_stats_vec_get_var (stats, 1) + 2.0 * ncm_stats_vec_get_cov (stats, 0, 1));

        const gdouble msd = sqrt (var / (i + 101.0)) / mean;

        if (msd < 1.0e-2)
          break;
      }
    }

    amise += mean;

    ncm_stats_vec_free (stats);
    ncm_vector_free (x1);
    ncm_vector_free (x2);
    ncm_rng_clear (&rng);
  }

  if (self->print_fit)
    ncm_message ("# over-smooth: % 22.15g, amise = % 22.15g\n",
                 self->over_smooth, amise);

  return amise;
}

static gdouble
_ncm_stats_dist_minimize_obj (NcmStatsDist *sd, gdouble (*objective) (const gsl_vector *, void *))
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  gdouble s                        = 0.1;
  gdouble lnos                     = log (self->over_smooth);
  NcmVector *x                     = ncm_vector_new_data_static (&lnos, 1, 1);
  NcmVector *ss                    = ncm_vector_new_data_static (&s, 1, 1);
  gsl_multimin_function minex_func;

  minex_func.n      = 1;
  minex_func.f      = objective;
  minex_func.params = sd;

  gsl_multimin_fminimizer_set (
    self->fmin,
    &minex_func,
    ncm_vector_gsl (x),
    ncm_vector_gsl (ss));

  {
    gint iter = 0;
    gint status;

    do {
      iter++;
      status = gsl_multimin_fminimizer_iterate (self->fmin);

      if (status)
        break;

      status = gsl_multimin_test_size (self->fmin->size, 1.0e-3);
    } while (status == GSL_CONTINUE && iter < 1000);

    if (self->print_fit)
      printf ("# iter: %d, over-smooth: % 22.15g, m2lnp = % 22.15g, gsl status (%d)\n",
              iter, self->over_smooth, self->fmin->fval, status);
  }

  ncm_vector_free (x);
  ncm_vector_free (ss);

  return gsl_multimin_fminimizer_minimum (self->fmin);
}

/*
 * Chooses the kernel together with the over-smooth factor, using the same held-out
 * objective. The Student-t kernel is (1 + chi2 / nu)^(-(nu + d) / 2), which is the Cauchy
 * kernel at nu = 1 and tends to the Gaussian one as nu grows, so the kernel is not a
 * discrete choice but the single continuous parameter nu. What is fitted is therefore
 * (ln over_smooth, ln nu) jointly, by the same objective and on the same held-out points.
 *
 * The range of nu is bounded on both sides. Center shrinkage needs a finite kernel
 * covariance, nu / (nu - 2), so nu > 2 there. At the other end the kernel is the Gaussian
 * one to better than 1e-4 by nu ~ 1e4, while the (1 + chi2 / nu) form starts losing
 * precision well above that, so 1e4 is the ceiling and the Gaussian kernel itself is
 * installed whenever the fit reaches it.
 */

#define NCM_STATS_DIST_AUTO_KERNEL_NU_MAX (1.0e4)

typedef struct _NcmStatsDistKernelFit
{
  NcmStatsDist *sd;

  gdouble (*objective) (const gsl_vector *, void *);

  gdouble nu_min;
  guint neval;
} NcmStatsDistKernelFit;

static gdouble
_ncm_stats_dist_m2lnp_kernel (const gsl_vector *v, void *params)
{
  NcmStatsDistKernelFit *fit       = (NcmStatsDistKernelFit *) params;
  NcmStatsDist *sd                 = fit->sd;
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  const gdouble nu                 = GSL_MIN (GSL_MAX (exp (gsl_vector_get (v, 1)), fit->nu_min),
                                              NCM_STATS_DIST_AUTO_KERNEL_NU_MAX);
  gsl_vector_const_view os_view = gsl_vector_const_subvector (v, 0, 1);

  ncm_stats_dist_kernel_st_set_nu (NCM_STATS_DIST_KERNEL_ST (self->kernel), nu);
  sd_class->update_kernel_norms (sd);

  /* The over-smooth part of the objective is unchanged; it reads element zero. */
  return fit->objective (&os_view.vector, sd);
}

static gdouble
_ncm_stats_dist_m2lnp_kernel_nlopt (guint n, const gdouble *x, gdouble *grad, gpointer params)
{
  NcmVector *v               = ncm_vector_new_data_static ((gdouble *) x, n, 1);
  NcmStatsDistKernelFit *fit = (NcmStatsDistKernelFit *) params;
  const gdouble res          = _ncm_stats_dist_m2lnp_kernel (ncm_vector_gsl (v), params);

  fit->neval++;
  g_assert (grad == NULL);
  ncm_vector_free (v);

  return res;
}

/*
 * One (ln over_smooth, ln nu) fit, at the center-shrinkage setting the caller has left in
 * place. Returns the objective it reached and leaves the object holding that solution.
 */
static gdouble
_ncm_stats_dist_fit_kernel (NcmStatsDist *sd, gdouble (*objective) (const gsl_vector *, void *),
                            NcmStatsDistKernelST *st, gdouble *os_best, gdouble *nu_best)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  const gdouble nu_min             = self->center_shrink ? 2.5 : 1.0;
  const gdouble nu_start           = GSL_MAX (nu_min, 10.0);
  NcmStatsDistKernelFit fit        = {sd, objective, nu_min, 0};
  gdouble p[2]                     = {log (self->over_smooth), log (nu_start)};
  gdouble step[2]                  = {0.1, 0.5};
  gdouble lb[2]                    = {log (1.0e-2), log (nu_min)};
  gdouble ub[2]                    = {log (2.0e1), log (NCM_STATS_DIST_AUTO_KERNEL_NU_MAX)};
  nlopt_opt opt                    = nlopt_create (NLOPT_LN_BOBYQA, 2);
  gdouble minf                     = 0.0;
  gint ret;

  ncm_stats_dist_kernel_st_set_nu (st, nu_start);

  /*
   * BOBYQA builds a quadratic model of a smooth objective and takes the bounds natively,
   * which keeps nu inside the range where the kernel is both defined and numerically
   * sound, without the objective having to clamp its own argument. The simplex methods
   * can terminate on the nu lower bound instead of the interior minimum. Comparison
   * against Nelder-Mead, Subplex and COBYLA: dev-notes/apes_center_shrink/auto_tuning.md.
   */
  nlopt_set_lower_bounds (opt, lb);
  nlopt_set_upper_bounds (opt, ub);
  nlopt_set_initial_step (opt, step);
  nlopt_set_xtol_rel (opt, 1.0e-3);
  nlopt_set_maxeval (opt, 1000);
  nlopt_set_min_objective (opt, &_ncm_stats_dist_m2lnp_kernel_nlopt, &fit);

  ret = nlopt_optimize (opt, p, &minf);

  if (ret < 0) /* LCOV_EXCL_BR_LINE */
    g_warning ("_ncm_stats_dist_select_kernel: nlopt failed with code %d, "
               "keeping the configuration it reached.", ret);

  /* Re-evaluating at the minimizer leaves the object holding the chosen configuration. */
  {
    NcmVector *x_best = ncm_vector_new_data_static (p, 2, 1);

    _ncm_stats_dist_m2lnp_kernel (ncm_vector_gsl (x_best), &fit);
    ncm_vector_free (x_best);
  }

  if (self->print_fit)
    printf ("# center shrinkage: %s, over-smooth: % 22.15g, nu: % 22.15g, m2lnp = % 22.15g, "
            "neval = %u, nlopt status (%d)\n",
            self->center_shrink ? "on" : "off", self->over_smooth,
            ncm_stats_dist_kernel_st_get_nu (st), minf, fit.neval, ret);

  nlopt_destroy (opt);

  os_best[0] = self->over_smooth;
  nu_best[0] = ncm_stats_dist_kernel_st_get_nu (st);

  return minf;
}

/*
 * Center shrinkage is not fitted here; it stays the caller's choice. The held-out
 * objective is a Kullback-Leibler criterion, and KL(pi || p~) penalizes a proposal that is
 * too narrow far more than one that is too wide, so it does not see the cost a mixture
 * wider by (1 + kappa h^2 s^2) imposes on the Metropolis-Hastings acceptance. Fitting it
 * needs an objective that penalizes both directions, such as the held-out acceptance
 * estimate in dev-notes/apes_center_shrink/auto_tuning.md.
 */
static void
_ncm_stats_dist_select_kernel (NcmStatsDist *sd, gdouble (*objective) (const gsl_vector *, void *))
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistKernelST *st         = ncm_stats_dist_kernel_st_new (self->d, 10.0);
  NcmStatsDistKernel *orig_kernel  = ncm_stats_dist_kernel_ref (self->kernel);
  gdouble best_os                  = self->over_smooth;
  gdouble best_nu                  = 10.0;

  ncm_stats_dist_kernel_clear (&self->kernel);
  self->kernel = NCM_STATS_DIST_KERNEL (st);

  _ncm_stats_dist_fit_kernel (sd, objective, st, &best_os, &best_nu);

  self->over_smooth = best_os;

  if (best_nu >= 0.99 * NCM_STATS_DIST_AUTO_KERNEL_NU_MAX)
  {
    NcmStatsDistKernelGauss *gauss = ncm_stats_dist_kernel_gauss_new (self->d);

    ncm_stats_dist_kernel_clear (&self->kernel);
    self->kernel = NCM_STATS_DIST_KERNEL (gauss);
  }
  else
  {
    ncm_stats_dist_kernel_st_set_nu (st, best_nu);
  }

  sd_class->update_kernel_norms (sd);
  _ncm_stats_dist_set_href (sd, _ncm_stats_dist_calc_href (sd));

  ncm_stats_dist_kernel_free (orig_kernel);
}

static void
_ncm_stats_dist_prepare (NcmStatsDist *sd)
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  switch (self->cv_type)
  {
    case NCM_STATS_DIST_CV_LOO:

      self->n_obs     = self->sample_array->len;
      self->n_kernels = self->sample_array->len;

      if ((self->n_obs != self->alloc_n_obs) || (self->n_kernels != self->alloc_n_kernels))
      {
        ncm_matrix_clear (&self->IM);
        ncm_vector_clear (&self->f);
        ncm_vector_clear (&self->f1);

        self->IM = ncm_matrix_new (self->n_obs, self->n_kernels);
        self->f  = ncm_vector_new (self->n_obs);
        self->f1 = ncm_vector_new (self->n_obs);

        ncm_vector_set_all (self->f1, 1.0);

        self->alloc_n_obs     = self->n_obs;
        self->alloc_n_kernels = self->n_kernels;
        self->alloc_subs      = FALSE;
      }

      break;
    case NCM_STATS_DIST_CV_NONE:
      self->n_obs     = self->sample_array->len;
      self->n_kernels = self->sample_array->len;
      break;
    case NCM_STATS_DIST_CV_SPLIT:
    case NCM_STATS_DIST_CV_SPLIT_NOFIT:
    case NCM_STATS_DIST_CV_SPLIT_ACCEPT:
      self->n_obs     = self->sample_array->len;
      self->n_kernels = ceil (self->sample_array->len * self->split_frac);
      break;
    case NCM_STATS_DIST_CV_LOO_M2LNP:
      self->n_obs     = self->sample_array->len;
      self->n_kernels = self->sample_array->len;
      break;
    default: /* LCOV_EXCL_BR_LINE */
      g_assert_not_reached ();
      break;
  }

  if (self->n_obs <= self->d)
    g_error ("_ncm_stats_dist_prepare: the sample is too small.");

  sd_class->prepare_kernel (sd, self->sample_array);

  if ((self->weights == NULL) ||
      (self->n_kernels != ncm_vector_len (self->weights)))
  {
    ncm_vector_clear (&self->weights);
    ncm_vector_clear (&self->wcum);

    self->weights    = ncm_vector_new (self->n_kernels);
    self->wcum       = ncm_vector_new (self->n_kernels + 1);
    self->alloc_subs = FALSE;
  }

  _ncm_stats_dist_set_href (sd, _ncm_stats_dist_calc_href (sd));

  ncm_vector_set_all (self->weights, 1.0 / (1.0 * self->n_kernels));
  self->wcum_ready = FALSE;

  switch (self->cv_type)
  {
    case NCM_STATS_DIST_CV_NONE:
    case NCM_STATS_DIST_CV_SPLIT:
      break;
    case NCM_STATS_DIST_CV_SPLIT_NOFIT:

      if (self->auto_kernel)
        _ncm_stats_dist_select_kernel (sd, &_ncm_stats_dist_m2lnp);
      else
        _ncm_stats_dist_minimize_obj (sd, &_ncm_stats_dist_m2lnp);

      break;
    case NCM_STATS_DIST_CV_LOO:

      if ((NCM_IS_STATS_DIST_KDE (sd)) && (NCM_IS_STATS_DIST_KERNEL_GAUSS (self->kernel)))
        _ncm_stats_dist_minimize_obj (sd, &_ncm_stats_dist_amise_kde_gauss);
      else
        _ncm_stats_dist_minimize_obj (sd, &_ncm_stats_dist_amise);

      break;
    case NCM_STATS_DIST_CV_SPLIT_ACCEPT:

      if (self->cv_m2lnL_sample == NULL)
        g_error ("_ncm_stats_dist_prepare: NCM_STATS_DIST_CV_SPLIT_ACCEPT needs the sample's -2ln(L), "
                 "use ncm_stats_dist_prepare_interp().");

      if (self->n_kernels >= self->n_obs)
        g_error ("_ncm_stats_dist_prepare: NCM_STATS_DIST_CV_SPLIT_ACCEPT needs held-out points, "
                 "split-frac %g leaves none of %u.", self->split_frac, self->n_obs);

      if (self->auto_kernel)
        _ncm_stats_dist_select_kernel (sd, &_ncm_stats_dist_accept);
      else
        _ncm_stats_dist_minimize_obj (sd, &_ncm_stats_dist_accept);

      break;
    case NCM_STATS_DIST_CV_LOO_M2LNP:

      if (self->auto_kernel)
        _ncm_stats_dist_select_kernel (sd, &_ncm_stats_dist_loo_m2lnp);
      else
        _ncm_stats_dist_minimize_obj (sd, &_ncm_stats_dist_loo_m2lnp);

      break;
    default: /* LCOV_EXCL_BR_LINE */
      g_assert_not_reached ();
      break;
  }
}

static void
_ncm_stats_dist_compute_IM_full (NcmStatsDist *sd)
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  gint i;

  sd_class->compute_IM (sd, self->IM);

  #pragma omp parallel for if (self->use_threads)

  for (i = 0; i < self->n_obs; i++)
    ncm_matrix_mul_row (self->IM, i, 1.0 / ncm_vector_get (self->f, i));
}

typedef struct _NcmStatsDistEval
{
  NcmStatsDist *sd;
  NcmStatsDistPrivate * const self;
  NcmStatsDistClass *sd_class;
  NcmVector *residuals;
  NcmVector *m2lnp;
} NcmStatsDistEval;

static void
_ncm_stats_dist_prepare_interp_fit_nnls_f (gdouble *p, gdouble *hx, gint m, gint n, gpointer adata)
{
  NcmStatsDistEval *eval = adata;
  NcmVector *f           = ncm_vector_new_data_static (hx, n, 1);
  gdouble rnorm          = 0.0;
  gint i;

  g_assert (eval->self->n_obs == (guint) n);

  eval->self->over_smooth = exp (p[0]);
  _ncm_stats_dist_set_href (eval->sd, _ncm_stats_dist_calc_href (eval->sd));

  _ncm_stats_dist_compute_IM_full (eval->sd);
  rnorm = NCM_NNLS_SOLVE (eval->self->nnls, eval->self->sub_IM, eval->self->sub_x, eval->self->f1);

  #pragma omp parallel for if (eval->self->use_threads)

  for (i = 0; i < eval->self->n_obs; i++)
  {
    NcmVector *x_i         = g_ptr_array_index (eval->self->sample_array, i);
    const gdouble m2lnpt_i = ncm_vector_get (eval->m2lnp, i) - eval->self->min_m2lnp;
    const gdouble m2lnpi_i = ncm_stats_dist_eval_m2lnp (eval->sd, x_i);

    /*printf ("%d %d % 22.15g % 22.15g\n", eval->self->n, i, m2lnpt_i, m2lnpi_i);*/
    /*ncm_vector_set (f, i, sqrt (fabs (m2lnpt_i - m2lnpi_i))); */
    ncm_vector_set (f, i, expm1 (-0.5 * (m2lnpi_i - m2lnpt_i)));
  }

  if (eval->self->print_fit)
    ncm_message ("# over-smooth: % 22.15g, rnorm = % 22.15g, fnorm = % 22.15g\n",
                 eval->self->over_smooth, rnorm, ncm_vector_dnrm2 (f));

  ncm_vector_free (f);
}

static void
_ncm_stats_dist_alloc_nnls (NcmStatsDist *sd, const guint nrows, const guint ncols)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  if ((self->nnls == NULL) ||
      ((ncm_nnls_get_nrows (self->nnls) != nrows) || (ncm_nnls_get_ncols (self->nnls) != ncols)))
  {
    ncm_nnls_clear (&self->nnls);
    self->nnls = ncm_nnls_new (nrows, ncols);
    ncm_nnls_set_umethod (self->nnls, NCM_NNLS_UMETHOD_NORMAL);

    self->alloc_subs = FALSE;
  }

  if (!self->alloc_subs)
  {
    ncm_matrix_clear (&self->sub_IM);
    ncm_vector_clear (&self->sub_x);

    self->sub_IM = ncm_matrix_get_submatrix (self->IM, 0, 0, nrows, ncols);
    self->sub_x  = ncm_vector_get_subvector (self->weights, 0, ncols);

    self->alloc_subs = TRUE;
  }
}

static void
_ncm_stats_dist_prepare_interp (NcmStatsDist *sd, NcmVector *m2lnp)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  _ncm_stats_dist_prepare (sd);

  g_assert_cmpuint (ncm_vector_len (m2lnp), ==, self->n_obs);
  {
    NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);
    NcmStatsDistEval eval       = {sd, self, sd_class, NULL, m2lnp};
    const gdouble dbl_limit     = 2.0;
    guint i;

    /*
     * Evaluating the right-hand-side
     */
    self->min_m2lnp = GSL_POSINF;
    self->max_m2lnp = GSL_NEGINF;

    for (i = 0; i < self->n_kernels; i++)
    {
      const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);

      self->min_m2lnp = MIN (self->min_m2lnp, m2lnp_i);
      self->max_m2lnp = MAX (self->max_m2lnp, m2lnp_i);
    }

    if (self->max_m2lnp - self->min_m2lnp > -2.0 * dbl_limit * GSL_LOG_DBL_EPSILON)
    {
      guint n_cut = 0;

      /*
       * Everything in this block is indexed by kernel: the minimum and maximum above,
       * the loop below and the weights the indices end up addressing. Sorting all n_obs
       * entries would write n_obs indices into an array holding n_kernels of them, which
       * overflows whenever a cross-validation splits the sample.
       */
      g_array_set_size (self->m2lnp_sort, self->n_kernels);
      gsl_sort_index (&g_array_index (self->m2lnp_sort, size_t, 0),
                      ncm_vector_data (m2lnp),
                      ncm_vector_stride (m2lnp),
                      self->n_kernels);

      for (i = 0; i < self->n_kernels; i++)
      {
        guint p               = g_array_index (self->m2lnp_sort, size_t, i);
        const gdouble m2lnp_p = ncm_vector_get (m2lnp, p);

        if (m2lnp_p - self->min_m2lnp > -2.0 * dbl_limit * GSL_LOG_DBL_EPSILON)
        {
          n_cut = i;
          break;
        }
      }

      /* printf ("n_cut %d n_kernels %d\n", n_cut, self->n_kernels); */

      /* Too many points falling outside using normal KDE using 10%
       * of the weight for the points falling outside and 90% for
       * the points falling inside.
       */
      if (n_cut < (guint) (0.5 * self->n_obs))
      {
        ncm_vector_set_all (self->weights, 0.1 / (self->n_kernels - n_cut));

        for (i = 0; i < n_cut; i++)
        {
          guint p = g_array_index (self->m2lnp_sort, size_t, i);

          ncm_vector_set (self->weights, p, 0.9 / n_cut);
        }

        return;
      }

      {
        /*
         * n_cut counts kernels, this branch rebuilds the whole sample, so it needs its
         * own count over the observations. The two coincide only when no
         * cross-validation splits the sample.
         */
        GPtrArray *sample_array_cut = g_ptr_array_new ();
        NcmVector *m2lnp_cut;
        guint n_keep = 0;
        guint j      = 0;

        for (i = 0; i < self->n_obs; i++)
        {
          if (ncm_vector_get (m2lnp, i) - self->min_m2lnp <= -2.0 * dbl_limit * GSL_LOG_DBL_EPSILON)
            n_keep++;
        }

        m2lnp_cut = ncm_vector_new (n_keep);

        for (i = 0; i < self->n_obs; i++)
        {
          const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);

          if (m2lnp_i - self->min_m2lnp <= -2.0 * dbl_limit * GSL_LOG_DBL_EPSILON)
          {
            ncm_vector_set (m2lnp_cut, j++, m2lnp_i);
            g_ptr_array_add (sample_array_cut,
                             ncm_vector_ref (g_ptr_array_index (self->sample_array, i)));
          }
        }

        g_assert_cmpuint (j, ==, n_keep);

        g_ptr_array_set_size (self->sample_array, 0);

        for (i = 0; i < n_keep; i++)
        {
          g_ptr_array_add (self->sample_array, g_ptr_array_index (sample_array_cut, i));
        }

        g_ptr_array_unref (sample_array_cut);

        ncm_stats_dist_prepare_interp (sd, m2lnp_cut);

        ncm_vector_free (m2lnp_cut);
      }

      return;
    }

    /*
     * Preparing allocations
     */
    if ((self->n_obs != self->alloc_n_obs) || (self->n_kernels != self->alloc_n_kernels))
    {
      ncm_matrix_clear (&self->IM);
      ncm_vector_clear (&self->f);
      ncm_vector_clear (&self->f1);

      self->IM = ncm_matrix_new (self->n_obs, self->n_kernels);
      self->f  = ncm_vector_new (self->n_obs);
      self->f1 = ncm_vector_new (self->n_obs);

      ncm_vector_set_all (self->f1, 1.0);

      self->alloc_n_obs     = self->n_obs;
      self->alloc_n_kernels = self->n_kernels;
      self->alloc_subs      = FALSE;
    }

    ncm_vector_set_zero (self->weights);

    for (i = 0; i < self->n_obs; i++)
    {
      const gdouble m2lnp_i = ncm_vector_get (m2lnp, i);

      ncm_vector_set (self->f, i, exp (-0.5 * (m2lnp_i - self->min_m2lnp)));
    }

    if (self->n_kernels > 20000)
      g_warning ("_ncm_stats_dist_prepare_interp: very large system n = %u!", self->n_kernels);

    if (self->uniform_weights)
    {
      /* The bandwidth (and kernel) fit already ran in _ncm_stats_dist_prepare (); the
       * weights stay 1 / n_kernels and the NNLS fit is skipped. */
      ncm_vector_set_all (self->weights, 1.0 / (1.0 * self->n_kernels));
      self->rnorm = 0.0;

      return;
    }

    switch (self->cv_type)
    {
      case NCM_STATS_DIST_CV_SPLIT:
      {
        gdouble info[LM_INFO_SZ];
        gdouble opts[LM_OPTS_SZ];
        gdouble cov, ln_os, rnorm0;

        _ncm_stats_dist_alloc_nnls (sd, self->n_obs, self->n_kernels);

        if (self->levmar_n != self->n_obs)
        {
          g_clear_pointer (&self->levmar_workz, g_free);

          self->levmar_workz = g_new0 (gdouble, LM_DIF_WORKSZ (self->d, self->n_obs));
          self->levmar_n     = self->n_obs;
        }

        opts[0] = LM_INIT_MU;
        opts[1] = 1.0e-7;
        opts[2] = 1.0e-7;
        opts[3] = 1.0e-10;
        opts[4] = LM_DIFF_DELTA;

        ln_os = log (self->over_smooth);

        _ncm_stats_dist_compute_IM_full (sd);
        rnorm0 = NCM_NNLS_SOLVE (self->nnls, self->sub_IM, self->sub_x, self->f1);

        for (i = 0; i < 10; i++)
        {
          const gdouble ln_os_try = ncm_rng_gaussian_gen (self->rng, ln_os, 0.5);
          gdouble rnorm_try;

          self->over_smooth = exp (ln_os_try);
          _ncm_stats_dist_set_href (sd, _ncm_stats_dist_calc_href (sd));

          _ncm_stats_dist_compute_IM_full (sd);
          rnorm_try = NCM_NNLS_SOLVE (self->nnls, self->sub_IM, self->sub_x, self->f1);

          if (rnorm_try < rnorm0)
          {
            ln_os  = ln_os_try;
            rnorm0 = rnorm_try;
          }
        }

        dlevmar_dif (&_ncm_stats_dist_prepare_interp_fit_nnls_f,
                     &ln_os, NULL, 1, self->n_obs,
                     10000, opts, info, self->levmar_workz, &cov, &eval);

        self->over_smooth = exp (ln_os);
        _ncm_stats_dist_set_href (sd, _ncm_stats_dist_calc_href (sd));

        _ncm_stats_dist_compute_IM_full (sd);
        self->rnorm = NCM_NNLS_SOLVE (self->nnls, self->sub_IM, self->sub_x, self->f1);
      }
      break;
      case NCM_STATS_DIST_CV_SPLIT_NOFIT:
      case NCM_STATS_DIST_CV_SPLIT_ACCEPT:
      case NCM_STATS_DIST_CV_LOO:
      case NCM_STATS_DIST_CV_LOO_M2LNP:
      case NCM_STATS_DIST_CV_NONE:
        _ncm_stats_dist_alloc_nnls (sd, self->n_obs, self->n_kernels);
        _ncm_stats_dist_compute_IM_full (sd);
        self->rnorm = NCM_NNLS_SOLVE (self->nnls, self->sub_IM, self->sub_x, self->f1);
        break;
      default: /* LCOV_EXCL_BR_LINE */
        g_assert_not_reached ();
        break;
    }
  }

  {
    const gdouble total_weight = ncm_vector_sum_cpts (self->weights);

    if (!(gsl_finite (total_weight) && (total_weight > 0.0)))
    {
      /* A degenerate fit (all weights zero or not finite) falls back to the plain
       * kernel density estimate, which is still a valid proposal. */
      g_warning ("_ncm_stats_dist_prepare_interp: interpolation returned no positive weights (sum = %g), using uniform weights.", total_weight);
      ncm_vector_set_all (self->weights, 1.0 / (1.0 * self->n_kernels));
    }
    else
    {
      ncm_vector_scale (self->weights, 1.0 / total_weight);
    }
  }
}

static void
_ncm_stats_dist_reset (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  g_ptr_array_set_size (self->sample_array, 0);
}

/**
 * ncm_stats_dist_ref:
 * @sd: a #NcmStatsDist
 *
 * Increases the reference count of @sd.
 *
 * Returns: (transfer full): @sd.
 */
NcmStatsDist *
ncm_stats_dist_ref (NcmStatsDist *sd)
{
  return g_object_ref (sd);
}

/**
 * ncm_stats_dist_free:
 * @sd: a #NcmStatsDist
 *
 * Decreases the reference count of @sd.
 *
 */
void
ncm_stats_dist_free (NcmStatsDist *sd)
{
  g_object_unref (sd);
}

/**
 * ncm_stats_dist_clear:
 * @sd: a #NcmStatsDist
 *
 * Decreases the reference count of *@sd and sets the pointer *@sd to NULL.
 *
 */
void
ncm_stats_dist_clear (NcmStatsDist **sd)
{
  g_clear_object (sd);
}

/**
 * ncm_stats_dist_set_kernel:
 * @sd: a #NcmStatsDist
 * @sdk: a #NcmStatsDistKernel
 *
 * Sets the kernel to be used in the interpolation.
 * The different types of kernels are: the gaussian kernel and the student-t kernel,
 * which are under the file names ncm_stats_dist_kernel_gauss.c and ncm_stats_dist_kernel_st.c.
 */
void
ncm_stats_dist_set_kernel (NcmStatsDist *sd, NcmStatsDistKernel *sdk)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  ncm_stats_dist_kernel_clear (&self->kernel);
  self->kernel = ncm_stats_dist_kernel_ref (sdk);

  NCM_STATS_DIST_GET_CLASS (sd)->set_dim (sd, ncm_stats_dist_kernel_get_dim (sdk));
}

/**
 * ncm_stats_dist_peek_kernel:
 * @sd: a #NcmStatsDist
 *
 * Gets the kernel to be used in the interpolation.
 *
 * Returns: (transfer none): current #NcmStatsDistKernel used.
 */
NcmStatsDistKernel *
ncm_stats_dist_peek_kernel (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->kernel;
}

/**
 * ncm_stats_dist_get_kernel:
 * @sd: a #NcmStatsDist
 *
 * Gets the kernel to be used in the interpolation.
 *
 * Returns: (transfer full): current #NcmStatsDistKernel used.
 */
NcmStatsDistKernel *
ncm_stats_dist_get_kernel (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return ncm_stats_dist_kernel_ref (self->kernel);
}

/**
 * ncm_stats_dist_get_dim:
 * @sd: a #NcmStatsDist
 *
 * Returns: an int d, the dimension of the sample space, which is the same dimension of the used kernel.
 */
guint
ncm_stats_dist_get_dim (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->d;
}

/**
 * ncm_stats_dist_get_sample_size:
 * @sd: a #NcmStatsDist
 *
 * After the prepare call, this function returns the size of the sample used in the
 * interpolation.
 *
 * Returns: the size of the sample used.
 */
guint
ncm_stats_dist_get_sample_size (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->n_obs;
}

/**
 * ncm_stats_dist_get_n_kernels:
 * @sd: a #NcmStatsDist
 *
 * After the prepare call, this function returns the number of kernels used in the
 * interpolation.
 *
 * Returns: the number of kernels used.
 */
guint
ncm_stats_dist_get_n_kernels (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->n_kernels;
}

/**
 * ncm_stats_dist_get_href:
 * @sd: a #NcmStatsDist
 *
 * Gets the bandwidth applied to the kernels in the last preparation, the same one
 * used by ncm_stats_dist_get_lnnorm() and ncm_stats_dist_get_Ki(). It is the nominal
 * bandwidth the subclass derives from #NcmStatsDist:over-smooth, see the class
 * description, times the center shrinkage factor
 * ncm_stats_dist_get_center_shrink_factor(), which is one unless
 * #NcmStatsDist:center-shrink is set. It is zero before the first preparation.
 *
 * Returns: a double h, the currently used @href.
 */
gdouble
ncm_stats_dist_get_href (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->href;
}

/**
 * ncm_stats_dist_set_over_smooth:
 * @sd: a #NcmStatsDist
 * @over_smooth: the over-smooth factor
 *
 * Sets the over-smooth factor to @over_smooth.
 *
 */
void
ncm_stats_dist_set_over_smooth (NcmStatsDist *sd, const gdouble over_smooth)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->over_smooth = over_smooth;
}

/**
 * ncm_stats_dist_get_over_smooth:
 * @sd: a #NcmStatsDist
 *
 * Returns: a double os, the over-smooth factor.
 */
gdouble
ncm_stats_dist_get_over_smooth (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->over_smooth;
}

/**
 * ncm_stats_dist_set_split_frac:
 * @sd: a #NcmStatsDist
 * @split_frac: the over-smooth factor
 *
 * Sets cross-correlation split fraction to @split_frac.
 * This method shall be used when the cv_type is the cv_split.
 * The split fraction determines the fraction of sample points
 * that will be left out to use the cross validation method.
 *
 */
void
ncm_stats_dist_set_split_frac (NcmStatsDist *sd, const gdouble split_frac)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  g_assert_cmpfloat (split_frac, >=, 0.01);
  g_assert_cmpfloat (split_frac, <=, 1.0);

  self->split_frac = split_frac;
}

/**
 * ncm_stats_dist_get_split_frac:
 * @sd: a #NcmStatsDist
 *
 * Returns: a double @split_frac, the cross-correlation split fraction.
 */
gdouble
ncm_stats_dist_get_split_frac (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->split_frac;
}

/**
 * ncm_stats_dist_set_center_shrink:
 * @sd: a #NcmStatsDist
 * @center_shrink: whether to shrink the kernel centers toward the sample mean
 *
 * Enables or disables center shrinkage, see the class description. Takes effect
 * at the next call to ncm_stats_dist_prepare() or ncm_stats_dist_prepare_interp().
 *
 */
void
ncm_stats_dist_set_center_shrink (NcmStatsDist *sd, const gboolean center_shrink)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->center_shrink = center_shrink;
}

/**
 * ncm_stats_dist_get_center_shrink:
 * @sd: a #NcmStatsDist
 *
 * Returns: whether center shrinkage is enabled.
 */
gboolean
ncm_stats_dist_get_center_shrink (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->center_shrink;
}

/**
 * ncm_stats_dist_set_auto_kernel:
 * @sd: a #NcmStatsDist
 * @auto_kernel: whether to choose the kernel with the bandwidth
 *
 * Sets whether the kernel is chosen together with the over-smooth factor, see
 * #NcmStatsDist:auto-kernel.
 *
 */
void
ncm_stats_dist_set_auto_kernel (NcmStatsDist *sd, gboolean auto_kernel)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->auto_kernel = auto_kernel;
}

/**
 * ncm_stats_dist_get_auto_kernel:
 * @sd: a #NcmStatsDist
 *
 * Returns: whether the kernel is chosen together with the over-smooth factor.
 */
gboolean
ncm_stats_dist_get_auto_kernel (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->auto_kernel;
}

/**
 * ncm_stats_dist_get_center_shrink_factor:
 * @sd: a #NcmStatsDist
 *
 * Gets the scalar part $a = \det(A)^{1/d}$ of the center shrinkage transform $A$
 * applied in the last preparation, see the class description. The bandwidth in use is
 * $a$ times the nominal one. It is $1$ when center shrinkage is disabled, and for
 * #NcmStatsDistKDE with the sample covariance $A = a\,I$ exactly.
 *
 * Returns: the center shrinkage factor.
 */
gdouble
ncm_stats_dist_get_center_shrink_factor (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->center_a;
}

/**
 * ncm_stats_dist_peek_center_shrink_matrix:
 * @sd: a #NcmStatsDist
 *
 * Gets the $d \times d$ center shrinkage transform $A$ applied in the last
 * preparation: the centers are $c_i = \mu + A (x_i - \mu)$ and the kernel scale
 * matrices $\hat A \Sigma_i \hat A^T$ with $\hat A = A / a$, see the class
 * description. It is the identity when center shrinkage is disabled, and %NULL before
 * the first preparation.
 *
 * Returns: (transfer none) (nullable): the transform $A$.
 */
NcmMatrix *
ncm_stats_dist_peek_center_shrink_matrix (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->center_A;
}

/**
 * ncm_stats_dist_set_uniform_weights:
 * @sd: a #NcmStatsDist
 * @uniform_weights: whether to keep uniform weights in ncm_stats_dist_prepare_interp()
 *
 * Sets #NcmStatsDist:uniform-weights. Takes effect at the next preparation.
 *
 */
void
ncm_stats_dist_set_uniform_weights (NcmStatsDist *sd, const gboolean uniform_weights)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->uniform_weights = uniform_weights;
}

/**
 * ncm_stats_dist_get_uniform_weights:
 * @sd: a #NcmStatsDist
 *
 * Returns: #NcmStatsDist:uniform-weights.
 */
gboolean
ncm_stats_dist_get_uniform_weights (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->uniform_weights;
}

/**
 * ncm_stats_dist_set_defensive_frac:
 * @sd: a #NcmStatsDist
 * @frac: weight of the wide component, in $[0, 1]$
 *
 * Sets #NcmStatsDist:defensive-frac. Takes effect at the next preparation.
 *
 */
void
ncm_stats_dist_set_defensive_frac (NcmStatsDist *sd, const gdouble frac)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  g_assert_cmpfloat (frac, >=, 0.0);
  g_assert_cmpfloat (frac, <=, 1.0);
  self->defensive_frac = frac;
}

/**
 * ncm_stats_dist_get_defensive_frac:
 * @sd: a #NcmStatsDist
 *
 * Returns: #NcmStatsDist:defensive-frac.
 */
gdouble
ncm_stats_dist_get_defensive_frac (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->defensive_frac;
}

/**
 * ncm_stats_dist_set_defensive_scale:
 * @sd: a #NcmStatsDist
 * @scale: covariance factor of the wide component
 *
 * Sets #NcmStatsDist:defensive-scale. Takes effect at the next preparation.
 *
 */
void
ncm_stats_dist_set_defensive_scale (NcmStatsDist *sd, const gdouble scale)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  g_assert_cmpfloat (scale, >, 0.0);
  self->defensive_scale = scale;
}

/**
 * ncm_stats_dist_get_defensive_scale:
 * @sd: a #NcmStatsDist
 *
 * Returns: #NcmStatsDist:defensive-scale.
 */
gdouble
ncm_stats_dist_get_defensive_scale (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->defensive_scale;
}

/**
 * ncm_stats_dist_set_defensive_nu:
 * @sd: a #NcmStatsDist
 * @nu: degrees of freedom of the wide component
 *
 * Sets #NcmStatsDist:defensive-nu. Takes effect at the next preparation.
 *
 */
void
ncm_stats_dist_set_defensive_nu (NcmStatsDist *sd, const gdouble nu)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  g_assert_cmpfloat (nu, >, 0.0);
  self->defensive_nu = nu;
}

/**
 * ncm_stats_dist_get_defensive_nu:
 * @sd: a #NcmStatsDist
 *
 * Returns: #NcmStatsDist:defensive-nu.
 */
gdouble
ncm_stats_dist_get_defensive_nu (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->defensive_nu;
}

/**
 * ncm_stats_dist_set_print_fit:
 * @sd: a #NcmStatsDist
 * @print_fit: a boolean
 *
 * Whether to print steps during the fitting process.
 *
 */
void
ncm_stats_dist_set_print_fit (NcmStatsDist *sd, const gboolean print_fit)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->print_fit = print_fit;
}

/**
 * ncm_stats_dist_get_print_fit:
 * @sd: a #NcmStatsDist
 *
 * Returns: Whether it is going to print steps during the fitting process.
 */
gboolean
ncm_stats_dist_get_print_fit (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->print_fit;
}

/**
 * ncm_stats_dist_set_cv_type:
 * @sd: a #NcmStatsDist
 * @cv_type: a #NcmStatsDistCV
 *
 * Sets the cross-validation method to @cv_type.
 * If the selected method is none, all the sample points
 * will be used to compute the interpolation. If the cv_type is the cv_split,
 * a split fraction of the points are randomly excluded and the interpolation
 * is computed to a best fit of the remaining sample points,
 * which leads to a more point independent interpolation.
 *
 */
void
ncm_stats_dist_set_cv_type (NcmStatsDist *sd, const NcmStatsDistCV cv_type)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->cv_type = cv_type;
}

/**
 * ncm_stats_dist_get_cv_type:
 * @sd: a #NcmStatsDist
 *
 * Returns: a string @cv_type, current cross-validation method used.
 */
NcmStatsDistCV
ncm_stats_dist_get_cv_type (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->cv_type;
}

/**
 * ncm_stats_dist_set_use_threads:
 * @sd: a #NcmStatsDist
 * @use_threads: whether to use threads
 *
 * Sets whether to use OpenMP threads during the computation.
 *
 */
void
ncm_stats_dist_set_use_threads (NcmStatsDist *sd, const gboolean use_threads)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->use_threads = use_threads;
}

/**
 * ncm_stats_dist_get_use_threads:
 * @sd: a #NcmStatsDist
 *
 * Returns: whether to use OpenMP threads during the computation.
 */
gboolean
ncm_stats_dist_get_use_threads (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->use_threads;
}

/**
 * ncm_stats_dist_prepare_kernel: (virtual prepare_kernel)
 * @sd: a #NcmStatsDist
 * @sample_array: (element-type NcmVector): an array of #NcmVector
 *
 * Prepares the object for computations of the individuals kernels
 * and is usually part of ncm_stats_dist_prepare() and is should not
 * be called directly.
 *
 * This virtual method does not have a default implementation and
 * must be defined by the descendants.
 *
 */
void
ncm_stats_dist_prepare_kernel (NcmStatsDist *sd, GPtrArray *sample_array)
{
  NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);

  sd_class->prepare_kernel (sd, sample_array);
}

/**
 * ncm_stats_dist_prepare: (virtual prepare)
 * @sd: a #NcmStatsDist
 *
 * Prepares the object for calculations. This function prepares
 * the weight matrix and sets all the weights to 1.0/sample size.
 * It also calls the kernel_prepare function, implemented by a child,
 * and calls the get_href function.
 */
void
ncm_stats_dist_prepare (NcmStatsDist *sd)
{
  NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);

  sd_class->prepare (sd);
}

/**
 * ncm_stats_dist_prepare_interp: (virtual prepare_interp)
 * @sd: a #NcmStatsDist
 * @m2lnp: a #NcmVector containing the distribution values that will be used to compute the interpolation function.
 *
 * Prepares the object for calculations. Using the distribution values
 * at the sample points. This function calls the prepare function and
 * prepares the needed objects to compute the least squares problem.
 * The interpolation matrix IM is prepared by a child object and called in this function.
 * Then, depending on the cross validation method, the function solves the least squares problem using the ncm_nnls object.
 */
void
ncm_stats_dist_prepare_interp (NcmStatsDist *sd, NcmVector *m2lnp)
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  /* The sample's -2ln(L), for the cross-validation objectives that need it
   * (NCM_STATS_DIST_CV_SPLIT_ACCEPT); the cut branch re-enters here with the cut vector. */
  ncm_vector_clear (&self->cv_m2lnL_sample);
  self->cv_m2lnL_sample = ncm_vector_ref (m2lnp);

  sd_class->prepare_interp (sd, m2lnp);

  ncm_vector_clear (&self->cv_m2lnL_sample);
}

/**
 * ncm_stats_dist_eval:
 * @sd: a #NcmStatsDist
 * @x: a #NcmVector
 *
 * Evaluate the distribution at $\vec{x}=$@x. The method ncm_stats_dist_eval_m2lnp()
 * can be used to avoid underflow.
 *
 * Returns: $P(\vec{x})$.
 */
gdouble
ncm_stats_dist_eval (NcmStatsDist *sd, NcmVector *x)
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  const gdouble p                  = sd_class->eval_weights (sd, self->weights, x);

  if (self->defensive_frac <= 0.0)
    return p;

  return (1.0 - self->defensive_frac) * p + self->defensive_frac * exp (-0.5 * _ncm_stats_dist_defensive_m2lnK (sd, x));
}

/**
 * ncm_stats_dist_eval_m2lnp:
 * @sd: a #NcmStatsDist
 * @x: a #NcmVector
 *
 * Evaluate the distribution at $\vec{x}=$@x. This method is more
 * stable than ncm_stats_dist_eval() since it avoids underflows
 * and overflows.
 *
 * Returns: $P(\vec{x})$.
 */
gdouble
ncm_stats_dist_eval_m2lnp (NcmStatsDist *sd, NcmVector *x)
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  const gdouble m2lnp              = sd_class->eval_weights_m2lnp (sd, self->weights, x);

  if (self->defensive_frac <= 0.0)
    return m2lnp;

  {
    /* -2 ln [(1 - eps) p + eps K], summed in the log to keep the far tail. */
    const gdouble a = log1p (-self->defensive_frac) - 0.5 * m2lnp;
    const gdouble b = log (self->defensive_frac) - 0.5 * _ncm_stats_dist_defensive_m2lnK (sd, x);
    const gdouble m = GSL_MAX (a, b);

    return -2.0 * (m + log (exp (a - m) + exp (b - m)));
  }
}

/**
 * ncm_stats_dist_kernel_choose:
 * @sd: a #NcmStatsDist
 * @rng: a #NcmRNG
 *
 * Using the pseudo-random number generator @rng chooses
 * a random kernel based on the computed weights.
 *
 */
guint
ncm_stats_dist_kernel_choose (NcmStatsDist *sd, NcmRNG *rng)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  guint i;

  if (!self->wcum_ready)
  {
    gdouble cum = 0.0;

    ncm_vector_set (self->wcum, 0, cum);

    for (i = 0; i < self->n_kernels; i++)
    {
      cum += ncm_vector_get (self->weights, i);
      ncm_vector_set (self->wcum, i + 1, cum);
    }

    ncm_vector_scale (self->wcum, 1.0 / cum);
    self->wcum_ready = TRUE;
  }

  {
    const gdouble p = ncm_rng_uniform_gen (rng, 0.0, 1.0);
    gint ilo        = 0;
    gint ihi        = self->n_kernels;

    while (ihi > ilo + 1)
    {
      gint mi = (ihi + ilo) / 2;

      if (ncm_vector_fast_get (self->wcum, mi) > p)
        ihi = mi;
      else
        ilo = mi;
    }

    i = ilo;
  }

  return i;
}

/**
 * ncm_stats_dist_sample:
 * @sd: a #NcmStatsDist
 * @x: a #NcmVector
 * @rng: a #NcmRNG
 *
 * Using the pseudo-random number generator @rng generates a
 * point from the distribution and copy it to @x.
 *
 */
void
ncm_stats_dist_sample (NcmStatsDist *sd, NcmVector *x, NcmRNG *rng)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  if ((self->defensive_frac > 0.0) && (ncm_rng_uniform_gen (rng, 0.0, 1.0) < self->defensive_frac))
  {
    ncm_stats_dist_kernel_sample (self->defensive_kernel, self->defensive_decomp, 1.0, self->center_mean, x, rng);

    return;
  }

  {
    const gint i     = ncm_stats_dist_kernel_choose (sd, rng);
    NcmVector *x_i   = g_ptr_array_index (self->center_array, i);
    NcmMatrix *cov_U = ncm_stats_dist_peek_cov_decomp (sd, i);

    ncm_stats_dist_kernel_sample (self->kernel, cov_U, self->href, x_i, x, rng);
  }
}

/* This is an internal function used to sample antithetic variates
 * to improve the variance of the Monte Carlo integration.
 */
void
ncm_stats_dist_sample2 (NcmStatsDist *sd, NcmVector *x1, NcmVector *x2, NcmRNG *rng)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  const gint i                     = ncm_stats_dist_kernel_choose (sd, rng);
  const gint o_i                   = g_array_index (self->m2lnp_sort, size_t, i);
  NcmVector *x_i                   = g_ptr_array_index (self->center_array, o_i);
  NcmMatrix *cov_U_i               = ncm_stats_dist_peek_cov_decomp (sd, o_i);

  /* Each point takes the wide component independently with probability epsilon. */
  if ((self->defensive_frac > 0.0) && (ncm_rng_uniform_gen (rng, 0.0, 1.0) < self->defensive_frac))
    ncm_stats_dist_kernel_sample (self->defensive_kernel, self->defensive_decomp, 1.0, self->center_mean, x1, rng);
  else
    ncm_stats_dist_kernel_sample (self->kernel, cov_U_i, self->href, x_i, x1, rng);

  {
    const gint j       = self->sample_array->len - 1 - i;
    const gint o_j     = g_array_index (self->m2lnp_sort, size_t, j);
    NcmVector *x_j     = g_ptr_array_index (self->center_array, o_j);
    NcmMatrix *cov_U_j = ncm_stats_dist_peek_cov_decomp (sd, o_j);

    if ((self->defensive_frac > 0.0) && (ncm_rng_uniform_gen (rng, 0.0, 1.0) < self->defensive_frac))
      ncm_stats_dist_kernel_sample (self->defensive_kernel, self->defensive_decomp, 1.0, self->center_mean, x2, rng);
    else
      ncm_stats_dist_kernel_sample (self->kernel, cov_U_j, self->href, x_j, x2, rng);
  }
}

/**
 * ncm_stats_dist_get_rnorm:
 * @sd: a #NcmStatsDist
 *
 * Gets the value of the last $\chi^2$ fit obtained
 * when computing the interpolation through
 * ncm_stats_dist_prepare_interp().
 *
 * Returns: a double, the value of the $\chi^2$.
 */
gdouble
ncm_stats_dist_get_rnorm (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->rnorm * self->rnorm;
}

/**
 * ncm_stats_dist_add_obs:
 * @sd: a #NcmStatsDist
 * @y: a #NcmVector
 *
 * Adds a new point @y to the sample with weight 1.0.
 * This function must be called to insert an initial sample into the object, so the interpolation can be computed.
 *
 */
void
ncm_stats_dist_add_obs (NcmStatsDist *sd, NcmVector *x)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  g_ptr_array_add (self->sample_array, ncm_vector_dup (x));
}

/**
 * ncm_stats_dist_peek_sample_array:
 * @sd: a #NcmStatsDist
 *
 * Returns: (transfer none) (element-type NcmVector): current sample array.
 */
GPtrArray *
ncm_stats_dist_peek_sample_array (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->sample_array;
}

/**
 * ncm_stats_dist_peek_center_array:
 * @sd: a #NcmStatsDist
 *
 * Gets the kernel centers $c_i$ used in the last preparation, one per kernel. They
 * coincide with the first #NcmStatsDist:N sample points unless center shrinkage is
 * enabled, see the class description.
 *
 * Returns: (transfer none) (element-type NcmVector): current center array.
 */
GPtrArray *
ncm_stats_dist_peek_center_array (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->center_array;
}

/**
 * ncm_stats_dist_peek_cov_decomp: (virtual peek_cov_decomp)
 * @sd: a #NcmStatsDist
 * @i: kernel index
 *
 * Gets the covariance matrix associated with the @i-th kernel.
 *
 * Returns: (transfer none): Cholesky decomposition of the @i-th covariance matrix.
 */
NcmMatrix *
ncm_stats_dist_peek_cov_decomp (NcmStatsDist *sd, guint i)
{
  NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);

  return sd_class->peek_cov_decomp (sd, i);
}

/**
 * ncm_stats_dist_peek_full_cov_decomp: (virtual peek_full_cov_decomp)
 * @sd: a #NcmStatsDist
 *
 * Gets the full covariance matrix decomposition. This is a the Cholesky decomposition
 * of the covariance matrix of the whole sample.
 *
 * Returns: (transfer none): full covariance matrix decomposition.
 */
NcmMatrix *
ncm_stats_dist_peek_full_cov_decomp (NcmStatsDist *sd)
{
  NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);

  return sd_class->peek_full_cov_decomp (sd);
}

/**
 * ncm_stats_dist_peek_full_cov: (virtual peek_full_cov)
 * @sd: a #NcmStatsDist
 *
 * Gets the full covariance matrix of the whole sample.
 *
 * Returns: (transfer none): full covariance matrix.
 */
NcmMatrix *
ncm_stats_dist_peek_full_cov (NcmStatsDist *sd)
{
  NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);

  return sd_class->peek_full_cov (sd);
}

/**
 * ncm_stats_dist_get_lnnorm: (virtual get_lnnorm)
 * @sd: a #NcmStatsDist
 * @i: kernel index
 *
 * Gets the logarithm of the @i-th kernel normalization.
 *
 * Returns: $\ln (N_i)$.
 */
gdouble
ncm_stats_dist_get_lnnorm (NcmStatsDist *sd, guint i)
{
  NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);

  return sd_class->get_lnnorm (sd, i);
}

/**
 * ncm_stats_dist_peek_weights:
 * @sd: a #NcmStatsDist
 *
 * Returns: (transfer none): current kernel weights vector.
 */
NcmVector *
ncm_stats_dist_peek_weights (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->weights;
}

/**
 * ncm_stats_dist_reset: (virtual reset)
 * @sd: a #NcmStatsDist
 *
 * Reset the object discarding all added points.
 *
 */
void
ncm_stats_dist_reset (NcmStatsDist *sd)
{
  NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);

  sd_class->reset (sd);
}

/**
 * ncm_stats_dist_get_Ki:
 * @sd: a #NcmStatsDist
 * @i: kernel index
 * @y_i: (out callee-allocates) (transfer full): kernel location
 * @cov_i: (out callee-allocates) (transfer full): kernel covariance U
 * @n_i: (out): kernel normalization
 * @w_i: (out): kernel weight
 *
 * Return all information about the @i-th kernel.
 *
 */
void
ncm_stats_dist_get_Ki (NcmStatsDist *sd, const guint i, NcmVector **y_i, NcmMatrix **cov_i, gdouble *n_i, gdouble *w_i)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  NcmMatrix *cov_decomp            = ncm_stats_dist_peek_cov_decomp (sd, i);
  const gdouble lnnorm             = ncm_stats_dist_get_lnnorm (sd, i);
  const gdouble href               = self->href;

  /* Kernels, not observations: the two differ when a split cross-validation is used. */
  g_assert (i < ncm_stats_dist_get_n_kernels (sd));

  y_i[0]   = ncm_vector_dup (g_ptr_array_index (self->center_array, i));
  cov_i[0] = ncm_matrix_dup (cov_decomp);
  n_i[0]   = exp (lnnorm);
  w_i[0]   = ncm_vector_get (self->weights, i);

  ncm_matrix_triang_to_sym (cov_decomp, 'U', TRUE, cov_i[0]);

  ncm_matrix_scale (cov_i[0], href * href);
}

