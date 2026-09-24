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
 * Add the sample points with ncm_stats_dist_add_obs(), then call ncm_stats_dist_prepare():
 * without the sample's $-2\ln g(x_i)$ it weights every kernel equally, with them it solves
 * a non-negative least-squares problem for the weights. Either way the object is then
 * ready for ncm_stats_dist_eval() and ncm_stats_dist_sample().
 *
 * Nothing else has to be set. #NcmStatsDist:over-smooth sets the bandwidth,
 * #NcmStatsDist:CV-type the cross-validation that can fit it and
 * #NcmStatsDist:split-frac the out-of-sample fraction; #NcmStatsDist:center-shrink
 * contracts the mixture so its covariance matches the sample, and
 * #NcmStatsDist:auto-kernel fits the tail of the #NcmStatsDistKernelST the object was
 * built with along with the bandwidth, which requires a cross-validation that fits it.
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
#include <gsl/gsl_sort.h>
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

/* A d-vector of scratch per thread, for the evaluation paths. */
static gpointer
_ncm_stats_dist_dx_new (gpointer userdata)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (NCM_STATS_DIST (userdata));

  return ncm_vector_new (self->d);
}

static void
ncm_stats_dist_init (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->kernel              = NULL;
  self->sample_array        = g_ptr_array_new ();
  self->weights             = NULL;
  self->wcum                = NULL;
  self->wcum_ready          = FALSE;
  self->print_fit           = FALSE;
  self->over_smooth         = 0.0;
  self->cv_type             = NCM_STATS_DIST_CV_LEN;
  self->use_threads         = FALSE;
  self->split_frac          = 0.0;
  self->m2lnL_min           = 0.0;
  self->href                = 0.0;
  self->rnorm               = 0.0;
  self->auto_kernel         = FALSE;
  self->center_array        = g_ptr_array_new ();
  self->sample_mean         = NULL;
  self->sample_decomp       = NULL;
  self->kernel_cov          = NULL;
  self->shrink.on           = FALSE;
  self->shrink.A            = NULL;
  self->shrink.Ahat         = NULL;
  self->shrink.scale        = 1.0;
  self->shrink.is_isotropic = TRUE;
  self->shrink.eigval       = NULL;
  self->shrink.diag         = NULL;
  self->shrink.W            = NULL;
  self->shrink.Winv         = NULL;
  self->shrink.tmp          = NULL;
  self->shrink.ws           = ncm_lapack_ws_new ();
  self->refactor_M          = NULL;
  self->refactor_B          = NULL;
  self->defensive_frac      = 0.0;
  self->defensive_scale     = 4.0;
  self->defensive_nu        = 3.0;
  self->defensive_kernel    = NULL;
  self->defensive_decomp    = NULL;
  self->defensive_lnnorm    = 0.0;
  self->mp_dx               = ncm_memory_pool_new (&_ncm_stats_dist_dx_new, sd, (GDestroyNotify) ncm_vector_free);
  self->n_obs               = 0;
  self->n_kernels           = 0;
  self->d                   = 0;
  self->nnls                = NULL;
  self->IM                  = NULL;
  self->target              = NULL;
  self->ones                = NULL;
  self->cv_m2lnp            = NULL;
  self->cv_x                = g_ptr_array_new ();
  self->m2lnL               = NULL;
  self->cv_w                = NULL;
  self->uniform_weights     = FALSE;
  self->fit_weights         = FALSE;
  self->cut_weights         = FALSE;
  self->kernel_order        = g_array_new (FALSE, FALSE, sizeof (size_t));
  self->kernel_density      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  self->amise_x1            = NULL;
  self->amise_x2            = NULL;
  self->amise_stats         = ncm_stats_vec_new (2, NCM_STATS_VEC_COV, FALSE);
  self->amise_rng           = ncm_rng_seeded_new (NULL, 0);

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

static void _ncm_stats_dist_shrink_clear (NcmStatsDistShrink *shrink);

static void
_ncm_stats_dist_dispose (GObject *object)
{
  NcmStatsDist *sd                 = NCM_STATS_DIST (object);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  ncm_stats_dist_kernel_clear (&self->kernel);

  g_clear_pointer (&self->sample_array, g_ptr_array_unref); /* LCOV_EXCL_BR_LINE */
  g_clear_pointer (&self->center_array, g_ptr_array_unref); /* LCOV_EXCL_BR_LINE */
  ncm_vector_clear (&self->sample_mean);
  ncm_matrix_clear (&self->sample_decomp);
  ncm_matrix_clear (&self->kernel_cov);
  _ncm_stats_dist_shrink_clear (&self->shrink);
  ncm_matrix_clear (&self->refactor_M);
  ncm_matrix_clear (&self->refactor_B);

  if (self->mp_dx != NULL)
  {
    ncm_memory_pool_free (self->mp_dx, TRUE);
    self->mp_dx = NULL;
  }

  ncm_stats_dist_kernel_clear (&self->defensive_kernel);
  ncm_matrix_clear (&self->defensive_decomp);
  ncm_vector_clear (&self->weights);
  ncm_vector_clear (&self->wcum);

  ncm_nnls_clear (&self->nnls);

  ncm_matrix_clear (&self->IM);
  ncm_vector_clear (&self->target);
  ncm_vector_clear (&self->ones);
  ncm_vector_clear (&self->cv_m2lnp);
  ncm_vector_clear (&self->m2lnL);
  ncm_vector_clear (&self->cv_w);
  g_clear_pointer (&self->cv_x, g_ptr_array_unref);

  g_clear_pointer (&self->kernel_order, g_array_unref); /* LCOV_EXCL_BR_LINE */
  g_clear_pointer (&self->kernel_density, g_array_unref);
  ncm_vector_clear (&self->amise_x1);
  ncm_vector_clear (&self->amise_x2);
  ncm_stats_vec_clear (&self->amise_stats);
  ncm_rng_clear (&self->amise_rng);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_parent_class)->dispose (object);
}

static void
_ncm_stats_dist_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_stats_dist_parent_class)->finalize (object);
}

static void _ncm_stats_dist_set_dim (NcmStatsDist *sd, const guint dim);

/* LCOV_EXCL_START these should be overwritten and never executed */

static void
_ncm_stats_dist_prepare_shapes (NcmStatsDist *sd, GPtrArray *sample_array)
{
  g_error ("method prepare_shapes not implemented by %s.", G_OBJECT_TYPE_NAME (sd));
}

static void
_ncm_stats_dist_prepare_kernels (NcmStatsDist *sd)
{
  g_error ("method prepare_kernels not implemented by %s.", G_OBJECT_TYPE_NAME (sd));
}

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

static gdouble
_ncm_stats_dist_bandwidth (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  return self->over_smooth * ncm_stats_dist_kernel_get_rot_bandwidth (self->kernel, self->n_kernels);
}

static void _ncm_stats_dist_reset (NcmStatsDist *sd);

/* LCOV_EXCL_STOP */

static void
_ncm_stats_dist_eval_weights_m2lnp_vec (NcmStatsDist *sd, NcmVector *weights, GPtrArray *x_a, NcmVector *m2lnp)
{
  NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);
  guint i;

  for (i = 0; i < x_a->len; i++)
    ncm_vector_set (m2lnp, i, sd_class->eval_weights_m2lnp (sd, weights, g_ptr_array_index (x_a, i)));
}

/*
 * Leave-one-out companion of the batched evaluator: point @i of @x_a is the centre of
 * kernel @i, and the density there is wanted without that kernel. The default drops the
 * weight of one kernel at a time and evaluates point by point; a subclass that sweeps the
 * kernels in a batch overrides this and masks the one index instead, which is the same
 * arithmetic over a sweep it does once for every point.
 */
static void
_ncm_stats_dist_eval_weights_m2lnp_loo (NcmStatsDist *sd, NcmVector *weights, GPtrArray *x_a, NcmVector *m2lnp)
{
  NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);
  NcmVector *loo_weights      = ncm_vector_dup (weights);
  guint i;

  g_assert_cmpuint (ncm_vector_len (m2lnp), >=, x_a->len);

  for (i = 0; i < x_a->len; i++)
  {
    const gdouble w_i = ncm_vector_get (weights, i);

    ncm_vector_set (loo_weights, i, 0.0);
    ncm_vector_set (m2lnp, i, sd_class->eval_weights_m2lnp (sd, loo_weights, g_ptr_array_index (x_a, i)));
    ncm_vector_set (loo_weights, i, w_i);
  }

  ncm_vector_free (loo_weights);
}

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
   * out-of-sample objective. Requires a cross-validation that fits the bandwidth, currently
   * #NCM_STATS_DIST_CV_SPLIT_M2LNP; it is ignored otherwise. The kernel is a
   * #NcmStatsDistKernelST the object was built with -- any other kernel is refused at
   * prepare -- whose degrees of freedom $\nu$ are fitted in place, jointly with the
   * over-smooth factor, over $\nu \in [\nu_\mathrm{min}, 10^4]$, with
   * $\nu_\mathrm{min} = 2.5$ when #NcmStatsDist:center-shrink is set and $1$
   * otherwise; at the upper bound the kernel is Gaussian to $10^{-4}$. See the class
   * description. Default: FALSE.
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
   * Whether ncm_stats_dist_prepare() keeps uniform kernel weights instead of
   * fitting them by NNLS. The sample's $-2\ln L$ is then used only by the bandwidth
   * objectives that need it and by the outlier cut. Default: FALSE.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_UNIFORM_WEIGHTS,
                                   g_param_spec_boolean ("uniform-weights",
                                                         NULL,
                                                         "Keep uniform kernel weights instead of the NNLS fit when the sample -2lnL is given",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));


  klass->set_dim                = &_ncm_stats_dist_set_dim;
  klass->bandwidth              = &_ncm_stats_dist_bandwidth;
  klass->prepare_shapes         = &_ncm_stats_dist_prepare_shapes;
  klass->prepare_kernels        = &_ncm_stats_dist_prepare_kernels;
  klass->compute_IM             = &_ncm_stats_dist_compute_IM;
  klass->amise                  = NULL;
  klass->peek_cov_decomp        = &_ncm_stats_dist_peek_cov_decomp;
  klass->peek_full_cov_decomp   = &_ncm_stats_dist_peek_full_cov_decomp;
  klass->peek_full_cov          = &_ncm_stats_dist_peek_full_cov;
  klass->get_lnnorm             = &_ncm_stats_dist_get_lnnorm;
  klass->eval_weights           = &_ncm_stats_dist_eval_weights;
  klass->eval_weights_m2lnp     = &_ncm_stats_dist_eval_weights_m2lnp;
  klass->eval_weights_m2lnp_vec = &_ncm_stats_dist_eval_weights_m2lnp_vec;
  klass->eval_weights_m2lnp_loo = &_ncm_stats_dist_eval_weights_m2lnp_loo;
  klass->reset                  = &_ncm_stats_dist_reset;
}

static void
_ncm_stats_dist_shrink_clear (NcmStatsDistShrink *shrink)
{
  ncm_matrix_clear (&shrink->A);
  ncm_matrix_clear (&shrink->Ahat);
  ncm_vector_clear (&shrink->eigval);
  ncm_vector_clear (&shrink->diag);
  ncm_matrix_clear (&shrink->W);
  ncm_matrix_clear (&shrink->Winv);
  ncm_matrix_clear (&shrink->tmp);
  ncm_lapack_ws_clear (&shrink->ws);
}

static void
_ncm_stats_dist_shrink_alloc (NcmStatsDistShrink *shrink, const guint dim)
{
  NcmLapackWS *ws = shrink->ws;

  shrink->ws = NULL;
  _ncm_stats_dist_shrink_clear (shrink);
  shrink->ws = ws;

  shrink->A      = ncm_matrix_new (dim, dim);
  shrink->Ahat   = ncm_matrix_new (dim, dim);
  shrink->eigval = ncm_vector_new (dim);
  shrink->diag   = ncm_vector_new (dim);
  shrink->W      = ncm_matrix_new (dim, dim);
  shrink->Winv   = ncm_matrix_new (dim, dim);
  shrink->tmp    = ncm_matrix_new (dim, dim);

  ncm_matrix_set_identity (shrink->A);
  ncm_matrix_set_identity (shrink->Ahat);
  shrink->scale        = 1.0;
  shrink->is_isotropic = TRUE;
}

static void
_ncm_stats_dist_set_dim (NcmStatsDist *sd, const guint dim)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->d = dim;
  g_ptr_array_set_size (self->center_array, 0);

  ncm_vector_clear (&self->sample_mean);
  ncm_matrix_clear (&self->sample_decomp);
  ncm_matrix_clear (&self->kernel_cov);

  self->sample_mean   = ncm_vector_new (dim);
  self->sample_decomp = ncm_matrix_new (dim, dim);
  self->kernel_cov    = ncm_matrix_new (dim, dim);

  _ncm_stats_dist_shrink_alloc (&self->shrink, dim);

  ncm_matrix_clear (&self->refactor_M);
  ncm_matrix_clear (&self->refactor_B);
  self->refactor_M = ncm_matrix_new (dim, dim);
  self->refactor_B = ncm_matrix_new (dim, dim);

  ncm_matrix_clear (&self->defensive_decomp);
  self->defensive_decomp = ncm_matrix_new (dim, dim);
  ncm_memory_pool_empty (self->mp_dx, TRUE);

  ncm_vector_clear (&self->amise_x1);
  ncm_vector_clear (&self->amise_x2);
  self->amise_x1 = ncm_vector_new (dim);
  self->amise_x2 = ncm_vector_new (dim);
}

static void _ncm_stats_dist_split (NcmStatsDist *sd);
static void _ncm_stats_dist_sample_mean (NcmStatsDist *sd);
static void _ncm_stats_dist_update_defensive (NcmStatsDist *sd);
static void _ncm_stats_dist_shrink_prepare_shapes (NcmStatsDistShrink *shrink, NcmMatrix *UC, NcmMatrix *kernel_cov);
static void _ncm_stats_dist_fit_bandwidth (NcmStatsDist *sd);
static void _ncm_stats_dist_fit_weights (NcmStatsDist *sd);
static void _ncm_stats_dist_cut_weights (NcmStatsDist *sd);

static void
_ncm_stats_dist_prepare (NcmStatsDist *sd)
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  _ncm_stats_dist_split (sd);

  sd_class->prepare_shapes (sd, self->sample_array);
  _ncm_stats_dist_sample_mean (sd);
  _ncm_stats_dist_update_defensive (sd);

  if (self->shrink.on)
    _ncm_stats_dist_shrink_prepare_shapes (&self->shrink, self->sample_decomp, self->kernel_cov);

  ncm_vector_set_all (self->weights, 1.0 / (1.0 * self->n_kernels));
  self->wcum_ready = FALSE;

  _ncm_stats_dist_fit_bandwidth (sd);

  if (self->fit_weights)
    _ncm_stats_dist_fit_weights (sd);
  else if (self->cut_weights)
    _ncm_stats_dist_cut_weights (sd);
}

/*
 * Observations whose -2lnL sits more than this above the sample minimum have a density
 * below DBL_EPSILON^2 relative to the best point: nothing the weight fit can represent, so
 * the split step drops them before anything is built.
 */
#define NCM_STATS_DIST_M2LNL_RANGE (-4.0 * GSL_LOG_DBL_EPSILON)

static guint _ncm_stats_dist_n_kernels (NcmStatsDist *sd, const guint n);

/*
 * Step 1 of a prepare: trim, split, guard, and (re)make everything sized by n_obs or
 * n_kernels -- this is the only place that does. Reads the sample's -2lnL from
 * m2lnL when the caller gave one.
 */
static void
_ncm_stats_dist_split (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  guint i;

  self->fit_weights = FALSE;
  self->cut_weights = FALSE;
  self->rnorm       = 0.0;

  if (self->m2lnL != NULL)
  {
    NcmVector *m2lnL   = self->m2lnL;
    const guint n      = self->sample_array->len;
    const guint n_kern = _ncm_stats_dist_n_kernels (sd, n);
    guint n_keep       = 0;
    guint n_cut        = 0;

    g_assert_cmpuint (ncm_vector_len (m2lnL), ==, n);

    self->m2lnL_min = ncm_vector_get_min (m2lnL);

    /* n_keep counts observations, n_cut the kernel centres among them (the first n_kern). */
    for (i = 0; i < n; i++)
    {
      if (ncm_vector_get (m2lnL, i) - self->m2lnL_min <= NCM_STATS_DIST_M2LNL_RANGE)
      {
        n_keep++;

        if (i < n_kern)
          n_cut++;
      }
    }

    if (n_cut < 0.5 * n)
    {
      /*
       * Fewer than half the kernels are within NCM_STATS_DIST_M2LNL_RANGE of the best
       * point. Cutting the rest would leave too little to build from, so every point stays
       * a kernel and the weights are set to 0.9/0.1 instead (_ncm_stats_dist_cut_weights),
       * which favours the good points while keeping support where the others are. No
       * weight fit: its rows would be divided by an underflowing density.
       */
      self->cut_weights = TRUE;
    }
    else
    {
      /*
       * At least half the kernels are within range. The observations beyond it leave the
       * sample -- as kernels and as cross-validation points -- and everything is built
       * from what remains.
       */
      if (n_keep < n)
      {
        GPtrArray *kept       = g_ptr_array_new ();
        NcmVector *m2lnL_keep = ncm_vector_new (n_keep);
        guint j               = 0;

        for (i = 0; i < n; i++)
        {
          const gdouble m2lnL_i = ncm_vector_get (m2lnL, i);

          if (m2lnL_i - self->m2lnL_min <= NCM_STATS_DIST_M2LNL_RANGE)
          {
            ncm_vector_set (m2lnL_keep, j++, m2lnL_i);
            g_ptr_array_add (kept, ncm_vector_ref (g_ptr_array_index (self->sample_array, i)));
          }
        }

        g_ptr_array_set_size (self->sample_array, 0);

        for (i = 0; i < n_keep; i++)
          g_ptr_array_add (self->sample_array, g_ptr_array_index (kept, i));

        g_ptr_array_unref (kept);
        ncm_vector_clear (&self->m2lnL);
        self->m2lnL = m2lnL_keep;
      }

      self->fit_weights = !self->uniform_weights;
    }
  }

  self->n_obs     = self->sample_array->len;
  self->n_kernels = _ncm_stats_dist_n_kernels (sd, self->n_obs);

  if ((self->n_kernels >= self->n_obs) &&
      ((self->cv_type == NCM_STATS_DIST_CV_SPLIT_M2LNP) || (self->cv_type == NCM_STATS_DIST_CV_SPLIT_ACCEPT)))
    g_error ("_ncm_stats_dist_prepare: a sample split needs out-of-sample points, "
             "split-frac %g leaves none of %u.", self->split_frac, self->n_obs);

  if (self->n_obs <= self->d)
    g_error ("_ncm_stats_dist_prepare: the sample is too small.");

  if ((self->cv_type == NCM_STATS_DIST_CV_SPLIT_ACCEPT) && (self->m2lnL == NULL))
    g_error ("_ncm_stats_dist_prepare: NCM_STATS_DIST_CV_SPLIT_ACCEPT needs the sample's -2ln(L).");

  if ((self->weights == NULL) || (ncm_vector_len (self->weights) != self->n_kernels))
  {
    ncm_vector_clear (&self->weights);
    ncm_vector_clear (&self->wcum);
    ncm_vector_clear (&self->cv_w);

    self->weights = ncm_vector_new (self->n_kernels);
    self->wcum    = ncm_vector_new (self->n_kernels + 1);
    self->cv_w    = ncm_vector_new (self->n_kernels);
  }

  if ((self->cv_m2lnp == NULL) || (ncm_vector_len (self->cv_m2lnp) != self->n_obs))
  {
    ncm_vector_clear (&self->cv_m2lnp);
    self->cv_m2lnp = ncm_vector_new (self->n_obs);
  }

  {
    const guint cur_len = self->center_array->len;

    g_ptr_array_set_size (self->center_array, self->n_kernels);

    for (i = cur_len; i < self->n_kernels; i++)
      g_ptr_array_index (self->center_array, i) = ncm_vector_new (self->d);
  }

  /* The interpolation matrix only for what reads it: the AMISE objective and the weight fit. */
  if ((self->cv_type == NCM_STATS_DIST_CV_LOO) || self->fit_weights)
  {
    if ((self->IM == NULL) || (ncm_matrix_nrows (self->IM) != self->n_obs) || (ncm_matrix_ncols (self->IM) != self->n_kernels))
    {
      ncm_matrix_clear (&self->IM);
      ncm_vector_clear (&self->target);
      ncm_vector_clear (&self->ones);

      self->IM     = ncm_matrix_new (self->n_obs, self->n_kernels);
      self->target = ncm_vector_new (self->n_obs);
      self->ones   = ncm_vector_new (self->n_obs);

      ncm_vector_set_all (self->ones, 1.0);
    }
  }

  if (self->fit_weights)
  {
    if ((self->nnls == NULL) || (ncm_nnls_get_nrows (self->nnls) != self->n_obs) || (ncm_nnls_get_ncols (self->nnls) != self->n_kernels))
    {
      ncm_nnls_clear (&self->nnls);
      self->nnls = ncm_nnls_new (self->n_obs, self->n_kernels);
      ncm_nnls_set_umethod (self->nnls, NCM_NNLS_UMETHOD_NORMAL);
    }
  }
}

/* How many of the first n observations are kernel centres under the current cross-validation. */
static guint
_ncm_stats_dist_n_kernels (NcmStatsDist *sd, const guint n)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  switch (self->cv_type)
  {
    case NCM_STATS_DIST_CV_NONE:
    case NCM_STATS_DIST_CV_LOO:
    case NCM_STATS_DIST_CV_LOO_M2LNP:

      return n;

    case NCM_STATS_DIST_CV_SPLIT_M2LNP:
    case NCM_STATS_DIST_CV_SPLIT_ACCEPT:

      return ceil (n * self->split_frac);

    default: /* LCOV_EXCL_BR_LINE */
      g_assert_not_reached ();

      return 0;
  }
}

/* The mean of the kernel centres' sample points, m in c_i = m + A (x_i - m). */
static void
_ncm_stats_dist_sample_mean (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  guint i;

  ncm_vector_set_zero (self->sample_mean);

  for (i = 0; i < self->n_kernels; i++)
    ncm_vector_add (self->sample_mean, g_ptr_array_index (self->sample_array, i));

  ncm_vector_scale (self->sample_mean, 1.0 / (1.0 * self->n_kernels));
}

/*
 * The wide component: a Student-t centered on the sample mean with c times the sample
 * covariance as scale matrix, weight epsilon in the proposal. Built once per prepare, after
 * the shapes: it depends only on the sample mean and covariance factor.
 */
static void
_ncm_stats_dist_update_defensive (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  if (self->defensive_frac <= 0.0)
    return;

  if ((self->defensive_kernel == NULL) || (ncm_stats_dist_kernel_get_dim (self->defensive_kernel) != self->d))
  {
    ncm_stats_dist_kernel_clear (&self->defensive_kernel);
    self->defensive_kernel = NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_st_new (self->d, self->defensive_nu));
  }
  else
  {
    ncm_stats_dist_kernel_st_set_nu (NCM_STATS_DIST_KERNEL_ST (self->defensive_kernel), self->defensive_nu);
  }

  ncm_matrix_memcpy (self->defensive_decomp, self->sample_decomp);
  ncm_matrix_zero_triangle (self->defensive_decomp, 'U');
  ncm_matrix_scale (self->defensive_decomp, sqrt (self->defensive_scale));
  self->defensive_lnnorm = ncm_stats_dist_kernel_get_lnnorm (self->defensive_kernel, self->defensive_decomp);
}

/*
 * The part of the shrinkage that does not move during a fit. Writing C = L L^T and
 * G = L^{-1} <Sigma> L^{-T} = V diag(s) V^T, the transform that matches the mixture
 * covariance to C is
 *
 *   A = L V diag[(1 + kappa h^2 s)^{-1/2}] V^T L^{-1},
 *
 * so V and s, and with them L V and its inverse, depend only on C and <Sigma> -- both
 * fixed by prepare_shapes(). Only the d diagonal entries follow the bandwidth and the
 * kernel, which _ncm_stats_dist_set_href() then updates.
 *
 * Two properties come with this form. <Sigma> is positive semi-definite, so every
 * s >= 0 and 1 + kappa h^2 s >= 1: the diagonal can never be singular and there is
 * nothing to rescue. And A comes out independent of which factor of C is used -- send
 * L -> L Q for orthogonal Q and the Qs cancel -- where the Cholesky-only construction
 * A = U_C^T U_M^{-T} does not: Cholesky does not commute with permutation, so that one
 * changed with the order of the parameters.
 */
static void
_ncm_stats_dist_shrink_prepare_shapes (NcmStatsDistShrink *shrink, NcmMatrix *UC, NcmMatrix *kernel_cov)
{
  const guint d = ncm_matrix_nrows (UC); /* C = UC^T UC, so L = UC^T */
  NcmMatrix *G  = shrink->tmp;           /* scratch until A is built */
  gint ret;

  /* G = L^{-1} <Sigma> L^{-T} = UC^{-T} <Sigma> UC^{-1}, both solves reading only UC's
   * upper triangle. */
  ncm_matrix_memcpy (G, kernel_cov);
  ncm_matrix_dtrsm (G, 'L', 'U', 'T', 1.0, UC);
  ncm_matrix_dtrsm (G, 'R', 'U', 'N', 1.0, UC);

  /* G is symmetric, so the column-major view LAPACK takes of it is the same matrix. It
   * returns the eigenvectors in that view's columns, which are this matrix's rows: G
   * holds V^T afterwards. */
  ret = ncm_lapack_dsyevd ('V', 'U', d, ncm_matrix_data (G), ncm_matrix_tda (G),
                           ncm_vector_data (shrink->eigval), shrink->ws);

  if (ret != 0)
    g_error ("_ncm_stats_dist_shrink_prepare_shapes: eigendecomposition of the whitened mean kernel "
             "covariance failed with code %d.", ret);

  /* W = L V = UC^T V */
  ncm_matrix_transpose_memcpy (shrink->W, G);
  ncm_matrix_dtrmm (shrink->W, 'L', 'U', 'T', 1.0, UC);

  /* Winv = (L V)^{-1} = V^T L^{-1} = V^T UC^{-T} */
  ncm_matrix_memcpy (shrink->Winv, G);
  ncm_matrix_dtrsm (shrink->Winv, 'R', 'U', 'T', 1.0, UC);
}

static void _ncm_stats_dist_set_href (NcmStatsDist *sd, const gdouble href);
static gdouble _ncm_stats_dist_calc_href (NcmStatsDist *sd);

/* A bandwidth objective: scores the mixture as currently built (over_smooth, nu applied). */
typedef gdouble (*NcmStatsDistObjective) (NcmStatsDist *sd);

static void _ncm_stats_dist_fit (NcmStatsDist *sd, NcmStatsDistObjective objective);
static void _ncm_stats_dist_fit_over_smooth (NcmStatsDist *sd, NcmStatsDistObjective objective);
static gdouble _ncm_stats_dist_m2lnp (NcmStatsDist *sd);
static gdouble _ncm_stats_dist_accept (NcmStatsDist *sd);
static gdouble _ncm_stats_dist_loo_m2lnp (NcmStatsDist *sd);

/*
 * Step 5 of a prepare: the bandwidth. Every objective opens by setting over-smooth and
 * calling _ncm_stats_dist_set_href() itself, so nothing is set before the fit starts:
 * doing so would apply centre shrinkage, and refactor every covariance, against a kernel
 * the fit is about to replace. Only NCM_STATS_DIST_CV_NONE sets it here.
 */
static void
_ncm_stats_dist_fit_bandwidth (NcmStatsDist *sd)
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  switch (self->cv_type)
  {
    case NCM_STATS_DIST_CV_NONE:
      _ncm_stats_dist_set_href (sd, _ncm_stats_dist_calc_href (sd));
      break;
    case NCM_STATS_DIST_CV_SPLIT_M2LNP:
      _ncm_stats_dist_fit (sd, &_ncm_stats_dist_m2lnp);
      break;
    case NCM_STATS_DIST_CV_SPLIT_ACCEPT:
      _ncm_stats_dist_fit (sd, &_ncm_stats_dist_accept);
      break;
    case NCM_STATS_DIST_CV_LOO_M2LNP:
      _ncm_stats_dist_fit (sd, &_ncm_stats_dist_loo_m2lnp);
      break;
    case NCM_STATS_DIST_CV_LOO:
      /* A subclass may know its AMISE in closed form; the Monte Carlo estimate otherwise. */
      _ncm_stats_dist_fit (sd, (sd_class->amise != NULL) ? sd_class->amise : &_ncm_stats_dist_amise);
      break;
    default: /* LCOV_EXCL_BR_LINE */
      g_assert_not_reached ();
      break;
  }
}

/*
 * Sets the bandwidth and recomputes the kernel centers. With center shrinkage the
 * centers are c_i = mu + A (x_i - mu) and the kernel scale matrices Ahat Sigma_i Ahat^T,
 * where A = a Ahat, det Ahat = 1, and A solves A (C + kappa h^2 <Sigma>) A^T = C for the
 * sample covariance C = U_C^T U_C and the mean kernel scale matrix <Sigma>, both handed
 * over by the subclass in prepare_shapes(). The bandwidth stored and applied is a h.
 * Without center shrinkage A is the identity. Every place that changes href must go
 * through here so that centers, factors and bandwidth stay consistent.
 */
static void _ncm_stats_dist_shrink_prepare_kernels (NcmStatsDistShrink *shrink, const gdouble kappa, const gdouble h);

static void
_ncm_stats_dist_set_href (NcmStatsDist *sd, const gdouble href)
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  gdouble kappa                    = 1.0;
  guint i;

  g_assert_cmpuint (self->sample_array->len, >=, self->n_kernels);

  if (self->shrink.on)
  {
    kappa = ncm_stats_dist_kernel_get_var_factor (self->kernel);

    if (!gsl_finite (kappa))
      g_error ("_ncm_stats_dist_set_href: center shrinkage requires a kernel with a finite "
               "covariance, but %s has none (a Student-t kernel needs nu > 2). Use a "
               "Gaussian kernel or a larger nu, or disable center-shrink.",
               G_OBJECT_TYPE_NAME (self->kernel));
  }

  _ncm_stats_dist_shrink_prepare_kernels (&self->shrink, kappa, href);
  self->href = self->shrink.scale * href;

  if (self->shrink.on)
  {
    NcmVector **dx_ptr = ncm_memory_pool_get (self->mp_dx);
    NcmVector *dx      = *dx_ptr;

    for (i = 0; i < self->n_kernels; i++)
    {
      NcmVector *x_i = g_ptr_array_index (self->sample_array, i);
      NcmVector *c_i = g_ptr_array_index (self->center_array, i);

      ncm_vector_memcpy (dx, x_i);
      ncm_vector_sub (dx, self->sample_mean);
      ncm_matrix_update_vector (self->shrink.A, 'N', 1.0, dx, 0.0, c_i);
      ncm_vector_add (c_i, self->sample_mean);
    }

    ncm_memory_pool_return (dx_ptr);
  }
  else
  {
    for (i = 0; i < self->n_kernels; i++)
      ncm_vector_memcpy (g_ptr_array_index (self->center_array, i), g_ptr_array_index (self->sample_array, i));
  }

  sd_class->prepare_kernels (sd);
}

/*
 * The kernel stage of the shrinkage: from the basis and the bandwidth, A, its scale
 * a = det(A)^{1/d} and Ahat = A / a. Off, or with a basis that leaves A the identity to
 * rounding, the transform is the identity and is_isotropic says so, which lets the
 * subclasses skip refactoring their kernels.
 */
static void
_ncm_stats_dist_shrink_prepare_kernels (NcmStatsDistShrink *shrink, const gdouble kappa, const gdouble h)
{
  const guint d = ncm_matrix_nrows (shrink->A);

  if (!shrink->on)
  {
    ncm_matrix_set_identity (shrink->A);
    ncm_matrix_set_identity (shrink->Ahat);
    shrink->scale        = 1.0;
    shrink->is_isotropic = TRUE;

    return;
  }

  {
    const gdouble lambda = kappa * h * h;
    gdouble lnden        = 0.0;
    guint p;

    for (p = 0; p < d; p++)
    {
      const gdouble den = 1.0 + lambda * ncm_vector_get (shrink->eigval, p);

      ncm_vector_set (shrink->diag, p, 1.0 / sqrt (den));
      lnden += log (den);
    }

    /* A = W diag Winv, and det A = prod diag, so a = det(A)^(1/d) in closed form. */
    ncm_matrix_memcpy (shrink->tmp, shrink->W);
    ncm_matrix_scale_cols (shrink->tmp, shrink->diag);
    ncm_matrix_dgemm (shrink->A, 'N', 'N', 1.0, shrink->tmp, shrink->Winv, 0.0);

    shrink->scale = exp (-0.5 * lnden / d);
    ncm_matrix_memcpy (shrink->Ahat, shrink->A);
    ncm_matrix_scale (shrink->Ahat, 1.0 / shrink->scale);
    shrink->is_isotropic = ncm_matrix_is_identity (shrink->Ahat, 1.0e-12);

    if (shrink->is_isotropic)
      ncm_matrix_set_identity (shrink->Ahat);
  }
}

/*
 * Nominal bandwidth, computed by the subclass from over_smooth and the kernel rule of
 * thumb. The bandwidth actually applied to the kernels is set by
 * _ncm_stats_dist_set_href() and read back with ncm_stats_dist_get_href().
 */
static gdouble
_ncm_stats_dist_calc_href (NcmStatsDist *sd)
{
  return NCM_STATS_DIST_GET_CLASS (sd)->bandwidth (sd);
}

static void _ncm_stats_dist_fit_kernel (NcmStatsDist *sd, NcmStatsDistObjective objective, NcmStatsDistKernelST *st);

static void
_ncm_stats_dist_fit (NcmStatsDist *sd, NcmStatsDistObjective objective)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  if (!self->auto_kernel)
  {
    _ncm_stats_dist_fit_over_smooth (sd, objective);

    return;
  }

  if (!NCM_IS_STATS_DIST_KERNEL_ST (self->kernel))
    g_error ("_ncm_stats_dist_fit: auto-kernel fits the degrees of freedom of a Student-t kernel "
             "in place and needs an NcmStatsDistKernelST; the object was built with %s.",
             G_OBJECT_TYPE_NAME (self->kernel));

  _ncm_stats_dist_fit_kernel (sd, objective, NCM_STATS_DIST_KERNEL_ST (self->kernel));
}

/*
 * Chooses the kernel together with the over-smooth factor, using the same out-of-sample
 * objective. The Student-t kernel is (1 + chi2 / nu)^(-(nu + d) / 2), which is the Cauchy
 * kernel at nu = 1 and tends to the Gaussian one as nu grows, so the kernel is not a
 * discrete choice but the single continuous parameter nu. What is fitted is therefore
 * (ln over_smooth, ln nu) jointly, by the same objective and on the same out-of-sample points.
 *
 * The range of nu is bounded on both sides. Center shrinkage needs a finite kernel
 * covariance, nu / (nu - 2), so nu > 2 there. At the other end the kernel is the Gaussian
 * one to better than 1e-4 by nu ~ 1e4, while the (1 + chi2 / nu) form starts losing
 * precision well above that, so 1e4 is the ceiling; a fit that reaches it leaves the
 * Student-t kernel there, Gaussian to that accuracy. The kernel the object was built with
 * is tuned in place: nothing is created or replaced during a prepare.
 */

#define NCM_STATS_DIST_AUTO_KERNEL_NU_MAX (1.0e4)

/* The over-smooth range every fit searches, in ln. */
#define NCM_STATS_DIST_LN_OS_MIN (log (1.0e-2))
#define NCM_STATS_DIST_LN_OS_MAX (log (2.0e1))

static gdouble _ncm_stats_dist_fit_run (NcmStatsDist *sd, NcmStatsDistObjective objective, const guint n,
                                        gdouble *p, const gdouble *lb, const gdouble *ub, const gdouble *step);

/*
 * The (ln over_smooth, ln nu) fit, at the center-shrinkage setting the caller has left in
 * place. Both start from the previous fit, as over_smooth always did: the evaluations go
 * into the walk from the start to the optimum, not into the refinement, and a warm start
 * cuts them from 49 to 14 per prepare at d = 50 with the same optimum. Leaves the object
 * holding the solution.
 */
static void
_ncm_stats_dist_fit_kernel (NcmStatsDist *sd, NcmStatsDistObjective objective, NcmStatsDistKernelST *st)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  const gdouble nu_min             = self->shrink.on ? 2.5 : 1.0;
  const gdouble nu_start           = CLAMP (ncm_stats_dist_kernel_st_get_nu (st), nu_min, NCM_STATS_DIST_AUTO_KERNEL_NU_MAX);
  gdouble p[2]                     = {CLAMP (log (self->over_smooth), NCM_STATS_DIST_LN_OS_MIN, NCM_STATS_DIST_LN_OS_MAX), log (nu_start)};
  const gdouble step[2]            = {0.1, 0.5};
  const gdouble lb[2]              = {NCM_STATS_DIST_LN_OS_MIN, log (nu_min)};
  const gdouble ub[2]              = {NCM_STATS_DIST_LN_OS_MAX, log (NCM_STATS_DIST_AUTO_KERNEL_NU_MAX)};

  ncm_stats_dist_kernel_st_set_nu (st, nu_start);
  _ncm_stats_dist_fit_run (sd, objective, 2, p, lb, ub, step);
}

typedef struct _NcmStatsDistFit
{
  NcmStatsDist *sd;
  NcmStatsDistObjective objective;
  guint neval;
} NcmStatsDistFit;

static gdouble _ncm_stats_dist_fit_f (guint n, const gdouble *x, gdouble *grad, gpointer data);

/*
 * Minimizes @objective over ln over_smooth (@n = 1) or over (ln over_smooth, ln nu)
 * (@n = 2, the kernel must be a Student-t) inside [@lb, @ub], starting at @p. BOBYQA
 * builds a quadratic model of a smooth objective and takes the bounds natively, which
 * keeps nu inside the range where the kernel is both defined and numerically sound
 * without the objective having to clamp its own argument; the simplex methods could
 * terminate on the nu lower bound instead of the interior minimum. Comparison against
 * Nelder-Mead, Subplex and COBYLA: dev-notes/apes_center_shrink/auto_tuning.md.
 *
 * Leaves the object holding the minimizer -- over_smooth, nu, and everything
 * _ncm_stats_dist_set_href() derives from them -- and @p updated to it. Returns the minimum.
 */
static gdouble
_ncm_stats_dist_fit_run (NcmStatsDist *sd, NcmStatsDistObjective objective, const guint n,
                         gdouble *p, const gdouble *lb, const gdouble *ub, const gdouble *step)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistFit fit              = {sd, objective, 0};
  nlopt_opt opt                    = nlopt_create (NLOPT_LN_BOBYQA, n);
  gdouble minf                     = 0.0;
  gint ret;

  nlopt_set_lower_bounds (opt, lb);
  nlopt_set_upper_bounds (opt, ub);
  nlopt_set_initial_step (opt, step);

  nlopt_set_xtol_rel (opt, 1.0e-3);
  nlopt_set_maxeval (opt, 1000);
  nlopt_set_min_objective (opt, &_ncm_stats_dist_fit_f, &fit);

  ret = nlopt_optimize (opt, p, &minf);

  if (ret < 0) /* LCOV_EXCL_BR_LINE */
    g_warning ("_ncm_stats_dist_fit_run: nlopt failed with code %d, keeping the configuration it reached.", ret);

  /* Re-evaluating at the minimizer leaves the object holding the chosen configuration. */
  _ncm_stats_dist_fit_f (n, p, NULL, &fit);

  if (self->print_fit)
  {
    if (n == 2)
      printf ("# center shrinkage: %s, over-smooth: % 22.15g, nu: % 22.15g, objective = % 22.15g, "
              "neval = %u, nlopt status (%d)\n",
              self->shrink.on ? "on" : "off", self->over_smooth, exp (p[1]), minf, fit.neval, ret);
    else
      printf ("# over-smooth: % 22.15g, objective = % 22.15g, neval = %u, nlopt status (%d)\n",
              self->over_smooth, minf, fit.neval, ret);
  }

  nlopt_destroy (opt);

  return minf;
}

/* The nlopt objective: applies the trial hyperparameters, then scores. */
static gdouble
_ncm_stats_dist_fit_f (guint n, const gdouble *x, gdouble *grad, gpointer data)
{
  NcmStatsDistFit *fit             = data;
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (fit->sd);

  g_assert (grad == NULL);

  self->over_smooth = exp (x[0]);

  if (n == 2)
    ncm_stats_dist_kernel_st_set_nu (NCM_STATS_DIST_KERNEL_ST (self->kernel), exp (x[1]));

  _ncm_stats_dist_set_href (fit->sd, _ncm_stats_dist_calc_href (fit->sd));
  fit->neval++;

  return fit->objective (fit->sd);
}

/* The bandwidth alone, the kernel as it is. */
static void
_ncm_stats_dist_fit_over_smooth (NcmStatsDist *sd, NcmStatsDistObjective objective)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  gdouble p[1]                     = {CLAMP (log (self->over_smooth), NCM_STATS_DIST_LN_OS_MIN, NCM_STATS_DIST_LN_OS_MAX)};
  const gdouble step[1]            = {0.1};
  const gdouble lb[1]              = {NCM_STATS_DIST_LN_OS_MIN};
  const gdouble ub[1]              = {NCM_STATS_DIST_LN_OS_MAX};

  _ncm_stats_dist_fit_run (sd, objective, 1, p, lb, ub, step);
}

static gdouble
_ncm_stats_dist_m2lnp (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  gdouble m2lnp                    = 0.0;
  gint i;

  /* Per-point values first, summed serially afterwards: a shared accumulator inside the
   * evaluation was a data race, and even a reduction would round in thread order. */
  /*
   * The out-of-sample points go in one batch rather than one at a time. They are measured
   * against the same kernels with the same weights, so a single sweep reads each kernel
   * covariance once for the whole batch instead of once per point, which is where this
   * objective spent its time.
   */
  g_ptr_array_set_size (self->cv_x, 0);

  for (i = self->n_kernels; i < (gint) self->n_obs; i++)
    g_ptr_array_add (self->cv_x, g_ptr_array_index (self->sample_array, i));

  {
    NcmVector *m2lnp_oos = ncm_vector_get_subvector (self->cv_m2lnp, self->n_kernels,
                                                     self->n_obs - self->n_kernels);

    ncm_stats_dist_eval_m2lnp_vec (sd, self->cv_x, m2lnp_oos);
    ncm_vector_free (m2lnp_oos);
  }

  for (i = self->n_kernels; i < (gint) self->n_obs; i++)
    m2lnp += ncm_vector_get (self->cv_m2lnp, i);

  if (self->print_fit)
    ncm_message ("# over-smooth: % 22.15g, m2lnp = % 22.15g\n",
                 self->over_smooth, m2lnp);

  return m2lnp;
}

static inline gdouble _ncm_stats_dist_logaddexp (const gdouble a, const gdouble b);

/*
 * Out-of-sample estimate of the independence-sampler acceptance (auto_tuning.md). With
 * r = ln q - ln pi at the kernel points (k) and at the out-of-sample points (j), the mean over j
 * of E_{x' ~ q}[min (1, exp (r_j - r_{x'}))] is estimated by self-normalized importance
 * sampling with the kernel points as draws from pi and weights proportional to exp (r_k).
 * Returns -ln of the estimate. Needs the sample's -2ln(L) (m2lnL).
 */
gdouble
_ncm_stats_dist_accept (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  const gint nk                    = self->n_kernels;
  const gint no                    = self->n_obs;
  gdouble lse_k                    = GSL_NEGINF;
  gdouble acc                      = 0.0;
  gint i, j;

  /* r_i = ln q (x_i) - ln pi (x_i) = (m2lnL_i - m2lnq_i) / 2 at every sample point. The
   * sample points are the whole batch, so this is one call: fanning out over the scalar
   * evaluator instead would re-sweep the kernels once per point, and the threading now
   * lives inside the batched evaluator. */
  ncm_stats_dist_eval_m2lnp_vec (sd, self->sample_array, self->cv_m2lnp);

  for (i = 0; i < no; i++)
    ncm_vector_set (self->cv_m2lnp, i, 0.5 * (ncm_vector_get (self->m2lnL, i) - ncm_vector_get (self->cv_m2lnp, i)));

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

static inline gdouble
_ncm_stats_dist_logaddexp (const gdouble a, const gdouble b)
{
  if (a == GSL_NEGINF)
    return b;

  if (b == GSL_NEGINF)
    return a;

  return GSL_MAX (a, b) + log1p (exp (-fabs (a - b)));
}

/*
 * Leave-one-out likelihood cross-validation, -2 sum_i ln q_{-i} (x_i): every point is a
 * kernel and q_{-i} is the mixture with weight i zeroed and the rest rescaled.
 */
gdouble
_ncm_stats_dist_loo_m2lnp (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  const gint n                     = self->n_kernels;
  gdouble m2lnp                    = 0.0;
  gint i;

  /*
   * The kernel centres go in one batch. Each point drops its own kernel, which is one
   * masked index of a sweep that is otherwise the same for all of them, so the kernels are
   * swept once here instead of once per point.
   */
  g_ptr_array_set_size (self->cv_x, 0);

  for (i = 0; i < n; i++)
    g_ptr_array_add (self->cv_x, g_ptr_array_index (self->sample_array, i));

  {
    NcmVector *m2lnp_loo = ncm_vector_get_subvector (self->cv_m2lnp, 0, n);

    sd_class->eval_weights_m2lnp_loo (sd, self->weights, self->cv_x, m2lnp_loo);
    ncm_vector_free (m2lnp_loo);
  }

  /* Renormalization of the weights that are left, one term per point. */
  for (i = 0; i < n; i++)
    m2lnp += ncm_vector_get (self->cv_m2lnp, i) + 2.0 * log1p (-ncm_vector_get (self->weights, i));

  if (self->print_fit)
    ncm_message ("# over-smooth: % 22.15g, loo m2lnp = % 22.15g\n", self->over_smooth, m2lnp);

  return m2lnp;
}

static void _ncm_stats_dist_sample2 (NcmStatsDist *sd, NcmVector *x1, NcmVector *x2, NcmRNG *rng);

/* The AMISE objective by Monte Carlo over antithetic pairs; shared with the subclasses. */
gdouble
_ncm_stats_dist_amise (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  gdouble amise                    = 0.0;
  guint i, j;

  sd_class->compute_IM (sd, self->IM);

  g_array_set_size (self->kernel_order, self->n_kernels);
  g_array_set_size (self->kernel_density, self->n_kernels);

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

    g_array_index (self->kernel_density, gdouble, i) = (row_sum + ncm_matrix_get (self->IM, i, i)) / self->n_kernels;
  }

  gsl_sort_index (&g_array_index (self->kernel_order, size_t, 0),
                  (gdouble *) self->kernel_density->data, 1, self->kernel_density->len);

  {
    NcmRNG *rng        = self->amise_rng;
    NcmVector *x1      = self->amise_x1;
    NcmVector *x2      = self->amise_x2;
    NcmStatsVec *stats = self->amise_stats;
    guint max_iter     = 100000000;
    gdouble mean;
    gdouble p1, p2;

    /* The same draws every evaluation, so the objective is a deterministic function of h. */
    ncm_rng_set_seed (rng, 0);
    ncm_stats_vec_reset (stats, TRUE);

    for (i = 0; i < 100; i++)
    {
      _ncm_stats_dist_sample2 (sd, x1, x2, rng);
      p1 = ncm_stats_dist_eval (sd, x1);
      p2 = ncm_stats_dist_eval (sd, x2);

      ncm_stats_vec_set (stats, 0, p1);
      ncm_stats_vec_set (stats, 1, p2);
      ncm_stats_vec_update (stats);
    }

    for (i = 0; i < max_iter; i++)
    {
      _ncm_stats_dist_sample2 (sd, x1, x2, rng);
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
  }

  if (self->print_fit)
    ncm_message ("# over-smooth: % 22.15g, amise = % 22.15g\n",
                 self->over_smooth, amise);

  return amise;
}

/* This is an internal function used to sample antithetic variates
 * to improve the variance of the Monte Carlo integration.
 */
static void
_ncm_stats_dist_sample2 (NcmStatsDist *sd, NcmVector *x1, NcmVector *x2, NcmRNG *rng)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  const gint i                     = ncm_stats_dist_kernel_choose (sd, rng);
  const gint o_i                   = g_array_index (self->kernel_order, size_t, i);
  NcmVector *x_i                   = g_ptr_array_index (self->center_array, o_i);
  NcmMatrix *cov_U_i               = ncm_stats_dist_peek_cov_decomp (sd, o_i);

  /* Each point takes the wide component independently with probability epsilon. */
  if ((self->defensive_frac > 0.0) && (ncm_rng_uniform_gen (rng, 0.0, 1.0) < self->defensive_frac))
    ncm_stats_dist_kernel_sample (self->defensive_kernel, self->defensive_decomp, 1.0, self->sample_mean, x1, rng);
  else
    ncm_stats_dist_kernel_sample (self->kernel, cov_U_i, self->href, x_i, x1, rng);

  {
    const gint j       = self->sample_array->len - 1 - i;
    const gint o_j     = g_array_index (self->kernel_order, size_t, j);
    NcmVector *x_j     = g_ptr_array_index (self->center_array, o_j);
    NcmMatrix *cov_U_j = ncm_stats_dist_peek_cov_decomp (sd, o_j);

    if ((self->defensive_frac > 0.0) && (ncm_rng_uniform_gen (rng, 0.0, 1.0) < self->defensive_frac))
      ncm_stats_dist_kernel_sample (self->defensive_kernel, self->defensive_decomp, 1.0, self->sample_mean, x2, rng);
    else
      ncm_stats_dist_kernel_sample (self->kernel, cov_U_j, self->href, x_j, x2, rng);
  }
}

static void _ncm_stats_dist_compute_IM_full (NcmStatsDist *sd);

/*
 * Step 6 of a prepare: the kernel weights by non-negative least squares, so that the
 * mixture interpolates the sample's density. Always last -- shrinkage is computed for
 * the uniform-weight mixture and this is a refinement on top of it.
 */
static void
_ncm_stats_dist_fit_weights (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  guint i;

  for (i = 0; i < self->n_obs; i++)
    ncm_vector_set (self->target, i, exp (-0.5 * (ncm_vector_get (self->m2lnL, i) - self->m2lnL_min)));

  if (self->n_kernels > 20000)
    g_warning ("_ncm_stats_dist_fit_weights: very large system n = %u!", self->n_kernels);

  ncm_vector_set_zero (self->weights);
  _ncm_stats_dist_compute_IM_full (sd);
  self->rnorm = NCM_NNLS_SOLVE (self->nnls, self->IM, self->weights, self->ones);

  {
    const gdouble total_weight = ncm_vector_sum_cpts (self->weights);

    if (!(gsl_finite (total_weight) && (total_weight > 0.0)))
    {
      /* A degenerate fit (all weights zero or not finite) falls back to the plain
       * kernel density estimate, which is still a valid proposal. */
      g_warning ("_ncm_stats_dist_fit_weights: the NNLS fit returned no positive weights (sum = %g), using uniform weights.", total_weight);
      ncm_vector_set_all (self->weights, 1.0 / (1.0 * self->n_kernels));
    }
    else
    {
      ncm_vector_scale (self->weights, 1.0 / total_weight);
    }
  }
}

static void
_ncm_stats_dist_compute_IM_full (NcmStatsDist *sd)
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  guint i;

  sd_class->compute_IM (sd, self->IM);

  /* #pragma omp parallel for if (self->use_threads) */

  for (i = 0; i < self->n_obs; i++)
    ncm_matrix_mul_row (self->IM, i, 1.0 / ncm_vector_get (self->target, i));
}

/*
 * The weights when more than half the kernels sit beyond NCM_STATS_DIST_M2LNL_RANGE of the
 * best point: 0.9 spread over the kernels within it, 0.1 over the rest, so the proposal
 * favours the good points while keeping support where the others are. Applies whether or
 * not the weight fit is enabled: that fit cannot run on rows whose density underflows.
 */
static void
_ncm_stats_dist_cut_weights (NcmStatsDist *sd)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  guint n_near                     = 0;
  guint i;

  for (i = 0; i < self->n_kernels; i++)
  {
    if (ncm_vector_get (self->m2lnL, i) - self->m2lnL_min <= NCM_STATS_DIST_M2LNL_RANGE)
      n_near++;
  }

  if ((n_near == 0) || (n_near == self->n_kernels))
    return;

  for (i = 0; i < self->n_kernels; i++)
  {
    const gboolean near = (ncm_vector_get (self->m2lnL, i) - self->m2lnL_min <= NCM_STATS_DIST_M2LNL_RANGE);

    ncm_vector_set (self->weights, i, near ? 0.9 / n_near : 0.1 / (self->n_kernels - n_near));
  }
}

/*
 * The family's one Cholesky: the upper factor of @cov into @decomp. A covariance that is
 * not positive definite to rounding is repaired with nearPD and reported, since it means
 * the sample is degenerate somewhere; one that nearPD cannot repair is fatal. @what names
 * the matrix in both messages.
 */
void
_ncm_stats_dist_cholesky (NcmMatrix *decomp, const NcmMatrix *cov, const guint maxiter, const gchar *what)
{
  gboolean repaired;
  const gint ret = ncm_matrix_cholesky_decomp_nearPD (cov, decomp, 'U', maxiter, &repaired);

  if (ret != 0)
    g_error ("_ncm_stats_dist_cholesky: %s is not positive definite and nearPD could not repair it in %u "
             "iterations. The sample is degenerate: a parameter that does not vary, or coincident points.",
             what, maxiter);

  if (repaired)
    g_warning ("_ncm_stats_dist_cholesky: %s was not positive definite to rounding; using the nearest "
               "positive definite matrix instead.", what);
}

/* Used by the subclasses; declared in ncm_stats_dist_private.h. */
#define NCM_STATS_DIST_CENTER_NEARPD_MAXITER (200)

void
_ncm_stats_dist_refactor_decomp (NcmStatsDist *sd, NcmMatrix *U0, NcmMatrix *U)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  if (self->shrink.is_isotropic)
  {
    ncm_matrix_memcpy (U, U0);

    return;
  }

  {
    NcmMatrix *M = self->refactor_M;
    NcmMatrix *B = self->refactor_B;

    /* M = U0 Ahat^T; only U0's upper triangle is read. */
    ncm_matrix_transpose_memcpy (M, self->shrink.Ahat);
    ncm_matrix_dtrmm (M, 'L', 'U', 'N', 1.0, U0);

    ncm_matrix_dgemm (B, 'T', 'N', 1.0, M, M, 0.0); /* B = M^T M = Ahat Sigma Ahat^T */
    _ncm_stats_dist_cholesky (U, B, NCM_STATS_DIST_CENTER_NEARPD_MAXITER, "the shrunk kernel covariance");
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
 * at the next call to ncm_stats_dist_prepare() or ncm_stats_dist_prepare().
 *
 */
void
ncm_stats_dist_set_center_shrink (NcmStatsDist *sd, const gboolean center_shrink)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  self->shrink.on = center_shrink;
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

  return self->shrink.on;
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

  return self->shrink.scale;
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

  return self->shrink.A;
}

/**
 * ncm_stats_dist_set_uniform_weights:
 * @sd: a #NcmStatsDist
 * @uniform_weights: whether to keep uniform weights in ncm_stats_dist_prepare()
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
 * ncm_stats_dist_prepare_shapes: (virtual prepare_shapes)
 * @sd: a #NcmStatsDist
 * @sample_array: (element-type NcmVector): an array of #NcmVector
 *
 * Runs the first stage of a prepare alone: the per-kernel covariance structures the
 * subclass builds from the sample, before any bandwidth is applied. For inspecting those
 * structures; ncm_stats_dist_prepare() runs it as part of the whole pipeline, and only
 * after that is the object ready for ncm_stats_dist_eval().
 *
 * This virtual method has no default implementation.
 */
void
ncm_stats_dist_prepare_shapes (NcmStatsDist *sd, GPtrArray *sample_array)
{
  NcmStatsDistClass *sd_class = NCM_STATS_DIST_GET_CLASS (sd);

  sd_class->prepare_shapes (sd, sample_array);
}

/**
 * ncm_stats_dist_prepare:
 * @sd: a #NcmStatsDist
 * @m2lnL: (nullable): the sample's $-2\ln L$, one entry per observation added with
 * ncm_stats_dist_add_obs(), or %NULL
 *
 * Builds the estimator from the sample: the kernel shapes, the bandwidth (fitted by the
 * cross-validation chosen with #NcmStatsDist:CV-type) and the weights. Leaves the object
 * ready for ncm_stats_dist_eval() and ncm_stats_dist_sample().
 *
 * With @m2lnL the weights are fitted by non-negative least squares so that the mixture
 * interpolates the sample's density (unless #NcmStatsDist:uniform-weights is set), and
 * the cross-validations that score against the target density become available
 * (#NCM_STATS_DIST_CV_SPLIT_ACCEPT). Observations whose $-2\ln L$ lies more than
 * $4 |\ln \epsilon|$ above the minimum have a density that underflows. While more than
 * half the kernels are such points every observation stays a kernel and the weights are
 * 0.9 spread over the kernels within that range and 0.1 over the rest; otherwise those
 * observations are dropped from the sample before anything is built. Without @m2lnL the
 * weights are uniform, $1 / n$.
 */
void
ncm_stats_dist_prepare (NcmStatsDist *sd, NcmVector *m2lnL)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

  ncm_vector_clear (&self->m2lnL);

  if (m2lnL != NULL)
    self->m2lnL = ncm_vector_ref (m2lnL);

  _ncm_stats_dist_prepare (sd);

  ncm_vector_clear (&self->m2lnL);
}

static gdouble _ncm_stats_dist_defensive_m2lnK (NcmStatsDist *sd, NcmVector *x);

static gdouble
_ncm_stats_dist_defensive_mix (NcmStatsDist *sd, NcmVector *x, const gdouble m2lnp)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);

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

/* -2 ln of the wide component's density at x. */
static gdouble
_ncm_stats_dist_defensive_m2lnK (NcmStatsDist *sd, NcmVector *x)
{
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  NcmVector **dx_ptr               = ncm_memory_pool_get (self->mp_dx);
  NcmVector *dx                    = *dx_ptr;
  gdouble chi2, m2lnK;

  ncm_vector_memcpy (dx, x);
  ncm_vector_sub (dx, self->sample_mean);
  ncm_matrix_dtrsv (self->defensive_decomp, 'U', 'T', dx);
  chi2  = ncm_vector_dot (dx, dx);
  m2lnK = 2.0 * self->defensive_lnnorm - 2.0 * log (ncm_stats_dist_kernel_eval_unnorm (self->defensive_kernel, chi2));
  ncm_memory_pool_return (dx_ptr);

  return m2lnK;
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

  return _ncm_stats_dist_defensive_mix (sd, x, sd_class->eval_weights_m2lnp (sd, self->weights, x));
}

/**
 * ncm_stats_dist_eval_m2lnp_vec:
 * @sd: a #NcmStatsDist
 * @x_a: (element-type NcmVector): an array of #NcmVector
 * @m2lnp: a #NcmVector holding the result
 *
 * Evaluates the distribution at every point of @x_a, writing $-2\ln P(\vec{x}_i)$ into
 * element $i$ of @m2lnp. Equivalent to calling ncm_stats_dist_eval_m2lnp() once per
 * point, but subclasses may share work across the batch, so the two need not agree to
 * the last bit.
 *
 */
void
ncm_stats_dist_eval_m2lnp_vec (NcmStatsDist *sd, GPtrArray *x_a, NcmVector *m2lnp)
{
  NcmStatsDistClass *sd_class      = NCM_STATS_DIST_GET_CLASS (sd);
  NcmStatsDistPrivate * const self = ncm_stats_dist_get_instance_private (sd);
  guint i;

  g_assert_cmpuint (ncm_vector_len (m2lnp), >=, x_a->len);

  sd_class->eval_weights_m2lnp_vec (sd, self->weights, x_a, m2lnp);

  if (self->defensive_frac <= 0.0)
    return;

  for (i = 0; i < x_a->len; i++)
  {
    NcmVector *x_i = g_ptr_array_index (x_a, i);

    ncm_vector_set (m2lnp, i, _ncm_stats_dist_defensive_mix (sd, x_i, ncm_vector_get (m2lnp, i)));
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
    ncm_stats_dist_kernel_sample (self->defensive_kernel, self->defensive_decomp, 1.0, self->sample_mean, x, rng);

    return;
  }

  {
    const gint i     = ncm_stats_dist_kernel_choose (sd, rng);
    NcmVector *x_i   = g_ptr_array_index (self->center_array, i);
    NcmMatrix *cov_U = ncm_stats_dist_peek_cov_decomp (sd, i);

    ncm_stats_dist_kernel_sample (self->kernel, cov_U, self->href, x_i, x, rng);
  }
}

/**
 * ncm_stats_dist_get_rnorm:
 * @sd: a #NcmStatsDist
 *
 * Gets the value of the last $\chi^2$ fit obtained
 * when computing the interpolation through
 * ncm_stats_dist_prepare().
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

