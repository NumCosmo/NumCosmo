/***************************************************************************
 *            ncm_fit_esmcmc_walker_apes.c
 *
 *  Sat October 27 13:08:13 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_apes.c
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
 * NcmFitESMCMCWalkerAPES:
 *
 * Ensemble sampler Markov Chain Monte Carlo walker - apes move.
 *
 * Implementing apes move walker for #NcmFitESMCMC.
 *
 * This object implements the Approximate Posterior Ensemble Sample (APES) step proposal
 * for a walker. This proposal was developed by Sandro Dias Pinto Vitenti and
 * implemented in this library. Below there is a description of the proposal.
 *
 * The APES proposal consists of using radial basis interpolation to generate an
 * interpolant $\tilde{\pi}$ from a target distribution $\pi$ and use this interpolant
 * to propose new points for the walker. By using a distribution $\tilde{\pi}$ that
 * resembles the original target distribution, the APES proposal generates samples that
 * converge faster to the target distribution and are more independent when compared to
 * other step proposals.
 *
 * The APES step is implemented as follows: suppose that there are $L$ walkers. They are
 * divided into two blocks $L_1$ and $L_2$, containing the first and the second half of
 * the walkers respectively. When proposing new points $Y$ for the walkers in the $L_1$
 * block, we use the points in the $L_2$ block to generate an interpolant
 * $\tilde{\pi}_{L_2}$ and then propose points $Y \sim \tilde{\pi}_{L_2}$ for the $L_1$
 * block. These points are accepted or rejected based on an acceptance probability
 * $A(Y|X)$, and after the points of the first block are updated, we do the same
 * procedure for the $L_2$ block using the $L_1$ block. This procedure can be seen in
 * the pseudocode below.
 *
 * ![apes_sketch](apes.png)
 *
 * The user must provide the input the values: @nwalkers, @nparams, @method, @k\_@type,
 * @over\_@smooth$ and @use\_@interp - ncm\_fit\_esmcmc\_walker\_apes\_new\_full(). The
 * user can also initialize the object with: @nwalkers, @nparams -
 * ncm\_fit\_esmcmc\_walker\_apes\_new() and let the remaining parameters as default,
 * which are defined in the properties of the class. For more information about the
 * algorithm, check the explanation below.
 *
 *    - This object shall be used in the #NcmFitESMCMC class to generate a Monte Carlo
 *      Markov Chain using an ensemble sampler. To see an example of its implementation,
 *      check the file example\_rosenbrock.py in NumCosmo/examples.
 *
 *    - Regarding the radial basis interpolation method is implemented, check the
 *      #NcmStatsDist class.
 *
 *    - Regarding the types of kernel used in the interpolation method as the radial
 *      basis function, check the #NcmStatsDistKernel class.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_esmcmc_walker.h"
#include "math/ncm_fit_esmcmc_walker_apes.h"

#include "math/ncm_c.h"
#include "math/ncm_fit_esmcmc.h"
#include "math/ncm_stats_dist_vkde.h"
#include "math/ncm_stats_dist_kernel_st.h"
#include "math/ncm_stats_dist_kernel_gauss.h"

#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */


enum
{
  PROP_0,
  PROP_METHOD,
  PROP_K_TYPE,
  PROP_OVER_SMOOTH,
  PROP_SHRINK,
  PROP_RANDOM_WALK_PROB,
  PROP_RANDOM_WALK_SCALE,
  PROP_USE_INTERP,
  PROP_USE_THREADS,
};

typedef struct _NcmFitESMCMCWalkerAPESRandomWalk
{
  NcmVector *std;
  NcmVector *lb;
  NcmVector *ub;
} NcmFitESMCMCWalkerAPESRandomWalk;

typedef struct _NcmFitESMCMCWalkerAPESPrivate
{
  guint size;
  guint size_2;
  guint nparams;
  guint a_size;
  guint a_nparams;
  guint mk;
  NcmVector *m2lnp_star;
  NcmVector *m2lnp_cur;
  gchar *desc;
  NcmStatsDist *sd0;
  NcmStatsDist *sd1;
  GPtrArray *thetastar;
  NcmVector *m2lnL_s0;
  NcmVector *m2lnL_s1;
  NcmFitESMCMCWalkerAPESMethod method;
  NcmFitESMCMCWalkerAPESKType k_type;
  gdouble over_smooth;
  gdouble shrink;            /* Shrink factor for weight computation */
  gdouble random_walk_prob;  /* Probability of random walk step */
  gdouble random_walk_scale; /* Scale of the random walk step */
  NcmFitESMCMCWalkerAPESRandomWalk rw0;
  NcmFitESMCMCWalkerAPESRandomWalk rw1;
  gboolean use_interp;
  gboolean use_threads;
  gboolean constructed;
  guint exploration;
} NcmFitESMCMCWalkerAPESPrivate;

struct _NcmFitESMCMCWalkerAPES
{
  NcmFitESMCMCWalker parent_instance;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmFitESMCMCWalkerAPES, ncm_fit_esmcmc_walker_apes, NCM_TYPE_FIT_ESMCMC_WALKER)

#define __MK(method, k_type) (method + (k_type << 8))

static void
ncm_fit_esmcmc_walker_apes_init (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->size              = 0;
  self->size_2            = 0;
  self->nparams           = 0;
  self->a_size            = 0;
  self->a_nparams         = 0;
  self->mk                = -1;
  self->m2lnp_star        = NULL;
  self->m2lnp_cur         = NULL;
  self->desc              = NULL;
  self->sd0               = NULL;
  self->sd1               = NULL;
  self->thetastar         = g_ptr_array_new ();
  self->m2lnL_s0          = NULL;
  self->m2lnL_s1          = NULL;
  self->method            = NCM_FIT_ESMCMC_WALKER_APES_METHOD_LEN;
  self->k_type            = NCM_FIT_ESMCMC_WALKER_APES_KTYPE_LEN;
  self->over_smooth       = 0.0;
  self->shrink            = 0.0;
  self->random_walk_prob  = 0.0;
  self->random_walk_scale = 0.0;
  self->use_interp        = FALSE;
  self->use_threads       = FALSE;
  self->constructed       = FALSE;
  self->exploration       = 0;

  self->rw0.std = NULL;
  self->rw0.lb  = NULL;
  self->rw0.ub  = NULL;
  self->rw1.std = NULL;
  self->rw1.lb  = NULL;
  self->rw1.ub  = NULL;

  g_ptr_array_set_free_func (self->thetastar, (GDestroyNotify) ncm_vector_free);
}

static void
_ncm_fit_esmcmc_walker_apes_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerAPES *apes = NCM_FIT_ESMCMC_WALKER_APES (object);

  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_APES (object));

  switch (prop_id)
  {
    case PROP_METHOD:
      ncm_fit_esmcmc_walker_apes_set_method (apes, g_value_get_enum (value));
      break;
    case PROP_K_TYPE:
      ncm_fit_esmcmc_walker_apes_set_k_type (apes, g_value_get_enum (value));
      break;
    case PROP_OVER_SMOOTH:
      ncm_fit_esmcmc_walker_apes_set_over_smooth (apes, g_value_get_double (value));
      break;
    case PROP_SHRINK:
      ncm_fit_esmcmc_walker_apes_set_shrink (apes, g_value_get_double (value));
      break;
    case PROP_RANDOM_WALK_PROB:
      ncm_fit_esmcmc_walker_apes_set_random_walk_prob (apes, g_value_get_double (value));
      break;
    case PROP_RANDOM_WALK_SCALE:
      ncm_fit_esmcmc_walker_apes_set_random_walk_scale (apes, g_value_get_double (value));
      break;
    case PROP_USE_INTERP:
      ncm_fit_esmcmc_walker_apes_use_interp (apes, g_value_get_boolean (value));
      break;
    case PROP_USE_THREADS:
      ncm_fit_esmcmc_walker_apes_set_use_threads (apes, g_value_get_boolean (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_fit_esmcmc_walker_apes_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerAPES *apes = NCM_FIT_ESMCMC_WALKER_APES (object);

  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_APES (object));

  switch (prop_id)
  {
    case PROP_METHOD:
      g_value_set_enum (value, ncm_fit_esmcmc_walker_apes_get_method (apes));
      break;
    case PROP_K_TYPE:
      g_value_set_enum (value, ncm_fit_esmcmc_walker_apes_get_k_type (apes));
      break;
    case PROP_OVER_SMOOTH:
      g_value_set_double (value, ncm_fit_esmcmc_walker_apes_get_over_smooth (apes));
      break;
    case PROP_SHRINK:
      g_value_set_double (value, ncm_fit_esmcmc_walker_apes_get_shrink (apes));
      break;
    case PROP_RANDOM_WALK_PROB:
      g_value_set_double (value, ncm_fit_esmcmc_walker_apes_get_random_walk_prob (apes));
      break;
    case PROP_RANDOM_WALK_SCALE:
      g_value_set_double (value, ncm_fit_esmcmc_walker_apes_get_random_walk_scale (apes));
      break;
    case PROP_USE_INTERP:
      g_value_set_boolean (value, ncm_fit_esmcmc_walker_apes_interp (apes));
      break;
    case PROP_USE_THREADS:
      g_value_set_boolean (value, ncm_fit_esmcmc_walker_apes_get_use_threads (apes));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void _ncm_fit_esmcmc_walker_apes_set_sys (NcmFitESMCMCWalker *walker);

static void
_ncm_fit_esmcmc_walker_apes_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_apes_parent_class)->constructed (object);
  {
    NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (object);
    NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

    self->constructed = TRUE;
    _ncm_fit_esmcmc_walker_apes_set_sys (NCM_FIT_ESMCMC_WALKER (object));
  }
}

static void
_ncm_fit_esmcmc_walker_apes_dispose (GObject *object)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (object);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  ncm_vector_clear (&self->m2lnp_star);
  ncm_vector_clear (&self->m2lnp_cur);
  ncm_vector_clear (&self->m2lnL_s0);
  ncm_vector_clear (&self->m2lnL_s1);

  ncm_stats_dist_clear (&self->sd0);
  ncm_stats_dist_clear (&self->sd1);

  ncm_vector_clear (&self->rw0.std);
  ncm_vector_clear (&self->rw0.lb);
  ncm_vector_clear (&self->rw0.ub);
  ncm_vector_clear (&self->rw1.std);
  ncm_vector_clear (&self->rw1.lb);
  ncm_vector_clear (&self->rw1.ub);

  g_clear_pointer (&self->thetastar, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_apes_parent_class)->dispose (object);
}

static void
_ncm_fit_esmcmc_walker_apes_finalize (GObject *object)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (object);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  g_clear_pointer (&self->desc, g_free);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_apes_parent_class)->finalize (object);
}

static void _ncm_fit_esmcmc_walker_apes_set_size (NcmFitESMCMCWalker *walker, guint size);
static guint _ncm_fit_esmcmc_walker_apes_get_size (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_apes_set_nparams (NcmFitESMCMCWalker *walker, guint nparams);
static guint _ncm_fit_esmcmc_walker_apes_get_nparams (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_apes_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng);
static void _ncm_fit_esmcmc_walker_apes_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static gdouble _ncm_fit_esmcmc_walker_apes_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star);
static gdouble _ncm_fit_esmcmc_walker_apes_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static void _ncm_fit_esmcmc_walker_apes_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf);
static const gchar *_ncm_fit_esmcmc_walker_apes_desc (NcmFitESMCMCWalker *walker);

static void
ncm_fit_esmcmc_walker_apes_class_init (NcmFitESMCMCWalkerAPESClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcmFitESMCMCWalkerClass *walker_class = NCM_FIT_ESMCMC_WALKER_CLASS (klass);

  object_class->set_property = &_ncm_fit_esmcmc_walker_apes_set_property;
  object_class->get_property = &_ncm_fit_esmcmc_walker_apes_get_property;
  object_class->constructed  = &_ncm_fit_esmcmc_walker_apes_constructed;
  object_class->dispose      = &_ncm_fit_esmcmc_walker_apes_dispose;
  object_class->finalize     = &_ncm_fit_esmcmc_walker_apes_finalize;

  /**
   * NcmFitESMCMCWalkerAPES:method:
   *
   * Method used in posterior approximation.
   * This property can be set to one of the #NcmFitESMCMCWalkerAPESMethod values.
   * The default value is #NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_METHOD,
                                   g_param_spec_enum ("method",
                                                      NULL,
                                                      "Method used in posterior approximation",
                                                      NCM_TYPE_FIT_ESMCMC_WALKER_APES_METHOD, NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:k-type:
   *
   * Kernel used in posterior approximation. This property can be set to one of the
   * #NcmFitESMCMCWalkerAPESKType values. The default value is
   * #NCM_FIT_ESMCMC_WALKER_APES_KTYPE_CAUCHY.
   */
  g_object_class_install_property (object_class,
                                   PROP_K_TYPE,
                                   g_param_spec_enum ("kernel-type",
                                                      NULL,
                                                      "Kernel used in posterior approximation",
                                                      NCM_TYPE_FIT_ESMCMC_WALKER_APES_KTYPE, NCM_FIT_ESMCMC_WALKER_APES_KTYPE_CAUCHY,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:over-smooth:
   *
   * Over-smooth parameter used to adjust kernel bandwidth. The default value is 1.0.
   */
  g_object_class_install_property (object_class,
                                   PROP_OVER_SMOOTH,
                                   g_param_spec_double ("over-smooth",
                                                        NULL,
                                                        "Over-smooth parameter used to adjust kernel bandwidth",
                                                        1.0e-10, 1.0e10, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:shrink:
   *
   * Shrink factor for weight computation. This property defines the shrink factor
   * used in the weight computation for the APES proposal. The default value is 0.01.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SHRINK,
                                   g_param_spec_double ("shrink",
                                                        NULL,
                                                        "Shrink factor for weight computation",
                                                        0.0, 1.0, 0.01,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:random-walk-prob:
   *
   * Probability of random walk step. This property defines the probability of a random
   * walk step being taken when proposing new points for the walkers. The default value
   * is 0.02, meaning a random walk step will be taken with probability of 2%.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RANDOM_WALK_PROB,
                                   g_param_spec_double ("random-walk-prob",
                                                        NULL,
                                                        "Probability of random walk step",
                                                        0.0001, 1.0, 0.02,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:random-walk-scale:
   *
   * Scale factor for the random walk step used in proposal generation. This property
   * defines the standard deviation of the random walk proposal as a fraction of the
   * empirical standard deviation computed from the current half-ensemble (i.e., the
   * half not being updated). The default value is 0.25, meaning the random walk step
   * will have a standard deviation equal to 25% of that empirical value.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_RANDOM_WALK_SCALE,
                                   g_param_spec_double ("random-walk-scale",
                                                        NULL,
                                                        "Scale of the random walk step",
                                                        0.01, 1.0, 0.25,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:use-interp:
   *
   * Whether to use interpolation to build the posterior approximation. This property
   * defines whether the walker will use interpolation to build the posterior
   * approximation. The default value is TRUE, meaning interpolation will be used. If
   * set to FALSE, the walker will not use interpolation.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_INTERP,
                                   g_param_spec_boolean ("use-interp",
                                                         NULL,
                                                         "Whether to use interpolation to build the posterior approximation",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:use-threads:
   *
   * Whether to use threads when building the posterior approximation. This property
   * defines whether the walker will use threads when building the posterior
   * approximation. The default value is FALSE, meaning threads will not be used. If set
   * to TRUE, the walker will use threads.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_USE_THREADS,
                                   g_param_spec_boolean ("use-threads",
                                                         NULL,
                                                         "Whether to use threads when building the posterior approximation",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  walker_class->set_size    = &_ncm_fit_esmcmc_walker_apes_set_size;
  walker_class->get_size    = &_ncm_fit_esmcmc_walker_apes_get_size;
  walker_class->set_nparams = &_ncm_fit_esmcmc_walker_apes_set_nparams;
  walker_class->get_nparams = &_ncm_fit_esmcmc_walker_apes_get_nparams;
  walker_class->setup       = &_ncm_fit_esmcmc_walker_apes_setup;
  walker_class->step        = &_ncm_fit_esmcmc_walker_apes_step;
  walker_class->prob        = &_ncm_fit_esmcmc_walker_apes_prob;
  walker_class->prob_norm   = &_ncm_fit_esmcmc_walker_apes_prob_norm;
  walker_class->clean       = &_ncm_fit_esmcmc_walker_apes_clean;
  walker_class->desc        = &_ncm_fit_esmcmc_walker_apes_desc;
}

static void
_ncm_fit_esmcmc_walker_apes_vkde_check_sizes (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  guint cov_estimates0 = ncm_stats_dist_vkde_get_local_frac (NCM_STATS_DIST_VKDE (self->sd0)) * self->size_2;
  guint cov_estimates1 = ncm_stats_dist_vkde_get_local_frac (NCM_STATS_DIST_VKDE (self->sd1)) * self->size_2;

  if (cov_estimates0 < 2)
    g_error ("Number of walkers per block (%d) is too low for the current dimension (%d).\n"
             "\tToo few points (%d) to estimate local covariances.",
             self->size_2, self->nparams, cov_estimates0);

  if (cov_estimates1 < 2)
    g_error ("Number of walkers per block (%d) is too low for the current dimension (%d).\n"
             "\tToo few points (%d) to estimate local covariances.",
             self->size_2, self->nparams, cov_estimates1);
}

static void
_ncm_fit_esmcmc_walker_apes_set_sys (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if ((self->size != self->a_size) ||
      (self->nparams != self->a_nparams) ||
      (self->mk != __MK (self->method, self->k_type)))
  {
    guint i;

    self->a_size    = self->size;
    self->a_nparams = self->nparams;
    self->mk        = __MK (self->method, self->k_type);

    ncm_stats_dist_clear (&self->sd0);
    ncm_stats_dist_clear (&self->sd1);
    ncm_vector_clear (&self->m2lnp_star);
    ncm_vector_clear (&self->m2lnp_cur);
    ncm_vector_clear (&self->m2lnL_s0);
    ncm_vector_clear (&self->m2lnL_s1);

    g_ptr_array_set_size (self->thetastar, 0);

    g_assert (self->size % 2 == 0);
    self->size_2 = self->size / 2;

    self->m2lnp_star = ncm_vector_new (self->size);
    self->m2lnp_cur  = ncm_vector_new (self->size);
    self->m2lnL_s0   = ncm_vector_new (self->size_2);
    self->m2lnL_s1   = ncm_vector_new (self->size_2);

    {
      NcmStatsDistKernel *kernel;

      switch (self->k_type)
      {
        case NCM_FIT_ESMCMC_WALKER_APES_KTYPE_CAUCHY:
          kernel = NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_st_new (self->nparams, 1.0));
          break;
        case NCM_FIT_ESMCMC_WALKER_APES_KTYPE_ST3:
          kernel = NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_st_new (self->nparams, 3.0));
          break;
        case NCM_FIT_ESMCMC_WALKER_APES_KTYPE_GAUSS:
          kernel = NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_gauss_new (self->nparams));
          break;
        default:
          g_assert_not_reached ();
          break;
      }

      switch (self->method)
      {
        case NCM_FIT_ESMCMC_WALKER_APES_METHOD_KDE:
        {
          self->sd0 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, NCM_STATS_DIST_CV_NONE));
          self->sd1 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, NCM_STATS_DIST_CV_NONE));
          break;
        }
        case NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE:
        {
          self->sd0 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, NCM_STATS_DIST_CV_NONE));
          self->sd1 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, NCM_STATS_DIST_CV_NONE));
          _ncm_fit_esmcmc_walker_apes_vkde_check_sizes (walker);
          break;
        }
        default:
          g_assert_not_reached ();
          break;
      }

      ncm_stats_dist_kernel_free (kernel);
    }

    ncm_stats_dist_set_over_smooth (self->sd0, self->over_smooth);
    ncm_stats_dist_set_over_smooth (self->sd1, self->over_smooth);

    ncm_stats_dist_set_shrink (self->sd0, self->shrink);
    ncm_stats_dist_set_shrink (self->sd1, self->shrink);

    ncm_stats_dist_set_use_threads (self->sd0, self->use_threads);
    ncm_stats_dist_set_use_threads (self->sd1, self->use_threads);

    for (i = 0; i < self->size; i++)
    {
      NcmVector *thetastar_i = ncm_vector_new (self->nparams);

      g_ptr_array_add (self->thetastar, thetastar_i);
    }
  }
}

static void
_ncm_fit_esmcmc_walker_apes_set_size (NcmFitESMCMCWalker *walker, guint size)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  g_assert_cmpuint (size, >, 0);
  self->size = size;

  if (self->constructed)
    _ncm_fit_esmcmc_walker_apes_set_sys (walker);
}

static guint
_ncm_fit_esmcmc_walker_apes_get_size (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->size;
}

static void
_ncm_fit_esmcmc_walker_apes_set_nparams (NcmFitESMCMCWalker *walker, guint nparams)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  g_assert_cmpuint (nparams, >, 0);
  self->nparams = nparams;

  if (self->constructed)
    _ncm_fit_esmcmc_walker_apes_set_sys (walker);
}

static guint
_ncm_fit_esmcmc_walker_apes_get_nparams (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->nparams;
}

static void
_ncm_fit_esmcmc_walker_apes_prepare_random_walk (NcmFitESMCMCWalker *walker, NcmStatsDist *sd, NcmFitESMCMCWalkerAPESRandomWalk *random_walk, NcmMSet *mset)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if (self->random_walk_prob > 0.0)
  {
    NcmMatrix *cov = ncm_stats_dist_peek_full_cov (sd);
    guint i;

    if ((random_walk->std == NULL) || (ncm_vector_len (random_walk->std) != self->nparams))
    {
      ncm_vector_clear (&random_walk->std);
      ncm_vector_clear (&random_walk->lb);
      ncm_vector_clear (&random_walk->ub);

      random_walk->std = ncm_vector_new (self->nparams);
      random_walk->lb  = ncm_vector_new (self->nparams);
      random_walk->ub  = ncm_vector_new (self->nparams);
    }

    for (i = 0; i < self->nparams; i++)
    {
      const gdouble var = ncm_matrix_get (cov, i, i);
      const gdouble lb  = ncm_mset_fparam_get_lower_bound (mset, i);
      const gdouble ub  = ncm_mset_fparam_get_upper_bound (mset, i);

      ncm_vector_set (random_walk->lb, i, lb);
      ncm_vector_set (random_walk->ub, i, ub);

      /*
       * The standard deviation is the square root of the variance. If the variance is
       * non-positive, then the standard deviation is undefined.
       */
      if (var <= 0.0)
        g_error ("Invalid covariance matrix: diagonal element %d is non-positive.", i);

      /*
       * 0.25 is a scaling factor for the standard deviation. TODO: Make this a
       * parameter.
       */
      ncm_vector_set (random_walk->std, i, sqrt (var) * 0.25);
    }
  }
}

static void
_ncm_fit_esmcmc_walker_apes_random_walk_sample (NcmFitESMCMCWalker *walker, NcmStatsDist *sd, NcmFitESMCMCWalkerAPESRandomWalk *random_walk, const NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);
  guint i;

  for (i = 0; i < self->nparams; i++)
  {
    const gdouble lb      = ncm_vector_fast_get (random_walk->lb, i);
    const gdouble ub      = ncm_vector_fast_get (random_walk->ub, i);
    const gdouble std     = ncm_vector_fast_get (random_walk->std, i);
    const gdouble theta_i = ncm_vector_get (theta, i);
    gdouble x;

    do {
      x = ncm_rng_gaussian_gen (rng, theta_i, std);
    } while ((x < lb) || (x > ub));

    ncm_vector_set (thetastar, i, x);
  }
}

static void
_ncm_fit_esmcmc_walker_apes_sample (NcmFitESMCMCWalker *walker, NcmStatsDist *sd, NcmMSet *mset, NcmFitESMCMCWalkerAPESRandomWalk *random_walk, const NcmVector *theta, NcmVector *thetastar, NcmRNG *rng)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if (self->random_walk_prob == 0.0)
  {
    do {
      ncm_stats_dist_sample (sd, thetastar, rng);
    } while (!ncm_mset_fparam_valid_bounds (mset, thetastar));
  }
  else
  {
    do {
      if (ncm_rng_uniform01_pos_gen (rng) < self->random_walk_prob)
        _ncm_fit_esmcmc_walker_apes_random_walk_sample (walker, sd, random_walk, theta, thetastar, rng);
      else
        ncm_stats_dist_sample (sd, thetastar, rng);

      /* Ensure the sampled point is within bounds */
    } while (!ncm_mset_fparam_valid_bounds (mset, thetastar));
  }
}

static void
_ncm_fit_esmcmc_walker_apes_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);
  const gdouble T                            = 1.0;
  guint i;

  if (ki < self->size_2)
  {
    ncm_stats_dist_reset (self->sd0);

    for (i = self->size_2; i < self->size; i++)
    {
      ncm_vector_set (self->m2lnL_s0, i - self->size_2, (1.0 / T) * ncm_vector_get (g_ptr_array_index (m2lnL, i), 0));
    }

    for (i = self->size_2; i < self->size; i++)
    {
      NcmVector *theta_i = g_ptr_array_index (theta, i);

      ncm_vector_set (self->m2lnL_s0, i - self->size_2, ncm_vector_get (g_ptr_array_index (m2lnL, i), 0));
      ncm_stats_dist_add_obs (self->sd0, theta_i);
    }

    if (self->use_interp)
      ncm_stats_dist_prepare_interp (self->sd0, self->m2lnL_s0);
    else
      ncm_stats_dist_prepare (self->sd0);

    _ncm_fit_esmcmc_walker_apes_prepare_random_walk (walker, self->sd0, &self->rw0, mset);

    for (i = ki; i < self->size_2; i++)
    {
      NcmVector *theta_i     = g_ptr_array_index (theta, i);
      NcmVector *thetastar_i = g_ptr_array_index (self->thetastar, i);

      _ncm_fit_esmcmc_walker_apes_sample (walker, self->sd0, mset, &self->rw0, theta_i, thetastar_i, rng);
    }
  }

  if (kf > self->size_2)
  {
    ncm_stats_dist_reset (self->sd1);

    for (i = 0; i < self->size_2; i++)
    {
      ncm_vector_set (self->m2lnL_s1, i, (1.0 / T) * ncm_vector_get (g_ptr_array_index (m2lnL, i), 0));
    }

    for (i = 0; i < self->size_2; i++)
    {
      NcmVector *theta_i = g_ptr_array_index (theta, i);

      ncm_vector_set (self->m2lnL_s1, i, ncm_vector_get (g_ptr_array_index (m2lnL, i), 0));
      ncm_stats_dist_add_obs (self->sd1, theta_i);
    }

    if (self->use_interp)
      ncm_stats_dist_prepare_interp (self->sd1, self->m2lnL_s1);
    else
      ncm_stats_dist_prepare (self->sd1);

    _ncm_fit_esmcmc_walker_apes_prepare_random_walk (walker, self->sd1, &self->rw1, mset);

    for (i = self->size_2; i < kf; i++)
    {
      NcmVector *theta_i     = g_ptr_array_index (theta, i);
      NcmVector *thetastar_i = g_ptr_array_index (self->thetastar, i);

      _ncm_fit_esmcmc_walker_apes_sample (walker, self->sd1, mset, &self->rw1, theta_i, thetastar_i, rng);
    }
  }

  if (self->exploration > 0)
    self->exploration--;
}

static gdouble
_ncm_fit_esmcmc_walker_apes_transition_prob (NcmFitESMCMCWalker *walker, NcmStatsDist *sd, NcmFitESMCMCWalkerAPESRandomWalk *random_walk, const NcmVector *theta, NcmVector *thetastar)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);
  gdouble m2lnp, sign;

  if (self->random_walk_prob > 0.0)
  {
    gdouble m2lnp_sd = ncm_stats_dist_eval_m2lnp (sd, thetastar);
    gdouble m2lnp_rw = 0.0;
    guint i;

    for (i = 0; i < self->nparams; i++)
    {
      const gdouble lb          = ncm_vector_fast_get (random_walk->lb, i);
      const gdouble ub          = ncm_vector_fast_get (random_walk->ub, i);
      const gdouble std         = ncm_vector_fast_get (random_walk->std, i);
      const gdouble theta_i     = ncm_vector_get (theta, i);
      const gdouble thetastar_i = ncm_vector_get (thetastar, i);
      const gdouble ln_norm     = 0.5 * ncm_c_ln2pi () + log (std) + ncm_util_log_gaussian_integral (lb, ub, theta_i, std, &sign);

      m2lnp_rw += gsl_pow_2 ((thetastar_i - theta_i) / std) + 2.0 * ln_norm;
    }

    m2lnp_rw += -2.0 * log (self->random_walk_prob);
    m2lnp_sd += -2.0 * log1p (-self->random_walk_prob);

    if (m2lnp_sd < m2lnp_rw)
      m2lnp = m2lnp_sd - 2.0 * log1p (exp (-0.5 * (m2lnp_rw - m2lnp_sd)));
    else
      m2lnp = m2lnp_rw - 2.0 * log1p (exp (-0.5 * (m2lnp_sd - m2lnp_rw)));
  }
  else
  {
    m2lnp = ncm_stats_dist_eval_m2lnp (sd, thetastar);
  }


  return m2lnp;
}

static void
_ncm_fit_esmcmc_walker_apes_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);
  NcmVector *theta_k                         = g_ptr_array_index (theta, k);

  ncm_vector_memcpy (thetastar, g_ptr_array_index (self->thetastar, k));

  if (k < self->size_2)
  {
    const gdouble m2lnapes_star = _ncm_fit_esmcmc_walker_apes_transition_prob (walker, self->sd0, &self->rw0, theta_k, thetastar);
    const gdouble m2lnapes_cur  = _ncm_fit_esmcmc_walker_apes_transition_prob (walker, self->sd0, &self->rw0, thetastar, theta_k);

    g_assert (gsl_finite (m2lnapes_star) && gsl_finite (m2lnapes_cur));

    ncm_vector_set (self->m2lnp_star, k, m2lnapes_star);
    ncm_vector_set (self->m2lnp_cur,  k, m2lnapes_cur);
  }

  if (k >= self->size_2)
  {
    const gdouble m2lnapes_star = _ncm_fit_esmcmc_walker_apes_transition_prob (walker, self->sd1, &self->rw1, theta_k, thetastar);
    const gdouble m2lnapes_cur  = _ncm_fit_esmcmc_walker_apes_transition_prob (walker, self->sd1, &self->rw1, thetastar, theta_k);

    g_assert (gsl_finite (m2lnapes_star) && gsl_finite (m2lnapes_cur));

    ncm_vector_set (self->m2lnp_star, k, m2lnapes_star);
    ncm_vector_set (self->m2lnp_cur,  k, m2lnapes_cur);
  }
}

static gdouble
_ncm_fit_esmcmc_walker_apes_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);
  const gdouble m2lnp_star                   = ncm_vector_get (self->m2lnp_star, k);
  const gdouble m2lnp_cur                    = ncm_vector_get (self->m2lnp_cur, k);

  if (self->exploration)
    return exp (-0.5 * (m2lnL_star - m2lnL_cur));
  else
    return exp (-0.5 * ((m2lnL_star - m2lnp_star) - (m2lnL_cur - m2lnp_cur)));
}

static gdouble
_ncm_fit_esmcmc_walker_apes_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);
  const gdouble m2lnp_star                   = ncm_vector_get (self->m2lnp_star, k);
  const gdouble m2lnp_cur                    = ncm_vector_get (self->m2lnp_cur, k);

  if (self->exploration)
    return 0.0;
  else
    return -0.5 * (m2lnp_cur - m2lnp_star);
}

static void
_ncm_fit_esmcmc_walker_apes_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  /*NcmFitESMCMCWalkerAPES *apes = NCM_FIT_ESMCMC_WALKER_APES (walker);*/
  /*NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);*/

  /* Nothing to do. */
}

const gchar *
_ncm_fit_esmcmc_walker_apes_desc (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  gchar *kernel, *method;

  switch (self->method)
  {
    case NCM_FIT_ESMCMC_WALKER_APES_METHOD_KDE:
      method = g_strdup ("KDE");
      break;
    case NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE:
      method = g_strdup ("VKDE");
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  {
    gchar *tmp = method;

    switch (ncm_stats_dist_kde_get_cov_type (NCM_STATS_DIST_KDE (self->sd0)))
    {
      case NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE:
        method = g_strdup (method);
        break;
      case NCM_STATS_DIST_KDE_COV_TYPE_FIXED:
        method = g_strdup_printf ("Fixed-%s", method);
        break;
      case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG:
        method = g_strdup_printf ("RobustDiag-%s", method);
        break;
      case NCM_STATS_DIST_KDE_COV_TYPE_ROBUST:
        method = g_strdup_printf ("Robust-%s", method);
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    g_free (tmp);
  }

  if (self->use_interp)
  {
    gchar *tmp = method;

    method = g_strdup_printf ("Interp-%s", method);
    g_free (tmp);
  }

  switch (self->k_type)
  {
    case NCM_FIT_ESMCMC_WALKER_APES_KTYPE_CAUCHY:
      kernel = "Cauchy";
      break;
    case NCM_FIT_ESMCMC_WALKER_APES_KTYPE_ST3:
      kernel = "ST3";
      break;
    case NCM_FIT_ESMCMC_WALKER_APES_KTYPE_GAUSS:
      kernel = "Gauss";
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  g_clear_pointer (&self->desc, g_free);
  self->desc = g_strdup_printf ("APES-Move:%s:%s", method, kernel);
  g_free (method);

  return self->desc;
}

/**
 * ncm_fit_esmcmc_walker_apes_new:
 * @nwalkers: number of walkers
 * @nparams: number of parameters
 *
 * Creates a new #NcmFitESMCMCWalkerAPES to be used
 * with @nwalkers.
 *
 * Returns: (transfer full): a new #NcmFitESMCMCWalkerAPES.
 */
NcmFitESMCMCWalkerAPES *
ncm_fit_esmcmc_walker_apes_new (guint nwalkers, guint nparams)
{
  NcmFitESMCMCWalkerAPES *apes = g_object_new (NCM_TYPE_FIT_ESMCMC_WALKER_APES,
                                               "size", nwalkers,
                                               "nparams", nparams,
                                               NULL);

  return apes;
}

/**
 * ncm_fit_esmcmc_walker_apes_new_full:
 * @nwalkers: number of walkers
 * @nparams: number of parameters
 * @method: a #NcmFitESMCMCWalkerAPESMethod
 * @k_type: a #NcmFitESMCMCWalkerAPESKType
 * @over_smooth: a double
 * @use_interp: a boolean
 *
 * Creates a new #NcmFitESMCMCWalkerAPES to be used with @nwalkers,
 * interpolation method @method, kernel @kernel and over-smooth parameter
 * @over_smooth. If @use_interp is TRUE computes the approximation
 * interpolating the computed likelihood values, otherwise, use standard
 * kernel density estimation.
 *
 * Returns: (transfer full): a new #NcmFitESMCMCWalkerAPES.
 */
NcmFitESMCMCWalkerAPES *
ncm_fit_esmcmc_walker_apes_new_full (guint nwalkers, guint nparams, NcmFitESMCMCWalkerAPESMethod method, NcmFitESMCMCWalkerAPESKType k_type, gdouble over_smooth, gboolean use_interp)
{
  NcmFitESMCMCWalkerAPES *apes = g_object_new (NCM_TYPE_FIT_ESMCMC_WALKER_APES,
                                               "size",        nwalkers,
                                               "nparams",     nparams,
                                               "method",      method,
                                               "kernel-type", k_type,
                                               "over-smooth", over_smooth,
                                               "use-interp",  use_interp,
                                               NULL);

  return apes;
}

/**
 * ncm_fit_esmcmc_walker_apes_ref:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Increases the reference count of @apes atomically.
 *
 * Returns: (transfer full): @apes.
 */
NcmFitESMCMCWalkerAPES *
ncm_fit_esmcmc_walker_apes_ref (NcmFitESMCMCWalkerAPES *apes)
{
  return g_object_ref (apes);
}

/**
 * ncm_fit_esmcmc_walker_apes_free:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Decreases the reference count of @apes atomically.
 *
 */
void
ncm_fit_esmcmc_walker_apes_free (NcmFitESMCMCWalkerAPES *apes)
{
  g_object_unref (apes);
}

/**
 * ncm_fit_esmcmc_walker_apes_clear:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Decreases the reference count of *@apes atomically and sets the pointer *@apes to null.
 *
 */
void
ncm_fit_esmcmc_walker_apes_clear (NcmFitESMCMCWalkerAPES **apes)
{
  g_clear_object (apes);
}

/**
 * ncm_fit_esmcmc_walker_apes_set_method:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @method: a #NcmFitESMCMCWalkerAPESMethod
 *
 * Sets the estimation method to be used when building the
 * posterior approximations.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_method (NcmFitESMCMCWalkerAPES *apes, NcmFitESMCMCWalkerAPESMethod method)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if (method >= NCM_FIT_ESMCMC_WALKER_APES_METHOD_LEN)
    g_error ("ncm_fit_esmcmc_walker_apes_set_method: invalid method `%d'.", method);

  self->method = method;

  if (self->constructed)
    _ncm_fit_esmcmc_walker_apes_set_sys (NCM_FIT_ESMCMC_WALKER (apes));
}

/**
 * ncm_fit_esmcmc_walker_apes_set_k_type:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @k_type: a #NcmFitESMCMCWalkerAPESKType
 *
 * Sets the kernel to be used when building the
 * posterior approximations.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_k_type (NcmFitESMCMCWalkerAPES *apes, NcmFitESMCMCWalkerAPESKType k_type)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if (k_type >= NCM_FIT_ESMCMC_WALKER_APES_KTYPE_LEN)
    g_error ("ncm_fit_esmcmc_walker_apes_set_method: invalid method `%d'.", k_type);

  self->k_type = k_type;

  if (self->constructed)
    _ncm_fit_esmcmc_walker_apes_set_sys (NCM_FIT_ESMCMC_WALKER (apes));
}

/**
 * ncm_fit_esmcmc_walker_apes_set_over_smooth:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @os: a double
 *
 * Sets the over smooth parameter to adjust the interpolation
 * bandwidth.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_over_smooth (NcmFitESMCMCWalkerAPES *apes, const gdouble os)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->over_smooth = os;

  if (self->constructed)
  {
    ncm_stats_dist_set_over_smooth (self->sd0, self->over_smooth);
    ncm_stats_dist_set_over_smooth (self->sd1, self->over_smooth);
  }
}

/**
 * ncm_fit_esmcmc_walker_apes_set_shrink:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @shrink: a double
 *
 * Sets the shrink parameter to adjust the kernel weights. The value must be between 0.0
 * and 1.0. See ncm_stats_dist_set_shrink() for more details.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_shrink (NcmFitESMCMCWalkerAPES *apes, const gdouble shrink)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if ((shrink < 0.0) || (shrink > 1.0))
    g_error ("ncm_fit_esmcmc_walker_apes_set_shrink: invalid shrink `%f'.", shrink);

  self->shrink = shrink;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_random_walk_prob:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @prob: a double
 *
 * Sets the probability of performing a random walk step. The value must be between 0.0
 * and 1.0.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_random_walk_prob (NcmFitESMCMCWalkerAPES *apes, const gdouble prob)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if ((prob < 0.0) || (prob > 1.0))
    g_error ("ncm_fit_esmcmc_walker_apes_set_random_walk_prob: invalid probability `%f'.", prob);

  self->random_walk_prob = prob;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_random_walk_scale:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @scale: a double
 *
 * Sets the scale factor for the random walk step. The value must be greater than 0.0.
 * This factor multiplies the standard deviation used in the random walk proposal.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_random_walk_scale (NcmFitESMCMCWalkerAPES *apes, const gdouble scale)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if (scale <= 0.0)
    g_error ("ncm_fit_esmcmc_walker_apes_set_random_walk_scale: invalid scale `%f'.", scale);

  self->random_walk_scale = scale;
}

/**
 * ncm_fit_esmcmc_walker_apes_get_method:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Gets the currently used estimation method.
 *
 * Returns: currently used estimation method #NcmFitESMCMCWalkerAPESMethod.
 */
NcmFitESMCMCWalkerAPESMethod
ncm_fit_esmcmc_walker_apes_get_method (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->method;
}

/**
 * ncm_fit_esmcmc_walker_apes_get_k_type:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Gets the currently used kernel.
 *
 * Returns: currently used kernel #NcmFitESMCMCWalkerAPESKType.
 */
NcmFitESMCMCWalkerAPESKType
ncm_fit_esmcmc_walker_apes_get_k_type (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->k_type;
}

/**
 * ncm_fit_esmcmc_walker_apes_get_over_smooth:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Gets the currently used over-smooth parameter.
 *
 * Returns: currently used over-smooth.
 */
gdouble
ncm_fit_esmcmc_walker_apes_get_over_smooth (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->over_smooth;
}

/**
 * ncm_fit_esmcmc_walker_apes_get_shrink:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Gets the currently used shrink parameter.
 *
 * Returns: currently used shrink.
 */
gdouble
ncm_fit_esmcmc_walker_apes_get_shrink (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->shrink;
}

/**
 * ncm_fit_esmcmc_walker_apes_get_random_walk_prob:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Gets the currently used random walk probability.
 *
 * Returns: currently used random walk probability.
 */
gdouble
ncm_fit_esmcmc_walker_apes_get_random_walk_prob (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->random_walk_prob;
}

/**
 * ncm_fit_esmcmc_walker_apes_get_random_walk_scale:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Gets the currently used random walk scale.
 *
 * Returns: currently used random walk scale.
 */
gdouble
ncm_fit_esmcmc_walker_apes_get_random_walk_scale (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->random_walk_scale;
}

/**
 * ncm_fit_esmcmc_walker_apes_use_interp:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @use_interp: whether to use interpolation of the posterior
 *
 * Sets whether to use interpolation of the posterior approximation (@use_interp == TRUE)
 * or kernel density estimate (@use_interp == FALSE).
 *
 */
void
ncm_fit_esmcmc_walker_apes_use_interp (NcmFitESMCMCWalkerAPES *apes, gboolean use_interp)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->use_interp = use_interp;
}

/**
 * ncm_fit_esmcmc_walker_apes_interp:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: whether interpolation is being used for posterior approximation.
 */
gboolean
ncm_fit_esmcmc_walker_apes_interp (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->use_interp;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_use_threads:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @use_threads: whether to use threads
 *
 * Sets whether to use threads for building the posterior
 * approximation.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_use_threads (NcmFitESMCMCWalkerAPES *apes, gboolean use_threads)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->use_threads = use_threads;

  if (self->constructed)
  {
    ncm_stats_dist_set_use_threads (self->sd0, self->use_threads);
    ncm_stats_dist_set_use_threads (self->sd1, self->use_threads);
  }
}

/**
 * ncm_fit_esmcmc_walker_apes_get_use_threads:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: whether threads are being used for building the posterior
 * approximation.
 */
gboolean
ncm_fit_esmcmc_walker_apes_get_use_threads (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);
  gboolean use_threads0                      = ncm_stats_dist_get_use_threads (self->sd0);
  gboolean use_threads1                      = ncm_stats_dist_get_use_threads (self->sd1);

  g_assert (self->use_threads == use_threads0);
  g_assert (use_threads0 == use_threads1);

  return use_threads0;
}

/**
 * ncm_fit_esmcmc_walker_apes_peek_sds:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @sd0: (out) (transfer none): a #NcmStatsDist
 * @sd1: (out) (transfer none): a #NcmStatsDist
 *
 * Peeks the currently used #NcmStatsDist objects.
 *
 */
void
ncm_fit_esmcmc_walker_apes_peek_sds (NcmFitESMCMCWalkerAPES *apes, NcmStatsDist **sd0, NcmStatsDist **sd1)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  g_assert (self->constructed);

  sd0[0] = self->sd0;
  sd1[0] = self->sd1;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_local_frac:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @local_frac: a double determining the local fraction to use in VKDE.
 *
 * Sets the local fraction to use in VKDE.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_local_frac (NcmFitESMCMCWalkerAPES *apes, gdouble local_frac)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if (self->method != NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE)
    g_error ("ncm_fit_esmcmc_walker_apes_set_local_frac: cannot set local fraction for a non-VKDE method.");

  ncm_stats_dist_vkde_set_local_frac (NCM_STATS_DIST_VKDE (self->sd0), local_frac);
  ncm_stats_dist_vkde_set_local_frac (NCM_STATS_DIST_VKDE (self->sd1), local_frac);

  _ncm_fit_esmcmc_walker_apes_vkde_check_sizes (NCM_FIT_ESMCMC_WALKER (apes));
}

/**
 * ncm_fit_esmcmc_walker_apes_set_cov_fixed_from_mset:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @mset: a #NcmMSet to get the covariance from.
 *
 * Sets the fixed covariance to the KDE interpolation using
 * the scales set into @mset.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_cov_fixed_from_mset (NcmFitESMCMCWalkerAPES *apes, NcmMSet *mset)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);
  NcmMatrix *cov_fixed                       = ncm_matrix_new (self->nparams, self->nparams);
  guint i;

  ncm_matrix_set_identity (cov_fixed);

  for (i = 0; i < self->nparams; i++)
  {
    const gdouble scale = ncm_mset_fparam_get_scale (mset, i);

    ncm_matrix_set (cov_fixed, i, i, scale * scale);
  }

  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd0), NCM_STATS_DIST_KDE_COV_TYPE_FIXED);
  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd1), NCM_STATS_DIST_KDE_COV_TYPE_FIXED);

  ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (self->sd0), cov_fixed);
  ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (self->sd1), cov_fixed);

  ncm_matrix_free (cov_fixed);
}

/**
 * ncm_fit_esmcmc_walker_apes_set_cov_robust_diag:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Sets the fixed covariance to the KDE interpolation using
 * robust estimates of scale.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_cov_robust_diag (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd0), NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG);
  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd1), NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG);
}

/**
 * ncm_fit_esmcmc_walker_apes_set_cov_robust:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Sets the fixed covariance to the KDE interpolation using
 * robust estimates of scale.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_cov_robust (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd0), NCM_STATS_DIST_KDE_COV_TYPE_ROBUST);
  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd1), NCM_STATS_DIST_KDE_COV_TYPE_ROBUST);
}

/**
 * ncm_fit_esmcmc_walker_apes_set_exploration:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @exploration: a guint
 *
 * Sets the exploration parameter to be used when building the
 * posterior approximation. During the exploration phase, the
 * new samples are accepted considering only the posterior.
 * This makes the exploration phase to be more efficient to find
 * the global maximum of the posterior, but this phase should be
 * discarded in the final analysis.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_exploration (NcmFitESMCMCWalkerAPES *apes, guint exploration)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->exploration = exploration;
}

