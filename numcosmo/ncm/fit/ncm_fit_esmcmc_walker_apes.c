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

#include "ncm/fit/ncm_fit_esmcmc_walker.h"
#include "ncm/fit/ncm_fit_esmcmc_walker_apes.h"

#include "ncm/core/ncm_c.h"
#include "ncm/fit/ncm_fit_esmcmc.h"
#include "ncm/stats/ncm_stats_dist_kde.h"
#include "ncm/stats/ncm_stats_dist_vkde.h"
#include "ncm/stats/ncm_stats_dist_kernel_st.h"
#include "ncm/stats/ncm_stats_dist_kernel_gauss.h"

#include "ncm_enum_types.h"

#ifndef NUMCOSMO_GIR_SCAN
#endif /* NUMCOSMO_GIR_SCAN */


enum
{
  PROP_0,
  PROP_METHOD,
  PROP_K_TYPE,
  PROP_OVER_SMOOTH,
  PROP_USE_INTERP,
  PROP_USE_THREADS,
  PROP_CENTER_SHRINK,
  PROP_DEFENSIVE_FRAC,
  PROP_DEFENSIVE_SCALE,
  PROP_DEFENSIVE_NU,
  PROP_VKDE_POINTS_PER_DIM,
  PROP_UNIFORM_WEIGHTS,
  PROP_CV_TYPE,
  PROP_SPLIT_FRAC,
  PROP_EXPLORATION,
  PROP_EXPLORATION_QRATIO_FLOOR,
  PROP_EXPLORATION_STOP_AFTER,
};

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
  gboolean use_interp;
  gboolean use_threads;
  gboolean center_shrink;
  gdouble defensive_frac;
  gdouble defensive_scale;
  gdouble defensive_nu;
  gdouble vkde_points_per_dim;
  gboolean uniform_weights;
  gdouble local_frac;
  NcmStatsDistCV cv_type;
  gdouble split_frac;
  NcmStatsDistKDECovType cov_type;
  NcmMatrix *cov_fixed;
  gboolean constructed;
  guint exploration;                /* cap on the exploration phase, iterations; 0: no cap */
  gdouble exploration_qratio_floor; /* floor of q(x)/q(x') while exploring; 0: posterior-only acceptance */
  guint exploration_stop_after;     /* iterations without a clipped acceptance that end the phase */
  gboolean exploring;               /* phase armed */
  guint expl_iters;                 /* iterations elapsed in the phase */
  guint unclipped_iters;            /* consecutive iterations in which no acceptance was clipped */
  GArray *clipped;                  /* per-walker: acceptance modified in the current iteration */
  gboolean last_markovian;          /* the last completed iteration used the exact acceptance everywhere */
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

  self->size                     = 0;
  self->size_2                   = 0;
  self->nparams                  = 0;
  self->a_size                   = 0;
  self->a_nparams                = 0;
  self->mk                       = -1;
  self->m2lnp_star               = NULL;
  self->m2lnp_cur                = NULL;
  self->desc                     = NULL;
  self->sd0                      = NULL;
  self->sd1                      = NULL;
  self->thetastar                = g_ptr_array_new ();
  self->m2lnL_s0                 = NULL;
  self->m2lnL_s1                 = NULL;
  self->method                   = NCM_FIT_ESMCMC_WALKER_APES_METHOD_LEN;
  self->k_type                   = NCM_FIT_ESMCMC_WALKER_APES_KTYPE_LEN;
  self->over_smooth              = 0.0;
  self->use_interp               = FALSE;
  self->use_threads              = FALSE;
  self->center_shrink            = FALSE;
  self->defensive_frac           = 0.0;
  self->defensive_scale          = 4.0;
  self->defensive_nu             = 3.0;
  self->vkde_points_per_dim      = 0.0;
  self->uniform_weights          = FALSE;
  self->local_frac               = 0.0;
  self->cv_type                  = NCM_STATS_DIST_CV_NONE;
  self->split_frac               = 0.0;
  self->cov_type                 = NCM_STATS_DIST_KDE_COV_TYPE_SAMPLE;
  self->cov_fixed                = NULL;
  self->constructed              = FALSE;
  self->exploration              = 0;
  self->exploration_qratio_floor = 0.0;
  self->exploration_stop_after   = 10;
  self->exploring                = FALSE;
  self->expl_iters               = 0;
  self->unclipped_iters          = 0;
  self->clipped                  = g_array_new (FALSE, TRUE, sizeof (gboolean));
  self->last_markovian           = TRUE;


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
    case PROP_USE_INTERP:
      ncm_fit_esmcmc_walker_apes_use_interp (apes, g_value_get_boolean (value));
      break;
    case PROP_USE_THREADS:
      ncm_fit_esmcmc_walker_apes_set_use_threads (apes, g_value_get_boolean (value));
      break;
    case PROP_CENTER_SHRINK:
      ncm_fit_esmcmc_walker_apes_set_center_shrink (apes, g_value_get_boolean (value));
      break;
    case PROP_DEFENSIVE_FRAC:
      ncm_fit_esmcmc_walker_apes_set_defensive_frac (apes, g_value_get_double (value));
      break;
    case PROP_DEFENSIVE_SCALE:
      ncm_fit_esmcmc_walker_apes_set_defensive_scale (apes, g_value_get_double (value));
      break;
    case PROP_DEFENSIVE_NU:
      ncm_fit_esmcmc_walker_apes_set_defensive_nu (apes, g_value_get_double (value));
      break;
    case PROP_VKDE_POINTS_PER_DIM:
      ncm_fit_esmcmc_walker_apes_set_vkde_points_per_dim (apes, g_value_get_double (value));
      break;
    case PROP_UNIFORM_WEIGHTS:
      ncm_fit_esmcmc_walker_apes_set_uniform_weights (apes, g_value_get_boolean (value));
      break;
    case PROP_CV_TYPE:
      ncm_fit_esmcmc_walker_apes_set_cv_type (apes, g_value_get_enum (value));
      break;
    case PROP_SPLIT_FRAC:
      ncm_fit_esmcmc_walker_apes_set_split_frac (apes, g_value_get_double (value));
      break;
    case PROP_EXPLORATION:
      ncm_fit_esmcmc_walker_apes_set_exploration (apes, g_value_get_uint (value));
      break;
    case PROP_EXPLORATION_QRATIO_FLOOR:
      ncm_fit_esmcmc_walker_apes_set_exploration_qratio_floor (apes, g_value_get_double (value));
      break;
    case PROP_EXPLORATION_STOP_AFTER:
      ncm_fit_esmcmc_walker_apes_set_exploration_stop_after (apes, g_value_get_uint (value));
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
    case PROP_USE_INTERP:
      g_value_set_boolean (value, ncm_fit_esmcmc_walker_apes_interp (apes));
      break;
    case PROP_USE_THREADS:
      g_value_set_boolean (value, ncm_fit_esmcmc_walker_apes_get_use_threads (apes));
      break;
    case PROP_CENTER_SHRINK:
      g_value_set_boolean (value, ncm_fit_esmcmc_walker_apes_get_center_shrink (apes));
      break;
    case PROP_DEFENSIVE_FRAC:
      g_value_set_double (value, ncm_fit_esmcmc_walker_apes_get_defensive_frac (apes));
      break;
    case PROP_DEFENSIVE_SCALE:
      g_value_set_double (value, ncm_fit_esmcmc_walker_apes_get_defensive_scale (apes));
      break;
    case PROP_DEFENSIVE_NU:
      g_value_set_double (value, ncm_fit_esmcmc_walker_apes_get_defensive_nu (apes));
      break;
    case PROP_VKDE_POINTS_PER_DIM:
      g_value_set_double (value, ncm_fit_esmcmc_walker_apes_get_vkde_points_per_dim (apes));
      break;
    case PROP_UNIFORM_WEIGHTS:
      g_value_set_boolean (value, ncm_fit_esmcmc_walker_apes_get_uniform_weights (apes));
      break;
    case PROP_CV_TYPE:
      g_value_set_enum (value, ncm_fit_esmcmc_walker_apes_get_cv_type (apes));
      break;
    case PROP_SPLIT_FRAC:
      g_value_set_double (value, ncm_fit_esmcmc_walker_apes_get_split_frac (apes));
      break;
    case PROP_EXPLORATION:
      g_value_set_uint (value, ncm_fit_esmcmc_walker_apes_get_exploration (apes));
      break;
    case PROP_EXPLORATION_QRATIO_FLOOR:
      g_value_set_double (value, ncm_fit_esmcmc_walker_apes_get_exploration_qratio_floor (apes));
      break;
    case PROP_EXPLORATION_STOP_AFTER:
      g_value_set_uint (value, ncm_fit_esmcmc_walker_apes_get_exploration_stop_after (apes));
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
  ncm_matrix_clear (&self->cov_fixed);


  g_clear_pointer (&self->thetastar, g_ptr_array_unref);

  g_clear_pointer (&self->clipped, g_array_unref);

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
static void _ncm_fit_esmcmc_walker_apes_start_run (NcmFitESMCMCWalker *walker, gboolean initial, guint exploration_done);
static void _ncm_fit_esmcmc_walker_apes_end_run (NcmFitESMCMCWalker *walker);
static gboolean _ncm_fit_esmcmc_walker_apes_is_markovian (NcmFitESMCMCWalker *walker);
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
   * #NcmFitESMCMCWalkerAPESKType values. The default,
   * #NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO, fits the kernel together with the bandwidth
   * and so needs a #NcmFitESMCMCWalkerAPES:cv-type that fits it.
   */
  g_object_class_install_property (object_class,
                                   PROP_K_TYPE,
                                   g_param_spec_enum ("kernel-type",
                                                      NULL,
                                                      "Kernel used in posterior approximation",
                                                      NCM_TYPE_FIT_ESMCMC_WALKER_APES_KTYPE, NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO,
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

  /**
   * NcmFitESMCMCWalkerAPES:center-shrink:
   *
   * Whether to shrink the kernel centres of the posterior approximation toward the
   * ensemble mean, so that the covariance of the approximation equals the covariance
   * of the half-ensemble it is built from for any value of
   * #NcmFitESMCMCWalkerAPES:over-smooth. See #NcmStatsDist:center-shrink.
   *
   * On by default. The gain grows with dimension: on a Gaussian target at 5200
   * walkers it leaves the autocorrelation time at 18 against 60 to 78 without it, and
   * the acceptance at 0.14 against 0.005. It needs a kernel with a finite covariance,
   * so not the Cauchy one.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_CENTER_SHRINK,
                                   g_param_spec_boolean ("center-shrink",
                                                         NULL,
                                                         "Whether to shrink the kernel centres toward the ensemble mean",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:defensive-frac:
   *
   * Weight of the wide Student-t component in the proposal, see
   * #NcmStatsDist:defensive-frac. Default: 0.
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
   * NcmFitESMCMCWalkerAPES:defensive-scale:
   *
   * Covariance factor of the wide component, see #NcmStatsDist:defensive-scale.
   * Default: 4.
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
   * NcmFitESMCMCWalkerAPES:defensive-nu:
   *
   * Degrees of freedom of the wide component, see #NcmStatsDist:defensive-nu.
   * Default: 3.
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
   * NcmFitESMCMCWalkerAPES:vkde-points-per-dim:
   *
   * See #NcmStatsDistVKDE:points-per-dim; ignored for the KDE method. Default: 12.
   * Values around 25 do about as well, so roughly 12 to 25 is the useful range; 0
   * falls back to #NcmFitESMCMCWalkerAPES:local-frac.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_VKDE_POINTS_PER_DIM,
                                   g_param_spec_double ("vkde-points-per-dim",
                                                        NULL,
                                                        "Nearest neighbors per dimension for the VKDE local covariances (0: local fraction)",
                                                        0.0, 1.0e6, 12.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:uniform-weights:
   *
   * See #NcmStatsDist:uniform-weights: with interpolation on, keep uniform kernel weights
   * (the bandwidth fit still runs) instead of the NNLS fit. On by default: it costs a few
   * per cent of autocorrelation time and removes the NNLS solve, which is serial and
   * quadratic in the walker count and so caps how many walkers a high-dimensional
   * problem can afford.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_UNIFORM_WEIGHTS,
                                   g_param_spec_boolean ("uniform-weights",
                                                         NULL,
                                                         "Uniform kernel weights instead of the NNLS fit",
                                                         TRUE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:cv-type:
   *
   * The cross-validation used to choose the over-smooth factor. With
   * #NCM_STATS_DIST_CV_SPLIT_M2LNP the approximation is built from a fraction
   * #NcmFitESMCMCWalkerAPES:split-frac of the block and the over-smooth factor is
   * chosen on the remaining points, which are an independent sample from the same
   * block. This is the default; #NCM_STATS_DIST_CV_NONE uses the over-smooth factor as
   * given.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_CV_TYPE,
                                   g_param_spec_enum ("cv-type",
                                                      NULL,
                                                      "Cross-validation used to choose the over-smooth factor",
                                                      NCM_TYPE_STATS_DIST_CV, NCM_STATS_DIST_CV_SPLIT_M2LNP,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:split-frac:
   *
   * The fraction of the block used as kernel centres when the cross-validation
   * splits the sample. Default: 0.8. Zero keeps whatever #NcmStatsDist uses.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SPLIT_FRAC,
                                   g_param_spec_double ("split-frac",
                                                        NULL,
                                                        "Fraction of the block used as kernel centres",
                                                        0.0, 1.0, 0.8,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:exploration:
   *
   * Length cap of the exploration phase, in iterations (ensemble steps). With
   * #NcmFitESMCMCWalkerAPES:exploration-qratio-floor at 0 the phase accepts by the
   * posterior ratio alone and lasts exactly this many iterations; with a positive floor it
   * ends earlier, after #NcmFitESMCMCWalkerAPES:exploration-stop-after consecutive iterations
   * in which no acceptance was clipped.
   * 0: no cap (and, with the floor at 0, no exploration). The phase is armed only when a
   * run starts from its initial ensemble, never on a continuation, and is never re-armed
   * within a run; the catalog's #NcmMSetCatalog:markovian-id records where it ended.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_EXPLORATION,
                                   g_param_spec_uint ("exploration",
                                                      NULL,
                                                      "Exploration phase length cap in iterations",
                                                      0, G_MAXUINT, 0,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:exploration-qratio-floor:
   *
   * During the exploration phase the proposal-density ratio $q(x)/q(x')$ in the
   * acceptance is clipped from below at this value, so a walker sitting where the
   * proposal has almost no mass (a straggler) can leave; moves whose ratio is above the
   * floor keep the exact Metropolis-Hastings acceptance. An iteration in which any walker
   * was clipped is not Markovian. 0: no clipping (the phase, if any, uses the posterior
   * ratio alone); 1: the posterior ratio alone for every blocked move.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_EXPLORATION_QRATIO_FLOOR,
                                   g_param_spec_double ("exploration-qratio-floor",
                                                        NULL,
                                                        "Floor of q(x)/q(x') during the exploration phase",
                                                        0.0, 1.0, 0.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmFitESMCMCWalkerAPES:exploration-stop-after:
   *
   * Number of consecutive iterations without any clipped acceptance after which the
   * exploration phase ends (the clip is disarmed for the rest of the run).
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_EXPLORATION_STOP_AFTER,
                                   g_param_spec_uint ("exploration-stop-after",
                                                      NULL,
                                                      "Consecutive iterations with no clipped acceptance that end the phase",
                                                      1, G_MAXUINT, 10,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  walker_class->set_size     = &_ncm_fit_esmcmc_walker_apes_set_size;
  walker_class->get_size     = &_ncm_fit_esmcmc_walker_apes_get_size;
  walker_class->set_nparams  = &_ncm_fit_esmcmc_walker_apes_set_nparams;
  walker_class->get_nparams  = &_ncm_fit_esmcmc_walker_apes_get_nparams;
  walker_class->setup        = &_ncm_fit_esmcmc_walker_apes_setup;
  walker_class->step         = &_ncm_fit_esmcmc_walker_apes_step;
  walker_class->prob         = &_ncm_fit_esmcmc_walker_apes_prob;
  walker_class->prob_norm    = &_ncm_fit_esmcmc_walker_apes_prob_norm;
  walker_class->clean        = &_ncm_fit_esmcmc_walker_apes_clean;
  walker_class->desc         = &_ncm_fit_esmcmc_walker_apes_desc;
  walker_class->start_run    = &_ncm_fit_esmcmc_walker_apes_start_run;
  walker_class->end_run      = &_ncm_fit_esmcmc_walker_apes_end_run;
  walker_class->is_markovian = &_ncm_fit_esmcmc_walker_apes_is_markovian;
}

static void
_ncm_fit_esmcmc_walker_apes_vkde_check_sizes (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  guint cov_estimates0 = ncm_stats_dist_vkde_get_n_neighbors (NCM_STATS_DIST_VKDE (self->sd0), self->size_2);
  guint cov_estimates1 = ncm_stats_dist_vkde_get_n_neighbors (NCM_STATS_DIST_VKDE (self->sd1), self->size_2);
  const gdouble ppd    = ncm_stats_dist_vkde_get_points_per_dim (NCM_STATS_DIST_VKDE (self->sd0));

  /* With points-per-dim set, the same count that a local covariance needs is the least
   * the half-ensemble must have for the global one: fewer walkers than that and the run
   * does not start. */
  if ((ppd > 0.0) && (self->size_2 < (guint) ceil (ppd * self->nparams)))
    g_error ("Number of walkers per block (%u) is below points-per-dim x dimension (%g x %u = %u): "
             "not enough walkers to estimate a %u x %u covariance. Increase nwalkers or lower vkde-points-per-dim.",
             self->size_2, ppd, self->nparams, (guint) ceil (ppd * self->nparams), self->nparams, self->nparams);

  if (cov_estimates0 < 2)
    g_error ("Number of walkers per block (%d) is too low for the current dimension (%d).\n"
             "\tToo few points (%d) to estimate local covariances.",
             self->size_2, self->nparams, cov_estimates0);

  if (cov_estimates1 < 2)
    g_error ("Number of walkers per block (%d) is too low for the current dimension (%d).\n"
             "\tToo few points (%d) to estimate local covariances.",
             self->size_2, self->nparams, cov_estimates1);
}

/* Centre shrinkage matches the covariance of the approximation to the ensemble
 * covariance, which the Cauchy kernel does not have. Refuse the combination rather than
 * silently producing a mismatched proposal, or silently dropping a default the caller
 * may be relying on. Since shrinkage is on by default, selecting the Cauchy kernel means
 * turning shrinkage off first. */
static void
_ncm_fit_esmcmc_walker_apes_check_center_shrink (NcmFitESMCMCWalkerAPESPrivate * const self)
{
  if (self->center_shrink && (self->k_type == NCM_FIT_ESMCMC_WALKER_APES_KTYPE_CAUCHY))
    g_error ("Centre shrinkage is on (the default) and the Cauchy kernel has no covariance for it to match.\n"
             "\tShrinkage makes the proposal reproduce the ensemble covariance, which Cauchy does not define: "
             "use the ST3 or GAUSS kernel, or turn centre shrinkage off before setting the kernel type.");
}

/* Only the cross-validations that fit the bandwidth can fit the kernel with it; under any
 * other the kernel would silently stay at its unfitted seed. */
static void
_ncm_fit_esmcmc_walker_apes_check_auto_kernel (NcmFitESMCMCWalkerAPESPrivate * const self)
{
  if (self->k_type != NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO)
    return;

  switch (self->cv_type)
  {
    case NCM_STATS_DIST_CV_SPLIT_M2LNP:
    case NCM_STATS_DIST_CV_SPLIT_ACCEPT:
    case NCM_STATS_DIST_CV_LOO_M2LNP:
      break;
    default:
      g_error ("Kernel selection is set to auto but the cross-validation does not fit the bandwidth.\n"
               "\tAuto chooses the kernel by the same out-of-sample objective that fits the bandwidth: "
               "use one that fits it, split-m2lnp for instance, or name a kernel instead of auto.");
      break;
  }
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
        case NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO:
          /* Where the kernel fit starts from, so the estimator is built at the point the
           * cross-validation will move away from rather than at an unrelated one. */
          kernel = NCM_STATS_DIST_KERNEL (ncm_stats_dist_kernel_st_new (self->nparams, 10.0));
          break;
        default:
          g_assert_not_reached ();
          break;
      }

      switch (self->method)
      {
        case NCM_FIT_ESMCMC_WALKER_APES_METHOD_KDE:
        {
          self->sd0 = NCM_STATS_DIST (ncm_stats_dist_kde_new (kernel, self->cv_type));
          self->sd1 = NCM_STATS_DIST (ncm_stats_dist_kde_new (kernel, self->cv_type));
          break;
        }
        case NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE:
        {
          self->sd0 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, self->cv_type));
          self->sd1 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, self->cv_type));
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

    ncm_stats_dist_set_cv_type (self->sd0, self->cv_type);
    ncm_stats_dist_set_cv_type (self->sd1, self->cv_type);

    if (self->split_frac > 0.0)
    {
      ncm_stats_dist_set_split_frac (self->sd0, self->split_frac);
      ncm_stats_dist_set_split_frac (self->sd1, self->split_frac);
    }

    ncm_stats_dist_set_auto_kernel (self->sd0, self->k_type == NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO);
    ncm_stats_dist_set_auto_kernel (self->sd1, self->k_type == NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO);

    ncm_stats_dist_set_use_threads (self->sd0, self->use_threads);
    ncm_stats_dist_set_use_threads (self->sd1, self->use_threads);

    _ncm_fit_esmcmc_walker_apes_check_center_shrink (self);
    _ncm_fit_esmcmc_walker_apes_check_auto_kernel (self);
    ncm_stats_dist_set_center_shrink (self->sd0, self->center_shrink);
    ncm_stats_dist_set_center_shrink (self->sd1, self->center_shrink);
    ncm_stats_dist_set_defensive_frac (self->sd0, self->defensive_frac);
    ncm_stats_dist_set_defensive_frac (self->sd1, self->defensive_frac);
    ncm_stats_dist_set_defensive_scale (self->sd0, self->defensive_scale);
    ncm_stats_dist_set_defensive_scale (self->sd1, self->defensive_scale);
    ncm_stats_dist_set_defensive_nu (self->sd0, self->defensive_nu);
    ncm_stats_dist_set_defensive_nu (self->sd1, self->defensive_nu);
    ncm_stats_dist_set_uniform_weights (self->sd0, self->uniform_weights);
    ncm_stats_dist_set_uniform_weights (self->sd1, self->uniform_weights);

    /* The objects above have just been created, so every setting that lives inside
     * them has to be applied again; otherwise changing the method or the kernel
     * would silently reset whatever the caller had configured. */
    if ((self->local_frac > 0.0) && (self->method == NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE))
    {
      ncm_stats_dist_vkde_set_local_frac (NCM_STATS_DIST_VKDE (self->sd0), self->local_frac);
      ncm_stats_dist_vkde_set_local_frac (NCM_STATS_DIST_VKDE (self->sd1), self->local_frac);
      _ncm_fit_esmcmc_walker_apes_vkde_check_sizes (walker);
    }

    if (self->method == NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE)
    {
      ncm_stats_dist_vkde_set_points_per_dim (NCM_STATS_DIST_VKDE (self->sd0), self->vkde_points_per_dim);
      ncm_stats_dist_vkde_set_points_per_dim (NCM_STATS_DIST_VKDE (self->sd1), self->vkde_points_per_dim);
      _ncm_fit_esmcmc_walker_apes_vkde_check_sizes (walker);
    }

    ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd0), self->cov_type);
    ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd1), self->cov_type);

    if (self->cov_fixed != NULL)
    {
      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (self->sd0), self->cov_fixed);
      ncm_stats_dist_kde_set_cov_fixed (NCM_STATS_DIST_KDE (self->sd1), self->cov_fixed);
    }

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
_ncm_fit_esmcmc_walker_apes_sample (NcmStatsDist *sd, NcmMSet *mset, NcmVector *thetastar, NcmRNG *rng)
{
  /* Redraw the whole proposal until it falls inside the parameter box: the proposal
   * density is then q(x) / Z with a Z common to every walker of the half, which cancels
   * in the acceptance ratio. */
  do {
    ncm_stats_dist_sample (sd, thetastar, rng);
  } while (!ncm_mset_fparam_valid_bounds (mset, thetastar));
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
      ncm_stats_dist_prepare (self->sd0, self->m2lnL_s0);
    else
      ncm_stats_dist_prepare (self->sd0, NULL);

    for (i = ki; i < self->size_2; i++)
    {
      NcmVector *thetastar_i = g_ptr_array_index (self->thetastar, i);

      _ncm_fit_esmcmc_walker_apes_sample (self->sd0, mset, thetastar_i, rng);
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
      ncm_stats_dist_prepare (self->sd1, self->m2lnL_s1);
    else
      ncm_stats_dist_prepare (self->sd1, NULL);

    for (i = self->size_2; i < kf; i++)
    {
      NcmVector *thetastar_i = g_ptr_array_index (self->thetastar, i);

      _ncm_fit_esmcmc_walker_apes_sample (self->sd1, mset, thetastar_i, rng);
    }
  }
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
    const gdouble m2lnapes_star = ncm_stats_dist_eval_m2lnp (self->sd0, thetastar);
    const gdouble m2lnapes_cur  = ncm_stats_dist_eval_m2lnp (self->sd0, theta_k);

    g_assert (gsl_finite (m2lnapes_star) && gsl_finite (m2lnapes_cur));

    ncm_vector_set (self->m2lnp_star, k, m2lnapes_star);
    ncm_vector_set (self->m2lnp_cur,  k, m2lnapes_cur);
  }

  if (k >= self->size_2)
  {
    const gdouble m2lnapes_star = ncm_stats_dist_eval_m2lnp (self->sd1, thetastar);
    const gdouble m2lnapes_cur  = ncm_stats_dist_eval_m2lnp (self->sd1, theta_k);

    g_assert (gsl_finite (m2lnapes_star) && gsl_finite (m2lnapes_cur));

    ncm_vector_set (self->m2lnp_star, k, m2lnapes_star);
    ncm_vector_set (self->m2lnp_cur,  k, m2lnapes_cur);
  }
}

/*
 * ln [q(x)/q(x')] for walker k, clipped from below during the exploration phase; records the
 * clip per walker (distinct index per thread, no shared state). With the floor at 0 the
 * phase uses the posterior ratio alone, i.e. the q ratio is replaced by 1 for every walker.
 */
static gdouble
_ncm_fit_esmcmc_walker_apes_lnqratio (NcmFitESMCMCWalkerAPESPrivate * const self, guint k)
{
  const gdouble m2lnp_star = ncm_vector_get (self->m2lnp_star, k);
  const gdouble m2lnp_cur  = ncm_vector_get (self->m2lnp_cur, k);
  const gdouble lnqratio   = -0.5 * (m2lnp_cur - m2lnp_star);

  if (!self->exploring)
    return lnqratio;

  if (self->exploration_qratio_floor == 0.0)
  {
    g_array_index (self->clipped, gboolean, k) = TRUE;

    return 0.0;
  }
  else
  {
    const gdouble ln_floor = log (self->exploration_qratio_floor);

    if (lnqratio < ln_floor)
    {
      g_array_index (self->clipped, gboolean, k) = TRUE;

      return ln_floor;
    }

    return lnqratio;
  }
}

static gdouble
_ncm_fit_esmcmc_walker_apes_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return exp (-0.5 * (m2lnL_star - m2lnL_cur) + _ncm_fit_esmcmc_walker_apes_lnqratio (self, k));
}

static gdouble
_ncm_fit_esmcmc_walker_apes_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return _ncm_fit_esmcmc_walker_apes_lnqratio (self, k);
}

static void
_ncm_fit_esmcmc_walker_apes_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  /* clean () is called once per iteration, after both half-ensembles; that is the
   * iteration boundary where the exploration phase is accounted. */
  if (!self->exploring)
  {
    self->last_markovian = TRUE;

    return;
  }
  else
  {
    gboolean any = FALSE;
    guint k;

    for (k = 0; k < self->clipped->len; k++)
    {
      any                                        = any || g_array_index (self->clipped, gboolean, k);
      g_array_index (self->clipped, gboolean, k) = FALSE;
    }

    self->last_markovian = !any;
    self->expl_iters++;
    self->unclipped_iters = any ? 0 : self->unclipped_iters + 1;

    if ((self->exploration > 0) && (self->expl_iters >= self->exploration))
      self->exploring = FALSE;
    else if ((self->exploration_qratio_floor > 0.0) && (self->unclipped_iters >= self->exploration_stop_after))
      self->exploring = FALSE;
  }
}

static void
_ncm_fit_esmcmc_walker_apes_start_run (NcmFitESMCMCWalker *walker, gboolean initial, guint exploration_done)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);
  const gboolean configured                  = (self->exploration > 0) || (self->exploration_qratio_floor > 0.0);
  const gboolean capped_out                  = (self->exploration > 0) && (exploration_done >= self->exploration);

  /* The phase exists to remove the stragglers of the initial ensemble: armed only while
   * the chain has no Markovian rows (a fresh start, or a resume of an interrupted phase),
   * never once Markovian rows exist, and never re-armed within a run. */
  self->exploring       = initial && configured && !capped_out;
  self->expl_iters      = exploration_done;
  self->unclipped_iters = 0;
  self->last_markovian  = TRUE;

  g_array_set_size (self->clipped, 0);
  g_array_set_size (self->clipped, self->size);
}

static void
_ncm_fit_esmcmc_walker_apes_end_run (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->exploring = FALSE;
}

static gboolean
_ncm_fit_esmcmc_walker_apes_is_markovian (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->last_markovian;
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

  if (self->center_shrink)
  {
    gchar *tmp = method;

    method = g_strdup_printf ("Shrink-%s", method);
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
    case NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO:
      kernel = "Auto";
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  {
    /* The settings that change the proposal but are invisible in the name above: which
     * weights the interpolation carries, how the local neighborhood of VKDE is chosen, the
     * cross-validation of the bandwidth, and the over-smoothing. A catalog or a log naming
     * only the structure cannot say which of two runs it came from. */
    const gchar *cv = NULL;
    gchar *opts;

    switch (self->cv_type)
    {
      case NCM_STATS_DIST_CV_NONE:
        cv = "none";
        break;
      case NCM_STATS_DIST_CV_SPLIT_M2LNP:
        cv = "split-m2lnp";
        break;
      case NCM_STATS_DIST_CV_LOO:
        cv = "loo";
        break;
      case NCM_STATS_DIST_CV_SPLIT_ACCEPT:
        cv = "split-accept";
        break;
      case NCM_STATS_DIST_CV_LOO_M2LNP:
        cv = "loo-m2lnp";
        break;
      default:
        g_assert_not_reached ();
        break;
    }

    opts = g_strdup_printf ("os=%g:cv=%s", self->over_smooth, cv);

    /* Weights are only fitted when the kernels are interpolated; without interpolation they
     * are uniform by construction and saying so would suggest a choice was made. */
    if (self->use_interp)
    {
      gchar *tmp = opts;

      opts = g_strdup_printf ("%s:%s", self->uniform_weights ? "unif" : "nnls", tmp);
      g_free (tmp);
    }

    /* The local neighborhood is a VKDE notion: points per dimension when it is positive,
     * otherwise a fraction of the ensemble. */
    if (self->method == NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE)
    {
      gchar *tmp = opts;

      if (self->vkde_points_per_dim > 0.0)
        opts = g_strdup_printf ("ppd=%g:%s", self->vkde_points_per_dim, tmp);
      else
        opts = g_strdup_printf ("lf=%g:%s", self->local_frac, tmp);

      g_free (tmp);
    }

    g_clear_pointer (&self->desc, g_free);
    self->desc = g_strdup_printf ("APES-Move:%s:%s:%s", method, kernel, opts);
    g_free (opts);
  }

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
 * ncm_fit_esmcmc_walker_apes_set_center_shrink:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @center_shrink: whether to shrink the kernel centres toward the ensemble mean
 *
 * Sets whether the posterior approximations use centre shrinkage, see
 * #NcmFitESMCMCWalkerAPES:center-shrink.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_center_shrink (NcmFitESMCMCWalkerAPES *apes, gboolean center_shrink)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->center_shrink = center_shrink;

  if (self->constructed)
  {
    /* No check here: the kernel type may still be set afterwards, as the Python
     * helpers do. The combination is rejected when the estimators are rebuilt and,
     * failing that, by NcmStatsDist when it prepares. */
    ncm_stats_dist_set_center_shrink (self->sd0, self->center_shrink);
    ncm_stats_dist_set_center_shrink (self->sd1, self->center_shrink);
    ncm_stats_dist_set_defensive_frac (self->sd0, self->defensive_frac);
    ncm_stats_dist_set_defensive_frac (self->sd1, self->defensive_frac);
    ncm_stats_dist_set_defensive_scale (self->sd0, self->defensive_scale);
    ncm_stats_dist_set_defensive_scale (self->sd1, self->defensive_scale);
    ncm_stats_dist_set_defensive_nu (self->sd0, self->defensive_nu);
    ncm_stats_dist_set_defensive_nu (self->sd1, self->defensive_nu);
    ncm_stats_dist_set_uniform_weights (self->sd0, self->uniform_weights);
    ncm_stats_dist_set_uniform_weights (self->sd1, self->uniform_weights);
  }
}

/**
 * ncm_fit_esmcmc_walker_apes_get_center_shrink:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: whether the posterior approximations use centre shrinkage.
 */
gboolean
ncm_fit_esmcmc_walker_apes_get_center_shrink (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->center_shrink;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_defensive_frac:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @frac: weight of the wide Student-t component in the proposal
 *
 * Sets #NcmFitESMCMCWalkerAPES:defensive-frac, forwarded to both estimators.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_defensive_frac (NcmFitESMCMCWalkerAPES *apes, const gdouble frac)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->defensive_frac = frac;

  if (self->constructed)
  {
    ncm_stats_dist_set_defensive_frac (self->sd0, frac);
    ncm_stats_dist_set_defensive_frac (self->sd1, frac);
  }
}

/**
 * ncm_fit_esmcmc_walker_apes_get_defensive_frac:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: #NcmFitESMCMCWalkerAPES:defensive-frac.
 */
gdouble
ncm_fit_esmcmc_walker_apes_get_defensive_frac (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->defensive_frac;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_defensive_scale:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @scale: covariance factor of the wide component
 *
 * Sets #NcmFitESMCMCWalkerAPES:defensive-scale, forwarded to both estimators.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_defensive_scale (NcmFitESMCMCWalkerAPES *apes, const gdouble scale)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->defensive_scale = scale;

  if (self->constructed)
  {
    ncm_stats_dist_set_defensive_scale (self->sd0, scale);
    ncm_stats_dist_set_defensive_scale (self->sd1, scale);
  }
}

/**
 * ncm_fit_esmcmc_walker_apes_get_defensive_scale:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: #NcmFitESMCMCWalkerAPES:defensive-scale.
 */
gdouble
ncm_fit_esmcmc_walker_apes_get_defensive_scale (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->defensive_scale;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_defensive_nu:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @nu: degrees of freedom of the wide component
 *
 * Sets #NcmFitESMCMCWalkerAPES:defensive-nu, forwarded to both estimators.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_defensive_nu (NcmFitESMCMCWalkerAPES *apes, const gdouble nu)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->defensive_nu = nu;

  if (self->constructed)
  {
    ncm_stats_dist_set_defensive_nu (self->sd0, nu);
    ncm_stats_dist_set_defensive_nu (self->sd1, nu);
  }
}

/**
 * ncm_fit_esmcmc_walker_apes_get_defensive_nu:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: #NcmFitESMCMCWalkerAPES:defensive-nu.
 */
gdouble
ncm_fit_esmcmc_walker_apes_get_defensive_nu (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->defensive_nu;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_vkde_points_per_dim:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @points_per_dim: see #NcmStatsDistVKDE:points-per-dim
 *
 * Sets #NcmFitESMCMCWalkerAPES:vkde-points-per-dim, forwarded to both estimators when the
 * method is VKDE.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_vkde_points_per_dim (NcmFitESMCMCWalkerAPES *apes, const gdouble points_per_dim)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->vkde_points_per_dim = points_per_dim;

  if (self->constructed && (self->method == NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE))
  {
    ncm_stats_dist_vkde_set_points_per_dim (NCM_STATS_DIST_VKDE (self->sd0), points_per_dim);
    ncm_stats_dist_vkde_set_points_per_dim (NCM_STATS_DIST_VKDE (self->sd1), points_per_dim);
    _ncm_fit_esmcmc_walker_apes_vkde_check_sizes (NCM_FIT_ESMCMC_WALKER (apes));
  }
}

/**
 * ncm_fit_esmcmc_walker_apes_get_vkde_points_per_dim:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: #NcmFitESMCMCWalkerAPES:vkde-points-per-dim.
 */
gdouble
ncm_fit_esmcmc_walker_apes_get_vkde_points_per_dim (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->vkde_points_per_dim;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_uniform_weights:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @uniform_weights: see #NcmStatsDist:uniform-weights
 *
 * Sets #NcmFitESMCMCWalkerAPES:uniform-weights, forwarded to both estimators.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_uniform_weights (NcmFitESMCMCWalkerAPES *apes, gboolean uniform_weights)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->uniform_weights = uniform_weights;

  if (self->constructed)
  {
    ncm_stats_dist_set_uniform_weights (self->sd0, uniform_weights);
    ncm_stats_dist_set_uniform_weights (self->sd1, uniform_weights);
  }
}

/**
 * ncm_fit_esmcmc_walker_apes_get_uniform_weights:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: #NcmFitESMCMCWalkerAPES:uniform-weights.
 */
gboolean
ncm_fit_esmcmc_walker_apes_get_uniform_weights (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->uniform_weights;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_cv_type:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @cv_type: a #NcmStatsDistCV
 *
 * Sets the cross-validation used to choose the over-smooth factor, see
 * #NcmFitESMCMCWalkerAPES:cv-type.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_cv_type (NcmFitESMCMCWalkerAPES *apes, NcmStatsDistCV cv_type)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->cv_type = cv_type;

  if (self->constructed)
  {
    ncm_stats_dist_set_cv_type (self->sd0, self->cv_type);
    ncm_stats_dist_set_cv_type (self->sd1, self->cv_type);
  }
}

/**
 * ncm_fit_esmcmc_walker_apes_get_cv_type:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: the cross-validation used to choose the over-smooth factor.
 */
NcmStatsDistCV
ncm_fit_esmcmc_walker_apes_get_cv_type (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->cv_type;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_split_frac:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @split_frac: the fraction of the block used as kernel centres
 *
 * Sets the fraction of the block used as kernel centres when the cross-validation
 * splits the sample, see #NcmFitESMCMCWalkerAPES:split-frac. Zero keeps the
 * #NcmStatsDist default.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_split_frac (NcmFitESMCMCWalkerAPES *apes, const gdouble split_frac)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if ((split_frac > 0.0) && (split_frac < 0.01))
    g_error ("ncm_fit_esmcmc_walker_apes_set_split_frac: invalid split fraction `%f', "
             "use zero to keep the NcmStatsDist default or a value in [0.01, 1].", split_frac);

  self->split_frac = split_frac;

  if (self->constructed && (split_frac > 0.0))
  {
    ncm_stats_dist_set_split_frac (self->sd0, split_frac);
    ncm_stats_dist_set_split_frac (self->sd1, split_frac);
  }
}

/**
 * ncm_fit_esmcmc_walker_apes_get_split_frac:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: the fraction of the block used as kernel centres.
 */
gdouble
ncm_fit_esmcmc_walker_apes_get_split_frac (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->split_frac;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_auto_kernel:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @auto_kernel: whether to choose the kernel with the bandwidth
 *
 * Sets whether the kernel is chosen together with the over-smooth factor, see
 * #NcmFitESMCMCWalkerAPES:auto-kernel.
 *
 * Deprecated: 0.24.0: set #NcmFitESMCMCWalkerAPES:kernel-type to
 * #NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO instead. Turning it off here leaves the Gaussian
 * kernel, since the previous choice is not recorded anywhere.
 */
void
ncm_fit_esmcmc_walker_apes_set_auto_kernel (NcmFitESMCMCWalkerAPES *apes, gboolean auto_kernel)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if (auto_kernel)
    ncm_fit_esmcmc_walker_apes_set_k_type (apes, NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO);
  else if (self->k_type == NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO)
    ncm_fit_esmcmc_walker_apes_set_k_type (apes, NCM_FIT_ESMCMC_WALKER_APES_KTYPE_GAUSS);
}

/**
 * ncm_fit_esmcmc_walker_apes_get_auto_kernel:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: whether the kernel is chosen together with the over-smooth factor.
 *
 * Deprecated: 0.24.0: read #NcmFitESMCMCWalkerAPES:kernel-type instead.
 */
gboolean
ncm_fit_esmcmc_walker_apes_get_auto_kernel (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->k_type == NCM_FIT_ESMCMC_WALKER_APES_KTYPE_AUTO;
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

  self->local_frac = local_frac;

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

  self->cov_type = NCM_STATS_DIST_KDE_COV_TYPE_FIXED;
  ncm_matrix_clear (&self->cov_fixed);
  self->cov_fixed = ncm_matrix_ref (cov_fixed);

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

  self->cov_type = NCM_STATS_DIST_KDE_COV_TYPE_ROBUST_DIAG;

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

  self->cov_type = NCM_STATS_DIST_KDE_COV_TYPE_ROBUST;

  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd0), NCM_STATS_DIST_KDE_COV_TYPE_ROBUST);
  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd1), NCM_STATS_DIST_KDE_COV_TYPE_ROBUST);
}

/**
 * ncm_fit_esmcmc_walker_apes_set_exploration:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @exploration: length cap of the exploration phase, in iterations (0: no cap)
 *
 * Sets #NcmFitESMCMCWalkerAPES:exploration. The rows of the phase are not part of the
 * Markovian chain; #NcmMSetCatalog:markovian-id records where it ended.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_exploration (NcmFitESMCMCWalkerAPES *apes, guint exploration)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->exploration = exploration;
}

/**
 * ncm_fit_esmcmc_walker_apes_get_exploration:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: #NcmFitESMCMCWalkerAPES:exploration.
 */
guint
ncm_fit_esmcmc_walker_apes_get_exploration (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->exploration;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_exploration_qratio_floor:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @qratio_floor: floor of $q(x)/q(x')$ during the exploration phase, in [0, 1]
 *
 * Sets #NcmFitESMCMCWalkerAPES:exploration-qratio-floor.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_exploration_qratio_floor (NcmFitESMCMCWalkerAPES *apes, const gdouble qratio_floor)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if ((qratio_floor < 0.0) || (qratio_floor > 1.0))
    g_error ("ncm_fit_esmcmc_walker_apes_set_exploration_qratio_floor: floor %g outside [0, 1].", qratio_floor);

  self->exploration_qratio_floor = qratio_floor;
}

/**
 * ncm_fit_esmcmc_walker_apes_get_exploration_qratio_floor:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: #NcmFitESMCMCWalkerAPES:exploration-qratio-floor.
 */
gdouble
ncm_fit_esmcmc_walker_apes_get_exploration_qratio_floor (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->exploration_qratio_floor;
}

/**
 * ncm_fit_esmcmc_walker_apes_set_exploration_stop_after:
 * @apes: a #NcmFitESMCMCWalkerAPES
 * @stop_after: consecutive iterations with no clipped acceptance that end the exploration phase (>= 1)
 *
 * Sets #NcmFitESMCMCWalkerAPES:exploration-stop-after.
 *
 */
void
ncm_fit_esmcmc_walker_apes_set_exploration_stop_after (NcmFitESMCMCWalkerAPES *apes, guint stop_after)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  if (stop_after == 0)
    g_error ("ncm_fit_esmcmc_walker_apes_set_exploration_stop_after: stop_after must be at least 1.");

  self->exploration_stop_after = stop_after;
}

/**
 * ncm_fit_esmcmc_walker_apes_get_exploration_stop_after:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: #NcmFitESMCMCWalkerAPES:exploration-stop-after.
 */
guint
ncm_fit_esmcmc_walker_apes_get_exploration_stop_after (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->exploration_stop_after;
}

/**
 * ncm_fit_esmcmc_walker_apes_is_exploring:
 * @apes: a #NcmFitESMCMCWalkerAPES
 *
 * Returns: whether the exploration phase is armed (between a fresh start of a run and the
 * end of the phase).
 */
gboolean
ncm_fit_esmcmc_walker_apes_is_exploring (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  return self->exploring;
}

