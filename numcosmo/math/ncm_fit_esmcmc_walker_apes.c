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
 * SECTION:ncm_fit_esmcmc_walker_apes
 * @title: NcmFitESMCMCWalkerAPES
 * @short_description: Ensemble sampler Markov Chain Monte Carlo walker - apes move.
 *
 * Implementing apes move walker for #NcmFitESMCMC (affine invariant).
 *
 * This object implements the Approximate Posterior Ensemble Sample (APES) step proposal for a walker.
 * This proposal was developed by Sandro Dias Pinto Vitenti and implemented in this library. Below there is a description of the proposal.
 *
 * The APES proposal consists of using radial basis interpolation to generate an interpolant $\tilde{\pi}$
 * from a target distribution $\pi$ and use this interpolant to propose new points for the walker.
 * By using a distribution $\tilde{\pi}$ that resembles the original target distribution, the APES proposal
 * generates samples that converge faster to the target distribution and are more independent when compared to other step proposals.
 *
 * The APES step is implemented as follows: suppose that there are $L$ walkers. They are divided into two blocks $L_1$ and $L_2$,
 * containing the first and the second half of the walkers respectively. When proposing new points $Y$ for the walkers in the $L_1$ block,
 * we use the points in the $L_2$ block to generate an interpolant $\tilde{\pi}_{L_2}$ and then propose points $Y \sim \tilde{\pi}_{L_2}$
 * for the $L_1$ block. These points are accepted or rejected based on an acceptance probability $A(Y|X)$, and after the points of the first
 * block are updated, we do the same procedure for the $L_2$ block using the $L_1$ block. This procedure can be seen in the pseudocode below.
 *
 * ![apes_sketch](apes.png)
 *
 * The user must provide the input the values: @nwalkers, @nparams, @method, @k\_@type, @over\_@smooth$ and @use\_@interp - ncm\_fit\_esmcmc\_walker\_apes\_new\_full().
 * The user can also initialize the object with: @nwalkers, @nparams - ncm\_fit\_esmcmc\_walker\_apes\_new() and let the remaining parameters as default,
 * which are defined in the properties of the class.
 * For more information about the algorithm, check the explanation below.
 *
 *		- This object shall be used in the #NcmFitESMCMC class to generate a Monte Carlo Markov Chain using an ensemble sampler.
 *                To see an example of its implementation, check the file example\_rosenbrock.py in NumCosmo/examples.
 *
 *		- Regarding the radial basis interpolation method is implemented, check the #NcmStatsDist class.
 *
 *		- Regarding the types of kernel used in the interpolation method as the radial basis function, check the #NcmStatsDistKernel class.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_esmcmc_walker.h"
#include "math/ncm_fit_esmcmc_walker_apes.h"

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
  PROP_USE_INTERP,
};

struct _NcmFitESMCMCWalkerAPESPrivate
{
  guint size;
  guint size_2;
  guint nparams;
  guint a_size;
  guint a_nparams;
  gint mk;
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
  gboolean constructed;
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmFitESMCMCWalkerAPES, ncm_fit_esmcmc_walker_apes, NCM_TYPE_FIT_ESMCMC_WALKER);

#define __MK(method, k_type) (method + (k_type << 8))

static void
ncm_fit_esmcmc_walker_apes_init (NcmFitESMCMCWalkerAPES *apes)
{
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv = ncm_fit_esmcmc_walker_apes_get_instance_private (apes);

  self->size        = 0;
  self->size_2      = 0;
  self->nparams     = 0;
  self->a_size      = 0;
  self->a_nparams   = 0;
  self->mk          = -1;
  self->m2lnp_star  = NULL;
  self->m2lnp_cur   = NULL;
  self->desc        = NULL;
  self->sd0         = NULL;
  self->sd1         = NULL;
  self->thetastar   = g_ptr_array_new ();
  self->m2lnL_s0    = NULL;
  self->m2lnL_s1    = NULL;
  self->method      = NCM_FIT_ESMCMC_WALKER_APES_METHOD_LEN;
  self->k_type      = NCM_FIT_ESMCMC_WALKER_APES_KTYPE_LEN;
  self->over_smooth = 0.0;
  self->use_interp  = FALSE;
  self->constructed = FALSE;

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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
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
    NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

    self->constructed = TRUE;
    _ncm_fit_esmcmc_walker_apes_set_sys (NCM_FIT_ESMCMC_WALKER (object));
  }
}

static void
_ncm_fit_esmcmc_walker_apes_dispose (GObject *object)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (object);
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

  ncm_vector_clear (&self->m2lnp_star);
  ncm_vector_clear (&self->m2lnp_cur);
  ncm_vector_clear (&self->m2lnL_s0);
  ncm_vector_clear (&self->m2lnL_s1);

  ncm_stats_dist_clear (&self->sd0);
  ncm_stats_dist_clear (&self->sd1);

  g_clear_pointer (&self->thetastar, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_apes_parent_class)->dispose (object);
}

static void
_ncm_fit_esmcmc_walker_apes_finalize (GObject *object)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (object);
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

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

  g_object_class_install_property (object_class,
                                   PROP_METHOD,
                                   g_param_spec_enum ("method",
                                                      NULL,
                                                      "Method used in posterior approximation",
                                                      NCM_TYPE_FIT_ESMCMC_WALKER_APES_METHOD, NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_K_TYPE,
                                   g_param_spec_enum ("kernel-type",
                                                      NULL,
                                                      "Kernel used in posterior approximation",
                                                      NCM_TYPE_FIT_ESMCMC_WALKER_APES_KTYPE, NCM_FIT_ESMCMC_WALKER_APES_KTYPE_CAUCHY,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_OVER_SMOOTH,
                                   g_param_spec_double ("over-smooth",
                                                        NULL,
                                                        "Over-smooth parameter used to adjust kernel bandwidth",
                                                        1.0e-10, 1.0e10, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  g_object_class_install_property (object_class,
                                   PROP_USE_INTERP,
                                   g_param_spec_boolean ("use-interp",
                                                         NULL,
                                                         "Whether to use interpolation to build the posterior approximation",
                                                         TRUE,
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
_ncm_fit_esmcmc_walker_apes_set_sys (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

  if ((self->size != self->a_size) ||
      (self->nparams != self->a_nparams) ||
      (self->mk != __MK (self->method, self->k_type)))
  {
    gint i;

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
          self->sd0 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, NCM_STATS_DIST_CV_NONE));
          self->sd1 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, NCM_STATS_DIST_CV_NONE));
          break;
        case NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE:
          self->sd0 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, NCM_STATS_DIST_CV_NONE));
          self->sd1 = NCM_STATS_DIST (ncm_stats_dist_vkde_new (kernel, NCM_STATS_DIST_CV_NONE));
          break;
        default:
          g_assert_not_reached ();
          break;
      }

      ncm_stats_dist_kernel_free (kernel);
    }

    ncm_stats_dist_set_over_smooth (self->sd0, self->over_smooth);
    ncm_stats_dist_set_over_smooth (self->sd1, self->over_smooth);

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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

  g_assert_cmpuint (size, >, 0);
  self->size = size;

  if (self->constructed)
    _ncm_fit_esmcmc_walker_apes_set_sys (walker);
}

static guint
_ncm_fit_esmcmc_walker_apes_get_size (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

  return self->size;
}

static void
_ncm_fit_esmcmc_walker_apes_set_nparams (NcmFitESMCMCWalker *walker, guint nparams)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

  g_assert_cmpuint (nparams, >, 0);
  self->nparams = nparams;

  if (self->constructed)
    _ncm_fit_esmcmc_walker_apes_set_sys (walker);
}

static guint
_ncm_fit_esmcmc_walker_apes_get_nparams (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

  return self->nparams;
}

static void
_ncm_fit_esmcmc_walker_apes_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;
  const gdouble T                            = 1.0;
  gint i;

  if (ki < self->size_2)
  {
    ncm_stats_dist_reset (self->sd0);

    for (i = self->size_2; i < self->size; i++)
    {
      /*NcmVector *theta_i = g_ptr_array_index (theta, i);*/
      ncm_vector_set (self->m2lnL_s0, i - self->size_2, (1.0 / T) * ncm_vector_get (g_ptr_array_index (m2lnL, i), 0));
      /*ncm_stats_dist_add_obs (NCM_STATS_DIST_ND_VBK (self->dndg0), theta_i);*/
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

    for (i = ki; i < self->size_2; i++)
    {
      NcmVector *thetastar_i = g_ptr_array_index (self->thetastar, i);

      do {
        ncm_stats_dist_sample (self->sd0, thetastar_i, rng);
      } while (!ncm_mset_fparam_valid_bounds (mset, thetastar_i));

      /*ncm_vector_log_vals (thetastar_i, "TS: ", "%12.5g", TRUE);*/
    }
  }

  if (kf > self->size_2)
  {
    ncm_stats_dist_reset (self->sd1);

    for (i = 0; i < self->size_2; i++)
    {
      /*NcmVector *theta_i = g_ptr_array_index (theta, i);*/
      ncm_vector_set (self->m2lnL_s1, i, (1.0 / T) * ncm_vector_get (g_ptr_array_index (m2lnL, i), 0));
      /*ncm_stats_dist_add_obs (self->dndg1, theta_i);*/
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

    for (i = self->size_2; i < kf; i++)
    {
      NcmVector *thetastar_i = g_ptr_array_index (self->thetastar, i);

      do {
        ncm_stats_dist_sample (self->sd1, thetastar_i, rng);
      } while (!ncm_mset_fparam_valid_bounds (mset, thetastar_i));

      /*ncm_vector_log_vals (thetastar_i, "TS: ", "%12.5g", TRUE);*/
    }
  }
}

static void
_ncm_fit_esmcmc_walker_apes_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;
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

  /*ncm_vector_log_vals (theta_k,   "    THETA: ", "% 22.15g", TRUE);*/
  /*ncm_vector_log_vals (thetastar, "THETASTAR: ", "% 22.15g", TRUE);*/
}

static gdouble
_ncm_fit_esmcmc_walker_apes_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;
  const gdouble m2lnp_star                   = ncm_vector_get (self->m2lnp_star, k);
  const gdouble m2lnp_cur                    = ncm_vector_get (self->m2lnp_cur, k);

/*
 *  if ((fabs (m2lnL_star) > 1.0e5) || (fabs (m2lnL_cur) > 1.0e5))
 *  {
 *   printf ("AAA m2lnL_star % 22.15g m2lnp_star % 22.15g m2lnL_cur % 22.15g m2lnp_cur % 22.15g L cur->star: %12.5g p star->cur: %12.5g | T %12.5g\n",
 *       m2lnL_star, m2lnp_star, m2lnL_cur, m2lnp_cur,
 *       exp (- 0.5 * (m2lnL_star - m2lnL_cur)), exp (- 0.5 * (m2lnp_cur - m2lnp_star)),
 *       MIN (exp (- 0.5 * ((m2lnL_star - m2lnp_star) - (m2lnL_cur - m2lnp_cur))), 1.0));
 *  }
 */
/*
 *  ncm_message ("AAA Dchi % 12.5f (% 12.5f - % 12.5f) DAchi % 12.5f (% 12.5f - % 12.5f) L cur->star: %12.5g p star->cur: %12.5g | ",
 *     m2lnL_star - m2lnL_cur, m2lnL_star, m2lnL_cur, m2lnp_star - m2lnp_cur, m2lnp_star, m2lnp_cur,
 *     exp (- 0.5 * (m2lnL_star - m2lnL_cur)), exp (- 0.5 * (m2lnp_cur - m2lnp_star)));
 */

  return exp (-0.5 * ((m2lnL_star - m2lnp_star) - (m2lnL_cur - m2lnp_cur)));
}

static gdouble
_ncm_fit_esmcmc_walker_apes_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;
  const gdouble m2lnp_star                   = ncm_vector_get (self->m2lnp_star, k);
  const gdouble m2lnp_cur                    = ncm_vector_get (self->m2lnp_cur, k);

  return -0.5 * (m2lnp_cur - m2lnp_star);
}

static void
_ncm_fit_esmcmc_walker_apes_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  /*NcmFitESMCMCWalkerAPES *apes = NCM_FIT_ESMCMC_WALKER_APES (walker);*/
  /*NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;*/

  /* Nothing to do. */
}

const gchar *
_ncm_fit_esmcmc_walker_apes_desc (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerAPES *apes               = NCM_FIT_ESMCMC_WALKER_APES (walker);
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

  return self->use_interp;
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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

  if (self->method != NCM_FIT_ESMCMC_WALKER_APES_METHOD_VKDE)
    g_error ("ncm_fit_esmcmc_walker_apes_set_local_frac: cannot set local fraction for a non-VKDE method.");

  ncm_stats_dist_vkde_set_local_frac (NCM_STATS_DIST_VKDE (self->sd0), local_frac);
  ncm_stats_dist_vkde_set_local_frac (NCM_STATS_DIST_VKDE (self->sd1), local_frac);
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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;
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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

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
  NcmFitESMCMCWalkerAPESPrivate * const self = apes->priv;

  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd0), NCM_STATS_DIST_KDE_COV_TYPE_ROBUST);
  ncm_stats_dist_kde_set_cov_type (NCM_STATS_DIST_KDE (self->sd1), NCM_STATS_DIST_KDE_COV_TYPE_ROBUST);
}

