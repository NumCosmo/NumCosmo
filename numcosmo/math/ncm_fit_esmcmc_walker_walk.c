/***************************************************************************
 *            ncm_fit_esmcmc_walker_walk.c
 *
 *  Tue March 29 10:41:49 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_walk.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_fit_esmcmc_walker_walk
 * @title: NcmFitESMCMCWalkerWalk
 * @short_description: Ensemble sampler Markov Chain Monte Carlo walker - walk move.
 * @stability: Unstable
 *
 * Implementing walk move walker for #NcmFitESMCMC (affine invariant).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_esmcmc_walker.h"
#include "math/ncm_fit_esmcmc_walker_walk.h"

#include "math/ncm_fit_esmcmc.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_SCALE
};

struct _NcmFitESMCMCWalkerWalk
{
  /*< private >*/
  NcmFitESMCMCWalker parent_instance;
  guint size;
  guint size_2;
  guint nparams;
  gdouble a;
  gdouble sqrt_nparams;
  NcmMatrix *z;
  GPtrArray *thetabar;
  GArray *indices;
  GArray *numbers;
};

G_DEFINE_TYPE (NcmFitESMCMCWalkerWalk, ncm_fit_esmcmc_walker_walk, NCM_TYPE_FIT_ESMCMC_WALKER);

static void
ncm_fit_esmcmc_walker_walk_init (NcmFitESMCMCWalkerWalk *walk)
{
  walk->size         = 0;
  walk->size_2       = 0;
  walk->nparams      = 0;
  walk->a            = 0.0;
  walk->sqrt_nparams = 0.0;
  walk->z            = NULL;
  walk->thetabar     = g_ptr_array_new ();
  walk->indices      = g_array_new (TRUE, TRUE, sizeof (guint));
  walk->numbers      = g_array_new (TRUE, TRUE, sizeof (guint));

  g_ptr_array_set_free_func (walk->thetabar, (GDestroyNotify) ncm_vector_free);
}

static void
ncm_fit_esmcmc_walker_walk_dispose (GObject *object)
{
  NcmFitESMCMCWalkerWalk *walk = NCM_FIT_ESMCMC_WALKER_WALK (object);

  ncm_matrix_clear (&walk->z);

  g_clear_pointer (&walk->thetabar, g_ptr_array_unref);
  g_clear_pointer (&walk->indices, g_array_unref);
  g_clear_pointer (&walk->numbers, g_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_walk_parent_class)->dispose (object);
}

static void
ncm_fit_esmcmc_walker_walk_finalize (GObject *object)
{
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_walk_parent_class)->finalize (object);
}

static void
ncm_fit_esmcmc_walker_walk_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerWalk *walk = NCM_FIT_ESMCMC_WALKER_WALK (object);

  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_WALK (object));

  switch (prop_id)
  {
    case PROP_SCALE:
      ncm_fit_esmcmc_walker_walk_set_scale (walk, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_esmcmc_walker_walk_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerWalk *walk = NCM_FIT_ESMCMC_WALKER_WALK (object);

  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_WALK (object));

  switch (prop_id)
  {
    case PROP_SCALE:
      g_value_set_double (value, ncm_fit_esmcmc_walker_walk_get_scale (walk));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void _ncm_fit_esmcmc_walker_walk_set_size (NcmFitESMCMCWalker *walker, guint size);
static guint _ncm_fit_esmcmc_walker_walk_get_size (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_walk_set_nparams (NcmFitESMCMCWalker *walker, guint nparams);
static guint _ncm_fit_esmcmc_walker_walk_get_nparams (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_walk_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng);
static void _ncm_fit_esmcmc_walker_walk_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static gdouble _ncm_fit_esmcmc_walker_walk_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star);
static gdouble _ncm_fit_esmcmc_walker_walk_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static void _ncm_fit_esmcmc_walker_walk_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf);
static const gchar *_ncm_fit_esmcmc_walker_walk_desc (NcmFitESMCMCWalker *walker);

static void
ncm_fit_esmcmc_walker_walk_class_init (NcmFitESMCMCWalkerWalkClass *klass)
{
  GObjectClass *object_class            = G_OBJECT_CLASS (klass);
  NcmFitESMCMCWalkerClass *walker_class = NCM_FIT_ESMCMC_WALKER_CLASS (klass);

  object_class->set_property = ncm_fit_esmcmc_walker_walk_set_property;
  object_class->get_property = ncm_fit_esmcmc_walker_walk_get_property;
  object_class->dispose      = ncm_fit_esmcmc_walker_walk_dispose;
  object_class->finalize     = ncm_fit_esmcmc_walker_walk_finalize;

  g_object_class_install_property (object_class,
                                   PROP_SCALE,
                                   g_param_spec_double ("scale",
                                                        NULL,
                                                        "Walk scale a",
                                                        1.0e-2, G_MAXDOUBLE, 0.2,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  walker_class->set_size    = &_ncm_fit_esmcmc_walker_walk_set_size;
  walker_class->get_size    = &_ncm_fit_esmcmc_walker_walk_get_size;
  walker_class->set_nparams = &_ncm_fit_esmcmc_walker_walk_set_nparams;
  walker_class->get_nparams = &_ncm_fit_esmcmc_walker_walk_get_nparams;
  walker_class->setup       = &_ncm_fit_esmcmc_walker_walk_setup;
  walker_class->step        = &_ncm_fit_esmcmc_walker_walk_step;
  walker_class->prob        = &_ncm_fit_esmcmc_walker_walk_prob;
  walker_class->prob_norm   = &_ncm_fit_esmcmc_walker_walk_prob_norm;
  walker_class->clean       = &_ncm_fit_esmcmc_walker_walk_clean;
  walker_class->desc        = &_ncm_fit_esmcmc_walker_walk_desc;
}

static void
_ncm_fit_esmcmc_walker_walk_set_sys (NcmFitESMCMCWalker *walker, guint size, guint nparams)
{
  NcmFitESMCMCWalkerWalk *walk = NCM_FIT_ESMCMC_WALKER_WALK (walker);

  g_assert_cmpuint (size, >, 0);
  g_assert_cmpuint (nparams, >, 0);

  if ((walk->size != size) || (walk->nparams != nparams))
  {
    guint i;

    ncm_matrix_clear (&walk->z);

    g_assert (size % 2 == 0);

    walk->z = ncm_matrix_new (size, nparams);

    g_array_set_size (walk->indices, size * nparams);

    g_array_set_size (walk->numbers, size);

    for (i = 0; i < size; i++)
      g_array_index (walk->numbers, guint, i) = i;

    g_ptr_array_set_size (walk->thetabar, 0);

    for (i = 0; i < size; i++)
    {
      NcmVector *thetabar_i = ncm_vector_new (nparams);

      g_ptr_array_add (walk->thetabar, thetabar_i);
    }

    walk->size   = size;
    walk->size_2 = size / 2;

    walk->size         = size;
    walk->nparams      = nparams;
    walk->sqrt_nparams = sqrt (1.0 * walk->nparams);
  }
}

static void
_ncm_fit_esmcmc_walker_walk_set_size (NcmFitESMCMCWalker *walker, guint size)
{
  NcmFitESMCMCWalkerWalk *walk = NCM_FIT_ESMCMC_WALKER_WALK (walker);

  g_assert_cmpuint (size, >, 0);

  if (walk->nparams != 0)
    _ncm_fit_esmcmc_walker_walk_set_sys (walker, size, walk->nparams);
  else
    walk->size = size;
}

static guint
_ncm_fit_esmcmc_walker_walk_get_size (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerWalk *walk = NCM_FIT_ESMCMC_WALKER_WALK (walker);

  return walk->size;
}

static void
_ncm_fit_esmcmc_walker_walk_set_nparams (NcmFitESMCMCWalker *walker, guint nparams)
{
  NcmFitESMCMCWalkerWalk *walk = NCM_FIT_ESMCMC_WALKER_WALK (walker);

  g_assert_cmpuint (nparams, >, 0);

  if (walk->size != 0)
    _ncm_fit_esmcmc_walker_walk_set_sys (walker, walk->size, nparams);
  else
    walk->nparams = nparams;
}

static guint
_ncm_fit_esmcmc_walker_walk_get_nparams (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerWalk *walk = NCM_FIT_ESMCMC_WALKER_WALK (walker);

  return walk->nparams;
}

static void
_ncm_fit_esmcmc_walker_walk_setup (NcmFitESMCMCWalker *walker, NcmMSet *mset, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng)
{
  NcmFitESMCMCWalkerWalk *walk = NCM_FIT_ESMCMC_WALKER_WALK (walker);
  guint k;

  for (k = ki; k < kf; k++)
  {
    const guint subensemble = (k < walk->size_2) ? walk->size_2 : 0;
    guint pi;

    for (pi = 0; pi < walk->nparams; pi++)
    {
      const gdouble zi = gsl_ran_ugaussian (rng->r);

      ncm_matrix_set (walk->z, k, pi, zi);
    }

    gsl_ran_choose (rng->r,
                    &g_array_index (walk->indices, guint, k * walk->nparams),
                    walk->nparams,
                    &g_array_index (walk->numbers, guint, subensemble),
                    walk->size_2,
                    g_array_get_element_size (walk->numbers));

/*
 *   for (pi = 0; pi < walk->nparams; pi++)
 *   {
 *     printf ("Choosing for walker %u, walker %u | nparams %u size_2 %u subensemble %u | numbers[subensemble + %u] = %u\n",
 *             k, g_array_index (walk->indices, guint, k * walk->nparams + pi),
 *             walk->nparams, walk->size_2, subensemble,
 *             pi, g_array_index (walk->numbers, guint, subensemble + pi));
 *   }
 */
  }
}

static void
_ncm_fit_esmcmc_walker_walk_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerWalk *walk = NCM_FIT_ESMCMC_WALKER_WALK (walker);
  NcmVector *theta_k           = g_ptr_array_index (theta, k);
  NcmVector *thetabar_k        = g_ptr_array_index (walk->thetabar, k);
  guint i;

  ncm_vector_set_zero (thetabar_k);

  for (i = 0; i < walk->nparams; i++)
  {
    const guint j      = g_array_index (walk->indices, guint, k * walk->nparams + i);
    NcmVector *theta_j = g_ptr_array_index (theta, j);

    /*printf ("walker %u using walker %u to calculate mean\n", k, j);*/

    ncm_vector_add (thetabar_k, theta_j);
  }

  ncm_vector_scale (thetabar_k, 1.0 / (1.0 * walk->nparams));

  ncm_vector_memcpy (thetastar, theta_k);

  for (i = 0; i < walk->nparams; i++)
  {
    const gdouble z    = ncm_matrix_get (walk->z, k, i);
    const guint j      = g_array_index (walk->indices, guint, k * walk->nparams + i);
    NcmVector *theta_j = g_ptr_array_index (theta, j);
    guint m;

    for (m = 0; m < walk->nparams; m++)
    {
      const gdouble thetabar_k_i = ncm_vector_get (thetabar_k, i);
      const gdouble theta_j_i    = ncm_vector_get (theta_j, i);

      ncm_vector_addto (thetastar, i, walk->a * z * (theta_j_i - thetabar_k_i) / walk->sqrt_nparams);
    }
  }
}

static gdouble
_ncm_fit_esmcmc_walker_walk_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  return exp ((m2lnL_cur - m2lnL_star) * 0.5);
}

static gdouble
_ncm_fit_esmcmc_walker_walk_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  return 0.0;
}

static void
_ncm_fit_esmcmc_walker_walk_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  /* Nothing to do. */
}

const gchar *
_ncm_fit_esmcmc_walker_walk_desc (NcmFitESMCMCWalker *walker)
{
  return "Walk-Move";
}

/**
 * ncm_fit_esmcmc_walker_walk_new:
 * @nwalkers: number of walkers
 *
 * Creates a new #NcmFitESMCMCWalkerWalk to be used
 * with @nwalkers.
 *
 * Returns: (transfer full): a new #NcmFitESMCMCWalkerWalk.
 */
NcmFitESMCMCWalkerWalk *
ncm_fit_esmcmc_walker_walk_new (guint nwalkers)
{
  NcmFitESMCMCWalkerWalk *walk = g_object_new (NCM_TYPE_FIT_ESMCMC_WALKER_WALK,
                                               "size", nwalkers,
                                               NULL);

  return walk;
}

/**
 * ncm_fit_esmcmc_walker_walk_set_scale:
 * @walk: a #NcmFitESMCMCWalkerWalk
 * @a: new scale $a > 1$
 *
 * Sets the value of the scale $a > 1$.
 *
 */
void
ncm_fit_esmcmc_walker_walk_set_scale (NcmFitESMCMCWalkerWalk *walk, const gdouble a)
{
  g_assert_cmpfloat (a, >, 1.0e-2);
  walk->a = a;
}

/**
 * ncm_fit_esmcmc_walker_walk_get_scale:
 * @walk: a #NcmFitESMCMCWalkerWalk
 *
 * Gets the value of the scale $a > 1$.
 *
 * Returns: current value of $a$.
 */
gdouble
ncm_fit_esmcmc_walker_walk_get_scale (NcmFitESMCMCWalkerWalk *walk)
{
  return walk->a;
}

