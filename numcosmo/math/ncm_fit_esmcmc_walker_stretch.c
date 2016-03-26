/***************************************************************************
 *            ncm_fit_esmcmc_walker_stretch.c
 *
 *  Wed March 16 15:53:28 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * ncm_fit_esmcmc_walker_stretch.c
 * Copyright (C) 2016 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:ncm_fit_esmcmc_walker_stretch
 * @title: NcmFitESMCMCWalkerStretch
 * @short_description: Ensemble sampler Markov Chain Monte Carlo walker - stretch move.
 *
 * Implementing stretch move walker for #NcmFitESMCMC (affine invariant).
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_fit_esmcmc_walker.h"
#include "math/ncm_fit_esmcmc_walker_stretch.h"

#include "math/ncm_fit_esmcmc.h"

#include <gsl/gsl_math.h>

enum
{
  PROP_0,
  PROP_SCALE
};

G_DEFINE_TYPE (NcmFitESMCMCWalkerStretch, ncm_fit_esmcmc_walker_stretch, NCM_TYPE_FIT_ESMCMC_WALKER);

static void
ncm_fit_esmcmc_walker_stretch_init (NcmFitESMCMCWalkerStretch *stretch)
{
  stretch->size    = 0;
  stretch->size_2  = 0;
  stretch->a       = 0.0;
  stretch->z       = NULL;
  stretch->indices = g_array_new (TRUE, TRUE, sizeof (guint));
}

static void
ncm_fit_esmcmc_walker_stretch_dispose (GObject *object)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (object);

  ncm_vector_clear (&stretch->z);
  g_clear_pointer (&stretch->indices, g_array_unref);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_stretch_parent_class)->dispose (object);
}

static void
ncm_fit_esmcmc_walker_stretch_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_stretch_parent_class)->finalize (object);
}

static void
ncm_fit_esmcmc_walker_stretch_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (object);
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_STRETCH (object));

  switch (prop_id)
  {
    case PROP_SCALE:
      ncm_fit_esmcmc_walker_stretch_set_scale (stretch, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_fit_esmcmc_walker_stretch_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (object);
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_STRETCH (object));

  switch (prop_id)
  {
    case PROP_SCALE:
      g_value_set_double (value, ncm_fit_esmcmc_walker_stretch_get_scale (stretch));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void _ncm_fit_esmcmc_walker_stretch_set_size (NcmFitESMCMCWalker *walker, guint size);
static guint _ncm_fit_esmcmc_walker_stretch_get_size (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_stretch_setup (NcmFitESMCMCWalker *walker, GPtrArray *theta, guint ki, guint kf, NcmRNG *rng);
static void _ncm_fit_esmcmc_walker_stretch_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, NcmVector *thetastar, guint k);
static gdouble _ncm_fit_esmcmc_walker_stretch_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star);
static void _ncm_fit_esmcmc_walker_stretch_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf);
static const gchar *_ncm_fit_esmcmc_walker_stretch_desc (NcmFitESMCMCWalker *walker);

static void
ncm_fit_esmcmc_walker_stretch_class_init (NcmFitESMCMCWalkerStretchClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmFitESMCMCWalkerClass *walker_class = NCM_FIT_ESMCMC_WALKER_CLASS (klass);

  object_class->set_property = ncm_fit_esmcmc_walker_stretch_set_property;
  object_class->get_property = ncm_fit_esmcmc_walker_stretch_get_property;
  object_class->dispose      = ncm_fit_esmcmc_walker_stretch_dispose;
  object_class->finalize     = ncm_fit_esmcmc_walker_stretch_finalize;

  g_object_class_install_property (object_class,
                                   PROP_SCALE,
                                   g_param_spec_double ("scale",
                                                        NULL,
                                                        "Strech scale a",
                                                        1.1, G_MAXDOUBLE, 2.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  walker_class->set_size = &_ncm_fit_esmcmc_walker_stretch_set_size;
  walker_class->get_size = &_ncm_fit_esmcmc_walker_stretch_get_size;
  walker_class->setup    = &_ncm_fit_esmcmc_walker_stretch_setup;
  walker_class->step     = &_ncm_fit_esmcmc_walker_stretch_step;
  walker_class->prob     = &_ncm_fit_esmcmc_walker_stretch_prob;
  walker_class->clean    = &_ncm_fit_esmcmc_walker_stretch_clean;
  walker_class->desc     = &_ncm_fit_esmcmc_walker_stretch_desc;
}

static void 
_ncm_fit_esmcmc_walker_stretch_set_size (NcmFitESMCMCWalker *walker, guint size)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);
  if (stretch->size != size)
  {
    ncm_vector_clear (&stretch->z);

    g_assert (size % 2 == 0);

    if (size > 0)
      stretch->z = ncm_vector_new (size);

    g_array_set_size (stretch->indices, size);
    stretch->size   = size;
    stretch->size_2 = size / 2;
  }
}

static guint 
_ncm_fit_esmcmc_walker_stretch_get_size (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);

  return stretch->size;
}

static void 
_ncm_fit_esmcmc_walker_stretch_setup (NcmFitESMCMCWalker *walker, GPtrArray *theta, guint ki, guint kf, NcmRNG *rng)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);
  /*const guint len = ncm_vector_len (g_ptr_array_index (theta, 0));*/
  /*const gint nvar = len - NCM_FIT_ESMCMC_NADD_VALS;*/
  guint k;

  for (k = ki; k < kf; k++)
  {
    const guint subensemble = (k < stretch->size_2) ? stretch->size_2 : 0;
    const gdouble u         = gsl_rng_uniform (rng->r);
    const gdouble z         = gsl_pow_2 (1.0 + (stretch->a - 1.0) * u) / stretch->a;
    /*const gdouble z         = pow (1.0 + (pow (stretch->a, nvar) - 1.0) * u, 2.0 / nvar) / stretch->a;*/
    const guint j           = gsl_rng_uniform_int (rng->r, stretch->size_2) + subensemble;

    /*printf ("# Walker %u using z = % 20.15g and other walker %u to move.\n", k, z, j);*/
    
    ncm_vector_set (stretch->z, k, z);
    g_array_index (stretch->indices, guint, k) = j;
  }
}

static void 
_ncm_fit_esmcmc_walker_stretch_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);
  const guint j   = g_array_index (stretch->indices, guint, k);
  const gdouble z = ncm_vector_get (stretch->z, k);
  const guint len = ncm_vector_len (thetastar);
  NcmVector *theta_k = g_ptr_array_index (theta, k);
  NcmVector *theta_j = g_ptr_array_index (theta, j);
  guint i;

  for (i = NCM_FIT_ESMCMC_NADD_VALS; i < len; i++)
  {
    const gdouble theta_j_i = ncm_vector_get (theta_j, i);
    const gdouble theta_k_i = ncm_vector_get (theta_k, i);
    
    ncm_vector_set (thetastar, i, theta_j_i + z * (theta_k_i - theta_j_i));
  }
}

static gdouble 
_ncm_fit_esmcmc_walker_stretch_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);
  const guint len = ncm_vector_len (thetastar);
  const gdouble z = ncm_vector_get (stretch->z, k);
  const gint nvar = len - NCM_FIT_ESMCMC_NADD_VALS;

  return pow (z, nvar - 1.0) * exp ((m2lnL_cur - m2lnL_star) * 0.5);
  /*return exp ((m2lnL_cur - m2lnL_star) * 0.5);*/
}

static void 
_ncm_fit_esmcmc_walker_stretch_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  /* Nothing to do. */
}

const gchar *
_ncm_fit_esmcmc_walker_stretch_desc (NcmFitESMCMCWalker *walker)
{
  return "Stretch-Move";
}

/**
 * ncm_fit_esmcmc_walker_stretch_new:
 * @nwalkers: number of walkers
 * 
 * Creates a new #NcmFitESMCMCWalkerStretch to be used
 * with @nwalkers.
 * 
 * Returns: (transfer full): a new #NcmFitESMCMCWalkerStretch.
 */
NcmFitESMCMCWalkerStretch *
ncm_fit_esmcmc_walker_stretch_new (guint nwalkers)
{
  NcmFitESMCMCWalkerStretch *stretch = g_object_new (NCM_TYPE_FIT_ESMCMC_WALKER_STRETCH,
                                                     "size", nwalkers,
                                                     NULL);

  return stretch;
}

/**
 * ncm_fit_esmcmc_walker_stretch_set_scale:
 * @stretch: a #NcmFitESMCMCWalkerStretch
 * @a: new scale $a > 1$
 * 
 * Sets the value of the scale $a > 1$.
 * 
 */
void
ncm_fit_esmcmc_walker_stretch_set_scale (NcmFitESMCMCWalkerStretch *stretch, const gdouble a)
{
  g_assert_cmpfloat (a, >, 1.0);
  stretch->a = a;
}

/**
 * ncm_fit_esmcmc_walker_stretch_get_scale:
 * @stretch: a #NcmFitESMCMCWalkerStretch
 * 
 * Gets the value of the scale $a > 1$.
 * 
 * Returns: current value of $a$.
 */
gdouble
ncm_fit_esmcmc_walker_stretch_get_scale (NcmFitESMCMCWalkerStretch *stretch)
{
  return stretch->a;
}
