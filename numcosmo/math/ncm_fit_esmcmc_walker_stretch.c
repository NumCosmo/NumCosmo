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
 * @stability: Stable
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

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#endif /* NUMCOSMO_GIR_SCAN */

enum
{
  PROP_0,
  PROP_SCALE,
  PROP_MULTI
};

G_DEFINE_TYPE (NcmFitESMCMCWalkerStretch, ncm_fit_esmcmc_walker_stretch, NCM_TYPE_FIT_ESMCMC_WALKER);

static void
ncm_fit_esmcmc_walker_stretch_init (NcmFitESMCMCWalkerStretch *stretch)
{
  stretch->size     = 0;
  stretch->size_2   = 0;
  stretch->nparams  = 0;
  stretch->a        = 0.0;
  stretch->z        = NULL;
  stretch->box      = NULL;
  stretch->norm_box = NULL;
  stretch->use_box  = g_array_new (TRUE, TRUE, sizeof (gboolean));
  stretch->indices  = g_array_new (TRUE, TRUE, sizeof (guint));
  stretch->numbers  = g_array_new (TRUE, TRUE, sizeof (guint));
  stretch->multi    = FALSE;
  stretch->desc     = NULL;
}

static void
_ncm_fit_esmcmc_walker_stretch_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (object);
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_STRETCH (object));

  switch (prop_id)
  {
    case PROP_SCALE:
      ncm_fit_esmcmc_walker_stretch_set_scale (stretch, g_value_get_double (value));
      break;
    case PROP_MULTI:
      ncm_fit_esmcmc_walker_stretch_multi (stretch, g_value_get_boolean (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_esmcmc_walker_stretch_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (object);
  g_return_if_fail (NCM_IS_FIT_ESMCMC_WALKER_STRETCH (object));

  switch (prop_id)
  {
    case PROP_SCALE:
      g_value_set_double (value, ncm_fit_esmcmc_walker_stretch_get_scale (stretch));
      break;
    case PROP_MULTI:
      g_value_set_boolean (value, stretch->multi);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_fit_esmcmc_walker_stretch_dispose (GObject *object)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (object);

  ncm_vector_clear (&stretch->norm_box);

  ncm_matrix_clear (&stretch->z);
  ncm_matrix_clear (&stretch->box);

  g_clear_pointer (&stretch->use_box, g_array_unref);
  g_clear_pointer (&stretch->indices, g_array_unref);
  g_clear_pointer (&stretch->numbers, g_array_unref);

  g_clear_pointer (&stretch->desc, g_free);
  
  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_stretch_parent_class)->dispose (object);
}

static void
_ncm_fit_esmcmc_walker_stretch_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_fit_esmcmc_walker_stretch_parent_class)->finalize (object);
}

static void _ncm_fit_esmcmc_walker_stretch_set_size (NcmFitESMCMCWalker *walker, guint size);
static guint _ncm_fit_esmcmc_walker_stretch_get_size (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_stretch_set_nparams (NcmFitESMCMCWalker *walker, guint nparams);
static guint _ncm_fit_esmcmc_walker_stretch_get_nparams (NcmFitESMCMCWalker *walker);
static void _ncm_fit_esmcmc_walker_stretch_setup (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng);
static void _ncm_fit_esmcmc_walker_stretch_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static gdouble _ncm_fit_esmcmc_walker_stretch_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star);
static gdouble _ncm_fit_esmcmc_walker_stretch_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k);
static void _ncm_fit_esmcmc_walker_stretch_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf);
static const gchar *_ncm_fit_esmcmc_walker_stretch_desc (NcmFitESMCMCWalker *walker);

static void
ncm_fit_esmcmc_walker_stretch_class_init (NcmFitESMCMCWalkerStretchClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmFitESMCMCWalkerClass *walker_class = NCM_FIT_ESMCMC_WALKER_CLASS (klass);

  object_class->set_property = &_ncm_fit_esmcmc_walker_stretch_set_property;
  object_class->get_property = &_ncm_fit_esmcmc_walker_stretch_get_property;
  object_class->dispose      = &_ncm_fit_esmcmc_walker_stretch_dispose;
  object_class->finalize     = &_ncm_fit_esmcmc_walker_stretch_finalize;

  g_object_class_install_property (object_class,
                                   PROP_SCALE,
                                   g_param_spec_double ("scale",
                                                        NULL,
                                                        "Strech scale a",
                                                        1.1, G_MAXDOUBLE, 2.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_MULTI,
                                   g_param_spec_boolean ("multi-stretch",
                                                         NULL,
                                                         "Whether it should use multiple stretchs per step",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  
  walker_class->set_size    = &_ncm_fit_esmcmc_walker_stretch_set_size;
  walker_class->get_size    = &_ncm_fit_esmcmc_walker_stretch_get_size;
  walker_class->set_nparams = &_ncm_fit_esmcmc_walker_stretch_set_nparams;
  walker_class->get_nparams = &_ncm_fit_esmcmc_walker_stretch_get_nparams;
  walker_class->setup       = &_ncm_fit_esmcmc_walker_stretch_setup;
  walker_class->step        = &_ncm_fit_esmcmc_walker_stretch_step;
  walker_class->prob        = &_ncm_fit_esmcmc_walker_stretch_prob;
  walker_class->prob_norm   = &_ncm_fit_esmcmc_walker_stretch_prob_norm;
  walker_class->clean       = &_ncm_fit_esmcmc_walker_stretch_clean;
  walker_class->desc        = &_ncm_fit_esmcmc_walker_stretch_desc;
}

static void 
_ncm_fit_esmcmc_walker_stretch_set_sys (NcmFitESMCMCWalker *walker, guint size, guint nparams)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);
  g_assert_cmpuint (size, >, 0);
  g_assert_cmpuint (nparams, >, 0);
  
  if (stretch->size != size || stretch->nparams != nparams)
  {
    guint i;

    ncm_vector_clear (&stretch->norm_box);
    ncm_matrix_clear (&stretch->z);
    ncm_matrix_clear (&stretch->box);

    g_assert (size % 2 == 0);

    stretch->z        = ncm_matrix_new (size, nparams);

    stretch->norm_box = ncm_vector_new (size);
    stretch->box      = ncm_matrix_new (nparams, 2);

    g_array_set_size (stretch->use_box, nparams);

    g_array_set_size (stretch->indices, size * nparams);

    g_array_set_size (stretch->numbers, size);
    for (i = 0; i < size; i++)
      g_array_index (stretch->numbers, guint, i) = i;
    
    stretch->size    = size;
    stretch->size_2  = size / 2;
    stretch->nparams = nparams;
  }
}

static void 
_ncm_fit_esmcmc_walker_stretch_set_size (NcmFitESMCMCWalker *walker, guint size)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);

  g_assert_cmpuint (size, >, 0);
  
  if (stretch->nparams != 0)
    _ncm_fit_esmcmc_walker_stretch_set_sys (walker, size, stretch->nparams);
  else
    stretch->size = size;  
}

static guint 
_ncm_fit_esmcmc_walker_stretch_get_size (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);

  return stretch->size;
}

static void 
_ncm_fit_esmcmc_walker_stretch_set_nparams (NcmFitESMCMCWalker *walker, guint nparams)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);

  g_assert_cmpuint (nparams, >, 0);
  
  if (stretch->size != 0)
    _ncm_fit_esmcmc_walker_stretch_set_sys (walker, stretch->size, nparams);
  else
    stretch->nparams = nparams;
}

static guint 
_ncm_fit_esmcmc_walker_stretch_get_nparams (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);

  return stretch->nparams;
}

static void 
_ncm_fit_esmcmc_walker_stretch_setup (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, guint ki, guint kf, NcmRNG *rng)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);
  guint k;

  if (!stretch->multi)
  {
    for (k = ki; k < kf; k++)
    {
      const guint subensemble = (k < stretch->size_2) ? stretch->size_2 : 0;
      const gdouble u         = gsl_rng_uniform (rng->r);
      const gdouble z         = gsl_pow_2 (1.0 + (stretch->a - 1.0) * u) / stretch->a;
      const guint j           = gsl_rng_uniform_int (rng->r, stretch->size_2) + subensemble;

      /*printf ("# Walker %u using z = % 20.15g and other walker %u to move.\n", k, z, j);*/

      ncm_matrix_set (stretch->z, k, 0, z);
      g_array_index (stretch->indices, guint, k) = j;
    }
  }
  else
  {
    for (k = ki; k < kf; k++)
    {
      const guint subensemble = (k < stretch->size_2) ? stretch->size_2 : 0;
      guint pi;
      for (pi = 0; pi < stretch->nparams; pi++)
      {
        const gdouble u = gsl_rng_uniform (rng->r);
        const gdouble z = gsl_pow_2 (1.0 + (stretch->a - 1.0) * u) / stretch->a;
        /*const guint j   = gsl_rng_uniform_int (rng->r, stretch->size_2) + subensemble;*/

        ncm_matrix_set (stretch->z, k, pi, z);
        /*g_array_index (stretch->indices, guint, k * stretch->nparams + pi) = j;*/
      }
      gsl_ran_choose (rng->r, 
                      &g_array_index (stretch->indices, guint, k * stretch->nparams), 
                      stretch->nparams, 
                      &g_array_index (stretch->numbers, guint, subensemble), 
                      stretch->size_2,
                      g_array_get_element_size (stretch->numbers));
    }
  }
}

static gdouble 
_ncm_fit_esmcmc_walker_stretch_theta_to_x (NcmFitESMCMCWalkerStretch *stretch, guint i, const gdouble theta_i)
{
  const gdouble lb          = ncm_matrix_get (stretch->box, i, 0);
  const gdouble ub          = ncm_matrix_get (stretch->box, i, 1);
  const gdouble bsize       = ub - lb;
  const gdouble tb          = theta_i - lb;
  const gdouble maxatan     = - 0.5 * log (GSL_DBL_EPSILON * 0.5);
  const gdouble twotb_bsize = 2.0 * tb / bsize;

  /*printf ("[% 22.15g % 22.15g]<% 22.15g> % 22.15g % 22.15g % 22.15g % 22.15g\n", lb, ub, theta_i, tb, maxatan, bsize, atanh (twotb_bsize - 1.0));*/
  
  if (twotb_bsize < GSL_DBL_EPSILON)             
    return -maxatan;
  else if (tb >= bsize)
    return maxatan;
  else
    return atanh (twotb_bsize - 1.0);
}

static gdouble 
_ncm_fit_esmcmc_walker_stretch_x_to_theta (NcmFitESMCMCWalkerStretch *stretch, guint i, const gdouble x_i)
{
  const gdouble lb     = ncm_matrix_get (stretch->box, i, 0);
  const gdouble ub     = ncm_matrix_get (stretch->box, i, 1);
  const gdouble bsize  = ub - lb;
  const gdouble middle = 0.5 * (ub + lb);

  return middle + bsize * tanh (x_i) * 0.5;
}

static gdouble 
_ncm_fit_esmcmc_walker_stretch_norm (NcmFitESMCMCWalkerStretch *stretch, guint i, const gdouble thetastar_i, const gdouble theta_i)
{
  const gdouble lb    = ncm_matrix_get (stretch->box, i, 0);
  const gdouble ub    = ncm_matrix_get (stretch->box, i, 1);
  const gdouble tb    = theta_i - lb;
  const gdouble bt    = ub - theta_i;
  const gdouble tsb   = thetastar_i - lb;
  const gdouble bts   = ub - thetastar_i;
  const gdouble norm  = tsb * bts / (tb * bt);
  
  return log (norm);
}

static void
_ncm_fit_esmcmc_walker_stretch_one_step (NcmFitESMCMCWalkerStretch *stretch, NcmVector *theta_k, NcmVector *theta_j, const gdouble z, NcmVector *thetastar)
{
  guint i;
  
  for (i = 0; i < stretch->nparams; i++)
  {
    const gdouble theta_k_i   = ncm_vector_get (theta_k, i);
    const gdouble theta_j_i   = ncm_vector_get (theta_j, i);

    if (!g_array_index (stretch->use_box, gboolean, i))
    {
      const gdouble thetastar_i = theta_j_i + z * (theta_k_i - theta_j_i);
      
      ncm_vector_set (thetastar, i, thetastar_i);
    }
    else
    {
      const gdouble x_k_i       = _ncm_fit_esmcmc_walker_stretch_theta_to_x (stretch, i, theta_k_i);
      const gdouble x_j_i       = _ncm_fit_esmcmc_walker_stretch_theta_to_x (stretch, i, theta_j_i);
      const gdouble xstar_i     = x_j_i + z * (x_k_i - x_j_i);
      const gdouble thetastar_i = _ncm_fit_esmcmc_walker_stretch_x_to_theta (stretch, i, xstar_i);
      
      ncm_vector_set (thetastar, i, thetastar_i);
    }
  }
}

static void
_ncm_fit_esmcmc_walker_stretch_step (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);
  NcmVector *theta_k = g_ptr_array_index (theta, k);
  guint i;

  ncm_vector_set (stretch->norm_box, k, 0.0);
  
  if (!stretch->multi)
  {
    const guint j      = g_array_index (stretch->indices, guint, k);
    const gdouble z    = ncm_matrix_get (stretch->z, k, 0);
    NcmVector *theta_j = g_ptr_array_index (theta, j);
    
    _ncm_fit_esmcmc_walker_stretch_one_step (stretch, theta_k, theta_j, z, thetastar);
  }
  else
  {
    NcmVector *theta_c = theta_k;
    guint si;

    /*ncm_vector_log_vals (theta_c, "before = ", "% 20.15g", TRUE);*/
    for (si = 0; si < stretch->nparams; si++)
    {
      const guint j      = g_array_index (stretch->indices, guint, k * stretch->nparams + si);
      const gdouble z    = ncm_matrix_get (stretch->z, k, si);
      NcmVector *theta_j = g_ptr_array_index (theta, j);
      /*printf ("Stretching k %u using j %u and z % 20.15g\n", k, j, z);*/
      _ncm_fit_esmcmc_walker_stretch_one_step (stretch, theta_c, theta_j, z, thetastar);
      theta_c = thetastar;
      /*ncm_vector_log_vals (theta_c, "during = ", "% 20.15g", TRUE);*/
    }
  }

  for (i = 0; i < stretch->nparams; i++)
  {
    if (g_array_index (stretch->use_box, gboolean, i))
    {
      const gdouble thetastar_i = ncm_vector_get (thetastar, i);
      const gdouble theta_k_i   = ncm_vector_get (theta_k, i);

      ncm_vector_addto (stretch->norm_box, k, 
                        _ncm_fit_esmcmc_walker_stretch_norm (stretch, i, thetastar_i, theta_k_i));
    }
  }
}

static gdouble 
_ncm_fit_esmcmc_walker_stretch_prob (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k, const gdouble m2lnL_cur, const gdouble m2lnL_star)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);
  
  if (!stretch->multi)
  {
    const gdouble z = ncm_matrix_get (stretch->z, k, 0);
    return pow (z, stretch->nparams - 1.0) * exp ((m2lnL_cur - m2lnL_star) * 0.5 + ncm_vector_get (stretch->norm_box, k));
  }
  else
  {
    return exp ((m2lnL_cur - m2lnL_star) * 0.5 + ncm_vector_get (stretch->norm_box, k));
  }
}

static gdouble 
_ncm_fit_esmcmc_walker_stretch_prob_norm (NcmFitESMCMCWalker *walker, GPtrArray *theta, GPtrArray *m2lnL, NcmVector *thetastar, guint k)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);
  
  if (!stretch->multi)
  {
    const gdouble z = ncm_matrix_get (stretch->z, k, 0);
    return (stretch->nparams - 1.0) * log (z) + ncm_vector_get (stretch->norm_box, k);
  }
  else
  {
    return ncm_vector_get (stretch->norm_box, k);
  }
}

static void 
_ncm_fit_esmcmc_walker_stretch_clean (NcmFitESMCMCWalker *walker, guint ki, guint kf)
{
  /* Nothing to do. */
}

const gchar *
_ncm_fit_esmcmc_walker_stretch_desc (NcmFitESMCMCWalker *walker)
{
  NcmFitESMCMCWalkerStretch *stretch = NCM_FIT_ESMCMC_WALKER_STRETCH (walker);

  g_clear_pointer (&stretch->desc, g_free);

  stretch->desc = g_strdup_printf ("Stretch-Move%s", stretch->multi ? "[multi-strecth]" : "");
  
  return stretch->desc;
}

/**
 * ncm_fit_esmcmc_walker_stretch_new:
 * @nwalkers: number of walkers
 * @nparams: number of parameters
 * 
 * Creates a new #NcmFitESMCMCWalkerStretch to be used
 * with @nwalkers.
 * 
 * Returns: (transfer full): a new #NcmFitESMCMCWalkerStretch.
 */
NcmFitESMCMCWalkerStretch *
ncm_fit_esmcmc_walker_stretch_new (guint nwalkers, guint nparams)
{
  NcmFitESMCMCWalkerStretch *stretch = g_object_new (NCM_TYPE_FIT_ESMCMC_WALKER_STRETCH,
                                                     "size", nwalkers,
                                                     "nparams", nparams,
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

/**
 * ncm_fit_esmcmc_walker_stretch_set_box:
 * @stretch: a #NcmFitESMCMCWalkerStretch
 * @n: parameter index
 * @lb: lower bound
 * @ub: upper bound
 * 
 * Sets box sampling for the @n-th parameter using @lb as lower bound
 * and @ub as upper bound.
 * 
 */
void 
ncm_fit_esmcmc_walker_stretch_set_box (NcmFitESMCMCWalkerStretch *stretch, guint n, const gdouble lb, const gdouble ub)
{
  g_assert_cmpuint (n, <, stretch->nparams);
  g_assert_cmpfloat (lb, <, ub);

  g_array_index (stretch->use_box, gboolean, n) = TRUE;
  ncm_matrix_set (stretch->box, n, 0, lb);
  ncm_matrix_set (stretch->box, n, 1, ub);
}

/**
 * ncm_fit_esmcmc_walker_stretch_set_box_mset:
 * @stretch: a #NcmFitESMCMCWalkerStretch
 * @mset: a #NcmMSet
 * 
 * Sets box sampling for the parameters using bounds from
 * @mset.
 * 
 */
void 
ncm_fit_esmcmc_walker_stretch_set_box_mset (NcmFitESMCMCWalkerStretch *stretch, NcmMSet *mset)
{
  const guint fparams_len = ncm_mset_fparams_len (mset);
  guint i;
  
  g_assert_cmpuint (stretch->nparams, ==, fparams_len);
  for (i = 0; i < fparams_len; i++)
  {
    ncm_fit_esmcmc_walker_stretch_set_box (stretch, i, 
                                           ncm_mset_fparam_get_lower_bound (mset, i),
                                           ncm_mset_fparam_get_upper_bound (mset, i));
  }  
}

/**
 * ncm_fit_esmcmc_walker_stretch_multi:
 * @stretch: a #NcmFitESMCMCWalkerStretch
 * @multi: a boolean
 * 
 * Sets whether it should use multi-stretchs in a single step.
 * 
 */
void 
ncm_fit_esmcmc_walker_stretch_multi (NcmFitESMCMCWalkerStretch *stretch, gboolean multi)
{
  stretch->multi = multi;
}
