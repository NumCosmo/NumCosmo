/***************************************************************************
 *            ncm_rng.c
 *
 *  Sat August 17 12:39:38 2013
 *  Copyright  2013  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * ncm_rng.c
 * Copyright (C) 2013 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * SECTION:ncm_rng
 * @title: NcmRNG
 * @short_description: Encapsulated GNU Scientific Library (GSL) random number generator with support for multhreading.
 * @stability: Stable
 * @include: numcosmo/math/ncm_rng.h
 *
 * This object encapsulates the [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) pseudo random number generator (PRNG).
 * Its main purpose is to add support for saving and loading state and multhreading.
 * For more information about the GSL routines see both links: [random number generation](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-generation)
 * and [random number distributions](https://www.gnu.org/software/gsl/doc/html/randist.html#random-number-distributions).
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_rng.h"
#include "math/ncm_cfg.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_randist.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _NcmRNGPrivate
{
  /*< private >*/
  GObject parent_instance;
  gsl_rng *r;
  gulong seed_val;
  gboolean seed_set;
  GMutex lock;
} NcmRNGPrivate;

struct _NcmRNGDiscrete
{
  /*< private >*/
  gsize n;
  gdouble *weights;
  gsl_ran_discrete_t *wran;
};

enum
{
  PROP_0,
  PROP_ALGO,
  PROP_STATE,
  PROP_SEED,
};

G_DEFINE_TYPE_WITH_PRIVATE (NcmRNG, ncm_rng, G_TYPE_OBJECT)
G_DEFINE_BOXED_TYPE (NcmRNGDiscrete, ncm_rng_discrete, ncm_rng_discrete_copy, ncm_rng_discrete_free)

static void
ncm_rng_init (NcmRNG *rng)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  self->r        = NULL;
  self->seed_val = 0;
  self->seed_set = FALSE;

  g_mutex_init (&self->lock);
}

static void
_ncm_rng_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_rng_parent_class)->constructed (object);
  {
    NcmRNG *rng                = NCM_RNG (object);
    NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

    if (!self->seed_set)
      ncm_rng_set_seed (rng, gsl_rng_default_seed);
  }
}

static void
_ncm_rng_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmRNG *rng = NCM_RNG (object);

  g_return_if_fail (NCM_IS_RNG (object));

  switch (prop_id)
  {
    case PROP_ALGO:
      ncm_rng_set_algo (rng, g_value_get_string (value));
      break;
    case PROP_STATE:
      ncm_rng_set_state (rng, g_value_get_string (value));
      break;
    case PROP_SEED:
      ncm_rng_set_seed (rng, g_value_get_ulong (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_rng_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmRNG *rng                = NCM_RNG (object);
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  g_return_if_fail (NCM_IS_RNG (object));

  switch (prop_id)
  {
    case PROP_ALGO:
      g_value_set_string (value, ncm_rng_get_algo (rng));
      break;
    case PROP_STATE:
      g_value_take_string (value, ncm_rng_get_state (rng));
      break;
    case PROP_SEED:
      g_value_set_ulong (value, self->seed_val);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_rng_finalize (GObject *object)
{
  NcmRNG *rng                = NCM_RNG (object);
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  g_clear_pointer (&self->r, gsl_rng_free);
  g_mutex_clear (&self->lock);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_rng_parent_class)->finalize (object);
}

static void
ncm_rng_class_init (NcmRNGClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->constructed  = &_ncm_rng_constructed;
  object_class->finalize     = &_ncm_rng_finalize;
  object_class->set_property = &_ncm_rng_set_property;
  object_class->get_property = &_ncm_rng_get_property;

  /**
   * NcmRNG:algorithm:
   *
   * The name of the pseudo random number algorithm to be used from [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/).
   * A list of the available algorithms can be find [here](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-generator-algorithms).
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_ALGO,
                                   g_param_spec_string ("algorithm",
                                                        NULL,
                                                        "Algorithm name",
                                                        gsl_rng_default->name,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmRNG:seed:
   *
   * Pseudo random number algorithm seed.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_SEED,
                                   g_param_spec_ulong ("seed",
                                                       NULL,
                                                       "Algorithm seed",
                                                       0, G_MAXULONG, 0,
                                                       G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmRNG:state:
   *
   * Pseudo random number algorithm state.
   *
   */
  g_object_class_install_property (object_class,
                                   PROP_STATE,
                                   g_param_spec_string ("state",
                                                        NULL,
                                                        "Algorithm state",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  /* Init the global gsl_rng variables */
  gsl_rng_env_setup ();

  klass->seed_gen  = g_rand_new ();
  klass->seed_hash = g_hash_table_new (g_direct_hash, g_direct_equal);
}

/**
 * ncm_rng_discrete_new:
 * @weights: (array length=n): array of weights
 * @n: number of elements in @weights
 *
 * Creates a new #NcmRNGDiscrete. This object is used to generate random numbers
 * from a discrete distribution for the given set of weights.
 *
 * Returns: (transfer full): a new #NcmRNGDiscrete.
 */
NcmRNGDiscrete *
ncm_rng_discrete_new (const gdouble *weights, const guint n)
{
  NcmRNGDiscrete *rng_discrete = g_slice_new (NcmRNGDiscrete);

  rng_discrete->n = n;
#if GLIB_CHECK_VERSION (2, 68, 0)
  rng_discrete->weights = g_memdup2 (weights, n * sizeof (gdouble));
#else
  rng_discrete->weights = g_memdup (weights, n * sizeof (gdouble));
#endif /* GLIB_CHECK_VERSION */
  rng_discrete->wran = gsl_ran_discrete_preproc (n, weights);

  return rng_discrete;
}

/**
 * ncm_rng_discrete_copy:
 * @rng: a #NcmRNGDiscrete
 *
 * Creates a copy of @rng_discrete.
 *
 * Returns: (transfer full): a copy of @rng_discrete.
 */
NcmRNGDiscrete *
ncm_rng_discrete_copy (NcmRNGDiscrete *rng)
{
  return ncm_rng_discrete_new (rng->weights, rng->n);
}

/**
 * ncm_rng_discrete_free:
 * @rng: a #NcmRNGDiscrete
 *
 * Frees the memory allocated by @rng.
 *
 */
void
ncm_rng_discrete_free (NcmRNGDiscrete *rng)
{
  g_clear_pointer (&rng->weights, g_free);
  g_clear_pointer (&rng->wran, gsl_ran_discrete_free);
  g_slice_free (NcmRNGDiscrete, rng);
}

/**
 * ncm_rng_new:
 * @algo: (allow-none): algorithm name
 *
 * Creates a new #NcmRNG using the algorithm @algo.
 * See the list of algorithms [here](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-generator-algorithms).
 * If @algo is NULL the default algorithm and seed are used.
 * See this [link](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-environment-variables) for more details.
 *
 * Returns: (transfer full): a new #NcmRNG.
 */
NcmRNG *
ncm_rng_new (const gchar *algo)
{
  NcmRNG *rng = g_object_new (NCM_TYPE_RNG,
                              "algorithm", algo,
                              NULL);

  return rng;
}

/**
 * ncm_rng_seeded_new:
 * @algo: (allow-none): algorithm name
 * @seed: seed used to initialize the PRNG
 *
 * Creates a new #NcmRNG using the algorithm @algo.
 * See the list of algorithms [here](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-generator-algorithms).
 * If @algo is NULL the default algorithm is used.
 * See this [link](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-environment-variables) for more details.
 *
 * Returns: (transfer full): a new #NcmRNG.
 */
NcmRNG *
ncm_rng_seeded_new (const gchar *algo, gulong seed)
{
  NcmRNG *rng = g_object_new (NCM_TYPE_RNG,
                              "algorithm", algo,
                              "seed", seed,
                              NULL);

  return rng;
}

/**
 * ncm_rng_ref:
 * @rng: a #NcmRNG
 *
 * Increases the reference count of @rng by one.
 *
 * Returns: (transfer full): @rng.
 */
NcmRNG *
ncm_rng_ref (NcmRNG *rng)
{
  return g_object_ref (rng);
}

/**
 * ncm_rng_free:
 * @rng: a #NcmRNG
 *
 * Decreases the reference count of @rng by one.
 *
 */
void
ncm_rng_free (NcmRNG *rng)
{
  g_object_unref (rng);
}

/**
 * ncm_rng_clear:
 * @rng: a #NcmRNG
 *
 * Decreases the reference count of *@rng by one and sets *@rng to NULL.
 *
 */
void
ncm_rng_clear (NcmRNG **rng)
{
  g_clear_object (rng);
}

/**
 * ncm_rng_lock:
 * @rng: a #NcmRNG
 *
 * Locks @rng.
 *
 */
void
ncm_rng_lock (NcmRNG *rng)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  g_mutex_lock (&self->lock);
}

/**
 * ncm_rng_unlock:
 * @rng: a #NcmRNG
 *
 * Unlocks @rng.
 *
 */
void
ncm_rng_unlock (NcmRNG *rng)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  g_mutex_unlock (&self->lock);
}

/**
 * ncm_rng_get_algo:
 * @rng: a #NcmRNG
 *
 * Gets the name of the algorithm.
 *
 * Returns: (transfer none): algorithm name.
 */
const gchar *
ncm_rng_get_algo (NcmRNG *rng)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_rng_name (self->r);
}

/**
 * ncm_rng_get_state:
 * @rng: a #NcmRNG
 *
 * Gets the state of the algorithm in [Base64](https://en.wikipedia.org/wiki/Base64).
 * It can be a very large string depending on the underlining state.
 *
 * Returns: (transfer full): algorithm state.
 */
gchar *
ncm_rng_get_state (NcmRNG *rng)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);
  gpointer state             = gsl_rng_state (self->r);
  gsize state_len            = gsl_rng_size (self->r);

  return g_base64_encode (state, state_len);
}

/**
 * ncm_rng_set_algo:
 * @rng: a #NcmRNG
 * @algo: algorithm name
 *
 * Sets the PRNG algorithm.
 *
 */
void
ncm_rng_set_algo (NcmRNG *rng, const gchar *algo)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);
  const gsl_rng_type *type;
  gboolean found = FALSE;

  if (algo != NULL)
  {
    const gsl_rng_type **t;
    const gsl_rng_type **t0;

    t0 = gsl_rng_types_setup ();

    for (t = t0; *t != 0; t++)
    {
      if (strcmp ((*t)->name, algo) == 0)
      {
        found = TRUE;
        break;
      }
    }

    if (!found)
      g_error ("ncm_rng_set_algo: cannot find algorithm %s.", algo);

    type = *t;
  }
  else
  {
    type = gsl_rng_default;
  }

  if (self->r == NULL)
  {
    self->r = gsl_rng_alloc (type);
  }
  else if (strcmp (gsl_rng_name (self->r), algo) != 0)
  {
    gsl_rng_free (self->r);
    self->r = gsl_rng_alloc (type);
  }
}

/**
 * ncm_rng_set_state:
 * @rng: a #NcmRNG
 * @state: algorithm state
 *
 * Sets the PRNG algorithm state.
 *
 */
void
ncm_rng_set_state (NcmRNG *rng, const gchar *state)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);
  gpointer state_ptr         = gsl_rng_state (self->r);
  gsize state_len            = gsl_rng_size (self->r);
  gsize state_dec_len        = 0;
  guchar *decoded_state      = g_base64_decode (state, &state_dec_len);

  g_assert_cmpuint (state_len, ==, state_dec_len);

  memcpy (state_ptr, decoded_state, state_len);

  g_free (decoded_state);
}

/**
 * ncm_rng_check_seed:
 * @rng: a #NcmRNG
 * @seed: seed for the PRNG
 *
 * Check if the seed was already used by any #NcmRNG.
 *
 * Returns: TRUE if @seed was never used and FALSE otherwise.
 */
gboolean
ncm_rng_check_seed (NcmRNG *rng, gulong seed)
{
  NcmRNGClass *rng_class     = NCM_RNG_GET_CLASS (rng);
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);
  gint seed_int              = seed;
  gpointer b                 = g_hash_table_lookup (rng_class->seed_hash, GINT_TO_POINTER (seed_int));

  return GPOINTER_TO_INT (b) == 0;
}

/**
 * ncm_rng_set_seed:
 * @rng: a #NcmRNG
 * @seed: seed for the PRNG
 *
 * Sets the PRNG algorithm seed.
 *
 */
void
ncm_rng_set_seed (NcmRNG *rng, gulong seed)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  self->seed_val = seed;

  if (self->r != NULL)
  {
    NcmRNGClass *rng_class = NCM_RNG_GET_CLASS (rng);
    gint seed_int          = seed;

    gsl_rng_set (self->r, seed);
    g_hash_table_insert (rng_class->seed_hash, GINT_TO_POINTER (seed_int), GINT_TO_POINTER (1));
    self->seed_set = TRUE;
  }
}

/**
 * ncm_rng_get_seed:
 * @rng: a #NcmRNG
 *
 * This functions returns the seed used to initialize the PRNG.
 *
 * Returns: @rng's @seed.
 *
 */
gulong
ncm_rng_get_seed (NcmRNG *rng)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return self->seed_val;
}

/**
 * ncm_rng_set_random_seed:
 * @rng: a #NcmRNG
 * @allow_colisions: a gboolean
 *
 * Sets the algorithm seed using a PRNG seeded by /dev/urandom (Unix/Linux)
 * or current time, when the first is not available (see #g_rand_new).
 * If @allow_colisions is FALSE this function will set the first unused seed generated.
 *
 */
void
ncm_rng_set_random_seed (NcmRNG *rng, gboolean allow_colisions)
{
  NcmRNGClass *rng_class = NCM_RNG_GET_CLASS (rng);
  gulong seed            = g_rand_int (rng_class->seed_gen) + 1;

  while (!ncm_rng_check_seed (rng, seed))
    seed = g_rand_int (rng_class->seed_gen) + 1;

  ncm_rng_set_seed (rng, seed);
}

static GHashTable *rng_table = NULL;

/**
 * ncm_rng_pool_get:
 * @name: a string
 *
 * Gets the #NcmRNG named @name from the pool.
 * If it doesn't exists, it creates one, add it to the pool and returns it.
 *
 * Returns: (transfer full): the #NcmRNG named @name.
 */
NcmRNG *
ncm_rng_pool_get (const gchar *name)
{
  NcmRNG *rng;

  G_LOCK_DEFINE_STATIC (create_lock);
  G_LOCK_DEFINE_STATIC (update_acess_lock);

  if (rng_table == NULL)
  {
    G_LOCK (create_lock);

    if (rng_table == NULL)
      rng_table = g_hash_table_new_full (g_str_hash, g_str_equal,
                                         &g_free, (GDestroyNotify) & ncm_rng_free);

    G_UNLOCK (create_lock);
  }

  G_LOCK (update_acess_lock);
  {
    rng = g_hash_table_lookup (rng_table, name);

    if (rng == NULL)
    {
      rng = ncm_rng_new (NULL);
      g_hash_table_insert (rng_table,
                           g_strdup (name),
                           ncm_rng_ref (rng));
    }
    else
    {
      ncm_rng_ref (rng);
    }
  }
  G_UNLOCK (update_acess_lock);

  return rng;
}

/**
 * ncm_rng_gen_ulong:
 * @rng: a #NcmRNG
 *
 * This function returns a random unsigned integer from the uniform distribution.
 *
 * Returns: a random unsigned long from the uniform distribution.
 */
gulong
ncm_rng_gen_ulong (NcmRNG *rng)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_rng_get (self->r);
}

/**
 * ncm_rng_uniform_int_gen:
 * @rng: a #NcmRNG
 * @n: upper limit
 *
 * This function returns a random number drawn from the
 * uniform distribution between zero and @n.
 *
 * Returns: a random number from the uniform distribution.
 */
gulong
ncm_rng_uniform_int_gen (NcmRNG *rng, gulong n)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_rng_uniform_int (self->r, n);
}

/**
 * ncm_rng_uniform01_gen:
 * @rng: a #NcmRNG
 *
 * This function returns a random number drawn from the
 * uniform distribution between zero and one $[0,1)$.
 *
 * Returns: a random number from the uniform distribution.
 */
gdouble
ncm_rng_uniform01_gen (NcmRNG *rng)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_rng_uniform (self->r);
}

/**
 * ncm_rng_uniform01_pos_gen:
 * @rng: a #NcmRNG
 *
 * This function returns a random number drawn from the
 * uniform distribution between zero and one $(0,1)$.
 *
 * Returns: a random number from the uniform distribution.
 */
gdouble
ncm_rng_uniform01_pos_gen (NcmRNG *rng)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_rng_uniform_pos (self->r);
}

/**
 * ncm_rng_uniform_gen:
 * @rng: a #NcmRNG
 * @xl: lower value
 * @xu: upper value
 *
 * This functions returns a random number drawn from the
 * uniform distribution between the values @xl and @xu.
 *
 * Returns: a random number from the uniform distribution.
 */
gdouble
ncm_rng_uniform_gen (NcmRNG *rng, const gdouble xl, const gdouble xu)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_flat (self->r, xl, xu);
}

/**
 * ncm_rng_gaussian_gen:
 * @rng: a #NcmRNG
 * @mu: mean
 * @sigma: standard deviation
 *
 * This function returns a random number drawn from the
 * [Gaussian distribution](https://en.wikipedia.org/wiki/Normal_distribution),
 * with mean @mu and standard deviation @sigma.
 *
 * Returns: a random number from the Gaussian distribution.
 */
gdouble
ncm_rng_gaussian_gen (NcmRNG *rng, const gdouble mu, const gdouble sigma)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_gaussian (self->r, sigma) + mu;
}

/**
 * ncm_rng_ugaussian_gen:
 * @rng: a #NcmRNG
 *
 * This function returns a random number drwan from the
 * [Gaussian distribution](https://en.wikipedia.org/wiki/Normal_distribution),
 * with mean zero and standard deviation one.
 * Equivalent as above but with @mean = 0 and @sigma = 1.
 *
 * Returns: a random number from the Gaussian distribution.
 */
gdouble
ncm_rng_ugaussian_gen (NcmRNG *rng)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_ugaussian (self->r);
}

/**
 * ncm_rng_gaussian_tail_gen:
 * @rng: a #NcmRNG
 * @a: positive lower limit
 * @sigma: standard deviation
 *
 * This function returns a random number drawn from the upper tail of the
 * [Gaussian distribution](https://en.wikipedia.org/wiki/Normal_distribution) with standard deviation @sigma.
 * The value returned is larger than the lower limit @a, which must be positive.
 *
 * Returns: a random number from the Gaussian distribution tail.
 */
gdouble
ncm_rng_gaussian_tail_gen (NcmRNG *rng, const gdouble a, const gdouble sigma)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_gaussian_tail (self->r, a, sigma);
}

/**
 * ncm_rng_exponential_gen:
 * @rng: a #NcmRNG
 * @mu: scale parameter
 *
 * This function returns a random number drawn from the
 * [exponential distribution](https://en.wikipedia.org/wiki/Exponential_distribution)
 * with scale parameter (mean) @mu.
 *
 * Returns: a random number from the exponential distribution.
 */
gdouble
ncm_rng_exponential_gen (NcmRNG *rng, const gdouble mu)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_exponential (self->r, mu);
}

/**
 * ncm_rng_laplace_gen:
 * @rng: a #NcmRNG
 * @a: width of the distribution
 *
 * This function returns a random number drawn from the
 * [Laplace distribution](https://en.wikipedia.org/wiki/Laplace_distribution)
 * with width @a.
 *
 * Returns: a random number from the Laplace distribution.
 */
gdouble
ncm_rng_laplace_gen (NcmRNG *rng, const gdouble a)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_laplace (self->r, a);
}

/**
 * ncm_rng_exppow_gen:
 * @rng: a #NcmRNG
 * @a: scale parameter
 * @b: exponent
 *
 * This function returns a random number drawn from the
 * [exponential power distribution](https://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1)
 * with scale parameter @a and exponent @b.
 *
 * Returns: a random number from the exponential power distribution.
 */
gdouble
ncm_rng_exppow_gen (NcmRNG *rng, const gdouble a, const gdouble b)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_exppow (self->r, a, b);
}

/**
 * ncm_rng_beta_gen:
 * @rng: a #NcmRNG
 * @a: shape parameter
 * @b: shape parameter
 *
 * This function returns a random number drawn from the
 * [beta distribution](https://en.wikipedia.org/wiki/Beta_distribution)
 * with shape parameters @a and @b. The shape parameters must be positive.
 *
 * Returns: a random number from the beta distribution.
 */
gdouble
ncm_rng_beta_gen (NcmRNG *rng, const gdouble a, const gdouble b)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_beta (self->r, a, b);
}

/**
 * ncm_rng_gamma_gen:
 * @rng: a #NcmRNG
 * @a: shape parameter
 * @b: scale parameter
 *
 * This function returns a random number drawn from the
 * [gamma distribution](https://en.wikipedia.org/wiki/Gamma_distribution)
 * with shape parameter @a and scale parameter @b.
 *
 * Returns: a random number from the gamma distribution.
 */
gdouble
ncm_rng_gamma_gen (NcmRNG *rng, const gdouble a, const gdouble b)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_gamma (self->r, a, b);
}

/**
 * ncm_rng_chisq_gen:
 * @rng: a #NcmRNG
 * @nu: degrees of freedom $\nu$
 *
 * This function returns a random number drawn from the
 * [Chi-square Distribution](https://en.wikipedia.org/wiki/Chi-square_distribution),
 * with $\nu$ degrees of freedom.
 * Returns: a random number from Chi-square distribution.
 */
gdouble
ncm_rng_chisq_gen (NcmRNG *rng, const gdouble nu)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_chisq (self->r, nu);
}

/**
 * ncm_rng_poisson_gen:
 * @rng: a #NcmRNG
 * @mu: degrees of freedom $\nu$
 *
 * This function returns a random number drawn from the Poisson distribution,
 * with frequency @mu.
 *
 * Returns: a random number from the Poisson distribution.
 */
gdouble
ncm_rng_poisson_gen (NcmRNG *rng, const gdouble mu)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_poisson (self->r, mu);
}

/**
 * ncm_rng_rayleigh_gen:
 * @rng: a #NcmRNG
 * @sigma: scale parameter
 *
 * This function returns a random number drawn from the
 * [Rayleigh distribution](https://en.wikipedia.org/wiki/Rayleigh_distribution)
 * with scale parameter @sigma.
 *
 * Returns: a random number from the Rayleigh distribution.
 */
gdouble
ncm_rng_rayleigh_gen (NcmRNG *rng, const gdouble sigma)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_rayleigh (self->r, sigma);
}

/**
 * ncm_rng_discrete_gen:
 * @rng: a #NcmRNG
 * @rng_discrete: a #NcmRNGDiscrete
 *
 * This function returns a random number drawn from the discrete distribution. The
 * weights must created using ncm_rng_discrete_new().
 *
 */
gsize
ncm_rng_discrete_gen (NcmRNG *rng, NcmRNGDiscrete *rng_discrete)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  return gsl_ran_discrete (self->r, rng_discrete->wran);
}

/**
 * ncm_rng_sample:
 * @rng: a #NcmRNG
 * @dest: an array of @k elements of size @size
 * @k: number of elements in @dest
 * @src: an array of @n elements of size @size
 * @n: number of elements in @src
 * @size: size of each element in @dest and @src
 *
 * This function fills the array @dest with @k elements from the array @src.
 * The elements are chosen randomly using the algorithm (sample with replecement) in
 * [gsl_ran_sample()](https://www.gnu.org/software/gsl/doc/html/randist.html#c.gsl_ran_sample).
 *
 */
void
ncm_rng_sample (NcmRNG *rng, void *dest, size_t k, void *src, size_t n, size_t size)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  gsl_ran_sample (self->r, dest, k, src, n, size);
}

/**
 * ncm_rng_choose:
 * @rng: a #NcmRNG
 * @dest: an array of @k elements of size @size
 * @k: number of elements in @dest
 * @src: an array of @n elements of size @size
 * @n: number of elements in @src
 * @size: size of each element in @dest and @src
 *
 * This function fills the array @dest with @k elements from the array @src.
 * The elements are chosen randomly using the algorithm (choose with replecement) in
 * [gsl_ran_choose()](https://www.gnu.org/software/gsl/doc/html/randist.html#c.gsl_ran_choose).
 *
 */
void
ncm_rng_choose (NcmRNG *rng, void *dest, size_t k, void *src, size_t n, size_t size)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  gsl_ran_choose (self->r, dest, k, src, n, size);
}

/**
 * ncm_rng_multinomial:
 * @rng: a #NcmRNG
 * @K: number of possible outcomes
 * @N: number of trials
 * @p: (array length=K) (element-type gdouble): array of probabilities
 * @n: (array length=K) (element-type guint): array of counts
 *
 * This function fills the array @n with @K elements using a multinomial distribution
 * defined by the array @p.
 *
 */
void
ncm_rng_multinomial (NcmRNG *rng, gsize K, guint N, const gdouble *p, guint *n)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  gsl_ran_multinomial (self->r, K, N, p, n);
}

/**
 * ncm_rng_bivariate_gaussian_gen:
 * @rng: a #NcmRNG
 * @sigma_x: standard deviation
 * @sigma_y: standard deviation
 * @rho: correlation coefficient
 * @x: (out): random number from the Bivariate Gaussian distribution
 * @y: (out): random number from the Bivariate Gaussian distribution
 *
 * This function returns a random number drawn from the
 * [Bivariate Gaussian distribution](https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Bivariate_case),
 * with standard deviations @sigma_x and @sigma_y and correlation coefficient @rho.
 * The correlation coefficient must be in the range $-1 \leq \rho \leq 1$.
 *
 */
void
ncm_rng_bivariate_gaussian_gen (NcmRNG *rng, const gdouble sigma_x, const gdouble sigma_y, const gdouble rho, gdouble *x, gdouble *y)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  gsl_ran_bivariate_gaussian (self->r, sigma_x, sigma_y, rho, x, y);
}

