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
 * NcmRNG:
 *
 * Pseudo-random number generator wrapping a GSL `gsl_rng`.
 *
 * Adds a serializable state (`NcmRNG:state`), a process-wide record of the seeds
 * already used, a pool of named generators and a mutex. The generation functions do
 * not lock; code sharing an #NcmRNG between threads brackets its calls with
 * ncm_rng_lock() and ncm_rng_unlock(). See the GSL documentation on
 * [random number generation](https://www.gnu.org/software/gsl/doc/html/rng.html) and
 * [random number distributions](https://www.gnu.org/software/gsl/doc/html/randist.html).
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_rng.h"
#include "ncm/core/ncm_cfg.h"

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
   * The GSL name of the algorithm, one of the
   * [GSL generators](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-generator-algorithms).
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
   * The last seed set.
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
   * The state of the generator encoded in Base64, see ncm_rng_get_state().
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
 * @weights: (array length=n): the weights
 * @n: number of weights
 *
 * Creates the lookup table for sampling the indexes $[0, n)$ with probabilities
 * proportional to @weights, for ncm_rng_discrete_gen().
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
 * Returns: (transfer full): a copy of @rng.
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
 * Frees @rng.
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
 * Creates a new #NcmRNG with the GSL algorithm @algo, or the GSL default if @algo is
 * %NULL. The seed is the GSL default seed. Both defaults are set by the
 * [GSL environment variables](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-environment-variables).
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
 * @seed: the seed
 *
 * Same as ncm_rng_new(), with seed @seed.
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
 * Decreases the reference count of *@rng by one and sets *@rng to %NULL.
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
 * Locks the mutex of @rng.
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
 * Unlocks the mutex of @rng.
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
 * Returns: (transfer none): the GSL name of the algorithm.
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
 * Encodes the state of the generator in Base64. Its length is proportional to the
 * size of the GSL state, which depends on the algorithm.
 *
 * Returns: (transfer full): the encoded state.
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
 * Sets the GSL algorithm, one of the
 * [GSL generators](https://www.gnu.org/software/gsl/doc/html/rng.html#random-number-generator-algorithms). Replacing the
 * algorithm allocates a new generator with the GSL default seed. Aborts if @algo is
 * not a GSL algorithm.
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
 * @state: a state from ncm_rng_get_state()
 *
 * Restores the state of the generator. @state must come from a generator with the
 * same algorithm; aborts if its length differs.
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
 * @seed: a seed
 *
 * Returns: whether no #NcmRNG in the process has been seeded with @seed.
 */
gboolean
ncm_rng_check_seed (NcmRNG *rng, gulong seed)
{
  NcmRNGClass *rng_class = NCM_RNG_GET_CLASS (rng);
  gint seed_int          = seed;
  gpointer b             = g_hash_table_lookup (rng_class->seed_hash, GINT_TO_POINTER (seed_int));

  return GPOINTER_TO_INT (b) == 0;
}

/**
 * ncm_rng_set_seed:
 * @rng: a #NcmRNG
 * @seed: the seed
 *
 * Seeds the generator with @seed and records it as used.
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
 * Returns: the last seed set.
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
 * @allow_colisions: whether a used seed is accepted
 *
 * Seeds the generator with a positive 32-bit seed drawn from a #GRand, itself seeded
 * from /dev/urandom or, if unavailable, the current time (see g_rand_new()). Unless
 * @allow_colisions is %TRUE, draws again until the seed passes
 * ncm_rng_check_seed().
 */
void
ncm_rng_set_random_seed (NcmRNG *rng, gboolean allow_colisions)
{
  NcmRNGClass *rng_class = NCM_RNG_GET_CLASS (rng);
  gulong seed            = g_rand_int (rng_class->seed_gen) + 1;

  if (!allow_colisions)
    while (!ncm_rng_check_seed (rng, seed))
      seed = g_rand_int (rng_class->seed_gen) + 1;


  ncm_rng_set_seed (rng, seed);
}

static GHashTable *rng_table = NULL;

/**
 * ncm_rng_pool_get:
 * @name: the name
 *
 * Returns the process-wide #NcmRNG named @name, creating it with ncm_rng_new() on
 * first use. Thread-safe.
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
 * Returns: a uniform integer between the minimum and maximum of the algorithm,
 * see gsl_rng_get().
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
 * @n: number of values
 *
 * Returns: a uniform integer in $[0, n - 1]$.
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
 * Returns: a uniform number in $[0, 1)$.
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
 * Returns: a uniform number in $(0, 1)$.
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
 * @xl: lower limit
 * @xu: upper limit
 *
 * Returns: a uniform number in $[x_l, x_u)$.
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
 * Returns: a Gaussian number with mean @mu and standard deviation @sigma.
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
 * Returns: a Gaussian number with zero mean and unit standard deviation.
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
 * @a: lower limit
 * @sigma: standard deviation
 *
 * Draws from the zero-mean Gaussian of standard deviation @sigma restricted to
 * $x > a$, with $a > 0$.
 *
 * Returns: the number drawn.
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
 * @mu: mean
 *
 * Draws from $p(x) = e^{-x/\mu}/\mu$, $x \geq 0$.
 *
 * Returns: the number drawn.
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
 * @a: width
 *
 * Draws from $p(x) = e^{-|x|/a}/(2a)$.
 *
 * Returns: the number drawn.
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
 * @a: scale
 * @b: exponent
 *
 * Draws from $p(x) = e^{-|x/a|^b} / [2a\,\Gamma(1 + 1/b)]$.
 *
 * Returns: the number drawn.
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
 * @a: first shape parameter
 * @b: second shape parameter
 *
 * Draws from $p(x) \propto x^{a-1} (1 - x)^{b-1}$, $0 \leq x \leq 1$, with $a, b > 0$.
 *
 * Returns: the number drawn.
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
 * @a: shape
 * @b: scale
 *
 * Draws from $p(x) = x^{a-1} e^{-x/b} / [\Gamma(a)\, b^a]$, $x > 0$.
 *
 * Returns: the number drawn.
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
 * Returns: a $\chi^2$ number with $\nu$ degrees of freedom.
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
 * @mu: mean
 *
 * Returns: a Poisson count with mean @mu.
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
 * @sigma: scale
 *
 * Draws from $p(x) = (x/\sigma^2)\, e^{-x^2/(2\sigma^2)}$, $x > 0$.
 *
 * Returns: the number drawn.
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
 * Returns: an index drawn with the probabilities of @rng_discrete.
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
 * @dest: array of @k elements
 * @k: number of elements of @dest
 * @src: array of @n elements
 * @n: number of elements of @src
 * @size: size in bytes of each element
 *
 * Fills @dest with @k elements of @src drawn with replacement, see
 * [gsl_ran_sample()](https://www.gnu.org/software/gsl/doc/html/randist.html#c.gsl_ran_sample).
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
 * @dest: array of @k elements
 * @k: number of elements of @dest
 * @src: array of @n elements
 * @n: number of elements of @src
 * @size: size in bytes of each element
 *
 * Fills @dest with @k distinct elements of @src drawn without replacement, in their
 * order in @src, with $k \leq n$, see
 * [gsl_ran_choose()](https://www.gnu.org/software/gsl/doc/html/randist.html#c.gsl_ran_choose).
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
 * @K: number of outcomes
 * @N: number of trials
 * @p: (array length=K) (element-type gdouble): probabilities
 * @n: (array length=K) (element-type guint): counts
 *
 * Fills @n with the counts of @N trials of a multinomial distribution with
 * probabilities proportional to @p.
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
 * @sigma_x: standard deviation of $x$
 * @sigma_y: standard deviation of $y$
 * @rho: correlation coefficient
 * @x: (out): the first component
 * @y: (out): the second component
 *
 * Draws a pair from the zero-mean bivariate Gaussian with standard deviations @sigma_x
 * and @sigma_y and correlation coefficient $-1 \leq \rho \leq 1$.
 */
void
ncm_rng_bivariate_gaussian_gen (NcmRNG *rng, const gdouble sigma_x, const gdouble sigma_y, const gdouble rho, gdouble *x, gdouble *y)
{
  NcmRNGPrivate * const self = ncm_rng_get_instance_private (rng);

  gsl_ran_bivariate_gaussian (self->r, sigma_x, sigma_y, rho, x, y);
}

