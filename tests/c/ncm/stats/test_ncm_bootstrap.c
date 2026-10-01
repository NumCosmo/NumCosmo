/***************************************************************************
 *            test_ncm_bootstrap.c
 *
 *  Tue September 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_bootstrap.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_BOOTSTRAP_SEED 123456
#define TEST_BOOTSTRAP_NDRAWS 10000

static void
test_ncm_bootstrap_sizes (void)
{
  NcmBootstrap *bstrap = ncm_bootstrap_new ();
  NcmBootstrap *sized  = ncm_bootstrap_sized_new (7);
  NcmBootstrap *full   = ncm_bootstrap_full_new (7, 3);

  g_assert_cmpuint (ncm_bootstrap_get_fsize (bstrap), ==, 0);
  g_assert_cmpuint (ncm_bootstrap_get_bsize (bstrap), ==, 0);
  g_assert_cmpuint (ncm_bootstrap_get_fsize (sized), ==, 7);
  g_assert_cmpuint (ncm_bootstrap_get_bsize (sized), ==, 7);
  g_assert_cmpuint (ncm_bootstrap_get_fsize (full), ==, 7);
  g_assert_cmpuint (ncm_bootstrap_get_bsize (full), ==, 3);
  g_assert_false (ncm_bootstrap_is_init (full));

  /* The full size does not change the bootstrap size */
  ncm_bootstrap_set_fsize (full, 11);
  g_assert_cmpuint (ncm_bootstrap_get_fsize (full), ==, 11);
  g_assert_cmpuint (ncm_bootstrap_get_bsize (full), ==, 3);

  ncm_bootstrap_set_bsize (full, 5);
  g_assert_cmpuint (ncm_bootstrap_get_fsize (full), ==, 11);
  g_assert_cmpuint (ncm_bootstrap_get_bsize (full), ==, 5);

  ncm_bootstrap_free (bstrap);
  ncm_bootstrap_free (sized);
  ncm_bootstrap_free (full);
}

/* A new size discards the realization, the same size keeps it */
static void
test_ncm_bootstrap_resize (void)
{
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, TEST_BOOTSTRAP_SEED);
  NcmBootstrap *bstrap = ncm_bootstrap_full_new (7, 3);

  ncm_bootstrap_resample (bstrap, rng);
  ncm_bootstrap_set_fsize (bstrap, 7);
  ncm_bootstrap_set_bsize (bstrap, 3);
  g_assert_true (ncm_bootstrap_is_init (bstrap));

  ncm_bootstrap_set_fsize (bstrap, 8);
  g_assert_false (ncm_bootstrap_is_init (bstrap));

  ncm_bootstrap_resample (bstrap, rng);
  ncm_bootstrap_set_bsize (bstrap, 4);
  g_assert_false (ncm_bootstrap_is_init (bstrap));

  ncm_bootstrap_free (bstrap);
  ncm_rng_free (rng);
}

/* Each index is drawn bsize / fsize times per realization on average */
static void
test_ncm_bootstrap_resample (void)
{
  const guint fsize    = 10;
  const guint bsize    = 20;
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, TEST_BOOTSTRAP_SEED);
  NcmBootstrap *bstrap = ncm_bootstrap_full_new (fsize, bsize);
  guint count[10]      = {0};
  guint i, n;

  ncm_bootstrap_resample (bstrap, rng);
  g_assert_true (ncm_bootstrap_is_init (bstrap));

  for (n = 0; n < TEST_BOOTSTRAP_NDRAWS; n++)
  {
    ncm_bootstrap_resample (bstrap, rng);

    for (i = 0; i < bsize; i++)
    {
      const guint k = ncm_bootstrap_get (bstrap, i);

      g_assert_cmpuint (k, <, fsize);
      count[k]++;
    }
  }

  {
    const gdouble p     = 1.0 / fsize;
    const gdouble mean  = TEST_BOOTSTRAP_NDRAWS * bsize * p;
    const gdouble sigma = sqrt (mean * (1.0 - p));

    for (i = 0; i < fsize; i++)
      ncm_assert_cmpdouble_e (count[i], ==, mean, 0.0, 5.0 * sigma);
  }

  ncm_bootstrap_free (bstrap);
  ncm_rng_free (rng);
}

/* Distinct indexes in increasing order, each chosen with probability bsize / fsize */
static void
test_ncm_bootstrap_remix (void)
{
  const guint fsize    = 10;
  const guint bsize    = 4;
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, TEST_BOOTSTRAP_SEED);
  NcmBootstrap *bstrap = ncm_bootstrap_full_new (fsize, bsize);
  guint count[10]      = {0};
  guint i, n;

  for (n = 0; n < TEST_BOOTSTRAP_NDRAWS; n++)
  {
    ncm_bootstrap_remix (bstrap, rng);

    for (i = 0; i < bsize; i++)
    {
      const guint k = ncm_bootstrap_get (bstrap, i);

      g_assert_cmpuint (k, <, fsize);

      if (i > 0)
        g_assert_cmpuint (ncm_bootstrap_get (bstrap, i - 1), <, k);

      count[k]++;
    }
  }

  {
    const gdouble p     = (gdouble) bsize / fsize;
    const gdouble mean  = TEST_BOOTSTRAP_NDRAWS * p;
    const gdouble sigma = sqrt (mean * (1.0 - p));

    for (i = 0; i < fsize; i++)
      ncm_assert_cmpdouble_e (count[i], ==, mean, 0.0, 5.0 * sigma);
  }

  /* bsize = fsize selects every index */
  ncm_bootstrap_set_bsize (bstrap, fsize);
  ncm_bootstrap_remix (bstrap, rng);

  for (i = 0; i < fsize; i++)
    g_assert_cmpuint (ncm_bootstrap_get (bstrap, i), ==, i);

  ncm_bootstrap_free (bstrap);
  ncm_rng_free (rng);
}

/* Pairs (index, count) in increasing index order, counts summing to bsize, realization unchanged */
static void
test_ncm_bootstrap_sortncomp (void)
{
  const guint fsize    = 8;
  const guint bsize    = 30;
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, TEST_BOOTSTRAP_SEED);
  NcmBootstrap *bstrap = ncm_bootstrap_full_new (fsize, bsize);
  guint before[30];
  guint count[8] = {0};
  guint total    = 0;
  GArray *pairs;
  guint i;

  ncm_bootstrap_resample (bstrap, rng);

  for (i = 0; i < bsize; i++)
  {
    before[i] = ncm_bootstrap_get (bstrap, i);
    count[before[i]]++;
  }

  pairs = ncm_bootstrap_get_sortncomp (bstrap);
  g_assert_cmpuint (pairs->len % 2, ==, 0);

  for (i = 0; i < pairs->len / 2; i++)
  {
    const guint k   = g_array_index (pairs, guint, 2 * i + 0);
    const guint c_k = g_array_index (pairs, guint, 2 * i + 1);

    if (i > 0)
      g_assert_cmpuint (g_array_index (pairs, guint, 2 * i - 2), <, k);

    g_assert_cmpuint (c_k, ==, count[k]);
    total += c_k;
  }

  g_assert_cmpuint (total, ==, bsize);

  for (i = 0; i < bsize; i++)
    g_assert_cmpuint (ncm_bootstrap_get (bstrap, i), ==, before[i]);

  g_array_unref (pairs);

  /* An empty realization gives no pairs */
  ncm_bootstrap_set_bsize (bstrap, 0);
  ncm_bootstrap_resample (bstrap, rng);
  pairs = ncm_bootstrap_get_sortncomp (bstrap);
  g_assert_cmpuint (pairs->len, ==, 0);
  g_array_unref (pairs);

  ncm_bootstrap_free (bstrap);
  ncm_rng_free (rng);
}

/* The realization survives a serialization round trip */
static void
test_ncm_bootstrap_serialize (void)
{
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, TEST_BOOTSTRAP_SEED);
  NcmBootstrap *bstrap = ncm_bootstrap_full_new (9, 6);
  NcmSerialize *ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmBootstrap *empty  = ncm_bootstrap_full_new (9, 6);
  NcmBootstrap *bstrap_dup;
  guint i;

  ncm_bootstrap_resample (bstrap, rng);
  bstrap_dup = NCM_BOOTSTRAP (ncm_serialize_dup_obj (ser, G_OBJECT (bstrap)));

  g_assert_true (ncm_bootstrap_is_init (bstrap_dup));
  g_assert_cmpuint (ncm_bootstrap_get_fsize (bstrap_dup), ==, 9);
  g_assert_cmpuint (ncm_bootstrap_get_bsize (bstrap_dup), ==, 6);

  for (i = 0; i < 6; i++)
    g_assert_cmpuint (ncm_bootstrap_get (bstrap_dup, i), ==, ncm_bootstrap_get (bstrap, i));

  ncm_bootstrap_free (bstrap_dup);

  /* Without a realization the copy has none */
  bstrap_dup = NCM_BOOTSTRAP (ncm_serialize_dup_obj (ser, G_OBJECT (empty)));
  g_assert_false (ncm_bootstrap_is_init (bstrap_dup));

  ncm_bootstrap_free (bstrap_dup);
  ncm_bootstrap_free (empty);
  ncm_bootstrap_free (bstrap);
  ncm_serialize_free (ser);
  ncm_rng_free (rng);
}

static void
test_ncm_bootstrap_realization_length_subprocess (void)
{
  NcmBootstrap *bstrap = ncm_bootstrap_sized_new (3);

  g_object_set (bstrap, "realization", g_variant_new_parsed ("[@u 0, 1]"), NULL);
}

static void
test_ncm_bootstrap_realization_index_subprocess (void)
{
  NcmBootstrap *bstrap = ncm_bootstrap_sized_new (3);

  g_object_set (bstrap, "realization", g_variant_new_parsed ("[@u 0, 1, 3]"), NULL);
}

static void
test_ncm_bootstrap_resample_empty_subprocess (void)
{
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, TEST_BOOTSTRAP_SEED);
  NcmBootstrap *bstrap = ncm_bootstrap_full_new (0, 3);

  ncm_bootstrap_resample (bstrap, rng);
}

static void
test_ncm_bootstrap_remix_too_large_subprocess (void)
{
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, TEST_BOOTSTRAP_SEED);
  NcmBootstrap *bstrap = ncm_bootstrap_full_new (3, 5);

  ncm_bootstrap_remix (bstrap, rng);
}

static void
test_ncm_bootstrap_sortncomp_no_realization_subprocess (void)
{
  NcmBootstrap *bstrap = ncm_bootstrap_sized_new (3);

  ncm_bootstrap_get_sortncomp (bstrap);
}

static void
test_ncm_bootstrap_traps (void)
{
  g_test_trap_subprocess ("/ncm/bootstrap/realization_length/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/bootstrap/realization_index/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/bootstrap/resample_empty/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/bootstrap/remix_too_large/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/bootstrap/sortncomp_no_realization/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/bootstrap/sizes", &test_ncm_bootstrap_sizes);
  g_test_add_func ("/ncm/bootstrap/resize", &test_ncm_bootstrap_resize);
  g_test_add_func ("/ncm/bootstrap/resample", &test_ncm_bootstrap_resample);
  g_test_add_func ("/ncm/bootstrap/remix", &test_ncm_bootstrap_remix);
  g_test_add_func ("/ncm/bootstrap/sortncomp", &test_ncm_bootstrap_sortncomp);
  g_test_add_func ("/ncm/bootstrap/serialize", &test_ncm_bootstrap_serialize);
  g_test_add_func ("/ncm/bootstrap/traps", &test_ncm_bootstrap_traps);
  g_test_add_func ("/ncm/bootstrap/realization_length/subprocess", &test_ncm_bootstrap_realization_length_subprocess);
  g_test_add_func ("/ncm/bootstrap/realization_index/subprocess", &test_ncm_bootstrap_realization_index_subprocess);
  g_test_add_func ("/ncm/bootstrap/resample_empty/subprocess", &test_ncm_bootstrap_resample_empty_subprocess);
  g_test_add_func ("/ncm/bootstrap/remix_too_large/subprocess", &test_ncm_bootstrap_remix_too_large_subprocess);
  g_test_add_func ("/ncm/bootstrap/sortncomp_no_realization/subprocess", &test_ncm_bootstrap_sortncomp_no_realization_subprocess);

  g_test_run ();
}

