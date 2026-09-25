/***************************************************************************
 *            test_ncm_function_cache.c
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_function_cache.c
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

static NcmFunctionCache *
_test_function_cache_new (void)
{
  NcmFunctionCache *cache = ncm_function_cache_new (2, 0.0, 1.0e-7);
  const gdouble x[]       = {1.0, 2.0, 4.0};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (x); i++)
  {
    NcmVector *p = ncm_vector_new (2);

    ncm_vector_set (p, 0, x[i]);
    ncm_vector_set (p, 1, x[i] * x[i]);
    ncm_function_cache_insert_vector (cache, x[i], p);
    ncm_vector_free (p);
  }

  return cache;
}

static void
test_ncm_function_cache_tolerances (void)
{
  NcmFunctionCache *cache = ncm_function_cache_new (3, 1.0e-3, 1.0e-5);
  guint dim;

  g_assert_cmpfloat (ncm_function_cache_get_abstol (cache), ==, 1.0e-3);
  g_assert_cmpfloat (ncm_function_cache_get_reltol (cache), ==, 1.0e-5);

  g_object_get (cache, "dimension", &dim, NULL);
  g_assert_cmpuint (dim, ==, 3);

  ncm_function_cache_free (cache);
}

static void
test_ncm_function_cache_get (void)
{
  NcmFunctionCache *cache = _test_function_cache_new ();
  NcmVector *v            = NULL;
  gdouble x               = 2.0;

  g_assert_true (ncm_function_cache_get (cache, &x, &v));
  g_assert_cmpfloat (ncm_vector_get (v, 0), ==, 2.0);
  g_assert_cmpfloat (ncm_vector_get (v, 1), ==, 4.0);

  x = 3.0;
  g_assert_false (ncm_function_cache_get (cache, &x, &v));

  ncm_function_cache_free (cache);
}

static void
test_ncm_function_cache_insert (void)
{
  NcmFunctionCache *cache = _test_function_cache_new ();
  NcmVector *p            = ncm_vector_new (2);
  NcmVector *v            = NULL;
  gdouble x               = 2.0;

  /* An argument already cached keeps its first value */
  ncm_vector_set_all (p, -1.0);
  ncm_function_cache_insert_vector (cache, 2.0, p);
  ncm_function_cache_insert (cache, 2.0, -1.0, -1.0);

  g_assert_true (ncm_function_cache_get (cache, &x, &v));
  g_assert_cmpfloat (ncm_vector_get (v, 0), ==, 2.0);
  g_assert_cmpfloat (ncm_vector_get (v, 1), ==, 4.0);

  ncm_function_cache_insert (cache, 8.0, 8.0, 64.0);
  x = 8.0;
  g_assert_true (ncm_function_cache_get (cache, &x, &v));
  g_assert_cmpfloat (ncm_vector_get (v, 0), ==, 8.0);
  g_assert_cmpfloat (ncm_vector_get (v, 1), ==, 64.0);

  ncm_vector_free (p);
  ncm_function_cache_free (cache);
}

static void
test_ncm_function_cache_get_near (void)
{
  NcmFunctionCache *cache = _test_function_cache_new ();

  /* {x, search type, found, x_c} */
  const struct
  {
    gdouble x;
    NcmFunctionCacheSearchType type;
    gboolean found;
    gdouble x_c;
  } cases[] = {
    {2.0, NCM_FUNCTION_CACHE_SEARCH_BOTH, TRUE, 2.0},
    {2.9, NCM_FUNCTION_CACHE_SEARCH_BOTH, TRUE, 2.0},
    {3.1, NCM_FUNCTION_CACHE_SEARCH_BOTH, TRUE, 4.0},
    {10.0, NCM_FUNCTION_CACHE_SEARCH_BOTH, TRUE, 4.0},
    {2.1, NCM_FUNCTION_CACHE_SEARCH_GT, TRUE, 4.0},
    {2.0, NCM_FUNCTION_CACHE_SEARCH_GT, TRUE, 2.0},
    {3.9, NCM_FUNCTION_CACHE_SEARCH_LT, TRUE, 2.0},
    {1.0, NCM_FUNCTION_CACHE_SEARCH_LT, TRUE, 1.0},
    {4.5, NCM_FUNCTION_CACHE_SEARCH_GT, FALSE, 0.0},
    {0.5, NCM_FUNCTION_CACHE_SEARCH_LT, FALSE, 0.0},

    /* Within reltol = 1e-7 of a cached argument, on the wrong side of it */
    {2.0 * (1.0 + 1.0e-9), NCM_FUNCTION_CACHE_SEARCH_GT, TRUE, 2.0},
    {2.0 * (1.0 - 1.0e-9), NCM_FUNCTION_CACHE_SEARCH_LT, TRUE, 2.0},
    {2.0 * (1.0 + 1.0e-6), NCM_FUNCTION_CACHE_SEARCH_GT, TRUE, 4.0},
  };

  guint i;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    NcmVector *v = NULL;
    gdouble x_c  = -1.0;

    g_assert_true (ncm_function_cache_get_near (cache, cases[i].x, &x_c, &v, cases[i].type) == cases[i].found);

    if (cases[i].found)
    {
      g_assert_cmpfloat (x_c, ==, cases[i].x_c);
      g_assert_cmpfloat (ncm_vector_get (v, 0), ==, cases[i].x_c);
      ncm_vector_free (v);
    }
    else
    {
      g_assert_null (v);
      g_assert_cmpfloat (x_c, ==, -1.0);
    }
  }

  ncm_function_cache_free (cache);
}

static void
test_ncm_function_cache_empty (void)
{
  NcmFunctionCache *cache = _test_function_cache_new ();
  NcmVector *v            = NULL;
  gdouble x               = 2.0;
  gdouble x_c;

  ncm_function_cache_empty_cache (cache);
  g_assert_false (ncm_function_cache_get (cache, &x, &v));
  g_assert_false (ncm_function_cache_get_near (cache, 2.0, &x_c, &v, NCM_FUNCTION_CACHE_SEARCH_BOTH));
  g_assert_null (v);

  ncm_function_cache_insert (cache, 2.0, 5.0, 6.0);
  g_assert_true (ncm_function_cache_get (cache, &x, &v));
  g_assert_cmpfloat (ncm_vector_get (v, 0), ==, 5.0);

  ncm_function_cache_clear (&cache);
  g_assert_null (cache);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/function_cache/tolerances", &test_ncm_function_cache_tolerances);
  g_test_add_func ("/ncm/function_cache/get", &test_ncm_function_cache_get);
  g_test_add_func ("/ncm/function_cache/insert", &test_ncm_function_cache_insert);
  g_test_add_func ("/ncm/function_cache/get_near", &test_ncm_function_cache_get_near);
  g_test_add_func ("/ncm/function_cache/empty", &test_ncm_function_cache_empty);

  g_test_run ();
}

