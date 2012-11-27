/***************************************************************************
 *            test_ncm_vector.c
 *
 *  Tue April 03 16:02:26 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

void test_ncm_vector_new (void);
void test_ncm_vector_new_gsl (void);
void test_ncm_vector_new_array (void);
void test_ncm_vector_new_data_slice (void);
void test_ncm_vector_new_data_malloc (void);
void test_ncm_vector_new_data_static (void);
void test_ncm_vector_new_data_const (void);
void test_ncm_vector_operations (void);
void test_ncm_vector_free (void);
void test_ncm_vector_subvector (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/numcosmo/ncm_vector/new", &test_ncm_vector_new);
  g_test_add_func ("/numcosmo/ncm_vector/new_gsl", &test_ncm_vector_new_gsl);
  g_test_add_func ("/numcosmo/ncm_vector/new_array", &test_ncm_vector_new_array);
  g_test_add_func ("/numcosmo/ncm_vector/new_data_slice", &test_ncm_vector_new_data_slice);
  g_test_add_func ("/numcosmo/ncm_vector/new_data_malloc", &test_ncm_vector_new_data_malloc);
  g_test_add_func ("/numcosmo/ncm_vector/new_data_static", &test_ncm_vector_new_data_static);
  g_test_add_func ("/numcosmo/ncm_vector/new_data_const", &test_ncm_vector_new_data_const);
  g_test_add_func ("/numcosmo/ncm_vector/operations", &test_ncm_vector_operations);
  g_test_add_func ("/numcosmo/ncm_vector/subvector", &test_ncm_vector_subvector);
  g_test_add_func ("/numcosmo/ncm_vector/free", &test_ncm_vector_free);

  g_test_run ();
}

static NcmVector *v;

#define _NCM_VECTOR_TEST_SIZE 100


void
test_ncm_vector_new_sanity (NcmVector *vv)
{
  guint i;

  g_assert (NCM_IS_VECTOR (vv));
  g_assert (g_object_is_floating (vv));
  ncm_vector_ref (vv);
  g_assert (!g_object_is_floating (vv));

  for (i = 0; i < 10 * _NCM_VECTOR_TEST_SIZE; i++)
  {
    const guint n = g_test_rand_int_range (0, _NCM_VECTOR_TEST_SIZE - 1);
    const gdouble d = g_test_rand_double ();
    ncm_vector_set (vv, n, d);
    g_assert_cmpfloat (ncm_vector_get (vv, n), ==, d);
  }
}

void
test_ncm_vector_new (void)
{
  v = ncm_vector_new (_NCM_VECTOR_TEST_SIZE);
  test_ncm_vector_new_sanity (v);
}

void
test_ncm_vector_operations (void)
{
  NcmVector *cv = ncm_vector_dup (v);
  guint i;
  test_ncm_vector_new_sanity (cv);

  for (i = 0; i < 10 * _NCM_VECTOR_TEST_SIZE; i++)
  {
    const guint n = g_test_rand_int_range (0, _NCM_VECTOR_TEST_SIZE - 1);
    const gdouble d = g_test_rand_double ();
    ncm_vector_set (v, n, d);
    g_assert_cmpfloat (ncm_vector_get (v, n), ==, d);
  }

  for (i = 0; i < 10 * _NCM_VECTOR_TEST_SIZE; i++)
  {
    guint n;
    gdouble *d;

    n = g_test_rand_int_range (0, _NCM_VECTOR_TEST_SIZE - 1);
    d = ncm_vector_ptr (v, n);
    (*d) *= g_test_rand_double ();

    g_assert_cmpfloat (ncm_vector_get (v, n), ==, *d);
  }

  for (i = 0; i < 10 * _NCM_VECTOR_TEST_SIZE; i++)
  {
    guint n;
    gdouble d, d1 = g_test_rand_double ();

    n = g_test_rand_int_range (0, _NCM_VECTOR_TEST_SIZE - 1);
    d = ncm_vector_get (v, n);
    d += d1;

    ncm_vector_addto (v, n, d1);
    g_assert_cmpfloat (ncm_vector_get (v, n), ==, d);
  }

  for (i = 0; i < 10 * _NCM_VECTOR_TEST_SIZE; i++)
  {
    guint n;
    gdouble d, d1 = g_test_rand_double ();

    n = g_test_rand_int_range (0, _NCM_VECTOR_TEST_SIZE - 1);
    d = ncm_vector_get (v, n);
    d -= d1;

    ncm_vector_subfrom (v, n, d1);
    g_assert_cmpfloat (ncm_vector_get (v, n), ==, d);
  }

  for (i = 0; i < 10 * _NCM_VECTOR_TEST_SIZE; i++)
  {
    guint n;
    gdouble d1 = g_test_rand_double ();

    n = g_test_rand_int_range (0, _NCM_VECTOR_TEST_SIZE - 1);
    ncm_vector_set_all (v, d1);
    g_assert_cmpfloat (ncm_vector_get (v, n), ==, d1);
  }

  for (i = 0; i < 10 * _NCM_VECTOR_TEST_SIZE; i++)
  {
    guint n;
    gdouble d, d1 = g_test_rand_double ();

    n = g_test_rand_int_range (0, _NCM_VECTOR_TEST_SIZE - 1);
    d = ncm_vector_get (v, n);
    ncm_vector_scale (v, d1);
    g_assert_cmpfloat (ncm_vector_get (v, n), ==, d * d1);
  }

  for (i = 0; i < 10 * _NCM_VECTOR_TEST_SIZE; i++)
  {
    guint n;
    gdouble d1 = g_test_rand_double (),
    d2 = g_test_rand_double ();

    n = g_test_rand_int_range (0, _NCM_VECTOR_TEST_SIZE - 1);

    ncm_vector_set_all (v, d1);
    ncm_vector_set_all (cv, d2);

    ncm_vector_div (v, cv);

    g_assert_cmpfloat (ncm_vector_get (v, n), ==, d1 / d2);
  }

  for (i = 0; i < 10 * _NCM_VECTOR_TEST_SIZE; i++)
  {
    guint n;
    gdouble d1 = g_test_rand_double (),
    d2 = g_test_rand_double ();

    n = g_test_rand_int_range (0, _NCM_VECTOR_TEST_SIZE - 1);

    ncm_vector_set_all (v, d1);
    ncm_vector_set_all (cv, d2);

    ncm_vector_add (v, cv);

    g_assert_cmpfloat (ncm_vector_get (v, n), ==, d1 + d2);
  }

  for (i = 0; i < 10 * _NCM_VECTOR_TEST_SIZE; i++)
  {
    guint n;
    gdouble d1 = g_test_rand_double (),
    d2 = g_test_rand_double ();

    n = g_test_rand_int_range (0, _NCM_VECTOR_TEST_SIZE - 1);

    ncm_vector_set_all (v, d1);
    ncm_vector_set_all (cv, d2);

    ncm_vector_sub (v, cv);

    g_assert_cmpfloat (ncm_vector_get (v, n), ==, d1 - d2);
  }

  ncm_vector_set_zero (v);
  for (i = 0; i < _NCM_VECTOR_TEST_SIZE; i++)
    g_assert_cmpfloat (ncm_vector_get (v, i), ==, 0.0);

  ncm_vector_memcpy (v, cv);
  for (i = 0; i < _NCM_VECTOR_TEST_SIZE; i++)
    g_assert_cmpfloat (ncm_vector_get (v, i), ==, ncm_vector_get (cv, i));

  ncm_vector_set_zero (v);
  ncm_vector_memcpy2 (v, cv, 2, 0, _NCM_VECTOR_TEST_SIZE - 2);

  g_assert_cmpfloat (ncm_vector_get (v, 0), ==, 0.0);
  g_assert_cmpfloat (ncm_vector_get (v, 1), ==, 0.0);
  for (i = 0; i < _NCM_VECTOR_TEST_SIZE - 2; i++)
    g_assert_cmpfloat (ncm_vector_get (v, i + 2), ==, ncm_vector_get (cv, i));
}

void
test_ncm_vector_new_gsl (void)
{
  gsl_vector *gv = gsl_vector_alloc (_NCM_VECTOR_TEST_SIZE);
  NcmVector *vv = ncm_vector_new_gsl (gv);
  test_ncm_vector_new_sanity (vv);
  g_assert_cmpfloat (ncm_vector_len (vv), ==, gv->size);

  ncm_vector_free (vv);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_vector_free (vv);
    exit (0);
  }
  g_test_trap_assert_failed ();

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    g_assert (vv->gv == NULL);
    exit (0);
  }
  g_test_trap_assert_passed ();
}

void
test_ncm_vector_new_array (void)
{
  GArray *ga = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), _NCM_VECTOR_TEST_SIZE);
  NcmVector *vv;
  g_array_set_size (ga, _NCM_VECTOR_TEST_SIZE);
  vv = ncm_vector_new_array (ga);
  test_ncm_vector_new_sanity (vv);
  g_array_unref (ga);

  g_assert_cmpuint (ncm_vector_len (vv), ==, ga->len);

  g_assert (ga == ncm_vector_get_array (vv));
  g_array_unref (ga);

  ncm_vector_free (vv);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_vector_free (vv);
    exit (0);
  }
  g_test_trap_assert_failed ();

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    g_array_unref (vv->a);
    exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_vector_new_data_slice (void)
{
  gdouble *d = g_slice_alloc (_NCM_VECTOR_TEST_SIZE * sizeof (gdouble));
  NcmVector *vv;
  vv = ncm_vector_new_data_slice (d, _NCM_VECTOR_TEST_SIZE, 1);
  test_ncm_vector_new_sanity (vv);

  g_assert_cmpuint (ncm_vector_len (vv), ==, _NCM_VECTOR_TEST_SIZE);

  ncm_vector_free (vv);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_vector_free (vv);
    exit (0);
  }
  g_test_trap_assert_failed ();

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    g_assert (NCM_VECTOR_DATA (vv) == NULL);
    exit (0);
  }
  g_test_trap_assert_passed ();
}

void
test_ncm_vector_new_data_malloc (void)
{
  gdouble *d = g_malloc (_NCM_VECTOR_TEST_SIZE * sizeof (gdouble));
  NcmVector *vv;
  vv = ncm_vector_new_data_malloc (d, _NCM_VECTOR_TEST_SIZE, 1);
  test_ncm_vector_new_sanity (vv);

  g_assert_cmpuint (ncm_vector_len (vv), ==, _NCM_VECTOR_TEST_SIZE);

  ncm_vector_free (vv);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_vector_free (vv);
    exit (0);
  }
  g_test_trap_assert_failed ();

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    g_assert (NCM_VECTOR_DATA (vv) == NULL);
    exit (0);
  }
  g_test_trap_assert_passed ();
}

void
test_ncm_vector_new_data_static (void)
{
  gdouble d[_NCM_VECTOR_TEST_SIZE];
  NcmVector *vv;
  vv = ncm_vector_new_data_static (d, _NCM_VECTOR_TEST_SIZE, 1);
  test_ncm_vector_new_sanity (vv);

  g_assert_cmpuint (ncm_vector_len (vv), ==, _NCM_VECTOR_TEST_SIZE);

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_vector_free (vv);
    exit (0);
  }
  g_test_trap_assert_passed ();
}

void
test_ncm_vector_new_data_const (void)
{
  const gdouble d[_NCM_VECTOR_TEST_SIZE] = {0.0, };
  const NcmVector *vv;
  vv = ncm_vector_new_data_const (d, _NCM_VECTOR_TEST_SIZE, 1);

  g_assert_cmpuint (ncm_vector_len (vv), ==, _NCM_VECTOR_TEST_SIZE);

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_vector_const_free (vv);
    exit (0);
  }
  g_test_trap_assert_passed ();
}

void
test_ncm_vector_subvector (void)
{
  NcmVector *sv = ncm_vector_get_subvector (v, 1, _NCM_VECTOR_TEST_SIZE - 1);
  guint ntests = 20 * _NCM_VECTOR_TEST_SIZE;

  g_assert_cmpuint (ncm_vector_len (sv), ==, (_NCM_VECTOR_TEST_SIZE - 1));

  while (ntests--)
  {
    guint i = g_test_rand_int_range (0, _NCM_VECTOR_TEST_SIZE - 1);
    ncm_vector_set (sv, i, g_test_rand_double ());
    g_assert_cmpfloat (ncm_vector_get (sv, i), ==, ncm_vector_get (v, i + 1));
  }

  g_assert (NCM_IS_VECTOR (sv));
  g_assert (g_object_is_floating (sv));
  ncm_vector_ref (sv);
  g_assert (!g_object_is_floating (sv));

  ncm_vector_free (v);
  g_assert (G_IS_OBJECT (v));
  ncm_vector_ref (v);

  ncm_vector_free (sv);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_vector_free (sv);
    exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_vector_free (void)
{
  ncm_vector_free (v);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_vector_free (v);
    exit (0);
  }
  g_test_trap_assert_failed ();
}
