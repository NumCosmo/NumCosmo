/***************************************************************************
 *            test_ncm_dtuple.c
 *
 *  Fri December 08 10:28:16 2023
 *  Copyright  2023  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_dtuple.c
 * Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

typedef struct _TestNcmDTuple2
{
  NcmDTuple2 *dt2;
} TestNcmDTuple2;

void test_ncm_dtuple2_new (TestNcmDTuple2 *test, gconstpointer pdata);
void test_ncm_dtuple2_free (TestNcmDTuple2 *test, gconstpointer pdata);

void test_ncm_dtuple2_copy (TestNcmDTuple2 *test, gconstpointer pdata);
void test_ncm_dtuple2_static (TestNcmDTuple2 *test, gconstpointer pdata);

typedef struct _TestNcmDTuple3
{
  NcmDTuple3 *dt3;
} TestNcmDTuple3;

void test_ncm_dtuple3_new (TestNcmDTuple3 *test, gconstpointer pdata);
void test_ncm_dtuple3_free (TestNcmDTuple3 *test, gconstpointer pdata);

void test_ncm_dtuple3_copy (TestNcmDTuple3 *test, gconstpointer pdata);
void test_ncm_dtuple3_static (TestNcmDTuple3 *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/dtuple2/copy", TestNcmDTuple2, NULL,
              &test_ncm_dtuple2_new,
              &test_ncm_dtuple2_copy,
              &test_ncm_dtuple2_free);

  g_test_add ("/ncm/dtuple2/static", TestNcmDTuple2, NULL,
              &test_ncm_dtuple2_new,
              &test_ncm_dtuple2_static,
              &test_ncm_dtuple2_free);

  g_test_add ("/ncm/dtuple3/copy", TestNcmDTuple3, NULL,
              &test_ncm_dtuple3_new,
              &test_ncm_dtuple3_copy,
              &test_ncm_dtuple3_free);

  g_test_add ("/ncm/dtuple3/static", TestNcmDTuple3, NULL,
              &test_ncm_dtuple3_new,
              &test_ncm_dtuple3_static,
              &test_ncm_dtuple3_free);

  g_test_run ();
}

void
test_ncm_dtuple2_new (TestNcmDTuple2 *test, gconstpointer pdata)
{
  const gdouble elem0 = g_test_rand_double ();
  const gdouble elem1 = g_test_rand_double ();

  test->dt2 = ncm_dtuple2_new (elem0, elem1);

  g_assert_cmpfloat (test->dt2->elements[0], ==, elem0);
  g_assert_cmpfloat (test->dt2->elements[1], ==, elem1);
}

void
test_ncm_dtuple2_free (TestNcmDTuple2 *test, gconstpointer pdata)
{
  ncm_dtuple2_free (test->dt2);
}

void
test_ncm_dtuple2_copy (TestNcmDTuple2 *test, gconstpointer pdata)
{
  NcmDTuple2 *dt2_copy = ncm_dtuple2_copy (test->dt2);

  g_assert_cmpfloat (test->dt2->elements[0], ==, dt2_copy->elements[0]);
  g_assert_cmpfloat (test->dt2->elements[1], ==, dt2_copy->elements[1]);

  ncm_dtuple2_free (dt2_copy);
}

void
test_ncm_dtuple2_static (TestNcmDTuple2 *test, gconstpointer pdata)
{
  NcmDTuple2 dt2 = NCM_DTUPLE2_STATIC_INIT (test->dt2->elements[0], test->dt2->elements[1]);

  g_assert_cmpfloat (dt2.elements[0], ==, test->dt2->elements[0]);
  g_assert_cmpfloat (dt2.elements[1], ==, test->dt2->elements[1]);
}

void
test_ncm_dtuple3_new (TestNcmDTuple3 *test, gconstpointer pdata)
{
  const gdouble elem0 = g_test_rand_double ();
  const gdouble elem1 = g_test_rand_double ();
  const gdouble elem2 = g_test_rand_double ();

  test->dt3 = ncm_dtuple3_new (elem0, elem1, elem2);

  g_assert_cmpfloat (test->dt3->elements[0], ==, elem0);
  g_assert_cmpfloat (test->dt3->elements[1], ==, elem1);
  g_assert_cmpfloat (test->dt3->elements[2], ==, elem2);
}

void
test_ncm_dtuple3_free (TestNcmDTuple3 *test, gconstpointer pdata)
{
  ncm_dtuple3_free (test->dt3);
}

void
test_ncm_dtuple3_copy (TestNcmDTuple3 *test, gconstpointer pdata)
{
  NcmDTuple3 *dt3_copy = ncm_dtuple3_copy (test->dt3);

  g_assert_cmpfloat (test->dt3->elements[0], ==, dt3_copy->elements[0]);
  g_assert_cmpfloat (test->dt3->elements[1], ==, dt3_copy->elements[1]);
  g_assert_cmpfloat (test->dt3->elements[2], ==, dt3_copy->elements[2]);

  ncm_dtuple3_free (dt3_copy);
}

void
test_ncm_dtuple3_static (TestNcmDTuple3 *test, gconstpointer pdata)
{
  NcmDTuple3 dt3 = NCM_DTUPLE3_STATIC_INIT (test->dt3->elements[0], test->dt3->elements[1], test->dt3->elements[2]);

  g_assert_cmpfloat (dt3.elements[0], ==, test->dt3->elements[0]);
  g_assert_cmpfloat (dt3.elements[1], ==, test->dt3->elements[1]);
  g_assert_cmpfloat (dt3.elements[2], ==, test->dt3->elements[2]);
}

