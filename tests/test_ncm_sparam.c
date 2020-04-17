/***************************************************************************
 *            ncm_sparam.c
 *
 *  Fri March 23 14:26:42 2012
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

typedef struct _TestNcmSparam
{
  NcmSParam *p;
  guint ntests;
} TestNcmSparam;

void test_ncm_sparam_new (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_sparam_free (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_sparam_setget_lower_bound (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_sparam_setget_upper_bound (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_sparam_setget_scale (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_sparam_setget_abstol (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_sparam_setget_default_value (TestNcmSparam *test, gconstpointer pdata);

void test_ncm_sparam_traps (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_sparam_invalid_lower_bound (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_sparam_invalid_upper_bound (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_sparam_invalid_scale (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_sparam_invalid_abstol (TestNcmSparam *test, gconstpointer pdata);
void test_ncm_sparam_invalid_default_value (TestNcmSparam *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/sparam/setget/lower_bound", TestNcmSparam, NULL,
              &test_ncm_sparam_new,
              &test_ncm_sparam_setget_lower_bound,
              &test_ncm_sparam_free);

  g_test_add ("/ncm/sparam/setget/upper_bound", TestNcmSparam, NULL,
              &test_ncm_sparam_new,
              &test_ncm_sparam_setget_upper_bound,
              &test_ncm_sparam_free);

  g_test_add ("/ncm/sparam/setget/scale", TestNcmSparam, NULL,
              &test_ncm_sparam_new,
              &test_ncm_sparam_setget_scale,
              &test_ncm_sparam_free);

  g_test_add ("/ncm/sparam/setget/abstol", TestNcmSparam, NULL,
              &test_ncm_sparam_new,
              &test_ncm_sparam_setget_abstol,
              &test_ncm_sparam_free);

  g_test_add ("/ncm/sparam/setget/default_value", TestNcmSparam, NULL,
              &test_ncm_sparam_new,
              &test_ncm_sparam_setget_default_value,
              &test_ncm_sparam_free);

  g_test_add ("/ncm/sparam/traps", TestNcmSparam, NULL,
              &test_ncm_sparam_new,
              &test_ncm_sparam_traps,
              &test_ncm_sparam_free);
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_add ("/ncm/sparam/invalid/lower_bound/subprocess", TestNcmSparam, NULL,
              &test_ncm_sparam_new,
              &test_ncm_sparam_invalid_lower_bound,
              &test_ncm_sparam_free);

  g_test_add ("/ncm/sparam/invalid/upper_bound/subprocess", TestNcmSparam, NULL,
              &test_ncm_sparam_new,
              &test_ncm_sparam_invalid_upper_bound,
              &test_ncm_sparam_free);

  g_test_add ("/ncm/sparam/invalid/scale/subprocess", TestNcmSparam, NULL,
              &test_ncm_sparam_new,
              &test_ncm_sparam_invalid_scale,
              &test_ncm_sparam_free);

  g_test_add ("/ncm/sparam/invalid/abstol/subprocess", TestNcmSparam, NULL,
              &test_ncm_sparam_new,
              &test_ncm_sparam_invalid_abstol,
              &test_ncm_sparam_free);

  g_test_add ("/ncm/sparam/invalid/default_value/subprocess", TestNcmSparam, NULL,
              &test_ncm_sparam_new,
              &test_ncm_sparam_invalid_default_value,
              &test_ncm_sparam_free);
#endif
  g_test_run ();
}

void
test_ncm_sparam_new (TestNcmSparam *test, gconstpointer pdata)
{
  gchar *name = "t1";
  gchar *symbol = "t2^2";
  gdouble lower_bound = -10.0;
  gdouble upper_bound = M_PI;
  gdouble scale = 1e-3;
  gdouble abstol = 1e-5;
  gdouble default_val = (upper_bound + lower_bound) * 0.5;
  NcmParamType ftype = NCM_PARAM_TYPE_FIXED;
  NcmSParam *p;

  test->ntests = 1000;
  p = test->p = ncm_sparam_new (name, symbol, lower_bound, upper_bound, scale, abstol, default_val, ftype);

  g_assert_true (p != NULL);
  g_assert_true (NCM_IS_SPARAM (p));

  g_assert_true (p->name != name);
  g_assert_true (p->symbol != symbol);

  g_assert_cmpstr (p->name, ==, name);
  g_assert_cmpstr (p->symbol, ==, symbol);

  ncm_assert_cmpdouble (p->lower_bound, ==, lower_bound);
  ncm_assert_cmpdouble (p->upper_bound, ==, upper_bound);
  ncm_assert_cmpdouble (p->scale, ==, scale);
  ncm_assert_cmpdouble (p->abstol, ==, abstol);
}

void
test_ncm_sparam_free (TestNcmSparam *test, gconstpointer pdata)
{
  NcmSParam *p = test->p;
  NCM_TEST_FREE (ncm_sparam_free, p);
}

void
test_ncm_sparam_setget_lower_bound (TestNcmSparam *test, gconstpointer pdata)
{
  NcmSParam *p = test->p;
  guint ntests = test->ntests;

  g_assert_true (p != NULL);
  g_assert_true (NCM_IS_SPARAM (p));

  while (ntests--)
  {
    const gdouble lb = g_test_rand_double_range (-G_MAXDOUBLE, ncm_sparam_get_upper_bound (p));
    const gdouble old_lb = ncm_sparam_get_lower_bound (p);
    gdouble new_lb = 0.0;

    ncm_assert_cmpdouble (ncm_sparam_get_lower_bound (p), ==, old_lb);

    ncm_sparam_set_lower_bound (p, lb);
    ncm_assert_cmpdouble (ncm_sparam_get_lower_bound (p), ==, lb);

    g_object_get (p, "lower-bound", &new_lb, NULL);
    ncm_assert_cmpdouble (new_lb, ==, lb);

    g_object_set (p, "lower-bound", old_lb, NULL);
    ncm_assert_cmpdouble (ncm_sparam_get_lower_bound (p), ==, old_lb);
  }
}

void
test_ncm_sparam_setget_upper_bound (TestNcmSparam *test, gconstpointer pdata)
{
  NcmSParam *p = test->p;
  guint ntests = test->ntests;

  g_assert_true (p != NULL);
  g_assert_true (NCM_IS_SPARAM (p));

  while (ntests--)
  {
    const gdouble ub = g_test_rand_double_range (ncm_sparam_get_lower_bound (p), G_MAXDOUBLE);
    const gdouble old_ub = ncm_sparam_get_upper_bound (p);
    gdouble new_ub = 0.0;

    ncm_assert_cmpdouble (ncm_sparam_get_upper_bound (p), ==, old_ub);
    ncm_sparam_set_upper_bound (p, ub);
    ncm_assert_cmpdouble (ncm_sparam_get_upper_bound (p), ==, ub);

    g_object_get (p, "upper-bound", &new_ub, NULL);
    ncm_assert_cmpdouble (new_ub, ==, ub);

    g_object_set (p, "upper-bound", old_ub, NULL);
    ncm_assert_cmpdouble (ncm_sparam_get_upper_bound (p), ==, old_ub);
  }
}

void
test_ncm_sparam_setget_scale (TestNcmSparam *test, gconstpointer pdata)
{
  NcmSParam *p = test->p;
  guint ntests = test->ntests;

  g_assert_true (p != NULL);
  g_assert_true (NCM_IS_SPARAM (p));

  while (ntests--)
  {
    const gdouble scale = fabs (g_test_rand_double ());
    const gdouble old_scale = ncm_sparam_get_scale (p);
    gdouble new_scale = 0.0;

    ncm_assert_cmpdouble (ncm_sparam_get_scale (p), ==, old_scale);
    ncm_sparam_set_scale (p, scale);
    ncm_assert_cmpdouble (ncm_sparam_get_scale (p), ==, scale);

    g_object_get (p, "scale", &new_scale, NULL);
    ncm_assert_cmpdouble (new_scale, ==, scale);

    g_object_set (p, "scale", old_scale, NULL);
    ncm_assert_cmpdouble (ncm_sparam_get_scale (p), ==, old_scale);
  }
}

void
test_ncm_sparam_setget_abstol (TestNcmSparam *test, gconstpointer pdata)
{
  NcmSParam *p = test->p;
  guint ntests = test->ntests;

  g_assert_true (p != NULL);
  g_assert_true (NCM_IS_SPARAM (p));

  while (ntests--)
  {
    const gdouble abstol = fabs (g_test_rand_double ());
    const gdouble old_abstol = ncm_sparam_get_absolute_tolerance (p);
    gdouble new_abstol = 0.0;

    ncm_assert_cmpdouble (ncm_sparam_get_absolute_tolerance (p), ==, old_abstol);
    ncm_sparam_set_absolute_tolerance (p, abstol);
    ncm_assert_cmpdouble (ncm_sparam_get_absolute_tolerance (p), ==, abstol);

    g_object_get (p, "absolute-tolerance", &new_abstol, NULL);
    ncm_assert_cmpdouble (new_abstol, ==, abstol);

    g_object_set (p, "absolute-tolerance", old_abstol, NULL);
    ncm_assert_cmpdouble (ncm_sparam_get_absolute_tolerance (p), ==, old_abstol);
  }
}

void
test_ncm_sparam_setget_default_value (TestNcmSparam *test, gconstpointer pdata)
{
  NcmSParam *p = test->p;
  guint ntests = test->ntests;

  g_assert_true (p != NULL);
  g_assert_true (NCM_IS_SPARAM (p));

  while (ntests--)
  {
    const gdouble default_value = g_test_rand_double ();
    const gdouble old_default_value = ncm_sparam_get_default_value (p);
    gdouble new_default_value = 0.0;

    ncm_assert_cmpdouble (ncm_sparam_get_default_value (p), ==, old_default_value);
    ncm_sparam_set_default_value (p, default_value);
    ncm_assert_cmpdouble (ncm_sparam_get_default_value (p), ==, default_value);

    g_object_get (p, "default-value", &new_default_value, NULL);
    ncm_assert_cmpdouble (new_default_value, ==, default_value);

    g_object_set (p, "default-value", old_default_value, NULL);
    ncm_assert_cmpdouble (ncm_sparam_get_default_value (p), ==, old_default_value);
  }
}

void
test_ncm_sparam_traps (TestNcmSparam *test, gconstpointer pdata)
{
#if GLIB_CHECK_VERSION(2,38,0)
  g_test_trap_subprocess ("/ncm/sparam/invalid/lower_bound/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/sparam/invalid/upper_bound/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/sparam/invalid/scale/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/sparam/invalid/abstol/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/sparam/invalid/default_value/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}


void
test_ncm_sparam_invalid_lower_bound (TestNcmSparam *test, gconstpointer pdata)
{
  NcmSParam *p = test->p;
  const gdouble ub = ncm_sparam_get_upper_bound (p);
  ncm_sparam_set_lower_bound (p, (ub > 0) ? (ub * 1.1) : (ub * 0.9));
}

void
test_ncm_sparam_invalid_upper_bound (TestNcmSparam *test, gconstpointer pdata)
{
  NcmSParam *p = test->p;
  const gdouble lb = ncm_sparam_get_lower_bound (p);
  ncm_sparam_set_upper_bound (p, (lb > 0) ? (lb * 0.9) : (lb * 1.1));
}

void
test_ncm_sparam_invalid_scale (TestNcmSparam *test, gconstpointer pdata)
{
  NcmSParam *p = test->p;
  ncm_sparam_set_scale (p, -fabs (g_test_rand_double ()));
}

void
test_ncm_sparam_invalid_abstol (TestNcmSparam *test, gconstpointer pdata)
{
  NcmSParam *p = test->p;
  ncm_sparam_set_absolute_tolerance (p, -fabs (g_test_rand_double ()));
}

void
test_ncm_sparam_invalid_default_value (TestNcmSparam *test, gconstpointer pdata)
{
  g_assert_not_reached ();
}
