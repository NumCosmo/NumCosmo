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

void test_ncm_sparam_new (void);
void test_ncm_sparam_free (void);
void test_ncm_sparam_setget_lower_bound (void);
void test_ncm_sparam_setget_upper_bound (void);
void test_ncm_sparam_setget_scale (void);
void test_ncm_sparam_setget_abstol (void);
void test_ncm_sparam_setget_default_value (void);

void test_ncm_sparam_invalid_lower_bound (void);
void test_ncm_sparam_invalid_upper_bound (void);
void test_ncm_sparam_invalid_scale (void);
void test_ncm_sparam_invalid_abstol (void);
void test_ncm_sparam_invalid_default_value (void);

#define _NCM_SPARAM_NTESTS 1000

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/numcosmo/ncm_sparam/new", &test_ncm_sparam_new);

  g_test_add_func ("/numcosmo/ncm_sparam/setget/lower_bound", &test_ncm_sparam_setget_lower_bound);
  g_test_add_func ("/numcosmo/ncm_sparam/setget/upper_bound", &test_ncm_sparam_setget_upper_bound);
  g_test_add_func ("/numcosmo/ncm_sparam/setget/scale", &test_ncm_sparam_setget_scale);
  g_test_add_func ("/numcosmo/ncm_sparam/setget/abstol", &test_ncm_sparam_setget_abstol);
  g_test_add_func ("/numcosmo/ncm_sparam/setget/default_value", &test_ncm_sparam_setget_default_value);

  g_test_add_func ("/numcosmo/ncm_sparam/invalid/lower_bound", &test_ncm_sparam_invalid_lower_bound);
  g_test_add_func ("/numcosmo/ncm_sparam/invalid/upper_bound", &test_ncm_sparam_invalid_upper_bound);
  g_test_add_func ("/numcosmo/ncm_sparam/invalid/scale", &test_ncm_sparam_invalid_scale);
  g_test_add_func ("/numcosmo/ncm_sparam/invalid/abstol", &test_ncm_sparam_invalid_abstol);
  g_test_add_func ("/numcosmo/ncm_sparam/invalid/default_value", &test_ncm_sparam_invalid_default_value);

  g_test_add_func ("/numcosmo/ncm_sparam/free", &test_ncm_sparam_free);

  g_test_run ();
}

static NcmSParam *p;

void
test_ncm_sparam_new (void)
{
  gchar *name = "t1";
  gchar *symbol = "t2^2";
  gdouble lower_bound = -10.0;
  gdouble upper_bound = M_PI;
  gdouble scale = 1e-3;
  gdouble abstol = 1e-5;
  gdouble default_val = (upper_bound + lower_bound) * 0.5;
  NcmParamType ftype = NCM_PARAM_TYPE_FIXED;

  p = ncm_sparam_new (name, symbol, lower_bound, upper_bound, scale, abstol, default_val, ftype);
  g_assert (p != NULL);
  g_assert (NCM_IS_SPARAM (p));

  g_assert (p->name != name);
  g_assert (p->symbol != symbol);

  g_assert_cmpstr (p->name, ==, name);
  g_assert_cmpstr (p->symbol, ==, symbol);
  ncm_assert_cmpdouble (p->lower_bound, ==, lower_bound);
  ncm_assert_cmpdouble (p->upper_bound, ==, upper_bound);
  ncm_assert_cmpdouble (p->scale, ==, scale);
  ncm_assert_cmpdouble (p->abstol, ==, abstol);
}


void
test_ncm_sparam_free (void)
{
  ncm_sparam_free (p);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_sparam_free (p);
    exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_sparam_setget_lower_bound (void)
{
  guint ntests = _NCM_SPARAM_NTESTS;

  g_assert (p != NULL);
  g_assert (NCM_IS_SPARAM (p));

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
test_ncm_sparam_setget_upper_bound (void)
{
  guint ntests = _NCM_SPARAM_NTESTS;

  g_assert (p != NULL);
  g_assert (NCM_IS_SPARAM (p));

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
test_ncm_sparam_setget_scale (void)
{
  guint ntests = _NCM_SPARAM_NTESTS;

  g_assert (p != NULL);
  g_assert (NCM_IS_SPARAM (p));

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
test_ncm_sparam_setget_abstol (void)
{
  guint ntests = _NCM_SPARAM_NTESTS;

  g_assert (p != NULL);
  g_assert (NCM_IS_SPARAM (p));

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
test_ncm_sparam_setget_default_value (void)
{
  guint ntests = _NCM_SPARAM_NTESTS;

  g_assert (p != NULL);
  g_assert (NCM_IS_SPARAM (p));

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
test_ncm_sparam_invalid_lower_bound (void)
{
  const gdouble ub = ncm_sparam_get_upper_bound (p);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_sparam_set_lower_bound (p, (ub > 0) ? (ub * 1.1) : (ub * 0.9));
    exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_sparam_invalid_upper_bound (void)
{
  const gdouble lb = ncm_sparam_get_lower_bound (p);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_sparam_set_upper_bound (p, (lb > 0) ? (lb * 0.9) : (lb * 1.1));
    exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_sparam_invalid_scale (void)
{
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_sparam_set_scale (p, -fabs (g_test_rand_double ()));
    exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_sparam_invalid_abstol (void)
{
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_sparam_set_absolute_tolerance (p, -fabs (g_test_rand_double ()));
    exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_sparam_invalid_default_value (void)
{
}
