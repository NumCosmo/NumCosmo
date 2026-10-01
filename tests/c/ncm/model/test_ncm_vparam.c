/***************************************************************************
 *            test_ncm_vparam.c
 *
 *  Mon September 28 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

static NcmVParam *
_test_vparam_new (const guint len)
{
  return ncm_vparam_full_new (len, "w", "w", -1.0, 1.0, 0.1, 0.0, 0.5, NCM_PARAM_TYPE_FIXED);
}

/* Components copy the default, named name_i and {symbol}_i; set_len grows and shrinks. */
static void
test_ncm_vparam_components (void)
{
  NcmVParam *vp = _test_vparam_new (3);
  guint i, default_len;

  g_assert_cmpuint (ncm_vparam_get_len (vp), ==, 3);
  g_assert_cmpuint (ncm_vparam_len (vp), ==, 3);

  for (i = 0; i < 3; i++)
  {
    gchar *name   = g_strdup_printf ("w_%u", i);
    gchar *symbol = g_strdup_printf ("{w}_%u", i);

    g_assert_cmpstr (ncm_sparam_name (ncm_vparam_peek_sparam (vp, i)), ==, name);
    g_assert_cmpstr (ncm_sparam_symbol (ncm_vparam_peek_sparam (vp, i)), ==, symbol);
    g_assert_cmpfloat (ncm_vparam_get_lower_bound (vp, i), ==, -1.0);
    g_assert_cmpfloat (ncm_vparam_get_upper_bound (vp, i), ==, 1.0);
    g_assert_cmpfloat (ncm_vparam_get_scale (vp, i), ==, 0.1);
    g_assert_cmpfloat (ncm_vparam_get_absolute_tolerance (vp, i), ==, 0.0);
    g_assert_cmpfloat (ncm_vparam_get_default_value (vp, i), ==, 0.5);
    g_assert_cmpint (ncm_vparam_get_fit_type (vp, i), ==, NCM_PARAM_TYPE_FIXED);

    g_free (name);
    g_free (symbol);
  }

  ncm_vparam_set_len (vp, 5);
  g_assert_cmpuint (ncm_vparam_get_len (vp), ==, 5);
  g_assert_cmpstr (ncm_sparam_name (ncm_vparam_peek_sparam (vp, 4)), ==, "w_4");

  ncm_vparam_set_len (vp, 2);
  g_assert_cmpuint (ncm_vparam_get_len (vp), ==, 2);
  g_object_get (vp, "default-len", &default_len, NULL);
  g_assert_cmpuint (default_len, ==, 2);

  ncm_vparam_free (vp);
}

/* The per-component setters change that component only. */
static void
test_ncm_vparam_setters (void)
{
  NcmVParam *vp = _test_vparam_new (2);

  ncm_vparam_set_lower_bound (vp, 1, -2.0);
  ncm_vparam_set_upper_bound (vp, 1, 2.0);
  ncm_vparam_set_scale (vp, 1, 0.3);
  ncm_vparam_set_absolute_tolerance (vp, 1, 1.0e-4);
  ncm_vparam_set_default_value (vp, 1, 0.25);
  ncm_vparam_set_fit_type (vp, 1, NCM_PARAM_TYPE_FREE);

  g_assert_cmpfloat (ncm_vparam_get_lower_bound (vp, 1), ==, -2.0);
  g_assert_cmpfloat (ncm_vparam_get_upper_bound (vp, 1), ==, 2.0);
  g_assert_cmpfloat (ncm_vparam_get_scale (vp, 1), ==, 0.3);
  g_assert_cmpfloat (ncm_vparam_get_absolute_tolerance (vp, 1), ==, 1.0e-4);
  g_assert_cmpfloat (ncm_vparam_get_default_value (vp, 1), ==, 0.25);
  g_assert_cmpint (ncm_vparam_get_fit_type (vp, 1), ==, NCM_PARAM_TYPE_FREE);

  g_assert_cmpfloat (ncm_vparam_get_lower_bound (vp, 0), ==, -1.0);
  g_assert_cmpfloat (ncm_vparam_get_scale (vp, 0), ==, 0.1);
  g_assert_cmpint (ncm_vparam_get_fit_type (vp, 0), ==, NCM_PARAM_TYPE_FIXED);

  ncm_vparam_free (vp);
}

/*
 * set_sparam takes a reference, also to the component already held (released
 * before the fix, before being referenced); set_sparam_full builds a new one; a
 * copy is independent of the original.
 */
static void
test_ncm_vparam_set_sparam_copy (void)
{
  NcmVParam *vp = _test_vparam_new (2);
  NcmSParam *sp = ncm_sparam_new ("x", "x", 0.0, 2.0, 0.5, 0.0, 1.0, NCM_PARAM_TYPE_FREE);
  NcmVParam *cp;
  NcmSParam *got;

  ncm_vparam_set_sparam (vp, 0, ncm_vparam_peek_sparam (vp, 0));
  g_assert_cmpstr (ncm_sparam_name (ncm_vparam_peek_sparam (vp, 0)), ==, "w_0");

  ncm_vparam_set_sparam (vp, 0, sp);
  ncm_sparam_free (sp);
  g_assert_cmpstr (ncm_sparam_name (ncm_vparam_peek_sparam (vp, 0)), ==, "x");

  ncm_vparam_set_sparam_full (vp, 1, "y", "y", -3.0, 3.0, 0.2, 0.0, 0.0, NCM_PARAM_TYPE_FREE);
  g_assert_cmpstr (ncm_sparam_name (ncm_vparam_peek_sparam (vp, 1)), ==, "y");
  g_assert_cmpfloat (ncm_vparam_get_upper_bound (vp, 1), ==, 3.0);

  got = ncm_vparam_get_sparam (vp, 1);
  g_assert_true (got == ncm_vparam_peek_sparam (vp, 1));
  ncm_sparam_free (got);

  cp = ncm_vparam_copy (vp);
  g_assert_cmpuint (ncm_vparam_get_len (cp), ==, 2);
  g_assert_true (ncm_vparam_peek_sparam (cp, 0) != ncm_vparam_peek_sparam (vp, 0));
  g_assert_cmpstr (ncm_sparam_name (ncm_vparam_peek_sparam (cp, 0)), ==, "x");
  g_assert_cmpstr (ncm_sparam_name (ncm_vparam_peek_sparam (cp, 1)), ==, "y");

  ncm_vparam_set_upper_bound (cp, 1, 5.0);
  g_assert_cmpfloat (ncm_vparam_get_upper_bound (vp, 1), ==, 3.0);

  ncm_vparam_free (cp);
  ncm_vparam_free (vp);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/vparam/components", &test_ncm_vparam_components);
  g_test_add_func ("/ncm/vparam/setters", &test_ncm_vparam_setters);
  g_test_add_func ("/ncm/vparam/set_sparam_copy", &test_ncm_vparam_set_sparam_copy);

  g_test_run ();
}

