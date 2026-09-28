/***************************************************************************
 *            test_ncm_powspec_spline2d.c
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

/*
 * The table is ln P = f(z, y), y = ln k, cubic in z and in y, which the bicubic
 * not-a-knot spline reproduces up to rounding, on z in [0.5, 3] and k in [1e-4, 10].
 */
#define TEST_ZMIN (0.5)
#define TEST_ZMAX (3.0)
#define TEST_NZ 26
#define TEST_KMIN (1.0e-4)
#define TEST_KMAX (10.0)
#define TEST_NK 60

static gdouble
_f (const gdouble z, const gdouble y)
{
  return log (2.0e4) + 0.96 * y - 0.02 * y * y - 0.001 * y * y * y - 0.8 * z + 0.1 * z * z - 0.01 * z * z * z + 0.05 * z * y;
}

static gdouble
_f_y (const gdouble z, const gdouble y)
{
  return 0.96 - 0.04 * y - 0.003 * y * y + 0.05 * z;
}

static gdouble
_f_z (const gdouble z, const gdouble y)
{
  return -0.8 + 0.2 * z - 0.03 * z * z + 0.05 * y;
}

static NcmSpline2d *
_test_table_new (void)
{
  NcmVector *zv  = ncm_vector_new (TEST_NZ);
  NcmVector *yv  = ncm_vector_new (TEST_NK);
  NcmMatrix *lnP = ncm_matrix_new (TEST_NK, TEST_NZ);
  NcmSpline2d *s2d;
  guint i, j;

  for (j = 0; j < TEST_NZ; j++)
    ncm_vector_set (zv, j, TEST_ZMIN + (TEST_ZMAX - TEST_ZMIN) * j / (TEST_NZ - 1.0));

  for (i = 0; i < TEST_NK; i++)
    ncm_vector_set (yv, i, log (TEST_KMIN) + log (TEST_KMAX / TEST_KMIN) * i / (TEST_NK - 1.0));

  for (i = 0; i < TEST_NK; i++)
    for (j = 0; j < TEST_NZ; j++)
      ncm_matrix_set (lnP, i, j, _f (ncm_vector_get (zv, j), ncm_vector_get (yv, i)));

  s2d = NCM_SPLINE2D (ncm_spline2d_bicubic_notaknot_new ());
  ncm_spline2d_set (s2d, zv, yv, lnP, TRUE);

  ncm_vector_free (zv);
  ncm_vector_free (yv);
  ncm_matrix_free (lnP);

  return s2d;
}

static NcmPowspec *
_test_powspec_new (void)
{
  NcmSpline2d *s2d = _test_table_new ();
  NcmPowspec *ps   = NCM_POWSPEC (ncm_powspec_spline2d_new (s2d));

  ncm_spline2d_free (s2d);
  ncm_powspec_prepare (ps, NULL);

  return ps;
}

/* The table's knots set the ranges. */
static void
test_ncm_powspec_spline2d_ranges (void)
{
  NcmPowspec *ps = _test_powspec_new ();
  guint Nz, Nk;

  g_assert_cmpfloat (ncm_powspec_get_zi (ps), ==, TEST_ZMIN);
  g_assert_cmpfloat (ncm_powspec_get_zf (ps), ==, TEST_ZMAX);
  {
    NcmVector *yv = ncm_spline2d_peek_yv (ncm_powspec_spline2d_peek_spline2d (NCM_POWSPEC_SPLINE2D (ps)));

    g_assert_cmpfloat (ncm_powspec_get_kmin (ps), ==, exp (ncm_vector_get (yv, 0)));
    g_assert_cmpfloat (ncm_powspec_get_kmax (ps), ==, exp (ncm_vector_get (yv, TEST_NK - 1)));
  }

  ncm_powspec_get_nknots (ps, &Nz, &Nk);
  g_assert_cmpuint (Nz, ==, TEST_NZ);
  g_assert_cmpuint (Nk, ==, TEST_NK);

  ncm_powspec_free (ps);
}

/* Inside the table P, dP/dz and dP/dk reproduce f; measured 5e-15, 2.3e-13 and 2.8e-14. */
static void
test_ncm_powspec_spline2d_inside (void)
{
  NcmPowspec *ps = _test_powspec_new ();
  guint i, j;

  for (i = 0; i <= 16; i++)
  {
    const gdouble z = TEST_ZMIN + (TEST_ZMAX - TEST_ZMIN) * i / 16.0;

    for (j = 1; j < 100; j++)
    {
      const gdouble y = log (TEST_KMIN) + log (TEST_KMAX / TEST_KMIN) * j / 100.0;
      const gdouble k = exp (y);
      const gdouble P = exp (_f (z, y));

      ncm_assert_cmpdouble_e (ncm_powspec_eval (ps, NULL, z, k), ==, P, 2.0e-14, 0.0);
      ncm_assert_cmpdouble_e (ncm_powspec_deriv_z (ps, NULL, z, k), ==, P * _f_z (z, y), 1.0e-12, 0.0);
      ncm_assert_cmpdouble_e (ncm_powspec_deriv_k (ps, NULL, z, k), ==, P * _f_y (z, y) / k, 1.0e-13, 0.0);
    }
  }

  ncm_powspec_free (ps);
}

/*
 * Outside [kmin, kmax], ln P = f(z, e) + f_y(z, e) (y - e) with e the nearest end, over
 * four decades on each side; measured 2.1e-13 on P and dP/dk and 8.9e-11 on dP/dz.
 */
static void
test_ncm_powspec_spline2d_continuation (void)
{
  NcmPowspec *ps = _test_powspec_new ();
  guint i, j;

  for (i = 0; i <= 16; i++)
  {
    const gdouble z = TEST_ZMIN + (TEST_ZMAX - TEST_ZMIN) * i / 16.0;

    for (j = 0; j < 16; j++)
    {
      const gdouble k = (j < 8) ? TEST_KMIN * pow (1.0e-4, (8.0 - j) / 8.0) : TEST_KMAX *pow (1.0e3, (j - 7.0) / 8.0);

      const gdouble y    = log (k);
      const gdouble e    = log ((j < 8) ? TEST_KMIN : TEST_KMAX);
      const gdouble P    = exp (_f (z, e) + _f_y (z, e) * (y - e));
      const gdouble dlnz = _f_z (z, e) + 0.05 * (y - e);

      ncm_assert_cmpdouble_e (ncm_powspec_eval (ps, NULL, z, k), ==, P, 1.0e-12, 0.0);
      ncm_assert_cmpdouble_e (ncm_powspec_deriv_z (ps, NULL, z, k), ==, P * dlnz, 3.0e-10, 0.0);
      ncm_assert_cmpdouble_e (ncm_powspec_deriv_k (ps, NULL, z, k), ==, P * _f_y (z, e) / k, 1.0e-12, 0.0);
    }
  }

  /* The slope is continuous across both ends; measured jump 2.8e-13. */
  ncm_assert_cmpdouble_e (ncm_powspec_deriv_k (ps, NULL, 1.2, TEST_KMIN * (1.0 - 1.0e-12)), ==,
                          ncm_powspec_deriv_k (ps, NULL, 1.2, TEST_KMIN * (1.0 + 1.0e-12)), 1.0e-12, 0.0);
  ncm_assert_cmpdouble_e (ncm_powspec_deriv_k (ps, NULL, 1.2, TEST_KMAX * (1.0 - 1.0e-12)), ==,
                          ncm_powspec_deriv_k (ps, NULL, 1.2, TEST_KMAX * (1.0 + 1.0e-12)), 1.0e-12, 0.0);

  ncm_powspec_free (ps);
}

/* eval_vec is eval at each k; a serialized copy evaluates identically. */
static void
test_ncm_powspec_spline2d_eval_vec_dup (void)
{
  NcmPowspec *ps    = _test_powspec_new ();
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmPowspec *dup   = NCM_POWSPEC (ncm_serialize_dup_obj (ser, G_OBJECT (ps)));
  NcmVector *k      = ncm_vector_new (30);
  NcmVector *Pk     = ncm_vector_new (30);
  guint i;

  ncm_powspec_prepare (dup, NULL);

  for (i = 0; i < 30; i++)
    ncm_vector_set (k, i, 1.0e-6 * pow (1.0e10, i / 29.0));

  ncm_powspec_eval_vec (ps, NULL, 1.7, k, Pk);

  for (i = 0; i < 30; i++)
  {
    const gdouble ki = ncm_vector_get (k, i);

    g_assert_cmpfloat (ncm_vector_get (Pk, i), ==, ncm_powspec_eval (ps, NULL, 1.7, ki));
    g_assert_cmpfloat (ncm_powspec_eval (dup, NULL, 1.7, ki), ==, ncm_powspec_eval (ps, NULL, 1.7, ki));
  }

  ncm_vector_free (k);
  ncm_vector_free (Pk);
  ncm_powspec_free (dup);
  ncm_serialize_free (ser);
  ncm_powspec_free (ps);
}

/* get_spline_2d is the base default, a spline of P on (z, k); measured 1.7e-9 at the default reltol. */
static void
test_ncm_powspec_spline2d_get_spline_2d (void)
{
  NcmPowspec *ps   = _test_powspec_new ();
  NcmSpline2d *s2d = ncm_powspec_get_spline_2d (ps, NULL);
  guint i, j;

  for (i = 0; i <= 16; i++)
  {
    const gdouble z = TEST_ZMIN + (TEST_ZMAX - TEST_ZMIN) * i / 16.0;

    for (j = 1; j < 100; j++)
    {
      const gdouble k = TEST_KMIN * pow (TEST_KMAX / TEST_KMIN, j / 100.0);

      ncm_assert_cmpdouble_e (ncm_spline2d_eval (s2d, z, k), ==, ncm_powspec_eval (ps, NULL, z, k), ncm_powspec_get_reltol_spline (ps), 0.0);
    }
  }

  ncm_spline2d_free (s2d);
  ncm_powspec_free (ps);
}

/* Setting the table the object already holds, with no other reference to it. */
static void
test_ncm_powspec_spline2d_set_same (void)
{
  NcmPowspec *ps             = _test_powspec_new ();
  NcmPowspecSpline2d *ps_s2d = NCM_POWSPEC_SPLINE2D (ps);

  ncm_powspec_spline2d_set_spline2d (ps_s2d, ncm_powspec_spline2d_peek_spline2d (ps_s2d));
  ncm_powspec_prepare (ps, NULL);

  ncm_assert_cmpdouble_e (ncm_powspec_eval (ps, NULL, 1.0, 0.1), ==, exp (_f (1.0, log (0.1))), 2.0e-14, 0.0);

  ncm_powspec_free (ps);
}

static void
test_ncm_powspec_spline2d_eval_outside_z (void)
{
  g_test_trap_subprocess ("/ncm/powspec_spline2d/eval_outside_z/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*is outside the table*");
}

static void
test_ncm_powspec_spline2d_eval_outside_z_subprocess (void)
{
  NcmPowspec *ps = _test_powspec_new ();

  ncm_powspec_eval (ps, NULL, TEST_ZMAX + 1.0e-10, 0.1);
}

static void
test_ncm_powspec_spline2d_prepare_outside_z (void)
{
  g_test_trap_subprocess ("/ncm/powspec_spline2d/prepare_outside_z/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*is not inside the table*");
}

static void
test_ncm_powspec_spline2d_prepare_outside_z_subprocess (void)
{
  NcmPowspec *ps = _test_powspec_new ();

  ncm_powspec_require_zi (ps, 0.1);
  ncm_powspec_prepare (ps, NULL);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/powspec_spline2d/ranges", &test_ncm_powspec_spline2d_ranges);
  g_test_add_func ("/ncm/powspec_spline2d/inside", &test_ncm_powspec_spline2d_inside);
  g_test_add_func ("/ncm/powspec_spline2d/continuation", &test_ncm_powspec_spline2d_continuation);
  g_test_add_func ("/ncm/powspec_spline2d/eval_vec_dup", &test_ncm_powspec_spline2d_eval_vec_dup);
  g_test_add_func ("/ncm/powspec_spline2d/get_spline_2d", &test_ncm_powspec_spline2d_get_spline_2d);
  g_test_add_func ("/ncm/powspec_spline2d/set_same", &test_ncm_powspec_spline2d_set_same);
  g_test_add_func ("/ncm/powspec_spline2d/eval_outside_z", &test_ncm_powspec_spline2d_eval_outside_z);
  g_test_add_func ("/ncm/powspec_spline2d/eval_outside_z/subprocess", &test_ncm_powspec_spline2d_eval_outside_z_subprocess);
  g_test_add_func ("/ncm/powspec_spline2d/prepare_outside_z", &test_ncm_powspec_spline2d_prepare_outside_z);
  g_test_add_func ("/ncm/powspec_spline2d/prepare_outside_z/subprocess", &test_ncm_powspec_spline2d_prepare_outside_z_subprocess);

  g_test_run ();
}

