/***************************************************************************
 *            test_nc_hicosmo_de_wspline.c
 *
 *  Wed Oct 21 15:10:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_hicosmo_de_wspline.c
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

#define TEST_NKNOTS 8
#define TEST_ZF 2.33

typedef struct _TestNcHICosmoDEWSpline
{
  NcHICosmoDEWSpline *wspline;
  NcDistance *dist;
} TestNcHICosmoDEWSpline;

/* Omega_k values at which the computed sound horizon is compared, flat and both curvatures. */
static const gdouble _test_Omega_k[] = {0.0, 0.05, -0.05};

void test_nc_hicosmo_de_wspline_new (TestNcHICosmoDEWSpline *test, gconstpointer pdata);
void test_nc_hicosmo_de_wspline_free (TestNcHICosmoDEWSpline *test, gconstpointer pdata);

void test_nc_hicosmo_de_wspline_as_drag_computed (TestNcHICosmoDEWSpline *test, gconstpointer pdata);
void test_nc_hicosmo_de_wspline_as_drag_param (TestNcHICosmoDEWSpline *test, gconstpointer pdata);
void test_nc_hicosmo_de_wspline_as_drag_serialization (TestNcHICosmoDEWSpline *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/hicosmo_de_wspline/as_drag/computed", TestNcHICosmoDEWSpline, NULL,
              &test_nc_hicosmo_de_wspline_new,
              &test_nc_hicosmo_de_wspline_as_drag_computed,
              &test_nc_hicosmo_de_wspline_free);
  g_test_add ("/nc/hicosmo_de_wspline/as_drag/param", TestNcHICosmoDEWSpline, NULL,
              &test_nc_hicosmo_de_wspline_new,
              &test_nc_hicosmo_de_wspline_as_drag_param,
              &test_nc_hicosmo_de_wspline_free);
  g_test_add ("/nc/hicosmo_de_wspline/as_drag/serialization", TestNcHICosmoDEWSpline, NULL,
              &test_nc_hicosmo_de_wspline_new,
              &test_nc_hicosmo_de_wspline_as_drag_serialization,
              &test_nc_hicosmo_de_wspline_free);

  g_test_run ();
}

/* A w = -1 spline, so that the computed sound horizon has an XCDM counterpart. */
void
test_nc_hicosmo_de_wspline_new (TestNcHICosmoDEWSpline *test, gconstpointer pdata)
{
  guint i;

  test->wspline = nc_hicosmo_de_wspline_new (TEST_NKNOTS, TEST_ZF);
  test->dist    = nc_distance_new (3.0);

  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (test->wspline), NULL);

  for (i = 0; i < TEST_NKNOTS; i++)
    ncm_model_orig_vparam_set (NCM_MODEL (test->wspline), NC_HICOSMO_DE_WSPLINE_W, i, -1.0);
}

void
test_nc_hicosmo_de_wspline_free (TestNcHICosmoDEWSpline *test, gconstpointer pdata)
{
  nc_distance_free (test->dist);
  NCM_TEST_FREE (nc_hicosmo_free, NC_HICOSMO (test->wspline));
}

/* By default the sound horizon at drag is computed, equal to NcDistance on the same model as XCDM. */
void
test_nc_hicosmo_de_wspline_as_drag_computed (TestNcHICosmoDEWSpline *test, gconstpointer pdata)
{
  NcHICosmo *cosmo      = NC_HICOSMO (test->wspline);
  NcHICosmo *xcdm       = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  const gchar *shared[] = {"H0", "Omegac", "Omegab", "Tgamma0", "ENnu"};
  guint i, j;

  g_assert_false (nc_hicosmo_de_wspline_get_use_as_drag (test->wspline));
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (xcdm), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (xcdm), "w", -1.0, NULL);

  for (j = 0; j < G_N_ELEMENTS (shared); j++)
    ncm_model_param_set_by_name (NCM_MODEL (xcdm), shared[j], ncm_model_param_get_by_name (NCM_MODEL (cosmo), shared[j], NULL), NULL);

  for (i = 0; i < G_N_ELEMENTS (_test_Omega_k); i++)
  {
    gdouble r_wspline, r_xcdm;

    ncm_model_param_set_by_name (NCM_MODEL (cosmo), "Omegak", _test_Omega_k[i], NULL);
    ncm_model_param_set_by_name (NCM_MODEL (xcdm), "Omegak", _test_Omega_k[i], NULL);

    nc_distance_prepare (test->dist, cosmo);
    r_wspline = nc_distance_r_zd (test->dist, cosmo);
    nc_distance_prepare (test->dist, xcdm);
    r_xcdm = nc_distance_r_zd (test->dist, xcdm);

    ncm_assert_cmpdouble_e (r_wspline, ==, r_xcdm, 1.0e-12, 0.0);
  }

  nc_hicosmo_free (xcdm);
}

/* With use-as-drag set the sound horizon at drag is the asdrag parameter, and unsetting it restores the computation. */
void
test_nc_hicosmo_de_wspline_as_drag_param (TestNcHICosmoDEWSpline *test, gconstpointer pdata)
{
  NcHICosmo *cosmo = NC_HICOSMO (test->wspline);
  gdouble r_computed;

  nc_distance_prepare (test->dist, cosmo);
  r_computed = nc_distance_r_zd (test->dist, cosmo);

  nc_hicosmo_de_wspline_set_use_as_drag (test->wspline, TRUE);
  ncm_model_param_set_by_name (NCM_MODEL (cosmo), "asdrag", 0.0331, NULL);
  nc_distance_prepare_if_needed (test->dist, cosmo);

  g_assert_true (nc_hicosmo_de_wspline_get_use_as_drag (test->wspline));
  ncm_assert_cmpdouble (nc_distance_r_zd (test->dist, cosmo), ==, 0.0331);
  g_assert_cmpfloat (fabs (r_computed - 0.0331), >, 1.0e-4);

  nc_hicosmo_de_wspline_set_use_as_drag (test->wspline, FALSE);
  nc_distance_prepare_if_needed (test->dist, cosmo);
  ncm_assert_cmpdouble (nc_distance_r_zd (test->dist, cosmo), ==, r_computed);
}

/* use-as-drag and asdrag survive serialization. */
void
test_nc_hicosmo_de_wspline_as_drag_serialization (TestNcHICosmoDEWSpline *test, gconstpointer pdata)
{
  NcHICosmoDEWSpline *dup;

  nc_hicosmo_de_wspline_set_use_as_drag (test->wspline, TRUE);
  ncm_model_param_set_by_name (NCM_MODEL (test->wspline), "asdrag", 0.0331, NULL);

  dup = NC_HICOSMO_DE_WSPLINE (ncm_serialize_global_dup_obj (G_OBJECT (test->wspline)));

  g_assert_true (nc_hicosmo_de_wspline_get_use_as_drag (dup));
  ncm_assert_cmpdouble (ncm_model_param_get_by_name (NCM_MODEL (dup), "asdrag", NULL), ==, 0.0331);
  ncm_assert_cmpdouble (nc_hicosmo_as_drag (NC_HICOSMO (dup)), ==, 0.0331);

  NCM_TEST_FREE (nc_hicosmo_free, NC_HICOSMO (dup));
}

