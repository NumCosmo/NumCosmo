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
void test_nc_hicosmo_de_wspline_lp_ls (TestNcHICosmoDEWSpline *test, gconstpointer pdata);
void test_nc_hicosmo_de_wspline_lp_ls_wrong_dim (void);
void test_nc_hicosmo_de_wspline_lp_ls_wrong_dim_subprocess (void);
void test_nc_hicosmo_de_wspline_lp_ls_bad_p_subprocess (void);

static void _test_nc_hicosmo_de_wspline_ls_check (NcHICosmoDEWSpline *model, NcmMSet *mset, NcmSpline *weight, const guint p, const gdouble sigma);
static void _test_nc_hicosmo_de_wspline_curved (NcHICosmoDEWSpline *wspline);

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
  g_test_add ("/nc/hicosmo_de_wspline/lp_ls", TestNcHICosmoDEWSpline, NULL,
              &test_nc_hicosmo_de_wspline_new,
              &test_nc_hicosmo_de_wspline_lp_ls,
              &test_nc_hicosmo_de_wspline_free);
  g_test_add_func ("/nc/hicosmo_de_wspline/lp_ls/wrong_dim", &test_nc_hicosmo_de_wspline_lp_ls_wrong_dim);
  g_test_add_func ("/nc/hicosmo_de_wspline/lp_ls/wrong_dim/subprocess", &test_nc_hicosmo_de_wspline_lp_ls_wrong_dim_subprocess);
  g_test_add_func ("/nc/hicosmo_de_wspline/lp_ls/bad_p/subprocess", &test_nc_hicosmo_de_wspline_lp_ls_bad_p_subprocess);

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

/* By default NcDistance computes the sound horizon at drag, as for the same model as XCDM. */
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

static void
_test_nc_hicosmo_de_wspline_ls_check (NcHICosmoDEWSpline *model, NcmMSet *mset, NcmSpline *weight, const guint p, const gdouble sigma)
{
  const gdouble x      = p;
  const guint len      = nc_hicosmo_de_wspline_d2_lp_residuals_len (model, weight, p);
  NcmMSetFuncList *fls = ncm_mset_func_list_new_dim ((weight == NULL) ? "NcHICosmoDEWSpline:lp_w2_ls" : "NcHICosmoDEWSpline:wlp_w2_ls", G_OBJECT (weight), len);
  NcmMSetFuncList *fsc = ncm_mset_func_list_new ((weight == NULL) ? "NcHICosmoDEWSpline:lp_w2" : "NcHICosmoDEWSpline:wlp_w2", G_OBJECT (weight));
  NcmPrior *pls        = NCM_PRIOR (ncm_prior_gauss_func_new (NCM_MSET_FUNC (fls), 0.0, sigma, x));
  NcmPrior *psc        = NCM_PRIOR (ncm_prior_gauss_func_new (NCM_MSET_FUNC (fsc), 0.0, sigma, x));
  GArray *f            = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), len);
  const gdouble norm   = (weight == NULL) ? nc_hicosmo_de_wspline_lp_norm (model, NCM_SPLINE_CURVATURE_D2, p) : nc_hicosmo_de_wspline_weighted_lp_norm (model, NCM_SPLINE_CURVATURE_D2, p, weight);
  gdouble x_f          = x;
  gdouble f2           = 0.0;
  gdouble prior_f2     = 0.0;
  guint k;

  g_assert_cmpuint (ncm_mset_func_get_dim (NCM_MSET_FUNC (pls)), ==, len);
  g_array_set_size (f, len);

  ncm_mset_func_eval (NCM_MSET_FUNC (fls), mset, &x_f, &g_array_index (f, gdouble, 0));

  for (k = 0; k < len; k++)
    f2 += gsl_pow_2 (g_array_index (f, gdouble, k));

  /* The residuals reproduce N_p^2: to the 1e-10 of the unweighted norm, and to the 2e-10 of
   * the weighted norm, whose numerator and denominator are integrals to 1e-10 each. */
  ncm_assert_cmpdouble_e (f2, ==, norm * norm, (weight == NULL) ? 1.0e-10 : 2.0e-10, 0.0);

  /* The prior on the residuals is the Gaussian prior on N_p, (N_p / sigma)^2. */
  ncm_mset_func_eval (NCM_MSET_FUNC (pls), mset, NULL, &g_array_index (f, gdouble, 0));

  for (k = 0; k < len; k++)
    prior_f2 += gsl_pow_2 (g_array_index (f, gdouble, k));

  ncm_assert_cmpdouble_e (prior_f2, ==, gsl_pow_2 (ncm_mset_func_eval0 (NCM_MSET_FUNC (psc), mset)),
                          (weight == NULL) ? 1.0e-10 : 2.0e-10, 0.0);

  g_array_unref (f);
  ncm_prior_free (pls);
  ncm_prior_free (psc);
  ncm_mset_func_list_free (fls);
  ncm_mset_func_list_free (fsc);
}

/* A curved w(alpha): w_i = -1 + 0.3 sin (1.7 i) + 0.05 i. */
static void
_test_nc_hicosmo_de_wspline_curved (NcHICosmoDEWSpline *wspline)
{
  guint i;

  for (i = 0; i < TEST_NKNOTS; i++)
    ncm_model_orig_vparam_set (NCM_MODEL (wspline), NC_HICOSMO_DE_WSPLINE_W, i, -1.0 + 0.3 * sin (1.7 * i) + 0.05 * i);
}

/*
 * The residual functions lp_w2_ls and wlp_w2_ls reproduce the squares of lp_w2 and wlp_w2 for
 * p = 2 and 10, and a Gaussian prior on them is the Gaussian prior on the norm.
 */
void
test_nc_hicosmo_de_wspline_lp_ls (TestNcHICosmoDEWSpline *test, gconstpointer pdata)
{
  const guint ps[]  = {2, 10};
  NcmMSet *mset     = ncm_mset_new (test->wspline, NULL, NULL);
  const gdouble a_f = log1p (TEST_ZF);
  NcmVector *xv     = ncm_vector_new (9);
  NcmVector *yv     = ncm_vector_new (9);
  NcmSpline *weight;
  guint i;

  _test_nc_hicosmo_de_wspline_curved (test->wspline);

  for (i = 0; i < 9; i++)
  {
    const gdouble a = a_f * i / 8.0;

    ncm_vector_set (xv, i, a);
    ncm_vector_set (yv, i, 1.0 + 0.5 * a + 0.2 * a * a);
  }

  weight = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (xv, yv, TRUE));

  for (i = 0; i < G_N_ELEMENTS (ps); i++)
  {
    _test_nc_hicosmo_de_wspline_ls_check (test->wspline, mset, NULL, ps[i], 0.7);
    _test_nc_hicosmo_de_wspline_ls_check (test->wspline, mset, weight, ps[i], 0.7);
  }

  ncm_spline_free (weight);
  ncm_vector_free (xv);
  ncm_vector_free (yv);
  ncm_mset_free (mset);
}

void
test_nc_hicosmo_de_wspline_lp_ls_wrong_dim (void)
{
  g_test_trap_subprocess ("/nc/hicosmo_de_wspline/lp_ls/wrong_dim/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*NcHICosmoDEWSpline:lp_w2_ls: the function has dimension 15, the model has 14 residuals at p = 2*");

  g_test_trap_subprocess ("/nc/hicosmo_de_wspline/lp_ls/bad_p/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*NcHICosmoDEWSpline:lp_w2_ls: the order p = 2.5 must be an even positive integer*");
}

void
test_nc_hicosmo_de_wspline_lp_ls_wrong_dim_subprocess (void)
{
  NcHICosmoDEWSpline *wspline = nc_hicosmo_de_wspline_new (TEST_NKNOTS, TEST_ZF);
  NcmMSet *mset               = ncm_mset_new (wspline, NULL, NULL);
  NcmMSetFuncList *f          = ncm_mset_func_list_new_dim ("NcHICosmoDEWSpline:lp_w2_ls", NULL, 15);
  gdouble x                   = 2.0;
  gdouble res[15];

  ncm_mset_func_eval (NCM_MSET_FUNC (f), mset, &x, res);
}

void
test_nc_hicosmo_de_wspline_lp_ls_bad_p_subprocess (void)
{
  NcHICosmoDEWSpline *wspline = nc_hicosmo_de_wspline_new (TEST_NKNOTS, TEST_ZF);
  NcmMSet *mset               = ncm_mset_new (wspline, NULL, NULL);
  NcmMSetFuncList *f          = ncm_mset_func_list_new_dim ("NcHICosmoDEWSpline:lp_w2_ls", NULL, 14);
  gdouble x                   = 2.5;
  gdouble res[14];

  ncm_mset_func_eval (NCM_MSET_FUNC (f), mset, &x, res);
}

