/***************************************************************************
 *            test_nc_halo_bias.c
 *
 *  Mon Jan 16/02 10:04:16 2023
 *  Copyright  2023  Eduardo José Barroso
 *  <eduardojsbarroso@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Eduardo José Barroso 2023 <eduardojsbarroso@gmail.com>
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

typedef struct _TestNcHaloBias TestNcHaloBias;

struct _TestNcHaloBias
{
  NcHaloBiasPS *ps;
  NcHaloBiasSTEllip *bste;
  NcHaloBiasSTSpher *bsts;
  NcHaloBiasTinker *bt;
  NcHICosmo *cosmo;
  NcHaloMassFunction *mfp;
  gdouble z;
  gdouble lnM;
  gdouble R1, R2, R3;
  gint ntests;
  gboolean (*rng_params) (TestNcHaloBias *test);
};

void test_nc_halo_bias_new (TestNcHaloBias *test, gconstpointer pdata);
void test_nc_halo_bias_eval (TestNcHaloBias *test, gconstpointer pdata);
void test_nc_halo_bias_free (TestNcHaloBias *test, gconstpointer pdata);
void test_nc_halo_bias_set_get (TestNcHaloBias *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/halo_bias/eval", TestNcHaloBias, NULL,
              &test_nc_halo_bias_new,
              &test_nc_halo_bias_eval,
              &test_nc_halo_bias_free);

  g_test_add ("/nc/halo_bias/set_get", TestNcHaloBias, NULL,
              &test_nc_halo_bias_new,
              &test_nc_halo_bias_set_get,
              &test_nc_halo_bias_free);
  g_test_add ("/nc/halo_bias/integrand", TestNcHaloBias, NULL,
              &test_nc_halo_bias_new,
              &test_nc_halo_bias_set_get,
              &test_nc_halo_bias_free);


  g_test_run ();
}

void
test_nc_halo_bias_free (TestNcHaloBias *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_halo_bias_free, NC_HALO_BIAS (test->ps));



  ncm_model_free (NCM_MODEL (test->cosmo));

  nc_halo_mass_function_free (test->mfp);
}

void
test_nc_halo_bias_new (TestNcHaloBias *test, gconstpointer pdata)
{
  NcHICosmo *cosmo         = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist         = nc_distance_new (3.0);
  NcTransferFunc *tf       = NC_TRANSFER_FUNC (ncm_serialize_global_from_string ("NcTransferFuncEH"));
  NcPowspecML *ps_ml       = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  NcmPowspecFilter *psf    = ncm_powspec_filter_new (NCM_POWSPEC (ps_ml), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  NcMultiplicityFunc *mulf = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_tinker_new_full (NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL, 500.0));
  NcHaloMassFunction *mfp  = nc_halo_mass_function_new (dist, psf, mulf);
  NcHaloBiasPS *ps         = nc_halo_bias_ps_new (mfp);
  NcHaloBiasSTEllip *bste  = nc_halo_bias_st_ellip_new (mfp);
  NcHaloBiasSTSpher *bsts  = nc_halo_bias_st_spher_new (mfp);
  NcHaloBiasTinker *bt     = nc_halo_bias_tinker_new (mfp);


  test->cosmo = cosmo;
  test->bsts  = bsts;
  test->bste  = bste;
  test->bt    = bt;
  test->ps    = ps;
  test->mfp   = mfp;
  test->z     = 1.0;
  test->lnM   = 32.23;

  {
    NcHaloBiasPS *ps_test = nc_halo_bias_ps_new_full (mfp, 1.0);

    g_assert_true (NC_IS_HALO_BIAS_PS (ps_test));
    NCM_TEST_FREE (nc_halo_bias_free, NC_HALO_BIAS (ps_test));
  }
  {
    NcHaloBiasTinker *bt_test = nc_halo_bias_tinker_new_full (mfp, 1.0, 2.0, 3.0, 4.0);

    g_assert_true (NC_IS_HALO_BIAS_TINKER (bt_test));
    NCM_TEST_FREE (nc_halo_bias_free, NC_HALO_BIAS (bt_test));
  }
  {
    NcHaloBiasSTEllip *bste_test = nc_halo_bias_st_ellip_new_full (mfp, 1.0, 2.0, 3.0, 4.0);

    g_assert_true (NC_IS_HALO_BIAS_ST_ELLIP (bste_test));
    NCM_TEST_FREE (nc_halo_bias_free, NC_HALO_BIAS (bste_test));
  }
  {
    NcHaloBiasSTSpher *bsts_test = nc_halo_bias_st_spher_new_full (mfp, 1.0, 2.0, 3.0);

    g_assert_true (NC_IS_HALO_BIAS_ST_SPHER (bsts_test));
    NCM_TEST_FREE (nc_halo_bias_free, NC_HALO_BIAS (bsts_test));
  }

  g_assert_true (NC_IS_HALO_BIAS_ST_SPHER (bsts));
  g_assert_true (NC_IS_HALO_BIAS_ST_ELLIP (bste));
  g_assert_true (NC_IS_HALO_BIAS_TINKER (bt));
  g_assert_true (NC_IS_HALO_BIAS_PS (ps));


  nc_multiplicity_func_free (mulf);
  ncm_powspec_filter_free (psf);
  nc_powspec_ml_free (ps_ml);
  nc_transfer_func_free (tf);
  nc_distance_free (dist);
}

void
test_nc_halo_bias_set_get (TestNcHaloBias *test, gconstpointer pdata)
{
  gdouble delta_c = g_test_rand_double_range (1.0, 10.0);
  gdouble a       = g_test_rand_double_range (1.0, 10.0);
  gdouble b       = g_test_rand_double_range (1.0, 10.0);
  gdouble B       = g_test_rand_double_range (1.0, 10.0);
  gdouble c       = g_test_rand_double_range (1.0, 10.0);
  gdouble p       = nc_halo_bias_st_spher_get_p (test->bsts);

  {
    nc_halo_bias_ps_set_delta_c (test->ps, delta_c);
    g_assert_cmpfloat (nc_halo_bias_ps_get_delta_c (test->ps), ==, delta_c);
  }

  {
    nc_halo_bias_tinker_set_delta_c (test->bt, delta_c);
    nc_halo_bias_tinker_set_B (test->bt, B);
    nc_halo_bias_tinker_set_b (test->bt, b);
    nc_halo_bias_tinker_set_c (test->bt, c);

    g_assert_cmpfloat (nc_halo_bias_tinker_get_delta_c (test->bt), ==, delta_c);
    g_assert_cmpfloat (nc_halo_bias_tinker_get_B (test->bt), ==, B);
    g_assert_cmpfloat (nc_halo_bias_tinker_get_b (test->bt), ==, b);
    g_assert_cmpfloat (nc_halo_bias_tinker_get_c (test->bt), ==, c);
  }

  {
    nc_halo_bias_st_spher_set_delta_c (test->bsts, delta_c);
    nc_halo_bias_st_spher_set_a (test->bsts, a);
    nc_halo_bias_st_spher_set_p (test->bsts, p);

    g_assert_cmpfloat (nc_halo_bias_st_spher_get_delta_c (test->bsts), ==, delta_c);
    g_assert_cmpfloat (nc_halo_bias_st_spher_get_a (test->bsts), ==, a);
    g_assert_cmpfloat (nc_halo_bias_st_spher_get_p (test->bsts), ==, p);
  }

  {
    nc_halo_bias_st_ellip_set_delta_c (test->bste, delta_c);
    nc_halo_bias_st_ellip_set_a (test->bste, a);
    nc_halo_bias_st_ellip_set_b (test->bste, b);
    nc_halo_bias_st_ellip_set_c (test->bste, c);

    g_assert_cmpfloat (nc_halo_bias_st_ellip_get_delta_c (test->bste), ==, delta_c);
    g_assert_cmpfloat (nc_halo_bias_st_ellip_get_a (test->bste), ==, a);
    g_assert_cmpfloat (nc_halo_bias_st_ellip_get_b (test->bste), ==, b);
    g_assert_cmpfloat (nc_halo_bias_st_ellip_get_c (test->bste), ==, c);
  }
}

void
test_nc_halo_bias_eval (TestNcHaloBias *test, gconstpointer pdata)
{
  gdouble eval_ps, eval_bste, eval_bsts, eval_bt;

  eval_ps   = nc_halo_bias_eval (NC_HALO_BIAS (test->ps), test->cosmo, 0.0, test->z , test->lnM);
  eval_bste = nc_halo_bias_eval (NC_HALO_BIAS (test->bste), test->cosmo, 5.0, test->z , test->lnM);
  eval_bsts = nc_halo_bias_eval (NC_HALO_BIAS (test->bsts), test->cosmo, GSL_POSINF, test->z , test->lnM);
  eval_bt   = nc_halo_bias_eval (NC_HALO_BIAS (test->bt), test->cosmo, GSL_POSINF, test->z , test->lnM);

  g_assert_cmpfloat (eval_ps, ==, GSL_POSINF);
  g_assert_cmpfloat (eval_bt, ==, 1.0);
  g_assert_true (gsl_finite (eval_bste));

  {
    gdouble delta_c  = nc_halo_bias_st_spher_get_delta_c (test->bsts);
    gdouble p        = nc_halo_bias_st_spher_get_p (test->bsts);
    gdouble bias_inf = 1.0 - 1.0 / delta_c + 2.0 * p / delta_c;

    g_assert_cmpfloat (eval_bsts - bias_inf, >, -1.0e-8);
  }
}

void
test_nc_halo_bias_integrand (TestNcHaloBias *test, gconstpointer pdata)
{
  gdouble lnM            = g_test_rand_double_range (12.0, 15.0);
  gdouble bias_integrand = nc_halo_bias_integrand (NC_HALO_BIAS (test->ps), test->cosmo, lnM, test->z);
  gdouble d2n_dzdlnM     = nc_halo_mass_function_d2n_dzdlnM (test->mfp, test->cosmo, lnM, test->z);
  gdouble sigma          = nc_halo_mass_function_sigma_lnM (test->mfp, test->cosmo, lnM, test->z);
  gdouble bias           = nc_halo_bias_eval (NC_HALO_BIAS (test->ps), test->cosmo, sigma, test->z, test->lnM);


  g_assert_cmpfloat (d2n_dzdlnM * bias, ==, bias_integrand);
}

