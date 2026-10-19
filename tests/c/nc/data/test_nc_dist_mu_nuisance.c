/***************************************************************************
 *            test_nc_dist_mu_nuisance.c
 *
 *  Sun Oct 18 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_dist_mu_nuisance.c
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

typedef struct _TestNcDistMuNuisance
{
  NcDistMuNuisance *nuisance;
} TestNcDistMuNuisance;

void test_nc_dist_mu_nuisance_new (TestNcDistMuNuisance *test, gconstpointer pdata);
void test_nc_dist_mu_nuisance_free (TestNcDistMuNuisance *test, gconstpointer pdata);

void test_nc_dist_mu_nuisance_params (TestNcDistMuNuisance *test, gconstpointer pdata);
void test_nc_dist_mu_nuisance_sigma_mu_v (TestNcDistMuNuisance *test, gconstpointer pdata);
void test_nc_dist_mu_nuisance_serialize (TestNcDistMuNuisance *test, gconstpointer pdata);
void test_nc_dist_mu_nuisance_sigma_v_per_sample (TestNcDistMuNuisance *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/dist_mu_nuisance/params", TestNcDistMuNuisance, NULL,
              &test_nc_dist_mu_nuisance_new,
              &test_nc_dist_mu_nuisance_params,
              &test_nc_dist_mu_nuisance_free);

  g_test_add ("/nc/dist_mu_nuisance/sigma_mu_v", TestNcDistMuNuisance, NULL,
              &test_nc_dist_mu_nuisance_new,
              &test_nc_dist_mu_nuisance_sigma_mu_v,
              &test_nc_dist_mu_nuisance_free);

  g_test_add ("/nc/dist_mu_nuisance/serialize", TestNcDistMuNuisance, NULL,
              &test_nc_dist_mu_nuisance_new,
              &test_nc_dist_mu_nuisance_serialize,
              &test_nc_dist_mu_nuisance_free);

  g_test_add ("/nc/dist_mu_nuisance/sigma_v_per_sample", TestNcDistMuNuisance, NULL,
              &test_nc_dist_mu_nuisance_new,
              &test_nc_dist_mu_nuisance_sigma_v_per_sample,
              &test_nc_dist_mu_nuisance_free);

  g_test_run ();
}

void
test_nc_dist_mu_nuisance_new (TestNcDistMuNuisance *test, gconstpointer pdata)
{
  test->nuisance = nc_dist_mu_nuisance_new (3);

  g_assert_true (NC_IS_DIST_MU_NUISANCE (test->nuisance));
}

void
test_nc_dist_mu_nuisance_free (TestNcDistMuNuisance *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_dist_mu_nuisance_free, test->nuisance);
}

void
test_nc_dist_mu_nuisance_params (TestNcDistMuNuisance *test, gconstpointer pdata)
{
  NcmModel *model = NCM_MODEL (test->nuisance);
  guint s;

  g_assert_cmpuint (nc_dist_mu_nuisance_get_nsamples (test->nuisance), ==, 3);
  ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_sigma_v (test->nuisance), ==, NC_DIST_MU_NUISANCE_DEFAULT_SIGMA_V);

  for (s = 0; s < 3; s++)
    ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_dmu (test->nuisance, s), ==, NC_DIST_MU_NUISANCE_DEFAULT_DMU);

  ncm_model_orig_param_set (model, NC_DIST_MU_NUISANCE_SIGMA_V, 250.0);
  ncm_model_orig_vparam_set (model, NC_DIST_MU_NUISANCE_DMU, 2, -0.15);

  ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_sigma_v (test->nuisance), ==, 250.0);
  ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_dmu (test->nuisance, 2), ==, -0.15);
  ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_dmu (test->nuisance, 0), ==, 0.0);

  /* A model without offsets */
  {
    NcDistMuNuisance *none = nc_dist_mu_nuisance_new (0);

    g_assert_cmpuint (nc_dist_mu_nuisance_get_nsamples (none), ==, 0);
    nc_dist_mu_nuisance_free (none);
  }
}

void
test_nc_dist_mu_nuisance_sigma_mu_v (TestNcDistMuNuisance *test, gconstpointer pdata)
{
  NcmModel *model = NCM_MODEL (test->nuisance);
  guint j;

  ncm_assert_cmpdouble (nc_dist_mu_nuisance_sigma_mu_v (test->nuisance, 0.01), ==, 0.0);

  ncm_model_orig_param_set (model, NC_DIST_MU_NUISANCE_SIGMA_V, 300.0);

  for (j = 1; j <= 5; j++)
  {
    const gdouble z   = 0.02 * j;
    const gdouble ref = 5.0 / M_LN10 * (300.0 / (1.0e-3 * ncm_c_c ())) * (1.0 + z) / (z * (1.0 + 0.5 * z));

    ncm_assert_cmpdouble_e (nc_dist_mu_nuisance_sigma_mu_v (test->nuisance, z), ==, ref, 1.0e-15, 0.0);
  }

  /* 300 km/s at z = 0.01 is about 0.21 mag */
  g_assert_cmpfloat (nc_dist_mu_nuisance_sigma_mu_v (test->nuisance, 0.01), >, 0.2);
  g_assert_cmpfloat (nc_dist_mu_nuisance_sigma_mu_v (test->nuisance, 0.01), <, 0.23);
}

void
test_nc_dist_mu_nuisance_serialize (TestNcDistMuNuisance *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmModel *model   = NCM_MODEL (test->nuisance);
  NcDistMuNuisance *dup;

  ncm_model_orig_param_set (model, NC_DIST_MU_NUISANCE_SIGMA_V, 120.0);
  ncm_model_orig_vparam_set (model, NC_DIST_MU_NUISANCE_DMU, 1, 0.05);

  dup = NC_DIST_MU_NUISANCE (ncm_serialize_dup_obj (ser, G_OBJECT (test->nuisance)));

  g_assert_cmpuint (nc_dist_mu_nuisance_get_nsamples (dup), ==, 3);
  ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_sigma_v (dup), ==, 120.0);
  ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_dmu (dup, 1), ==, 0.05);

  nc_dist_mu_nuisance_free (dup);
  ncm_serialize_free (ser);
}

void
test_nc_dist_mu_nuisance_sigma_v_per_sample (TestNcDistMuNuisance *test, gconstpointer pdata)
{
  NcDistMuNuisance *per = nc_dist_mu_nuisance_new_full (3, TRUE);
  NcmModel *model       = NCM_MODEL (per);
  const gdouble sv[3]   = {150.0, 0.0, 420.0};
  const gdouble z       = 0.02;
  gint s;

  /* The plain constructor has no per-sample dispersion: every sample uses sigma_v */
  g_assert_cmpuint (nc_dist_mu_nuisance_get_nsigma_v (test->nuisance), ==, 0);
  ncm_model_orig_param_set (NCM_MODEL (test->nuisance), NC_DIST_MU_NUISANCE_SIGMA_V, 250.0);

  for (s = -1; s < 4; s++)
  {
    ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_sigma_v_sample (test->nuisance, s), ==, 250.0);
    ncm_assert_cmpdouble (nc_dist_mu_nuisance_sigma_mu_v_sample (test->nuisance, z, s), ==,
                          nc_dist_mu_nuisance_sigma_mu_v (test->nuisance, z));
  }

  /* With one per sample: the sample's value inside the range, sigma_v outside it */
  g_assert_cmpuint (nc_dist_mu_nuisance_get_nsamples (per), ==, 3);
  g_assert_cmpuint (nc_dist_mu_nuisance_get_nsigma_v (per), ==, 3);
  g_assert_cmpint (ncm_model_param_get_ftype (model, ncm_model_vparam_index (model, NC_DIST_MU_NUISANCE_SIGMA_V_S, 0)), ==, NCM_PARAM_TYPE_FIXED);
  ncm_model_orig_param_set (model, NC_DIST_MU_NUISANCE_SIGMA_V, 250.0);

  for (s = 0; s < 3; s++)
    ncm_model_orig_vparam_set (model, NC_DIST_MU_NUISANCE_SIGMA_V_S, s, sv[s]);

  for (s = 0; s < 3; s++)
  {
    const gdouble ref = 5.0 / M_LN10 * (sv[s] / (1.0e-3 * ncm_c_c ())) * (1.0 + z) / (z * (1.0 + 0.5 * z));

    ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_sigma_v_sample (per, s), ==, sv[s]);
    ncm_assert_cmpdouble_e (nc_dist_mu_nuisance_sigma_mu_v_sample (per, z, s), ==, ref, 1.0e-15, 0.0);
  }

  ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_sigma_v_sample (per, -1), ==, 250.0);
  ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_sigma_v_sample (per, 3), ==, 250.0);
  ncm_assert_cmpdouble (nc_dist_mu_nuisance_sigma_mu_v_sample (per, z, 7), ==, nc_dist_mu_nuisance_sigma_mu_v (per, z));

  /* Serialization keeps the per-sample values */
  {
    NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    NcDistMuNuisance *dup = NC_DIST_MU_NUISANCE (ncm_serialize_dup_obj (ser, G_OBJECT (per)));

    g_assert_cmpuint (nc_dist_mu_nuisance_get_nsigma_v (dup), ==, 3);
    ncm_assert_cmpdouble (nc_dist_mu_nuisance_get_sigma_v_sample (dup, 2), ==, 420.0);

    nc_dist_mu_nuisance_free (dup);
    ncm_serialize_free (ser);
  }

  nc_dist_mu_nuisance_free (per);
}

