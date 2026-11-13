/***************************************************************************
 *            test_nc_halo_mass_function.c
 *
 *  Fri Nov 13 10:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

typedef struct _TestNcHaloMassFunction
{
  NcHICosmo *cosmo;
  NcmPowspecFilter *psf;
  NcHaloMassFunction *mfp;
} TestNcHaloMassFunction;

typedef enum _TestNcHaloMassFunctionMulf
{
  TEST_NC_HALO_MASS_FUNCTION_MULF_TINKER_200M = 0,
  TEST_NC_HALO_MASS_FUNCTION_MULF_BOCQUET_200C,
  TEST_NC_HALO_MASS_FUNCTION_MULF_PS,
} TestNcHaloMassFunctionMulf;

static const TestNcHaloMassFunctionMulf mulf_tinker_200m  = TEST_NC_HALO_MASS_FUNCTION_MULF_TINKER_200M;
static const TestNcHaloMassFunctionMulf mulf_bocquet_200c = TEST_NC_HALO_MASS_FUNCTION_MULF_BOCQUET_200C;
static const TestNcHaloMassFunctionMulf mulf_ps           = TEST_NC_HALO_MASS_FUNCTION_MULF_PS;

static const gdouble test_lnM[] = {30.0, 32.0, 34.0};
static const gdouble test_z[]   = {0.0, 0.5, 1.0};

void test_nc_halo_mass_function_new (TestNcHaloMassFunction *test, gconstpointer pdata);
void test_nc_halo_mass_function_free (TestNcHaloMassFunction *test, gconstpointer pdata);
void test_nc_halo_mass_function_f_lnR (TestNcHaloMassFunction *test, gconstpointer pdata);
void test_nc_halo_mass_function_f_lnR_ps (TestNcHaloMassFunction *test, gconstpointer pdata);
void test_nc_halo_mass_function_dn_dlnR (TestNcHaloMassFunction *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/halo_mass_function/tinker_200m/f_lnR", TestNcHaloMassFunction, &mulf_tinker_200m,
              &test_nc_halo_mass_function_new,
              &test_nc_halo_mass_function_f_lnR,
              &test_nc_halo_mass_function_free);
  g_test_add ("/nc/halo_mass_function/bocquet_200c/f_lnR", TestNcHaloMassFunction, &mulf_bocquet_200c,
              &test_nc_halo_mass_function_new,
              &test_nc_halo_mass_function_f_lnR,
              &test_nc_halo_mass_function_free);
  g_test_add ("/nc/halo_mass_function/ps/f_lnR", TestNcHaloMassFunction, &mulf_ps,
              &test_nc_halo_mass_function_new,
              &test_nc_halo_mass_function_f_lnR_ps,
              &test_nc_halo_mass_function_free);
  g_test_add ("/nc/halo_mass_function/tinker_200m/dn_dlnR", TestNcHaloMassFunction, &mulf_tinker_200m,
              &test_nc_halo_mass_function_new,
              &test_nc_halo_mass_function_dn_dlnR,
              &test_nc_halo_mass_function_free);
  g_test_add ("/nc/halo_mass_function/bocquet_200c/dn_dlnR", TestNcHaloMassFunction, &mulf_bocquet_200c,
              &test_nc_halo_mass_function_new,
              &test_nc_halo_mass_function_dn_dlnR,
              &test_nc_halo_mass_function_free);

  g_test_run ();
}

void
test_nc_halo_mass_function_new (TestNcHaloMassFunction *test, gconstpointer pdata)
{
  const TestNcHaloMassFunctionMulf mulf_type = *((const TestNcHaloMassFunctionMulf *) pdata);
  NcHIPrimPowerLaw *prim                     = nc_hiprim_power_law_new ();
  NcHICosmo *cosmo                           = NC_HICOSMO (g_object_new (NC_TYPE_HICOSMO_DE_XCDM, "prim", prim, NULL));
  NcDistance *dist                           = nc_distance_new (3.0);
  NcTransferFunc *tf                         = NC_TRANSFER_FUNC (ncm_serialize_global_from_string ("NcTransferFuncEH"));
  NcPowspecML *ps_ml                         = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  NcmPowspecFilter *psf                      = ncm_powspec_filter_new (NCM_POWSPEC (ps_ml), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  NcMultiplicityFunc *mulf                   = NULL;

  switch (mulf_type)
  {
    case TEST_NC_HALO_MASS_FUNCTION_MULF_TINKER_200M:
      mulf = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_tinker_new_full (NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN, 200.0));
      break;
    case TEST_NC_HALO_MASS_FUNCTION_MULF_BOCQUET_200C:
      mulf = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_bocquet_new_full (NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL,
                                                                          NC_MULTIPLICITY_FUNC_BOCQUET_SIM_HYDRO,
                                                                          200.0));
      break;
    case TEST_NC_HALO_MASS_FUNCTION_MULF_PS:
      mulf = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_ps_new ());
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  test->cosmo = cosmo;
  test->psf   = ncm_powspec_filter_ref (psf);
  test->mfp   = nc_halo_mass_function_new (dist, psf, mulf);

  ncm_powspec_filter_prepare_if_needed (psf, NCM_MODEL (cosmo));

  nc_hiprim_free (NC_HIPRIM (prim));
  nc_multiplicity_func_free (mulf);
  ncm_powspec_filter_free (psf);
  nc_powspec_ml_free (ps_ml);
  nc_transfer_func_free (tf);
  nc_distance_free (dist);
}

void
test_nc_halo_mass_function_free (TestNcHaloMassFunction *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_halo_mass_function_free, test->mfp);
  ncm_powspec_filter_free (test->psf);
  ncm_model_free (NCM_MODEL (test->cosmo));
}

void
test_nc_halo_mass_function_f_lnR (TestNcHaloMassFunction *test, gconstpointer pdata)
{
  NcMultiplicityFunc *mulf = nc_halo_mass_function_peek_multiplicity_function (test->mfp);
  const gboolean has_cf    = nc_multiplicity_func_has_correction_factor (mulf);
  guint i, j;

  g_assert_true (has_cf == (*((const TestNcHaloMassFunctionMulf *) pdata) == TEST_NC_HALO_MASS_FUNCTION_MULF_BOCQUET_200C));

  for (i = 0; i < G_N_ELEMENTS (test_lnM); i++)
  {
    for (j = 0; j < G_N_ELEMENTS (test_z); j++)
    {
      const gdouble z      = test_z[j];
      const gdouble lnR    = nc_halo_mass_function_lnM_to_lnR (test->mfp, test->cosmo, test_lnM[i]);
      const gdouble lnM    = nc_halo_mass_function_lnR_to_lnM (test->mfp, test->cosmo, lnR);
      const gdouble sigma  = nc_halo_mass_function_sigma_lnR (test->mfp, test->cosmo, lnR, z);
      const gdouble f_eval = nc_multiplicity_func_eval (mulf, test->cosmo, sigma, lnR, z);
      const gdouble cf     = has_cf ? nc_multiplicity_func_correction_factor (mulf, test->cosmo, sigma, z, lnM) : 1.0;
      const gdouble f      = nc_halo_mass_function_f_lnR (test->mfp, test->cosmo, lnR, z);

      g_assert_cmpfloat (f, >, 0.0);
      ncm_assert_cmpdouble_e (f, ==, f_eval * cf, 1.0e-15, 0.0);

      if (has_cf)
        g_assert_cmpfloat (fabs (cf - 1.0), >, 1.0e-3);
    }
  }
}

void
test_nc_halo_mass_function_f_lnR_ps (TestNcHaloMassFunction *test, gconstpointer pdata)
{
  NcMultiplicityFunc *mulf = nc_halo_mass_function_peek_multiplicity_function (test->mfp);
  const gdouble delta_c    = nc_multiplicity_func_ps_get_delta_c (NC_MULTIPLICITY_FUNC_PS (mulf));
  guint i, j;

  g_assert_false (nc_multiplicity_func_has_correction_factor (mulf));

  for (i = 0; i < G_N_ELEMENTS (test_lnM); i++)
  {
    for (j = 0; j < G_N_ELEMENTS (test_z); j++)
    {
      const gdouble z     = test_z[j];
      const gdouble lnR   = nc_halo_mass_function_lnM_to_lnR (test->mfp, test->cosmo, test_lnM[i]);
      const gdouble sigma = ncm_powspec_filter_eval_sigma_lnr (test->psf, z, lnR);
      const gdouble nu    = delta_c / sigma;
      const gdouble f_PS  = sqrt (2.0 / M_PI) * nu * exp (-0.5 * nu * nu);
      const gdouble f     = nc_halo_mass_function_f_lnR (test->mfp, test->cosmo, lnR, z);

      ncm_assert_cmpdouble_e (f, ==, f_PS, 1.0e-15, 0.0);
    }
  }
}

void
test_nc_halo_mass_function_dn_dlnR (TestNcHaloMassFunction *test, gconstpointer pdata)
{
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (test_lnM); i++)
  {
    for (j = 0; j < G_N_ELEMENTS (test_z); j++)
    {
      const gdouble z           = test_z[j];
      const gdouble lnR         = nc_halo_mass_function_lnM_to_lnR (test->mfp, test->cosmo, test_lnM[i]);
      const gdouble V           = ncm_powspec_filter_volume_rm3 (test->psf) * exp (3.0 * lnR);
      const gdouble dlnvar_dlnR = ncm_powspec_filter_eval_dlnvar_dlnr (test->psf, z, lnR);
      const gdouble f           = nc_halo_mass_function_f_lnR (test->mfp, test->cosmo, lnR, z);
      const gdouble dn_dlnR     = nc_halo_mass_function_dn_dlnR (test->mfp, test->cosmo, lnR, z);

      g_assert_cmpfloat (dn_dlnR, >, 0.0);
      ncm_assert_cmpdouble_e (dn_dlnR, ==, -f * 0.5 * dlnvar_dlnR / V, 1.0e-14, 0.0);
    }
  }
}

