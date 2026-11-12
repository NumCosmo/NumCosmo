/***************************************************************************
 *            test_nc_galaxy_shape_factor_spectral_cert.c
 *
 *  Thu Nov 12 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_galaxy_shape_factor_spectral_cert.c
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

/*
 * Acceptance test of NcGalaxyShapeFactorSpectral: ln P at the default reltol against every
 * value of the certification grid of the Arb reference
 * tests/tools/nc_galaxy_shape_marginal_arb.c, downloaded from the data-file release (the
 * table is written by tests/tools/make_shape_marginal_arb_table.py).
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>

#define TEST_SPECTRAL_CERT_FILE "nc_galaxy_shape_marginal_arb_cert.bin"

/*
 * Largest |ln P - ln P_Arb| accepted at the default reltol (1e-5). The largest differences
 * measured on the 26,023 certified points are 6.9e-9 for (|eps_obs| - 1) / std_noise <= 8
 * and |g| != 1, 2.5e-10 for (|eps_obs| - 1) / std_noise > 8 and 2.9e-11 at |g| = 1.
 */
#define TEST_SPECTRAL_CERT_TOL (1.0e-7)

/* The columns of the table */
typedef enum _TestSpectralCertCol
{
  TEST_SPECTRAL_CERT_COL_CONV,
  TEST_SPECTRAL_CERT_COL_POP,
  TEST_SPECTRAL_CERT_COL_P1,
  TEST_SPECTRAL_CERT_COL_P2,
  TEST_SPECTRAL_CERT_COL_STD_NOISE,
  TEST_SPECTRAL_CERT_COL_G1,
  TEST_SPECTRAL_CERT_COL_G2,
  TEST_SPECTRAL_CERT_COL_E1,
  TEST_SPECTRAL_CERT_COL_E2,
  TEST_SPECTRAL_CERT_COL_LNP,
  TEST_SPECTRAL_CERT_COL_LEN,
} TestSpectralCertCol;

void test_nc_galaxy_shape_factor_spectral_cert (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/nc/galaxy_shape_factor_spectral/cert", &test_nc_galaxy_shape_factor_spectral_cert);

  g_test_run ();
}

/*
 * Every row of the certified table at the default reltol: the truncated Gaussian (column
 * pop 0, sigma = p1) and the Beta population (pop 1, alpha = p1, beta = p2) in both
 * conventions.
 */
void
test_nc_galaxy_shape_factor_spectral_cert (void)
{
  gchar *filename           = nc_data_download_release_file (TEST_SPECTRAL_CERT_FILE);
  NcmSerialize *ser         = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmObjDictStr *table      = ncm_serialize_dict_str_from_binfile (ser, filename);
  NcmMatrix *points         = NCM_MATRIX (ncm_obj_dict_str_peek (table, "points"));
  NcGalaxyShapePop *pops[2] = {
    NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ()),
    NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_beta_new ()),
  };
  NcGalaxyShapePopData *pop_data[2] = {
    nc_galaxy_shape_pop_data_new (pops[0]),
    nc_galaxy_shape_pop_data_new (pops[1]),
  };
  NcGalaxyShapeFactorSpectral *gsfs[2] = {
    nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE),
    nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET),
  };
  guint i;

  g_assert_nonnull (points);
  g_assert_cmpuint (ncm_matrix_ncols (points), ==, TEST_SPECTRAL_CERT_COL_LEN);
  g_assert_cmpuint (ncm_matrix_nrows (points), >, 0);

  for (i = 0; i < ncm_matrix_nrows (points); i++)
  {
    const NcGalaxyWLObsEllipConv conv = (NcGalaxyWLObsEllipConv) ncm_matrix_get (points, i, TEST_SPECTRAL_CERT_COL_CONV);
    const guint ip                    = (guint) ncm_matrix_get (points, i, TEST_SPECTRAL_CERT_COL_POP);
    const guint ic                    = (conv == NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE) ? 0 : 1;
    gdouble lnP;

    g_assert_cmpuint (ip, <, 2);
    g_assert_true ((conv == NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE) || (conv == NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET));

    if (ip == 0)
    {
      ncm_model_param_set_by_name (NCM_MODEL (pops[0]), "sigma", ncm_matrix_get (points, i, TEST_SPECTRAL_CERT_COL_P1), NULL);
    }
    else
    {
      ncm_model_param_set_by_name (NCM_MODEL (pops[1]), "alpha", ncm_matrix_get (points, i, TEST_SPECTRAL_CERT_COL_P1), NULL);
      ncm_model_param_set_by_name (NCM_MODEL (pops[1]), "beta", ncm_matrix_get (points, i, TEST_SPECTRAL_CERT_COL_P2), NULL);
    }

    nc_galaxy_shape_pop_prepare (pops[ip], pop_data[ip]);

    lnP = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs[ic], pops[ip], pop_data[ip],
                                                          ncm_matrix_get (points, i, TEST_SPECTRAL_CERT_COL_G1),
                                                          ncm_matrix_get (points, i, TEST_SPECTRAL_CERT_COL_G2),
                                                          ncm_matrix_get (points, i, TEST_SPECTRAL_CERT_COL_E1),
                                                          ncm_matrix_get (points, i, TEST_SPECTRAL_CERT_COL_E2),
                                                          ncm_matrix_get (points, i, TEST_SPECTRAL_CERT_COL_STD_NOISE));

    g_assert_cmpfloat_with_epsilon (lnP, ncm_matrix_get (points, i, TEST_SPECTRAL_CERT_COL_LNP), TEST_SPECTRAL_CERT_TOL);
  }

  for (i = 0; i < 2; i++)
  {
    nc_galaxy_shape_factor_spectral_free (gsfs[i]);
    nc_galaxy_shape_pop_data_unref (pop_data[i]);
    nc_galaxy_shape_pop_free (pops[i]);
  }

  ncm_obj_dict_str_unref (table);
  ncm_serialize_free (ser);
  g_free (filename);
}

