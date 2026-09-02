/***************************************************************************
 *            test_nc_recomb.c
 *
 *  Wed November 14 11:38:56 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

void test_nc_recomb_seager_new (void);
void test_nc_recomb_seager_wmap_zstar (void);
void test_nc_recomb_seager_Xe_ini (void);
void test_nc_recomb_cbe_seager_redshifts (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  /* ncm_cfg_enable_gsl_err_handler (); */

  g_test_add_func ("/nc/recomb/seager/new", &test_nc_recomb_seager_new);
  g_test_add_func ("/nc/recomb/seager/wmap/zstar", &test_nc_recomb_seager_wmap_zstar);
  g_test_add_func ("/nc/recomb/seager/wmap/Xe_ini", &test_nc_recomb_seager_Xe_ini);
  g_test_add_func ("/nc/recomb/cbe/seager/redshifts", &test_nc_recomb_cbe_seager_redshifts);

  g_test_run ();
}

void
test_nc_recomb_seager_new ()
{
  NcRecomb *recomb = NC_RECOMB (nc_recomb_seager_new ());

  g_assert_true (NC_IS_RECOMB (recomb));
  g_assert_true (NC_IS_RECOMB_SEAGER (recomb));

  NCM_TEST_FREE (nc_recomb_free, recomb);

  recomb = NC_RECOMB (nc_recomb_seager_new_full (1.0e-10, 2.2e9, 1.0e-5));
  g_assert_true (NC_IS_RECOMB (recomb));
  g_assert_true (NC_IS_RECOMB_SEAGER (recomb));

  ncm_assert_cmpdouble (recomb->init_frac, ==, 1.0e-10);
  ncm_assert_cmpdouble (recomb->zi, ==, 2.2e9);
  ncm_assert_cmpdouble (recomb->prec, ==, 1.0e-5);

  NCM_TEST_FREE (nc_recomb_free, recomb);

  recomb = NC_RECOMB (nc_recomb_seager_new_full (1.0e-11, NC_RECOMB_STARTING_X, 2.0e-7));

  g_assert_true (NC_IS_RECOMB (recomb));
  g_assert_true (NC_IS_RECOMB_SEAGER (recomb));

  ncm_assert_cmpdouble (recomb->prec, ==, 2.0e-7);

  NCM_TEST_FREE (nc_recomb_free, recomb);
}

void
test_nc_recomb_seager_wmap_zstar (void)
{
  NcRecomb *recomb = NC_RECOMB (nc_recomb_seager_new ());
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  nc_hicosmo_de_set_wmap5_params (NC_HICOSMO_DE (cosmo));
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W,  -1.0);

  nc_recomb_prepare_if_needed (recomb, cosmo);
  {
    const gdouble zstar = nc_recomb_get_tau_z (recomb, cosmo);

    ncm_assert_cmpdouble_e (zstar, ==, 1088.76, 1.0e-4, 0.0);
  }

  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_T_GAMMA0,  2.2250);
  nc_recomb_prepare_if_needed (recomb, cosmo);
  {
    const gdouble zstar = nc_recomb_get_tau_z (recomb, cosmo);

    ncm_assert_cmpdouble_e (zstar, ==, 1325.06, 1.0e-4, 0.0);
  }

  nc_hicosmo_free (cosmo);
  nc_recomb_free (recomb);
}

void
test_nc_recomb_seager_Xe_ini (void)
{
  NcRecomb *recomb = NC_RECOMB (nc_recomb_seager_new ());
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  nc_recomb_prepare_if_needed (recomb, cosmo);

  {
    const gdouble Xe_ini = nc_recomb_Xe (recomb, cosmo, recomb->lambdai);

    ncm_assert_cmpdouble_e (Xe_ini, ==, 1.0 + 2.0 * nc_hicosmo_XHe (cosmo), recomb->prec * 1.0e2, 0.0);
  }

  nc_hicosmo_free (cosmo);
  nc_recomb_free (recomb);
}

void
test_nc_recomb_cbe_seager_redshifts (void)
{
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim   = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new_full (reion, prim, NULL));
  NcRecomb *seager = NC_RECOMB (nc_recomb_seager_new ());
  NcRecomb *cbe    = NC_RECOMB (nc_recomb_cbe_new ());

  nc_recomb_prepare_if_needed (seager, cosmo);
  nc_recomb_prepare_if_needed (cbe, cosmo);

  /*
   * Both backends compute the four characteristic redshifts from their own optical
   * depth splines with the same definitions, so they agree to the accuracy of the two
   * recombination histories. CLASS publishes z_rec and z_d, which the CBE backend does
   * not use: z_rec maximizes the visibility function per unit conformal time, not per
   * unit lambda, and differs from v_tau_max_z by about 4 in z.
   */
  g_assert_true (gsl_finite (nc_recomb_get_v_tau_max_z (cbe, cosmo)));
  g_assert_true (gsl_finite (nc_recomb_get_tau_z (cbe, cosmo)));
  g_assert_true (gsl_finite (nc_recomb_get_tau_drag_z (cbe, cosmo)));
  g_assert_true (gsl_finite (nc_recomb_get_tau_cutoff_z (cbe, cosmo)));

  ncm_assert_cmpdouble_e (nc_recomb_get_v_tau_max_z (cbe, cosmo), ==,
                          nc_recomb_get_v_tau_max_z (seager, cosmo), 1.0e-4, 0.0);
  ncm_assert_cmpdouble_e (nc_recomb_get_tau_z (cbe, cosmo), ==,
                          nc_recomb_get_tau_z (seager, cosmo), 1.0e-4, 0.0);
  ncm_assert_cmpdouble_e (nc_recomb_get_tau_drag_z (cbe, cosmo), ==,
                          nc_recomb_get_tau_drag_z (seager, cosmo), 1.0e-4, 0.0);
  ncm_assert_cmpdouble_e (nc_recomb_get_tau_cutoff_z (cbe, cosmo), ==,
                          nc_recomb_get_tau_cutoff_z (seager, cosmo), 1.0e-3, 0.0);

  /* The lambda fields must match the redshifts they are derived from. */
  ncm_assert_cmpdouble_e (nc_recomb_get_v_tau_max_lambda (cbe, cosmo), ==,
                          -log1p (nc_recomb_get_v_tau_max_z (cbe, cosmo)), 1.0e-14, 0.0);
  ncm_assert_cmpdouble_e (nc_recomb_get_tau_lambda (cbe, cosmo), ==,
                          -log1p (nc_recomb_get_tau_z (cbe, cosmo)), 1.0e-14, 0.0);

  nc_hiprim_free (prim);
  nc_hireion_free (reion);
  nc_hicosmo_free (cosmo);
  nc_recomb_free (seager);
  nc_recomb_free (cbe);
}

