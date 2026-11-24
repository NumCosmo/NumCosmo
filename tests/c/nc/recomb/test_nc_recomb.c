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
void test_nc_recomb_v_tau_reion_min (void);
void test_nc_recomb_seager_Xe_panels (void);
void test_nc_recomb_seager_Xe_panels_tight (void);

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
  g_test_add_func ("/nc/recomb/v_tau_reion_min", &test_nc_recomb_v_tau_reion_min);
  g_test_add_func ("/nc/recomb/seager/Xe_panels", &test_nc_recomb_seager_Xe_panels);
  g_test_add_func ("/nc/recomb/seager/Xe_panels/tight", &test_nc_recomb_seager_Xe_panels_tight);

  g_test_run ();
}

static void
_test_nc_recomb_seager_Xe_panels_check (NcRecomb *recomb, NcHICosmo *cosmo, gboolean cmp_spline)
{
  NcmSpectralPanels *sp   = nc_recomb_seager_peek_Xe_recomb_panels (NC_RECOMB_SEAGER (recomb));
  const gdouble lambda_sw = -log (nc_recomb_HeII_ion_saha_x_by_HeIII_He (cosmo, recomb->init_frac));
  const guint n_panels    = ncm_spectral_panels_get_n_panels (sp);
  gboolean edge_at_switch = FALSE;
  gdouble a, b;
  guint i;

  g_assert_cmpuint (n_panels, >, 0);
  ncm_spectral_panels_get_range (sp, &a, &b);
  ncm_assert_cmpdouble_e (a, ==, recomb->lambdai, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (b, ==, 0.0, 0.0, 1.0e-15);

  /* the Saha-to-ODE switch is a breakpoint: a panel starts there */
  for (i = 0; i < n_panels; i++)
  {
    gdouble pa, pb;
    guint pN;

    ncm_spectral_panels_peek_panel (sp, i, &pa, &pb, &pN);

    if (fabs (pa - lambda_sw) <= 1.0e-14 * fabs (lambda_sw))
      edge_at_switch = TRUE;
  }

  g_assert_true (edge_at_switch);

  /* Saha phase: the expanded function (the Saha equilibrium), to the relative tolerance prec */
  for (i = 0; i < 200; i++)
  {
    const gdouble lambda = recomb->lambdai + (lambda_sw - recomb->lambdai) * (i + 0.5) / 200.0;
    const gdouble Xe     = ncm_spectral_panels_eval_row (sp, 0, lambda);

    ncm_assert_cmpdouble_e (Xe, ==, nc_recomb_equilibrium_Xe (recomb, cosmo, exp (-lambda)), recomb->prec, 0.0);
  }

  /* ODE phase: the CVODE solution, which the X_HII and X_HeII splines interpolate to about prec */
  if (cmp_spline)
  {
    for (i = 0; i < 200; i++)
    {
      const gdouble lambda = lambda_sw * (i + 0.5) / 200.0;
      const gdouble Xe     = ncm_spectral_panels_eval_row (sp, 0, lambda);

      ncm_assert_cmpdouble_e (Xe, ==, nc_recomb_XHII (recomb, cosmo, lambda) + nc_recomb_XHeII (recomb, cosmo, lambda), recomb->prec, 0.0);
    }
  }
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

    /* the Saha equilibrium, which is not 1 + 2 X_He at z_i = 1e12 (see nc_recomb_equilibrium_Xe) */
    ncm_assert_cmpdouble_e (Xe_ini, ==, nc_recomb_equilibrium_Xe (recomb, cosmo, exp (-recomb->lambdai)), recomb->prec, 0.0);
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

void
test_nc_recomb_v_tau_reion_min (void)
{
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim   = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new_full (reion, prim, NULL));
  NcRecomb *seager = NC_RECOMB (nc_recomb_seager_new ());
  NcRecomb *cbe    = NC_RECOMB (nc_recomb_cbe_new ());
  NcRecomb *recombs[2];
  guint i;

  recombs[0] = seager;
  recombs[1] = cbe;

  for (i = 0; i < 2; i++)
  {
    NcRecomb *recomb          = recombs[i];
    const gdouble lambda_min  = nc_recomb_get_v_tau_reion_min_lambda (recomb, cosmo);
    const gdouble z_min       = nc_recomb_get_v_tau_reion_min_z (recomb, cosmo);
    const gdouble lambda_max  = nc_recomb_get_v_tau_max_lambda (recomb, cosmo);
    const gdouble lambda_init = -log (nc_hireion_get_init_x (reion, cosmo));
    const gdouble log_v_min   = nc_recomb_log_v_tau (recomb, cosmo, lambda_min);
    const gdouble dlambda     = 1.0e-2;

    /*
     * The minimum lies in the valley between the recombination shell and the
     * reionization bump: below the reionization onset, far above the bump, with
     * the visibility much smaller than at the peak and rising on both sides.
     */
    g_assert_true (gsl_finite (lambda_min));
    g_assert_cmpfloat (lambda_min, >, lambda_max);
    g_assert_cmpfloat (lambda_min, >, lambda_init - 1.0);
    g_assert_cmpfloat (z_min, >, 5.0);
    g_assert_cmpfloat (z_min, <, 50.0);
    ncm_assert_cmpdouble_e (lambda_min, ==, -log1p (z_min), 1.0e-14, 0.0);

    g_assert_cmpfloat (log_v_min, <, nc_recomb_log_v_tau (recomb, cosmo, lambda_max) - 4.0 * M_LN10);
    g_assert_cmpfloat (log_v_min, <, nc_recomb_log_v_tau (recomb, cosmo, lambda_min - dlambda));
    g_assert_cmpfloat (log_v_min, <, nc_recomb_log_v_tau (recomb, cosmo, lambda_min + dlambda));

    /* The cached value is returned as is. */
    g_assert_cmpfloat (lambda_min, ==, nc_recomb_get_v_tau_reion_min_lambda (recomb, cosmo));
  }

  /* Both recombination histories share the ionization floor and the reionization model. */
  ncm_assert_cmpdouble_e (nc_recomb_get_v_tau_reion_min_z (cbe, cosmo), ==,
                          nc_recomb_get_v_tau_reion_min_z (seager, cosmo), 5.0e-2, 0.0);

  /* Without reionization the visibility keeps falling and the minimum is today. */
  {
    NcHIPrim *prim_noreion   = NC_HIPRIM (nc_hiprim_power_law_new ());
    NcHICosmo *cosmo_noreion = NC_HICOSMO (nc_hicosmo_de_xcdm_new_full (NULL, prim_noreion, NULL));
    NcRecomb *recomb         = NC_RECOMB (nc_recomb_seager_new ());

    g_assert_cmpfloat (nc_recomb_get_v_tau_reion_min_z (recomb, cosmo_noreion), ==, 0.0);

    nc_hiprim_free (prim_noreion);
    nc_hicosmo_free (cosmo_noreion);
    nc_recomb_free (recomb);
  }

  nc_hiprim_free (prim);
  nc_hireion_free (reion);
  nc_hicosmo_free (cosmo);
  nc_recomb_free (seager);
  nc_recomb_free (cbe);
}

void
test_nc_recomb_seager_Xe_panels (void)
{
  NcRecomb *recomb = NC_RECOMB (nc_recomb_seager_new ());
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  NcRecomb *tight = NC_RECOMB (nc_recomb_seager_new ());
  guint i;

  /* default prec 1e-7, measured: Saha 7.2e-14, spline 1.6e-9 */
  nc_recomb_prepare_if_needed (recomb, cosmo);
  _test_nc_recomb_seager_Xe_panels_check (recomb, cosmo, TRUE);

  /* against the panels of a prec 1e-11 solution over the whole range, measured: 2.1e-8 */
  g_object_set (tight, "prec", 1.0e-11, NULL);
  nc_recomb_prepare_if_needed (tight, cosmo);
  {
    NcmSpectralPanels *sp   = nc_recomb_seager_peek_Xe_recomb_panels (NC_RECOMB_SEAGER (recomb));
    NcmSpectralPanels *sp_t = nc_recomb_seager_peek_Xe_recomb_panels (NC_RECOMB_SEAGER (tight));

    for (i = 0; i < 400; i++)
    {
      const gdouble lambda = recomb->lambdai * (1.0 - (i + 0.5) / 400.0);

      ncm_assert_cmpdouble_e (ncm_spectral_panels_eval_row (sp, 0, lambda), ==,
                              ncm_spectral_panels_eval_row (sp_t, 0, lambda), recomb->prec, 0.0);
    }
  }

  nc_hicosmo_free (cosmo);
  nc_recomb_free (recomb);
  nc_recomb_free (tight);
}

void
test_nc_recomb_seager_Xe_panels_tight (void)
{
  NcRecomb *recomb = NC_RECOMB (nc_recomb_seager_new ());
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  g_object_set (recomb, "prec", 1.0e-11, NULL);

  /* the build must not abort below the CVODE tolerance (measured Saha 4.6e-15); the spline is
   * not accurate to 1e-11 */
  nc_recomb_prepare_if_needed (recomb, cosmo);
  _test_nc_recomb_seager_Xe_panels_check (recomb, cosmo, FALSE);

  nc_hicosmo_free (cosmo);
  nc_recomb_free (recomb);
}

