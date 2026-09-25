/***************************************************************************
 *            test_ncm_c.c
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_c.c
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

#include <gsl/gsl_sf_erf.h>

static void
test_ncm_c_math (void)
{
  ncm_assert_cmpdouble_e (ncm_c_pi (), ==, M_PI, 1.0e-16, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_sqrt_pi (), ==, sqrt (M_PI), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_sqrt_2pi (), ==, sqrt (2.0 * M_PI), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_sqrt_pi_2 (), ==, sqrt (M_PI / 2.0), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_sqrt_1_4pi (), ==, sqrt (1.0 / (4.0 * M_PI)), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_sqrt_3_4pi (), ==, sqrt (3.0 / (4.0 * M_PI)), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_ln2 (), ==, M_LN2, 1.0e-16, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_ln3 (), ==, log (3.0), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_lnpi (), ==, log (M_PI), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_lnpi_4 (), ==, log (M_PI) / 4.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_ln2pi (), ==, log (2.0 * M_PI), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_two_pi_2 (), ==, 2.0 * M_PI * M_PI, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_tan_1arcsec (), ==, tan (2.0 * M_PI / (360.0 * 60.0 * 60.0)), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_deg2_steradian (), ==, gsl_pow_2 (M_PI / 180.0), 1.0e-15, 0.0);

  /* The probability within n standard deviations is erf(n / sqrt(2)) */
  ncm_assert_cmpdouble_e (ncm_c_stats_1sigma (), ==, gsl_sf_erf (1.0 / M_SQRT2), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_stats_2sigma (), ==, gsl_sf_erf (2.0 / M_SQRT2), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_stats_3sigma (), ==, gsl_sf_erf (3.0 / M_SQRT2), 1.0e-15, 0.0);
}

static void
test_ncm_c_angles (void)
{
  const gdouble rs[] = {-7.0, -M_PI_2, 0.0, 1.0, M_PI + 0.1, 5.0 * M_PI + 0.2, 100.0};
  guint i;

  ncm_assert_cmpdouble_e (ncm_c_degree_to_radian (180.0), ==, M_PI, 1.0e-16, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_radian_to_degree (M_PI_2), ==, 90.0, 1.0e-16, 0.0);

  for (i = 0; i < G_N_ELEMENTS (rs); i++)
  {
    const gdouble r0 = ncm_c_radian_0_2pi (rs[i]);

    g_assert_cmpfloat (r0, >=, 0.0);
    g_assert_cmpfloat (r0, <, 2.0 * M_PI);
    ncm_assert_cmpdouble_e (sin (r0), ==, sin (rs[i]), 1.0e-13, 1.0e-13);
    g_assert_cmpfloat (ncm_c_sign_sin (rs[i]), ==, (sin (rs[i]) >= 0.0) ? 1.0 : -1.0);
  }
}

static void
test_ncm_c_derived (void)
{
  const gdouble kb = ncm_c_kb ();
  const gdouble h  = ncm_c_h ();
  const gdouble c  = ncm_c_c ();

  /* Astronomical units */
  ncm_assert_cmpdouble_e (ncm_c_pc (), ==, 648000.0 * ncm_c_au () / M_PI, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_Mpc (), ==, 1.0e6 * ncm_c_pc (), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_lightyear (), ==, c * 365.25 * 86400.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_lightyear_pc (), ==, ncm_c_lightyear () / ncm_c_pc (), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_Glightyear_Mpc (), ==, 1.0e9 * ncm_c_lightyear () / ncm_c_Mpc (), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_mass_solar (), ==, ncm_c_G_mass_solar () / ncm_c_G (), 1.0e-15, 0.0);

  /* CODATA relations, exact up to the rounding of the published values */
  ncm_assert_cmpdouble_e (ncm_c_hbar (), ==, h / (2.0 * M_PI), 1.0e-9, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_Ry (), ==, h * c * ncm_c_Rinf (), 1.0e-12, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_electric_constant () * ncm_c_magnetic_constant () * c * c, ==, 1.0, 1.0e-15, 0.0);
  /* The H-I 2p mean is weighted by the statistical weights 2J + 1 = 2, 4 */
  ncm_assert_cmpdouble_e (ncm_c_HI_Lyman_wn_2p_2Pmean (), ==, (ncm_c_HI_Lyman_wn_2p_2P0_5 () + 2.0 * ncm_c_HI_Lyman_wn_2p_2P3_5 ()) / 3.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_HI_ion_wn_2p_2Pmean (), ==, (ncm_c_HI_ion_wn_2p_2P0_5 () + 2.0 * ncm_c_HI_ion_wn_2p_2P3_5 ()) / 3.0, 1.0e-15, 0.0);
  /* RECFAST's Lyman-alpha wavenumber, L_H_alpha */
  ncm_assert_cmpdouble_e (ncm_c_HI_Lyman_wn_2p_2Pmean (), ==, 8.225916453e6, 2.0e-8, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_blackbody_energy_density (), ==, 8.0 * gsl_pow_5 (M_PI) * gsl_pow_4 (kb) / (15.0 * gsl_pow_3 (h * c)), 1.0e-9, 0.0);

  /* Thermal wavelength and wavenumber are inverses */
  ncm_assert_cmpdouble_e (ncm_c_thermal_wl_e () * ncm_c_thermal_wn_e (), ==, 1.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_thermal_wl_p () * ncm_c_thermal_wn_p (), ==, 1.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_thermal_wl_n () * ncm_c_thermal_wn_n (), ==, 1.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_thermal_wl_e (), ==, sqrt (2.0 * M_PI * gsl_pow_2 (ncm_c_hbar ()) / (ncm_c_mass_e () * kb)), 1.0e-15, 0.0);

  /* Critical density for H0 = 100 km/s/Mpc */
  {
    const gdouble H0 = 1.0e5 / ncm_c_Mpc ();

    ncm_assert_cmpdouble_e (ncm_c_crit_density_h2 (), ==, 3.0 * c * c * H0 * H0 / (8.0 * M_PI * ncm_c_G ()), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_c_crit_mass_density_h2 (), ==, ncm_c_crit_density_h2 () / (c * c), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (ncm_c_hubble_radius_hm1_Mpc (), ==, c / 1.0e5, 1.0e-15, 0.0);
  }

  /* Ionization wavenumbers from the ground level minus the transition */
  ncm_assert_cmpdouble_e (ncm_c_HI_ion_wn_2s_2S0_5 (), ==, ncm_c_HI_ion_wn_1s_2S0_5 () - ncm_c_HI_Lyman_wn_2s_2S0_5 (), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_HI_ion_E_1s_2S0_5 (), ==, ncm_c_HI_ion_wn_1s_2S0_5 () * h * c, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_c_HI_Lyman_wl_2s_2S0_5 () * ncm_c_HI_Lyman_wn_2s_2S0_5 (), ==, 1.0, 1.0e-15, 0.0);
}

static void
test_ncm_c_H_bind (void)
{
  /* {n, j, NIST ionization energy / hc [m^-1], relative envelope of the missing Lamb shift} */
  const gdouble levels[][4] = {
    {1.0, 0.5, 1.0967877174307e7, 3.0e-6},
    {2.0, 0.5, 1.0967877174307e7 - 8.22589543992821e6, 1.5e-6},
    {2.0, 0.5, 1.0967877174307e7 - 8.22589191133e6, 3.0e-8},
    {2.0, 1.5, 1.0967877174307e7 - 8.22592850014e6, 3.0e-8},
  };
  guint i;

  for (i = 0; i < G_N_ELEMENTS (levels); i++)
  {
    const gdouble E_nist = levels[i][2] * ncm_c_hc ();

    ncm_assert_cmpdouble_e (ncm_c_H_bind (levels[i][0], levels[i][1]), ==, E_nist, levels[i][3], 0.0);
  }
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/c/math", &test_ncm_c_math);
  g_test_add_func ("/ncm/c/angles", &test_ncm_c_angles);
  g_test_add_func ("/ncm/c/derived", &test_ncm_c_derived);
  g_test_add_func ("/ncm/c/H_bind", &test_ncm_c_H_bind);

  g_test_run ();
}

