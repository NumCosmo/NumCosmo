/***************************************************************************
 *            test_ncm_csq1d.c
 *
 *  Tue September 30 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_csq1d.c
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

#include "test_ncm_csq1d_bessel.h"

static void
test_ncm_csq1d_defaults (void)
{
  /*
   * A subclass implementing only xi, nu and F1 gets nu^2 = nu nu, m = e^xi / nu and
   * F2 = F1' / (2 nu) from NcmCSQ1D; the Bessel system has them in closed form.
   */
  const gdouble a          = 2.0;
  const gdouble k          = 1.3;
  TestCSQ1DBesselMin *bmin = test_csq1d_bessel_min_new (a, k, TRUE);
  NcmCSQ1D *csq1d          = NCM_CSQ1D (bmin);
  gdouble max_F2_err       = 0.0;
  guint i;

  for (i = 0; i < 50; i++)
  {
    const gdouble t = -pow (10.0, -1.0 + 5.0 * i / 49.0);

    ncm_assert_cmpdouble_e (ncm_csq1d_eval_nu2 (csq1d, NULL, t), ==, k * k, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (ncm_csq1d_eval_m (csq1d, NULL, t), ==, test_csq1d_bessel_m (a, k, -1.0, t), 1.0e-13, 0.0);

    max_F2_err = GSL_MAX (max_F2_err, fabs (ncm_csq1d_eval_F2 (csq1d, NULL, t) / test_csq1d_bessel_F2 (a, k, -1.0, t) - 1.0));
  }

  /* Measured 3.0e-13 over t in [-1e4, -0.1]. */
  g_assert_cmpfloat (max_F2_err, <, 3.0e-12);

  ncm_csq1d_free (csq1d);
}

/* Points (alpha, gamma) covering the branches of the distance: small, large alpha, large
 * gamma differences. */
static const gdouble _test_ag[][2] =
{
  {0.0, 0.0}, {0.3, -0.2}, {-0.7, 0.4}, {1.5, 0.1}, {-2.5, -1.3}, {0.05, 3.0}, {4.0, -4.0},
};

/* The hyperbolic distance from the upper half-plane coordinates. */
static gdouble
_test_dist_half_plane (NcmCSQ1DState *s0, NcmCSQ1DState *s1)
{
  gdouble x0, lny0, x1, lny1;

  ncm_csq1d_state_get_poincare_half_plane (s0, &x0, &lny0);
  ncm_csq1d_state_get_poincare_half_plane (s1, &x1, &lny1);

  {
    const gdouble y0 = exp (lny0);
    const gdouble y1 = exp (lny1);

    return acosh (1.0 + (gsl_pow_2 (x0 - x1) + gsl_pow_2 (y0 - y1)) / (2.0 * y0 * y1));
  }
}

/* The hyperbolic distance from the Poincare disc coordinates. */
static gdouble
_test_dist_disc (NcmCSQ1DState *s0, NcmCSQ1DState *s1)
{
  gdouble x0, y0, x1, y1;

  ncm_csq1d_state_get_poincare_disc (s0, &x0, &y0);
  ncm_csq1d_state_get_poincare_disc (s1, &x1, &y1);

  return acosh (1.0 + 2.0 * (gsl_pow_2 (x0 - x1) + gsl_pow_2 (y0 - y1)) / ((1.0 - x0 * x0 - y0 * y0) * (1.0 - x1 * x1 - y1 * y1)));
}

static void
test_ncm_csq1d_state_maps (void)
{
  /*
   * The (chi, U+-) maps invert each other, with chi = sinh alpha and
   * U+- = ln cosh alpha +- gamma. The complex structure J has unit determinant,
   * e^U- = J11 and e^U+ = J22. The vector (phi, P_phi) rebuilds J as
   * J11 = 2 |phi|^2, J22 = 2 |P_phi|^2, J12 = 2 Re (phi P_phi^*), with Wronskian
   * phi P_phi^* - phi^* P_phi = i, in the phase where phi is real and positive.
   */
  NcmCSQ1DState *s = ncm_csq1d_state_new ();
  NcmCSQ1DState *r = ncm_csq1d_state_new ();
  guint i;

  for (i = 0; i < G_N_ELEMENTS (_test_ag); i++)
  {
    const gdouble alpha = _test_ag[i][0];
    const gdouble gamma = _test_ag[i][1];
    gdouble chi, Up, Um, a, g, J11, J12, J22, phi[2], Pphi[2];

    ncm_csq1d_state_set_ag (s, NCM_CSQ1D_FRAME_ORIG, 1.0, alpha, gamma);

    ncm_csq1d_state_get_up (s, &chi, &Up);
    ncm_assert_cmpdouble_e (chi, ==, sinh (alpha), 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (Up, ==, log (cosh (alpha)) + gamma, 1.0e-14, 1.0e-15);
    ncm_csq1d_state_set_up (r, NCM_CSQ1D_FRAME_ORIG, 1.0, chi, Up);
    ncm_csq1d_state_get_ag (r, &a, &g);
    ncm_assert_cmpdouble_e (a, ==, alpha, 1.0e-14, 1.0e-15);
    ncm_assert_cmpdouble_e (g, ==, gamma, 1.0e-14, 1.0e-14);

    ncm_csq1d_state_get_um (s, &chi, &Um);
    ncm_assert_cmpdouble_e (Um, ==, log (cosh (alpha)) - gamma, 1.0e-14, 1.0e-15);
    ncm_csq1d_state_set_um (r, NCM_CSQ1D_FRAME_ORIG, 1.0, chi, Um);
    ncm_csq1d_state_get_ag (r, &a, &g);
    ncm_assert_cmpdouble_e (a, ==, alpha, 1.0e-14, 1.0e-15);
    ncm_assert_cmpdouble_e (g, ==, gamma, 1.0e-14, 1.0e-14);

    ncm_csq1d_state_get_J (s, &J11, &J12, &J22);
    ncm_assert_cmpdouble_e (J11 * J22 - J12 * J12, ==, 1.0, 1.0e-12, 0.0);
    ncm_assert_cmpdouble_e (J11, ==, exp (Um), 1.0e-13, 0.0);
    ncm_assert_cmpdouble_e (J22, ==, exp (Up), 1.0e-13, 0.0);

    ncm_csq1d_state_get_phi_Pphi (s, phi, Pphi);
    ncm_assert_cmpdouble_e (2.0 * (phi[0] * phi[0] + phi[1] * phi[1]), ==, J11, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (2.0 * (Pphi[0] * Pphi[0] + Pphi[1] * Pphi[1]), ==, J22, 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (2.0 * (phi[0] * Pphi[0] + phi[1] * Pphi[1]), ==, J12, 1.0e-14, 1.0e-15);
    ncm_assert_cmpdouble_e (2.0 * (phi[1] * Pphi[0] - phi[0] * Pphi[1]), ==, 1.0, 1.0e-14, 0.0);

    /* phi is real and positive. */
    g_assert_cmpfloat (phi[1], ==, 0.0);
    g_assert_cmpfloat (phi[0], >, 0.0);
  }

  ncm_csq1d_state_free (s);
  ncm_csq1d_state_free (r);
}

static void
test_ncm_csq1d_state_distance (void)
{
  /* The distance agrees with those of the half-plane and disc coordinates, is
   * symmetric, and vanishes between equal points. */
  NcmCSQ1DState *s0 = ncm_csq1d_state_new ();
  NcmCSQ1DState *s1 = ncm_csq1d_state_new ();
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (_test_ag); i++)
  {
    ncm_csq1d_state_set_ag (s0, NCM_CSQ1D_FRAME_ORIG, 1.0, _test_ag[i][0], _test_ag[i][1]);
    g_assert_cmpfloat (ncm_csq1d_state_compute_distance (s0, s0), ==, 0.0);

    for (j = 0; j < G_N_ELEMENTS (_test_ag); j++)
    {
      gdouble d;

      if (i == j)
        continue;

      ncm_csq1d_state_set_ag (s1, NCM_CSQ1D_FRAME_ORIG, 1.0, _test_ag[j][0], _test_ag[j][1]);
      d = ncm_csq1d_state_compute_distance (s0, s1);

      ncm_assert_cmpdouble_e (ncm_csq1d_state_compute_distance (s1, s0), ==, d, 1.0e-14, 0.0);
      /* Measured 6.7e-16 and 6.9e-15. */
      ncm_assert_cmpdouble_e (_test_dist_half_plane (s0, s1), ==, d, 1.0e-14, 0.0);
      ncm_assert_cmpdouble_e (_test_dist_disc (s0, s1), ==, d, 1.0e-13, 0.0);
    }
  }

  ncm_csq1d_state_free (s0);
  ncm_csq1d_state_free (s1);
}

static void
test_ncm_csq1d_state_distance_exact (void)
{
  /*
   * Along alpha at fixed gamma the distance is |delta alpha|; along gamma at fixed
   * alpha it is 2 asinh (cosh (alpha) sinh (delta gamma / 2)). Both down to tiny
   * separations and at large alpha, where cancellations lose them.
   */
  const gdouble alphas[] = {0.0, 0.5, 1.5, 5.0, 50.0, 300.0};
  const gdouble deltas[] = {1.0e-12, 1.0e-8, 1.0e-4, 0.1, 3.0};
  NcmCSQ1DState *s0      = ncm_csq1d_state_new ();
  NcmCSQ1DState *s1      = ncm_csq1d_state_new ();
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (alphas); i++)
  {
    const gdouble alpha = alphas[i];
    const gdouble gamma = -0.4;

    ncm_csq1d_state_set_ag (s0, NCM_CSQ1D_FRAME_ORIG, 1.0, alpha, gamma);
    g_assert_cmpfloat (ncm_csq1d_state_compute_distance (s0, s0), ==, 0.0);

    for (j = 0; j < G_N_ELEMENTS (deltas); j++)
    {
      const gdouble d = deltas[j];

      ncm_csq1d_state_set_ag (s1, NCM_CSQ1D_FRAME_ORIG, 1.0, alpha + d, gamma);
      ncm_assert_cmpdouble_e (ncm_csq1d_state_compute_distance (s0, s1), ==, fabs ((alpha + d) - alpha), 1.0e-12, 0.0);

      ncm_csq1d_state_set_ag (s1, NCM_CSQ1D_FRAME_ORIG, 1.0, alpha, gamma + d);
      ncm_assert_cmpdouble_e (ncm_csq1d_state_compute_distance (s0, s1), ==,
                              2.0 * asinh (cosh (alpha) * sinh (0.5 * ((gamma + d) - gamma))), 1.0e-12, 0.0);
    }
  }

  ncm_csq1d_state_free (s0);
  ncm_csq1d_state_free (s1);
}

static void
test_ncm_csq1d_state_circle (void)
{
  /* The circle of radius r around a point lies at distance r from it; the angle 0
   * moves alpha by r at fixed gamma. */
  const gdouble radii[] = {1.0e-3, 0.3, 2.0};
  NcmCSQ1DState *s      = ncm_csq1d_state_new ();
  NcmCSQ1DState *c      = ncm_csq1d_state_new ();
  guint i, j, l;

  for (i = 0; i < G_N_ELEMENTS (_test_ag); i++)
  {
    ncm_csq1d_state_set_ag (s, NCM_CSQ1D_FRAME_ORIG, 1.0, _test_ag[i][0], _test_ag[i][1]);

    for (j = 0; j < G_N_ELEMENTS (radii); j++)
    {
      gdouble a, g;

      ncm_csq1d_state_get_circle (s, radii[j], 0.0, c);
      ncm_csq1d_state_get_ag (c, &a, &g);
      ncm_assert_cmpdouble_e (a, ==, _test_ag[i][0] + radii[j], 1.0e-13, 1.0e-15);
      ncm_assert_cmpdouble_e (g, ==, _test_ag[i][1], 1.0e-14, 1.0e-15);

      for (l = 0; l < 12; l++)
      {
        ncm_csq1d_state_get_circle (s, radii[j], 2.0 * M_PI * l / 12.0, c);
        /* Measured 3.6e-12, at the smallest radius. */
        ncm_assert_cmpdouble_e (ncm_csq1d_state_compute_distance (s, c), ==, radii[j], 3.0e-11, 0.0);
      }
    }
  }

  ncm_csq1d_state_free (s);
  ncm_csq1d_state_free (c);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/csq1d/defaults", &test_ncm_csq1d_defaults);
  g_test_add_func ("/ncm/csq1d/state/maps", &test_ncm_csq1d_state_maps);
  g_test_add_func ("/ncm/csq1d/state/distance", &test_ncm_csq1d_state_distance);
  g_test_add_func ("/ncm/csq1d/state/distance/exact", &test_ncm_csq1d_state_distance_exact);
  g_test_add_func ("/ncm/csq1d/state/circle", &test_ncm_csq1d_state_circle);

  g_test_run ();
}

