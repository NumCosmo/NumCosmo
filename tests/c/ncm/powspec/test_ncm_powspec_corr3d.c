/***************************************************************************
 *            test_ncm_powspec_corr3d.c
 *
 *  Mon September 28 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

/* NcmPowspecAnalytic (BBKS, LCDM growth) on k in [1e-6, 1e2] and z in [0, 2]. */
static NcmPowspec *
_test_powspec_new (void)
{
  NcmPowspec *ps = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                          NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));

  ncm_powspec_set_kmin (ps, 1.0e-6);
  ncm_powspec_set_kmax (ps, 1.0e2);
  ncm_powspec_set_zi (ps, 0.0);
  ncm_powspec_set_zf (ps, 2.0);

  return ps;
}

static NcmPowspecCorr3d *
_test_corr3d_new (NcmPowspec *ps, const gdouble reltol)
{
  NcmPowspecCorr3d *psc = ncm_powspec_corr3d_new (ps);

  ncm_powspec_corr3d_set_reltol (psc, reltol);
  ncm_powspec_corr3d_prepare (psc, NULL);

  return psc;
}

/* Largest |xi| over 400 points of [r_lo, r_hi] at z = 0. */
static gdouble
_test_peak (NcmPowspecCorr3d *psc, const gdouble r_lo, const gdouble r_hi)
{
  gdouble peak = 0.0;
  guint i;

  for (i = 0; i < 400; i++)
    peak = GSL_MAX (peak, fabs (ncm_powspec_corr3d_eval_xi (psc, 0.0, r_lo * pow (r_hi / r_lo, i / 399.0))));

  return peak;
}

/*
 * The grid at reltol 1e-5 against one at 1e-9, on their common r range, relative
 * to the peak of |xi|: within reltol (measured 1.4e-6).
 */
static void
test_ncm_powspec_corr3d_calibration (void)
{
  NcmPowspec *ps        = _test_powspec_new ();
  NcmPowspecCorr3d *psc = _test_corr3d_new (ps, 1.0e-5);
  NcmPowspecCorr3d *ref = _test_corr3d_new (ps, 1.0e-9);
  const gdouble r_lo    = GSL_MAX (ncm_powspec_corr3d_get_r_min (psc), ncm_powspec_corr3d_get_r_min (ref));
  const gdouble r_hi    = GSL_MIN (ncm_powspec_corr3d_get_r_max (psc), ncm_powspec_corr3d_get_r_max (ref));
  const gdouble peak    = _test_peak (ref, r_lo, r_hi);
  gdouble err           = 0.0;
  guint i;

  for (i = 0; i < 400; i++)
  {
    const gdouble r = r_lo * pow (r_hi / r_lo, i / 399.0);

    err = GSL_MAX (err, fabs (ncm_powspec_corr3d_eval_xi (psc, 0.0, r) - ncm_powspec_corr3d_eval_xi (ref, 0.0, r)));
  }

  g_assert_cmpfloat (err / peak, <, 1.0e-5);

  ncm_powspec_corr3d_free (psc);
  ncm_powspec_corr3d_free (ref);
  ncm_powspec_free (ps);
}

/*
 * Against the quadrature over the table, ncm_powspec_corr3d(), for r >= 100 / kmax,
 * where the continuation past the table contributes little: measured 2.8e-5 of the
 * peak of |xi|. Below that the two differ by design (2.9e-3 at r = 0.1).
 */
static void
test_ncm_powspec_corr3d_quadrature (void)
{
  NcmPowspec *ps        = _test_powspec_new ();
  NcmPowspecCorr3d *psc = _test_corr3d_new (ps, 1.0e-9);
  const gdouble peak    = _test_peak (psc, ncm_powspec_corr3d_get_r_min (psc), ncm_powspec_corr3d_get_r_max (psc));
  guint i;

  for (i = 0; i < 30; i++)
  {
    const gdouble r = 1.0 * pow (1.0e3, i / 29.0);

    ncm_assert_cmpdouble_e (ncm_powspec_corr3d_eval_xi (psc, 0.0, r), ==, ncm_powspec_corr3d (ps, NULL, 1.0e-8, 0.0, r), 0.0, 5.0e-5 * peak);
  }

  ncm_powspec_corr3d_free (psc);
  ncm_powspec_free (ps);
}

/* xi(r, z) = xi(r, 0) D(z)^2 for the separable spectrum; within reltol-z 1e-6 (measured 4.2e-8). */
static void
test_ncm_powspec_corr3d_redshift (void)
{
  NcmPowspec *ps        = _test_powspec_new ();
  NcmPowspecCorr3d *psc = _test_corr3d_new (ps, 1.0e-7);
  const gdouble rs[]    = {1.0, 10.0, 100.0};
  guint i, j;

  for (i = 0; i <= 22; i++)
  {
    const gdouble z  = 2.0 * i / 22.0;
    const gdouble D2 = gsl_pow_2 (ncm_powspec_analytic_eval_growth (NCM_POWSPEC_ANALYTIC (ps), z));

    for (j = 0; j < G_N_ELEMENTS (rs); j++)
      ncm_assert_cmpdouble_e (ncm_powspec_corr3d_eval_xi (psc, z, rs[j]), ==, ncm_powspec_corr3d_eval_xi (psc, 0.0, rs[j]) * D2,
                              ncm_powspec_corr3d_get_reltol_z (psc), 0.0);
  }

  ncm_powspec_corr3d_free (psc);
  ncm_powspec_free (ps);
}

/* A new reltol recalibrates through prepare_if_needed; preparing again keeps the values;
 * the reported ends evaluate. */
static void
test_ncm_powspec_corr3d_reprepare (void)
{
  NcmPowspec *ps        = _test_powspec_new ();
  NcmPowspecCorr3d *psc = _test_corr3d_new (ps, 1.0e-3);
  const gdouble r_min   = ncm_powspec_corr3d_get_r_min (psc);
  const gdouble xi_10   = ncm_powspec_corr3d_eval_xi (psc, 0.5, 10.0);

  ncm_powspec_corr3d_prepare (psc, NULL);
  g_assert_cmpfloat (ncm_powspec_corr3d_eval_xi (psc, 0.5, 10.0), ==, xi_10);

  /* The ends reported are inside the grid. */
  g_assert_true (isfinite (ncm_powspec_corr3d_eval_xi (psc, 0.0, ncm_powspec_corr3d_get_r_min (psc))));
  g_assert_true (isfinite (ncm_powspec_corr3d_eval_xi (psc, 2.0, ncm_powspec_corr3d_get_r_max (psc))));

  ncm_powspec_corr3d_set_reltol (psc, 1.0e-9);
  g_assert_cmpfloat (ncm_powspec_corr3d_get_reltol (psc), ==, 1.0e-9);
  ncm_powspec_corr3d_prepare_if_needed (psc, NULL);
  g_assert_cmpfloat (ncm_powspec_corr3d_get_r_min (psc), !=, r_min);

  ncm_powspec_corr3d_free (psc);
  ncm_powspec_free (ps);
}

static void
test_ncm_powspec_corr3d_outside (void)
{
  g_test_trap_subprocess ("/ncm/powspec_corr3d/outside/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*is outside the grid*");
}

static void
test_ncm_powspec_corr3d_outside_subprocess (void)
{
  NcmPowspec *ps        = _test_powspec_new ();
  NcmPowspecCorr3d *psc = _test_corr3d_new (ps, 1.0e-3);

  ncm_powspec_corr3d_eval_xi (psc, 0.0, ncm_powspec_corr3d_get_r_max (psc) * 1.01);
}

static void
test_ncm_powspec_corr3d_unprepared (void)
{
  g_test_trap_subprocess ("/ncm/powspec_corr3d/unprepared/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*exists only after*");
}

static void
test_ncm_powspec_corr3d_unprepared_subprocess (void)
{
  NcmPowspec *ps        = _test_powspec_new ();
  NcmPowspecCorr3d *psc = ncm_powspec_corr3d_new (ps);

  ncm_powspec_corr3d_get_r_min (psc);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/powspec_corr3d/calibration", &test_ncm_powspec_corr3d_calibration);
  g_test_add_func ("/ncm/powspec_corr3d/quadrature", &test_ncm_powspec_corr3d_quadrature);
  g_test_add_func ("/ncm/powspec_corr3d/redshift", &test_ncm_powspec_corr3d_redshift);
  g_test_add_func ("/ncm/powspec_corr3d/reprepare", &test_ncm_powspec_corr3d_reprepare);
  g_test_add_func ("/ncm/powspec_corr3d/outside", &test_ncm_powspec_corr3d_outside);
  g_test_add_func ("/ncm/powspec_corr3d/outside/subprocess", &test_ncm_powspec_corr3d_outside_subprocess);
  g_test_add_func ("/ncm/powspec_corr3d/unprepared", &test_ncm_powspec_corr3d_unprepared);
  g_test_add_func ("/ncm/powspec_corr3d/unprepared/subprocess", &test_ncm_powspec_corr3d_unprepared_subprocess);

  g_test_run ();
}

