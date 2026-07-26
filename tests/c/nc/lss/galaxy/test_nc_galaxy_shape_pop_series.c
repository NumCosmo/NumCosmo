/***************************************************************************
 *            test_nc_galaxy_shape_pop_series.c
 *
 *  Wed Jul 15 2026
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

/* Checks NcGalaxyShapePop's eval_p_rho2_g_series implementations (Gauss and
 * Beta) directly against their own exact, non-series eval_p(): despite its
 * name (shared with the eval_p_rho2/eval_p_rho2 vfunc pair), the series
 * eval_p_rho2_g_series composes with is x(g)=|chi_I(chi_L,g)|^2 itself (see
 * its own doc comment in nc_galaxy_shape_pop.h) -- the same variable eval_p()
 * takes directly, not the disc-compactified rho^2=x/(1-x) that eval_p_rho2()
 * uses. So for a synthetic x(g) truncated power series and small real g
 * inside the truncation radius, eval_p_rho2_g_series(x_series)
 * Horner-evaluated at g must match eval_p(x(g)) to O(g^(order+1)) -- the
 * population composition layer alone, decoupled from the shear-map series
 * (nc_wl_ellipticity_series.c, checked independently against
 * nc_wl_ellipticity.h in its own test file). */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <complex.h>
#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_ORDER 8
#define TEST_TOL   1.0e-9

/* A synthetic x(g) = r0 + r1*g + r2*g^2 (r0>0, higher terms zero): every
 * coefficient a plain real scalar (harmonic 0 only), matching how
 * x(g)=|chi_I(chi_L,g)|^2 always looks in the real pipeline. */
static NcmLaurentSeriesTPS *
_make_x_series (gdouble r0, gdouble r1, gdouble r2)
{
  NcmLaurentSeriesTPS *tps = ncm_laurent_series_tps_new (TEST_ORDER);
  const gdouble coeffs[3]  = {r0, r1, r2};
  guint n;

  for (n = 0; n <= TEST_ORDER; n++)
  {
    NcmLaurentSeries *slot = ncm_laurent_series_tps_get (tps, n);
    const gdouble c        = (n < 3) ? coeffs[n] : 0.0;

    ncm_laurent_series_set_single_into (slot, 0, c);
  }

  return tps;
}

/* Every coefficient here is harmonic-0 only, so evaluating at w=1 just
 * reads off the (real) scalar value -- no angular dependence to worry
 * about, unlike the shear-map series. */
static gdouble
_eval_real (NcmLaurentSeriesTPS *tps, gdouble g)
{
  return creal (ncm_laurent_series_tps_eval (tps, 1.0, g));
}

static const gdouble test_gs[] = {0.0, 0.01, -0.015, 0.02};

void test_nc_galaxy_shape_pop_gauss_series_matches_exact (void);
void test_nc_galaxy_shape_pop_beta_series_matches_exact (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/nc/galaxy_shape_pop_series/gauss/matches_exact",
                   &test_nc_galaxy_shape_pop_gauss_series_matches_exact);
  g_test_add_func ("/nc/galaxy_shape_pop_series/beta/matches_exact",
                   &test_nc_galaxy_shape_pop_beta_series_matches_exact);

  return g_test_run ();
}

void
test_nc_galaxy_shape_pop_gauss_series_matches_exact (void)
{
  NcGalaxyShapePopGauss *pop = nc_galaxy_shape_pop_gauss_new ();
  NcGalaxyShapePopData *data = nc_galaxy_shape_pop_data_new (NC_GALAXY_SHAPE_POP (pop));
  NcmLaurentSeriesTPS *x     = _make_x_series (0.05, 0.10, -0.03);
  NcmLaurentSeriesTPS *out   = ncm_laurent_series_tps_new (TEST_ORDER);
  guint i;

  ncm_model_param_set_by_name (NCM_MODEL (pop), "sigma", 0.25, NULL);
  nc_galaxy_shape_pop_prepare (NC_GALAXY_SHAPE_POP (pop), data);

  nc_galaxy_shape_pop_eval_p_rho2_g_series (NC_GALAXY_SHAPE_POP (pop), data, x, out);

  for (i = 0; i < G_N_ELEMENTS (test_gs); i++)
  {
    const gdouble g          = test_gs[i];
    const gdouble x_g        = _eval_real (x, g);
    const gdouble series_val = _eval_real (out, g);
    const gdouble exact_val  = nc_galaxy_shape_pop_eval_p (NC_GALAXY_SHAPE_POP (pop), data, x_g);

    g_assert_cmpfloat (fabs (series_val - exact_val), <, TEST_TOL);
  }

  ncm_laurent_series_tps_unref (x);
  ncm_laurent_series_tps_unref (out);
  nc_galaxy_shape_pop_data_unref (data);
  nc_galaxy_shape_pop_gauss_free (pop);
}

void
test_nc_galaxy_shape_pop_beta_series_matches_exact (void)
{
  NcGalaxyShapePopBeta *pop  = nc_galaxy_shape_pop_beta_new ();
  NcGalaxyShapePopData *data = nc_galaxy_shape_pop_data_new (NC_GALAXY_SHAPE_POP (pop));
  NcmLaurentSeriesTPS *x     = _make_x_series (0.08, 0.06, 0.02);
  NcmLaurentSeriesTPS *out   = ncm_laurent_series_tps_new (TEST_ORDER);
  guint i;

  /* alpha=0.9, beta=4.1 (the class's own pre->=1-bound default, still
   * directly settable here): a non-integer alpha (alpha<1) deliberately
   * exercises the "real, non-integer exponent" path
   * ncm_laurent_series_tps_pow() exists for. */
  ncm_model_param_set_by_name (NCM_MODEL (pop), "alpha", 0.9, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (pop), "beta", 4.1, NULL);
  nc_galaxy_shape_pop_prepare (NC_GALAXY_SHAPE_POP (pop), data);

  nc_galaxy_shape_pop_eval_p_rho2_g_series (NC_GALAXY_SHAPE_POP (pop), data, x, out);

  for (i = 0; i < G_N_ELEMENTS (test_gs); i++)
  {
    const gdouble g          = test_gs[i];
    const gdouble x_g        = _eval_real (x, g);
    const gdouble series_val = _eval_real (out, g);
    const gdouble exact_val  = nc_galaxy_shape_pop_eval_p (NC_GALAXY_SHAPE_POP (pop), data, x_g);

    g_assert_cmpfloat (fabs (series_val - exact_val), <, TEST_TOL);
  }

  ncm_laurent_series_tps_unref (x);
  ncm_laurent_series_tps_unref (out);
  nc_galaxy_shape_pop_data_unref (data);
  nc_galaxy_shape_pop_beta_free (pop);
}

