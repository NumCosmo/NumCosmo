/***************************************************************************
 *            test_nc_galaxy_shape_intrinsic_mode.c
 *
 *  Mon Jul 6 2026
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

/* Regression guard for NcGalaxyShapeIntrinsicMode's closed-form mode-finders
 * (nc_galaxy_shape_intrinsic_mode_find_trace_det()/_trace()) against the
 * original, convention-agnostic finite-difference one
 * (nc_galaxy_shape_intrinsic_mode_find()). Not introspectable (this whole
 * API is NUMCOSMO_GIR_SCAN-guarded), so this has to be a C-level test
 * rather than a Python one -- see the class docs in
 * nc_galaxy_shape_intrinsic_mode.c for the derivation (a Mobius/Blaschke
 * circle-image closed form for TRACE_DET's theta-profile, a rotation-
 * gauge-fixed quartic root for TRACE's, both paired with an
 * envelope-theorem closed form for the rho-Newton derivatives).
 *
 * IMPORTANT ASYMMETRY, found during development: TRACE_DET's closed-form
 * theta-profile is the EXACT global argmin at every rho (a provable
 * geometric fact, see the circle-distance derivation), so it agrees with
 * the finite-difference baseline almost everywhere and is tested with a
 * tight tolerance below. TRACE's theta-profile is also exact, but because
 * it's a genuinely broader/noisier convention in practice, the old nested
 * finite-difference Newton occasionally converges to the wrong local
 * optimum (verified against a 2000x2000 brute-force grid: in every
 * observed disagreement, the CLOSED FORM matched the true global mode and
 * the finite-difference one did not) -- so TRACE's test asserts the
 * closed form is never meaningfully WORSE (lower ln_peak) than the
 * finite-difference one, rather than asserting close agreement, which
 * would be the wrong invariant to enforce here. */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>
#include <numcosmo/nc/lss/galaxy/nc_galaxy_shape_intrinsic_mode.h>

#include <complex.h>
#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_math.h>

/* Tolerances set from an empirical 5000-trial sweep during development
 * (max relative Laplace-value difference ~6e-5, max NAN-mismatch rate
 * ~1/5000, both isolated to broad/high-noise cases already documented as
 * degrading the Laplace approximation itself) -- kept an order of magnitude
 * looser here to absorb compiler/platform variation without becoming a
 * flaky test. */
#define TEST_N_TRIALS          2000
#define TEST_MAX_REL_LAPLACE   1.0e-3
#define TEST_MAX_NAN_MISMATCH_FRAC 0.01

static gdouble
test_uniform (GRand *rng, gdouble a, gdouble b)
{
  return a + (b - a) * g_rand_double (rng);
}

static void
test_closed_form_matches_finite_difference (void)
{
  GRand *rng                       = g_rand_new_with_seed (20260706);
  NcGalaxyShapePopGauss *pop_gauss = nc_galaxy_shape_pop_gauss_new ();
  NcGalaxyShapePopBeta *pop_beta   = nc_galaxy_shape_pop_beta_new ();
  guint n_nan_mismatch             = 0;
  gdouble max_rel_laplace          = 0.0;
  guint i;

  ncm_model_orig_param_set (NCM_MODEL (pop_gauss), 0, 0.28);
  ncm_model_orig_param_set (NCM_MODEL (pop_beta), 0, 0.35);
  ncm_model_orig_param_set (NCM_MODEL (pop_beta), 1, 6.0);

  for (i = 0; i < TEST_N_TRIALS; i++)
  {
    NcGalaxyShapePop *pop = (i % 2 == 0) ?
                            NC_GALAXY_SHAPE_POP (pop_gauss) : NC_GALAXY_SHAPE_POP (pop_beta);
    NcGalaxyShapePopData *pop_data = nc_galaxy_shape_pop_data_new (pop);
    const complex double g         = test_uniform (rng, -0.4, 0.4) + I * test_uniform (rng, -0.4, 0.4);
    const complex double eps_obs   = test_uniform (rng, -0.5, 0.5) + I * test_uniform (rng, -0.5, 0.5);
    const gdouble std_noise        = test_uniform (rng, 0.03, 0.35);
    NcGalaxyShapeIntrinsicMode mode_old, mode_new;
    gboolean old_ok, new_ok;

    nc_galaxy_shape_pop_prepare (pop, pop_data);

    nc_galaxy_shape_intrinsic_mode_find (&nc_wl_ellipticity_apply_shear_trace_det,
                                         &nc_wl_ellipticity_apply_shear_inv_trace_det,
                                         pop, pop_data, g, eps_obs, std_noise, &mode_old);
    nc_galaxy_shape_intrinsic_mode_find_trace_det (pop, pop_data, g, eps_obs, std_noise, &mode_new);

    old_ok = gsl_finite (nc_galaxy_shape_intrinsic_mode_laplace (&mode_old));
    new_ok = gsl_finite (nc_galaxy_shape_intrinsic_mode_laplace (&mode_new));

    if (old_ok != new_ok)
    {
      n_nan_mismatch++;
    }
    else if (old_ok && new_ok)
    {
      const gdouble lap_old = nc_galaxy_shape_intrinsic_mode_laplace (&mode_old);
      const gdouble lap_new = nc_galaxy_shape_intrinsic_mode_laplace (&mode_new);
      const gdouble rel     = fabs (lap_old - lap_new) / fabs (lap_old);

      max_rel_laplace = MAX (max_rel_laplace, rel);
      g_assert_cmpfloat (rel, <, TEST_MAX_REL_LAPLACE);
    }

    nc_galaxy_shape_pop_data_unref (pop_data);
  }

  g_assert_cmpfloat ((gdouble) n_nan_mismatch / TEST_N_TRIALS, <, TEST_MAX_NAN_MISMATCH_FRAC);

  nc_galaxy_shape_pop_gauss_free (pop_gauss);
  nc_galaxy_shape_pop_beta_free (pop_beta);
  g_rand_free (rng);
}

/* g=0: the shear map is the identity for both TRACE and TRACE_DET, so the
 * closest-point-on-circle construction degenerates to the trivial case
 * (center 0, radius rho) -- theta_hat = arg(eps_obs) exactly, independent
 * of rho and of the population. An easy, independently-derivable sanity
 * check beyond "agrees with the finite-difference method". */
static void
test_zero_shear_theta_matches_arg_eps_obs (void)
{
  NcGalaxyShapePopGauss *pop     = nc_galaxy_shape_pop_gauss_new ();
  NcGalaxyShapePopData *pop_data = nc_galaxy_shape_pop_data_new (NC_GALAXY_SHAPE_POP (pop));
  const complex double g         = 0.0;
  const complex double eps_obs   = 0.25 - 0.18 * I;
  NcGalaxyShapeIntrinsicMode mode;

  ncm_model_orig_param_set (NCM_MODEL (pop), 0, 0.3);
  nc_galaxy_shape_pop_prepare (NC_GALAXY_SHAPE_POP (pop), pop_data);

  nc_galaxy_shape_intrinsic_mode_find_trace_det (NC_GALAXY_SHAPE_POP (pop), pop_data, g, eps_obs, 0.2, &mode);

  g_assert_cmpfloat (fabs (remainder (mode.theta - carg (eps_obs), 2.0 * M_PI)), <, 1.0e-8);

  nc_galaxy_shape_pop_data_unref (pop_data);
  nc_galaxy_shape_pop_gauss_free (pop);
}

/* TRACE never finds a meaningfully worse mode than the finite-difference
 * baseline -- see the file-level comment on this asymmetry. A small
 * negative tolerance absorbs the two methods occasionally landing on
 * bit-different but equally-good points. */
#define TEST_TRACE_MIN_LN_PEAK_SLACK -1.0e-6

static void
test_trace_closed_form_never_worse_than_finite_difference (void)
{
  GRand *rng                       = g_rand_new_with_seed (20260706);
  NcGalaxyShapePopGauss *pop_gauss = nc_galaxy_shape_pop_gauss_new ();
  NcGalaxyShapePopBeta *pop_beta   = nc_galaxy_shape_pop_beta_new ();
  guint n_nan_mismatch             = 0;
  guint i;

  ncm_model_orig_param_set (NCM_MODEL (pop_gauss), 0, 0.28);
  ncm_model_orig_param_set (NCM_MODEL (pop_beta), 0, 0.35);
  ncm_model_orig_param_set (NCM_MODEL (pop_beta), 1, 6.0);

  for (i = 0; i < TEST_N_TRIALS; i++)
  {
    NcGalaxyShapePop *pop = (i % 2 == 0) ?
                            NC_GALAXY_SHAPE_POP (pop_gauss) : NC_GALAXY_SHAPE_POP (pop_beta);
    NcGalaxyShapePopData *pop_data = nc_galaxy_shape_pop_data_new (pop);
    const complex double g         = test_uniform (rng, -0.4, 0.4) + I * test_uniform (rng, -0.4, 0.4);
    const complex double eps_obs   = test_uniform (rng, -0.5, 0.5) + I * test_uniform (rng, -0.5, 0.5);
    const gdouble std_noise        = test_uniform (rng, 0.03, 0.35);
    NcGalaxyShapeIntrinsicMode mode_old, mode_new;
    gboolean old_ok, new_ok;

    nc_galaxy_shape_pop_prepare (pop, pop_data);

    nc_galaxy_shape_intrinsic_mode_find (&nc_wl_ellipticity_apply_shear_trace,
                                         &nc_wl_ellipticity_apply_shear_inv_trace,
                                         pop, pop_data, g, eps_obs, std_noise, &mode_old);
    nc_galaxy_shape_intrinsic_mode_find_trace (pop, pop_data, g, eps_obs, std_noise, &mode_new);

    old_ok = gsl_finite (nc_galaxy_shape_intrinsic_mode_laplace (&mode_old));
    new_ok = gsl_finite (nc_galaxy_shape_intrinsic_mode_laplace (&mode_new));

    if (old_ok != new_ok)
      n_nan_mismatch++;
    else if (old_ok && new_ok)
      g_assert_cmpfloat (mode_new.ln_peak - mode_old.ln_peak, >, TEST_TRACE_MIN_LN_PEAK_SLACK);

    nc_galaxy_shape_pop_data_unref (pop_data);
  }

  g_assert_cmpfloat ((gdouble) n_nan_mismatch / TEST_N_TRIALS, <, TEST_MAX_NAN_MISMATCH_FRAC);

  nc_galaxy_shape_pop_gauss_free (pop_gauss);
  nc_galaxy_shape_pop_beta_free (pop_beta);
  g_rand_free (rng);
}

static void
test_trace_zero_shear_theta_matches_arg_eps_obs (void)
{
  NcGalaxyShapePopGauss *pop     = nc_galaxy_shape_pop_gauss_new ();
  NcGalaxyShapePopData *pop_data = nc_galaxy_shape_pop_data_new (NC_GALAXY_SHAPE_POP (pop));
  const complex double g         = 0.0;
  const complex double eps_obs   = 0.25 - 0.18 * I;
  NcGalaxyShapeIntrinsicMode mode;

  ncm_model_orig_param_set (NCM_MODEL (pop), 0, 0.3);
  nc_galaxy_shape_pop_prepare (NC_GALAXY_SHAPE_POP (pop), pop_data);

  nc_galaxy_shape_intrinsic_mode_find_trace (NC_GALAXY_SHAPE_POP (pop), pop_data, g, eps_obs, 0.2, &mode);

  g_assert_cmpfloat (fabs (remainder (mode.theta - carg (eps_obs), 2.0 * M_PI)), <, 1.0e-8);

  nc_galaxy_shape_pop_data_unref (pop_data);
  nc_galaxy_shape_pop_gauss_free (pop);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/nc/galaxy_shape_intrinsic_mode/trace_det/closed_form_matches_finite_difference",
                   &test_closed_form_matches_finite_difference);
  g_test_add_func ("/nc/galaxy_shape_intrinsic_mode/trace_det/zero_shear_theta_matches_arg_eps_obs",
                   &test_zero_shear_theta_matches_arg_eps_obs);
  g_test_add_func ("/nc/galaxy_shape_intrinsic_mode/trace/closed_form_never_worse_than_finite_difference",
                   &test_trace_closed_form_never_worse_than_finite_difference);
  g_test_add_func ("/nc/galaxy_shape_intrinsic_mode/trace/zero_shear_theta_matches_arg_eps_obs",
                   &test_trace_zero_shear_theta_matches_arg_eps_obs);

  return g_test_run ();
}

