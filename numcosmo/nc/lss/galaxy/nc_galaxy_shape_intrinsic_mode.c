/***************************************************************************
 *            nc_galaxy_shape_intrinsic_mode.c
 *
 *  Fri Jul 3 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * nc_galaxy_shape_intrinsic_mode.c
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
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:nc_galaxy_shape_intrinsic_mode
 *
 * Finds the joint mode (rho, theta) of
 * $P_\mathrm{pop}(\chi_I)\,N_2(\epsilon_\mathrm{obs}-f_g(\chi_I);\sigma_\mathrm{noise}^2)$
 * over the disc $\chi_I=\rho e^{i\theta}$, $\rho\in[0,1)$, plus the 2x2
 * Hessian of its log in Cartesian $\chi_I$ coordinates at that point.
 *
 * The search is nested rather than a blind joint 2D minimization: since
 * $P_\mathrm{pop}(\chi_I)=P_\mathrm{pop}(\rho^2)$ has no $\theta$ dependence,
 * the argmin over $\theta$ at ANY fixed $\rho$ comes entirely from the noise
 * term (a fixed 3-step Newton refinement on
 * $|\epsilon_\mathrm{obs}-f_g(\rho e^{i\theta})|^2$, the same scheme already
 * used by #NcGalaxyShapeFactorQuad's Beta-mode hint, generalized here to
 * arbitrary $\rho$, not just the population's own marginal mode). A blind
 * joint 2D minimizer cannot exploit this: it would needlessly re-evaluate
 * $P_\mathrm{pop}$ (a special-function call for e.g. #NcGalaxyShapePopBeta)
 * at every stencil point even when only $\theta$ moved. Nesting the search
 * instead means $P_\mathrm{pop}$ is only ever evaluated along the outer
 * (cheap-to-profile) $\rho$ direction.
 *
 * Both Newton loops are a small fixed number of steps with a concavity/
 * curvature guard (fall back to the current estimate rather than diverge),
 * matching the project's existing convention for these one-shot,
 * once-per-evaluation hint computations (see
 * nc_galaxy_shape_factor_quad.c's own theta-only refiner) -- this is a peak
 * finder for a single, generally well-behaved (unimodal) bump, not a general-
 * purpose robust optimizer.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/lss/galaxy/nc_galaxy_shape_intrinsic_mode.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_math.h>
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct _ModeCtx
{
  complex double (*apply_shear) (complex double, complex double);
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *pop_data;
  complex double g;
  complex double eps_obs;
  gdouble noise_var;
} ModeCtx;

static gdouble
_ln_integrand (const ModeCtx *ctx, const gdouble rho, const gdouble theta)
{
  const complex double chi_i    = rho * cexp (I * theta);
  const gdouble x_i             = rho * rho;
  const gdouble p_pop           = nc_galaxy_shape_pop_eval_p (ctx->pop, ctx->pop_data, x_i) / M_PI;
  const complex double eps_true = ctx->apply_shear (ctx->g, chi_i);
  const complex double delta    = ctx->eps_obs - eps_true;
  const gdouble noise           = exp (-(gsl_pow_2 (creal (delta)) + gsl_pow_2 (cimag (delta))) / (2.0 * ctx->noise_var)) / (2.0 * M_PI * ctx->noise_var);
  const gdouble val             = p_pop * noise;

  return (val > 0.0) && isfinite (val) ? log (val) : -G_MAXDOUBLE;
}

/* argmin over theta of |eps_obs - f_g(rho e^{i theta})|^2 at fixed rho: the
 * whole theta-dependence of the integrand, since P_pop is theta-independent
 * (see the file doc). */
static gdouble
_theta_hat (const ModeCtx *ctx, const gdouble rho, const gdouble theta0)
{
  const gdouble h = 1.0e-4;
  gdouble theta   = theta0;
  guint i;

  for (i = 0; i < 3; i++)
  {
    const complex double delta_p = ctx->eps_obs - ctx->apply_shear (ctx->g, rho * cexp (I * (theta + h)));
    const complex double delta_0 = ctx->eps_obs - ctx->apply_shear (ctx->g, rho * cexp (I * theta));
    const complex double delta_m = ctx->eps_obs - ctx->apply_shear (ctx->g, rho * cexp (I * (theta - h)));
    const gdouble fp             = gsl_pow_2 (creal (delta_p)) + gsl_pow_2 (cimag (delta_p));
    const gdouble f0             = gsl_pow_2 (creal (delta_0)) + gsl_pow_2 (cimag (delta_0));
    const gdouble fm             = gsl_pow_2 (creal (delta_m)) + gsl_pow_2 (cimag (delta_m));
    const gdouble grad           = (fp - fm) / (2.0 * h);
    const gdouble hess           = (fp - 2.0 * f0 + fm) / gsl_pow_2 (h);

    if (hess <= 0.0)
      break;

    theta -= grad / hess;
  }

  return theta;
}

/* ln P_pop(rho^2) + ln N(...) profiled over its own theta_hat(rho). */
static gdouble
_profile_h (const ModeCtx *ctx, const gdouble rho, const gdouble theta_seed, gdouble *theta_out)
{
  const gdouble theta = (rho > 0.0) ? _theta_hat (ctx, rho, theta_seed) : theta_seed;

  if (theta_out != NULL)
    *theta_out = theta;

  return _ln_integrand (ctx, rho, theta);
}

/* Outer 1D Newton search over rho in (0,1), nesting _profile_h() at every
 * trial point. */
static gdouble
_rho_hat (const ModeCtx *ctx, const gdouble rho0, const gdouble theta_seed, gdouble *theta_out)
{
  const gdouble h = 1.0e-4;
  gdouble rho     = rho0;
  guint i;

  for (i = 0; i < 6; i++)
  {
    const gdouble rho_p = MIN (rho + h, 1.0 - 1.0e-6);
    const gdouble rho_m = MAX (rho - h, 1.0e-6);
    const gdouble fp    = _profile_h (ctx, rho_p, theta_seed, NULL);
    const gdouble f0    = _profile_h (ctx, rho, theta_seed, NULL);
    const gdouble fm    = _profile_h (ctx, rho_m, theta_seed, NULL);
    const gdouble grad  = (fp - fm) / (rho_p - rho_m);
    const gdouble hess  = (fp - 2.0 * f0 + fm) / gsl_pow_2 (0.5 * (rho_p - rho_m));

    if (hess >= 0.0) /* not concave: not a local max along this direction */
      break;

    rho -= grad / hess;
    rho  = CLAMP (rho, 1.0e-6, 1.0 - 1.0e-6);
  }

  /* Only used to resolve the theta matching the final rho; the profile's own
   * ln-integrand value is not what this function returns. */
  (void) _profile_h (ctx, rho, theta_seed, theta_out);

  return rho;
}

/**
 * nc_galaxy_shape_intrinsic_mode_find:
 * @apply_shear: convention-dispatched forward shear map
 * @apply_shear_inv: convention-dispatched inverse shear map (used only for
 * the initial guess)
 * @pop: a #NcGalaxyShapePop
 * @pop_data: its resolved #NcGalaxyShapePopData
 * @g: reduced shear
 * @eps_obs: observed ellipticity
 * @std_noise: noise standard deviation
 * @mode: (out): the found joint mode and its log-integrand Hessian
 *
 * Seeds the search at the naive noiseless intrinsic preimage
 * $\chi_{I,0}=f_g^{-1}(\epsilon_\mathrm{obs})$ (already a good starting point
 * for the noise-only theta direction, see nc_galaxy_shape_factor_quad.c), then
 * nests the theta profile inside a 1D Newton search over rho.
 */
void
nc_galaxy_shape_intrinsic_mode_find (complex double ( *apply_shear ) (complex double, complex double),
                                     complex double (*apply_shear_inv) (complex double, complex double),
                                     NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                     complex double g, complex double eps_obs, gdouble std_noise,
                                     NcGalaxyShapeIntrinsicMode *mode)
{
  ModeCtx ctx                      = {apply_shear, pop, pop_data, g, eps_obs, gsl_pow_2 (std_noise)};
  const complex double chi_i_naive = apply_shear_inv (g, eps_obs);
  const gdouble theta0             = carg (chi_i_naive);
  const gdouble rho0               = CLAMP (cabs (chi_i_naive), 1.0e-3, 1.0 - 1.0e-3);
  gdouble theta_star;
  gdouble x1, x2;
  const gdouble h = 1.0e-4;

  mode->rho     = _rho_hat (&ctx, rho0, theta0, &theta_star);
  mode->theta   = theta_star;
  mode->ln_peak = _ln_integrand (&ctx, mode->rho, mode->theta);

  x1 = mode->rho * cos (mode->theta);
  x2 = mode->rho * sin (mode->theta);

  {
    /* Full 2x2 Hessian of ln(integrand) in Cartesian (chi_I_1, chi_I_2),
     * central differences, no theta/rho decoupling assumed (unlike the mode
     * search itself, the curvature at the peak mixes both directions in
     * general). */
    const gdouble f00 = _ln_integrand (&ctx, hypot (x1, x2), atan2 (x2, x1));
    gdouble fpp, fmm, fpm, fmp, fp0, fm0, f0p, f0m;

    #define LN_AT(a, b) _ln_integrand (&ctx, hypot ((a), (b)), atan2 ((b), (a)))
    fp0 = LN_AT (x1 + h, x2);
    fm0 = LN_AT (x1 - h, x2);
    f0p = LN_AT (x1, x2 + h);
    f0m = LN_AT (x1, x2 - h);
    fpp = LN_AT (x1 + h, x2 + h);
    fmm = LN_AT (x1 - h, x2 - h);
    fpm = LN_AT (x1 + h, x2 - h);
    fmp = LN_AT (x1 - h, x2 + h);
    #undef LN_AT

    mode->hxx = (fp0 - 2.0 * f00 + fm0) / gsl_pow_2 (h);
    mode->hyy = (f0p - 2.0 * f00 + f0m) / gsl_pow_2 (h);
    mode->hxy = (fpp - fpm - fmp + fmm) / (4.0 * gsl_pow_2 (h));
  }
}

/**
 * nc_galaxy_shape_intrinsic_mode_laplace:
 * @mode: a #NcGalaxyShapeIntrinsicMode, as filled by nc_galaxy_shape_intrinsic_mode_find()
 *
 * Laplace approximation of $\int d^2\chi_I\,\exp(\ln\mathrm{integrand})$
 * around @mode: $2\pi/\sqrt{\det(-H)}\,\exp(\ln_\mathrm{peak})$.
 *
 * Returns: the Laplace estimate, or NAN if the mode is not a proper local
 * maximum (non-positive-definite $-H$), which can happen for a very broad,
 * nearly-flat, or boundary-hugging population where a single local Gaussian
 * is not a meaningful description of the peak.
 */
gdouble
nc_galaxy_shape_intrinsic_mode_laplace (const NcGalaxyShapeIntrinsicMode *mode)
{
  const gdouble det_neg_h   = mode->hxx * mode->hyy - gsl_pow_2 (mode->hxy);
  const gdouble trace_neg_h = -(mode->hxx + mode->hyy);

  if ((det_neg_h <= 0.0) || (trace_neg_h <= 0.0))
    return NAN;

  return (2.0 * M_PI / sqrt (det_neg_h)) * exp (mode->ln_peak);
}

