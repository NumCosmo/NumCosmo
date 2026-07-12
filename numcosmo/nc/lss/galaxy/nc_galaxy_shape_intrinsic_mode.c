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
 * NcGalaxyShapeIntrinsicMode:
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
#include "nc/lss/wl/nc_wl_ellipticity.h"
#include "ncm/algebra/ncm_poly_roots.h"

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

/* Full 2x2 Hessian of ln(integrand) in Cartesian (chi_I_1, chi_I_2) at
 * (rho,theta), central differences -- shared by every mode-finding entry
 * point regardless of how (rho,theta) itself was found, since this part
 * has no convention-specific closed form (yet). */
static void
_hessian_at (const ModeCtx *ctx, const gdouble rho, const gdouble theta, NcGalaxyShapeIntrinsicMode *mode)
{
  const gdouble h   = 1.0e-4;
  const gdouble x1  = rho * cos (theta);
  const gdouble x2  = rho * sin (theta);
  const gdouble f00 = _ln_integrand (ctx, hypot (x1, x2), atan2 (x2, x1));
  gdouble fpp, fmm, fpm, fmp, fp0, fm0, f0p, f0m;

  #define LN_AT(a, b) _ln_integrand (ctx, hypot ((a), (b)), atan2 ((b), (a)))
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

/* --- TRACE_DET closed-form theta-profile and rho-Newton -----------------
 *
 * For the TRACE_DET convention, f_g is the Blaschke/Mobius disk automorphism
 * nc_wl_ellipticity_apply_shear_trace_det_c(g,chi)=(chi+g)/(1+conj(g)chi).
 * As theta ranges over a full turn at fixed rho, f_g(rho e^{i theta}) traces
 * an EXACT circle in the eps_obs plane -- center c(rho), radius r_rho(rho)
 * below (derived by inverting the map and imposing |chi|=rho; verified
 * against a brute-force theta grid search to machine precision). So
 *   min_theta |eps_obs - f_g(rho e^{i theta})|^2 = (d(rho)-r_rho(rho))^2,
 *   d(rho) = |eps_obs - c(rho)|,
 * closed-form, no search over theta at all. Differentiating this profile
 * w.r.t. rho (also closed-form, verified against finite differences) turns
 * the whole outer Newton search into one with no apply_shear/cexp calls
 * anywhere in the loop -- apply_shear_inv is used only once at the end, to
 * recover theta_hat from the closest point on the circle. */

/* mu(rho), A(rho)=d(rho)^2, r_rho(rho) and their rho-derivatives, in terms
 * of gamma=|g|^2, beta=Re(eps_obs * conj(g)), chiO2=|eps_obs|^2 (fixed per
 * galaxy). CSE'd directly from the symbolic derivation, verified against
 * finite differences of A/r themselves. */
static void
_trace_det_A_r (const gdouble rho, const gdouble gamma, const gdouble beta, const gdouble chiO2,
                gdouble *A, gdouble *Ap, gdouble *App, gdouble *r_rho, gdouble *rp, gdouble *rpp)
{
  const gdouble x0  = rho * rho;
  const gdouble x1  = gamma * x0;
  const gdouble x2  = x1 - 1.0;
  const gdouble x3  = 1.0 / x2;
  const gdouble x4  = x0 - 1.0;
  const gdouble x5  = beta * x4;
  const gdouble x6  = x4 * x4;
  const gdouble x7  = x2 * x2;
  const gdouble x8  = 1.0 / x7;
  const gdouble x9  = x2 * x7;
  const gdouble x10 = 1.0 / x9;
  const gdouble x11 = gamma * gamma;
  const gdouble x12 = 4.0 * x0;
  const gdouble x13 = gamma - 1.0;
  const gdouble x14 = rho * x13;

  *A   = chiO2 + gamma * x6 * x8 - 2.0 * x3 * x5;
  *Ap  = -4.0 * rho * x10 * (beta * x7 - gamma * x2 * x4 * (beta + 1.0) + x11 * x6);
  *App = -4.0 * (beta * x9 - 6.0 * x11 * gamma * x0 * x6 - gamma * x7 * (beta * x12 + 3.0 * x0 + x5 - 1.0)
                 + x11 * x2 * (8.0 * x0 * x4 + x12 * x5 + x6)) / (x7 * x7);
  *r_rho = x14 * x3;
  *rp    = -x13 * x8 * (x1 + 1.0);
  *rpp   = 2.0 * gamma * x10 * x14 * (x1 + 3.0);
}

/* Q(rho) = (d(rho)-r_rho(rho))^2 and its two rho-derivatives, from A/Ap/App
 * and r_rho/rp/rpp above (d = sqrt(A)); both derivatives verified against
 * finite differences of Q itself. */
static void
_trace_det_Q (const gdouble rho, const gdouble gamma, const gdouble beta, const gdouble chiO2,
              gdouble *Q, gdouble *Qp, gdouble *Qpp)
{
  gdouble A, Ap, App, r_rho, rp, rpp;

  _trace_det_A_r (rho, gamma, beta, chiO2, &A, &Ap, &App, &r_rho, &rp, &rpp);

  {
    const gdouble d = sqrt (fmax (A, 0.0));

    *Q   = gsl_pow_2 (d - r_rho);
    *Qp  = Ap * (1.0 - r_rho / d) - 2.0 * rp * (d - r_rho);
    *Qpp = 2.0 * gsl_pow_2 (rp) - 2.0 * rp * Ap / d + r_rho * gsl_pow_2 (Ap) / (2.0 * gsl_pow_3 (d))
           + App * (d - r_rho) / d - 2.0 * rpp * (d - r_rho);
  }
}

/* ln(P_pop(x)/pi), x=rho^2 -- same normalization _ln_integrand uses for its
 * own p_pop term, so mode->ln_peak stays exactly consistent with the
 * existing finite-difference path. */
static gdouble
_ln_p_pop (const ModeCtx *ctx, const gdouble x)
{
  const gdouble p = nc_galaxy_shape_pop_eval_p (ctx->pop, ctx->pop_data, x) / M_PI;

  return (p > 0.0) && isfinite (p) ? log (p) : -G_MAXDOUBLE;
}

/* Outer 1D Newton search over rho, TRACE_DET closed form: h(rho) =
 * ln P_pop(rho^2) - Q(rho)/(2 sigma_n^2), h'/h'' both closed-form (Q part)
 * plus a finite-difference L'(x)/L''(x) for the population term (a plain
 * scalar special-function call, no apply_shear involved -- cheap either
 * way, and not worth a closed form of its own yet). */
static gdouble
_rho_hat_trace_det (const ModeCtx *ctx, const gdouble gamma, const gdouble beta, const gdouble chiO2, const gdouble rho0)
{
  gdouble rho = rho0;
  guint i;

  for (i = 0; i < 6; i++)
  {
    const gdouble x0      = rho * rho;
    const gdouble hx      = fmin (1.0e-5, 0.5 * x0);
    const gdouble Lp      = _ln_p_pop (ctx, x0 + hx);
    const gdouble L0      = _ln_p_pop (ctx, x0);
    const gdouble Lm      = _ln_p_pop (ctx, x0 - hx);
    const gdouble Lprime  = (Lp - Lm) / (2.0 * hx);
    const gdouble Lprime2 = (Lp - 2.0 * L0 + Lm) / gsl_pow_2 (hx);
    gdouble Q, Qp, Qpp, hp, hpp;

    _trace_det_Q (rho, gamma, beta, chiO2, &Q, &Qp, &Qpp);

    hp  = 2.0 * rho * Lprime - Qp / (2.0 * ctx->noise_var);
    hpp = 2.0 * Lprime + 4.0 * x0 * Lprime2 - Qpp / (2.0 * ctx->noise_var);

    if (hpp >= 0.0) /* not concave: not a local max along this direction */
      break;

    rho -= hp / hpp;
    rho  = CLAMP (rho, 1.0e-6, 1.0 - 1.0e-6);
  }

  return rho;
}

/**
 * nc_galaxy_shape_intrinsic_mode_find_trace_det:
 * @pop: a #NcGalaxyShapePop
 * @pop_data: its resolved #NcGalaxyShapePopData
 * @g: reduced shear
 * @eps_obs: observed ellipticity (TRACE_DET convention)
 * @std_noise: noise standard deviation
 * @mode: (out): the found joint mode and its log-integrand Hessian
 *
 * TRACE_DET-specific closed-form replacement for
 * nc_galaxy_shape_intrinsic_mode_find(): the theta-profile at any fixed rho,
 * and both derivatives of the profiled objective w.r.t. rho, are exact
 * closed forms (see the file-local comment above), so the whole rho-Newton
 * search runs with no apply_shear/apply_shear_inv calls in its loop at
 * all -- apply_shear_inv is called exactly once, at the end, to recover
 * theta_hat from the closest point on the circle. The final Hessian is
 * still the shared finite-difference _hessian_at() (not yet closed-form).
 */
void
nc_galaxy_shape_intrinsic_mode_find_trace_det (NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                               complex double g, complex double eps_obs, gdouble std_noise,
                                               NcGalaxyShapeIntrinsicMode *mode)
{
  ModeCtx ctx                      = {&nc_wl_ellipticity_apply_shear_trace_det_c, pop, pop_data, g, eps_obs, gsl_pow_2 (std_noise)};
  const gdouble gamma              = creal (g * conj (g));
  const gdouble beta               = creal (eps_obs * conj (g));
  const gdouble chiO2              = creal (eps_obs * conj (eps_obs));
  const complex double chi_i_naive = nc_wl_ellipticity_apply_shear_inv_trace_det_c (g, eps_obs);
  const gdouble rho0               = CLAMP (cabs (chi_i_naive), 1.0e-3, 1.0 - 1.0e-3);
  complex double c, w_star, chi_i_star;
  gdouble mu_rho, r_rho, d;

  mode->rho = _rho_hat_trace_det (&ctx, gamma, beta, chiO2, rho0);

  /* Recover theta_hat: closest point on the circle (center c, radius
   * r_rho) to eps_obs, then invert through the Blaschke map. */
  mu_rho = (1.0 - gsl_pow_2 (mode->rho)) / (1.0 - gamma * gsl_pow_2 (mode->rho));
  r_rho  = mode->rho * (1.0 - gamma) / (1.0 - gamma * gsl_pow_2 (mode->rho));
  c      = g * mu_rho;
  d      = fmax (cabs (eps_obs - c), 1.0e-300);
  w_star = c + r_rho * (eps_obs - c) / d;

  chi_i_star    = nc_wl_ellipticity_apply_shear_inv_trace_det_c (g, w_star);
  mode->theta   = carg (chi_i_star);
  mode->ln_peak = _ln_integrand (&ctx, mode->rho, mode->theta);

  _hessian_at (&ctx, mode->rho, mode->theta, mode);
}

/* --- TRACE closed-form theta-profile and rho-Newton ----------------------
 *
 * Unlike TRACE_DET, TRACE's forward shear map does not send circles to
 * circles (verified numerically: a least-squares circle fit to the image
 * of |chi_I|=rho has ~10-30% residuals, not noise), so there is no fully
 * closed-form theta-profile here. Two facts still make this closed-form,
 * just via a polynomial root instead of a single formula:
 *
 * 1. Rotation covariance: nc_wl_ellipticity_apply_shear_trace_c(g,chi)
 *    satisfies f(chi*e^{i a}, g*e^{i a}) = f(chi,g)*e^{i a} for any real a
 *    (verified numerically). So rotating g to be real and non-negative
 *    (a = -arg(g)) and applying the SAME rotation to eps_obs leaves the
 *    whole |eps_obs - f(chi_I,g)|^2 objective unchanged -- only theta
 *    itself needs unrotating at the very end. This gauge-fixing is exact,
 *    not an approximation.
 * 2. With g real, minimizing |eps_obs_rot - f(rho e^{i theta},g)|^2 over
 *    theta reduces (Weierstrass substitution t=tan(theta/4), clearing
 *    denominators) to finding the real roots of an explicit QUARTIC in t
 *    (not a sextic, as the fully general-complex-g case would give --
 *    exploiting the rotation is what buys the lower degree) -- verified
 *    against a brute-force theta grid search.
 *
 * Once theta_hat(rho) is found this way, Q(rho), Q'(rho), Q''(rho) use the
 * exact same envelope-theorem trick as TRACE_DET (see the file doc): since
 * theta_hat is by construction a stationary point of D2(rho,theta) in
 * theta, Q'(rho) is just the PARTIAL rho-derivative of D2 at theta_hat (no
 * need to differentiate through the root at all), and Q''(rho) is the
 * standard profile-curvature correction. These partials are taken directly
 * on the ORIGINAL (unrotated) g/eps_obs -- rotation is only needed to find
 * theta_hat cheaply, not for the derivatives themselves. Both facts
 * verified against finite differences of a precision-refined brute-force
 * theta_hat(rho). */

/* D2(rho,theta)=|eps_obs-f(rho e^{i theta};g)|^2 and its rho/theta partials,
 * general complex g -- CSE'd directly from the symbolic derivation. */
static void
_trace_D2_and_partials (const gdouble rho, const gdouble theta, const gdouble g_r, const gdouble g_i,
                        const gdouble cO_r, const gdouble cO_i,
                        gdouble *D2, gdouble *D2_rho, gdouble *D2_theta,
                        gdouble *D2_rr, gdouble *D2_rt, gdouble *D2_tt)
{
  const gdouble x0  = g_i * g_i;
  const gdouble x1  = g_r * g_r;
  const gdouble x2  = 2.0 * g_i;
  const gdouble x3  = sin (theta);
  const gdouble x4  = rho * x3;
  const gdouble x5  = x2 * x4;
  const gdouble x6  = 2.0 * g_r;
  const gdouble x7  = cos (theta);
  const gdouble x8  = rho * x7;
  const gdouble x9  = x6 * x8;
  const gdouble x10 = x0 + x1 + x5 + x9 + 1.0;
  const gdouble x11 = x10 * x10;
  const gdouble x12 = 1.0 / x11;
  const gdouble x13 = cO_r * cO_r;
  const gdouble x14 = g_i * x4 * x6 - x0 * x8 + x1 * x8 + x6 + x8;
  const gdouble x15 = cO_r * x10;
  const gdouble x16 = 2.0 * x14;
  const gdouble x17 = cO_i * x0 + cO_i * x1 + cO_i * x5 + cO_i * x9 + cO_i - g_i * x9 - x0 * x4 + x1 * x4 - x2 - x4;
  const gdouble x18 = x11 * x13 + x14 * x14 - x15 * x16 + x17 * x17;
  const gdouble x19 = x12 * x18;
  const gdouble x20 = g_i * x3;
  const gdouble x21 = g_r * x7 + x20;
  const gdouble x22 = 2.0 * x21;
  const gdouble x23 = 1.0 / x10;
  const gdouble x24 = x18 * x23;
  const gdouble x25 = x10 * x13;
  const gdouble x26 = x22 * x25;
  const gdouble x27 = cO_r * x16;
  const gdouble x28 = x21 * x27;
  const gdouble x29 = -x0 * x7 + x1 * x7 + x20 * x6;
  const gdouble x30 = x29 + x7;
  const gdouble x31 = x15 * x30;
  const gdouble x32 = x14 * x30;
  const gdouble x33 = -x3;
  const gdouble x34 = x1 * x3;
  const gdouble x35 = cO_i * x3;
  const gdouble x36 = cO_i * x7;
  const gdouble x37 = g_i * x7;
  const gdouble x38 = x37 * x6;
  const gdouble x39 = x0 * x3;
  const gdouble x40 = x2 * x35 + x33 + x34 + x36 * x6 - x38 - x39;
  const gdouble x41 = x17 * x40;
  const gdouble x42 = x26 - x28 - x31 + x32 + x41;
  const gdouble x43 = 2.0 * x12;
  const gdouble x44 = -g_r * x3 + x37;
  const gdouble x45 = 2.0 * x44;
  const gdouble x46 = x33 - x34 + x38 + x39;
  const gdouble x47 = x2 * x36 + x29 - x35 * x6 - x7;
  const gdouble x48 = x14 * x46 - x15 * x46 + x17 * x47 + x25 * x45 - x27 * x44;
  const gdouble x49 = -x24 * x45 + x48;
  const gdouble x50 = rho * x43;
  const gdouble x51 = x21 * x21;
  const gdouble x52 = 4.0 * x13;
  const gdouble x53 = 4.0 * x21;
  const gdouble x54 = 12.0 * x19;
  const gdouble x55 = x23 * x42;
  const gdouble x56 = rho * x30;
  const gdouble x57 = rho * x44;
  const gdouble x58 = 4.0 * x57;
  const gdouble x59 = x23 * x48;
  const gdouble x60 = cO_r * x46;
  const gdouble x61 = x21 * x57;
  const gdouble x62 = rho * x44 * x44;

  *D2       = x19;
  *D2_rho   = x43 * (-x22 * x24 + x42);
  *D2_theta = x49 * x50;
  *D2_rr    = x43 * (-cO_r * x30 * x53 - 8.0 * x21 * x55 + x30 * x30 + x40 * x40 + x51 * x52 + x51 * x54);
  *D2_rt    = x43 * (-cO_r * x45 * x56 - rho * x22 * x60 + rho * x40 * x47 - rho * x53 * x59 + x46 * x56 + x49 + x52 * x61 + x54 * x61 - x55 * x58);
  *D2_tt    = x50 * (rho * x46 * x46 + rho * x47 * x47 + 2.0 * x24 * (x21 + 6.0 * x23 * x62) - x26 + x28 + x31 - x32 - x41 + x52 * x62 - 8.0 * x57 * x59 - x58 * x60);
}

/* D2 value only (no partials), for cheap evaluation at each candidate root
 * inside the theta-profile search below -- computing the full Hessian-
 * feeding partials for every candidate (up to 4 roots x 3 branches = 12
 * evaluations per rho-trial) was needlessly expensive when only the value
 * is needed to pick the best one; the full _trace_D2_and_partials() above
 * is only called once, on the winning candidate. */
static gdouble
_trace_D2_only (const gdouble rho, const gdouble theta, const gdouble g_r, const gdouble g_i,
                const gdouble cO_r, const gdouble cO_i)
{
  const gdouble x0  = g_i * g_i;
  const gdouble x1  = g_r * g_r;
  const gdouble x2  = 2.0 * g_i;
  const gdouble x3  = rho * sin (theta);
  const gdouble x4  = x2 * x3;
  const gdouble x5  = 2.0 * g_r;
  const gdouble x6  = rho * cos (theta);
  const gdouble x7  = x5 * x6;
  const gdouble x8  = x0 + x1 + x4 + x7 + 1.0;
  const gdouble x9  = x8 * x8;
  const gdouble x10 = g_i * x3 * x5 - x0 * x6 + x1 * x6 + x5 + x6;
  const gdouble x11 = cO_i * x0 + cO_i * x1 + cO_i * x4 + cO_i * x7 + cO_i - g_i * x7 - x0 * x3 + x1 * x3 - x2 - x3;

  return (gsl_pow_2 (cO_r) * x9 - 2.0 * cO_r * x10 * x8 + gsl_pow_2 (x10) + gsl_pow_2 (x11)) / x9;
}

/* Quartic-in-t (t=tan(theta_rot/4)) coefficients for the stationarity
 * condition of D2 in the rotated (g real >= 0) frame, ascending order
 * (a[0]..a[4], matching gsl_poly_complex_solve's convention) -- CSE'd
 * directly from the symbolic derivation. chiO/lam here are the ROTATED
 * eps_obs's magnitude/phase (magnitude is rotation-invariant; lam=2*arg of
 * the rotated eps_obs). */
static void
_trace_quartic_coeffs (const gdouble rho, const gdouble lam, const gdouble g, const gdouble chiO,
                       gdouble a[5])
{
  const gdouble x0  = gsl_pow_5 (g);
  const gdouble x1  = 4.0 * rho;
  const gdouble x2  = g * x1;
  const gdouble x3  = gsl_pow_4 (g);
  const gdouble x4  = gsl_pow_6 (g);
  const gdouble x5  = g * g;
  const gdouble x6  = rho * rho;
  const gdouble x7  = 4.0 * x6;
  const gdouble x8  = x3 * x7 + x3 + x4 - x5 * x7 - x5 - 1.0;
  const gdouble x9  = 0.5 * lam;
  const gdouble x10 = chiO * x1 * sin (x9);
  const gdouble x11 = gsl_pow_3 (g);
  const gdouble x12 = 2.0 * g;
  const gdouble x13 = 2.0 * x0;
  const gdouble x14 = cos (x9);
  const gdouble x15 = 4.0 * x11;
  const gdouble x16 = x15 * x6;
  const gdouble x17 = chiO * x14;
  const gdouble x18 = x17 * x5;
  const gdouble x19 = x17 * x3;
  const gdouble x20 = rho * x17;
  const gdouble x21 = x12 * x20 + x13 * x20 - x15 * x20;
  const gdouble x22 = 8.0 * rho;

  a[4] = x10 * (4.0 * rho * x0 - x2 - x8);
  a[3] = x22 * (chiO * x14 * x4 + chiO * x14 + 2.0 * g * x6 + 2.0 * x0 * x6 + 4.0 * x11 - x12 - x13 - x16 - x18 - x19 - x21);
  a[2] = 0.0;
  a[1] = x22 * (x12 * x6 - x12 + x13 * x6 - x13 + x15 - x16 + x17 * x4 + x17 - x18 - x19 + x21);
  a[0] = x10 * (x0 * x1 - x2 + x8);
}

/* Real-root-finding itself (bracket-and-bisect on the derivative chain) has
 * no physics content and lives in ncm_poly_roots_real_quartic_or_lower()
 * (numcosmo/ncm/algebra/ncm_poly_roots.c), with its own independent tests.
 * See that file's docs for the derivation; faster than
 * gsl_poly_complex_solve()'s general eigenvalue approach. */

/* Solve the theta-profile stationarity quartic and return the best
 * theta_hat (in the ROTATED frame, i.e. still needs -alpha applied by the
 * caller) via D2 evaluated at each real root. */
static gdouble
_trace_theta_hat_rot (const gdouble rho, const gdouble lam,
                      const gdouble g_rot, const gdouble chiO_mag, const gdouble chiOr_r, const gdouble chiOr_i)
{
  gdouble a[5];
  gdouble best_val = G_MAXDOUBLE, best_theta = 0.0;
  gint k;

  _trace_quartic_coeffs (rho, lam, g_rot, chiO_mag, a);

  {
    gdouble z[4];
    gint n_roots = ncm_poly_roots_real_quartic_or_lower (a, z);

    for (k = 0; k < n_roots; k++)
    {
      const gdouble zr = z[k];
      gint branch;

      for (branch = -1; branch <= 1; branch++)
      {
        /* t=tan(theta_B/4) where theta_B is the Weierstrass-substitution
         * variable from the symbolic derivation; the actual angle is
         * theta_B/2, NOT theta_B itself -- i.e. 2*atan(t)+branch*pi, not
         * 4*atan(t)+branch*2pi. See docs/theory/wl_shape_factor_history.md
         * for the off-by-factor-of-2 mistake this guards against. */
        const gdouble theta_rot = 2.0 * atan (zr) + branch * M_PI;
        const gdouble D2        = _trace_D2_only (rho, theta_rot, g_rot, 0.0, chiOr_r, chiOr_i);

        if (D2 < best_val)
        {
          best_val   = D2;
          best_theta = theta_rot;
        }
      }
    }
  }

  return best_theta;
}

/* Outer 1D Newton search over rho, TRACE closed form: same structure as
 * _rho_hat_trace_det, but theta_hat(rho) comes from the quartic above
 * (in a rotated frame) and Q/Q'/Q'' come from the envelope theorem applied
 * to _trace_D2_and_partials on the ORIGINAL (unrotated) g/eps_obs. */
static gdouble
_rho_hat_trace (const ModeCtx *ctx,
                const gdouble alpha, const gdouble g_rot, const gdouble chiOr_r, const gdouble chiOr_i,
                const gdouble lam, const gdouble g_r, const gdouble g_i, const gdouble cO_r, const gdouble cO_i,
                const gdouble rho0, gdouble *theta_out)
{
  const gdouble chiO_mag = hypot (chiOr_r, chiOr_i); /* rotation-invariant, same as cabs(eps_obs) */
  gdouble rho            = rho0;
  gdouble theta_true     = 0.0;
  guint i;

  for (i = 0; i < 6; i++)
  {
    const gdouble x0        = rho * rho;
    const gdouble hx        = fmin (1.0e-5, 0.5 * x0);
    const gdouble Lp        = _ln_p_pop (ctx, x0 + hx);
    const gdouble L0        = _ln_p_pop (ctx, x0);
    const gdouble Lm        = _ln_p_pop (ctx, x0 - hx);
    const gdouble Lprime    = (Lp - Lm) / (2.0 * hx);
    const gdouble Lprime2   = (Lp - 2.0 * L0 + Lm) / gsl_pow_2 (hx);
    const gdouble theta_rot = _trace_theta_hat_rot (rho, lam, g_rot, chiO_mag, chiOr_r, chiOr_i);
    gdouble D2, D2_rho, D2_theta, D2_rr, D2_rt, D2_tt, Qp, Qpp, hp, hpp;

    theta_true = theta_rot - alpha;

    _trace_D2_and_partials (rho, theta_true, g_r, g_i, cO_r, cO_i,
                            &D2, &D2_rho, &D2_theta, &D2_rr, &D2_rt, &D2_tt);

    Qp  = D2_rho;
    Qpp = D2_rr - gsl_pow_2 (D2_rt) / D2_tt;

    hp  = 2.0 * rho * Lprime - Qp / (2.0 * ctx->noise_var);
    hpp = 2.0 * Lprime + 4.0 * x0 * Lprime2 - Qpp / (2.0 * ctx->noise_var);

    if (hpp >= 0.0)
      break;

    rho -= hp / hpp;
    rho  = CLAMP (rho, 1.0e-6, 1.0 - 1.0e-6);
  }

  /* theta_true above matches whatever rho the loop last USED, not the rho
   * it produced on its final update (or the pre-loop rho0 if the very
   * first concavity check already failed) -- re-evaluate theta_hat once
   * more at the actual final rho, exactly like _rho_hat's own trailing
   * _profile_h() call does for the finite-difference path. */
  theta_true = _trace_theta_hat_rot (rho, lam, g_rot, chiO_mag, chiOr_r, chiOr_i) - alpha;
  *theta_out = theta_true;

  return rho;
}

/**
 * nc_galaxy_shape_intrinsic_mode_find_trace:
 * @pop: a #NcGalaxyShapePop
 * @pop_data: its resolved #NcGalaxyShapePopData
 * @g: reduced shear
 * @eps_obs: observed ellipticity (TRACE convention)
 * @std_noise: noise standard deviation
 * @mode: (out): the found joint mode and its log-integrand Hessian
 *
 * TRACE-specific closed-form replacement for
 * nc_galaxy_shape_intrinsic_mode_find(): theta_hat(rho) comes from a quartic
 * root (in a rotation-gauge-fixed frame, see the file-local comment above)
 * instead of a finite-difference Newton search, and both rho-derivatives of
 * the profiled objective are closed-form via the envelope theorem -- no
 * apply_shear/apply_shear_inv calls anywhere in the rho-Newton loop.
 * apply_shear_inv is not needed at all here (the quartic already gives
 * theta_hat directly, unlike TRACE_DET's circle-inversion step). The final
 * Hessian is still the shared finite-difference _hessian_at().
 */
void
nc_galaxy_shape_intrinsic_mode_find_trace (NcGalaxyShapePop *pop, NcGalaxyShapePopData *pop_data,
                                           complex double g, complex double eps_obs, gdouble std_noise,
                                           NcGalaxyShapeIntrinsicMode *mode)
{
  ModeCtx ctx                      = {&nc_wl_ellipticity_apply_shear_trace_c, pop, pop_data, g, eps_obs, gsl_pow_2 (std_noise)};
  const gdouble alpha              = -carg (g);
  const gdouble g_rot              = cabs (g);
  const complex double eps_obs_rot = eps_obs * cexp (I * alpha);
  const gdouble lam                = 2.0 * carg (eps_obs_rot);
  const complex double chi_i_naive = nc_wl_ellipticity_apply_shear_inv_trace_c (g, eps_obs);
  const gdouble rho0               = CLAMP (cabs (chi_i_naive), 1.0e-3, 1.0 - 1.0e-3);
  gdouble theta_star;

  mode->rho = _rho_hat_trace (&ctx, alpha, g_rot, creal (eps_obs_rot), cimag (eps_obs_rot), lam,
                              creal (g), cimag (g), creal (eps_obs), cimag (eps_obs), rho0, &theta_star);
  mode->theta = theta_star;

  mode->ln_peak = _ln_integrand (&ctx, mode->rho, mode->theta);

  _hessian_at (&ctx, mode->rho, mode->theta, mode);
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

  mode->rho     = _rho_hat (&ctx, rho0, theta0, &theta_star);
  mode->theta   = theta_star;
  mode->ln_peak = _ln_integrand (&ctx, mode->rho, mode->theta);

  _hessian_at (&ctx, mode->rho, mode->theta, mode);
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

