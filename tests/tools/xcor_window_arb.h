/***************************************************************************
 *            xcor_window_arb.h
 *
 *  Fri August 28 12:00:00 2026
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

/*
 * The analytic xcor windows in Arb ball arithmetic, shared by the reference
 * generators in this directory.
 *
 *   I_ell(k) = int W(chi) g(chi, k) j_ell(k chi) dchi
 *
 * with chi in Mpc, W normalized so that int W dchi = 1 over its truncated
 * support, and g the optional scale-dependent growth (1 unless kdep is on).
 *
 * Header-only and all-static: these tools are standalone programs compiled one
 * at a time, so this is the whole build integration.
 *
 * The wavenumber is carried as an `acb_t`, not a double. That is what lets an
 * outer integrator evaluate I_ell on a *ball* of k, which nc_xcor_kquad_arb.c
 * needs and nc_xcor_kernel_analytic_arb.c does not -- the latter simply sets a
 * point ball.
 */

#ifndef _XCOR_WINDOW_ARB_H
#define _XCOR_WINDOW_ARB_H

#if defined (__has_include )
#if __has_include (<flint/acb.h>)
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#include <flint/acb_calc.h>
#else
#include <acb.h>
#include <acb_hypgeom.h>
#include <acb_calc.h>
#endif
#else
#include <acb.h>
#include <acb_hypgeom.h>
#include <acb_calc.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_BUMPS 8
#define MAX_PANELS 16

/* Inner panels: how many oscillation periods of j_ell (k chi) each spans, and
 * the ceiling on how many that may produce. */
#define WINDOW_PANEL_PERIODS 4.0
#define MAX_CHI_PANELS 8192

typedef enum
{
  SHAPE_GAUSS,
  SHAPE_TOPHAT,
  SHAPE_TOPHAT_SMOOTH,
  SHAPE_STUDENT_T,
  SHAPE_POWER_EXP,
  SHAPE_LENSING,
  SHAPE_MULTI,
  SHAPE_LEN
} Shape;

static const char *shape_names[SHAPE_LEN] = {
  "gauss", "tophat", "tophat_smooth", "student_t", "power_exp", "lensing", "multi"
};

typedef struct
{
  Shape shape;
  long ell;
  acb_t k;

  /* gauss / student_t / multi */
  double chi_mean, chi_sigma, n_sigma;
  /* tophat / tophat_smooth / power_exp */
  double chi_lower, chi_upper;
  /* student_t */
  double chi_scale, nu, n_scale;
  /* power_exp */
  double alpha, beta;
  /* lensing */
  double chi_source_lower, chi_source_upper;

  /* multi. Bumps whose n_sigma supports meet form one group; the window is
   * the sum of a group's bumps over the group's whole stretch and zero in the
   * gaps between groups, as nc_xcor_kernel_analytic_multi.c defines it. */
  int n_bumps;
  double mu[MAX_BUMPS], sg[MAX_BUMPS], wt[MAX_BUMPS];
  int n_groups;
  double grp_lo[MAX_BUMPS], grp_hi[MAX_BUMPS];

  /* Scale-dependent growth, off unless kdep_on. Matches
   * _nc_xcor_kernel_radial_kdep_growth_eval(). */
  int kdep_on;
  double kdep_amplitude, kdep_k_transition, kdep_chi_ref;

  /* support, filled by shape_support () */
  double chi_min, chi_max;
  /* internal breakpoints where the shape is not analytic */
  int n_break;
  double breaks[MAX_PANELS];

  /* set while integrating: 0 = bare window, 1 = window times j_ell */
  int with_bessel;
} Par;

/* Index of the multi group holding chi, or -1 in a gap. */
static int
g_of (const Par *p, double c)
{
  int g;

  for (g = 0; g < p->n_groups; g++)
    if ((c >= p->grp_lo[g]) && (c <= p->grp_hi[g]))
      return g;

  return -1;
}

static void
par_init (Par *p)
{
  memset (p, 0, sizeof (Par));
  acb_init (p->k);
}

static void
par_clear (Par *p)
{
  acb_clear (p->k);
}

/*
 * The unnormalized window, on a complex ball.
 *
 * Returns non-analytic (an indeterminate ball) when @order is 1 and the input
 * ball may cross a branch cut, which is how acb_calc_integrate is told to
 * subdivide rather than trust the value. Only `lensing` and `power_exp` have a
 * cut at all -- both take a log or a real power of chi, analytic for Re chi > 0
 * -- so the guard is a positivity test on the real part.
 */
static int
window_u (acb_t out, const acb_t chi, Par *p, slong order, slong prec)
{
  acb_t t, u, v, w;
  int i;

  acb_init (t);
  acb_init (u);
  acb_init (v);
  acb_init (w);

  /* Analyticity guard for the shapes that take log/pow of chi. */
  if ((order > 0) &&
      ((p->shape == SHAPE_LENSING) || (p->shape == SHAPE_POWER_EXP)))
  {
    arb_t re;

    arb_init (re);
    acb_get_real (re, chi);

    if (!arb_is_positive (re))
    {
      acb_indeterminate (out);
      arb_clear (re);
      acb_clear (t);
      acb_clear (u);
      acb_clear (v);
      acb_clear (w);

      return 0;
    }

    arb_clear (re);
  }

  switch (p->shape)
  {
    case SHAPE_GAUSS:
      acb_set_d (t, p->chi_mean);
      acb_sub (u, chi, t, prec);
      acb_set_d (t, p->chi_sigma);
      acb_div (u, u, t, prec);
      acb_sqr (u, u, prec);
      acb_div_si (u, u, -2, prec);
      acb_exp (out, u, prec);
      break;

    case SHAPE_TOPHAT:
      acb_one (out);
      break;

    case SHAPE_TOPHAT_SMOOTH:
      /* erf ((chi_upper - chi)/(sqrt2 sigma)) - erf ((chi_lower - chi)/(sqrt2 sigma)) */
      acb_set_d (w, M_SQRT2 * p->chi_sigma);

      acb_set_d (t, p->chi_upper);
      acb_sub (u, t, chi, prec);
      acb_div (u, u, w, prec);
      acb_hypgeom_erf (u, u, prec);

      acb_set_d (t, p->chi_lower);
      acb_sub (v, t, chi, prec);
      acb_div (v, v, w, prec);
      acb_hypgeom_erf (v, v, prec);

      acb_sub (out, u, v, prec);
      break;

    case SHAPE_STUDENT_T:
      /* (1 + t^2/nu)^(-(nu+1)/2),  t = (chi - mean)/scale */
      acb_set_d (t, p->chi_mean);
      acb_sub (u, chi, t, prec);
      acb_set_d (t, p->chi_scale);
      acb_div (u, u, t, prec);
      acb_sqr (u, u, prec);
      acb_set_d (t, p->nu);
      acb_div (u, u, t, prec);
      acb_add_si (u, u, 1, prec);
      acb_set_d (v, -0.5 * (p->nu + 1.0));
      acb_pow (out, u, v, prec);
      break;

    case SHAPE_POWER_EXP:
      /* x^alpha exp (-x^beta),  x = chi / chi_scale */
      acb_set_d (t, p->chi_scale);
      acb_div (u, chi, t, prec);
      acb_set_d (v, p->alpha);
      acb_pow (w, u, v, prec);
      acb_set_d (v, p->beta);
      acb_pow (u, u, v, prec);
      acb_neg (u, u);
      acb_exp (u, u, prec);
      acb_mul (out, w, u, prec);
      break;

    case SHAPE_LENSING:
    {
      /* chi ((b - m) - chi log (b/m)) / (b - a), with m = max (chi, a).
       * The two branches meet at chi = a, where the shape has a kink; the
       * caller panels on it, so each call sees only one branch. */
      const double a = p->chi_source_lower;
      const double b = p->chi_source_upper;
      arb_t re, av;
      int upper;

      arb_init (re);
      arb_init (av);
      acb_get_real (re, chi);
      arb_set_d (av, a);
      upper = arb_gt (re, av);
      arb_clear (re);
      arb_clear (av);

      if (upper)
        acb_set (v, chi);  /* m = chi */
      else
        acb_set_d (v, a);  /* m = a   */

      acb_set_d (t, b);
      acb_div (u, t, v, prec); /* b / m */
      acb_log (u, u, prec);
      acb_mul (u, u, chi, prec); /* chi log (b/m) */

      acb_sub (w, t, v, prec); /* b - m */
      acb_sub (w, w, u, prec);
      acb_mul (w, w, chi, prec);

      acb_set_d (t, b - a);
      acb_div (out, w, t, prec);
      break;
    }

    case SHAPE_MULTI:
    {
      /* Zero in a gap between groups; inside a group, every bump of that group
       * across the group's whole stretch. The caller places panel edges at the
       * gap edges, so any one call lies wholly on one side of them, and the
       * midpoint of chi is enough to select the group. */
      double c;
      int g;

      /* Through a local copy: an interior pointer into @chi makes GCC 16 size
       * the ball as its midpoint alone and warn on every later read of it. */
      acb_set (v, chi);
      c = arf_get_d (arb_midref (acb_realref (v)), ARF_RND_NEAR);
      g = g_of (p, c);

      acb_zero (out);

      for (i = 0; (g >= 0) && (i < p->n_bumps); i++)
      {
        const double lo_i = fmax (0.0, p->mu[i] - p->n_sigma * p->sg[i]);

        if ((lo_i < p->grp_lo[g]) || (lo_i > p->grp_hi[g]))
          continue;

        acb_set_d (t, p->mu[i]);
        acb_sub (u, chi, t, prec);
        acb_set_d (t, p->sg[i]);
        acb_div (u, u, t, prec);
        acb_sqr (u, u, prec);
        acb_div_si (u, u, -2, prec);
        acb_exp (u, u, prec);
        acb_set_d (t, p->wt[i]);
        acb_mul (u, u, t, prec);
        acb_add (out, out, u, prec);
      }

      break;
    }

    default:
      fprintf (stderr, "window_u: bad shape\n");
      exit (1);
  }

  acb_clear (t);
  acb_clear (u);
  acb_clear (v);
  acb_clear (w);

  return 0;
}

/*
 * exp (-alpha * sigma * (chi_ref - chi) / chi_ref), sigma = x^2/(1 + x^2),
 * x = k / k_transition. Multiplies the radial integrand, not the window: the
 * library normalizes W alone and applies the growth alongside sqrt (P).
 */
static void
kdep_growth (acb_t out, const acb_t chi, const Par *p, slong prec)
{
  acb_t x, s, t;

  if (!p->kdep_on)
  {
    acb_one (out);

    return;
  }

  acb_init (x);
  acb_init (s);
  acb_init (t);

  acb_set_d (t, p->kdep_k_transition);
  acb_div (x, p->k, t, prec);
  acb_sqr (x, x, prec); /* x^2 */
  acb_add_si (s, x, 1, prec);
  acb_div (s, x, s, prec); /* sigma */

  acb_set_d (t, p->kdep_chi_ref);
  acb_sub (x, t, chi, prec);
  acb_div (x, x, t, prec); /* (chi_ref - chi)/chi_ref */

  acb_mul (x, x, s, prec);
  acb_set_d (t, -p->kdep_amplitude);
  acb_mul (x, x, t, prec);
  acb_exp (out, x, prec);

  acb_clear (x);
  acb_clear (s);
  acb_clear (t);
}

/*
 * j_ell (z), by whichever form is well conditioned at this z.
 *
 * Near the origin: the entire form
 *
 *   j_ell (z) = sqrt (pi) / 2^(ell+1) * z^ell * 0F1~ (; ell + 3/2; -z^2/4)
 *
 * The textbook sqrt (pi / 2z) J_{ell+1/2} (z) has a removable singularity at
 * z = 0, and a window whose support reaches the observer puts z = 0 inside the
 * domain, where the ball arithmetic sees 1/0 and subdivides without end.
 *
 * Away from it: the textbook form. The 0F1 series needs about z^2/4 terms and
 * loses about that many bits to cancellation, which is affordable at z ~ 100
 * and is not a method at all beyond it -- the C_ell integrand reaches
 * z = k chi ~ 4400 at ell = 200, where it would want 5 million bits. Arb
 * evaluates J_nu by asymptotic expansion out there, at no precision penalty.
 *
 * The switch is on a certified lower bound of |z|, so a ball that straddles
 * the threshold takes the entire form and stays correct either way.
 */
static void
sph_bessel (acb_t out, const acb_t z, long ell, slong prec)
{
  acb_t nu, J, t;
  arb_t az;
  arf_t lb;
  int far;

  arb_init (az);
  arf_init (lb);
  acb_abs (az, z, prec);
  arb_get_lbound_arf (lb, az, prec);
  far = arf_cmp_d (lb, 4.0) > 0;
  arb_clear (az);
  arf_clear (lb);

  acb_init (nu);
  acb_init (J);
  acb_init (t);

  if (far)
  {
    /* sqrt (pi / (2 z)) J_{ell + 1/2} (z) */
    acb_set_si (nu, 2 * ell + 1);
    acb_div_si (nu, nu, 2, prec);
    acb_hypgeom_bessel_j (J, nu, z, prec);
    acb_const_pi (t, prec);
    acb_div (t, t, z, prec);
    acb_div_si (t, t, 2, prec);
    acb_sqrt (t, t, prec);
    acb_mul (out, J, t, prec);
  }
  else
  {
    acb_set_si (nu, 2 * ell + 3);
    acb_div_si (nu, nu, 2, prec); /* ell + 3/2 */
    acb_sqr (J, z, prec);
    acb_div_si (J, J, -4, prec); /* -z^2/4 */
    acb_hypgeom_0f1 (J, nu, J, 1, prec);
    acb_pow_ui (t, z, (ulong) ell, prec);
    acb_mul (J, J, t, prec);
    acb_const_pi (t, prec);
    acb_sqrt (t, t, prec);
    acb_mul (J, J, t, prec);
    acb_mul_2exp_si (J, J, -(ell + 1));
    acb_set (out, J);
  }

  acb_clear (nu);
  acb_clear (J);
  acb_clear (t);
}

/* The integrand handed to acb_calc_integrate: W_u, optionally times g j_ell. */
static int
window_integrand (acb_ptr out, const acb_t chi, void *param, slong order, slong prec)
{
  Par *p = (Par *) param;
  acb_t W, z, J, g;

  if (order > 1)
    return 0;

  acb_init (W);
  window_u (W, chi, p, order, prec);

  if (!p->with_bessel)
  {
    acb_set (out, W);
    acb_clear (W);

    return 0;
  }

  acb_init (z);
  acb_init (J);
  acb_init (g);

  acb_mul (z, chi, p->k, prec);
  sph_bessel (J, z, p->ell, prec);
  kdep_growth (g, chi, p, prec);

  acb_mul (out, W, J, prec);
  acb_mul (out, out, g, prec);

  acb_clear (W);
  acb_clear (z);
  acb_clear (J);
  acb_clear (g);

  return 0;
}

/*
 * Cap on the inner integrator's work, 0 for Arb's default.
 *
 * Only a nested caller sets this. When k is a wide ball the integrand is
 * uncertain no matter how far chi is subdivided -- the uncertainty lives in k
 * -- so an uncapped run spends its whole budget and returns a wide ball
 * anyway. Capped, it returns that ball at once and the *outer* integrator
 * subdivides k, which is the axis that helps.
 */
static slong window_eval_limit = 0;

/* Sum of the certified integral over the shape's panels, at one precision. */
static void
integrate_panels (acb_t res, Par *p, slong prec)
{
  double edges[MAX_CHI_PANELS];
  acb_calc_integrate_opt_t opt;
  acb_t A, B, part;
  mag_t tol;
  int n = 0, i;

  edges[n++] = p->chi_min;

  for (i = 0; i < p->n_break; i++)
    if ((p->breaks[i] > p->chi_min) && (p->breaks[i] < p->chi_max))
      edges[n++] = p->breaks[i];

  edges[n++] = p->chi_max;

  /*
   * Subdivide on the oscillation scale of j_ell (k chi), which is 2 pi / k in
   * chi. Handed a single interval spanning hundreds of periods, the integrator
   * finds them by bisection and the cost of the whole nested integral grows as
   * the square of the phase k chi_max -- 3e4 for a Gaussian at ell = 2 against
   * 2.4e7 for a broad lensing window at ell = 10, which is the difference
   * between five seconds and not finishing.
   *
   * Only the width matters, so this is a refinement of the edge list above and
   * leaves the shape's own breakpoints where they are.
   */
  if (p->with_bessel)
  {
    arb_t ak;
    arf_t kub;
    double kmax;

    arb_init (ak);
    arf_init (kub);
    acb_abs (ak, p->k, prec);
    arb_get_ubound_arf (kub, ak, prec);
    kmax = arf_get_d (kub, ARF_RND_UP);
    arb_clear (ak);
    arf_clear (kub);

    if ((kmax > 0.0) && isfinite (kmax))
    {
      const double width = WINDOW_PANEL_PERIODS * 2.0 * M_PI / kmax;
      double refined[MAX_CHI_PANELS];
      int m = 0, j;

      for (j = 0; j + 1 < n && m + 2 < MAX_CHI_PANELS; j++)
      {
        double x = edges[j];

        refined[m++] = x;

        while ((edges[j + 1] - x > width) && (m + 2 < MAX_CHI_PANELS))
        {
          x           += width;
          refined[m++] = x;
        }
      }

      refined[m++] = edges[n - 1];

      if (m <= MAX_CHI_PANELS)
      {
        for (j = 0; j < m; j++)
          edges[j] = refined[j];

        n = m;
      }
    }
  }

  acb_calc_integrate_opt_init (opt);
  opt->eval_limit = window_eval_limit;
  acb_init (A);
  acb_init (B);
  acb_init (part);
  mag_init (tol);
  mag_set_ui_2exp_si (tol, 1, -prec);

  acb_zero (res);

  for (i = 0; i + 1 < n; i++)
  {
    acb_set_d (A, edges[i]);
    acb_set_d (B, edges[i + 1]);
    acb_calc_integrate (part, window_integrand, p, A, B, prec, tol, opt, prec);

    /* A capped run leaves an indeterminate ball, and one infinity poisons
     * every enclosing cell of a nested caller. Fall back to the trivial
     * enclosure (B - A) f([A, B]): crude, but finite, which is what lets the
     * caller subdivide instead. */
    if ((window_eval_limit > 0) && !acb_is_finite (part))
    {
      acb_t span, f;

      acb_init (span);
      acb_init (f);
      acb_union (span, A, B, prec);
      window_integrand (f, span, p, 0, prec);
      acb_sub (span, B, A, prec);
      acb_mul (part, f, span, prec);
      acb_clear (span);
      acb_clear (f);
    }

    acb_add (res, res, part, prec);
  }

  acb_clear (A);
  acb_clear (B);
  acb_clear (part);
  mag_clear (tol);
}

/* Recompute at doubling precision until the relative radius clears @target. */
static slong
certified (acb_t res, Par *p, double target)
{
  slong prec;

  for (prec = 128; prec <= 8192; prec *= 2)
  {
    double r, m;

    integrate_panels (res, p, prec);

    if (acb_is_zero (res))
      return prec;

    r = mag_get_d (arb_radref (acb_realref (res)));
    m = arf_get_d (arb_midref (acb_realref (res)), ARF_RND_NEAR);

    if ((m != 0.0) && (r / fabs (m) < target))
      return prec;
  }

  fprintf (stderr, "certified: did not reach %g (shape %s ell %ld)\n",
           target, shape_names[p->shape], p->ell);
  exit (1);
}

/* Support and internal breakpoints, matching the library's constructors. */
static void
shape_support (Par *p)
{
  int i;

  p->n_break = 0;

  switch (p->shape)
  {
    case SHAPE_GAUSS:
      p->chi_min = fmax (0.0, p->chi_mean - p->n_sigma * p->chi_sigma);
      p->chi_max = p->chi_mean + p->n_sigma * p->chi_sigma;
      break;

    case SHAPE_TOPHAT:
      p->chi_min = p->chi_lower;
      p->chi_max = p->chi_upper;
      break;

    case SHAPE_TOPHAT_SMOOTH:
      p->chi_min = fmax (0.0, p->chi_lower - p->n_sigma * p->chi_sigma);
      p->chi_max = p->chi_upper + p->n_sigma * p->chi_sigma;
      break;

    case SHAPE_STUDENT_T:
      p->chi_min = fmax (0.0, p->chi_mean - p->n_scale * p->chi_scale);
      p->chi_max = p->chi_mean + p->n_scale * p->chi_scale;
      break;

    case SHAPE_POWER_EXP:
      p->chi_min = p->chi_lower;
      p->chi_max = p->chi_upper;
      break;

    case SHAPE_LENSING:
      p->chi_min              = p->chi_lower;
      p->chi_max              = p->chi_source_upper;
      p->breaks[p->n_break++] = p->chi_source_lower;
      break;

    case SHAPE_MULTI:
    {
      /* Merge the bumps' supports into groups, as the library does: sort by
       * lower edge, and a bump whose support starts inside the current group
       * extends it. The window is zero between groups, so the gap edges are
       * breakpoints. */
      double lo[MAX_BUMPS], hi[MAX_BUMPS];
      int j;

      for (i = 0; i < p->n_bumps; i++)
      {
        lo[i] = fmax (0.0, p->mu[i] - p->n_sigma * p->sg[i]);
        hi[i] = p->mu[i] + p->n_sigma * p->sg[i];
      }

      for (i = 1; i < p->n_bumps; i++)
        for (j = i; (j > 0) && (lo[j] < lo[j - 1]); j--)
        {
          double s;

          s         = lo[j];
          lo[j]     = lo[j - 1];
          lo[j - 1] = s;
          s         = hi[j];
          hi[j]     = hi[j - 1];
          hi[j - 1] = s;
        }

      p->n_groups = 0;

      for (i = 0; i < p->n_bumps; i++)
      {
        if ((p->n_groups > 0) && (lo[i] <= p->grp_hi[p->n_groups - 1]))
        {
          p->grp_hi[p->n_groups - 1] = fmax (p->grp_hi[p->n_groups - 1], hi[i]);
        }
        else
        {
          p->grp_lo[p->n_groups] = lo[i];
          p->grp_hi[p->n_groups] = hi[i];
          p->n_groups++;
        }
      }

      p->chi_min = p->grp_lo[0];
      p->chi_max = p->grp_hi[p->n_groups - 1];

      for (i = 1; i < p->n_groups; i++)
      {
        p->breaks[p->n_break++] = p->grp_hi[i - 1];
        p->breaks[p->n_break++] = p->grp_lo[i];
      }

      break;
    }

    default:
      fprintf (stderr, "shape_support: bad shape\n");
      exit (1);
  }
}

static Shape
parse_shape (const char *s)
{
  int i;

  for (i = 0; i < SHAPE_LEN; i++)
    if (strcmp (s, shape_names[i]) == 0)
      return (Shape) i;

  fprintf (stderr, "unknown shape '%s'\n", s);
  exit (1);
}

static void
parse_list (const char *s, double *out, int *n)
{
  char *copy = strdup (s);
  char *tok  = strtok (copy, ",");

  *n = 0;

  while ((tok != NULL) && (*n < MAX_BUMPS))
  {
    out[(*n)++] = atof (tok);
    tok         = strtok (NULL, ",");
  }

  free (copy);
}

#endif /* _XCOR_WINDOW_ARB_H */

