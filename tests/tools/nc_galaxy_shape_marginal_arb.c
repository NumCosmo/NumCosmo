/***************************************************************************
 *            nc_galaxy_shape_marginal_arb.c
 *
 *  Mon November 09 12:00:00 2026
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
 * Certified reference values for the weak-lensing shape marginal,
 *
 *   P(eps_obs | g) = INT_{|chi_I| < 1} d^2chi_I P_2D(chi_I) N_2(eps_obs - f_g(chi_I)),
 *   N_2(d) = exp (-|d|^2 / 2 s_n^2) / (2 pi s_n^2),
 *
 * with s_n the per-component noise dispersion (NcGalaxyShapeFactorData:std_noise) and
 * f_g the forward shear map of the TRACE or TRACE_DET convention, both written here from
 * their formulas (no NumCosmo code is used):
 *
 *   TRACE:              f = (chi + g (g chi^* + 2)) / (1 + |g|^2 + 2 Re (g chi^*)),
 *   TRACE_DET, |g|<=1:  f = (chi + g) / (1 + g^* chi),
 *   TRACE_DET, |g|>1:   f = (1 + g chi^*) / (chi^* + g^*).
 *
 * Populations, as the radial density p(r) of r = |chi_I|, P_2D = p(r) / (2 pi r):
 *
 *   gauss sigma:      P_2D = exp (-r^2 / 2 sigma^2) / Z,
 *                     Z = 2 pi sigma^2 (1 - exp (-1 / 2 sigma^2)),
 *   beta alpha beta:  p(r) = r^(alpha - 1) (1 - r)^(beta - 1) / B(alpha, beta).
 *
 * These are NcGalaxyShapePopGauss (truncated to the unit disc) and NcGalaxyShapePopBeta.
 * Both depend on |chi_I| only, which the reductions below use.
 *
 * Reductions, exact for populations that depend on |chi_I| only, with h = g / |g|^2:
 *
 *   |g| > 1: P(eps_obs | g) = P(eps_obs | h). For TRACE the map depends on g only
 *            through 2 g / (1 + |g|^2), which is the same for g and h. For TRACE_DET
 *            the strong map is f_g (chi) = f_h ((g / g^*) chi^*), the weak map at h after
 *            a reflection, which preserves the measure and |chi|.
 *   |g| = 1: every point of the open disc maps to g, so P = N_2(eps_obs - g).
 *
 * h is computed in ball arithmetic from the input g, so the result is certified for the
 * input. After the reductions only the weak maps with |g| < 1 are integrated.
 *
 * Two coordinate systems are available (--coords); both give the same integral.
 *
 * polar (default): chi_I = r e^{i theta}, d^2chi_I = r dr dtheta, theta inside and r
 *   outside. The integrand is p(r) N_2 / (2 pi), and the population factor, which
 *   depends on r only, multiplies the inner integral. As |g| -> 1 the pole of the map
 *   approaches the unit circle at the direction arg (g) + pi, and the integrand varies on
 *   the scale of its distance to the circle there.
 *
 * uv: with g rotated to the positive real axis, z = (1 + eps) / (1 - eps) maps the disc
 *   onto the right half-plane, eps = chi_I for TRACE_DET and eps = T^{-1} (chi_I),
 *   T (eps) = 2 eps / (1 + |eps|^2), for TRACE (f^TRACE_g = T f^TRACE_DET_g T^{-1}). In
 *   z the TRACE_DET map is z -> lambda z, lambda = (1 + |g|) / (1 - |g|), and with
 *   z = e^{u + i v}, v in (-pi/2, pi/2), it is the translation u -> u + ln lambda. The
 *   measure is d^2eps = 4 |z|^2 / |z + 1|^4 du dv, times |det T'| =
 *   4 (1 - |eps|^2) / (1 + |eps|^2)^3 for TRACE. No pole comes near the domain for any
 *   |g| < 1. u is integrated over [-U_left, U_right] and the two ends by the
 *   substitutions rho = e^u in (0, e^{-U_left}] and s = e^{-u} in (0, e^{-U_right}],
 *   in which the integrand is analytic up to rho = 0 and s = 0, each split geometrically
 *   toward 0 for a peak at a fixed point of the map. For Beta, two patches: the point
 *   chi_I = 0 is the interior point (u, v) = (0, 0), where |chi_I| = sqrt (x^2 + y^2) has a
 *   branch and P_2D ~ r^(alpha - 2); the square |u|, |v| <= R0 (--uv-r0) is integrated on
 *   four triangles with apex at (0, 0), (u, v) = R0 t d (xi), du dv = R0^2 t dt dxi, where
 *   |tanh (w/2)|^2 = |w|^2 G / 4, G = sinc (iw/2) sinc (iwbar/2) / (cosh (w/2) cosh (wbar/2))
 *   analytic and nonzero on the strip, makes |chi_I| = t M with M analytic, so the only
 *   non-analytic factor is t^(alpha - 1) at the end of the t interval; the rest of the strip
 *   is integrated as rectangles. On the circle the zero of 1 - |chi_I| is factored out
 *   exactly, 1 - |eps|^2 = 4 W cos v / |a + b|^2 with W = rho, e^u or s.
 *
 * Non-integer powers of a ball that contains 0, at the end of an integration interval, are
 * enclosed by the disc of radius (max |x|)^y when only an enclosure is requested
 * (pow_safe ()); acb_pow () would return a non-finite ball through log (0).
 *
 * A complex number is carried as its two real components, each an acb_t: on real
 * arguments this is the ordinary arithmetic, and off the real line it is the analytic
 * continuation of the integrand, which acb_calc_integrate needs. Taking conj () of an
 * acb_t instead would not be analytic. Poles give a division by a ball containing zero,
 * which returns a non-finite ball, and acb_pow_analytic and acb_sqrt_analytic report a
 * ball crossing a branch cut; both make the integrator subdivide.
 *
 * The integrand is evaluated as exp (L0 - L) with L = -ln of the unnormalized integrand
 * and L0 its minimum found in double precision; the shift is exact (L0 is a double
 * taken as given) and is added back to ln P at the end. With it the integrand is of order
 * one at its peak whatever the value of P, which can be far below the double range.
 *
 * The breakpoints only set the cost; the value and its radius are certified whatever
 * they are. They are geometric ladders around the peak of the integrand, around the
 * peak of the population and, in polar, around the point of the unit circle closest to
 * the pole of the map when it is within 0.5 of the circle.
 *
 * --radial computes the same quantity at g = 0 as a one-dimensional integral, using
 * INT dtheta N_2 = exp (-(r^2 + e^2) / 2 s_n^2) I_0 (r e / s_n^2) / s_n^2, e = |eps_obs|.
 * It is an independent check of the nested integration and accepts only g = 0.
 *
 * Input, one point per line on stdin ('#' starts a comment):
 *   conv pop par1 par2 std_noise g1 g2 e1 e2
 * with conv in {trace, trace_det}, and pop gauss (par1 = sigma, par2 ignored) or beta
 * (par1 = alpha, par2 = beta).
 *
 * Output, one TSV row per point: the inputs, ln P, the radius of ln P, the precision
 * that certified it (0 if --prec-max did not suffice), the number of evaluations of the
 * inner integrand and the wall time in seconds.
 *
 * Build (standalone):
 *   gcc -O2 -o nc_galaxy_shape_marginal_arb nc_galaxy_shape_marginal_arb.c \
 *       $(pkg-config --cflags --libs flint) -lm
 */

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

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_EDGES 256

/* uv: length in u of the ranges beyond each end of the core searched for L0. */
#define UV_END_SEARCH 60.0

/* Precision of the input-derived balls (h, ln lambda, the rotated eps_obs). */
#define INPUT_PREC 1024


typedef enum
{
  CONV_TRACE,
  CONV_TRACE_DET,
} Conv;

typedef enum
{
  POP_GAUSS,
  POP_BETA,
} Pop;

typedef enum
{
  COORDS_POLAR,
  COORDS_UV,
} Coords;

/*
 * Outer variable of the uv integral: rho = e^u (left end), u (core), s = e^{-u} (right end),
 * or t on one of the four triangles of the inner square (two patches, see setup_uv ()).
 */
typedef enum
{
  UV_RHO,
  UV_U,
  UV_S,
  UV_TRI,
} UVPiece;

typedef struct
{
  Conv conv;
  Pop pop;
  Coords coords;
  double par1, par2;
  double std_noise;
  double g1, g2, e1, e2;

  /* polar: the shear after the |g| > 1 reduction, as balls and rounded to double */
  acb_t gx, gy;
  double gd1, gd2;

  /* uv: ln lambda and eps_obs rotated by -arg g, as balls and rounded to double */
  acb_t lnlam, ex, ey;
  double lnlam_d, ex_d, ey_d;
  double U_left, U_right, U_margin;

  double L0;
  double r_edges[MAX_EDGES];
  int n_r;
  double th_edges[MAX_EDGES];
  int n_th;

  /* uv, two patches: half-side R0 of the inner square in w, and its edges */
  int patches;
  double R0;
  double thx_edges[MAX_EDGES];
  int n_thx;
  double t_edges[MAX_EDGES];
  int n_t;
  double xi_edges[MAX_EDGES];
  int n_xi;

  /* the current outer ball and, for uv, which outer variable it is */
  acb_t r;
  UVPiece piece;
  int tri;
  int excise;
  int outer_analytic;

  mag_t inner_tol;
  slong inner_eval_limit;
  long n_eval;
} Point;

/* ------------------------------------------------------------------------- */
/* Input reduction                                                            */
/* ------------------------------------------------------------------------- */

/*
 * Sets the balls of @p from its inputs: g, or h = g / |g|^2 when |g| > 1; for uv also
 * ln lambda = 2 atanh |g| and eps_obs rotated by -arg g. Returns 1 when |g| = 1 exactly
 * (|g|^2 of two doubles is computed exactly at INPUT_PREC), 0 otherwise.
 */
static int
reduce_input (Point *p)
{
  const slong prec = INPUT_PREC;
  arb_t g2, t, one;
  int unit = 0;

  arb_init (g2);
  arb_init (t);
  arb_init (one);
  arb_one (one);

  acb_set_d (p->gx, p->g1);
  acb_set_d (p->gy, p->g2);
  arb_set_d (g2, p->g1);
  arb_sqr (g2, g2, prec);
  arb_set_d (t, p->g2);
  arb_addmul (g2, t, t, prec);

  if (arb_is_one (g2))
  {
    unit = 1;
  }
  else if (arb_gt (g2, one))
  {
    acb_div_arb (p->gx, p->gx, g2, prec);
    acb_div_arb (p->gy, p->gy, g2, prec);
  }

  p->gd1 = arf_get_d (arb_midref (acb_realref (p->gx)), ARF_RND_NEAR);
  p->gd2 = arf_get_d (arb_midref (acb_realref (p->gy)), ARF_RND_NEAR);

  if ((p->coords == COORDS_UV) && !unit)
  {
    acb_t gm, c, s, e1, e2;

    acb_init (gm);
    acb_init (c);
    acb_init (s);
    acb_init (e1);
    acb_init (e2);

    /* |g|, then cos (arg g) = g_x / |g|, sin (arg g) = g_y / |g| */
    acb_sqr (gm, p->gx, prec);
    acb_addmul (gm, p->gy, p->gy, prec);
    acb_sqrt (gm, gm, prec);

    if (acb_is_zero (gm))
    {
      acb_one (c);
      acb_zero (s);
      acb_zero (p->lnlam);
    }
    else
    {
      acb_div (c, p->gx, gm, prec);
      acb_div (s, p->gy, gm, prec);
      acb_atanh (p->lnlam, gm, prec);
      acb_mul_2exp_si (p->lnlam, p->lnlam, 1);
    }

    /* eps' = eps e^{-i arg g}: (e1 c + e2 s, e2 c - e1 s) */
    acb_set_d (e1, p->e1);
    acb_set_d (e2, p->e2);
    acb_mul (p->ex, e1, c, prec);
    acb_addmul (p->ex, e2, s, prec);
    acb_mul (p->ey, e2, c, prec);
    acb_submul (p->ey, e1, s, prec);

    p->lnlam_d = arf_get_d (arb_midref (acb_realref (p->lnlam)), ARF_RND_NEAR);
    p->ex_d    = arf_get_d (arb_midref (acb_realref (p->ex)), ARF_RND_NEAR);
    p->ey_d    = arf_get_d (arb_midref (acb_realref (p->ey)), ARF_RND_NEAR);

    acb_clear (gm);
    acb_clear (c);
    acb_clear (s);
    acb_clear (e1);
    acb_clear (e2);
  }

  arb_clear (g2);
  arb_clear (t);
  arb_clear (one);

  return unit;
}

/* ------------------------------------------------------------------------- */
/* Double precision: the peak of the integrand, for L0 and the breakpoints.   */
/* ------------------------------------------------------------------------- */

/* The weak map with the reduced shear, |g| < 1. */
static double complex
map_d (const Point *p, double complex chi)
{
  const double complex g = p->gd1 + I * p->gd2;
  const double g2        = p->gd1 * p->gd1 + p->gd2 * p->gd2;

  if (p->conv == CONV_TRACE)
  {
    const double complex num = chi + g * (g * conj (chi) + 2.0);
    const double den         = 1.0 + g2 + 2.0 * creal (g * conj (chi));

    return num / den;
  }

  return (chi + g) / (1.0 + conj (g) * chi);
}

/* -ln of the unnormalized 2D population density at r: r^2/2s^2, or -(a-2) ln r - (b-1) ln (1-r). */
static double
pop2d_L_d (const Point *p, double r)
{
  double L = 0.0;

  if (p->pop == POP_GAUSS)
    return r * r / (2.0 * p->par1 * p->par1);

  if (p->par1 != 2.0)
    L -= (p->par1 - 2.0) * log (r);

  if (p->par2 != 1.0)
    L -= (p->par2 - 1.0) * log1p (-r);

  return L;
}

/* polar: L(r, theta) = pop2d_L (r) - ln r + |eps - f_g (r e^{i theta})|^2 / 2 s_n^2. */
static double
L_polar_d (const Point *p, double r, double th)
{
  const double complex chi = r * cexp (I * th);
  const double complex d   = (p->e1 + I * p->e2) - map_d (p, chi);
  const double q           = creal (d) * creal (d) + cimag (d) * cimag (d);

  return pop2d_L_d (p, r) - log (r) + q / (2.0 * p->std_noise * p->std_noise);
}

/* uv: eps (u, v) = (z - 1) / (z + 1), z = e^{u + i v}, written to avoid overflow at large u. */
static double complex
eps_uv_d (double u, double v)
{
  if (u < 0.0)
  {
    const double complex a = exp (u) * cexp (I * v);

    return (a - 1.0) / (a + 1.0);
  }
  else
  {
    const double complex a = cexp (I * v);
    const double s         = exp (-u);

    return (a - s) / (a + s);
  }
}

static double complex
T_d (double complex e)
{
  return 2.0 * e / (1.0 + creal (e) * creal (e) + cimag (e) * cimag (e));
}

/* uv: L(u, v) = pop2d_L (|chi_I|) - ln m (u, v) + |eps' - chi_L|^2 / 2 s_n^2. */
static double
L_uv_d (const Point *p, double u, double v)
{
  const double complex eI  = eps_uv_d (u, v);
  const double complex eL  = eps_uv_d (u + p->lnlam_d, v);
  const double complex chI = (p->conv == CONV_TRACE) ? T_d (eI) : eI;
  const double complex chL = (p->conv == CONV_TRACE) ? T_d (eL) : eL;
  const double complex d   = (p->ex_d + I * p->ey_d) - chL;
  const double q           = creal (d) * creal (d) + cimag (d) * cimag (d);
  const double e2          = creal (eI) * creal (eI) + cimag (eI) * cimag (eI);
  /* m = 4 |z|^2 / |z + 1|^4 = |1 - eps^2|^2 / 4 in terms of eps, which avoids overflow */
  double lnm = 2.0 * log (cabs (1.0 - eI * eI)) - log (4.0);

  if (p->conv == CONV_TRACE)
    lnm += log (4.0) + log1p (-e2) - 3.0 * log1p (e2);

  return pop2d_L_d (p, cabs (chI)) - lnm + q / (2.0 * p->std_noise * p->std_noise);
}

typedef double (*LFunc) (const Point *p, double x, double y);

/*
 * Finds the minimum of @L over [x_lo, x_hi] x [y_lo, y_hi] by a grid search followed
 * by repeated zooming. The result is not required to be the exact minimum: it only sets
 * the exponent shift and the breakpoints.
 */
static void
find_peak (const Point *p, LFunc L, double x_lo, double x_hi, double y_lo, double y_hi,
           int nx, int ny, double *x_m, double *y_m, double *L_m)
{
  double best = INFINITY, bx = 0.5 * (x_lo + x_hi), by = 0.5 * (y_lo + y_hi);
  double hx = (x_hi - x_lo) / nx, hy = (y_hi - y_lo) / ny;
  int i, j, round;

  for (i = 0; i < nx; i++)
  {
    for (j = 0; j < ny; j++)
    {
      const double x = x_lo + (i + 0.5) * hx;
      const double y = y_lo + (j + 0.5) * hy;
      const double v = L (p, x, y);

      if (v < best)
      {
        best = v;
        bx   = x;
        by   = y;
      }
    }
  }

  for (round = 0; round < 12; round++)
  {
    const double x0 = bx, y0 = by;

    for (i = -8; i <= 8; i++)
    {
      for (j = -8; j <= 8; j++)
      {
        const double x = x0 + i * hx / 8.0;
        const double y = y0 + j * hy / 8.0;
        double v;

        if ((x <= x_lo) || (x >= x_hi) || (y <= y_lo) || (y >= y_hi))
          continue;

        v = L (p, x, y);

        if (v < best)
        {
          best = v;
          bx   = x;
          by   = y;
        }
      }
    }

    hx /= 4.0;
    hy /= 4.0;
  }

  *x_m = bx;
  *y_m = by;
  *L_m = best;
}

/*
 * Widths of the peak of @L at (x_m, y_m) from its second derivatives (and, in x, from
 * the first, for a peak against a boundary), clamped to [1e-7, x_max] and [1e-7, y_max].
 */
static void
peak_widths (const Point *p, LFunc L, double x_m, double y_m, double L_m, double hx,
             double x_max, double y_max, double *w_x, double *w_y)
{
  const double Lxx = (L (p, x_m + hx, y_m) - 2.0 * L_m + L (p, x_m - hx, y_m)) / (hx * hx);
  const double Lx  = (L (p, x_m + hx, y_m) - L (p, x_m - hx, y_m)) / (2.0 * hx);
  const double Lyy = (L (p, x_m, y_m + 1.0e-4) - 2.0 * L_m + L (p, x_m, y_m - 1.0e-4)) / 1.0e-8;
  double wx, wy;

  wx = (Lxx > 0.0) ? 1.0 / sqrt (Lxx) : 0.1;

  if (fabs (Lx) > 0.0)
    wx = fmin (wx, 1.0 / fabs (Lx));

  wy = (Lyy > 0.0) ? 1.0 / sqrt (Lyy) : 1.0;

  *w_x = fmin (fmax (wx, 1.0e-7), x_max);
  *w_y = fmin (fmax (wy, 1.0e-7), y_max);
}

static int
cmp_double (const void *a, const void *b)
{
  const double x = *(const double *) a, y = *(const double *) b;

  return (x > y) - (x < y);
}

/* Sorts @e, removes points closer than 1e-13 and points outside [lo, hi], keeps lo and hi. */
static int
finalize_edges (double *e, int n, double lo, double hi)
{
  int i, m = 0;

  e[n++] = lo;
  e[n++] = hi;
  qsort (e, n, sizeof (double), cmp_double);

  for (i = 0; i < n; i++)
  {
    if ((e[i] < lo) || (e[i] > hi))
      continue;

    if ((m > 0) && (e[i] - e[m - 1] < 1.0e-13))
      continue;

    e[m++] = e[i];
  }

  return m;
}

/* Appends c, and c +- w 2^k for k = 0, 1, ..., keeping only points inside (lo, hi). */
static int
add_ladder (double *e, int n, double c, double w, double lo, double hi)
{
  double s;

  if (!(w > 0.0))
    return n;

  if ((c > lo) && (c < hi) && (n < MAX_EDGES - 3))
    e[n++] = c;

  for (s = w; s < 2.0 * (hi - lo); s *= 2.0)
  {
    if (n >= MAX_EDGES - 4)
      break;

    if ((c - s > lo) && (c - s < hi))
      e[n++] = c - s;

    if ((c + s > lo) && (c + s < hi))
      e[n++] = c + s;
  }

  return n;
}

/*
 * polar: distance from the unit circle to the pole of the weak map (the zero of its
 * denominator closest to the disc), which lies in the direction arg (g) + pi:
 * (1 - |g|)^2 / 2|g| for TRACE, 1/|g| - 1 for TRACE_DET. Infinite at g = 0.
 */
static double
pole_distance (const Point *p)
{
  const double gm = hypot (p->gd1, p->gd2);

  if (gm == 0.0)
    return INFINITY;

  if (p->conv == CONV_TRACE)
    return (1.0 - gm) * (1.0 - gm) / (2.0 * gm);

  return 1.0 / gm - 1.0;
}

/*
 * polar: sets L0 and the r and theta breakpoints: ladders around the peak of the
 * integrand, around the peak of the population, and around the point of the unit
 * circle closest to the pole of the map. The theta range is [th_m - pi, th_m + pi], th_m
 * the angle of the peak; a ladder near one end is also placed at its copy 2 pi away.
 */
static void
setup_polar (Point *p)
{
  double r_m, th_m, L_m, w_r, w_t, d_pole;
  int k;

  find_peak (p, L_polar_d, 0.0, 1.0, 0.0, 2.0 * M_PI, 256, 256, &r_m, &th_m, &L_m);
  p->L0 = L_m;
  peak_widths (p, L_polar_d, r_m, th_m, L_m, 1.0e-4 * fmin (r_m, 1.0 - r_m), 0.5, M_PI, &w_r, &w_t);

  p->n_r  = 0;
  p->n_th = 0;
  p->n_r  = add_ladder (p->r_edges, p->n_r, r_m, 0.5 * w_r, 0.0, 1.0);
  p->n_th = add_ladder (p->th_edges, p->n_th, th_m, 0.5 * w_t, th_m - M_PI, th_m + M_PI);

  if (p->pop == POP_GAUSS)
  {
    p->n_r = add_ladder (p->r_edges, p->n_r, p->par1, 0.25 * p->par1, 0.0, 1.0);
  }
  else
  {
    const double a = p->par1, b = p->par2;
    const double sd = sqrt (a * b / ((a + b) * (a + b) * (a + b + 1.0)));

    if (a + b > 2.0)
      p->n_r = add_ladder (p->r_edges, p->n_r, (a - 1.0) / (a + b - 2.0), 0.25 * sd, 0.0, 1.0);
  }

  d_pole = pole_distance (p);

  if (d_pole < 0.5)
  {
    const double th_pole = atan2 (p->gd2, p->gd1) + M_PI;
    const double w       = fmax (0.5 * d_pole, 1.0e-12);

    p->n_r = add_ladder (p->r_edges, p->n_r, 1.0, w, 0.0, 1.0);

    for (k = -2; k <= 2; k++)
      p->n_th = add_ladder (p->th_edges, p->n_th, th_pole + 2.0 * M_PI * k, w,
                            th_m - M_PI, th_m + M_PI);
  }

  p->n_r  = finalize_edges (p->r_edges, p->n_r, 0.0, 1.0);
  p->n_th = finalize_edges (p->th_edges, p->n_th, th_m - M_PI, th_m + M_PI);
}

/*
 * uv: sets the core range [-U_left, U_right] of u, L0 and the u and v breakpoints:
 * ladders around the peak of the integrand and around (u, v) = (0, 0), where chi_I = 0.
 * The two ends are integrated in rho = e^u and s = e^{-u}, so the core range only has
 * to keep lambda rho and s below e^{-U_margin}: U_left = ln lambda + U_margin and
 * U_right = U_margin.
 */
static void
setup_uv (Point *p)
{
  double u_m, v_m, L_m, w_u, w_v;

  p->U_left  = p->lnlam_d + p->U_margin;
  p->U_right = p->U_margin;

  find_peak (p, L_uv_d, -p->U_left, p->U_right, -0.5 * M_PI, 0.5 * M_PI, 1024, 128, &u_m, &v_m, &L_m);
  p->L0 = L_m;
  peak_widths (p, L_uv_d, u_m, v_m, L_m, 1.0e-4, 4.0, 0.5 * M_PI, &w_u, &w_v);

  /*
   * For eps_obs near a fixed point of the map (chi = -+ g/|g|) the peak lies in an end
   * piece, beyond the core range; L0 is the smaller of the core and end minima.
   */
  {
    double ue, ve, Le;

    find_peak (p, L_uv_d, -p->U_left - UV_END_SEARCH, -p->U_left, -0.5 * M_PI, 0.5 * M_PI, 256, 64, &ue, &ve, &Le);
    p->L0 = fmin (p->L0, Le);
    find_peak (p, L_uv_d, p->U_right, p->U_right + UV_END_SEARCH, -0.5 * M_PI, 0.5 * M_PI, 256, 64, &ue, &ve, &Le);
    p->L0 = fmin (p->L0, Le);
  }

  p->n_r  = 0;
  p->n_th = 0;
  p->n_r  = add_ladder (p->r_edges, p->n_r, u_m, 0.5 * w_u, -p->U_left, p->U_right);
  p->n_th = add_ladder (p->th_edges, p->n_th, v_m, 0.5 * w_v, -0.5 * M_PI, 0.5 * M_PI);

  if (p->pop == POP_GAUSS)
  {
    p->n_r  = add_ladder (p->r_edges, p->n_r, 0.0, 0.25 * p->par1, -p->U_left, p->U_right);
    p->n_th = add_ladder (p->th_edges, p->n_th, 0.0, 0.25 * p->par1, -0.5 * M_PI, 0.5 * M_PI);
  }
  else
  {
    p->n_r  = add_ladder (p->r_edges, p->n_r, 0.0, 0.05, -p->U_left, p->U_right);
    p->n_th = add_ladder (p->th_edges, p->n_th, 0.0, 0.05, -0.5 * M_PI, 0.5 * M_PI);
  }

  /*
   * Two patches (Beta): the square |u|, |v| <= R0 around w = 0, where chi_I = 0, is
   * integrated on four triangles with apex at w = 0 (t in [0, 1], xi in [-1, 1]); the core
   * panels with |u| <= R0 then skip |v| < R0 in their v integral (thx_edges).
   */
  p->patches = (p->pop == POP_BETA);

  if (p->patches)
  {
    int i, k = 0;
    double tm = -1.0, xm = 0.0;

    p->r_edges[p->n_r++] = -p->R0;
    p->r_edges[p->n_r++] = p->R0;

    p->n_thx = 0;

    for (i = 0; i < p->n_th; i++)
      p->thx_edges[p->n_thx++] = p->th_edges[i];

    p->thx_edges[p->n_thx++] = -p->R0;
    p->thx_edges[p->n_thx++] = p->R0;
    p->n_thx                 = finalize_edges (p->thx_edges, p->n_thx, -0.5 * M_PI, 0.5 * M_PI);

    /* The peak's (t, xi) when it lies in the square, in the triangle containing it. */
    if ((fabs (u_m) < p->R0) && (fabs (v_m) < p->R0))
    {
      if (fabs (u_m) >= fabs (v_m))
      {
        tm = fabs (u_m) / p->R0;
        xm = v_m / fabs (u_m);
      }
      else
      {
        tm = fabs (v_m) / p->R0;
        xm = u_m / fabs (v_m);
      }
    }

    p->n_t = 0;

    for (k = 1; k <= 20; k++)
      p->t_edges[p->n_t++] = ldexp (1.0, -k);

    if (tm > 0.0)
      p->n_t = add_ladder (p->t_edges, p->n_t, tm, 0.5 * fmin (w_u, w_v) / p->R0, 0.0, 1.0);

    p->n_t = finalize_edges (p->t_edges, p->n_t, 0.0, 1.0);

    p->n_xi = 0;

    if (tm > 0.0)
      p->n_xi = add_ladder (p->xi_edges, p->n_xi, xm, 0.5 * fmin (w_u, w_v) / p->R0, -1.0, 1.0);

    p->n_xi = finalize_edges (p->xi_edges, p->n_xi, -1.0, 1.0);
  }

  p->n_r  = finalize_edges (p->r_edges, p->n_r, -p->U_left, p->U_right);
  p->n_th = finalize_edges (p->th_edges, p->n_th, -0.5 * M_PI, 0.5 * M_PI);
}

/* ------------------------------------------------------------------------- */
/* Ball arithmetic on pairs of components.                                   */
/* ------------------------------------------------------------------------- */

typedef struct
{
  acb_t x;
  acb_t y;
} V2;

static void
v2_init (V2 *a)
{
  acb_init (a->x);
  acb_init (a->y);
}

static void
v2_clear (V2 *a)
{
  acb_clear (a->x);
  acb_clear (a->y);
}

/* c = a b, as complex numbers with components (x, y). @c may not alias @a or @b. */
static void
v2_mul (V2 *c, const V2 *a, const V2 *b, slong prec)
{
  acb_t t;

  acb_init (t);
  acb_mul (c->x, a->x, b->x, prec);
  acb_mul (t, a->y, b->y, prec);
  acb_sub (c->x, c->x, t, prec);
  acb_mul (c->y, a->x, b->y, prec);
  acb_mul (t, a->y, b->x, prec);
  acb_add (c->y, c->y, t, prec);
  acb_clear (t);
}

/* n = |a|^2 = a_x^2 + a_y^2, the analytic continuation of the modulus squared. */
static void
v2_abs2 (acb_t n, const V2 *a, slong prec)
{
  acb_sqr (n, a->x, prec);
  acb_addmul (n, a->y, a->y, prec);
}

/* c = a / b = a b^* / |b|^2, with b^* = (b_x, -b_y). @c may not alias @a or @b. */
static void
v2_div (V2 *c, const V2 *a, const V2 *b, slong prec)
{
  V2 bc;
  acb_t n;

  v2_init (&bc);
  acb_init (n);
  acb_set (bc.x, b->x);
  acb_neg (bc.y, b->y);
  v2_mul (c, a, &bc, prec);
  v2_abs2 (n, b, prec);
  acb_div (c->x, c->x, n, prec);
  acb_div (c->y, c->y, n, prec);
  acb_clear (n);
  v2_clear (&bc);
}

/*
 * out = x^y (principal branch), y >= 0 a double. With @analytic set, a ball crossing the
 * branch cut gives a non-finite result (acb_pow_analytic), which makes the integrator
 * subdivide. Otherwise only an enclosure is needed, and a ball containing 0 is enclosed by
 * the disc of radius (max |x|)^y, valid for every branch; acb_pow () would return a
 * non-finite ball there, through log (0), for any y that is not an integer or 1/2.
 */
static void
pow_safe (acb_t out, const acb_t x, double y, int analytic, slong prec)
{
  acb_t e;

  acb_init (e);
  acb_set_d (e, y);

  if (!analytic && (y > 0.0) && acb_contains_zero (x))
  {
    arb_t m, ye;
    arf_t ub;
    mag_t rad;

    arb_init (m);
    arb_init (ye);
    arf_init (ub);
    mag_init (rad);
    acb_abs (m, x, prec);
    arb_get_ubound_arf (ub, m, prec);
    arb_set_arf (m, ub);
    arb_set_d (ye, y);
    arb_pow (m, m, ye, prec);
    arb_get_mag (rad, m);
    acb_zero (out);
    mag_set (arb_radref (acb_realref (out)), rad);
    mag_set (arb_radref (acb_imagref (out)), rad);
    arb_clear (m);
    arb_clear (ye);
    arf_clear (ub);
    mag_clear (rad);
  }
  else
  {
    acb_pow_analytic (out, x, e, analytic, prec);
  }

  acb_clear (e);
}

/* out = 2 s^2, computed from the double s with no rounding before it enters the ball. */
static void
two_sq_acb (acb_t out, double s, slong prec)
{
  acb_set_d (out, s);
  acb_sqr (out, out, prec);
  acb_mul_2exp_si (out, out, 1);
}

/* polar: f = f_g (chi), the weak map with the reduced shear of @p. */
static void
map_acb (V2 *f, const Point *p, const V2 *chi, slong prec)
{
  V2 g, gc, num, den, chic;

  v2_init (&g);
  v2_init (&gc);
  v2_init (&num);
  v2_init (&den);
  v2_init (&chic);

  acb_set (g.x, p->gx);
  acb_set (g.y, p->gy);
  acb_set (chic.x, chi->x);
  acb_neg (chic.y, chi->y);

  if (p->conv == CONV_TRACE)
  {
    acb_t d;

    acb_init (d);

    /* num = chi + g (g chi^* + 2), den = 1 + |g|^2 + 2 Re (g chi^*) */
    v2_mul (&gc, &g, &chic, prec);
    acb_mul_2exp_si (d, gc.x, 1);
    acb_add_si (gc.x, gc.x, 2, prec);
    v2_mul (&num, &g, &gc, prec);
    acb_add (num.x, num.x, chi->x, prec);
    acb_add (num.y, num.y, chi->y, prec);

    v2_abs2 (den.x, &g, prec);
    acb_add_si (den.x, den.x, 1, prec);
    acb_add (d, d, den.x, prec);
    acb_div (f->x, num.x, d, prec);
    acb_div (f->y, num.y, d, prec);

    acb_clear (d);
  }
  else
  {
    V2 gcj;

    v2_init (&gcj);

    /* (chi + g) / (1 + g^* chi) */
    acb_add (num.x, chi->x, g.x, prec);
    acb_add (num.y, chi->y, g.y, prec);
    acb_set (gcj.x, g.x);
    acb_neg (gcj.y, g.y);
    v2_mul (&den, &gcj, chi, prec);
    acb_add_si (den.x, den.x, 1, prec);
    v2_div (f, &num, &den, prec);

    v2_clear (&gcj);
  }

  v2_clear (&g);
  v2_clear (&gc);
  v2_clear (&num);
  v2_clear (&den);
  v2_clear (&chic);
}

/*
 * polar: the unnormalized radial population factor on a ball, r exp (-r^2 / 2 sigma^2)
 * or r^(alpha - 1) (1 - r)^(beta - 1). With @analytic set, a non-integer power whose
 * ball crosses the branch cut returns a non-finite ball.
 */
static void
pop_r_acb (acb_t out, const Point *p, const acb_t r, int analytic, slong prec)
{
  acb_t t, u;

  acb_init (t);
  acb_init (u);

  if (p->pop == POP_GAUSS)
  {
    acb_sqr (t, r, prec);
    two_sq_acb (u, p->par1, prec);
    acb_div (t, t, u, prec);
    acb_neg (t, t);
    acb_exp (t, t, prec);
    acb_mul (out, t, r, prec);
  }
  else
  {
    /* alpha - 1 and beta - 1 are exact in double for alpha, beta >= 1. */
    pow_safe (t, r, p->par1 - 1.0, analytic, prec);
    acb_sub_si (out, r, 1, prec);
    acb_neg (out, out);
    pow_safe (u, out, p->par2 - 1.0, analytic, prec);
    acb_mul (out, t, u, prec);
  }

  acb_clear (t);
  acb_clear (u);
}

/* out = exp (L0 - |e - f|^2 / 2 s_n^2), with e = (ex, ey). */
static void
noise_factor (acb_t out, const Point *p, const acb_t ex, const acb_t ey, const V2 *f, slong prec)
{
  acb_t q, t, d;

  acb_init (q);
  acb_init (t);
  acb_init (d);

  acb_sub (d, ex, f->x, prec);
  acb_sqr (q, d, prec);
  acb_sub (d, ey, f->y, prec);
  acb_addmul (q, d, d, prec);
  two_sq_acb (t, p->std_noise, prec);
  acb_div (q, q, t, prec);
  acb_set_d (t, p->L0);
  acb_sub (q, t, q, prec);
  acb_exp (out, q, prec);

  acb_clear (q);
  acb_clear (t);
  acb_clear (d);
}

/* polar: out = exp (L0 - |eps - f_g (r e^{i theta})|^2 / 2 s_n^2). */
static void
noise_polar (acb_t out, const Point *p, const acb_t r, const acb_t th, slong prec)
{
  V2 chi, f;
  acb_t c, s, ex, ey;

  v2_init (&chi);
  v2_init (&f);
  acb_init (c);
  acb_init (s);
  acb_init (ex);
  acb_init (ey);

  acb_sin_cos (s, c, th, prec);
  acb_mul (chi.x, r, c, prec);
  acb_mul (chi.y, r, s, prec);
  map_acb (&f, p, &chi, prec);
  acb_set_d (ex, p->e1);
  acb_set_d (ey, p->e2);
  noise_factor (out, p, ex, ey, &f, prec);

  v2_clear (&chi);
  v2_clear (&f);
  acb_clear (c);
  acb_clear (s);
  acb_clear (ex);
  acb_clear (ey);
}

/* eps = (a - b) / (a + b) for complex a and real b. */
static void
eps_from_ab (V2 *e, const V2 *a, const acb_t b, slong prec)
{
  V2 n, d;

  v2_init (&n);
  v2_init (&d);
  acb_sub (n.x, a->x, b, prec);
  acb_set (n.y, a->y);
  acb_add (d.x, a->x, b, prec);
  acb_set (d.y, a->y);
  v2_div (e, &n, &d, prec);
  v2_clear (&n);
  v2_clear (&d);
}

/* chi = T (eps) = 2 eps / (1 + |eps|^2), and jac = |det T'| = 4 (1 - |eps|^2) / (1 + |eps|^2)^3. */
static void
T_acb (V2 *chi, acb_t jac, const V2 *e, slong prec)
{
  acb_t n, d;

  acb_init (n);
  acb_init (d);
  v2_abs2 (n, e, prec);
  acb_add_si (d, n, 1, prec);
  acb_div (chi->x, e->x, d, prec);
  acb_div (chi->y, e->y, d, prec);
  acb_mul_2exp_si (chi->x, chi->x, 1);
  acb_mul_2exp_si (chi->y, chi->y, 1);

  if (jac != NULL)
  {
    acb_sub_si (n, n, 1, prec);
    acb_neg (n, n);
    acb_mul_2exp_si (n, n, 2);
    acb_pow_ui (d, d, 3, prec);
    acb_div (jac, n, d, prec);
  }

  acb_clear (n);
  acb_clear (d);
}

/*
 * uv: the unnormalized 2D population density at chi, exp (-|chi|^2 / 2 sigma^2) or
 * |chi|^(alpha - 2) (1 - |chi|)^(beta - 1), with branch crossings reported as non-finite.
 */
static void
pop2d_acb (acb_t out, const Point *p, const V2 *chi, slong prec)
{
  acb_t R2, t, u;

  acb_init (R2);
  acb_init (t);
  acb_init (u);
  v2_abs2 (R2, chi, prec);

  if (p->pop == POP_GAUSS)
  {
    two_sq_acb (t, p->par1, prec);
    acb_div (t, R2, t, prec);
    acb_neg (t, t);
    acb_exp (out, t, prec);
  }
  else
  {
    /* |chi|^(alpha - 2) = (R2)^((alpha - 2) / 2); (alpha - 2) / 2 is exact for alpha >= 2. */
    acb_set_d (u, 0.5 * (p->par1 - 2.0));
    acb_pow_analytic (t, R2, u, 1, prec);
    acb_sqrt_analytic (R2, R2, 1, prec);
    acb_sub_si (R2, R2, 1, prec);
    acb_neg (R2, R2);
    acb_set_d (u, p->par2 - 1.0);
    acb_pow_analytic (u, R2, u, 1, prec);
    acb_mul (out, t, u, prec);
  }

  acb_clear (R2);
  acb_clear (t);
  acb_clear (u);
}

/*
 * uv: the integrand exp (L0) P_2D (chi_I) N_2 (eps' - chi_L) times the measure, at the
 * outer ball x of the current piece and v. With a = rho e^{iv}, b = 1 (rho = x or e^x)
 * or a = e^{iv}, b = s (s = x), eps = (a - b) / (a + b), the measure per d(outer) dv is
 * 4 w / |a + b|^4 with w = rho, rho^2 or s, and lambda multiplies a for chi_L.
 */

/*
 * uv geometry: with a and b as in integrand_uv () and @w the factor of the measure for the
 * current outer variable, writes to @out exp (L0) N_2 (eps' - chi_L) times the measure
 * 4 w / |a + b|^4 (times |det T'| for TRACE) per d(outer) dv, and to @cI the intrinsic
 * shape chi_I. The population factor is not included.
 */
static void
uv_geometry (acb_t out, V2 *cI, V2 *eI_out, acb_t apb2, Point *p, V2 *a, const acb_t b, acb_t w, slong prec)
{
  V2 aL, eI, eL, cL;
  acb_t lam, t, jac;

  v2_init (&aL);
  v2_init (&eI);
  v2_init (&eL);
  v2_init (&cL);
  acb_init (lam);
  acb_init (t);
  acb_init (jac);

  acb_exp (lam, p->lnlam, prec);
  acb_mul (aL.x, a->x, lam, prec);
  acb_mul (aL.y, a->y, lam, prec);
  eps_from_ab (&eI, a, b, prec);
  eps_from_ab (&eL, &aL, b, prec);

  /* measure 4 w / |a + b|^4 */
  acb_add (a->x, a->x, b, prec);
  v2_abs2 (t, a, prec);

  if (apb2 != NULL)
    acb_set (apb2, t);

  if (eI_out != NULL)
  {
    acb_set (eI_out->x, eI.x);
    acb_set (eI_out->y, eI.y);
  }

  acb_sqr (t, t, prec);
  acb_div (w, w, t, prec);
  acb_mul_2exp_si (w, w, 2);

  if (p->conv == CONV_TRACE)
  {
    T_acb (cI, jac, &eI, prec);
    T_acb (&cL, NULL, &eL, prec);
    acb_mul (w, w, jac, prec);
  }
  else
  {
    acb_set (cI->x, eI.x);
    acb_set (cI->y, eI.y);
    acb_set (cL.x, eL.x);
    acb_set (cL.y, eL.y);
  }

  noise_factor (t, p, p->ex, p->ey, &cL, prec);
  acb_mul (out, w, t, prec);

  v2_clear (&aL);
  v2_clear (&eI);
  v2_clear (&eL);
  v2_clear (&cL);
  acb_clear (lam);
  acb_clear (t);
  acb_clear (jac);
}

/*
 * Two patches: the integrand on triangle p->tri of the inner square, at outer t and inner
 * xi, without the factor t^(alpha - 1) (Beta) or t (Gauss) that outer_f () applies. On the
 * triangle w = R0 t d (xi), d = (1, xi), (-1, xi), (xi, 1) or (xi, -1), du dv = R0^2 t dt dxi,
 * and |eps|^2 = |w|^2 G / 4 with G = sinc (i w/2) sinc (i wbar/2) / (cosh (w/2) cosh (wbar/2))
 * analytic and nonzero on the strip, so |chi_I| = t M with M analytic: M = (R0/2)
 * sqrt ((1 + xi^2) G) for TRACE_DET and 2 M_eps / (1 + t^2 M_eps^2) for TRACE. The Beta
 * factor is then t^(alpha - 2) M^(alpha - 2) (1 - t M)^(beta - 1).
 */
static void
integrand_tri (acb_t out, Point *p, const acb_t tt, const acb_t xi, slong prec)
{
  V2 a, cI;
  acb_t u, v, b, w, c, s, t, ww, wb, G, M;

  v2_init (&a);
  v2_init (&cI);
  acb_init (u);
  acb_init (v);
  acb_init (b);
  acb_init (w);
  acb_init (c);
  acb_init (s);
  acb_init (t);
  acb_init (ww);
  acb_init (wb);
  acb_init (G);
  acb_init (M);

  /* (u, v) = R0 t d (xi) */
  acb_set_d (t, p->R0);
  acb_mul (t, t, tt, prec);

  switch (p->tri)
  {
    case 0:
      acb_set (u, t);
      acb_mul (v, t, xi, prec);
      break;
    case 1:
      acb_neg (u, t);
      acb_mul (v, t, xi, prec);
      break;
    case 2:
      acb_mul (u, t, xi, prec);
      acb_set (v, t);
      break;
    default:
      acb_mul (u, t, xi, prec);
      acb_neg (v, t);
      break;
  }

  /* a = e^u (cos v, sin v), b = 1, measure factor e^{2u} per du dv */
  acb_sin_cos (s, c, v, prec);
  acb_exp (t, u, prec);
  acb_mul (a.x, t, c, prec);
  acb_mul (a.y, t, s, prec);
  acb_one (b);
  acb_sqr (w, t, prec);
  uv_geometry (out, &cI, NULL, NULL, p, &a, b, w, prec);

  /* du dv = R0^2 t dt dxi; the t goes to outer_f () */
  acb_set_d (t, p->R0 * p->R0);
  acb_mul (out, out, t, prec);

  if (p->pop == POP_GAUSS)
  {
    pop2d_acb (t, p, &cI, prec);
    acb_mul (out, out, t, prec);
  }
  else
  {
    /* ww = w/2, wb = wbar/2 with w = u + i v */
    acb_mul_onei (t, v);
    acb_add (ww, u, t, prec);
    acb_sub (wb, u, t, prec);
    acb_mul_2exp_si (ww, ww, -1);
    acb_mul_2exp_si (wb, wb, -1);

    /* G = sinc (i ww) sinc (i wb) / (cosh ww cosh wb) */
    acb_mul_onei (t, ww);
    acb_sinc (G, t, prec);
    acb_mul_onei (t, wb);
    acb_sinc (t, t, prec);
    acb_mul (G, G, t, prec);
    acb_cosh (t, ww, prec);
    acb_div (G, G, t, prec);
    acb_cosh (t, wb, prec);
    acb_div (G, G, t, prec);

    /* M_eps = (R0/2) sqrt ((1 + xi^2) G) */
    acb_sqr (t, xi, prec);
    acb_add_si (t, t, 1, prec);
    acb_mul (t, t, G, prec);
    acb_sqrt_analytic (M, t, 1, prec);
    acb_set_d (t, 0.5 * p->R0);
    acb_mul (M, M, t, prec);

    if (p->conv == CONV_TRACE)
    {
      /* M = 2 M_eps / (1 + t^2 M_eps^2) */
      acb_mul (t, tt, M, prec);
      acb_sqr (t, t, prec);
      acb_add_si (t, t, 1, prec);
      acb_div (M, M, t, prec);
      acb_mul_2exp_si (M, M, 1);
    }

    /* M^(alpha - 2) (1 - t M)^(beta - 1); alpha - 2 and beta - 1 are exact doubles */
    acb_set_d (t, p->par1 - 2.0);
    acb_pow_analytic (s, M, t, 1, prec);
    acb_mul (out, out, s, prec);
    acb_mul (t, tt, M, prec);
    acb_sub_si (t, t, 1, prec);
    acb_neg (t, t);
    acb_set_d (s, p->par2 - 1.0);
    acb_pow_analytic (t, t, s, 1, prec);
    acb_mul (out, out, t, prec);
  }

  v2_clear (&a);
  v2_clear (&cI);
  acb_clear (u);
  acb_clear (v);
  acb_clear (b);
  acb_clear (w);
  acb_clear (c);
  acb_clear (s);
  acb_clear (t);
  acb_clear (ww);
  acb_clear (wb);
  acb_clear (G);
  acb_clear (M);
}

/*
 * uv, Beta: |chi_I|^(alpha - 2) (1 - |chi_I|)^(beta - 1) at chi_I = @cI, with the zero of
 * 1 - |chi_I| on the circle (v = +-pi/2) factored out exactly. With W = rho, e^u or s and
 * c = cos v, 1 - |eps|^2 = 4 W c / |a + b|^2, so 1 - |chi_I| = (W c D)^m E with
 * D = 4 / (|a + b|^2 (1 + |eps|)), m = 1, E = 1 for TRACE_DET and m = 2,
 * E = 1 / (1 + |eps|^2) for TRACE (1 - |chi_I| = (1 - |eps|)^2 / (1 + |eps|^2)). The branch of
 * c^(m (beta - 1)) is at v = +-pi/2 for every u and uses the inner integrator's flag
 * @inner_an; that of W^(m (beta - 1)) is at W = 0 (the ends rho, s -> 0) and uses the outer
 * one's.
 */
static void
pop2d_beta_uv (acb_t out, const Point *p, const V2 *cI, const V2 *eI, const acb_t apb2,
               const acb_t W, const acb_t c, int inner_an, slong prec)
{
  const int m     = (p->conv == CONV_TRACE) ? 2 : 1;
  const double ex = m * (p->par2 - 1.0);
  acb_t R2, e2, t, u, D;

  acb_init (R2);
  acb_init (e2);
  acb_init (t);
  acb_init (u);
  acb_init (D);

  /* |chi_I|^(alpha - 2), away from chi_I = 0 */
  v2_abs2 (R2, cI, prec);
  acb_set_d (u, 0.5 * (p->par1 - 2.0));
  acb_pow_analytic (out, R2, u, 1, prec);

  /* D = 4 / (|a + b|^2 (1 + |eps|)) */
  v2_abs2 (e2, eI, prec);
  acb_sqrt_analytic (t, e2, 1, prec);
  acb_add_si (t, t, 1, prec);
  acb_mul (t, t, apb2, prec);
  acb_set_si (D, 4);
  acb_div (D, D, t, prec);

  pow_safe (t, c, ex, inner_an, prec);
  acb_mul (out, out, t, prec);

  if (p->piece == UV_U)
  {
    /* W = e^u: W^ex = e^{ex u}, analytic */
    acb_log (t, W, prec);
    acb_set_d (u, ex);
    acb_mul (t, t, u, prec);
    acb_exp (t, t, prec);
  }
  else
  {
    pow_safe (t, W, ex, p->outer_analytic, prec);
  }

  acb_mul (out, out, t, prec);
  acb_set_d (u, ex);
  acb_pow_analytic (t, D, u, 1, prec);
  acb_mul (out, out, t, prec);

  if (p->conv == CONV_TRACE)
  {
    /* E^(beta - 1) = (1 + |eps|^2)^(1 - beta) */
    acb_add_si (t, e2, 1, prec);
    acb_set_d (u, 1.0 - p->par2);
    acb_pow_analytic (t, t, u, 1, prec);
    acb_mul (out, out, t, prec);
  }

  acb_clear (R2);
  acb_clear (e2);
  acb_clear (t);
  acb_clear (u);
  acb_clear (D);
}

static void
integrand_uv (acb_t out, Point *p, const acb_t x, const acb_t v, int inner_an, slong prec)
{
  V2 a, cI, eI;
  acb_t b, w, c, s, t, W, apb2;

  v2_init (&a);
  v2_init (&cI);
  v2_init (&eI);
  acb_init (b);
  acb_init (w);
  acb_init (c);
  acb_init (s);
  acb_init (t);
  acb_init (W);
  acb_init (apb2);

  acb_sin_cos (s, c, v, prec);

  switch (p->piece)
  {
    case UV_RHO:
      acb_mul (a.x, x, c, prec);
      acb_mul (a.y, x, s, prec);
      acb_one (b);
      acb_set (w, x);
      acb_set (W, x);
      break;
    case UV_U:
      acb_exp (t, x, prec);
      acb_mul (a.x, t, c, prec);
      acb_mul (a.y, t, s, prec);
      acb_one (b);
      acb_sqr (w, t, prec);
      acb_set (W, t);
      break;
    case UV_S:
    default:
      acb_set (a.x, c);
      acb_set (a.y, s);
      acb_set (b, x);
      acb_set (w, x);
      acb_set (W, x);
      break;
  }

  uv_geometry (out, &cI, &eI, apb2, p, &a, b, w, prec);

  if (p->pop == POP_BETA)
    pop2d_beta_uv (t, p, &cI, &eI, apb2, W, c, inner_an, prec);
  else
    pop2d_acb (t, p, &cI, prec);

  acb_mul (out, out, t, prec);

  v2_clear (&a);
  v2_clear (&cI);
  v2_clear (&eI);
  acb_clear (W);
  acb_clear (apb2);
  acb_clear (b);
  acb_clear (w);
  acb_clear (c);
  acb_clear (s);
  acb_clear (t);
}

/* The inner integrand at the current outer ball p->r: theta (polar) or v (uv). */
static void
inner_eval (acb_t out, Point *p, const acb_t y, int analytic, slong prec)
{
  if (p->coords == COORDS_POLAR)
    noise_polar (out, p, p->r, y, prec);
  else if (p->piece == UV_TRI)
    integrand_tri (out, p, p->r, y, prec);
  else
    integrand_uv (out, p, p->r, y, analytic, prec);
}

static int
inner_f (acb_ptr out, const acb_t y, void *param, slong order, slong prec)
{
  Point *p = (Point *) param;

  if (order > 1)
    return 0;

  p->n_eval++;
  inner_eval (out, p, y, order != 0, prec);

  return 0;
}

/*
 * Adds to @res @sign times the enclosure (end - edge) f ([edge, end]) of the inner
 * integral from the double @edge to the exact @end, by the mean value theorem: sign 1
 * for an upper end, -1 for a lower end (the integral from @end to @edge).
 */
static void
add_inner_sliver (acb_t res, Point *p, double edge, const acb_t end, int sign, slong prec)
{
  acb_t e, span, f;

  acb_init (e);
  acb_init (span);
  acb_init (f);
  acb_set_d (e, edge);
  acb_union (span, e, end, prec);
  inner_eval (f, p, span, 0, prec);
  acb_sub (span, end, e, prec);
  acb_mul (f, f, span, prec);

  if (sign < 0)
    acb_sub (res, res, f, prec);
  else
    acb_add (res, res, f, prec);

  acb_clear (e);
  acb_clear (span);
  acb_clear (f);
}

/*
 * Integrates inner_f over the inner range at the current outer ball p->r: the panels
 * between the double edges, plus the slivers between the outer edges and the exact ends
 * of the range, [theta_0, theta_0 + 2 pi] in polar and [-pi/2, pi/2] in uv, which the
 * doubles cannot represent. A panel that reaches the evaluation cap returns a
 * non-finite ball; it is replaced by the enclosure (B - A) f ([A, B]), which is finite
 * whenever the integrand is, so that the outer integrator subdivides rather than
 * receiving an infinity.
 */
static void
inner_integral (acb_t res, Point *p, slong prec)
{
  acb_calc_integrate_opt_t opt;
  acb_t A, B, part;
  const double *edges = p->th_edges;
  int n               = p->n_th;
  int i;

  /* The triangles integrate xi over [-1, 1]; the excised core panels skip |v| < R0. */
  if ((p->coords == COORDS_UV) && (p->piece == UV_TRI))
  {
    edges = p->xi_edges;
    n     = p->n_xi;
  }
  else if ((p->coords == COORDS_UV) && p->excise)
  {
    edges = p->thx_edges;
    n     = p->n_thx;
  }

  acb_calc_integrate_opt_init (opt);
  opt->eval_limit = p->inner_eval_limit;
  acb_init (A);
  acb_init (B);
  acb_init (part);
  acb_zero (res);

  for (i = 0; i + 1 < n; i++)
  {
    if (p->excise && (edges[i] >= -p->R0) && (edges[i + 1] <= p->R0))
      continue;

    acb_set_d (A, edges[i]);
    acb_set_d (B, edges[i + 1]);
    acb_calc_integrate (part, inner_f, p, A, B, prec, p->inner_tol, opt, prec);

    if (!acb_is_finite (part))
    {
      acb_t span, f;

      acb_init (span);
      acb_init (f);
      acb_union (span, A, B, prec);
      inner_eval (f, p, span, 0, prec);
      acb_sub (span, B, A, prec);
      acb_mul (part, f, span, prec);
      acb_clear (span);
      acb_clear (f);
    }

    acb_add (res, res, part, prec);
  }

  acb_const_pi (B, prec);

  if (p->coords == COORDS_POLAR)
  {
    acb_mul_2exp_si (B, B, 1);
    acb_set_d (A, edges[0]);
    acb_add (B, B, A, prec);
    add_inner_sliver (res, p, edges[n - 1], B, 1, prec);
  }
  else if (p->piece != UV_TRI)
  {
    acb_mul_2exp_si (B, B, -1);
    acb_neg (A, B);
    add_inner_sliver (res, p, edges[n - 1], B, 1, prec);
    add_inner_sliver (res, p, edges[0], A, -1, prec);
  }

  acb_clear (A);
  acb_clear (B);
  acb_clear (part);
}

/* polar: pop_r (r) times the theta integral; uv: the v integral. */
static int
outer_f (acb_ptr out, const acb_t x, void *param, slong order, slong prec)
{
  Point *p = (Point *) param;

  if (order > 1)
    return 0;

  p->outer_analytic = (order != 0);

  if (p->coords == COORDS_POLAR)
  {
    acb_t w;

    acb_init (w);
    pop_r_acb (w, p, x, order != 0, prec);

    if (!acb_is_finite (w))
    {
      acb_indeterminate (out);
    }
    else
    {
      acb_set (p->r, x);
      inner_integral (out, p, prec);
      acb_mul (out, out, w, prec);
    }

    acb_clear (w);
  }
  else if (p->piece == UV_TRI)
  {
    acb_t w, e;

    /* t^(alpha - 1) (Beta) or t (Gauss), with the integrator's analyticity request */
    acb_init (w);
    acb_init (e);
    pow_safe (w, x, (p->pop == POP_BETA) ? p->par1 - 1.0 : 1.0, order != 0, prec);

    if (!acb_is_finite (w))
    {
      acb_indeterminate (out);
    }
    else
    {
      acb_set (p->r, x);
      inner_integral (out, p, prec);
      acb_mul (out, out, w, prec);
    }

    acb_clear (w);
    acb_clear (e);
  }
  else
  {
    acb_set (p->r, x);
    inner_integral (out, p, prec);
  }

  return 0;
}

/* g = 0 only: pop_r (r) exp (L0 - (r - e)^2/2s_n^2) 2 pi I0s (r e / s_n^2). */
static int
radial_f (acb_ptr out, const acb_t r, void *param, slong order, slong prec)
{
  Point *p = (Point *) param;
  acb_t q, t, z, nu, e;

  if (order > 1)
    return 0;

  p->n_eval++;
  acb_init (q);
  acb_init (t);
  acb_init (z);
  acb_init (nu);
  acb_init (e);

  /* e = |eps_obs| */
  acb_set_d (e, p->e1);
  acb_sqr (e, e, prec);
  acb_set_d (t, p->e2);
  acb_addmul (e, t, t, prec);
  acb_sqrt (e, e, prec);

  acb_sub (q, r, e, prec);
  acb_sqr (q, q, prec);
  two_sq_acb (t, p->std_noise, prec);
  acb_div (q, q, t, prec);
  acb_set_d (t, p->L0);
  acb_sub (q, t, q, prec);
  acb_exp (q, q, prec);

  /* z = 2 r e / (2 s_n^2) */
  two_sq_acb (t, p->std_noise, prec);
  acb_mul (z, r, e, prec);
  acb_mul_2exp_si (z, z, 1);
  acb_div (z, z, t, prec);
  acb_zero (nu);
  acb_hypgeom_bessel_i_scaled (t, nu, z, prec);
  acb_mul (q, q, t, prec);
  acb_const_pi (t, prec);
  acb_mul_2exp_si (t, t, 1);
  acb_mul (q, q, t, prec);

  pop_r_acb (t, p, r, order != 0, prec);
  acb_mul (out, q, t, prec);

  acb_clear (q);
  acb_clear (t);
  acb_clear (z);
  acb_clear (nu);
  acb_clear (e);

  return 0;
}

/* Integrates @f over [lo, hi] and adds the result to @I_res. */
static void
add_panel (acb_t I_res, acb_calc_func_t f, Point *p, double lo, double hi, const mag_t tol,
           const acb_calc_integrate_opt_t opt, slong prec)
{
  acb_t A, B, part;

  acb_init (A);
  acb_init (B);
  acb_init (part);
  acb_set_d (A, lo);
  acb_set_d (B, hi);
  acb_calc_integrate (part, f, p, A, B, prec, tol, opt, prec);
  acb_add (I_res, I_res, part, prec);
  acb_clear (A);
  acb_clear (B);
  acb_clear (part);
}

/*
 * uv: adds the integral of one end piece (p->piece, UV_RHO or UV_S) over [0, e^{-U}],
 * where the core range ends at u = -U or u = U. The double e_d = exp (-U) is not the
 * exact end e^{-U}, so the piece is integrated over [0, e_d] and the remainder over
 * [e_d, e^{-U}] is enclosed as (e^{-U} - e_d) f ([e_d, e^{-U}]), which contains it by the
 * mean value theorem.
 */
static void
add_end (acb_t I_res, Point *p, double U, const mag_t tol, const acb_calc_integrate_opt_t opt,
         slong prec)
{
  const double e_d = exp (-U);
  acb_t e, span, f;

  acb_init (e);
  acb_init (span);
  acb_init (f);

  {
    /* Geometric panels toward 0, where a peak at a fixed point of the map (u -> -+inf) lies. */
    mag_t etol;
    int k;

    mag_init (etol);
    mag_set (etol, tol);
    mag_mul_2exp_si (etol, etol, -6);
    add_panel (I_res, outer_f, p, 0.0, ldexp (e_d, -48), etol, opt, prec);

    for (k = 48; k >= 1; k--)
      add_panel (I_res, outer_f, p, ldexp (e_d, -k), ldexp (e_d, -k + 1), etol, opt, prec);

    mag_clear (etol);
  }

  acb_set_d (e, -U);
  acb_exp (e, e, prec);
  acb_set_d (f, e_d);
  acb_union (span, e, f, prec);
  outer_f (f, span, p, 0, prec);
  acb_set_d (span, e_d);
  acb_sub (span, e, span, prec);
  acb_mul (f, f, span, prec);
  acb_add (I_res, I_res, f, prec);

  acb_clear (e);
  acb_clear (span);
  acb_clear (f);
}

/*
 * The unnormalized integral I, over the outer panels (and, for uv, the two ends), at
 * working precision @prec, to an absolute tolerance @tol. The outer tolerance is @tol
 * shared over the outer panels; the inner one is that shared over the inner panels
 * (the outer measure is at most 1 in polar and the integrand of order one in uv).
 */
static void
integrate_I (acb_t I_res, Point *p, int radial, const mag_t tol, slong eval_limit, slong prec)
{
  acb_calc_integrate_opt_t opt;
  mag_t ptol;
  int i;

  acb_calc_integrate_opt_init (opt);
  opt->eval_limit = eval_limit;
  mag_init (ptol);
  acb_zero (I_res);

  mag_set (ptol, tol);
  mag_mul_2exp_si (ptol, ptol, -(slong) (log2 (p->n_r + 2) + 1));
  /* Computed in a local: GCC 16 warns on mag_mul_2exp_si applied to the member directly. */
  {
    mag_t itol;

    mag_init (itol);
    mag_mul_2exp_si (itol, ptol, -(slong) (log2 (p->n_th) + 1));
    mag_set (p->inner_tol, itol);
    mag_clear (itol);
  }

  if (p->coords == COORDS_POLAR)
  {
    for (i = 0; i + 1 < p->n_r; i++)
      add_panel (I_res, radial ? radial_f : outer_f, p, p->r_edges[i], p->r_edges[i + 1], ptol, opt, prec);
  }
  else
  {
    p->piece = UV_RHO;
    add_end (I_res, p, p->U_left, ptol, opt, prec);
    p->piece = UV_U;

    for (i = 0; i + 1 < p->n_r; i++)
    {
      p->excise = p->patches && (p->r_edges[i] >= -p->R0) && (p->r_edges[i + 1] <= p->R0);
      add_panel (I_res, outer_f, p, p->r_edges[i], p->r_edges[i + 1], ptol, opt, prec);
    }

    p->excise = 0;
    p->piece  = UV_S;
    add_end (I_res, p, p->U_right, ptol, opt, prec);

    if (p->patches)
    {
      int k;

      p->piece = UV_TRI;

      for (p->tri = 0; p->tri < 4; p->tri++)
        for (k = 0; k + 1 < p->n_t; k++)
          add_panel (I_res, outer_f, p, p->t_edges[k], p->t_edges[k + 1], ptol, opt, prec);
    }
  }

  mag_clear (ptol);
}

/*
 * ln Z of the 2D population: Z = 2 pi sigma^2 (1 - exp (-1/2 sigma^2)) or 2 pi B (alpha, beta).
 * In polar the integrand carries p(r) = 2 pi r P_2D, in uv P_2D itself; both use this Z.
 */
static void
pop_lnZ (arb_t lnZ, const Point *p, slong prec)
{
  arb_t t, u;

  arb_init (t);
  arb_init (u);

  if (p->pop == POP_GAUSS)
  {
    arb_set_d (u, p->par1);
    arb_sqr (u, u, prec);
    arb_mul_2exp_si (u, u, 1);
    arb_inv (t, u, prec);
    arb_neg (t, t);
    arb_expm1 (t, t, prec);
    arb_neg (t, t);
    arb_mul (t, t, u, prec);
    arb_const_pi (u, prec);
    arb_mul (t, t, u, prec);
    arb_log (lnZ, t, prec);
  }
  else
  {
    arb_set_d (t, p->par1);
    arb_lgamma (lnZ, t, prec);
    arb_set_d (t, p->par2);
    arb_lgamma (u, t, prec);
    arb_add (lnZ, lnZ, u, prec);
    arb_set_d (t, p->par1);
    arb_set_d (u, p->par2);
    arb_add (t, t, u, prec);
    arb_lgamma (u, t, prec);
    arb_sub (lnZ, lnZ, u, prec);
    arb_const_pi (u, prec);
    arb_mul_2exp_si (u, u, 1);
    arb_log (u, u, prec);
    arb_add (lnZ, lnZ, u, prec);
  }

  arb_clear (t);
  arb_clear (u);
}

/* ln (2 pi s_n^2) */
static void
noise_lnZ (arb_t out, const Point *p, slong prec)
{
  arb_t u;

  arb_init (u);
  arb_const_pi (u, prec);
  arb_set_d (out, p->std_noise);
  arb_sqr (out, out, prec);
  arb_mul_2exp_si (out, out, 1);
  arb_mul (out, out, u, prec);
  arb_log (out, out, prec);
  arb_clear (u);
}

/* |g| = 1: ln P = -|eps - g|^2 / 2 s_n^2 - ln (2 pi s_n^2), exact. */
static void
compute_unit_shear (arb_t lnP, const Point *p, slong prec)
{
  arb_t d, t, u;

  arb_init (d);
  arb_init (t);
  arb_init (u);
  arb_set_d (d, p->e1);
  arb_set_d (t, p->g1);
  arb_sub (d, d, t, prec);
  arb_sqr (d, d, prec);
  arb_set_d (t, p->e2);
  arb_set_d (u, p->g2);
  arb_sub (t, t, u, prec);
  arb_addmul (d, t, t, prec);
  arb_set_d (t, p->std_noise);
  arb_sqr (t, t, prec);
  arb_mul_2exp_si (t, t, 1);
  arb_div (d, d, t, prec);
  arb_neg (lnP, d);
  noise_lnZ (t, p, prec);
  arb_sub (lnP, lnP, t, prec);
  arb_clear (d);
  arb_clear (t);
  arb_clear (u);
}

/*
 * ln P for the point @p: ln I - L0 - ln Z_pop - ln (2 pi s_n^2), certified to a radius
 * below @target in ln P, escalating the working precision from 64 bits up to
 * @prec_max. Pass 1 (64 bits, loose tolerance) estimates I; pass 2 uses it as the
 * scale of the absolute tolerance. Returns the precision used, or 0 on failure.
 */
static slong
compute_point (arb_t lnP, Point *p, int radial, double target, slong prec_max, slong eval_limit)
{
  acb_t I_res;
  arb_t t;
  mag_t tol;
  slong prec, used = 0;
  double I_est;

  acb_init (I_res);
  arb_init (t);
  mag_init (tol);

  /*
   * Pass 1: the integrand is at most of order one at its peak if L0 is its minimum. If the
   * search missed a narrow peak, I is far from one (|ln I| > 64); L0 is then shifted by
   * ln I and pass 1 repeated, so that I_est is a finite double.
   */
  mag_set_d (tol, 1.0e-8);
  integrate_I (I_res, p, radial, tol, eval_limit, 64);

  if (acb_is_finite (I_res) && arb_is_positive (acb_realref (I_res)))
  {
    arb_log (t, acb_realref (I_res), 64);

    if (arf_cmpabs_2exp_si (arb_midref (t), 6) > 0)
    {
      p->L0 -= arf_get_d (arb_midref (t), ARF_RND_NEAR);
      integrate_I (I_res, p, radial, tol, eval_limit, 64);
    }
  }

  I_est = arf_get_d (arb_midref (acb_realref (I_res)), ARF_RND_NEAR);

  if (!(I_est > 0.0) || !acb_is_finite (I_res))
    I_est = 1.0e-8;

  for (prec = 64; prec <= prec_max; prec *= 2)
  {
    double rel;

    mag_set_d (tol, target * I_est / 8.0);
    integrate_I (I_res, p, radial, tol, eval_limit, prec);

    if (!acb_is_finite (I_res))
      continue;

    rel = mag_get_d (arb_radref (acb_realref (I_res))) /
          fabs (arf_get_d (arb_midref (acb_realref (I_res)), ARF_RND_NEAR));

    if (rel < target)
    {
      used = prec;
      break;
    }
  }

  if (used == 0)
    prec = prec_max;

  arb_log (lnP, acb_realref (I_res), prec);
  arb_set_d (t, p->L0);
  arb_sub (lnP, lnP, t, prec);
  pop_lnZ (t, p, prec);
  arb_sub (lnP, lnP, t, prec);
  noise_lnZ (t, p, prec);
  arb_sub (lnP, lnP, t, prec);

  acb_clear (I_res);
  arb_clear (t);
  mag_clear (tol);

  return used;
}

static int
parse_point (Point *p, const char *line)
{
  char conv[32], pop[32];
  int n;

  n = sscanf (line, "%31s %31s %lf %lf %lf %lf %lf %lf %lf", conv, pop, &p->par1, &p->par2,
              &p->std_noise, &p->g1, &p->g2, &p->e1, &p->e2);

  if (n != 9)
    return -1;

  if (!strcmp (conv, "trace"))
    p->conv = CONV_TRACE;
  else if (!strcmp (conv, "trace_det"))
    p->conv = CONV_TRACE_DET;
  else
    return -2;

  if (!strcmp (pop, "gauss"))
  {
    p->pop  = POP_GAUSS;
    p->par2 = 0.0;

    if (!(p->par1 > 0.0))
      return -4;
  }
  else if (!strcmp (pop, "beta"))
  {
    p->pop = POP_BETA;

    if (!(p->par1 >= 1.0) || !(p->par2 >= 1.0))
      return -4;
  }
  else
  {
    return -3;
  }

  if (!(p->std_noise > 0.0))
    return -5;

  return 0;
}

int
main (int argc, char **argv)
{
  double target = 1.0e-15;
  slong prec_max = 1024, eval_limit = 0, inner_eval_limit = 2000;
  int radial     = 0, verbose = 0, i;
  Coords coords = COORDS_POLAR;
  char line[1024];
  Point p;

  p.U_margin = 2.0;
  p.R0       = 0.375;
  p.excise   = 0;
  p.patches  = 0;

  for (i = 1; i < argc; i++)
  {
    const char *s = argv[i];

    if (!strncmp (s, "--target-rel=", 13))
    {
      target = atof (s + 13);
    }
    else if (!strncmp (s, "--prec-max=", 11))
    {
      prec_max = atol (s + 11);
    }
    else if (!strncmp (s, "--eval-limit=", 13))
    {
      eval_limit = atol (s + 13);
    }
    else if (!strncmp (s, "--inner-eval-limit=", 19))
    {
      inner_eval_limit = atol (s + 19);
    }
    else if (!strcmp (s, "--coords=polar"))
    {
      coords = COORDS_POLAR;
    }
    else if (!strcmp (s, "--coords=uv"))
    {
      coords = COORDS_UV;
    }
    else if (!strncmp (s, "--uv-margin=", 12))
    {
      p.U_margin = atof (s + 12);
    }
    else if (!strncmp (s, "--uv-r0=", 8))
    {
      p.R0 = atof (s + 8);
    }
    else if (!strcmp (s, "--radial"))
    {
      radial = 1;
    }
    else if (!strcmp (s, "--verbose"))
    {
      verbose = 1;
    }
    else
    {
      fprintf (stderr, "unknown argument '%s'\n", s);
      exit (1);
    }
  }

  if (radial && (coords != COORDS_POLAR))
  {
    fprintf (stderr, "--radial is a polar integral; it does not take --coords=uv\n");
    exit (1);
  }

  acb_init (p.r);
  acb_init (p.gx);
  acb_init (p.gy);
  acb_init (p.lnlam);
  acb_init (p.ex);
  acb_init (p.ey);
  mag_init (p.inner_tol);
  p.inner_eval_limit = inner_eval_limit;
  p.coords           = coords;

  printf ("# conv\tpop\tpar1\tpar2\tstd_noise\tg1\tg2\te1\te2\tlnP\tlnP_radius\tprec\tn_eval\tseconds\n");

  while (fgets (line, sizeof (line), stdin) != NULL)
  {
    struct timespec t0, t1;
    const char *s = line;
    arb_t lnP;
    slong used;
    char *str;
    int st;

    while (*s == ' ' || *s == '\t')
      s++;

    if ((*s == '#') || (*s == '\n') || (*s == '\0'))
      continue;

    st = parse_point (&p, s);

    if (st != 0)
    {
      fprintf (stderr, "cannot parse line (%d): %s", st, line);
      exit (1);
    }

    if (radial && ((p.g1 != 0.0) || (p.g2 != 0.0)))
    {
      fprintf (stderr, "--radial accepts only g = 0: %s", line);
      exit (1);
    }

    clock_gettime (CLOCK_MONOTONIC, &t0);
    p.n_eval = 0;
    arb_init (lnP);

    if (reduce_input (&p))
    {
      compute_unit_shear (lnP, &p, 256);
      used = 256;
    }
    else
    {
      if (coords == COORDS_POLAR)
        setup_polar (&p);
      else
        setup_uv (&p);

      if (verbose)
        fprintf (stderr, "L0 %.6g, %d outer panels, %d inner panels\n", p.L0, p.n_r - 1, p.n_th - 1);

      used = compute_point (lnP, &p, radial, target, prec_max, eval_limit);
    }

    clock_gettime (CLOCK_MONOTONIC, &t1);

    str = arb_get_str (lnP, 25, ARB_STR_MORE | ARB_STR_NO_RADIUS);
    printf ("%s\t%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%s\t%.3e\t%ld\t%ld\t%.3f\n",
            p.conv == CONV_TRACE ? "trace" : "trace_det", p.pop == POP_GAUSS ? "gauss" : "beta",
            p.par1, p.par2, p.std_noise, p.g1, p.g2, p.e1, p.e2, str,
            mag_get_d (arb_radref (lnP)), (long) used, p.n_eval,
            (t1.tv_sec - t0.tv_sec) + 1.0e-9 * (t1.tv_nsec - t0.tv_nsec));
    fflush (stdout);
    flint_free (str);
    arb_clear (lnP);
  }

  acb_clear (p.r);
  acb_clear (p.gx);
  acb_clear (p.gy);
  acb_clear (p.lnlam);
  acb_clear (p.ex);
  acb_clear (p.ey);
  mag_clear (p.inner_tol);
  flint_cleanup ();

  return 0;
}

