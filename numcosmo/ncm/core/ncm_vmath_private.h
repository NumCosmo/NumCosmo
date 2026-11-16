/***************************************************************************
 *            ncm_vmath_private.h
 *
 *  Sat November 14 2026
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * ncm_vmath_private.h
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

#ifndef _NCM_VMATH_PRIVATE_H_
#define _NCM_VMATH_PRIVATE_H_

#include <glib.h>
#include <float.h>
#include <math.h>
#include <string.h>

G_BEGIN_DECLS

/*
 * Elementary functions written so that a loop calling them vectorizes: no branch, every
 * candidate value is computed and the result is chosen by a bitwise blend. They are
 * static inline so that the compiler sees their body inside the caller's loop (an
 * omp simd loop, for instance). Without -ffast-math, as in NumCosmo, glibc declares its
 * vector variants (libmvec) only under __FAST_MATH__, so with GCC 14.3 and glibc 2.44 a
 * loop that calls exp() or log() stays scalar; other toolchains may differ (on the system
 * GCC of the Ubuntu CI image, a loop of log() and scalar calls of it differed in the last
 * bit). They are for the loops that dominate a computation, and libm is the reference for
 * the rest.
 *
 * ncm_vmath_select() is the blend. Under -ftrapping-math (the default) and without masked
 * instructions (SSE2 and AVX2), GCC 14.3 keeps a conditional expression on doubles as
 * control flow, which leaves the loop scalar, except the forms it recognizes as a minimum
 * or a maximum; with AVX-512 it vectorizes them.
 *
 * ncm_vmath_atan2() and ncm_vmath_sincos() are faster than libm with SSE2 already.
 * ncm_vmath_exp(), ncm_vmath_log() and ncm_vmath_pow_pos() are faster only with AVX2 or
 * AVX-512 (two-wide division and no fma with SSE2): there NCM_VMATH_HAVE_VECTOR_EXP_LOG
 * is defined, and
 * elsewhere they are exp(), log() and pow(). The code is chosen when it is compiled, with
 * no dispatch at run time: these three are faster only in a build for AVX2 or AVX-512
 * (-Dmarch-native=enabled, or an -march that has them). The default build and the conda
 * packages (-march=nocona) use libm for them. Results then depend on the build in the
 * last bits.
 */

/* x + ROUND - ROUND rounds x to an integer (ncm_vmath_exp()) and the blends rely on IEEE
 * semantics, which -ffast-math lets the compiler drop */
#ifdef __FAST_MATH__
#error "ncm_vmath_private.h needs IEEE arithmetic: do not build with -ffast-math"
#endif

/*
 * ncm_vmath_select:
 * @cond: a condition
 * @a: value for a true @cond
 * @b: value for a false @cond
 *
 * Returns: @cond ? @a : @b, as a bitwise blend of the two values.
 */
static inline gdouble
ncm_vmath_select (const gboolean cond, const gdouble a, const gdouble b)
{
  const guint64 mask = -(guint64) (cond != 0);
  guint64 ua, ub, ur;
  gdouble r;

  memcpy (&ua, &a, sizeof (ua));
  memcpy (&ub, &b, sizeof (ub));
  ur = (ua & mask) | (ub & ~mask);
  memcpy (&r, &ur, sizeof (r));

  return r;
}

/*
 * ncm_vmath_atan2:
 * @y: a value
 * @x: a value
 *
 * The arc tangent of @y / @x in $[-\pi, \pi]$, with the quadrant given by the signs, as
 * atan2(). It is the double-precision rational approximation of atan on $[0, 1]$ of
 * Cephes (Moshier), applied to min (|x|, |y|) / max (|x|, |y|) and folded back by selects.
 * It is within 2 ulp of atan2(), with the same values at zeros and when one argument is
 * infinite; with both arguments infinite it returns NaN instead of an odd multiple of
 * $\pi/4$, and a NaN argument gives NaN.
 *
 * Returns: atan2 (@y, @x).
 */
static inline gdouble
ncm_vmath_atan2 (const gdouble y, const gdouble x)
{
  const gdouble P0       = -8.750608600031904122785e-1;
  const gdouble P1       = -1.615753718733365076637e1;
  const gdouble P2       = -7.500855792314704667340e1;
  const gdouble P3       = -1.228866684490136173410e2;
  const gdouble P4       = -6.485021904942025371773e1;
  const gdouble Q0       = 2.485846490142306297962e1;
  const gdouble Q1       = 1.650270098316988542046e2;
  const gdouble Q2       = 4.328810604912902668951e2;
  const gdouble Q3       = 4.853903996359136964868e2;
  const gdouble Q4       = 1.945506571482613964425e2;
  const gdouble MOREBITS = 6.123233995736765886130e-17; /* pi / 2 - M_PI_2 */
  const gdouble ay       = fabs (y);
  const gdouble ax       = fabs (x);
  const gboolean x_major = (ax > ay);
  const gdouble mx       = ncm_vmath_select (x_major, ax, ay);
  const gdouble mn       = ncm_vmath_select (x_major, ay, ax);
  const gdouble t0       = mn / ncm_vmath_select (mx == 0.0, 1.0, mx);
  const gboolean big     = (t0 > 0.66);
  const gdouble t1       = (t0 - 1.0) / (t0 + 1.0);
  const gdouble t        = ncm_vmath_select (big, t1, t0);
  const gdouble z        = t * t;
  const gdouble p        = (((P0 * z + P1) * z + P2) * z + P3) * z + P4;
  const gdouble q        = ((((z + Q0) * z + Q1) * z + Q2) * z + Q3) * z + Q4;
  const gdouble r0       = t * (z * p / q) + t;
  const gdouble r1       = M_PI_4 + (r0 + 0.5 * MOREBITS);
  const gdouble r        = ncm_vmath_select (big, r1, r0);
  const gdouble a1       = (M_PI_2 - r) + MOREBITS;
  const gdouble a        = ncm_vmath_select (ay > ax, a1, r);
  const gdouble b1       = (M_PI - a) + 2.0 * MOREBITS;

  /* The sign of x, -0 included so that it gives pi as in atan2(): signbit() reads half
   * of the double and a 64-bit integer test needs SSE4.1, neither vectorizes with SSE2 */
  const gboolean x_neg = (copysign (1.0, x) < 0.0);

  return copysign (ncm_vmath_select (x_neg, b1, a), y);
}

/*
 * ncm_vmath_sincos:
 * @x: a value, $|x| \le 2^{20}$
 * @s: (out): sin (@x)
 * @c: (out): cos (@x)
 *
 * The sine and the cosine of @x. @x is reduced to $r = x - n \pi / 2$, $|r| \le \pi / 4$,
 * with n rounded by adding and subtracting $1.5 \times 2^{52}$ and $\pi / 2$ in three
 * parts (Cody and Waite), of which the first two have products with n that are exact;
 * sin r and cos r are the polynomials of Cephes (Moshier) on $[-\pi/4, \pi/4]$, and the
 * last two bits of n choose and sign them by selects. Up to $|x| = 10^5$ both are within
 * 2 ulp of sin() and cos() where their modulus exceeds $10^{-3}$ and within $2.3 \times
 * 10^{-16}$ of them elsewhere; sin (-0) is +0. Beyond $|x| = 2^{20}$ the reduction loses
 * digits; an infinite or NaN @x gives NaN.
 */
static inline void
ncm_vmath_sincos (const gdouble x, gdouble *s, gdouble *c)
{
  const gdouble ROUND = 6755399441055744.0;     /* 1.5 2^52: x + ROUND - ROUND rounds x to an integer */
  const gdouble TWOPI = 0.63661977236758134308; /* 2 / pi */
  const gdouble DP1   = 1.57079625129699707031e0;
  const gdouble DP2   = 7.54978941586159635335e-8;
  const gdouble DP3   = 5.39030285815811905290e-15;
  const gdouble S0    = 1.58962301576546568060e-10;
  const gdouble S1    = -2.50507477628578072866e-8;
  const gdouble S2    = 2.75573136213857245213e-6;
  const gdouble S3    = -1.98412698295895385996e-4;
  const gdouble S4    = 8.33333333332211858878e-3;
  const gdouble S5    = -1.66666666666666307295e-1;
  const gdouble C0    = -1.13585365213876817300e-11;
  const gdouble C1    = 2.08757008419747316778e-9;
  const gdouble C2    = -2.75573141792967388112e-7;
  const gdouble C3    = 2.48015872888517045348e-5;
  const gdouble C4    = -1.38888888888730564116e-3;
  const gdouble C5    = 4.16666666666665929218e-2;
  const gdouble tn    = x * TWOPI + ROUND;
  const gdouble n     = tn - ROUND;
  const gdouble r     = ((x - n * DP1) - n * DP2) - n * DP3;
  const gdouble z     = r * r;
  const gdouble sr    = r + r * z * (((((S0 * z + S1) * z + S2) * z + S3) * z + S4) * z + S5);
  const gdouble cr    = (1.0 - 0.5 * z) + z * z * (((((C0 * z + C1) * z + C2) * z + C3) * z + C4) * z + C5);
  guint64 q, us, uc, ns, nc;
  gdouble sv, cv;

  /* n mod 4 is in the last two bits of the mantissa of tn (1.5 2^52 is a multiple of 4) */
  memcpy (&q, &tn, sizeof (q));

  /* Quadrant q: (sin, cos) = (sr, cr), (cr, -sr), (-sr, -cr), (-cr, sr) */
  sv = ncm_vmath_select (q & 1, cr, sr);
  cv = ncm_vmath_select (q & 1, sr, cr);
  ns = (q & 2) << 62;
  nc = ((q + 1) & 2) << 62;

  memcpy (&us, &sv, sizeof (us));
  memcpy (&uc, &cv, sizeof (uc));
  us ^= ns;
  uc ^= nc;
  memcpy (s, &us, sizeof (us));
  memcpy (c, &uc, sizeof (uc));
}

#if defined (__AVX2__) || defined (__AVX512F__)
#define NCM_VMATH_HAVE_VECTOR_EXP_LOG 1
#endif

#ifdef NCM_VMATH_HAVE_VECTOR_EXP_LOG

/* 2^k for an integer k in [-1022, 1023]: k + 1023 is in the low bits of the mantissa of
 * 1.5 2^52 + k + 1023, and is shifted to the exponent field. */
static inline gdouble
_ncm_vmath_pow2i (const gdouble k)
{
  const gdouble t = k + (6755399441055744.0 + 1023.0);
  guint64 u;
  gdouble r;

  memcpy (&u, &t, sizeof (u));
  u <<= 52;
  memcpy (&r, &u, sizeof (r));

  return r;
}

#endif /* NCM_VMATH_HAVE_VECTOR_EXP_LOG */

/*
 * ncm_vmath_exp:
 * @x: a value
 *
 * The exponential of @x. With NCM_VMATH_HAVE_VECTOR_EXP_LOG it is the Pade form of Cephes
 * (Moshier) on $|r| \le \ln 2 / 2$, with $x = n \ln 2 + r$, n rounded by adding and
 * subtracting $1.5 \times 2^{52}$, and $2^n$ applied in two halves so that results below
 * the smallest normal are rounded once, as subnormals. It is within 2 ulp of exp(), 0 and
 * +inf beyond the range of double, and NaN for NaN. Otherwise it is exp().
 *
 * Returns: exp (@x).
 */
static inline gdouble
ncm_vmath_exp (const gdouble x)
{
#ifdef NCM_VMATH_HAVE_VECTOR_EXP_LOG
  const gdouble ROUND = 6755399441055744.0; /* 1.5 2^52: x + ROUND - ROUND rounds x to an integer */
  const gdouble LOG2E = 1.4426950408889634073599;
  const gdouble C1    = 6.93145751953125e-1;
  const gdouble C2    = 1.42860682030941723212e-6;
  const gdouble P0    = 1.26177193074810590878e-4;
  const gdouble P1    = 3.02994407707441961300e-2;
  const gdouble P2    = 9.99999999999999999910e-1;
  const gdouble Q0    = 3.00198505138664455042e-6;
  const gdouble Q1    = 2.52448340349684104192e-3;
  const gdouble Q2    = 2.27265548208155028766e-1;
  const gdouble Q3    = 2.00000000000000000009e0;
  const gdouble xl    = ncm_vmath_select (x < -746.0, -746.0, x);
  const gdouble xc    = ncm_vmath_select (x > 710.0, 710.0, xl);
  const gdouble n     = (xc * LOG2E + ROUND) - ROUND;
  const gdouble r     = (xc - n * C1) - n * C2;
  const gdouble rr    = r * r;
  const gdouble px    = r * ((P0 * rr + P1) * rr + P2);
  const gdouble qx    = ((Q0 * rr + Q1) * rr + Q2) * rr + Q3;
  const gdouble e     = 1.0 + 2.0 * (px / (qx - px));
  const gdouble h     = (0.5 * n + ROUND) - ROUND;

  return (e * _ncm_vmath_pow2i (h)) * _ncm_vmath_pow2i (n - h);

#else /* NCM_VMATH_HAVE_VECTOR_EXP_LOG */

  return exp (x);

#endif /* NCM_VMATH_HAVE_VECTOR_EXP_LOG */
}

/*
 * ncm_vmath_log:
 * @x: a value
 *
 * The natural logarithm of @x. With NCM_VMATH_HAVE_VECTOR_EXP_LOG it is the rational
 * approximation of Cephes (Moshier) on [sqrt(1/2) - 1, sqrt(2) - 1], with the exponent
 * and the mantissa of @x taken from its bits (the exponent turned into a double through
 * the mantissa of 2^52, since AVX2 has no conversion of 64-bit integers) and subnormals
 * scaled by 2^54 first. It is within 1 ulp of log(), -inf at 0, +inf at +inf, and NaN for
 * a negative @x or NaN. Otherwise it is log().
 *
 * Returns: log (@x).
 */
static inline gdouble
ncm_vmath_log (const gdouble x)
{
#ifdef NCM_VMATH_HAVE_VECTOR_EXP_LOG
  const gdouble P0    = 1.01875663804580931796e-4;
  const gdouble P1    = 4.97494994976747001425e-1;
  const gdouble P2    = 4.70579119878881725854e0;
  const gdouble P3    = 1.44989225341610930846e1;
  const gdouble P4    = 1.79368678507819816313e1;
  const gdouble P5    = 7.70838733755885391666e0;
  const gdouble Q0    = 1.12873587189167450590e1;
  const gdouble Q1    = 4.52279145837532221105e1;
  const gdouble Q2    = 8.29875266912776603211e1;
  const gdouble Q3    = 7.11544750618563894466e1;
  const gdouble Q4    = 2.31251620126765340583e1;
  const gdouble SQRTH = 0.70710678118654752440;
  const gboolean sub  = (x < DBL_MIN);
  const gdouble xs    = ncm_vmath_select (sub, x * 18014398509481984.0, x);
  guint64 u, eb, mb;
  gdouble te, m;

  memcpy (&u, &xs, sizeof (u));
  eb = 0x4330000000000000ULL | ((u >> 52) & 0x7ffULL);
  mb = (u & 0x000fffffffffffffULL) | 0x3fe0000000000000ULL;
  memcpy (&te, &eb, sizeof (te));
  memcpy (&m, &mb, sizeof (m));

  {
    /* x = 2^e m, m in [1/2, 1), moved to [sqrt(1/2), sqrt(2)) */
    const gdouble e0 = (te - 4503599627370496.0) - ncm_vmath_select (sub, 1022.0 + 54.0, 1022.0);
    const gboolean c = (m < SQRTH);
    const gdouble e  = ncm_vmath_select (c, e0 - 1.0, e0);
    const gdouble f  = ncm_vmath_select (c, (m + m) - 1.0, m - 1.0);
    const gdouble z  = f * f;
    const gdouble p  = ((((P0 * f + P1) * f + P2) * f + P3) * f + P4) * f + P5;
    const gdouble q  = ((((f + Q0) * f + Q1) * f + Q2) * f + Q3) * f + Q4;
    const gdouble y  = (f * (z * p / q) - e * 2.121944400546905827679e-4) - 0.5 * z;
    gdouble r        = (f + y) + e * 0.693359375;

    r = ncm_vmath_select (x == 0.0, -INFINITY, r);
    r = ncm_vmath_select (x < 0.0, NAN, r);
    r = ncm_vmath_select (x == INFINITY, INFINITY, r);
    r = ncm_vmath_select (x != x, x, r);

    return r;
  }

#else /* NCM_VMATH_HAVE_VECTOR_EXP_LOG */

  return log (x);

#endif /* NCM_VMATH_HAVE_VECTOR_EXP_LOG */
}

/*
 * ncm_vmath_pow_pos:
 * @x: a positive value
 * @y: an exponent
 *
 * @x to the power @y for @x > 0. With NCM_VMATH_HAVE_VECTOR_EXP_LOG it is
 * ncm_vmath_exp (@y ncm_vmath_log (@x)): the error of the logarithm is multiplied by @y,
 * so the relative error is within about (2 + |@y ln @x|) ulp. Otherwise it is pow(). On
 * both paths it is NaN for @x <= 0 or NaN.
 *
 * Returns: @x^@y.
 */
static inline gdouble
ncm_vmath_pow_pos (const gdouble x, const gdouble y)
{
#ifdef NCM_VMATH_HAVE_VECTOR_EXP_LOG

  return ncm_vmath_select (x > 0.0, ncm_vmath_exp (y * ncm_vmath_log (x)), NAN);

#else /* NCM_VMATH_HAVE_VECTOR_EXP_LOG */

  return (x > 0.0) ? pow (x, y) : NAN;

#endif /* NCM_VMATH_HAVE_VECTOR_EXP_LOG */
}

G_END_DECLS

#endif /* _NCM_VMATH_PRIVATE_H_ */

