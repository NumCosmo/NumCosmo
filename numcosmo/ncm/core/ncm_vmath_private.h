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
 * omp simd loop, for instance); a call to libm stops the vectorization of a loop when
 * -ffast-math is not used, as in NumCosmo, since glibc declares its vector variants only
 * under __FAST_MATH__. They are for hot loops, and libm is the reference for the rest.
 *
 * ncm_vmath_select() is the blend: GCC keeps a conditional expression on doubles as
 * control flow when floating-point exceptions trap (the default) and the target has no
 * masked instructions (SSE2 and AVX2), which leaves the loop scalar.
 *
 * ncm_vmath_atan2() is faster than libm with SSE2 already. ncm_vmath_exp(),
 * ncm_vmath_log() and ncm_vmath_pow_pos() are faster only with AVX2 or AVX-512 (two-wide
 * division and no fma with SSE2): there NCM_VMATH_HAVE_FAST_EXP is defined, and elsewhere
 * they are exp(), log() and pow().
 */

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

#if defined (__AVX2__) || defined (__AVX512F__)
#define NCM_VMATH_HAVE_FAST_EXP 1
#endif

#ifdef NCM_VMATH_HAVE_FAST_EXP

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

#endif /* NCM_VMATH_HAVE_FAST_EXP */

/*
 * ncm_vmath_exp:
 * @x: a value
 *
 * The exponential of @x. With NCM_VMATH_HAVE_FAST_EXP it is the Pade form of Cephes
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
#ifdef NCM_VMATH_HAVE_FAST_EXP
  const gdouble MAGIC = 6755399441055744.0;
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
  const gdouble n     = (xc * LOG2E + MAGIC) - MAGIC;
  const gdouble r     = (xc - n * C1) - n * C2;
  const gdouble rr    = r * r;
  const gdouble px    = r * ((P0 * rr + P1) * rr + P2);
  const gdouble qx    = ((Q0 * rr + Q1) * rr + Q2) * rr + Q3;
  const gdouble e     = 1.0 + 2.0 * (px / (qx - px));
  const gdouble h     = (0.5 * n + MAGIC) - MAGIC;

  return (e * _ncm_vmath_pow2i (h)) * _ncm_vmath_pow2i (n - h);

#else /* NCM_VMATH_HAVE_FAST_EXP */

  return exp (x);

#endif /* NCM_VMATH_HAVE_FAST_EXP */
}

/*
 * ncm_vmath_log:
 * @x: a value
 *
 * The natural logarithm of @x. With NCM_VMATH_HAVE_FAST_EXP it is the rational
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
#ifdef NCM_VMATH_HAVE_FAST_EXP
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

#else /* NCM_VMATH_HAVE_FAST_EXP */

  return log (x);

#endif /* NCM_VMATH_HAVE_FAST_EXP */
}

/*
 * ncm_vmath_pow_pos:
 * @x: a positive value
 * @y: an exponent
 *
 * @x to the power @y for @x > 0. With NCM_VMATH_HAVE_FAST_EXP it is
 * ncm_vmath_exp (@y ncm_vmath_log (@x)): the error of the logarithm is multiplied by @y,
 * so the relative error is within about (2 + |@y ln @x|) ulp. Otherwise it is pow().
 *
 * Returns: @x^@y.
 */
static inline gdouble
ncm_vmath_pow_pos (const gdouble x, const gdouble y)
{
#ifdef NCM_VMATH_HAVE_FAST_EXP

  return ncm_vmath_exp (y * ncm_vmath_log (x));

#else /* NCM_VMATH_HAVE_FAST_EXP */

  return pow (x, y);

#endif /* NCM_VMATH_HAVE_FAST_EXP */
}

G_END_DECLS

#endif /* _NCM_VMATH_PRIVATE_H_ */

