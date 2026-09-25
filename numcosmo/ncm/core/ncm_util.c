/***************************************************************************
 *            ncm_util.c
 *
 *  Tue Jun  5 00:21:11 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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

/**
 * NcmUtil:
 *
 * Utility functions and macros: elementary functions evaluated without
 * cancellation, Gaussian integrals, floating-point comparison, the #NcmComplex boxed
 * type, on-sky geometry, SUNDIALS CVODE helpers, #GError helpers, test assertions and
 * callback declarations.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_util.h"
#include "ncm/core/ncm_memory_pool.h"
#include "numcosmo/nc/background/nc_hicosmo.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_lambert.h>

#include <cvode/cvode.h>
#include <fftw3.h>
#ifdef HAVE_CFITSIO
#include <fitsio.h>
#endif /* HAVE_CFITSIO */
#if _POSIX_C_SOURCE >= 199309L
#include <time.h> /* for nanosleep */
#else
#include <unistd.h> /* for usleep */
#endif
#include <complex.h>
#endif /* NUMCOSMO_GIR_SCAN */

typedef struct
{
  mpz_t him1, him2, kim1, kim2, hi, ki, a;
} NcmCoarseDbl;

static gpointer
_besselj_bs_alloc (gpointer userdata)
{
  NcmCoarseDbl *cdbl = g_slice_new (NcmCoarseDbl);

  NCM_UNUSED (userdata);
  mpz_inits (cdbl->him1, cdbl->him2, cdbl->kim1, cdbl->kim2, cdbl->hi, cdbl->ki, cdbl->a, NULL);

  return cdbl;
}

static void
_besselj_bs_free (gpointer p)
{
  NcmCoarseDbl *cdbl = (NcmCoarseDbl *) p;

  mpz_clears (cdbl->him1, cdbl->him2, cdbl->kim1, cdbl->kim2, cdbl->hi, cdbl->ki, cdbl->a, NULL);

  g_slice_free (NcmCoarseDbl, cdbl);
}

NcmCoarseDbl **
_ncm_coarse_dbl_get_bs (void)
{
  G_LOCK_DEFINE_STATIC (create_lock);

  static NcmMemoryPool *mp = NULL;

  G_LOCK (create_lock);

  if (mp == NULL)
    mp = ncm_memory_pool_new (_besselj_bs_alloc, NULL, _besselj_bs_free);

  G_UNLOCK (create_lock);

  return ncm_memory_pool_get (mp);
}

/**
 * ncm_rational_coarse_double: (skip)
 * @x: a double
 * @q: the result
 *
 * Sets @q to the continued-fraction approximation of @x truncated when the next
 * correction is below the double precision. Aborts if the relative difference between
 * @q and @x exceeds $10^{-15}$.
 */
void
ncm_rational_coarse_double (gdouble x, mpq_t q)
{
  NcmCoarseDbl **cdbl_ptr = _ncm_coarse_dbl_get_bs ();
  NcmCoarseDbl *cdbl      = *cdbl_ptr;
  gint expo2              = 0;
  gdouble xo              = fabs (frexp (x, &expo2));
  gdouble xi              = xo;

  if (x == 0)
  {
    mpq_set_ui (q, 0, 1);
    ncm_memory_pool_return (cdbl_ptr);

    return;
  }

#define him1 cdbl->him1
#define him2 cdbl->him2
#define kim1 cdbl->kim1
#define kim2 cdbl->kim2
#define hi cdbl->hi
#define ki cdbl->ki
#define a cdbl->a

  mpz_set_ui (him1, 1);
  mpz_set_ui (him2, 0);
  mpz_set_ui (kim1, 0);
  mpz_set_ui (kim2, 1);
  mpz_set_ui (a, 1);

  while (TRUE)
  {
    mpz_set_d (a, xo);
    xo = 1.0 / (xo - mpz_get_d (a));
    mpz_mul (hi, him1, a);
    mpz_add (hi, hi, him2);
    mpz_mul (ki, kim1, a);
    mpz_add (ki, ki, kim2);

    mpz_swap (him2, him1);
    mpz_swap (him1, hi);
    mpz_swap (kim2, kim1);
    mpz_swap (kim1, ki);

    mpz_mul (ki, kim1, kim2);

    if (fabs (1.0 / (xi * mpz_get_d (ki))) < GSL_DBL_EPSILON)
    {
      mpz_set (mpq_numref (q), him2);
      mpz_set (mpq_denref (q), kim2);
      break;
    }

    if (!gsl_finite (xo))
    {
      mpz_set (mpq_numref (q), him1);
      mpz_set (mpq_denref (q), kim1);
      break;
    }
  }

  mpq_canonicalize (q);
  expo2 > 0 ? mpq_mul_2exp (q, q, expo2) : mpq_div_2exp (q, q, -expo2);

  if (GSL_SIGN (x) == -1)
    mpq_neg (q, q);

  if (fabs (mpq_get_d (q) / x - 1) > 1e-15)
  {
    mpfr_fprintf (stderr, "# Q = %Qd\n", q);
    g_error ("Wrong rational approximation for x = %.16g N(q) = %.16g [%.5e] | 2^(%d).", x, mpq_get_d (q), fabs (mpq_get_d (q) / x - 1), expo2);
  }

  ncm_memory_pool_return (cdbl_ptr);

#undef him1
#undef him2
#undef kim1
#undef kim2
#undef hi
#undef ki
#undef a

  return;
}

/**
 * ncm_mpz_inits: (skip)
 * @z: a #mpz_t
 * @...: a %NULL-terminated list of #mpz_t
 *
 * Initializes @z and every #mpz_t in the list.
 */
void
ncm_mpz_inits (mpz_t z, ...)
{
  va_list ap;
  mpz_ptr z1;

  mpz_init (z);
  va_start (ap, z);

  while ((z1 = va_arg (ap, mpz_ptr)) != NULL)
    mpz_init (z1);

  va_end (ap);
}

/**
 * ncm_mpz_clears: (skip)
 * @z: a #mpz_t
 * @...: a %NULL-terminated list of #mpz_t
 *
 * Clears @z and every #mpz_t in the list.
 */
void
ncm_mpz_clears (mpz_t z, ...)
{
  va_list ap;
  mpz_ptr z1;

  mpz_clear (z);
  va_start (ap, z);

  while ((z1 = va_arg (ap, mpz_ptr)) != NULL)
    mpz_clear (z1);

  va_end (ap);
}

/**
 * ncm_util_sqrt1px_m1:
 * @x: a real number $x > -1$
 *
 * Computes $\sqrt{1+x} - 1$ as $x / (\sqrt{1+x} + 1)$, without cancellation at
 * $x \approx 0$.
 *
 * Returns: $\sqrt{1+x} - 1$.
 */
/**
 * ncm_util_ln1pexpx:
 * @x: a real number $x$
 *
 * Computes $\ln(1 + e^x)$ without overflow; it returns $x$ when $e^{-x}$ is below the
 * double precision.
 *
 * Returns: $\ln(1 + e^x)$.
 */
/**
 * ncm_util_1pcosx:
 * @sinx: $\sin x$
 * @cosx: $\cos x$
 *
 * Computes $1 + \cos x$, as $\sin^2 x / (1 - \cos x)$ when $\cos x \leq -0.9$.
 *
 * Returns: $1 + \cos x$.
 */
/**
 * ncm_util_1mcosx:
 * @sinx: $\sin x$
 * @cosx: $\cos x$
 *
 * Computes $1 - \cos x$, as $\sin^2 x / (1 + \cos x)$ when $\cos x \geq 0.9$.
 *
 * Returns: $1 - \cos x$.
 */
/**
 * ncm_util_1psinx:
 * @sinx: $\sin x$
 * @cosx: $\cos x$
 *
 * Computes $1 + \sin x$, as $\cos^2 x / (1 - \sin x)$ when $\sin x \leq -0.9$.
 *
 * Returns: $1 + \sin x$.
 */
/**
 * ncm_util_1msinx:
 * @sinx: $\sin x$
 * @cosx: $\cos x$
 *
 * Computes $1 - \sin x$, as $\cos^2 x / (1 + \sin x)$ when $\sin x \geq 0.9$.
 *
 * Returns: $1 - \sin x$.
 */
/**
 * ncm_util_cos2x:
 * @sinx: $\sin x$
 * @cosx: $\cos x$
 *
 * Computes $\cos(2x)$ as $(\cos x - \sin x)(\cos x + \sin x)$.
 *
 * Returns: $\cos(2x)$.
 */

/**
 * ncm_cmpdbl:
 * @x: a double
 * @y: a double
 *
 * Returns: $|2(x - y)/(x + y)|$, or zero if $x = y$.
 */
gdouble
ncm_cmpdbl (const gdouble x, const gdouble y)
{
  if (x == y)
    return 0.0;
  else
    return fabs (2.0 * (x - y) / (x + y));
}

/**
 * ncm_exprel:
 * @x: a double
 *
 * Returns: $(e^x - 1)/x$.
 */
gdouble
ncm_exprel (const gdouble x)
{
  return gsl_sf_hyperg_1F1_int (1, 2, x);
}

/**
 * ncm_d1exprel:
 * @x: a double
 *
 * Returns: the first derivative of $(e^x - 1)/x$.
 */
gdouble
ncm_d1exprel (const gdouble x)
{
  return 0.5 * gsl_sf_hyperg_1F1_int (2, 3, x);
}

/**
 * ncm_d2exprel:
 * @x: a double
 *
 * Returns: the second derivative of $(e^x - 1)/x$.
 */
gdouble
ncm_d2exprel (const gdouble x)
{
  return gsl_sf_hyperg_1F1_int (3, 4, x) / 3.0;
}

/**
 * ncm_d3exprel:
 * @x: a double
 *
 * Returns: the third derivative of $(e^x - 1)/x$.
 */
gdouble
ncm_d3exprel (const gdouble x)
{
  return 0.25 * gsl_sf_hyperg_1F1_int (4, 5, x);
}

/**
 * ncm_util_sinh1:
 * @x: a double
 *
 * Computes $\sinh(x)/x$, from its Taylor series for $|x| < 0.9$.
 *
 * Returns: $\sinh(x)/x$.
 */
gdouble
ncm_util_sinh1 (const gdouble x)
{
  const gdouble cut = 9.0 / 10.0;

  if (fabs (x) < cut)
  {
    const gdouble x2 = x * x;
    gdouble d        = 1.0;
    gdouble x2n      = 1.0;
    gdouble p        = 2.0;
    gdouble res      = 0.0;
    gint n;

    for (n = 0; n < 9; n++)
    {
      res += x2n / d;

      d   *= (p * (p + 1.0));
      x2n *= x2;
      p   += 2.0;
    }

    return res;
  }
  else
  {
    return sinh (x) / x;
  }
}

/**
 * ncm_util_sinh3:
 * @x: a double
 *
 * Computes $[\sinh(x) - x] / (x^3/3!)$, from its Taylor series for $|x| < 0.9$.
 *
 * Returns: $[\sinh(x) - x] / (x^3/3!)$.
 */
gdouble
ncm_util_sinh3 (const gdouble x)
{
  const gdouble cut = 9.0 / 10.0;

  if (fabs (x) < cut)
  {
    const gdouble x2 = x * x;
    gdouble d        = 1.0;
    gdouble x2n      = 1.0;
    gdouble p        = 4.0;
    gdouble res      = 0.0;
    gint n;

    for (n = 0; n < 8; n++)
    {
      res += x2n / d;

      d   *= (p * (p + 1.0));
      x2n *= x2;
      p   += 2.0;
    }

    return res;
  }
  else
  {
    return (sinh (x) - x) * 6.0 / gsl_pow_3 (x);
  }
}

/**
 * ncm_util_sinhx_m_xcoshx_x3:
 * @x: a double
 *
 * Computes $[\sinh(x) - x\cosh(x)]/x^3$, without cancellation for $|x| < 0.9$.
 *
 * Returns: $[\sinh(x) - x\cosh(x)]/x^3$.
 */
gdouble
ncm_util_sinhx_m_xcoshx_x3 (const gdouble x)
{
  const gdouble cut = 9.0 / 10.0;
  const gdouble shx = sinh (x);

  if (fabs (x) < cut)
    return ncm_util_sinh3 (x) / 6.0 - gsl_pow_2 (ncm_util_sinh1 (x)) / (sqrt (1.0 + shx * shx) + 1.0);
  else
    return (shx - x * sqrt (1.0 + shx * shx)) / gsl_pow_3 (x);
}

/**
 * ncm_util_lambert_W0_ln:
 * @ln_y: the logarithm $\ln y$
 *
 * Computes the principal branch $W_0(y)$ of the Lambert $W$ function from $\ln y$,
 * so that $y$ may exceed the double range. Uses gsl_sf_lambert_W0() for
 * $\ln y < \ln(\mathrm{DBL\_MAX}) - 1$, otherwise Newton's method on
 * $W + \ln W = \ln y$. Aborts if Newton's method does not converge.
 *
 * Returns: $W_0(y)$.
 */
gdouble
ncm_util_lambert_W0_ln (const gdouble ln_y)
{
  /* gsl_sf_lambert_W0() fails as y approaches DBL_MAX */
  if (ln_y < GSL_LOG_DBL_MAX - 1.0)
  {
    return gsl_sf_lambert_W0 (exp (ln_y));
  }
  else
  {
    const guint max_iter = 100;
    gdouble W            = ln_y - log (ln_y);
    guint i;

    for (i = 0; i < max_iter; i++)
    {
      const gdouble dW = (W + log (W) - ln_y) / (1.0 + 1.0 / W);

      W -= dW;

      if (fabs (dW) <= GSL_DBL_EPSILON * W)
        return W;
    }

    g_error ("ncm_util_lambert_W0_ln: Newton's method did not converge for ln_y = %g.", ln_y);

    return GSL_NAN;
  }
}

/**
 * ncm_util_mln_1mIexpzA_1pIexpmzA:
 * @rho: $\rho$
 * @theta: $\theta$
 * @A: $A$
 * @rho1: (out): $\rho_1$
 * @theta1: (out): $\theta_1$
 *
 * Computes
 * $$z_1 = z - \ln\left(\frac{1 - i A e^{z}}{1 + i A e^{-z}}\right), \qquad z = \rho + i\theta,$$
 * from the series of the logarithm in $A$ when $e^{|\rho|}|A| < 0.1$, and sets
 * $z_1 = \rho_1 + i\theta_1$.
 */
void
ncm_util_mln_1mIexpzA_1pIexpmzA (const gdouble rho, const gdouble theta, const gdouble A, gdouble *rho1, gdouble *theta1)
{
  const double complex z = rho + I * theta;
  double complex zp;

  if (exp (fabs (rho)) * fabs (A) < 0.1)
  {
    const double complex z_p_ipi_2 = z + 0.5 * M_PI * I;
    const double complex T         = cexp (z_p_ipi_2);
    double complex Tn              = T;
    gdouble An                     = A;
    gint i;

    zp = 0.0;

    for (i = 0; ; i++)
    {
      const gdouble n   = i + 1.0;
      double complex dz = (1.0 / Tn - Tn) * An / n;

      An *= A;
      Tn *= T;

      zp += dz;

      if (cabs (dz / zp) < GSL_DBL_EPSILON * 1.0e-0)
        break;
    }

    zp = z - zp;
  }
  else
  {
    zp = z - clog ((1.0 - I * A * cexp (z)) / (1.0 + I * A * cexp (-z)));
  }

  rho1[0]   = creal (zp);
  theta1[0] = cimag (zp);
}

#define ERF_BOUND (3.0)

/**
 * ncm_util_normal_gaussian_integral:
 * @xl: lower limit
 * @xu: upper limit
 *
 * Computes $\int_{x_l}^{x_u} e^{-x^2/2}\,\mathrm{d}x / \sqrt{2\pi}$, using erfc() when both
 * limits are in the same tail, beyond $3\sqrt{2}$, to keep the relative precision. The
 * result is negative for $x_u < x_l$.
 *
 * Returns: the integral.
 */
gdouble
ncm_util_normal_gaussian_integral (const gdouble xl, const gdouble xu)
{
  if (xl == xu)
    return 0.0;

  {
    const gdouble sqrt_half = M_SQRT1_2;
    gdouble ul, uu, sign;

    if (xl < xu)
    {
      ul   = xl * sqrt_half;
      uu   = xu * sqrt_half;
      sign = +1.0;
    }
    else
    {
      ul   = xu * sqrt_half;
      uu   = xl * sqrt_half;
      sign = -1.0;
    }

    if (ul > ERF_BOUND)
    {
      /* Both limits in the right tail */
      const gdouble val = 0.5 * (erfc (ul) - erfc (uu));

      return sign * val;
    }
    else if (uu < -ERF_BOUND)
    {
      /* Both limits in the left tail */
      const gdouble val = 0.5 * (erfc (-uu) - erfc (-ul));

      return sign * val;
    }
    else
    {
      const gdouble val = 0.5 * (erf (uu) - erf (ul));

      return sign * val;
    }
  }
}

/**
 * ncm_util_gaussian_integral:
 * @xl: lower limit
 * @xu: upper limit
 * @mu: mean
 * @sigma: standard deviation
 *
 * Same as ncm_util_normal_gaussian_integral() for the Gaussian of mean @mu and standard
 * deviation @sigma.
 *
 * Returns: the integral.
 */
gdouble
ncm_util_gaussian_integral (const gdouble xl, const gdouble xu, const gdouble mu, const gdouble sigma)
{
  return ncm_util_normal_gaussian_integral ((xl - mu) / sigma, (xu - mu) / sigma);
}

/**
 * ncm_util_log_normal_gaussian_integral:
 * @xl: lower limit
 * @xu: upper limit
 * @sign: (out): the sign of the integral
 *
 * Computes the logarithm of the absolute value of ncm_util_normal_gaussian_integral(),
 * keeping the relative precision in both tails and when the integral is close to one.
 *
 * Returns: the logarithm of the absolute value of the integral, $-\infty$ if
 * $x_l = x_u$.
 */
gdouble
ncm_util_log_normal_gaussian_integral (const gdouble xl, const gdouble xu, gdouble *sign)
{
  if (xl == xu)
    return -INFINITY;

  {
    const gdouble sqrt_half = M_SQRT1_2;
    gdouble ul, uu;

    if (xl < xu)
    {
      ul    = xl * sqrt_half;
      uu    = xu * sqrt_half;
      *sign = 1.0;
    }
    else
    {
      ul    = xu * sqrt_half;
      uu    = xl * sqrt_half;
      *sign = -1.0;
    }

    if (ul > ERF_BOUND)
    {
      const gdouble val = 0.5 * (erfc (ul) - erfc (uu));

      return log (fabs (val));
    }
    else if (uu < -ERF_BOUND)
    {
      const gdouble val = 0.5 * (erfc (-ul) - erfc (-uu));

      return log (fabs (val));
    }
    else if ((uu > ERF_BOUND) && (ul < -ERF_BOUND))
    {
      const gdouble val = -0.5 * (erfc (uu) + erfc (-ul));

      return log1p (val);
    }
    else
    {
      const gdouble val = 0.5 * (erf (uu) - erf (ul));

      return log (fabs (val));
    }
  }
}

/**
 * ncm_util_log_gaussian_integral:
 * @xl: lower limit
 * @xu: upper limit
 * @mu: mean
 * @sigma: standard deviation
 * @sign: (out): the sign of the integral
 *
 * Same as ncm_util_log_normal_gaussian_integral() for the Gaussian of mean @mu and
 * standard deviation @sigma.
 *
 * Returns: the logarithm of the absolute value of the integral.
 */
gdouble
ncm_util_log_gaussian_integral (const gdouble xl, const gdouble xu, const gdouble mu, const gdouble sigma, gdouble *sign)
{
  const gdouble zl = (xl - mu) / sigma;
  const gdouble zu = (xu - mu) / sigma;

  return ncm_util_log_normal_gaussian_integral (zl, zu, sign);
}

/**
 * ncm_cmp:
 * @x: a double
 * @y: a double
 * @reltol: relative tolerance
 * @abstol: absolute tolerance
 *
 * Compares @x and @y, which are equal when
 * $|x - y| \leq \epsilon_\mathrm{rel} \max(|x|, |y|) + \epsilon_\mathrm{abs}$. If one of them is
 * zero, $\max(|x|, |y|)$ is replaced by one, so @reltol also acts as an absolute tolerance.
 *
 * Returns: $-1$ if $x < y$, $0$ if they are equal, $1$ if $x > y$.
 */
gint
ncm_cmp (gdouble x, gdouble y, const gdouble reltol, const gdouble abstol)
{
  if (G_UNLIKELY ((x == 0.0) && (y == 0.0)))
  {
    return 0;
  }
  else
  {
    const gdouble delta = (x - y);
    const gdouble abs_x = fabs (x);
    const gdouble abs_y = fabs (y);
    const gdouble mean  = G_UNLIKELY (x == 0.0 || y == 0.0) ? 1.0 : GSL_MAX (abs_x, abs_y);

    if (fabs (delta) <= reltol * mean + abstol)
      return 0;
    else
      return delta < 0 ? -1 : 1;
  }
}

void
_ncm_assertion_message_cmpdouble (const gchar *domain, const gchar *file, gint line, const gchar *func, const gchar *expr, gdouble arg1, const gchar *cmp, gdouble arg2, const gdouble reltol, const gdouble abstol)
{
  gchar *s = g_strdup_printf ("assertion failed (%s): (%.17g %s %.17g) (reltol %.17g diff_rel %.17g, abstol %.17g diff %.17g)",
                              expr, arg1, cmp, arg2, reltol,
                              fabs (arg1) > fabs (arg2) ? fabs (arg2 / arg1 - 1.0) : fabs (arg1 / arg2 - 1.0),
                              abstol,
                              fabs (arg1 - arg2));

  g_assertion_message (domain, file, line, func, s);
  g_free (s);
}

G_DEFINE_BOXED_TYPE (NcmComplex, ncm_complex, ncm_complex_dup, ncm_complex_free)

/**
 * ncm_complex_new:
 *
 * Allocates a complex number set to zero.
 *
 * Returns: (transfer full): a new #NcmComplex.
 */
NcmComplex *
ncm_complex_new ()
{
  return g_new0 (NcmComplex, 1);
}

/**
 * ncm_complex_dup:
 * @c: a #NcmComplex
 *
 * Returns: (transfer full): a newly allocated copy of @c.
 */
NcmComplex *
ncm_complex_dup (NcmComplex *c)
{
  NcmComplex *cc = ncm_complex_new ();

  *cc = *c;

  return cc;
}

/**
 * ncm_complex_free:
 * @c: a #NcmComplex
 *
 * Frees @c, which must come from ncm_complex_new() or ncm_complex_dup().
 */
void
ncm_complex_free (NcmComplex *c)
{
  g_free (c);
}

/**
 * ncm_complex_clear:
 * @c: a #NcmComplex
 *
 * Frees *@c, as ncm_complex_free(), and sets *@c to %NULL.
 */
void
ncm_complex_clear (NcmComplex **c)
{
  g_clear_pointer (c, g_free);
}

/**
 * ncm_complex_set:
 * @c: a #NcmComplex
 * @a: the real part $a$
 * @b: the imaginary part $b$
 *
 * Sets @c to $a + i b$.
 */
/**
 * ncm_complex_set_c: (skip)
 * @c: a #NcmComplex
 * @z: a complex double
 *
 * Sets @c to @z.
 */
/**
 * ncm_complex_set_zero:
 * @c: a #NcmComplex
 *
 * Sets @c to zero.
 */
/**
 * ncm_complex_Re:
 * @c: a #NcmComplex
 *
 * Returns: the real part of @c.
 */
/**
 * ncm_complex_Im:
 * @c: a #NcmComplex
 *
 * Returns: the imaginary part of @c.
 */
/**
 * ncm_complex_Abs:
 * @c: a #NcmComplex
 *
 * Returns: $|c|$.
 */
/**
 * ncm_complex_c: (skip)
 * @c: a #NcmComplex
 *
 * Returns: @c as a complex double.
 */

/**
 * ncm_complex_res_add_mul_real:
 * @c1: a #NcmComplex
 * @c2: a #NcmComplex
 * @v: a double
 *
 * Sets $c_1 \to c_1 + c_2 v$. @c1 and @c2 must not overlap.
 */
/**
 * ncm_complex_res_add_mul:
 * @c1: a #NcmComplex
 * @c2: a #NcmComplex
 * @c3: a #NcmComplex
 *
 * Sets $c_1 \to c_1 + c_2 c_3$. @c1 must not overlap @c2 or @c3.
 */

/**
 * ncm_complex_mul_real:
 * @c: a #NcmComplex
 * @v: a double
 *
 * Sets $c \to c v$.
 */
/**
 * ncm_complex_res_mul:
 * @c1: a #NcmComplex
 * @c2: a #NcmComplex
 *
 * Sets $c_1 \to c_1 c_2$. @c1 and @c2 must not overlap.
 */

/**
 * ncm_util_position_angle:
 * @ra1: right ascension of the first object, in degrees
 * @dec1: declination of the first object, in degrees
 * @ra2: right ascension of the second object, in degrees
 * @dec2: declination of the second object, in degrees
 *
 * Returns: the position angle of the second object seen from the first, East of
 * North, in radians.
 */

/**
 * ncm_util_great_circle_distance:
 * @ra1: right ascension of the first object, in degrees
 * @dec1: declination of the first object, in degrees
 * @ra2: right ascension of the second object, in degrees
 * @dec2: declination of the second object, in degrees
 *
 * Computes the angular separation from the Vincenty formula, see
 * [great-circle distance](https://en.wikipedia.org/wiki/Great-circle_distance).
 *
 * Returns: the angular separation, in degrees.
 */

/**
 * ncm_util_projected_radius:
 * @theta: angular separation, in radians
 * @d: distance
 *
 * Returns: $d \sin\theta$, the separation projected at distance @d, in the units of @d.
 */

/**
 * ncm_util_smooth_trans:
 * @f0: value before the transition
 * @f1: value after the transition
 * @z0: start of the transition
 * @dz: length of the transition
 * @z: the point
 *
 * Computes $f = f_0 \theta_0 + f_1 \theta_1$ with the logistic weights
 * $\theta_0 = 1/(1 + e^{g})$ and $\theta_1 = 1/(1 + e^{-g})$, where
 * $g = 72 (z - z_0 - \Delta z/2)/\Delta z$. At $z_0$ and $z_0 + \Delta z$ the weight of the
 * other value is $e^{-36} \approx 2 \times 10^{-16}$.
 *
 * Returns: $f$.
 */
/**
 * ncm_util_smooth_trans_get_theta:
 * @z0: start of the transition
 * @dz: length of the transition
 * @z: the point
 * @theta0: (out): the weight $\theta_0$
 * @theta1: (out): the weight $\theta_1$
 *
 * Computes the weights of ncm_util_smooth_trans().
 */

/**
 * ncm_util_cvode_check_flag:
 * @flagvalue: the value returned by a SUNDIALS function
 * @funcname: the name of that function
 * @opt: the kind of value
 *
 * Checks the value returned by the SUNDIALS function @funcname and logs a message on
 * failure. For @opt 0, @flagvalue is a pointer that must not be %NULL; for 1, a pointer
 * to an int flag that must not be negative; for 2, a pointer returned by a memory
 * allocation that must not be %NULL. Aborts for other values of @opt.
 *
 * Returns: whether the call succeeded.
 */
gboolean
ncm_util_cvode_check_flag (gpointer flagvalue, const gchar *funcname, gint opt)
{
  gint *errflag;

  switch (opt)
  {
    case 0:
    {
      if (flagvalue == NULL)
      {
        g_message ("\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);

        return FALSE;
      }

      break;
    }
    case 1:
    {
      errflag = (int *) flagvalue;

      if (*errflag < 0)
      {
        g_message ("\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);

        return FALSE;
      }

      break;
    }
    case 2:
    {
      if (flagvalue == NULL)
      {
        g_message ("\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);

        return FALSE;
      }

      break;
    }
    default:
      g_assert_not_reached ();
  }

  return TRUE;
}

/**
 * ncm_util_cvode_print_stats:
 * @cvode: a CVODE memory block
 *
 * Logs the integrator statistics of @cvode.
 *
 * Returns: %TRUE.
 */
gboolean
ncm_util_cvode_print_stats (gpointer cvode)
{
  glong nsteps, nfunceval, nlinsetups, njaceval, ndiffjaceval, nnonliniter,
        nconvfail, nerrortests, nrooteval;
  gint flag, qcurorder, qlastorder;
  gdouble hinused, hlast, hcur, tcur;

  flag = CVodeGetIntegratorStats (cvode, &nsteps, &nfunceval,
                                  &nlinsetups, &nerrortests, &qlastorder, &qcurorder,
                                  &hinused, &hlast, &hcur, &tcur);
  ncm_util_cvode_check_flag (&flag, "CVodeGetIntegratorStats", 1);


  flag = CVodeGetNumNonlinSolvIters (cvode, &nnonliniter);
  ncm_util_cvode_check_flag (&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails (cvode, &nconvfail);
  ncm_util_cvode_check_flag (&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVodeGetNumJacEvals (cvode, &njaceval);
  ncm_util_cvode_check_flag (&flag, "CVodeGetNumJacEvals", 1);

  flag = CVodeGetNumRhsEvals (cvode, &ndiffjaceval);
  ncm_util_cvode_check_flag (&flag, "CVodeGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals (cvode, &nrooteval);
  ncm_util_cvode_check_flag (&flag, "CVodeGetNumGEvals", 1);

  flag = CVodeGetLastOrder (cvode, &qlastorder);
  ncm_util_cvode_check_flag (&flag, "CVodeGetLastOrder", 1);

  flag = CVodeGetCurrentOrder (cvode, &qcurorder);
  ncm_util_cvode_check_flag (&flag, "CVodeGetCurrentOrder", 1);

  g_message ("# Final Statistics:\n");
  g_message ("# nerrortests = %-6ld | %.5e %.5e %.5e %.5e\n", nerrortests, hinused, hlast, hcur, tcur);
  g_message ("# nsteps = %-6ld nfunceval  = %-6ld nlinsetups = %-6ld njaceval = %-6ld ndiffjaceval = %ld\n",
             nsteps, nfunceval, nlinsetups, njaceval, ndiffjaceval);
  g_message ("# nnonliniter = %-6ld nconvfail = %-6ld nerrortests = %-6ld nrooteval = %-6ld qcurorder = %-6d qlastorder = %-6d\n",
             nnonliniter, nconvfail, nerrortests, nrooteval, qcurorder, qlastorder);

  return TRUE;
}

/**
 * ncm_util_basename_fits:
 * @fits_filename: a file name
 *
 * Removes a final ".fits" or ".fit" extension, in any case, from @fits_filename.
 *
 * Returns: (transfer full): @fits_filename without the extension.
 */
gchar *
ncm_util_basename_fits (const gchar *fits_filename)
{
  GError *error          = NULL;
  GRegex *fits_ext       = g_regex_new ("(.*)\\.[fF][iI][tT][sS]?$", 0, 0, &error);
  GMatchInfo *match_info = NULL;
  gchar *base_name       = NULL;

  if (g_regex_match (fits_ext, fits_filename, 0, &match_info))
    base_name = g_match_info_fetch (match_info, 1);
  else
    base_name = g_strdup (fits_filename);

  g_match_info_free (match_info);
  g_regex_unref (fits_ext);

  return base_name;
}

/**
 * ncm_util_function_params:
 * @func: a string "name" or "name(x1, x2, ...)"
 * @x: (out) (array length=len) (element-type double) (transfer full): the parameters
 * @len: (out): number of parameters
 *
 * Parses the function name and its numerical parameters from @func. Aborts if a
 * parameter is not a number.
 *
 * Returns: (transfer full) (nullable): the function name, or %NULL if @func does not
 * have that form.
 */
gchar *
ncm_util_function_params (const gchar *func, gdouble **x, guint *len)
{
  GError *error          = NULL;
  GRegex *func_params    = g_regex_new ("^([A-Za-z][A-Za-z0-9\\_\\:]*)\\s*(?:\\(\\s*([0-9\\.,eE\\-\\+\\s]*?)\\s*\\))?\\s*$", 0, 0, &error);
  GMatchInfo *match_info = NULL;
  gchar *func_name       = NULL;

  if (error != NULL)
    g_error ("ncm_util_function_params: `%s'.", error->message);

  func_name = NULL;
  *x        = NULL;
  *len      = 0;

  if (g_regex_match (func_params, func, 0, &match_info))
  {
    gchar *fpa = g_match_info_fetch (match_info, 2);

    func_name = g_match_info_fetch (match_info, 1);

    if (fpa != NULL)
    {
      gchar **xs = g_regex_split_simple ("\\s*,\\s*", fpa, 0, 0);

      *len = g_strv_length (xs);

      if (*len > 0)
      {
        guint i;

        *x = g_new0 (gdouble, *len);

        for (i = 0; i < *len; i++)
        {
          gchar *endptr = NULL;

          (*x)[i] = strtod (xs[i], &endptr);

          if ((endptr == xs[i]) || (strlen (endptr) > 0))
            g_error ("ncm_util_function_params: cannot identify double in string `%s'.", xs[i]);
        }
      }

      g_strfreev (xs);
      g_free (fpa);
    }
  }

  g_match_info_free (match_info);
  g_regex_unref (func_params);

  return func_name;
}

/**
 * ncm_util_fact_size:
 * @n: a positive integer
 *
 * Computes an integer $n_f \geq n$ whose prime factors are 2, 3, 5 and 7, a size for
 * which FFTW is efficient. $n_f$ is built by dividing out these factors and
 * incrementing when none divides, so it is not always the smallest such integer.
 *
 * Returns: $n_f$.
 */
gulong
ncm_util_fact_size (const gulong n)
{
  if (n == 1)
  {
    return 0;
  }
  else
  {
    const gulong r2 = n % 2;
    const gulong r3 = n % 3;
    const gulong r5 = n % 5;
    const gulong r7 = n % 7;
    gulong m        = 1;

    if (r2 == 0)
      m *= 2;

    if (r3 == 0)
      m *= 3;

    if (r5 == 0)
      m *= 5;

    if (r7 == 0)
      m *= 7;

    if (m != 1)
    {
      if (n / m == 1)
        return m;
      else
        return m * ncm_util_fact_size (n / m);
    }
    else
    {
      return ncm_util_fact_size (n + 1);
    }
  }
}

void
_ncm_util_set_destroyed (gpointer b)
{
  gboolean *destroyed = b;

  *destroyed = TRUE;
}

/**
 * ncm_util_sleep_ms:
 * @milliseconds: time in milliseconds
 *
 * Suspends the calling thread for @milliseconds.
 */
void
ncm_util_sleep_ms (gint milliseconds)
{
#if _POSIX_C_SOURCE >= 199309L
  struct timespec ts;

  ts.tv_sec  = milliseconds / 1000;
  ts.tv_nsec = (milliseconds % 1000) * 1000000;
  nanosleep (&ts, NULL);
#else

  if (milliseconds >= 1000)
    sleep (milliseconds / 1000);

  usleep ((milliseconds % 1000) * 1000);
#endif
}

/**
 * ncm_util_set_or_call_error:
 * @error: (nullable): a #GError location
 * @domain: the error domain
 * @code: the error code
 * @format: a printf format string
 * @...: arguments for @format
 *
 * Sets *@error to a new error with the formatted message, or aborts with that message
 * if @error is %NULL. Aborts if *@error is already set.
 */
void
ncm_util_set_or_call_error (GError **error, GQuark domain, gint code, const gchar *format, ...)
{
  va_list ap;
  gchar *message;

  va_start (ap, format);
  message = g_strdup_vprintf (format, ap);
  va_end (ap);

  if (error != NULL)
  {
    if (*error != NULL)
      g_error ("ncm_util_set_or_call_error: error piling up: (%s) over (%s)",
               message,
               (*error)->message
      );

    g_set_error (error, domain, code, "%s", message);
  }
  else
  {
    g_error ("%s", message);
  }

  g_free (message);
}

/**
 * ncm_util_forward_or_call_error:
 * @error: (nullable): a #GError location
 * @local_error: (nullable) (transfer full): an error
 * @format: a printf format string
 * @...: arguments for @format
 *
 * Does nothing if @local_error is %NULL. Otherwise moves @local_error to *@error with the
 * formatted message as a prefix, or aborts with both messages if @error is %NULL.
 */
void
ncm_util_forward_or_call_error (GError **error, GError *local_error, const gchar *format, ...)
{
  va_list ap;
  gchar *message;

  /* No error to forward */
  if (local_error == NULL)
    return;

  va_start (ap, format);
  message = g_strdup_vprintf (format, ap);
  va_end (ap);


  if (error != NULL)
    g_propagate_prefixed_error (error, local_error, "%s: ", message);
  else
    g_error ("%s: %s", message, local_error->message);

  g_free (message);
}

