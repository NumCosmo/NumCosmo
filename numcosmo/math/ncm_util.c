/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

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
 * SECTION:ncm_util
 * @title: NcmUtil
 * @short_description: Miscellaneous utilities.
 *
 * Miscellaneous utility functions, macros and objects.
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_util.h"
#include "math/ncm_memory_pool.h"
#include "numcosmo/nc_hicosmo.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_hyperg.h>

#include <cvode/cvode.h>
#ifdef HAVE_FFTW3
#include <fftw3.h>
#endif /* HAVE_FFTW3 */
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
 * ncm_rational_coarce_double: (skip)
 * @x: a double
 * @q: a #mpq_t to store the result
 *
 * Computes a rational approximation for @x and stores the result in @q.
 *
 */
void
ncm_rational_coarce_double (gdouble x, mpq_t q)
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

  /*  mpq_set_d (q, x); */
  /*  mpq_canonicalize (q); */
  /*  mpfr_printf ("# AUTO: %.15g %.15f %d | %Qd\n", x, xo, expo2, q); */
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

    /*mpfr_printf ("# ---- %.5e %Zu %.15g | %Zd/%Zd %Zd/%Zd\n", fabs(1.0/(xi * mpz_get_d (ki))), a, xo, him1, kim1, him2, kim2); */
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

  /*  mpfr_printf ("# MINE: %.15g %.15f %d | %Qd\n", x, xo, expo2, q); */
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
 * @z: a #mpz_t to initialize
 * @...: a null terminated list of #mpz_t to initialize
 *
 * Initializes @z and all the #mpz_t in the list.
 *
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
 * @z: a #mpz_t to clear
 * @...: a null terminated list of #mpz_t to clear
 *
 * Clears @z and all the #mpz_t in the list.
 *
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
 * @x: a real number $&gt;-1$
 *
 * Calculates $\sqrt{1+x}-1$ using the appropriated expression
 * to avoid round-off when $x \approx 0$.
 *
 * Returns: $\sqrt{1+x}-1$.
 */
/**
 * ncm_util_ln1pexpx:
 * @x: a real number $x$
 *
 * Calculates $\ln[1+\exp(x)]$.
 *
 * Returns: $\ln[1+\exp(x)]$.
 */
/**
 * ncm_util_1pcosx:
 * @sinx: a real number $\sin(x)$
 * @cosx: a real number $\cos(x)$
 *
 * Calculates $1 + \cos(x)$ using the appropriated taylor series when
 * $\cos(x) \approx -1$.
 *
 * Returns: $1 + \cos(x)$.
 */
/**
 * ncm_util_1mcosx:
 * @sinx: a real number $\sin(x)$
 * @cosx: a real number $\cos(x)$
 *
 * Calculates $1 - \cos(x)$ using the appropriated taylor series when
 * $\cos(x) \approx 1$.
 *
 * Returns: $1 - \cos(x)$.
 */
/**
 * ncm_util_1psinx:
 * @sinx: a real number $\sin(x)$
 * @cosx: a real number $\cos(x)$
 *
 * Calculates $1 + \sin(x)$ using the appropriated taylor series when
 * $\sin(x) \approx -1$.
 *
 * Returns: $1 + \sin(x)$.
 */
/**
 * ncm_util_1msinx:
 * @sinx: a real number $\sin(x)$
 * @cosx: a real number $\cos(x)$
 *
 * Calculates $1 - \sin(x)$ using the appropriated taylor series when
 * $\sin(x) \approx 1$.
 *
 * Returns: $1 - \sin(x)$.
 */
/**
 * ncm_util_cos2x:
 * @sinx: a real number $\sin(x)$
 * @cosx: a real number $\cos(x)$
 *
 * Calculates $\cos(2x)$ using the appropriated taylor series when
 * $\sin(x) \approx 1$.
 *
 * Returns: $1 - \sin(x)$.
 */

/**
 * ncm_cmpdbl:
 * @x: a double.
 * @y: a double.
 *
 * Compares @x and @y and returns the difference between them
 * relative to their mean.
 *
 * Returns: $\frac{2(x-y)}{x+y}$.
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
 * Computes the relative exponential $(\exp(x) - 1)/x$.
 *
 * Returns: $(\exp(x) - 1)/x$.
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
 * Computes the first derivative of the relative exponential $(\exp(x) - 1)/x$.
 *
 * Returns: first derivative of $(\exp(x) - 1)/x$.
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
 * Computes the second derivative of the relative exponential $(\exp(x) - 1)/x$.
 *
 * Returns: second derivative of $(\exp(x) - 1)/x$.
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
 * Computes the third derivative of the relative exponential $(\exp(x) - 1)/x$.
 *
 * Returns: third derivative of $(\exp(x) - 1)/x$.
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
 * Computes $\frac{\sinh(x)}{x}$. For small values of @x the taylor series is used.
 * For large values of @x the value is computed using the standard library function.
 *
 * Returns: $\frac{\sinh(x)}{x}$
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
 * Returns: $\frac{\sinh(x)-x}{x^3/3!}$
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
 * Returns: $\frac{\sinh(x)-x\cosh(x)}{x^3}$
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
 * ncm_util_mln_1mIexpzA_1pIexpmzA:
 * @rho: a double $\rho$
 * @theta: a double $\theta$
 * @A: a double $A$
 * @rho1: (out): a double $\rho_1$
 * @theta1: (out): a double $\theta_1$
 *
 * Computes $$z_1 = z - \ln\left(\frac{1-i e^{+z} A}{1+i e^{-z} A}\right),$$ where $z = \rho + i\theta$
 * and return the new $z_1 = \rho_1 + i\theta_1$ into @rho1 and $\theta1$.
 *
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

/**
 * ncm_cmp:
 * @x: a double
 * @y: a double
 * @reltol: relative precision
 * @abstol: the absolute precision
 *
 * Compare x and y and return -1 if x < y, 0 if x == y and 1 if x > y,
 * all comparisons are done with precision @reltol and @abstol.
 *
 * Returns: -1, 0, 1.
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
                              fabs (arg1 - arg2)
                             );

  g_assertion_message (domain, file, line, func, s);
  g_free (s);
}

G_DEFINE_BOXED_TYPE (NcmComplex, ncm_complex, ncm_complex_dup, ncm_complex_free)

/**
 * ncm_complex_new:
 *
 * Allocates a new complex number.
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
 * Allocates a new complex number and copy the contents of @c to it.
 *
 * Returns: (transfer full): a new #NcmComplex.
 */
NcmComplex *
ncm_complex_dup (NcmComplex *c)
{
  NcmComplex *cc = ncm_complex_new ();

  cc->z[0] = c->z[0];
  cc->z[1] = c->z[1];

  return cc;
}

/**
 * ncm_complex_free:
 * @c: a #NcmComplex
 *
 * Frees @c, it should not be used on a statically allocated NcmComplex.
 *
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
 * Frees *@c and sets *@c to NULL, it should not be used on a statically allocated NcmComplex.
 *
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
 * Sets @c to $a + I b$.
 *
 */
/**
 * ncm_complex_set_zero:
 * @c: a #NcmComplex
 *
 * Sets @c to $0 + I 0$.
 *
 */
/**
 * ncm_complex_Re:
 * @c: a #NcmComplex
 *
 * Returns the real part of @c.
 *
 * Returns: Re$(c)$.
 */
/**
 * ncm_complex_Im:
 * @c: a #NcmComplex
 *
 * Returns the imaginary part of @c.
 *
 * Returns: Im$(c)$.
 */

/**
 * ncm_complex_res_add_mul_real:
 * @c1: a #NcmComplex
 * @c2: a #NcmComplex
 * @v: a gdouble
 *
 * Computes @c1 = @c1 + @c2 * @v, assuming that
 * @c1 and @c2 are different.
 *
 */
/**
 * ncm_complex_res_add_mul:
 * @c1: a #NcmComplex
 * @c2: a #NcmComplex
 * @c3: #NcmComplex
 *
 * Computes @c1 = @c1 + @c2 * @c3, assuming that
 * @c1 and @c2 are different.
 *
 */

/**
 * ncm_complex_mul_real:
 * @c: a #NcmComplex
 * @v: a gdouble
 *
 * Computes @c1 = @c1 * @v.
 *
 */
/**
 * ncm_complex_res_mul:
 * @c1: a #NcmComplex
 * @c2: a #NcmComplex
 *
 * Computes @c1 = @c1 * @c2, assuming that
 * @c1 and @c2 are different.
 *
 */

/**
 * ncm_util_position_angle:
 * @ra1: Right ascension of object 1
 * @dec1: Declination of object 1
 * @ra2: Right ascension of object 2
 * @dec2: Declination of object 2
 *
 * Computes the on-sky position angle (East of North) between object1 (@ra1, @dec1) and object2 (@ra2, dec2).
 * The input coordinates ((@ra1, @dec1), (@ra2, @dec2)) must be given in decimal degrees.
 *
 * Returns: the position angle in radians
 */

/**
 * ncm_util_great_circle_distance:
 * @ra1: Right ascension of object 1
 * @dec1: Declination of object 1
 * @ra2: Right ascension of object 2
 * @dec2: Declination of object 2
 *
 * Compute the great circle distance (or separation, as defined in astropy) between poistion 1 (@ra1, @dec1) and position 2 (@ra2, @dec2).
 * See [Great-circle distance](https://en.wikipedia.org/wiki/Great-circle_distance), in particular the Vincenty equation (implemented here).
 * The input coordinates ((@ra1, @dec1), (@ra2, @dec2)) must be given in decimal degrees.
 *
 * Returns: the great circle distance in decimal degrees
 */

/**
 * ncm_util_projected_radius:
 * @theta: a gdouble in radians
 * @d: a gdouble in Mpc
 *
 * Converts the the angular separation `$\theta$' of a galaxy
 * at redshift `$z$' to the projected physical distance in Mpc.
 *
 * Returns: the physical distance in Mpc.
 */

/**
 * ncm_util_cvode_check_flag:
 * @flagvalue: pointer to flag value
 * @funcname: cvode function name
 * @opt: option
 *
 * Checks the CVode flag value and prints a message if an error occured.
 *
 * Returns: TRUE if no error occured, FALSE otherwise.
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
 * @cvode: a CVodeMem
 *
 * Prints the statistics of the CVodeMem object @cvode.
 *
 * Returns: TRUE.
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
 * @fits_filename: a fits filename
 *
 * Extracts the extension .fits or .fit from @fits_filename and returns
 * the prefix. If the extension is not found a copy of @fits_filename is
 * returned.
 *
 * Returns: (transfer full): prefix of @fits_filename.
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
 * @func: string representing function and its parameters
 * @x: (out) (array) (element-type double): the parameters or NULL if none found
 * @len: (out caller-allocates): number of parameters
 *
 * Extracts the function name and its numerical parameters.
 *
 * Returns: (transfer full): function name or NULL if it fails.
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
 * @n: a unsigned long integer
 *
 * Calculate the smallest factorization of @n such that
 * $n_f = 2^\mu \times 3^\nu \times 5^\alpha \times 7^\beta$
 * and $n_f \geq n$.
 *
 * This functions is useful to find a fft size such that fftw
 * can optimized it more easily.
 *
 * Returns: $n_f$
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
 * @milliseconds: sleep time in milliseconds
 *
 * Suspend the thread execution for @milliseconds.
 *
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

