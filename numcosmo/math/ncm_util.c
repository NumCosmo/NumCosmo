/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            ncm_util.c
 *
 *  Tue Jun  5 00:21:11 2007
 *  Copyright  2007  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */
#endif /* NUMCOSMO_GIR_SCAN */

/**
 * ncm_random_seed:
 *
 * FIXME
 *
 * Returns: FIXME
 */
gulong
ncm_random_seed ()
{
  guint seed;
  FILE *devrandom;

  if ( (devrandom = fopen("/dev/random","r")) == NULL )
  {
    g_warning ("Cannot open /dev/random, setting seed to 0\n");
    seed = 0;
  }
  else
  {
    if (fread (&seed, sizeof(seed), 1, devrandom) != 1)
      g_error ("ncm_random_seed: io error");
    fclose (devrandom);
  }

  return seed;
}

/**
 * ncm_smoothd:
 * @in: FIXME
 * @N: FIXME
 * @points: FIXME
 * @pass: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble *
ncm_smoothd (gdouble *in, size_t N, size_t points, size_t pass)
{
  gdouble *out = g_slice_alloc (sizeof(gdouble) * N);
  gdouble *temp = g_slice_alloc (sizeof(gdouble) * N);
  size_t i, j;
  g_assert (N > 0);
  g_assert (points > 0);
  g_assert (pass > 0);
  g_assert (N > points);

  for (j = 0; j < N; j++)
    temp[j] = out[j] = in[j];

  while (pass--)
  {
    for (i = points; i < N - points - 1; i++)
    {
      gdouble avg = 0.0;
      for (j = i - points; j <= points + i; j++)
        avg += temp[j];
      out[i] = avg / (2.0 * points + 1.0);
    }
    for (j = 0; j < N; j++)
      temp[j] = out[j];
    /*
     for (i = 0; i < N - points; i++)
     {
       gdouble avg = 0.0;
       for (j = i; j < points + i; j++)
       avg += temp[j];
       out[i] = avg / ((gdouble)points);
  }
  for (j = 0; j < N; j++)
  temp[j] = out[j];
  */
    /*
     for (i = 1; i < N - 1; i++)
     out[i] = (temp[i-1] + 2.0 * temp[i] + temp[i+1]) / 4.0;
     for (j = 0; j < N; j++)
     temp[j] = out[j];
     */
  }

  g_slice_free1 (sizeof(gdouble) * N, temp);
  return out;
}

/**
 * ncm_topology_comoving_a0_lss:
 * @n: FIXME
 * @alpha: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_topology_comoving_a0_lss (guint n, gdouble alpha)
{
  return atan (tan (M_PI / n) / cos (alpha));
}

/**
 * ncm_topology_sigma_comoving_a0_lss:
 * @n: FIXME
 * @alpha: FIXME
 * @sigma_alpha: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_topology_sigma_comoving_a0_lss (guint n, gdouble alpha, gdouble sigma_alpha)
{
  gdouble sigma = sin (alpha) * tan (M_PI/n) / (cos (alpha) * cos (alpha) + tan (M_PI / n) * tan (M_PI / n));
  sigma = sqrt (sigma * sigma * sigma_alpha * sigma_alpha);
  if (n == 2)
    sigma = sigma_alpha / alpha;
  return sigma;
}

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
  NcmCoarseDbl *cdbl = (NcmCoarseDbl *)p;

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
 * @x: FIXME
 * @q: FIXME
 *
 * FIXME
 *
 */
void
ncm_rational_coarce_double (gdouble x, mpq_t q)
{
  NcmCoarseDbl **cdbl_ptr = _ncm_coarse_dbl_get_bs ();
  NcmCoarseDbl *cdbl = *cdbl_ptr;
  gint expo2 = 0;
  gdouble xo = fabs(frexp (x, &expo2));
  gdouble xi = xo;

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

  //  mpq_set_d (q, x);
  //  mpq_canonicalize (q);
  //  mpfr_printf ("# AUTO: %.15g %.15f %d | %Qd\n", x, xo, expo2, q);
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
    //mpfr_printf ("# ---- %.5e %Zu %.15g | %Zd/%Zd %Zd/%Zd\n", fabs(1.0/(xi * mpz_get_d (ki))), a, xo, him1, kim1, him2, kim2);
    if (fabs(1.0/(xi * mpz_get_d (ki))) < GSL_DBL_EPSILON)
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
  if (GSL_SIGN(x) == -1)
    mpq_neg (q, q);

  //  mpfr_printf ("# MINE: %.15g %.15f %d | %Qd\n", x, xo, expo2, q);
  if (fabs(mpq_get_d (q) / x - 1) > 1e-15)
  {
    mpfr_fprintf (stderr, "# Q = %Qd\n", q);
    g_error ("Wrong rational approximation for x = %.16g N(q) = %.16g [%.5e] | 2^(%d).", x, mpq_get_d (q), fabs(mpq_get_d (q) / x - 1), expo2);
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

static gdouble
ncm_sphPlm_x_f (gdouble x, gpointer p)
{
  gint *params = p;
  gdouble val = fabs(gsl_sf_legendre_sphPlm (params[0], params[1], x));
  val = (val == 0.0) ? -300 : log (val);
  //  printf ("[%g] = %g\n", x, val);
  return (val+params[2]*M_LN10);
}

/**
 * ncm_sphPlm_x:
 * @l: FIXME
 * @m: FIXME
 * @order: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_sphPlm_x (gint l, gint m, gint order)
{
  gint status;
  gint iter = 0, max_iter = 1000000;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  gdouble x = 1.0, x0 = 0.1;
  gdouble prec = 1e-4, x1 = x;
  gint params[3];
  params[0] = l;
  params[1] = m;
  params[2] = order;

  F.function = &ncm_sphPlm_x_f;
  F.params = params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x0, x1);

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    if (status)
    {
      g_warning ("%s", gsl_strerror (status));
      gsl_root_fsolver_free (s);
      return GSL_NAN;
    }
    if (iter > 100)
      prec *= 10;

    x = gsl_root_fsolver_root (s);
    x0 = gsl_root_fsolver_x_lower (s);
    x1 = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x0, x1, 0, prec);

    if (!gsl_finite (ncm_sphPlm_x_f (x, params)))
    {
      g_debug ("Ops");
      x = GSL_NAN;
      break;
    }

    if ( (status == GSL_SUCCESS) && FALSE)
    {
      printf ("# Converged: %d, with prec = %g\n", iter, prec);
      //    }{
      printf ("#\t[%d]: (%g, %g):%g %g\n", iter, x0, x1, x, ncm_sphPlm_x_f (x, params));
      fflush (stdout);
    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);
  return x;
}

/**
 * ncm_mpfr_out_raw: (skip)
 * @stream: FIXME
 * @op: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gsize
ncm_mpfr_out_raw (FILE *stream, mpfr_t op)
{
  mpz_t m;
  mpz_init (m);
  gint64 expo = mpfr_get_z_exp (m, op);
  guint64 prec = mpfr_get_prec (op);
  gsize bout = 0;

  expo = GINT64_TO_BE (expo);
  prec = GUINT64_TO_BE (prec);

  bout += fwrite (&prec, 8, 1, stream);
  bout += fwrite (&expo, 8, 1, stream);
  bout += mpz_out_raw (stream, m);

  mpz_clear (m);
  return bout;
}

/**
 * ncm_mpfr_inp_raw: (skip)
 * @rop: FIXME
 * @stream: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gsize
ncm_mpfr_inp_raw (mpfr_t rop, FILE *stream)
{
  mpz_t m;
  mpz_init (m);
  gint64 expo;
  guint64 prec;
  gsize bin = 0;

  bin += fread (&prec, 8, 1, stream);
  bin += fread (&expo, 8, 1, stream);
  bin += mpz_inp_raw (m, stream);

  prec = GINT64_FROM_BE (prec);
  expo = GINT64_FROM_BE (expo);

  mpfr_set_prec (rop, prec);
  mpfr_set_z (rop, m, GMP_RNDN);
  if (expo != 0)
    mpfr_mul_2si (rop, rop, expo, GMP_RNDN);

  return bin;
}

/**
 * ncm_mpq_out_raw: (skip)
 * @f: FIXME
 * @q: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gsize
ncm_mpq_out_raw (FILE *f, mpq_t q)
{
  gsize t = 0;
  t += mpz_out_raw (f, mpq_numref (q));
  t += mpz_out_raw (f, mpq_denref (q));
  return t;
}

/**
 * ncm_mpq_inp_raw: (skip)
 * @q: FIXME
 * @f: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gsize
ncm_mpq_inp_raw (mpq_t q, FILE *f)
{
  gsize t = 0;
  t += mpz_inp_raw (mpq_numref (q), f);
  t += mpz_inp_raw (mpq_denref (q), f);
  return t;
}

/**
 * ncm_mpz_inits: (skip)
 * @z: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 */
void
ncm_mpz_inits (mpz_t z, ...)
{
  va_list ap;
  mpz_ptr z1;

  mpz_init (z);
  va_start(ap, z);
  while ((z1 = va_arg(ap, mpz_ptr)) != NULL)
    mpz_init (z1);
  va_end(ap);
}

/**
 * ncm_mpz_clears: (skip)
 * @z: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 */
void
ncm_mpz_clears (mpz_t z, ...)
{
  va_list ap;
  mpz_ptr z1;

  mpz_clear (z);
  va_start(ap, z);
  while ((z1 = va_arg(ap, mpz_ptr)) != NULL)
    mpz_clear (z1);
  va_end (ap);
}

/**
 * ncm_sum:
 * @d: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_sum (gdouble *d, gulong n)
{
  gdouble *part = g_slice_alloc (sizeof(gdouble) * n);
  gdouble res = 0.0;
  guint part_size = 0;
  guint i = 0;

  for (i = 0; i < n; i++)
  {
    guint j, k = 0;
    gdouble x = d[i];
    for (j = 0; j < part_size; j++)
    {
      const gdouble y = part[j];
      const gdouble hi = x + y;
      const gdouble lo = (fabs(x) < fabs(y)) ? (x - (hi - y)) : (y - (hi - x));
      if (lo != 0.0)
        part[k++] = lo;
      x = hi;
    }
    part[k] = x;
    part_size = k + 1;
  }

  for (i = 0; i < part_size; i++)
    res += part[i];

  g_slice_free1 (sizeof(gdouble) * n, part);

  return res;
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
 * Returns: FIXME
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
 * Returns: FIXME
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
 * Returns: FIXME
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
 * Returns: FIXME
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
 * Returns: FIXME
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
	{
		return ncm_util_sinh3 (x) / 6.0 - gsl_pow_2 (ncm_util_sinh1 (x)) / (sqrt (1.0 + shx * shx) + 1.0);
	}
	else
	{
		return (shx - x * sqrt (1.0 + shx * shx)) / gsl_pow_3 (x);
	}
}

/**
 * ncm_cmp:
 * @x: a double
 * @y: a double
 * @reltol: relative precision
 * @abstol: the absolute precision
 *
 * Compare x and y and return -1 if x < y, 0 if x == y and 1 if x > y,
 * all comparisons are done with precision @prec.
 *
 * Returns: -1, 0, 1.
 */
gint
ncm_cmp (gdouble x, gdouble y, const gdouble reltol, const gdouble abstol)
{
  if (G_UNLIKELY ((x == 0.0) && (y == 0.0)))
    return 0;
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

G_DEFINE_BOXED_TYPE (NcmComplex, ncm_complex, ncm_complex_dup, ncm_complex_free);

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
 * ncm_util_cvode_check_flag:
 * @flagvalue: FIXME
 * @funcname: FIXME
 * @opt: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
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
 * @cvode: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_util_cvode_print_stats (gpointer cvode)
{
  glong nsteps, nfunceval, nlinsetups, njaceval, ndiffjaceval, nnonliniter,
		nconvfail, nerrortests, nrooteval;
  gint flag, qcurorder, qlastorder;
	gdouble hinused, hlast, hcur, tcur;

	flag = CVodeGetIntegratorStats(cvode, &nsteps, &nfunceval,
	                               &nlinsetups, &nerrortests, &qlastorder, &qcurorder,
	                               &hinused, &hlast, &hcur, &tcur);
	ncm_util_cvode_check_flag (&flag, "CVodeGetIntegratorStats", 1);


  flag = CVodeGetNumNonlinSolvIters(cvode, &nnonliniter);
  ncm_util_cvode_check_flag (&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode, &nconvfail);
  ncm_util_cvode_check_flag (&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVodeGetNumJacEvals (cvode, &njaceval);
  ncm_util_cvode_check_flag (&flag, "CVodeGetNumJacEvals", 1);

  flag = CVodeGetNumRhsEvals (cvode, &ndiffjaceval);
  ncm_util_cvode_check_flag (&flag, "CVodeGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode, &nrooteval);
  ncm_util_cvode_check_flag (&flag, "CVodeGetNumGEvals", 1);

  flag = CVodeGetLastOrder(cvode, &qlastorder);
  ncm_util_cvode_check_flag (&flag, "CVodeGetLastOrder", 1);

  flag = CVodeGetCurrentOrder(cvode, &qcurorder);
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
  {
    base_name = g_match_info_fetch (match_info, 1);
  }
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
    gchar *fpa    = g_match_info_fetch (match_info, 2);
    func_name     = g_match_info_fetch (match_info, 1);

    if (fpa != NULL)
    {
      gchar **xs = g_regex_split_simple ("\\s*,\\s*", fpa, 0, 0);
      *len       = g_strv_length (xs);
      
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
 * ncm_util_print_bits:
 * @num: FIXME
 * 
 * FIXME
 * 
 */
void 
ncm_util_print_bits (guint64 num)
{
  gint bit;
  for(bit = 0; bit < (sizeof (guint64) * 8); bit++)
  {
    printf("%lu ", (glong)(num & 0x01));
    num = num >> 1;
  }
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
    return 0;
  else
  {
    const gulong r2 = n % 2;
    const gulong r3 = n % 3;
    const gulong r5 = n % 5;
    const gulong r7 = n % 7;
    gulong m = 1;

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
      return ncm_util_fact_size (n + 1);
  }
}


void
_ncm_util_set_destroyed (gpointer b)
{
  gboolean *destroyed = b;
  *destroyed = TRUE;
}

/********************************************************************
 *
 * File:          KolmogorovSmirnovDist.c
 * Environment:   ISO C99 or ANSI C89
 * Author:        Richard Simard
 * Organization:  DIRO, Université de Montréal
 * Date:          1 February 2012
 * Version        1.1

 * Copyright 1 march 2010 by Université de Montréal,
                             Richard Simard and Pierre L'Ecuyer
 =====================================================================

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 =====================================================================*/

#include <stdlib.h>

#define num_Pi     3.14159265358979323846 /* PI */
#define num_Ln2    0.69314718055994530941 /* log(2) */

/* For x close to 0 or 1, we use the exact formulae of Ruben-Gambino in all
   cases. For n <= NEXACT, we use exact algorithms: the Durbin matrix and
   the Pomeranz algorithms. For n > NEXACT, we use asymptotic methods
   except for x close to 0 where we still use the method of Durbin
   for n <= NKOLMO. For n > NKOLMO, we use asymptotic methods only and
   so the precision is less for x close to 0.
   We could increase the limit NKOLMO to 10^6 to get better precision
   for x close to 0, but at the price of a slower speed. */
#define NEXACT 500
#define NKOLMO 100000

/* The Durbin matrix algorithm for the Kolmogorov-Smirnov distribution */
static double DurbinMatrix (int n, double d);


/*========================================================================*/
#if 0

/* For ANSI C89 only, not for ISO C99 */
#define MAXI 50
#define EPSILON 1.0e-15

double log1p (double x)
{
   /* returns a value equivalent to log(1 + x) accurate also for small x. */
   if (fabs (x) > 0.1) {
      return log (1.0 + x);
   } else {
      double term = x;
      double sum = x;
      int s = 2;
      while ((fabs (term) > EPSILON * fabs (sum)) && (s < MAXI)) {
         term *= -x;
         sum += term / s;
         s++;
      }
      return sum;
   }
}

#undef MAXI
#undef EPSILON

#endif

/*========================================================================*/
#define MFACT 30

/* The natural logarithm of factorial n! for  0 <= n <= MFACT */
static double LnFactorial[MFACT + 1] = {
   0.,
   0.,
   0.6931471805599453,
   1.791759469228055,
   3.178053830347946,
   4.787491742782046,
   6.579251212010101,
   8.525161361065415,
   10.60460290274525,
   12.80182748008147,
   15.10441257307552,
   17.50230784587389,
   19.98721449566188,
   22.55216385312342,
   25.19122118273868,
   27.89927138384088,
   30.67186010608066,
   33.50507345013688,
   36.39544520803305,
   39.33988418719949,
   42.33561646075348,
   45.3801388984769,
   48.47118135183522,
   51.60667556776437,
   54.7847293981123,
   58.00360522298051,
   61.26170176100199,
   64.55753862700632,
   67.88974313718154,
   71.257038967168,
   74.65823634883016
};

/*------------------------------------------------------------------------*/

static double getLogFactorial (int n)
{
   /* Returns the natural logarithm of factorial n! */
   if (n <= MFACT) {
      return LnFactorial[n];

   } else {
      double x = (double) (n + 1);
      double y = 1.0 / (x * x);
      double z = ((-(5.95238095238E-4 * y) + 7.936500793651E-4) * y -
         2.7777777777778E-3) * y + 8.3333333333333E-2;
      z = ((x - 0.5) * log (x) - x) + 9.1893853320467E-1 + z / x;
      return z;
   }
}

/*------------------------------------------------------------------------*/

static double rapfac (int n)
{
   /* Computes n! / n^n */
   int i;
   double res = 1.0 / n;
   for (i = 2; i <= n; i++) {
      res *= (double) i / n;
   }
   return res;
}


/*========================================================================*/

static double **CreateMatrixD (int N, int M)
{
   int i;
   double **T2;

   T2 = (double **) malloc (N * sizeof (double *));
   T2[0] = (double *) malloc (N * M * sizeof (double));
   for (i = 1; i < N; i++)
      T2[i] = T2[0] + i * M;
   return T2;
}


static void DeleteMatrixD (double **T)
{
   free (T[0]);
   free (T);
}


/*========================================================================*/

static double KSPlusbarAsymp (int n, double x)
{
   /* Compute the probability of the KS+ distribution using an asymptotic
      formula */
   double t = (6.0 * n * x + 1);
   double z = t * t / (18.0 * n);
   double v = 1.0 - (2.0 * z * z - 4.0 * z - 1.0) / (18.0 * n);
   if (v <= 0.0)
      return 0.0;
   v = v * exp (-z);
   if (v >= 1.0)
      return 1.0;
   return v;
}


/*-------------------------------------------------------------------------*/

static double KSPlusbarUpper (int n, double x)
{
   /* Compute the probability of the KS+ distribution in the upper tail using
      Smirnov's stable formula */
   const double EPSILON = 1.0E-12;
   double q;
   double Sum = 0.0;
   double term;
   double t;
   double LogCom;
   double LOGJMAX;
   int j;
   int jdiv;
   int jmax = (int) (n * (1.0 - x));

   if (n > 200000)
      return KSPlusbarAsymp (n, x);

   /* Avoid log(0) for j = jmax and q ~ 1.0 */
   if ((1.0 - x - (double) jmax / n) <= 0.0)
      jmax--;

   if (n > 3000)
      jdiv = 2;
   else
      jdiv = 3;

   j = jmax / jdiv + 1;
   LogCom = getLogFactorial (n) - getLogFactorial (j) -
            getLogFactorial (n - j);
   LOGJMAX = LogCom;

   while (j <= jmax) {
      q = (double) j / n + x;
      term = LogCom + (j - 1) * log (q) + (n - j) * log1p (-q);
      t = exp (term);
      Sum += t;
      LogCom += log ((double) (n - j) / (j + 1));
      if (t <= Sum * EPSILON)
         break;
      j++;
   }

   j = jmax / jdiv;
   LogCom = LOGJMAX + log ((double) (j + 1) / (n - j));

   while (j > 0) {
      q = (double) j / n + x;
      term = LogCom + (j - 1) * log (q) + (n - j) * log1p (-q);
      t = exp (term);
      Sum += t;
      LogCom += log ((double) j / (n - j + 1));
      if (t <= Sum * EPSILON)
         break;
      j--;
   }

   Sum *= x;
   /* add the term j = 0 */
   Sum += exp (n * log1p (-x));
   return Sum;
}


/*========================================================================*/

static double Pelz (int n, double x)
{
   /* Approximating the Lower Tail-Areas of the Kolmogorov-Smirnov One-Sample
      Statistic,
      Wolfgang Pelz and I. J. Good,
      Journal of the Royal Statistical Society, Series B.
      Vol. 38, No. 2 (1976), pp. 152-156
   */

   const int JMAX = 20;
   const double EPS = 1.0e-10;
   const double C = 2.506628274631001;         /* sqrt(2*Pi) */
   const double C2 = 1.2533141373155001;       /* sqrt(Pi/2) */
   const double PI2 = num_Pi * num_Pi;
   const double PI4 = PI2 * PI2;
   const double RACN = sqrt ((double) n);
   const double z = RACN * x;
   const double z2 = z * z;
   const double z4 = z2 * z2;
   const double z6 = z4 * z2;
   const double w = PI2 / (2.0 * z * z);
   double ti, term, tom;
   double sum;
   int j;

   term = 1;
   j = 0;
   sum = 0;
   while (j <= JMAX && term > EPS * sum) {
      ti = j + 0.5;
      term = exp (-ti * ti * w);
      sum += term;
      j++;
   }
   sum *= C / z;

   term = 1;
   tom = 0;
   j = 0;
   while (j <= JMAX && fabs (term) > EPS * fabs (tom)) {
      ti = j + 0.5;
      term = (PI2 * ti * ti - z2) * exp (-ti * ti * w);
      tom += term;
      j++;
   }
   sum += tom * C2 / (RACN * 3.0 * z4);

   term = 1;
   tom = 0;
   j = 0;
   while (j <= JMAX && fabs (term) > EPS * fabs (tom)) {
      ti = j + 0.5;
      term = 6 * z6 + 2 * z4 + PI2 * (2 * z4 - 5 * z2) * ti * ti +
         PI4 * (1 - 2 * z2) * ti * ti * ti * ti;
      term *= exp (-ti * ti * w);
      tom += term;
      j++;
   }
   sum += tom * C2 / (n * 36.0 * z * z6);

   term = 1;
   tom = 0;
   j = 1;
   while (j <= JMAX && term > EPS * tom) {
      ti = j;
      term = PI2 * ti * ti * exp (-ti * ti * w);
      tom += term;
      j++;
   }
   sum -= tom * C2 / (n * 18.0 * z * z2);

   term = 1;
   tom = 0;
   j = 0;
   while (j <= JMAX && fabs (term) > EPS * fabs (tom)) {
      ti = j + 0.5;
      ti = ti * ti;
      term = -30 * z6 - 90 * z6 * z2 + PI2 * (135 * z4 - 96 * z6) * ti +
         PI4 * (212 * z4 - 60 * z2) * ti * ti + PI2 * PI4 * ti * ti * ti * (5 -
         30 * z2);
      term *= exp (-ti * w);
      tom += term;
      j++;
   }
   sum += tom * C2 / (RACN * n * 3240.0 * z4 * z6);

   term = 1;
   tom = 0;
   j = 1;
   while (j <= JMAX && fabs (term) > EPS * fabs (tom)) {
      ti = j * j;
      term = (3 * PI2 * ti * z2 - PI4 * ti * ti) * exp (-ti * w);
      tom += term;
      j++;
   }
   sum += tom * C2 / (RACN * n * 108.0 * z6);

   return sum;
}


/*=========================================================================*/

static void CalcFloorCeil (
   int n,                         /* sample size */
   double t,                      /* = nx */
   double *A,                     /* A_i */
   double *Atflo,                 /* floor (A_i - t) */
   double *Atcei                  /* ceiling (A_i + t) */
   )
{
   /* Precompute A_i, floors, and ceilings for limits of sums in the Pomeranz
      algorithm */
   int i;
   int ell = (int) t;             /* floor (t) */
   double z = t - ell;            /* t - floor (t) */
   double w = ceil (t) - t;

   if (z > 0.5) {
      for (i = 2; i <= 2 * n + 2; i += 2)
         Atflo[i] = i / 2 - 2 - ell;
      for (i = 1; i <= 2 * n + 2; i += 2)
         Atflo[i] = i / 2 - 1 - ell;

      for (i = 2; i <= 2 * n + 2; i += 2)
         Atcei[i] = i / 2 + ell;
      for (i = 1; i <= 2 * n + 2; i += 2)
         Atcei[i] = i / 2 + 1 + ell;

   } else if (z > 0.0) {
      for (i = 1; i <= 2 * n + 2; i++)
         Atflo[i] = i / 2 - 1 - ell;

      for (i = 2; i <= 2 * n + 2; i++)
         Atcei[i] = i / 2 + ell;
      Atcei[1] = 1 + ell;

   } else {                       /* z == 0 */
      for (i = 2; i <= 2 * n + 2; i += 2)
         Atflo[i] = i / 2 - 1 - ell;
      for (i = 1; i <= 2 * n + 2; i += 2)
         Atflo[i] = i / 2 - ell;

      for (i = 2; i <= 2 * n + 2; i += 2)
         Atcei[i] = i / 2 - 1 + ell;
      for (i = 1; i <= 2 * n + 2; i += 2)
         Atcei[i] = i / 2 + ell;
   }

   if (w < z)
      z = w;
   A[0] = A[1] = 0;
   A[2] = z;
   A[3] = 1 - A[2];
   for (i = 4; i <= 2 * n + 1; i++)
      A[i] = A[i - 2] + 1;
   A[2 * n + 2] = n;
}


/*========================================================================*/

static double Pomeranz (int n, double x)
{
   /* The Pomeranz algorithm to compute the KS distribution */
   const double EPS = 1.0e-15;
   const int ENO = 350;
   const double RENO = ldexp (1.0, ENO); /* for renormalization of V */
   int coreno;                    /* counter: how many renormalizations */
   const double t = n * x;
   double w, sum, minsum;
   int i, j, k, s;
   int r1, r2;                    /* Indices i and i-1 for V[i][] */
   int jlow, jup, klow, kup, kup0;
   double *A;
   double *Atflo;
   double *Atcei;
   double **V;
   double **H;                    /* = pow(w, j) / Factorial(j) */

   A = (double *) calloc ((size_t) (2 * n + 3), sizeof (double));
   Atflo = (double *) calloc ((size_t) (2 * n + 3), sizeof (double));
   Atcei = (double *) calloc ((size_t) (2 * n + 3), sizeof (double));
   V = (double **) CreateMatrixD (2, n + 2);
   H = (double **) CreateMatrixD (4, n + 2);

   CalcFloorCeil (n, t, A, Atflo, Atcei);

   for (j = 1; j <= n + 1; j++)
      V[0][j] = 0;
   for (j = 2; j <= n + 1; j++)
      V[1][j] = 0;
   V[1][1] = RENO;
   coreno = 1;

   /* Precompute H[][] = (A[j] - A[j-1]^k / k! for speed */
   H[0][0] = 1;
   w = 2.0 * A[2] / n;
   for (j = 1; j <= n + 1; j++)
      H[0][j] = w * H[0][j - 1] / j;

   H[1][0] = 1;
   w = (1.0 - 2.0 * A[2]) / n;
   for (j = 1; j <= n + 1; j++)
      H[1][j] = w * H[1][j - 1] / j;

   H[2][0] = 1;
   w = A[2] / n;
   for (j = 1; j <= n + 1; j++)
      H[2][j] = w * H[2][j - 1] / j;

   H[3][0] = 1;
   for (j = 1; j <= n + 1; j++)
      H[3][j] = 0;

   r1 = 0;
   r2 = 1;
   for (i = 2; i <= 2 * n + 2; i++) {
      jlow = 2 + (int) Atflo[i];
      if (jlow < 1)
         jlow = 1;
      jup = (int) Atcei[i];
      if (jup > n + 1)
         jup = n + 1;

      klow = 2 + (int) Atflo[i - 1];
      if (klow < 1)
         klow = 1;
      kup0 = (int) Atcei[i - 1];

      /* Find to which case it corresponds */
      w = (A[i] - A[i - 1]) / n;
      s = -1;
      for (j = 0; j < 4; j++) {
         if (fabs (w - H[j][1]) <= EPS) {
            s = j;
            break;
         }
      }
      /* assert (s >= 0, "Pomeranz: s < 0"); */

      minsum = RENO;
      r1 = (r1 + 1) & 1;          /* i - 1 */
      r2 = (r2 + 1) & 1;          /* i */

      for (j = jlow; j <= jup; j++) {
         kup = kup0;
         if (kup > j)
            kup = j;
         sum = 0;
         for (k = kup; k >= klow; k--)
            sum += V[r1][k] * H[s][j - k];
         V[r2][j] = sum;
         if (sum < minsum)
            minsum = sum;
      }

      if (minsum < 1.0e-280) {
         /* V is too small: renormalize to avoid underflow of probabilities */
         for (j = jlow; j <= jup; j++)
            V[r2][j] *= RENO;
         coreno++;                /* keep track of log of RENO */
      }
   }

   sum = V[r2][n + 1];
   free (A);
   free (Atflo);
   free (Atcei);
   DeleteMatrixD (H);
   DeleteMatrixD (V);
   w = getLogFactorial (n) - coreno * ENO * num_Ln2 + log (sum);
   if (w >= 0.)
      return 1.;
   return exp (w);
}


/*========================================================================*/

static double cdfSpecial (int n, double x)
{
   /* The KS distribution is known exactly for these cases */

   /* For nx^2 > 18, ncm_util_KSfbar(n, x) is smaller than 5e-16 */
   if ((n * x * x >= 18.0) || (x >= 1.0))
      return 1.0;

   if (x <= 0.5 / n)
      return 0.0;

   if (n == 1)
      return 2.0 * x - 1.0;

   if (x <= 1.0 / n) {
      double t = 2.0 * x * n - 1.0;
      double w;
      if (n <= NEXACT) {
         w = rapfac (n);
         return w * pow (t, (double) n);
      }
      w = getLogFactorial (n) + n * log (t / n);
      return exp (w);
   }

   if (x >= 1.0 - 1.0 / n) {
      return 1.0 - 2.0 * pow (1.0 - x, (double) n);
   }

   return -1.0;
}


/*========================================================================*/

/**
 * ncm_util_KScdf: (skip)
 * @n: FIXME
 * @x: FIXME
 *
 */
double ncm_util_KScdf (int n, double x)
{
   const double w = n * x * x;
   double u = cdfSpecial (n, x);
   if (u >= 0.0)
      return u;

   if (n <= NEXACT) {
      if (w < 0.754693)
         return DurbinMatrix (n, x);
      if (w < 4.0)
         return Pomeranz (n, x);
      return 1.0 - ncm_util_KSfbar (n, x);
   }

   if ((w * x * n <= 7.0) && (n <= NKOLMO))
      return DurbinMatrix (n, x);

   return Pelz (n, x);
}


/*=========================================================================*/

static double fbarSpecial (int n, double x)
{
   const double w = n * x * x;

   if ((w >= 370.0) || (x >= 1.0))
      return 0.0;
   if ((w <= 0.0274) || (x <= 0.5 / n))
      return 1.0;
   if (n == 1)
      return 2.0 - 2.0 * x;

   if (x <= 1.0 / n) {
      double z;
      double t = 2.0 * x * n - 1.0;
      if (n <= NEXACT) {
         z = rapfac (n);
         return 1.0 - z * pow (t, (double) n);
      }
      z = getLogFactorial (n) + n * log (t / n);
      return 1.0 - exp (z);
   }

   if (x >= 1.0 - 1.0 / n) {
      return 2.0 * pow (1.0 - x, (double) n);
   }
   return -1.0;
}


/*========================================================================*/

/**
 * ncm_util_KSfbar: (skip)
 * @n: FIXME
 * @x: FIXME
 *
 */
double ncm_util_KSfbar (int n, double x)
{
   const double w = n * x * x;
   double v = fbarSpecial (n, x);
   if (v >= 0.0)
      return v;

   if (n <= NEXACT) {
      if (w < 4.0)
         return 1.0 - ncm_util_KScdf (n, x);
      else
         return 2.0 * KSPlusbarUpper (n, x);
   }

   if (w >= 2.65)
      return 2.0 * KSPlusbarUpper (n, x);

   return 1.0 - ncm_util_KScdf (n, x);
}


/*=========================================================================

The following implements the Durbin matrix algorithm and was programmed by
G. Marsaglia, Wai Wan Tsang and Jingbo Wong.

I have made small modifications in their program. (Richard Simard)



=========================================================================*/

/*
 The C program to compute Kolmogorov's distribution

             K(n,d) = Prob(D_n < d),         where

      D_n = max(x_1-0/n,x_2-1/n...,x_n-(n-1)/n,1/n-x_1,2/n-x_2,...,n/n-x_n)

    with  x_1<x_2,...<x_n  a purported set of n independent uniform [0,1)
    random variables sorted into increasing order.
    See G. Marsaglia, Wai Wan Tsang and Jingbo Wong,
       J.Stat.Software, 8, 18, pp 1--4, (2003).
*/

#define NORM 1.0e140
#define INORM 1.0e-140
#define LOGNORM 140


/* Matrix product */
static void mMultiply (double *A, double *B, double *C, int m);

/* Matrix power */
static void mPower (double *A, int eA, double *V, int *eV, int m, int n);


static double DurbinMatrix (int n, double d)
{
   int k, m, i, j, g, eH, eQ;
   double h, s, *H, *Q;
   /* OMIT NEXT TWO LINES IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL */
#if 0
   s = d * d * n;
   if (s > 7.24 || (s > 3.76 && n > 99))
      return 1 - 2 * exp (-(2.000071 + .331 / sqrt (n) + 1.409 / n) * s);
#endif
   k = (int) (n * d) + 1;
   m = 2 * k - 1;
   h = k - n * d;
   H = (double *) malloc ((m * m) * sizeof (double));
   Q = (double *) malloc ((m * m) * sizeof (double));
   for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
         if (i - j + 1 < 0)
            H[i * m + j] = 0;
         else
            H[i * m + j] = 1;
   for (i = 0; i < m; i++) {
      H[i * m] -= pow (h, (double) (i + 1));
      H[(m - 1) * m + i] -= pow (h, (double) (m - i));
   }
   H[(m - 1) * m] += (2 * h - 1 > 0 ? pow (2 * h - 1, (double) m) : 0);
   for (i = 0; i < m; i++)
      for (j = 0; j < m; j++)
         if (i - j + 1 > 0)
            for (g = 1; g <= i - j + 1; g++)
               H[i * m + j] /= g;
   eH = 0;
   mPower (H, eH, Q, &eQ, m, n);
   s = Q[(k - 1) * m + k - 1];

   for (i = 1; i <= n; i++) {
      s = s * (double) i / n;
      if (s < INORM) {
         s *= NORM;
         eQ -= LOGNORM;
      }
   }
   s *= pow (10., (double) eQ);
   free (H);
   free (Q);
   return s;
}


static void mMultiply (double *A, double *B, double *C, int m)
{
   int i, j, k;
   double s;
   for (i = 0; i < m; i++)
      for (j = 0; j < m; j++) {
         s = 0.;
         for (k = 0; k < m; k++)
            s += A[i * m + k] * B[k * m + j];
         C[i * m + j] = s;
      }
}


static void renormalize (double *V, int m, int *p)
{
   int i;
   for (i = 0; i < m * m; i++)
      V[i] *= INORM;
   *p += LOGNORM;
}


static void mPower (double *A, int eA, double *V, int *eV, int m, int n)
{
   double *B;
   int eB, i;
   if (n == 1) {
      for (i = 0; i < m * m; i++)
         V[i] = A[i];
      *eV = eA;
      return;
   }
   mPower (A, eA, V, eV, m, n / 2);
   B = (double *) malloc ((m * m) * sizeof (double));
   mMultiply (V, V, B, m);
   eB = 2 * (*eV);
   if (B[(m / 2) * m + (m / 2)] > NORM)
      renormalize (B, m, &eB);

   if (n % 2 == 0) {
      for (i = 0; i < m * m; i++)
         V[i] = B[i];
      *eV = eB;
   } else {
      mMultiply (A, B, V, m);
      *eV = eA + eB;
   }

   if (V[(m / 2) * m + (m / 2)] > NORM)
      renormalize (V, m, eV);
   free (B);
}


/*=========================================================================*/
#if 0
#include <stdio.h>

int main (void)
{
   double x, y, z;
   const int K = 100;
   int n = 60;
   int j;
   printf ("n = %5d\n\n", n);
   printf ("      x                    cdf                     fbar\n");

   for (j = 0; j <= K; j++) {
      x = (double) j / K;
      y = ncm_util_KScdf (n, x);
      z = ncm_util_KSfbar (n, x);
      printf ("%8.3g     %22.15g      %22.15g\n", x, y, z);
   }
   return 0;
}
#endif

static double poly(const double *, int, double);

void
ncm_util_swilk (double *x, int n, double *w, double *pw, int *ifault)
{
    int nn2 = n / 2;
    double a[nn2 + 1]; /* 1-based */

/*	ALGORITHM AS R94 APPL. STATIST. (1995) vol.44, no.4, 547-551.

	Calculates the Shapiro-Wilk W test and its significance level
*/

    double small = 1e-19;

    /* polynomial coefficients */
    double g[2] = { -2.273,.459 };
    double c1[6] = { 0.,.221157,-.147981,-2.07119, 4.434685, -2.706056 };
    double c2[6] = { 0.,.042981,-.293762,-1.752461,5.682633, -3.582633 };
    double c3[4] = { .544,-.39978,.025054,-6.714e-4 };
    double c4[4] = { 1.3822,-.77857,.062767,-.0020322 };
    double c5[4] = { -1.5861,-.31082,-.083751,.0038915 };
    double c6[3] = { -.4803,-.082676,.0030302 };

    /* Local variables */
    int i, j, i1;

    double ssassx, summ2, ssumm2, gamma, range;
    double a1, a2, an, m, s, sa, xi, sx, xx, y, w1;
    double fac, asa, an25, ssa, sax, rsn, ssx, xsx;

    *pw = 1.;
    if (n < 3) { *ifault = 1; return;}

    an = (double) n;

    if (n == 3) {
	a[1] = 0.70710678;/* = sqrt(1/2) */
    } else {
	an25 = an + .25;
	summ2 = 0.;
	for (i = 1; i <= nn2; i++) {
    /* qnorm ((i - 0.375) / an25, 0.0, 1.0, 1, 0); */
	    a[i] = gsl_cdf_ugaussian_Pinv ((i - .375f) / an25);
	    double r__1 = a[i];
	    summ2 += r__1 * r__1;
	}
	summ2 *= 2.;
	ssumm2 = sqrt(summ2);
	rsn = 1. / sqrt(an);
	a1 = poly(c1, 6, rsn) - a[1] / ssumm2;

	/* Normalize a[] */
	if (n > 5) {
	    i1 = 3;
	    a2 = -a[2] / ssumm2 + poly(c2, 6, rsn);
	    fac = sqrt((summ2 - 2. * (a[1] * a[1]) - 2. * (a[2] * a[2]))
		       / (1. - 2. * (a1 * a1) - 2. * (a2 * a2)));
	    a[2] = a2;
	} else {
	    i1 = 2;
	    fac = sqrt((summ2 - 2. * (a[1] * a[1])) /
		       ( 1.  - 2. * (a1 * a1)));
	}
	a[1] = a1;
	for (i = i1; i <= nn2; i++) a[i] /= - fac;
    }

/*	Check for zero range */

    range = x[n - 1] - x[0];
    if (range < small) {*ifault = 6; return;}

/*	Check for correct sort order on range - scaled X */

    /* *ifault = 7; <-- a no-op, since it is changed below, in ANY CASE! */
    *ifault = 0;
    xx = x[0] / range;
    sx = xx;
    sa = -a[1];
    for (i = 1, j = n - 1; i < n; j--) {
	xi = x[i] / range;
	if (xx - xi > small) {
	    /* Fortran had:	 print *, "ANYTHING"
	     * but do NOT; it *does* happen with sorted x (on Intel GNU/linux 32bit):
	     *  shapiro.test(c(-1.7, -1,-1,-.73,-.61,-.5,-.24, .45,.62,.81,1))
	     */
	    *ifault = 7;
	}
	sx += xi;
	i++;
	if (i != j) sa += GSL_SIGN (i - j) * a[GSL_MIN (i, j)];
	xx = xi;
    }
    if (n > 5000) *ifault = 2;

/*	Calculate W statistic as squared correlation
	between data and coefficients */

    sa /= n;
    sx /= n;
    ssa = ssx = sax = 0.;
    for (i = 0, j = n - 1; i < n; i++, j--) {
	if (i != j) asa = GSL_SIGN (i - j) * a[1 + GSL_MIN (i, j)] - sa; else asa = -sa;
	xsx = x[i] / range - sx;
	ssa += asa * asa;
	ssx += xsx * xsx;
	sax += asa * xsx;
    }

/*	W1 equals (1-W) calculated to avoid excessive rounding error
	for W very near 1 (a potential problem in very large samples) */

    ssassx = sqrt(ssa * ssx);
    w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
    *w = 1. - w1;

/*	Calculate significance level for W */

    if (n == 3) {/* exact P value : */
	double pi6 = 1.90985931710274, /* = 6/pi */
	    stqr = 1.04719755119660; /* = asin(sqrt(3/4)) */
	*pw = pi6 * (asin(sqrt(*w)) - stqr);
	if(*pw < 0.) *pw = 0.;
	return;
    }
    y = log(w1);
    xx = log(an);
    if (n <= 11) {
	gamma = poly(g, 2, an);
	if (y >= gamma) {
	    *pw = 1e-99;/* an "obvious" value, was 'small' which was 1e-19f */
	    return;
	}
	y = -log(gamma - y);
	m = poly(c3, 4, an);
	s = exp(poly(c4, 4, an));
    } else {/* n >= 12 */
	m = poly(c5, 4, xx);
	s = exp(poly(c6, 3, xx));
    }
    /*DBG printf("c(w1=%g, w=%g, y=%g, m=%g, s=%g)\n",w1,*w,y,m,s); */

  /* pnorm (y, m, s, 0, 0) */
    *pw = gsl_cdf_gaussian_Q (y - m, s);

    return;
} /* swilk */

static double poly(const double *cc, int nord, double x)
{
/* Algorithm AS 181.2	Appl. Statist.	(1982) Vol. 31, No. 2

	Calculates the algebraic polynomial of order nord-1 with
	array of coefficients cc.  Zero order coefficient is cc(1) = cc[0]
*/
    double p, ret_val;

    ret_val = cc[0];
    if (nord > 1) {
	p = x * cc[nord-1];
	for (int j = nord - 2; j > 0; j--) p = (p + cc[j]) * x;
	ret_val += p;
    }
    return ret_val;
} /* poly */

