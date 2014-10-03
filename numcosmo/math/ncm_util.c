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
 * @title: Miscellaneous Utilities
 * @short_description: Miscellaneous Utilities
 *
 * Miscellaneous utility functions, macros and objects.
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_util.h"
#include "math/memory_pool.h"

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#ifdef NUMCOSMO_HAVE_CFITSIO
#include <fitsio.h>
#endif /* NUMCOSMO_HAVE_CFITSIO */

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
  sigma = sqrt(sigma * sigma * sigma_alpha * sigma_alpha);
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
_ncm_coarse_dbl_get_bs ()
{
  _NCM_STATIC_MUTEX_DECL (create_lock);
  static NcmMemoryPool *mp = NULL;

  _NCM_MUTEX_LOCK (&create_lock);
  if (mp == NULL)
    mp = ncm_memory_pool_new (_besselj_bs_alloc, NULL, _besselj_bs_free);
  _NCM_MUTEX_UNLOCK (&create_lock);

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
    g_error ("Wrong rational approximation for x = %.16g N(q) = %.16g [%.5e] | 2^(%d)\n", x, mpq_get_d (q), fabs(mpq_get_d (q) / x - 1), expo2);
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
 * ncm_sphPlm_test_theta:
 * @theta: FIXME
 * @lmax: FIXME
 * @lmin_data: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_sphPlm_test_theta (gdouble theta, gint lmax, gint *lmin_data)
{
  gdouble x = cos (theta);
  gdouble Plm_data[4096];
  gint m;
  gsl_vector_int_view lmin_view = gsl_vector_int_view_array (lmin_data, lmax+1);
  g_assert (lmax <= 4096);
  gsl_vector_int_set_all (&lmin_view.vector, 1e9);

  for (m = 0; m <= lmax; m++)
  {
    gint last = gsl_sf_legendre_array_size (lmax, m);
    gint l;
    gsl_sf_legendre_sphPlm_array (lmax, m, x, Plm_data);
    for (l = 0; l < last; l++)
      if (fabs(Plm_data[l]) > 1e-20)
    {
      lmin_data[m] = l+m;
      break;
    }
  }
  return 0.0;
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
 * ncm_numdiff_1: (skip)
 * @F: FIXME
 * @x: FIXME
 * @ho: FIXME
 * @err: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_numdiff_1 (gsl_function *F, const gdouble x, const gdouble ho, gdouble *err)
{
  const gint ntab = 20;
  const gdouble con = 2.0;
  const gdouble con2 = (con * con);
  const gdouble big = 1e300;
  const gdouble safe = 2.0;
  volatile gdouble temp = x + ho;
  const gdouble h = temp - x;
  gint i,j;
  gdouble errt, fac, hh, ans;
  gdouble a[ntab * ntab];
  //  guint count = 0;

  if (h == 0.0)
    g_error ("ncm_numdiff_1: Step h too small");

  hh = h;

  a[ntab * 0 + 0] = ans = (F->function (x + hh, F->params) - F->function (x - hh, F->params)) / (2.0*hh);
  //  count += 2;
  *err = big;

  for (i = 1; i < ntab; i++)
  {
    hh /= con;
    a[0 * ntab + i] = (F->function (x + hh, F->params) - F->function (x - hh, F->params)) / (2.0 * hh);
    //    count += 2;
    fac = con2;

    for (j = 1; j <= i; j++)
    {
      a[j * ntab + i] = (a[(j-1) * ntab + i] * fac - a[(j-1) * ntab + i - 1]) / (fac - 1.0);
      fac = con2 * fac;
      errt = GSL_MAX (fabs (a[j * ntab + i] - a[(j-1) * ntab + i]), fabs (a[j * ntab + i] - a[(j-1) * ntab + i - 1]));
      if (errt <= *err)
      {
        *err = errt;
        ans = a[j * ntab + i];
      }
    }
    if (fabs (a[i*ntab + i] - a[(i-1) * ntab + i - 1]) >= safe * (*err))
      break;
  }
  //  printf ("# Count %u\n", count);
  return ans;
}

#define NCM_NUMDIFF_RETRY_PP (0.1)
#define NCM_NUMDIFF_RETRY_N (8)

/**
 * ncm_numdiff_2_err: (skip)
 * @F: FIXME
 * @ofx: FIXME
 * @x: FIXME
 * @ho: FIXME
 * @err: FIXME
 * @ferr: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_numdiff_2_err (gsl_function *F, gdouble *ofx, const gdouble x, const gdouble ho, gdouble err, gdouble *ferr)
{
  gdouble herr;
  gdouble hans, lans;
  gdouble try_h = ho;
  gdouble try_again = FALSE;
  gdouble c = NCM_NUMDIFF_RETRY_PP;
  gint deteriorate = 0;

  lans = ncm_numdiff_2 (F, ofx, x, ho, ferr);
  //printf ("#D2 % 20.15g +/- % 20.15g [% 20.15g] | (% 20.15g, % 20.15g)\n", lans, lerr, ho, x, *ofx);
  if (fabs (ferr[0] / lans) < err)
    return lans;
  try_h = ho;

  do {
    try_h *= c;
    hans = ncm_numdiff_2 (F, ofx, x, try_h, &herr);
    //printf ("#   % 20.15g +/- % 20.15g [% 20.15g]<%d> lerr % 20.15g\n", hans, herr, try_h, deteriorate, lerr);
    if (fabs(herr) > fabs(ferr[0]))
    {
      if (deteriorate < NCM_NUMDIFF_RETRY_N)
      {
        deteriorate++;
        try_again = TRUE;
      }
      else if (c == NCM_NUMDIFF_RETRY_PP)
      {
        try_h = ho;
        c = 1.0 + NCM_NUMDIFF_RETRY_PP;
        try_again = TRUE;
        deteriorate = 0;
      }
      else
      {
        try_again = FALSE;
        lans = hans;
      }
    }
    else
    {
      try_again = (fabs(herr/hans) > fabs(err));
      lans = hans;
      ferr[0] = herr;
    }
  } while (try_again);

  return lans;
}

/**
 * ncm_numdiff_2: (skip)
 * @F: FIXME
 * @ofx: FIXME
 * @x: FIXME
 * @ho: FIXME
 * @err: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gdouble
ncm_numdiff_2 (gsl_function *F, gdouble *ofx, const gdouble x, const gdouble ho, gdouble *err)
{
  const gint ntab = 20;
  const gdouble con = 2.0;
  const gdouble con2 = (con * con);
  const gdouble big = 1e300;
  const gdouble safe = 2.0;
  volatile gdouble temp = x + ho;
  const gdouble h = temp - x;
  gint i,j;
  gdouble errt, fac, hh, ans;
  gdouble a[ntab * ntab];
  gdouble fx;
  //  guint count = 0;
  if (ofx == NULL)
    fx = F->function (x, F->params);
  else
    fx = *ofx;

  if (h == 0.0)
    g_error ("ncm_numdiff_2: Step h too small");

  hh = h;

  a[ntab * 0 + 0] = ans = (F->function (x + hh, F->params) - 2.0 * fx + F->function (x - hh, F->params)) / (hh * hh);
  //  count += 2;
  *err = big;

  for (i = 1; i < ntab; i++)
  {
    hh /= con;
    a[0 * ntab + i] = (F->function (x + hh, F->params) - 2.0 * fx + F->function (x - hh, F->params)) / (hh * hh);
    //    count += 2;
    fac = con2;

    for (j = 1; j <= i; j++)
    {
      a[j * ntab + i] = (a[(j-1) * ntab + i] * fac - a[(j-1) * ntab + i - 1]) / (fac - 1.0);
      fac = con2 * fac;
      errt = GSL_MAX (fabs (a[j * ntab + i] - a[(j-1) * ntab + i]), fabs (a[j * ntab + i] - a[(j-1) * ntab + i - 1]));
      if (errt <= *err)
      {
        *err = errt;
        ans = a[j * ntab + i];
      }
    }
    if (fabs (a[i*ntab + i] - a[(i-1) * ntab + i - 1]) >= safe * (*err))
      break;
  }

  return ans;
}

/**
 * ncm_sqrt1px_m1:
 * @x: a real number $&gt;-1$
 * 
 * Calculates $\sqrt{1+x}-1$ using the appropriated taylor series when 
 * $x \approx 1$.
 * 
 * Returns: $\sqrt{1+x}-1$.
 */
gdouble
ncm_sqrt1px_m1 (gdouble x)
{
  gdouble binfact = 1.0;
  gdouble res = 0;
  gdouble xn = 1;
  gint n = 0;

  if (!gsl_finite (x))
    return x;

  if (fabs(x) > 1e-1)
    return sqrt (1.0 + x) - 1.0;

  while (TRUE)
  {
    binfact *= 3.0 / (2.0 * (1.0 + n++)) - 1.0;
    xn *= x;
    res += binfact * xn;

    if (fabs(binfact * xn / res) < GSL_DBL_EPSILON)
      break;
  }

  return res;
}

/**
 * ncm_cmp:
 * @x: a double.
 * @y: a double.
 * @reltol: relative precision.
 * 
 * Compare x and y and return -1 if x < y, 0 if x == y and 1 if x > y,
 * all comparisons are done with precision @prec.
 * 
 * Returns: -1, 0, 1.
 */
gint 
ncm_cmp (gdouble x, gdouble y, gdouble reltol)
{
  if (G_UNLIKELY (x == 0 && y == 0))
    return 0;
  else
  {
    const gdouble delta = (x - y);
    const gdouble mean  = G_UNLIKELY (x == 0 || y == 0) ? 1.0 : 0.5 * fabs (x + y);
    if (fabs (delta / mean) < reltol)
      return 0;
    else
      return delta < 0 ? -1 : 1;
  }
}

void
_ncm_assertion_message_cmpdouble (const gchar *domain, const gchar *file, gint line, const gchar *func, const gchar *expr, gdouble arg1, const gchar *cmp, gdouble arg2)
{
  gchar *s = g_strdup_printf ("assertion failed (%s): (%.17g %s %.17g)", expr, arg1, cmp, arg2);
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
 * @c: a #NcmComplex.
 * 
 * Allocates a new complex number and copy the contents of @c to it.
 * 
 * Returns: (transfer full): a new #NcmComplex.
 */
NcmComplex *
ncm_complex_dup (NcmComplex *c)
{
  NcmComplex *cc = ncm_complex_new ();
printf ("Nhaca %p\n", c);
  cc->z = c->z;
  return cc;
}

/**
 * ncm_complex_free:
 * @c: a #NcmComplex.
 * 
 * Frees @c, it should not be used on a statically allocated NcmComplex.
 * 
 */
void 
ncm_complex_free (NcmComplex *c)
{
printf ("Nhoco %p\n", c);
  g_free (c);
}

/**
 * ncm_complex_clear:
 * @c: a #NcmComplex.
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
 * ncm_complex_Re:
 * @c: a #NcmComplex.
 * 
 * Returns the real part of @c.
 * 
 * Returns: Re$(c)$.
 */
gdouble 
ncm_complex_Re (NcmComplex *c)
{
  return creal (c->z);
}

/**
 * ncm_complex_Im:
 * @c: a #NcmComplex.
 * 
 * Returns the imaginary part of @c.
 * 
 * Returns: Im$(c)$.
 */
gdouble 
ncm_complex_Im (NcmComplex *c)
{
  printf ("Huga %p\n", c);
  return cimag (c->z);
}

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
ncm_util_cvode_check_flag (gpointer flagvalue, gchar *funcname, gint opt)
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

  flag = CVDlsGetNumJacEvals (cvode, &njaceval);
  ncm_util_cvode_check_flag (&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals (cvode, &ndiffjaceval);
  ncm_util_cvode_check_flag (&flag, "CVDlsGetNumRhsEvals", 1);  
    
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
