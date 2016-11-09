/***************************************************************************
 *            ncm_util.h
 *
 *  Mon Jul 16 18:02:22 2007
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

#ifndef _NCM_UTIL_H_
#define _NCM_UTIL_H_

#include <glib.h>
#include <glib-object.h>
#include <numcosmo/build_cfg.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_min.h>
#ifndef NUMCOSMO_GIR_SCAN
#include <complex.h>
#include <gmp.h>
#include <mpfr.h>
#endif

G_BEGIN_DECLS

gdouble *ncm_smoothd (gdouble *in, size_t N, size_t points, size_t pass);
gdouble ncm_topology_comoving_a0_lss (guint n, gdouble alpha);
gdouble ncm_topology_sigma_comoving_a0_lss (guint n, gdouble alpha, gdouble sigma_alpha);
gdouble ncm_sphPlm_x (gint l, gint m, gint order);
gdouble ncm_sphPlm_test_theta (gdouble theta, gint lmax, gint *lmin_data);
gdouble ncm_sum (gdouble *d, gulong n);
gdouble ncm_numdiff_1 (gsl_function *F, const gdouble x, const gdouble ho, gdouble *err);
gdouble ncm_numdiff_2 (gsl_function *F, gdouble *ofx, const gdouble x, const gdouble ho, gdouble *err);
gdouble ncm_numdiff_2_err (gsl_function *F, gdouble *ofx, const gdouble x, const gdouble ho, gdouble err, gdouble *ferr);
gdouble ncm_sqrt1px_m1 (gdouble x);
gdouble ncm_cmpdbl (const gdouble x, const gdouble y);
gdouble ncm_exprel (const gdouble x);
gdouble ncm_d1exprel (const gdouble x);
gdouble ncm_d2exprel (const gdouble x);
gdouble ncm_d3exprel (const gdouble x);

gsize ncm_mpfr_out_raw (FILE *stream, mpfr_t op);
gsize ncm_mpfr_inp_raw (mpfr_t rop, FILE *stream);
gsize ncm_mpq_out_raw (FILE *f, mpq_t q);
gsize ncm_mpq_inp_raw (mpq_t q, FILE *f);

gulong ncm_random_seed (void);

gint ncm_cmp (gdouble x, gdouble y, gdouble reltol);

void ncm_rational_coarce_double (gdouble x, mpq_t q);
void ncm_mpz_inits (mpz_t z, ...) G_GNUC_NULL_TERMINATED;
void ncm_mpz_clears (mpz_t z, ...) G_GNUC_NULL_TERMINATED;
void _ncm_assertion_message_cmpdouble (const gchar *domain, const gchar *file, gint line, const gchar *func, const gchar *expr, gdouble arg1, const gchar *cmp, gdouble arg2);

gboolean ncm_util_cvode_check_flag (gpointer flagvalue, const gchar *funcname, gint opt);
gboolean ncm_util_cvode_print_stats (gpointer cvode);

gchar *ncm_util_basename_fits (const gchar *fits_filename); 
gchar *ncm_util_function_params (const gchar *func, gdouble **x, guint *len);

typedef struct _NcmComplex NcmComplex;

struct _NcmComplex
{
  /*< private >*/
#ifndef NUMCOSMO_GIR_SCAN
  complex double z;
#else
  gdouble Rez, Imz;
#endif /* NUMCOSMO_GIR_SCAN */
};

GType ncm_complex_get_type (void) G_GNUC_CONST;

NcmComplex *ncm_complex_new (void);
NcmComplex *ncm_complex_ref (NcmComplex *c);
NcmComplex *ncm_complex_dup (NcmComplex *c);
void ncm_complex_free (NcmComplex *c);
void ncm_complex_clear (NcmComplex **c);

gdouble ncm_complex_Re (NcmComplex *c);
gdouble ncm_complex_Im (NcmComplex *c);

G_INLINE_FUNC gdouble ncm_util_smooth_trans (gdouble f0, gdouble f1, gdouble z0, gdouble dz, gdouble z);
G_INLINE_FUNC void ncm_util_smooth_trans_get_theta (gdouble z0, gdouble dz, gdouble z, gdouble *theta0, gdouble *theta1);

/* Macros */

#ifndef HAVE_EXP10
#define exp10(x) (exp ((x) * M_LN10))
#endif /* HAVE_EXP10 */

#ifndef NUMCOSMO_GIR_SCAN
#ifndef HAVE_SINCOS
G_INLINE_FUNC void sincos (gdouble x, gdouble *s, gdouble *c);
#endif
#endif

#define ncm_acb_get_complex(z) (arf_get_d (arb_midref (acb_realref (z)), ARF_RND_NEAR) + I * arf_get_d (arb_midref (acb_imagref (z)), ARF_RND_NEAR))

#define ncm_util_exp10(x) (exp ((x) * M_LN10))

#define NCM_GARRAY_MEMCPY(dest,src) \
G_STMT_START { \
g_assert_cmpuint ((src)->len, ==, (dest)->len); \
g_assert_cmpuint (g_array_get_element_size (src), ==, g_array_get_element_size (dest)); \
memcpy ((dest)->data, (src)->data, (src)->len * g_array_get_element_size (src)); \
} G_STMT_END

#define NCM_GARRAY_DUP(dest,src) \
G_STMT_START { \
dest = g_array_sized_new (FALSE, FALSE, g_array_get_element_size (src), (src)->len); \
g_array_set_size ((dest), (src)->len); \
memcpy ((dest)->data, (src)->data, (src)->len * g_array_get_element_size (src)); \
} G_STMT_END

#define ncm_assert_cmpdouble(n1,cmp,n2) \
do { \
  if (ncm_cmp ((n1), (n2), GSL_DBL_EPSILON) cmp 0) ; else \
    _ncm_assertion_message_cmpdouble (G_LOG_DOMAIN, __FILE__, __LINE__, G_STRFUNC, \
                                      #n1 " " #cmp " " #n2, (n1), #cmp, (n2)); \
} while (0)

#define ncm_assert_cmpdouble_e(n1,cmp,n2,epsilon) \
do { \
  if (ncm_cmp ((n1), (n2), (epsilon)) cmp 0) ; else \
    _ncm_assertion_message_cmpdouble (G_LOG_DOMAIN, __FILE__, __LINE__, G_STRFUNC, \
                                      #n1 " " #cmp " " #n2, (n1), #cmp, (n2)); \
} while (0)

#define NCM_RETURN_IF_INF(a) if (gsl_isinf(a)) return a

#define NCM_FLOOR_TRUNC(a,b) (floor ((b) * (a)) / (b))
#define NCM_CEIL_TRUNC(a,b) (ceil ((b) * (a)) / (b))
#define NCM_ROUND_TRUNC(a,b) (round ((b) * (a)) / (b))

#define NCM_TEST_GSL_RESULT(func,ret) if (ret != GSL_SUCCESS) g_error ("%s: %s", func, gsl_strerror (ret))

#define NCM_COMPLEX_INC_MUL_REAL_TEST(a,b,c) \
((fabs(GSL_REAL((b))*(c)/GSL_REAL((a))) < 1e-16) && (fabs(GSL_IMAG((b))*(c)/GSL_IMAG((a))) < 1e-16))

#define NCM_COMPLEX_INC_MUL_REAL(a,b,c) \
G_STMT_START { \
  GSL_REAL((a)) += GSL_REAL((b))*(c); \
  GSL_IMAG((a)) += GSL_IMAG((b))*(c); \
} G_STMT_END

#define NCM_COMPLEX_INC_MUL(a,b,c) \
G_STMT_START { \
  GSL_REAL((a)) += GSL_REAL((b)) * GSL_REAL((c)) - GSL_IMAG((b)) * GSL_IMAG((c)); \
  GSL_IMAG((a)) += GSL_REAL((b)) * GSL_IMAG((c)) + GSL_IMAG((b)) * GSL_REAL((c)); \
} G_STMT_END

#define NCM_COMPLEX_INC_MUL_MUL_REAL(a,b,c,d) \
G_STMT_START { \
  GSL_REAL((a)) += (GSL_REAL((b)) * GSL_REAL((c)) - GSL_IMAG((b)) * GSL_IMAG((c))) * (d); \
  GSL_IMAG((a)) += (GSL_REAL((b)) * GSL_IMAG((c)) + GSL_IMAG((b)) * GSL_REAL((c))) * (d); \
} G_STMT_END

#define NCM_COMPLEX_MUL_REAL(a,b,c) \
G_STMT_START { \
  GSL_REAL((a)) = GSL_REAL((b)) * (c); \
  GSL_IMAG((a)) = GSL_IMAG((b)) * (c); \
} G_STMT_END

#define NCM_COMPLEX_MUL(a,b) \
G_STMT_START { \
  gdouble temp = GSL_REAL((a)) * GSL_REAL((b)) - GSL_IMAG((a)) * GSL_IMAG((b)); \
  GSL_IMAG(a) = GSL_REAL((a)) * GSL_IMAG((b)) + GSL_IMAG((a)) * GSL_REAL((b)); \
  GSL_REAL(a) = temp; \
} G_STMT_END

#define NCM_COMPLEX_ADD(a,b) \
G_STMT_START { \
  GSL_REAL(a) += GSL_REAL((b)); \
  GSL_IMAG(a) += GSL_IMAG((b)); \
} G_STMT_END

#define NCM_COMPLEX_MUL_CONJUGATE(a,b) \
G_STMT_START { \
  gdouble temp = GSL_REAL((a)) * GSL_REAL((b)) + GSL_IMAG((a)) * GSL_IMAG((b)); \
  GSL_IMAG(a) = - GSL_REAL((a)) * GSL_IMAG((b)) + GSL_IMAG((a)) * GSL_REAL((b)); \
  GSL_REAL(a) = temp; \
} G_STMT_END

#define NCM_WRITE_INT32(_ff,_ii) G_STMT_START { gint32 _temp_i = GINT32_TO_BE ((_ii)); if (fwrite (&_temp_i, sizeof(gint32), (1), _ff) != 1) g_error ("NCM_WRITE_INT32: io error"); } G_STMT_END
#define NCM_WRITE_UINT32(_ff,_ii) G_STMT_START { guint32 _temp_i = GUINT32_TO_BE ((_ii)); if (fwrite (&_temp_i, sizeof(guint32), (1), _ff) != 1) g_error ("NCM_WRITE_UINT32: io error"); } G_STMT_END

#define NCM_WRITE_INT64(_ff,_ii) G_STMT_START { gint64 _temp_i = GINT64_TO_BE ((_ii)); if (fwrite (&_temp_i, sizeof(gint64), (1), _ff) != 1) g_error ("NCM_WRITE_INT64: io error"); } G_STMT_END
#define NCM_WRITE_UINT64(_ff,_ii) G_STMT_START { guint64 _temp_i = GUINT64_TO_BE ((_ii)); if (fwrite (&_temp_i, sizeof(guint64), (1), _ff) != 1) g_error ("NCM_WRITE_INT64: io error"); } G_STMT_END

#define NCM_WRITE_DOUBLE(_ff,_ii) G_STMT_START { NcmDoubleInt64 _iii; _iii.x = _ii; _iii.i = GINT64_TO_BE ((_iii.i)); if (fwrite (&_iii.i, sizeof(gint64), (1), _ff) != 1) g_error ("NCM_WRITE_DOUBLE: io error"); } G_STMT_END

#define NCM_READ_INT32(_ff,_ii) G_STMT_START { gint32 _temp_i; if (fread (&_temp_i, sizeof(gint32), (1), _ff) != 1) g_error ("NCM_READ_INT32: io error"); _ii = GINT32_FROM_BE (_temp_i); } G_STMT_END
#define NCM_READ_UINT32(_ff,_ii) G_STMT_START { guint32 _temp_i; if (fread (&_temp_i, sizeof(guint32), (1), _ff) != 1) g_error ("NCM_READ_UINT32: io error"); _ii = GUINT32_FROM_BE (_temp_i); } G_STMT_END

#define NCM_READ_INT64(_ff,_ii) G_STMT_START { gint64 _temp_i; if (fread (&_temp_i, sizeof(gint64), (1), _ff) != 1) g_error ("NCM_READ_INT64: io error"); _ii = GINT64_FROM_BE (_temp_i); } G_STMT_END
#define NCM_READ_UINT64(_ff,_ii) G_STMT_START { guint64 _temp_i; if (fread (&_temp_i, sizeof(guint64), (1), _ff) != 1) g_error ("NCM_READ_UINT64: io error"); _ii = GUINT64_FROM_BE (_temp_i); } G_STMT_END

#define NCM_READ_DOUBLE(_ff,_ii) G_STMT_START { NcmDoubleInt64 _iii; if (fread (&_iii.i, sizeof(gint64), (1), _ff) != 1) g_error ("NCM_READ_DOUBLE: io error"); _iii.i = GINT64_FROM_BE (_iii.i); _ii = _iii.x; } G_STMT_END

#define ncm_g_string_clear(s) G_STMT_START if (*(s) != NULL) { g_string_free (*(s), TRUE); *(s) = NULL; } G_STMT_END

#define NCM_UNUSED(x) (void)(x)

void _ncm_util_set_destroyed (gpointer b);

#define NCM_TEST_FREE(cmd,obj) \
G_STMT_START { \
  gboolean destroyed = FALSE; \
  g_object_set_data_full (G_OBJECT (obj), "test-destroy", &destroyed, _ncm_util_set_destroyed); \
  cmd (obj); \
  g_assert (destroyed); \
} G_STMT_END

/* Minumum version here is 2.38 but it segfault during tests so we start at 2.40. */
#if ((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 40))
#define NCM_TEST_FAIL(cmd) \
G_STMT_START { \
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR)) \
  { \
    cmd; \
    exit (0); \
  } \
  g_test_trap_assert_failed (); \
} G_STMT_END

#define NCM_TEST_PASS(cmd) \
G_STMT_START { \
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR)) \
  { \
    cmd; \
    exit (0); \
  } \
  g_test_trap_assert_passed (); \
} G_STMT_END

#else

#define NCM_TEST_FAIL(cmd) \
G_STMT_START { \
  if (g_test_subprocess ()) \
  { \
    cmd; \
    exit (0); \
  } \
  else \
  { \
    g_test_trap_subprocess (NULL, 0, 0); \
    g_test_trap_assert_failed (); \
  } \
} G_STMT_END

#define NCM_TEST_PASS(cmd) \
G_STMT_START { \
  if (g_test_subprocess ()) \
  { \
    cmd; \
    exit (0); \
  } \
  else \
  { \
    g_test_trap_subprocess (NULL, 0, 0); \
    g_test_trap_assert_passed (); \
  } \
} G_STMT_END
#endif /* !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 38)) */

#define NCM_CVODE_CHECK(chk,name,val,ret) \
G_STMT_START { \
  if (!ncm_util_cvode_check_flag (chk, name, val)) \
    return ret; \
} G_STMT_END

G_END_DECLS
#endif /* _NCM_UTIL_H_ */

#ifndef _NCM_UTIL_INLINE_H_
#define _NCM_UTIL_INLINE_H_
#ifdef NUMCOSMO_HAVE_INLINE

G_BEGIN_DECLS

G_INLINE_FUNC gdouble 
ncm_util_smooth_trans (gdouble f0, gdouble f1, gdouble z0, gdouble dz, gdouble z)
{
  const gdouble C0      = 18.0;
  const gdouble Delta   = 0.25 * dz / C0;
  const gdouble a       = - z0 - 0.5 * dz;
  const gdouble gz      = (z + a) / Delta;
  const gdouble exp_gz  = exp (gz);
  const gdouble exp_mgz = 1.0 / exp_gz;
  const gdouble f       = f0 / (1.0 + exp_gz) + f1 / (1.0 + exp_mgz);

  return f;
}

G_INLINE_FUNC void 
ncm_util_smooth_trans_get_theta (gdouble z0, gdouble dz, gdouble z, gdouble *theta0, gdouble *theta1)
{
  const gdouble C0      = 18.0;
  const gdouble Delta   = 0.25 * dz / C0;
  const gdouble a       = - z0 - 0.5 * dz;
  const gdouble gz      = (z + a) / Delta;
  const gdouble exp_gz  = exp (gz);
  const gdouble exp_mgz = 1.0 / exp_gz;
  theta0[0] = 1.0 / (1.0 + exp_gz);
  theta1[0] = 1.0 / (1.0 + exp_mgz);
}

#ifndef NUMCOSMO_GIR_SCAN
#ifndef HAVE_SINCOS
G_INLINE_FUNC void 
sincos (gdouble x, gdouble *s, gdouble *c)
{
  s[0] = sin (x);
  c[0] = cos (x);
}
#endif
#endif

G_END_DECLS

#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_UTIL_INLINE_H_ */

#ifndef KOLMOGOROVSMIRNOVDIST_H
#define KOLMOGOROVSMIRNOVDIST_H

#ifdef __cplusplus
extern "C" {
#endif

/********************************************************************
 *
 * File:          KolmogorovSmirnovDist.h
 * Environment:   ISO C99 or ANSI C89
 * Author:        Richard Simard
 * Organization:  DIRO, Université de Montréal
 * Date:          1 February 2012
 * Version        1.1
 *
 * Copyright March 2010 by Université de Montréal,
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
/*
 *
 * The Kolmogorov-Smirnov test statistic D_n is defined by
 *
 *        D_n = sup_x |F(x) - S_n(x)|
 *
 * where n is the sample size, F(x) is a completely specified theoretical
 * distribution, and S_n(x) is an empirical distribution function.
 *
 *
 * The function
 *
 *        double ncm_util_KScdf (int n, double x);
 *
 * computes the cumulative probability P[D_n <= x] of the 2-sided 1-sample
 * Kolmogorov-Smirnov distribution with sample size n at x.
 * It returns at least 13 decimal digits of precision for n <= 500,
 * at least 7 decimal digits of precision for 500 < n <= 100000,
 * and a few correct decimal digits for n > 100000.
 *
 */
double ncm_util_KScdf (int n, double x);

/*
 * The function
 *
 *        double ncm_util_KSfbar (int n, double x);
 *
 * computes the complementary cumulative probability P[D_n >= x] of the
 * 2-sided 1-sample Kolmogorov-Smirnov distribution with sample size n at x.
 * It returns at least 10 decimal digits of precision for n <= 500,
 * at least 6 decimal digits of precision for 500 < n <= 200000,
 * and a few correct decimal digits for n > 200000.
 *
 */
double ncm_util_KSfbar (int n, double x);

/*
 * NOTE:
 * The ISO C99 function log1p of the standard math library does not exist in
 * ANSI C89. Here, it is programmed explicitly in KolmogorovSmirnovDist.c.

 * For ANSI C89 compilers, change the preprocessor condition to make it
 * available.
 */

void ncm_util_swilk (double *x, int n, double *w, double *pw, int *ifault);

#ifdef __cplusplus
}
#endif

#endif
