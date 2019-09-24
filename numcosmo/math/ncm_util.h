/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */
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

#ifndef NUMCOSMO_GIR_SCAN
#include <gsl/gsl_rng.h>
#include <gsl/gsl_min.h>
#include <complex.h>
#include <gmp.h>
#include <mpfr.h>
#endif /* NUMCOSMO_GIR_SCAN */

G_BEGIN_DECLS

gdouble *ncm_smoothd (gdouble *in, size_t N, size_t points, size_t pass);
gdouble ncm_topology_comoving_a0_lss (guint n, gdouble alpha);
gdouble ncm_topology_sigma_comoving_a0_lss (guint n, gdouble alpha, gdouble sigma_alpha);
gdouble ncm_sphPlm_x (gint l, gint m, gint order);
gdouble ncm_sum (gdouble *d, gulong n);

NCM_INLINE gdouble ncm_util_sqrt1px_m1 (const gdouble x) G_GNUC_CONST;
NCM_INLINE gdouble ncm_util_ln1pexpx (const gdouble x) G_GNUC_CONST;
NCM_INLINE gdouble ncm_util_1pcosx (const gdouble sinx, const gdouble cosx) G_GNUC_CONST;
NCM_INLINE gdouble ncm_util_1mcosx (const gdouble sinx, const gdouble cosx) G_GNUC_CONST;
NCM_INLINE gdouble ncm_util_1psinx (const gdouble sinx, const gdouble cosx) G_GNUC_CONST;
NCM_INLINE gdouble ncm_util_1msinx (const gdouble sinx, const gdouble cosx) G_GNUC_CONST;
NCM_INLINE gdouble ncm_util_cos2x (const gdouble sinx, const gdouble cosx) G_GNUC_CONST;

gdouble ncm_cmpdbl (const gdouble x, const gdouble y) G_GNUC_CONST;
gdouble ncm_exprel (const gdouble x) G_GNUC_CONST;
gdouble ncm_d1exprel (const gdouble x) G_GNUC_CONST;
gdouble ncm_d2exprel (const gdouble x) G_GNUC_CONST;
gdouble ncm_d3exprel (const gdouble x) G_GNUC_CONST;

gdouble ncm_util_sinh1 (const gdouble x) G_GNUC_CONST;
gdouble ncm_util_sinh3 (const gdouble x) G_GNUC_CONST;

gdouble ncm_util_sinhx_m_xcoshx_x3 (const gdouble x) G_GNUC_CONST;

gsize ncm_mpfr_out_raw (FILE *stream, mpfr_t op);
gsize ncm_mpfr_inp_raw (mpfr_t rop, FILE *stream);
gsize ncm_mpq_out_raw (FILE *f, mpq_t q);
gsize ncm_mpq_inp_raw (mpq_t q, FILE *f);

gulong ncm_random_seed (void);

gint ncm_cmp (gdouble x, gdouble y, const gdouble reltol, const gdouble abstol);

void ncm_rational_coarce_double (gdouble x, mpq_t q);
void ncm_mpz_inits (mpz_t z, ...) G_GNUC_NULL_TERMINATED;
void ncm_mpz_clears (mpz_t z, ...) G_GNUC_NULL_TERMINATED;
void _ncm_assertion_message_cmpdouble (const gchar *domain, const gchar *file, gint line, const gchar *func, const gchar *expr, gdouble arg1, const gchar *cmp, gdouble arg2, const gdouble reltol, const gdouble abstol);

gboolean ncm_util_cvode_check_flag (gpointer flagvalue, const gchar *funcname, gint opt);
gboolean ncm_util_cvode_print_stats (gpointer cvode);

gchar *ncm_util_basename_fits (const gchar *fits_filename); 
gchar *ncm_util_function_params (const gchar *func, gdouble **x, guint *len);
void ncm_util_print_bits (guint64 num);

gulong ncm_util_fact_size (const gulong n);

typedef struct _NcmComplex NcmComplex;

struct _NcmComplex 
{
  gdouble z[2];
};

GType ncm_complex_get_type (void) G_GNUC_CONST;

NcmComplex *ncm_complex_new (void);
NcmComplex *ncm_complex_ref (NcmComplex *c);
NcmComplex *ncm_complex_dup (NcmComplex *c);
void ncm_complex_free (NcmComplex *c);
void ncm_complex_clear (NcmComplex **c);

NCM_INLINE void ncm_complex_set (NcmComplex *c, const gdouble a, const gdouble b);
NCM_INLINE void ncm_complex_set_zero (NcmComplex *c);

NCM_INLINE gdouble ncm_complex_Re (NcmComplex *c);
NCM_INLINE gdouble ncm_complex_Im (NcmComplex *c);

NCM_INLINE void ncm_complex_res_add_mul_real (NcmComplex * restrict c1, const NcmComplex * restrict c2, const gdouble v);
NCM_INLINE void ncm_complex_res_add_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2, const NcmComplex * restrict c3);

NCM_INLINE void ncm_complex_mul_real (NcmComplex *c, const gdouble v);
NCM_INLINE void ncm_complex_res_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2);

NCM_INLINE gdouble ncm_util_smooth_trans (gdouble f0, gdouble f1, gdouble z0, gdouble dz, gdouble z);
NCM_INLINE void ncm_util_smooth_trans_get_theta (gdouble z0, gdouble dz, gdouble z, gdouble *theta0, gdouble *theta1);

NCM_INLINE gdouble ncm_util_position_angle (gdouble ra1, gdouble dec1, gdouble ra2, gdouble dec2);
NCM_INLINE gdouble ncm_util_great_circle_distance (gdouble ra1, gdouble dec1, gdouble ra2, gdouble dec2);

/* Macros */

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
  if (ncm_cmp ((n1), (n2), GSL_DBL_EPSILON, 0.0) cmp 0) ; else \
    _ncm_assertion_message_cmpdouble (G_LOG_DOMAIN, __FILE__, __LINE__, G_STRFUNC, \
                                      #n1 " " #cmp " " #n2, (n1), #cmp, (n2), GSL_DBL_EPSILON, 0.0); \
} while (0)

#define ncm_assert_cmpdouble_e(n1,cmp,n2,epsilon,abstol) \
do { \
  if (ncm_cmp ((n1), (n2), (epsilon), (abstol)) cmp 0) ; else \
    _ncm_assertion_message_cmpdouble (G_LOG_DOMAIN, __FILE__, __LINE__, G_STRFUNC, \
                                      #n1 " " #cmp " " #n2, (n1), #cmp, (n2), (epsilon), (abstol)); \
} while (0)

#define NCM_RETURN_IF_INF(a) if (gsl_isinf(a)) return a

#define NCM_FLOOR_TRUNC(a,b) (floor ((b) * (a)) / (b))
#define NCM_CEIL_TRUNC(a,b) (ceil ((b) * (a)) / (b))
#define NCM_ROUND_TRUNC(a,b) (round ((b) * (a)) / (b))

#define NCM_TEST_GSL_RESULT(func,ret) if (ret != GSL_SUCCESS) g_error ("%s: %s", func, gsl_strerror (ret))

#define NCM_COMPLEX_ZERO {{0.0, 0.0}}
#define NCM_COMPLEX(p) ((NcmComplex *)(p))
#define NCM_COMPLEX_PTR(p) ((NcmComplex **)(p))

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
#if !GLIB_CHECK_VERSION(2,40,0)
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
#endif /* !GLIB_CHECK_VERSION(2,40,0) */

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
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

NCM_INLINE gdouble
ncm_util_sqrt1px_m1 (const gdouble x)
{
  return x / (sqrt (1.0 + x) + 1.0);
}

NCM_INLINE gdouble 
ncm_util_ln1pexpx (const gdouble x)
{
	if (x > - GSL_LOG_DBL_EPSILON)
		return x;
	else
	{
		if (x > 0.0)
			return x + log1p (exp (-x));
		else
			return log1p (exp (x));
	}
}

NCM_INLINE gdouble 
ncm_util_1pcosx (const gdouble sinx, const gdouble cosx)
{
  if (cosx > -0.9)
    return 1.0 + cosx;
  else
    return sinx * sinx / (1.0 - cosx);
}

NCM_INLINE gdouble 
ncm_util_1mcosx (const gdouble sinx, const gdouble cosx)
{
  if (cosx < 0.9)
    return 1.0 - cosx;
  else
    return sinx * sinx / (1.0 + cosx);
}

NCM_INLINE gdouble 
ncm_util_1psinx (const gdouble sinx, const gdouble cosx)
{
  if (sinx > -0.9)
    return 1.0 + sinx;
  else
    return cosx * cosx / (1.0 - sinx);
}

NCM_INLINE gdouble 
ncm_util_1msinx (const gdouble sinx, const gdouble cosx)
{
  if (sinx < 0.9)
    return 1.0 - sinx;
  else
    return cosx * cosx / (1.0 + sinx);
}

NCM_INLINE gdouble 
ncm_util_cos2x (const gdouble sinx, const gdouble cosx)
{
  return (cosx - sinx) * (cosx + sinx);
}

NCM_INLINE gdouble 
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

NCM_INLINE void 
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

NCM_INLINE gdouble 
ncm_util_position_angle (gdouble ra1, gdouble dec1, gdouble ra2, gdouble dec2)
{
  const gdouble deg2rad  = M_PI / 180.0; 
	const gdouble ra1_rad  = ra1 * deg2rad;  
  const gdouble dec1_rad = dec1 * deg2rad;
	const gdouble ra2_rad  = ra2 * deg2rad;
  const gdouble dec2_rad = dec2 * deg2rad;
  const gdouble raDelta  = ra2_rad - ra1_rad;
  const gdouble theta    = atan2 (sin(raDelta), cos(dec1_rad) * tan(dec2_rad) - sin(dec1_rad) * cos(raDelta)); 

	return theta;
}

NCM_INLINE gdouble 
ncm_util_great_circle_distance (gdouble ra1, gdouble dec1, gdouble ra2, gdouble dec2)
{
	const gdouble deg2rad  = M_PI / 180.0;
	const gdouble phi1     = dec1 * deg2rad; 
  const gdouble lam1     = ra1 * deg2rad;
  const gdouble phi2     = dec2 * deg2rad;;
  const gdouble lam2     = ra2 * deg2rad;;
  const gdouble deltaLam = fabs(lam1 - lam2);

	const gdouble cosphi1  = cos (phi1);
	const gdouble sinphi1  = sin (phi1);
	const gdouble cosphi2  = cos (phi2);
	const gdouble sinphi2  = sin (phi2);
	
	const gdouble cosdeltaLam  = cos(deltaLam);
	const gdouble sindeltaLam  = sin(deltaLam);	
	
  const gdouble n1      = gsl_pow_2 (cosphi2 * sindeltaLam);
  const gdouble n2      = gsl_pow_2 (cosphi1 * sinphi2 - sinphi1 * cosphi2 * cosdeltaLam);
	const gdouble num     = sqrt (n1 + n2);
  const gdouble d1      = sinphi1 * sinphi2;
  const gdouble d2      = cosphi1 * cosphi2 * cosdeltaLam;
  const gdouble denom   = d1 + d2;

	return atan2(num, denom) / deg2rad; 
}

/* NcmComplex methods */

NCM_INLINE void 
ncm_complex_set (NcmComplex *c, const gdouble a, const gdouble b)
{
  c->z[0] = a;
  c->z[1] = b;
}

NCM_INLINE void 
ncm_complex_set_zero (NcmComplex *c)
{
  c->z[0] = c->z[1] = 0.0;
}

NCM_INLINE gdouble
ncm_complex_Re (NcmComplex *c)
{
  return c->z[0];
}

NCM_INLINE gdouble
ncm_complex_Im (NcmComplex *c)
{
  return c->z[1];
}

NCM_INLINE void
ncm_complex_res_add_mul_real (NcmComplex * restrict c1, const NcmComplex * restrict c2, const gdouble v)
{
  c1->z[0] += c2->z[0] * v;
  c1->z[1] += c2->z[1] * v;
}

NCM_INLINE void
ncm_complex_res_add_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2, const NcmComplex * restrict c3)
{
  c1->z[0] += c2->z[0] * c3->z[0] - c2->z[1] * c3->z[1];
  c1->z[1] += c2->z[0] * c3->z[1] + c2->z[1] * c3->z[0];
}

NCM_INLINE void 
ncm_complex_mul_real (NcmComplex *c, const gdouble v)
{
  c->z[0] *= v;
  c->z[1] *= v;  
}

NCM_INLINE void 
ncm_complex_res_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2)
{
  const gdouble Re_c1 = c1->z[0] * c2->z[0] - c1->z[1] * c1->z[1];
  c1->z[1] = c1->z[0] * c2->z[1] + c1->z[1] * c2->z[0];
  c1->z[0] = Re_c1;
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_UTIL_INLINE_H_ */

#ifndef KOLMOGOROVSMIRNOVDIST_H
#define KOLMOGOROVSMIRNOVDIST_H
#ifndef __GTK_DOC_IGNORE__

G_BEGIN_DECLS

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

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* KOLMOGOROVSMIRNOVDIST_H */
