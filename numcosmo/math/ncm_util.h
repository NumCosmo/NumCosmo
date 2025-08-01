/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 2; tab-width: 2 -*-  */

/***************************************************************************
 *            ncm_util.h
 *
 *  Mon Jul 16 18:02:22 2007
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

void ncm_util_mln_1mIexpzA_1pIexpmzA (const gdouble rho, const gdouble theta, const gdouble A, gdouble *rho1, gdouble *theta1);

gdouble ncm_util_normal_gaussian_integral (const gdouble xl, const gdouble xu);
gdouble ncm_util_gaussian_integral (const gdouble xl, const gdouble xu, const gdouble mu, const gdouble sigma);
gdouble ncm_util_log_normal_gaussian_integral (const gdouble xl, const gdouble xu, gdouble *sign);
gdouble ncm_util_log_gaussian_integral (const gdouble xl, const gdouble xu, const gdouble mu, const gdouble sigma, gdouble *sign);

gint ncm_cmp (gdouble x, gdouble y, const gdouble reltol, const gdouble abstol);

void ncm_rational_coarse_double (gdouble x, mpq_t q);
void ncm_mpz_inits (mpz_t z, ...) G_GNUC_NULL_TERMINATED;
void ncm_mpz_clears (mpz_t z, ...) G_GNUC_NULL_TERMINATED;
void _ncm_assertion_message_cmpdouble (const gchar *domain, const gchar *file, gint line, const gchar *func, const gchar *expr, gdouble arg1, const gchar *cmp, gdouble arg2, const gdouble reltol, const gdouble abstol);

gboolean ncm_util_cvode_check_flag (gpointer flagvalue, const gchar *funcname, gint opt);
gboolean ncm_util_cvode_print_stats (gpointer cvode);

gchar *ncm_util_basename_fits (const gchar *fits_filename);
gchar *ncm_util_function_params (const gchar *func, gdouble **x, guint *len);

gulong ncm_util_fact_size (const gulong n);

void ncm_util_sleep_ms (gint milliseconds);

void ncm_util_set_or_call_error (GError **error, GQuark domain, gint code, const gchar *format, ...);
void ncm_util_forward_or_call_error (GError **error, GError *local_error, const gchar *format, ...);


#define NCM_UTIL_ON_ERROR_RETURN(error, body, value) \
        G_STMT_START {                               \
          if ((error != NULL) && (*error) != NULL)   \
          {                                          \
            body;                                    \
            return value;                            \
          }                                          \
        } G_STMT_END

#define NCM_UTIL_ON_ERROR_FORWARD(error, body, value, prefix, ...) \
        G_STMT_START {                                             \
          if ((error != NULL) && (*error) != NULL)                 \
          {                                                        \
            g_prefix_error (error, prefix, ## __VA_ARGS__);        \
            body;                                                  \
            return value;                                          \
          }                                                        \
        } G_STMT_END

#ifndef NUMCOSMO_GIR_SCAN
typedef complex double NcmComplex;
#else /* NUMCOSMO_GIR_SCAN */
typedef struct _NcmComplexShouldNeverAppear NcmComplex;
#endif /* NUMCOSMO_GIR_SCAN */

GType ncm_complex_get_type (void) G_GNUC_CONST;

NcmComplex *ncm_complex_new (void);
NcmComplex *ncm_complex_dup (NcmComplex *c);
void ncm_complex_free (NcmComplex *c);
void ncm_complex_clear (NcmComplex **c);

NCM_INLINE void ncm_complex_set (NcmComplex *c, const gdouble a, const gdouble b);
NCM_INLINE void ncm_complex_set_zero (NcmComplex *c);

NCM_INLINE gdouble ncm_complex_Re (const NcmComplex *c);
NCM_INLINE gdouble ncm_complex_Im (const NcmComplex *c);

#ifndef NUMCOSMO_GIR_SCAN
NCM_INLINE void ncm_complex_set_c (NcmComplex *c, const complex double z);
NCM_INLINE complex double ncm_complex_c (const NcmComplex *c);

#endif /* NUMCOSMO_GIR_SCAN */

NCM_INLINE void ncm_complex_res_add_mul_real (NcmComplex * restrict c1, const NcmComplex * restrict c2, const gdouble v);
NCM_INLINE void ncm_complex_res_add_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2, const NcmComplex * restrict c3);

NCM_INLINE void ncm_complex_mul_real (NcmComplex *c, const gdouble v);
NCM_INLINE void ncm_complex_res_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2);

NCM_INLINE gdouble ncm_util_smooth_trans (gdouble f0, gdouble f1, gdouble z0, gdouble dz, gdouble z);
NCM_INLINE void ncm_util_smooth_trans_get_theta (gdouble z0, gdouble dz, gdouble z, gdouble *theta0, gdouble *theta1);

NCM_INLINE gdouble ncm_util_position_angle (gdouble ra1, gdouble dec1, gdouble ra2, gdouble dec2);
NCM_INLINE gdouble ncm_util_great_circle_distance (gdouble ra1, gdouble dec1, gdouble ra2, gdouble dec2);
NCM_INLINE gdouble ncm_util_projected_radius (gdouble theta, gdouble d);

/* Macros */

#define ncm_acb_get_complex(z) (arf_get_d (arb_midref (acb_realref (z)), ARF_RND_NEAR) + I * arf_get_d (arb_midref (acb_imagref (z)), ARF_RND_NEAR))

#define ncm_util_exp10(x) (exp ((x) * M_LN10))

#define NCM_GARRAY_MEMCPY(dest, src)                                                              \
        G_STMT_START {                                                                            \
          g_assert_cmpuint ((src)->len, ==, (dest)->len);                                         \
          g_assert_cmpuint (g_array_get_element_size (src), ==, g_array_get_element_size (dest)); \
          memcpy ((dest)->data, (src)->data, (src)->len * g_array_get_element_size (src));        \
        } G_STMT_END

#define NCM_GARRAY_DUP(dest, src)                                                              \
        G_STMT_START {                                                                         \
          dest = g_array_sized_new (FALSE, FALSE, g_array_get_element_size (src), (src)->len); \
          g_array_set_size ((dest), (src)->len);                                               \
          memcpy ((dest)->data, (src)->data, (src)->len * g_array_get_element_size (src));     \
        } G_STMT_END

#define ncm_assert_cmpdouble(n1, cmp, n2)                                                                  \
        do {                                                                                               \
          if (ncm_cmp ((n1), (n2), GSL_DBL_EPSILON, 0.0) cmp 0); else                                      \
          _ncm_assertion_message_cmpdouble (G_LOG_DOMAIN, __FILE__, __LINE__, G_STRFUNC,                   \
                                            #n1 " " #cmp " " #n2, (n1), #cmp, (n2), GSL_DBL_EPSILON, 0.0); \
        } while (0)

#define ncm_assert_cmpdouble_e(n1, cmp, n2, epsilon, abstol)                                              \
        do {                                                                                              \
          if (ncm_cmp ((n1), (n2), (epsilon), (abstol)) cmp 0); else                                      \
          _ncm_assertion_message_cmpdouble (G_LOG_DOMAIN, __FILE__, __LINE__, G_STRFUNC,                  \
                                            #n1 " " #cmp " " #n2, (n1), #cmp, (n2), (epsilon), (abstol)); \
        } while (0)

#define NCM_TEST_GSL_RESULT(func, ret) \
        if (ret != GSL_SUCCESS) g_error ("%s: %s", func, gsl_strerror (ret))

#define NCM_COMPLEX_ZERO (0.0)
#define NCM_COMPLEX(p) ((NcmComplex *) (p))
#define NCM_COMPLEX_PTR(p) ((NcmComplex **) (p))
#define NCM_COMPLEX_INIT(z) (z)
#define NCM_COMPLEX_INIT_REAL(z) (z)

#define ncm_g_string_clear(s)                      \
        G_STMT_START                               \
        if (*(s) != NULL)                          \
        {                                          \
          g_string_free (*(s), TRUE); *(s) = NULL; \
        }                                          \
        G_STMT_END

#define NCM_UNUSED(x) (void) (x)

void _ncm_util_set_destroyed (gpointer b);

#define NCM_TEST_FREE(cmd, obj)                                                                         \
        G_STMT_START {                                                                                  \
          gboolean destroyed = FALSE;                                                                   \
          g_object_set_data_full (G_OBJECT (obj), "test-destroy", &destroyed, _ncm_util_set_destroyed); \
          cmd (obj);                                                                                    \
          g_assert (destroyed);                                                                         \
        } G_STMT_END

#define NCM_TEST_FAIL(cmd)                       \
        G_STMT_START {                           \
          if (g_test_subprocess ())              \
          {                                      \
            cmd;                                 \
            exit (0);                            \
          }                                      \
          else                                   \
          {                                      \
            g_test_trap_subprocess (NULL, 0, 0); \
            g_test_trap_assert_failed ();        \
          }                                      \
        } G_STMT_END

#define NCM_TEST_PASS(cmd)                       \
        G_STMT_START {                           \
          if (g_test_subprocess ())              \
          {                                      \
            cmd;                                 \
            exit (0);                            \
          }                                      \
          else                                   \
          {                                      \
            g_test_trap_subprocess (NULL, 0, 0); \
            g_test_trap_assert_passed ();        \
          }                                      \
        } G_STMT_END


#define NCM_CVODE_CHECK(chk, name, val, ret)               \
        G_STMT_START {                                     \
          if (!ncm_util_cvode_check_flag (chk, name, val)) \
          return ret;                                      \
        }                                                  \
        G_STMT_END



/* Simple Callback macros */

/**
 * NCM_UTIL_DECLARE_CALLBACK:
 * @CallBack: The name of the callback structure in camel case
 * @CALL_BACK: The name of the callback function in uppercase
 * @callback: The name of the callback function in lowercase
 * @ret: The return type of the callback function
 * @args_decl: The declaration of the arguments of the callback function
 *
 * This macro declares a callback structure and the functions to handle it. You must use
 * NCM_UTIL_CALLBACK_ARGS to declare the arguments of the callback function.
 * The argument @arg_decl can be empty.
 *
 */
#define NCM_UTIL_CALLBACK_ARGS(...) , ## __VA_ARGS__
#define NCM_UTIL_DECLARE_CALLBACK(CallBack, CALL_BACK, callback, ret, args_decl)         \
        typedef struct _ ## CallBack CallBack;                                           \
        G_GNUC_UNUSED static inline CallBack *CALL_BACK (gpointer callback_ptr) {        \
          return callback_ptr;                                                           \
        }                                                                                \
        typedef ret (*CallBack ## Func) (gpointer callback_data args_decl);              \
        typedef gpointer (*CallBack ## CopyData) (gpointer callback_data);               \
        typedef void (*CallBack ## FreeData) (gpointer callback_data);                   \
        typedef void (*CallBack ## PrepareData) (gpointer callback_data, NcmMSet *mset); \
        struct _ ## CallBack                                                             \
        {                                                                                \
          /*< private >*/                                                                \
          CallBack ## Func func;                                                         \
          CallBack ## FreeData callback_data_free;                                       \
          CallBack ## CopyData callback_data_copy;                                       \
          CallBack ## PrepareData callback_data_prepare;                                 \
          gpointer callback_data;                                                        \
        };                                                                               \
        GType callback ## _get_type (void) G_GNUC_CONST;                                 \
        CallBack *callback ## _new (CallBack ## Func func,                               \
                                    CallBack ## FreeData callback_data_free,             \
                                    CallBack ## CopyData callback_data_copy,             \
                                    CallBack ## PrepareData callback_data_prepare,       \
                                    gpointer callback_data);                             \
        CallBack *callback ## _copy (CallBack * callback_obj);                           \
        void callback ## _free (CallBack * callback_obj);                                \
        ret callback ## _eval (CallBack * callback_obj args_decl);                       \
        void callback ## _prepare (CallBack * callback_obj, NcmMSet * mset);

/**
 * NCM_UTIL_DEFINE_CALLBACK:
 * @CallBack: The name of the callback structure in camel case
 * @CALL_BACK: The name of the callback function in uppercase
 * @callback: The name of the callback function in lowercase
 * @ret: The return type of the callback function
 * @args_decl: The declaration of the arguments of the callback function
 * @args: The arguments of the callback function
 *
 * This macro defines the functions to handle the callback structure. You must use
 * NCM_UTIL_CALLBACK_ARGS to declare both the arguments declaration of the callback
 * function and the arguments of the function. They can be empty.
 *
 */
#define NCM_UTIL_DEFINE_CALLBACK(CallBack, CALL_BACK, callback, ret, args_decl, args)                               \
        G_DEFINE_BOXED_TYPE (CallBack, callback, callback ## _copy, callback ## _free)                              \
        CallBack *callback ## _new (CallBack ## Func func,                                                          \
                                    CallBack ## FreeData callback_data_free,                                        \
                                    CallBack ## CopyData callback_data_copy,                                        \
                                    CallBack ## PrepareData callback_data_prepare,                                  \
                                    gpointer callback_data)                                                         \
        {                                                                                                           \
          CallBack *callback_obj = g_new0 (CallBack, 1);                                                            \
          g_assert_nonnull (func);                                                                                  \
          g_assert_nonnull (callback_data_free);                                                                    \
          g_assert_nonnull (callback_data_copy);                                                                    \
          callback_obj->func                  = func;                                                               \
          callback_obj->callback_data_free    = callback_data_free;                                                 \
          callback_obj->callback_data_copy    = callback_data_copy;                                                 \
          callback_obj->callback_data         = callback_data;                                                      \
          callback_obj->callback_data_prepare = callback_data_prepare;                                              \
          return callback_obj;                                                                                      \
        }                                                                                                           \
        CallBack *                                                                                                  \
        callback ## _copy (CallBack * callback_obj)                                                                 \
        {                                                                                                           \
          CallBack *new_callback_obj = g_new0 (CallBack, 1);                                                        \
          new_callback_obj->func                  = callback_obj->func;                                             \
          new_callback_obj->callback_data_free    = callback_obj->callback_data_free;                               \
          new_callback_obj->callback_data_copy    = callback_obj->callback_data_copy;                               \
          new_callback_obj->callback_data_prepare = callback_obj->callback_data_prepare;                            \
          new_callback_obj->callback_data         = callback_obj->callback_data_copy (callback_obj->callback_data); \
          return new_callback_obj;                                                                                  \
        }                                                                                                           \
        void                                                                                                        \
        callback ## _free (CallBack * callback_obj)                                                                 \
        {                                                                                                           \
          callback_obj->callback_data_free (callback_obj->callback_data);                                           \
          g_free (callback_obj);                                                                                    \
        }                                                                                                           \
        ret                                                                                                         \
        callback ## _eval (CallBack * callback_obj args_decl)                                                       \
        {                                                                                                           \
          return callback_obj->func (callback_obj->callback_data args);                                             \
        }                                                                                                           \
        void                                                                                                        \
        callback ## _prepare (CallBack * callback_obj, NcmMSet * mset)                                              \
        {                                                                                                           \
          if (callback_obj->callback_data_prepare != NULL)                                                          \
          callback_obj->callback_data_prepare (callback_obj->callback_data, mset);                                  \
        }


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
  if (x > -GSL_LOG_DBL_EPSILON)
  {
    return x;
  }
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
  const gdouble a       = -z0 - 0.5 * dz;
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
  const gdouble a       = -z0 - 0.5 * dz;
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
  const gdouble theta    = atan2 (sin (raDelta), cos (dec1_rad) * tan (dec2_rad) - sin (dec1_rad) * cos (raDelta));

  return theta;
}

NCM_INLINE gdouble
ncm_util_great_circle_distance (gdouble ra1, gdouble dec1, gdouble ra2, gdouble dec2)
{
  const gdouble deg2rad     = M_PI / 180.0;
  const gdouble phi1        = dec1 * deg2rad;
  const gdouble lam1        = ra1 * deg2rad;
  const gdouble phi2        = dec2 * deg2rad;
  const gdouble lam2        = ra2 * deg2rad;
  const gdouble deltaLam    = fabs (lam1 - lam2);
  const gdouble cosphi1     = cos (phi1);
  const gdouble sinphi1     = sin (phi1);
  const gdouble cosphi2     = cos (phi2);
  const gdouble sinphi2     = sin (phi2);
  const gdouble cosdeltaLam = cos (deltaLam);
  const gdouble sindeltaLam = sin (deltaLam);
  const gdouble n1          = gsl_pow_2 (cosphi2 * sindeltaLam);
  const gdouble n2          = gsl_pow_2 (cosphi1 * sinphi2 - sinphi1 * cosphi2 * cosdeltaLam);
  const gdouble num         = sqrt (n1 + n2);
  const gdouble d1          = sinphi1 * sinphi2;
  const gdouble d2          = cosphi1 * cosphi2 * cosdeltaLam;
  const gdouble denom       = d1 + d2;

  return atan2 (num, denom) / deg2rad;
}

NCM_INLINE gdouble
ncm_util_projected_radius (gdouble theta, gdouble d)
{
  return d * sin (theta);
}

/* NcmComplex methods */

NCM_INLINE void
ncm_complex_set (NcmComplex *c, const gdouble a, const gdouble b)
{
  *c = a + I * b;
}

NCM_INLINE void
ncm_complex_set_zero (NcmComplex *c)
{
  *c = 0.0;
}

NCM_INLINE gdouble
ncm_complex_Re (const NcmComplex *c)
{
  return creal (*c);
}

NCM_INLINE gdouble
ncm_complex_Im (const NcmComplex *c)
{
  return cimag (*c);
}

#ifndef NUMCOSMO_GIR_SCAN

NCM_INLINE void
ncm_complex_set_c (NcmComplex *c, const complex double z)
{
  *c = z;
}

NCM_INLINE complex double
ncm_complex_c (const NcmComplex *c)
{
  return *c;
}

#endif /* NUMCOSMO_GIR_SCAN */

NCM_INLINE void
ncm_complex_res_add_mul_real (NcmComplex * restrict c1, const NcmComplex * restrict c2, const gdouble v)
{
  *c1 += (*c2) * v;
}

NCM_INLINE void
ncm_complex_res_add_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2, const NcmComplex * restrict c3)
{
  *c1 += (*c2) * (*c3);
}

NCM_INLINE void
ncm_complex_mul_real (NcmComplex *c, const gdouble v)
{
  *c *= v;
}

NCM_INLINE void
ncm_complex_res_mul (NcmComplex * restrict c1, const NcmComplex * restrict c2)
{
  *c1 *= *c2;
}

G_END_DECLS

#endif /* __GTK_DOC_IGNORE__ */
#endif /* NUMCOSMO_HAVE_INLINE */
#endif /* _NCM_UTIL_INLINE_H_ */

