/***************************************************************************
 *            ncm_spectral.c
 *
 *  Tue Feb 04 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/

/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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
 * SECTION:ncm_spectral
 * @title: NcmSpectral
 * @short_description: Spectral methods for function approximation
 *
 * Provides spectral methods for function approximation using Chebyshev
 * and Gegenbauer polynomials.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_spectral.h"
#include "math/ncm_cfg.h"

#include <math.h>
#ifndef NUMCOSMO_GIR_SCAN
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

struct _NcmSpectral
{
  /*< private >*/
  GObject parent_instance;
  /* Chebyshev coefficients computation cache */
  guint cheb_N_cached;     /* Cached N value */
  gdouble *cheb_f_vals;    /* Cached function values array */
  fftw_plan cheb_plan_r2r; /* Cached FFTW plan */
};

G_DEFINE_TYPE (NcmSpectral, ncm_spectral, G_TYPE_OBJECT)

static void
ncm_spectral_init (NcmSpectral *spectral)
{
  spectral->cheb_N_cached = 0;
  spectral->cheb_f_vals   = NULL;
  spectral->cheb_plan_r2r = NULL;
}

static void
ncm_spectral_finalize (GObject *object)
{
  NcmSpectral *spectral = NCM_SPECTRAL (object);

  if (spectral->cheb_plan_r2r != NULL)
  {
    fftw_destroy_plan (spectral->cheb_plan_r2r);
    spectral->cheb_plan_r2r = NULL;
  }

  if (spectral->cheb_f_vals != NULL)
  {
    fftw_free (spectral->cheb_f_vals);
    spectral->cheb_f_vals = NULL;
  }

  G_OBJECT_CLASS (ncm_spectral_parent_class)->finalize (object);
}

static void
ncm_spectral_class_init (NcmSpectralClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize = ncm_spectral_finalize;
}

/**
 * ncm_spectral_new:
 *
 * Creates a new #NcmSpectral object.
 *
 * Returns: (transfer full): a new #NcmSpectral
 */
NcmSpectral *
ncm_spectral_new (void)
{
  return g_object_new (NCM_TYPE_SPECTRAL, NULL);
}

/**
 * ncm_spectral_ref:
 * @spectral: a #NcmSpectral
 *
 * Increases the reference count of @spectral by one.
 *
 * Returns: (transfer full): @spectral
 */
NcmSpectral *
ncm_spectral_ref (NcmSpectral *spectral)
{
  return g_object_ref (spectral);
}

/**
 * ncm_spectral_free:
 * @spectral: a #NcmSpectral
 *
 * Decreases the reference count of @spectral by one. If the reference count
 * reaches zero, the object is freed.
 */
void
ncm_spectral_free (NcmSpectral *spectral)
{
  g_object_unref (spectral);
}

/**
 * ncm_spectral_clear:
 * @spectral: a #NcmSpectral
 *
 * If @spectral is not NULL, decreases the reference count of @spectral by one
 * and sets the pointer to NULL.
 */
void
ncm_spectral_clear (NcmSpectral **spectral)
{
  g_clear_object (spectral);
}

/**
 * ncm_spectral_compute_chebyshev_coeffs:
 * @spectral: a #NcmSpectral
 * @F: (scope call): function to evaluate
 * @a: left endpoint
 * @b: right endpoint
 * @coeffs: output vector
 * @user_data: user data for @F
 *
 * Computes Chebyshev coefficients of f(x) on [a,b] using FFTW DCT-I. The function is
 * sampled at Chebyshev nodes $x_k = (a+b)/2 - (b-a)/2\cos(k\pi/(N-1))$ and transformed
 * using a Type-I discrete cosine transform.
 *
 * The user must provide a pre-allocated #NcmVector @coeffs with length N, where N is
 * the number of Chebyshev coefficients to compute. The coefficients are stored in
 * @coeffs in increasing order, i.e., coeffs[0] corresponds to T_0, coeffs[1] to T_1,
 * ..., coeffs[N-1] to T_{N-1}. This method caches the FFTW plan and working arrays for
 * efficiency. They are only reallocated if N changes.
 *
 */
void
ncm_spectral_compute_chebyshev_coeffs (NcmSpectral *spectral, NcmSpectralF F, gdouble a, gdouble b, NcmVector *coeffs, gpointer user_data)
{
  const gdouble mid    = 0.5 * (a + b);
  const gdouble half_h = 0.5 * (b - a);
  const guint N        = ncm_vector_len (coeffs);
  guint i;

  /* Reallocate and replan if N has changed */
  if (spectral->cheb_N_cached != N)
  {
    /* Clean up old resources */
    if (spectral->cheb_plan_r2r != NULL)
    {
      fftw_destroy_plan (spectral->cheb_plan_r2r);
      spectral->cheb_plan_r2r = NULL;
    }

    if (spectral->cheb_f_vals != NULL)
    {
      fftw_free (spectral->cheb_f_vals);
      spectral->cheb_f_vals = NULL;
    }

    /* Allocate new resources */
    spectral->cheb_f_vals = fftw_malloc (sizeof (gdouble) * N);

    /* Create new FFTW plan */
    ncm_cfg_load_fftw_wisdom ("ncm_sbessel_ode_solver");
    ncm_cfg_lock_plan_fftw ();
    spectral->cheb_plan_r2r = fftw_plan_r2r_1d (N, spectral->cheb_f_vals, ncm_vector_data (coeffs),
                                                FFTW_REDFT00, ncm_cfg_get_fftw_default_flag ());
    ncm_cfg_unlock_plan_fftw ();
    ncm_cfg_save_fftw_wisdom ("ncm_sbessel_ode_solver");

    spectral->cheb_N_cached = N;
  }

  /* Sample function at Chebyshev nodes */
  for (i = 0; i < N; i++)
  {
    const gdouble theta = M_PI * i / (N - 1);
    const gdouble x     = mid + half_h * cos (theta);

    spectral->cheb_f_vals[i] = F (user_data, x);
  }

  /* Execute FFTW plan (need to update output pointer for this execution) */
  fftw_execute_r2r (spectral->cheb_plan_r2r, spectral->cheb_f_vals, ncm_vector_data (coeffs));

  /* Normalize coefficients */
  g_assert_cmpuint (ncm_vector_stride (coeffs), ==, 1);
  {
    gdouble *coeffs_data = ncm_vector_data (coeffs);

    coeffs_data[0]     = coeffs_data[0] / ((N - 1.0) * 2.0);
    coeffs_data[N - 1] = coeffs_data[N - 1] / ((N - 1.0) * 2.0);

    for (i = 1; i < N - 1; i++)
      coeffs_data[i] = coeffs_data[i] / (N - 1.0);
  }
}

/**
 * ncm_spectral_chebT_to_gegenbauer_alpha1:
 * @c: Chebyshev coefficients vector
 * @g: Gegenbauer $C^{(1)}_n$ coefficients vector (must have same length as @c, pre-allocated by caller)
 *
 * Converts Chebyshev $T_n$ coefficients to Gegenbauer $C^{(1)}_n$ coefficients ($\alpha=1$).
 * Uses the relationship: $T_n = \frac{1}{2}(C^{(1)}_n + C^{(1)}_{n-2})$ for $n \geq 2$.
 */
void
ncm_spectral_chebT_to_gegenbauer_alpha1 (NcmVector *c, NcmVector *g)
{
  const guint N = ncm_vector_len (c);
  guint i;

  g_assert_cmpuint (ncm_vector_len (g), ==, N);
  g_assert_cmpuint (ncm_vector_stride (c), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (g), ==, 1);

  if (N == 0)
    return;

  {
    const gdouble *c_data = ncm_vector_data (c);
    gdouble *g_data       = ncm_vector_data (g);

    memset (g_data, 0, N * sizeof (gdouble));

    /* n = 0 case */
    g_data[0] = c_data[0];

    if (N == 1)
      return;

    /* n = 1 case */
    g_data[1] = c_data[1] * 0.5;

    /* n >= 2 */
    for (i = 2; i < N; i++)
    {
      const gdouble ci = c_data[i];

      g_data[i]     += 0.5 * ci;
      g_data[i - 2] -= 0.5 * ci;
    }
  }
}

/**
 * ncm_spectral_chebT_to_gegenbauer_alpha2:
 * @c: Chebyshev coefficients vector
 * @g: Gegenbauer $C^{(2)}_k$ coefficients vector (must have same length as @c, pre-allocated by caller)
 *
 * Converts Chebyshev $T_n$ coefficients to Gegenbauer $C^{(2)}_k$ coefficients ($\alpha=2$).
 *
 * Uses the projection formula:
 * $$g_k = \frac{1}{2} c_0 \delta_{k,0} + \frac{c_k}{2(k+1)} - \frac{(k+2) c_{k+2}}{(k+1)(k+3)} + \frac{c_{k+4}}{2(k+3)}$$
 * where $f(x) = \sum_n c_n T_n(x)$.
 */
void
ncm_spectral_chebT_to_gegenbauer_alpha2 (NcmVector *c, NcmVector *g)
{
  const guint N = ncm_vector_len (c);
  guint k;

  g_assert_cmpuint (ncm_vector_len (g), ==, N);
  g_assert_cmpuint (ncm_vector_stride (c), ==, 1);
  g_assert_cmpuint (ncm_vector_stride (g), ==, 1);

  if (N == 0)
    return;

  {
    const gdouble *c_data = ncm_vector_data (c);
    gdouble *g_data       = ncm_vector_data (g);

    /* Zero output vector */
    memset (g_data, 0, N * sizeof (gdouble));

    /* Apply projection formula for each k */
    for (k = 0; k < N; k++)
    {
      const gdouble kd = (gdouble) k;
      gdouble gk       = 0.0;

      /* Special case: k=0 has additional 1/2 * c[0] contribution */
      if (k == 0)
        gk += 0.5 * c_data[0];

      /* First term: c[k] / (2*(k+1)) */
      gk += c_data[k] / (2.0 * (kd + 1.0));

      /* Second term: -(k+2) * c[k+2] / ((k+1)*(k+3)) */
      if (k + 2 < N)
        gk -= (kd + 2.0) * c_data[k + 2] / ((kd + 1.0) * (kd + 3.0));

      /* Third term: c[k+4] / (2*(k+3)) */
      if (k + 4 < N)
        gk += c_data[k + 4] / (2.0 * (kd + 3.0));

      g_data[k] = gk;
    }
  }
}

/**
 * ncm_spectral_gegenbauer_alpha1_eval:
 * @c: Gegenbauer $C^{(1)}_n$ coefficients vector
 * @x: point to evaluate
 *
 * Evaluates a Gegenbauer $C^{(1)}_n$ expansion at x using Clenshaw recurrence.
 * For $\alpha=1$, $C^{(1)}_n(x) = U_n(x)$ (Chebyshev polynomials of the second kind).
 *
 * Returns: the value of $\sum_{n=0}^{N-1} c_n C^{(1)}_n(x)$
 */
gdouble
ncm_spectral_gegenbauer_alpha1_eval (NcmVector *c, gdouble x)
{
  const guint N = ncm_vector_len (c);

  if (N == 0)
    return 0.0;

  g_assert_cmpuint (ncm_vector_stride (c), ==, 1);

  {
    const gdouble *c_data = ncm_vector_data (c);

    /* Endpoint handling: C_n^{(1)}(+/-1) = (n+1)*(+/-1)^n */
    if (fabs (x - 1.0) < 1e-15)
    {
      gdouble sum = 0.0;
      guint n;

      for (n = 0; n < N; n++)
        sum += c_data[n] * (gdouble) (n + 1);

      return sum;
    }

    if (fabs (x + 1.0) < 1e-15)
    {
      gdouble sum = 0.0;
      guint n;

      for (n = 0; n < N; n++)
        sum += c_data[n] * ((n & 1) ? -(gdouble) (n + 1) : (gdouble) (n + 1));

      return sum;
    }

    {
      /* Stable recurrence for interior x */
      gdouble Cnm1 = 1.0; /* U_0 */
      gdouble sum  = c_data[0] * Cnm1;

      if (N == 1)
        return sum;

      gdouble Cn = 2.0 * x; /* U_1 */

      sum += c_data[1] * Cn;

      for (guint n = 1; n < N - 1; n++)
      {
        gdouble Cnp1 = 2.0 * x * Cn - Cnm1; /* U_{n+1} */

        sum += c_data[n + 1] * Cnp1;
        Cnm1 = Cn;
        Cn   = Cnp1;
      }

      return sum;
    }
  }
}

/**
 * ncm_spectral_gegenbauer_alpha2_eval:
 * @c: Gegenbauer $C^{(2)}_n$ coefficients vector
 * @x: point to evaluate
 *
 * Evaluates a Gegenbauer $C^{(2)}_n$ expansion at x using Clenshaw recurrence.
 * For $\alpha=2$, the recurrence relation is:
 * $(n+1) C^{(2)}_{n+1}(x) = 2(n+2)x C^{(2)}_n(x) - (n+3) C^{(2)}_{n-1}(x)$
 *
 * Returns: the value of $\sum_{n=0}^{N-1} c_n C^{(2)}_n(x)$
 */
gdouble
ncm_spectral_gegenbauer_alpha2_eval (NcmVector *c, gdouble x)
{
  const guint N = ncm_vector_len (c);

  if (N == 0)
    return 0.0;

  g_assert_cmpuint (ncm_vector_stride (c), ==, 1);

  {
    const gdouble *c_data = ncm_vector_data (c);

    /* Endpoint handling: C_n^{(2)}(+/-1) = binom(n+3,3)*(+/-1)^n = ((n+1)*(n+2)*(n+3)/6)*(+/-1)^n */
    if (fabs (x - 1.0) < 1e-15)
    {
      gdouble sum = 0.0;
      guint n;

      for (n = 0; n < N; n++)
        sum += c_data[n] * (gdouble) ((n + 1) * (n + 2) * (n + 3)) / 6.0;

      return sum;
    }

    if (fabs (x + 1.0) < 1e-15)
    {
      gdouble sum = 0.0;
      guint n;

      for (n = 0; n < N; n++)
      {
        const gdouble val = (gdouble) ((n + 1) * (n + 2) * (n + 3)) / 6.0;

        sum += c_data[n] * ((n & 1) ? -val : val);
      }

      return sum;
    }

    {
      /* Stable recurrence for interior x */
      gdouble Cnm1 = 1.0; /* C_0^{(2)} = 1 */
      gdouble sum  = c_data[0] * Cnm1;

      if (N == 1)
        return sum;

      gdouble Cn = 4.0 * x; /* C_1^{(2)} = 4x */

      sum += c_data[1] * Cn;

      for (guint n = 1; n < N - 1; n++)
      {
        /* (n+1) C_{n+1}^{(2)} = 2(n+2)x C_n^{(2)} - (n+3) C_{n-1}^{(2)} */
        gdouble Cnp1 = (2.0 * (gdouble) (n + 2) * x * Cn - (gdouble) (n + 3) * Cnm1) / (gdouble) (n + 1);

        sum += c_data[n + 1] * Cnp1;
        Cnm1 = Cn;
        Cn   = Cnp1;
      }

      return sum;
    }
  }
}

/**
 * ncm_spectral_chebyshev_eval:
 * @a: Chebyshev coefficients vector
 * @t: point to evaluate in [-1,1]
 *
 * Evaluates a Chebyshev expansion $f(t) = \sum_{k=0}^{N-1} a_k T_k(t)$
 * using Clenshaw recurrence with Reinsch modification near endpoints.
 *
 * Returns: the value of the Chebyshev expansion at t
 */
gdouble
ncm_spectral_chebyshev_eval (NcmVector *a, gdouble t)
{
  const guint N           = ncm_vector_len (a);
  const gdouble *a_data   = ncm_vector_data (a);
  const gdouble threshold = 0.9;
  const gdouble eps       = 1e-15;

  if (N == 0)
    return 0.0;

  g_assert_cmpuint (ncm_vector_stride (a), ==, 1);

  if (N == 1)
    return a_data[0];

  /* Endpoint handling: T_k(+1) = 1, T_k(-1) = (-1)^k */
  if (fabs (t - 1.0) < eps)
  {
    gdouble sum = 0.0;
    guint k;

    for (k = 0; k < N; k++)
      sum += a_data[k];

    return sum;
  }

  if (fabs (t + 1.0) < eps)
  {
    gdouble sum = 0.0;
    guint k;

    for (k = 0; k < N; k++)
      sum += ((k & 1) ? -a_data[k] : a_data[k]);

    return sum;
  }

  if (fabs (t) < threshold)
  {
    /* Clenshaw recurrence for interior points */
    gdouble b_kplus1 = 0.0;
    gdouble b_kplus2 = 0.0;
    gdouble two_t    = t + t;
    gint k;

    for (k = (gint) N - 1; k >= 1; k--)
    {
      gdouble b_k = two_t * b_kplus1 - b_kplus2 + a_data[k];

      b_kplus2 = b_kplus1;
      b_kplus1 = b_k;
    }

    return t * b_kplus1 - b_kplus2 + a_data[0];
  }

  /* Near +1 : Reinsch modification */
  if (t > 0.0)
  {
    gdouble d_kplus1      = 0.0;
    gdouble e_kplus1      = 0.0;
    const gdouble tm1     = (t - 0.5) - 0.5;
    const gdouble two_tm1 = tm1 + tm1;

    for (gint k = (gint) N - 1; k >= 1; k--)
    {
      gdouble d_k = two_tm1 * e_kplus1 + d_kplus1 + a_data[k];
      gdouble e_k = d_k + e_kplus1;

      d_kplus1 = d_k;
      e_kplus1 = e_k;
    }

    return tm1 * e_kplus1 + d_kplus1 + a_data[0];
  }

  /* Near -1 : Reinsch modification */
  {
    gdouble d_kplus1      = 0.0;
    gdouble e_kplus1      = 0.0;
    const gdouble tp1     = (t + 0.5) + 0.5;
    const gdouble two_tp1 = tp1 + tp1;

    for (gint k = (gint) N - 1; k >= 1; k--)
    {
      gdouble d_k = two_tp1 * e_kplus1 - d_kplus1 + a_data[k];
      gdouble e_k = d_k - e_kplus1;

      d_kplus1 = d_k;
      e_kplus1 = e_k;
    }

    return tp1 * e_kplus1 - d_kplus1 + a_data[0];
  }
}

/**
 * ncm_spectral_chebyshev_deriv:
 * @a: Chebyshev coefficients vector (a_j multiplies T_j)
 * @t: point to evaluate in [-1,1]
 *
 * Evaluates the first derivative of a Chebyshev expansion at $t$.
 *
 * The Chebyshev series is
 * $$ f(t) = \sum_{j=0}^{N-1} a_j T_j(t) , $$
 * and its derivative can be written as
 * $$ f'(t) = \sum_{k=0}^{N-2} b_k T_k(t) . $$
 *
 * The derivative coefficients $b_k$ satisfy
 * $$ b_k = \sum_{j=k+1,k+3,\dots}^{N-1} 2 j a_j , \quad k \ge 1, $$
 * and
 * $$ b_0 = \sum_{j=1,3,5,\dots}^{N-1} j a_j . $$
 *
 * The derivative is evaluated using a fused backward recurrence and
 * the Clenshaw algorithm, without explicitly forming the coefficients $b_k$.
 *
 * Returns: the value of the derivative at $t$
 */
gdouble
ncm_spectral_chebyshev_deriv (NcmVector *a, gdouble t)
{
  const gint N = ncm_vector_len (a);

  if (N <= 1)
    return 0.0;

  g_assert_cmpuint (ncm_vector_stride (a), ==, 1);

  {
    const gdouble *a_data = ncm_vector_data (a);

    if (N == 2)
      return a_data[1];


    if (fabs (t - 1.0) < 1.0e-15)
    {
      /* ---- x = +1 ---- */
      gdouble d1  = 0.0;
      gdouble d2  = 0.0;
      gdouble sum = 0.0;

      for (gint k = N - 2; k >= 1; k--)
      {
        gdouble bk = d2 + 2.0 * (k + 1) * a_data[k + 1];

        d2   = d1;
        d1   = bk;
        sum += bk;
      }

      /* b0 has the 1/2 factor */
      gdouble b0 = 0.5 * (d2 + 2.0 * a_data[1]);

      sum += b0;

      return sum;
    }

    if (fabs (t + 1.0) < 1.0e-15)
    {
      /* ---- x = -1 ---- */
      gdouble d1  = 0.0;
      gdouble d2  = 0.0;
      gdouble sum = 0.0;

      /* start with (-1)^(N-2) */
      gdouble sign = ((N - 2) & 1) ? -1.0 : 1.0;

      for (gint k = N - 2; k >= 1; k--)
      {
        gdouble bk = d2 + 2.0 * (k + 1) * a_data[k + 1];

        d2 = d1;
        d1 = bk;

        sum += sign * bk;
        sign = -sign;
      }

      /* b0 has sign +1 */
      gdouble b0 = 0.5 * (d2 + 2.0 * a_data[1]);

      sum += b0;

      return sum;
    }

    {
      gdouble c1          = 0.0; /* Clenshaw state k+1 */
      gdouble c2          = 0.0; /* Clenshaw state k+2 */
      gdouble d1          = 0.0; /* recurrence helper */
      gdouble d2          = 0.0;
      const gdouble two_t = 2.0 * t;

      /* k = N-2 ... 1 */
      for (gint k = N - 2; k >= 1; k--)
      {
        /* build b[k] on the fly */
        gdouble bk = d2 + 2.0 * (k + 1) * a_data[k + 1];

        /* update derivative recurrence */
        d2 = d1;
        d1 = bk;

        /* Clenshaw step */
        gdouble c0 = two_t * c1 - c2 + bk;

        c2 = c1;
        c1 = c0;
      }

      {
        /* k = 0 needs the 1/2 factor */
        gdouble b0 = 0.5 * (d2 + 2.0 * a_data[1]);

        return t * c1 - c2 + b0;
      }
    }
  }
}

