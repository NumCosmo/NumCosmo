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

enum
{
  PROP_0,
  PROP_MAX_ORDER,
};

struct _NcmSpectral
{
  /*< private >*/
  GObject parent_instance;

  /* Adaptive refinement fields */
  guint max_order;       /* Maximum k: N_max = 2^max_order + 1 */
  gdouble *f_vals;       /* Function values array, size: 2^max_order + 1 */
  gdouble *coeffs_work;  /* Coefficients work array, size: 2^max_order + 1 */
  GArray *coeffs;        /* Coefficients array, size: 2^max_order + 1 */
  GPtrArray *cos_arrays; /* Precomputed cosines for each k level */
  GPtrArray *fftw_plans; /* FFTW plans for each k level */

  /* Legacy fields (backward compatibility) */
  guint cheb_N_cached;     /* Cached N value */
  gdouble *cheb_f_vals;    /* Cached function values array */
  gdouble *cheb_cos_vals;  /* Cached cosine values at Chebyshev nodes */
  fftw_plan cheb_plan_r2r; /* Cached FFTW plan */
};

G_DEFINE_TYPE (NcmSpectral, ncm_spectral, G_TYPE_OBJECT)

static void
ncm_spectral_init (NcmSpectral *spectral)
{
  spectral->max_order   = 16; /* Default: N_max = 1025 */
  spectral->f_vals      = NULL;
  spectral->coeffs_work = NULL;
  spectral->coeffs      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  spectral->cos_arrays  = NULL;
  spectral->fftw_plans  = NULL;

  spectral->cheb_N_cached = 0;
  spectral->cheb_f_vals   = NULL;
  spectral->cheb_cos_vals = NULL;
  spectral->cheb_plan_r2r = NULL;
}

static void
ncm_spectral_finalize (GObject *object)
{
  NcmSpectral *spectral = NCM_SPECTRAL (object);

  /* Clean up adaptive refinement resources */
  if (spectral->fftw_plans != NULL)
  {
    g_ptr_array_unref (spectral->fftw_plans);
    spectral->fftw_plans = NULL;
  }

  if (spectral->cos_arrays != NULL)
  {
    g_ptr_array_unref (spectral->cos_arrays);
    spectral->cos_arrays = NULL;
  }

  if (spectral->f_vals != NULL)
  {
    fftw_free (spectral->f_vals);
    spectral->f_vals = NULL;
  }

  if (spectral->coeffs_work != NULL)
  {
    fftw_free (spectral->coeffs_work);
    spectral->coeffs_work = NULL;
  }

  g_clear_pointer (&spectral->coeffs, g_array_unref);

  /* Clean up legacy resources */
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

  if (spectral->cheb_cos_vals != NULL)
  {
    g_free (spectral->cheb_cos_vals);
    spectral->cheb_cos_vals = NULL;
  }

  G_OBJECT_CLASS (ncm_spectral_parent_class)->finalize (object);
}

static void
ncm_spectral_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmSpectral *spectral = NCM_SPECTRAL (object);

  g_return_if_fail (NCM_IS_SPECTRAL (object));

  switch (prop_id)
  {
    case PROP_MAX_ORDER:
      ncm_spectral_set_max_order (spectral, g_value_get_uint (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_spectral_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmSpectral *spectral = NCM_SPECTRAL (object);

  g_return_if_fail (NCM_IS_SPECTRAL (object));

  switch (prop_id)
  {
    case PROP_MAX_ORDER:
      g_value_set_uint (value, ncm_spectral_get_max_order (spectral));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
ncm_spectral_class_init (NcmSpectralClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize     = ncm_spectral_finalize;
  object_class->set_property = ncm_spectral_set_property;
  object_class->get_property = ncm_spectral_get_property;

  /**
   * NcmSpectral:max-order:
   *
   * Maximum refinement level k for adaptive computations. The maximum
   * number of nodes is N_max = 2^max_order + 1.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_ORDER,
                                   g_param_spec_uint ("max-order",
                                                      NULL,
                                                      "Maximum refinement order",
                                                      1, G_MAXUINT, 16,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_spectral_new:
 *
 * Creates a new #NcmSpectral object with default maximum order (10).
 *
 * Returns: (transfer full): a new #NcmSpectral
 */
NcmSpectral *
ncm_spectral_new (void)
{
  return g_object_new (NCM_TYPE_SPECTRAL, NULL);
}

/**
 * ncm_spectral_new_with_max_order:
 * @max_order: maximum refinement level k (N_max = 2^@max_order + 1)
 *
 * Creates a new #NcmSpectral object with specified maximum order.
 *
 * Returns: (transfer full): a new #NcmSpectral
 */
NcmSpectral *
ncm_spectral_new_with_max_order (guint max_order)
{
  return g_object_new (NCM_TYPE_SPECTRAL,
                       "max-order", max_order,
                       NULL);
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
 * ncm_spectral_set_max_order:
 * @spectral: a #NcmSpectral
 * @max_order: maximum refinement level k (N_max = 2^@max_order + 1)
 *
 * Sets the maximum refinement order for adaptive computations.
 */
void
ncm_spectral_set_max_order (NcmSpectral *spectral, guint max_order)
{
  g_return_if_fail (NCM_IS_SPECTRAL (spectral));

  if (spectral->max_order != max_order)
  {
    spectral->max_order = max_order;

    /* Clear cached plans and arrays as they depend on max_order */
    if (spectral->fftw_plans != NULL)
    {
      g_ptr_array_unref (spectral->fftw_plans);
      spectral->fftw_plans = NULL;
    }

    if (spectral->cos_arrays != NULL)
    {
      g_ptr_array_unref (spectral->cos_arrays);
      spectral->cos_arrays = NULL;
    }

    if (spectral->f_vals != NULL)
    {
      fftw_free (spectral->f_vals);
      spectral->f_vals = NULL;
    }

    if (spectral->coeffs_work != NULL)
    {
      fftw_free (spectral->coeffs_work);
      spectral->coeffs_work = NULL;
    }
  }
}

/**
 * ncm_spectral_get_max_order:
 * @spectral: a #NcmSpectral
 *
 * Gets the maximum refinement order for adaptive computations.
 *
 * Returns: the maximum refinement order k
 */
guint
ncm_spectral_get_max_order (NcmSpectral *spectral)
{
  g_return_val_if_fail (NCM_IS_SPECTRAL (spectral), 0);

  return spectral->max_order;
}

/**
 * ncm_spectral_compute_chebyshev_coeffs:
 * @spectral: a #NcmSpectral
 * @F: (scope call): function to evaluate, receives x in [a,b]
 * @a: left endpoint of the interval
 * @b: right endpoint of the interval
 * @order: number of Chebyshev coefficients to compute
 * @coeffs: (out callee-allocates) (transfer full) (element-type gdouble): output array of coefficients
 * @user_data: user data for @F
 *
 * Computes Chebyshev coefficients of f(x) on [a,b] using FFTW DCT-I. The function @F
 * is sampled at Chebyshev nodes $x_k = (a+b)/2 - (b-a)/2\cos(k\pi/(N-1))$ which correspond
 * to the Chebyshev points $t_k = \cos(k\pi/(N-1))$ in $[-1,1]$ transformed to $[a,b]$.
 * The Chebyshev expansion is $f(x) = f(t) = \sum_{k=0}^{N-1} a_k T_k(t)$ where
 * $t = (2x - (a+b))/(b-a)$.
 *
 * If @coeffs points to NULL, allocates a new GArray of size @order. If @coeffs points
 * to an existing GArray, resizes it to @order. Through bindings, @coeffs always receives NULL.
 */
void
ncm_spectral_compute_chebyshev_coeffs (NcmSpectral *spectral, NcmSpectralF F, gdouble a, gdouble b, guint order, GArray **coeffs, gpointer user_data)
{
  const gdouble mid    = 0.5 * (a + b);
  const gdouble half_h = 0.5 * (b - a);
  const guint N        = order;
  guint i;

  if (*coeffs == NULL)
    *coeffs = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), N);

  g_array_set_size (*coeffs, N);

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

    if (spectral->cheb_cos_vals != NULL)
    {
      g_free (spectral->cheb_cos_vals);
      spectral->cheb_cos_vals = NULL;
    }

    /* Allocate new resources */
    spectral->cheb_f_vals   = fftw_malloc (sizeof (gdouble) * N);
    spectral->cheb_cos_vals = g_new (gdouble, N);

    /* Precompute cosine values at Chebyshev nodes */
    {
      const gdouble inv_Nm1 = 1.0 / (N - 1.0);
      const gdouble pi_Nm1  = M_PI * inv_Nm1;

      for (i = 0; i < N; i++)
        spectral->cheb_cos_vals[i] = cos (pi_Nm1 * i);
    }

    /* Create new FFTW plan */
    ncm_cfg_load_fftw_wisdom ("ncm_spectral");
    ncm_cfg_lock_plan_fftw ();
    spectral->cheb_plan_r2r = fftw_plan_r2r_1d (N, spectral->cheb_f_vals, (gdouble *) (*coeffs)->data,
                                                FFTW_REDFT00, ncm_cfg_get_fftw_default_flag ());
    ncm_cfg_unlock_plan_fftw ();
    ncm_cfg_save_fftw_wisdom ("ncm_spectral");

    spectral->cheb_N_cached = N;
  }

  /* Sample function at Chebyshev nodes using precomputed cosines */
  {
    gdouble * restrict f_vals       = spectral->cheb_f_vals;
    const gdouble * restrict c_vals = spectral->cheb_cos_vals;

    for (i = 0; i < N; i++)
    {
      const gdouble x = mid + half_h * c_vals[i];

      f_vals[i] = F (user_data, x);
    }
  }

  /* Execute FFTW plan */
  fftw_execute_r2r (spectral->cheb_plan_r2r, spectral->cheb_f_vals, (gdouble *) (*coeffs)->data);

  /* Normalize coefficients */
  {
    gdouble * restrict coeffs_data = (gdouble *) (*coeffs)->data;
    const gdouble inv_2Nm1         = 1.0 / (2.0 * (N - 1.0));
    const gdouble inv_Nm1          = 2.0 * inv_2Nm1;

    coeffs_data[0]     *= inv_2Nm1;
    coeffs_data[N - 1] *= inv_2Nm1;

    for (i = 1; i < N - 1; i++)
      coeffs_data[i] *= inv_Nm1;
  }
}

/* Helper functions for adaptive refinement */

static void
_ncm_spectral_prepare_plan_for_k (NcmSpectral *spectral, guint k)
{
  const guint N = (1 << k) + 1;
  guint j;

  /* Check if already prepared */
  if ((spectral->fftw_plans != NULL) && (k < spectral->fftw_plans->len) &&
      (g_ptr_array_index (spectral->fftw_plans, k) != NULL))
    return;

  /* Initialize arrays if needed */
  if (spectral->fftw_plans == NULL)
  {
    spectral->fftw_plans = g_ptr_array_new_with_free_func ((GDestroyNotify) fftw_destroy_plan);
    spectral->cos_arrays = g_ptr_array_new_with_free_func (g_free);
  }

  /* Ensure arrays are large enough */
  while (spectral->fftw_plans->len <= k)
  {
    g_ptr_array_add (spectral->fftw_plans, NULL);
    g_ptr_array_add (spectral->cos_arrays, NULL);
  }

  /* Allocate working arrays if needed */
  if (spectral->f_vals == NULL)
  {
    const guint N_max = (1 << spectral->max_order) + 1;

    spectral->f_vals      = fftw_malloc (sizeof (gdouble) * N_max);
    spectral->coeffs_work = fftw_malloc (sizeof (gdouble) * N_max);
  }

  /* Precompute Chebyshev-Lobatto cosines: cos(j*pi/2^k) */
  gdouble *cos_vals        = g_new (gdouble, N);
  const gdouble pi_over_2k = M_PI / (1 << k);

  for (j = 0; j < N; j++)
    cos_vals[j] = cos (j * pi_over_2k);

  g_ptr_array_index (spectral->cos_arrays, k) = cos_vals;

  /* Create out-of-place FFTW plan */
  ncm_cfg_load_fftw_wisdom ("ncm_spectral");
  ncm_cfg_lock_plan_fftw ();

  fftw_plan plan = fftw_plan_r2r_1d (N,
                                     spectral->f_vals,
                                     spectral->coeffs_work,
                                     FFTW_REDFT00,
                                     ncm_cfg_get_fftw_default_flag ());

  ncm_cfg_unlock_plan_fftw ();
  ncm_cfg_save_fftw_wisdom ("ncm_spectral");

  g_ptr_array_index (spectral->fftw_plans, k) = plan;
}

static void
_ncm_spectral_evaluate_all_nodes (NcmSpectral *spectral, NcmSpectralF F,
                                  gdouble a, gdouble b, guint k, gpointer user_data)
{
  const guint N           = (1 << k) + 1;
  const gdouble mid       = 0.5 * (a + b);
  const gdouble half_h    = 0.5 * (b - a);
  const gdouble *cos_vals = g_ptr_array_index (spectral->cos_arrays, k);
  guint j;

  for (j = 0; j < N; j++)
  {
    const gdouble x = mid + half_h * cos_vals[j];

    spectral->f_vals[j] = F (user_data, x);
  }
}

static void
_ncm_spectral_refine_to_k (NcmSpectral *spectral, NcmSpectralF F,
                           gdouble a, gdouble b, guint k_old, guint k_new,
                           gpointer user_data)
{
  const guint N_old       = (1 << k_old) + 1;
  const guint N_new       = (1 << k_new) + 1;
  const gdouble mid       = 0.5 * (a + b);
  const gdouble half_h    = 0.5 * (b - a);
  const gdouble *cos_vals = g_ptr_array_index (spectral->cos_arrays, k_new);
  gint j;
  guint jj;

  g_assert (k_new == k_old + 1);

  /* Move existing values to even positions (BACKWARD to avoid overwriting) */
  for (j = (gint) N_old - 1; j >= 0; j--)
    spectral->f_vals[2 * j] = spectral->f_vals[j];

  /* Compute new odd positions */
  for (jj = 1; jj < N_new; jj += 2)
  {
    const gdouble x = mid + half_h * cos_vals[jj];

    spectral->f_vals[jj] = F (user_data, x);
  }
}

static void
_ncm_spectral_normalize_coeffs (gdouble *coeffs_work, GArray *coeffs, guint N)
{
  const gdouble inv_2Nm1 = 1.0 / (2.0 * (N - 1.0));
  const gdouble inv_Nm1  = 2.0 * inv_2Nm1;
  gdouble *coeffs_data   = (gdouble *) coeffs->data;
  guint i;

  coeffs_data[0]     = coeffs_work[0] * inv_2Nm1;
  coeffs_data[N - 1] = coeffs_work[N - 1] * inv_2Nm1;

  for (i = 1; i < N - 1; i++)
  {
    coeffs_data[i] = coeffs_work[i] * inv_Nm1;
  }
}

static gboolean
_ncm_spectral_check_convergence (GArray *coeffs_2N, GArray *coeffs_N, gdouble tol)
{
  const gdouble *coeffs_2N_data = (gdouble *) coeffs_2N->data;
  const gdouble *coeffs_N_data  = (gdouble *) coeffs_N->data;
  gdouble norm2_diff            = 0.0;
  gdouble norm2_2N              = 0.0;
  guint i;

  for (i = 0; i < coeffs_N->len; i++)
  {
    const gdouble diff  = (coeffs_2N_data[i] - coeffs_N_data[i]);
    const gdouble diff2 = diff * diff;

    norm2_diff += diff2;
    norm2_2N   += coeffs_2N_data[i] * coeffs_2N_data[i];
  }

  if (norm2_diff < tol * tol * norm2_2N + 1.0e-100)
    return TRUE;

  return FALSE;
}

/**
 * ncm_spectral_compute_chebyshev_coeffs_adaptive:
 * @spectral: a #NcmSpectral
 * @F: (scope call): function to evaluate, receives x in [a,b]
 * @a: left endpoint of the interval
 * @b: right endpoint of the interval
 * @k_min: minimum refinement level (N_min = 2^@k_min + 1)
 * @tol: spectral convergence tolerance
 * @coeffs: (out callee-allocates) (transfer full) (element-type gdouble): output array of coefficients
 * @user_data: user data for @F
 *
 * Computes Chebyshev coefficients adaptively using nested Chebyshev-Lobatto nodes.
 * The function @F is evaluated at points x in [a,b]. Starts at level @k_min and refines
 * by doubling until spectral convergence is achieved or max_order is reached. Uses nested
 * nodes: only new odd nodes are computed at each refinement level.
 * The Chebyshev expansion is $f(x) = f(t) = \sum_{k=0}^{N-1} a_k T_k(t)$ where
 * $t = (2x - (a+b))/(b-a)$.
 *
 * If @coeffs points to NULL, allocates a new GArray. If @coeffs points to an existing
 * GArray, resizes it as needed. Through bindings, @coeffs always receives NULL.
 */
void
ncm_spectral_compute_chebyshev_coeffs_adaptive (NcmSpectral *spectral, NcmSpectralF F,
                                                gdouble a, gdouble b, guint k_min,
                                                gdouble tol, GArray **coeffs, gpointer user_data)
{
  guint k = k_min;
  GArray *c_tmp, *c_final;

  g_assert (k_min <= spectral->max_order);

  if (*coeffs == NULL)
    *coeffs = g_array_new (FALSE, FALSE, sizeof (gdouble));

  /* Initial evaluation at k_min */
  _ncm_spectral_prepare_plan_for_k (spectral, k);
  _ncm_spectral_evaluate_all_nodes (spectral, F, a, b, k, user_data);

  c_tmp   = spectral->coeffs;
  c_final = *coeffs;

  /* Transform using N and store in coeffs_work */
  {
    const guint N  = (1 << k) + 1;
    fftw_plan plan = g_ptr_array_index (spectral->fftw_plans, k);

    /* Transform f_vals -> coeffs_work */
    fftw_execute (plan);

    g_array_set_size (c_tmp, N);
    _ncm_spectral_normalize_coeffs (spectral->coeffs_work, c_tmp, N);
  }

  while (k < spectral->max_order)
  {
    /* Transform using 2N and store in coeffs */
    _ncm_spectral_prepare_plan_for_k (spectral, k + 1);
    _ncm_spectral_refine_to_k (spectral, F, a, b, k, k + 1, user_data);
    k++;
    {
      const guint N  = (1 << k) + 1;
      fftw_plan plan = g_ptr_array_index (spectral->fftw_plans, k);

      /* Transform f_vals -> coeffs_work */
      fftw_execute (plan);

      g_array_set_size (c_final, N);
      _ncm_spectral_normalize_coeffs (spectral->coeffs_work, c_final, N);
    }

    /* Check convergence using pre-computed e_total */
    if (_ncm_spectral_check_convergence (c_final, c_tmp, tol))
      break;

    /* Swap c_tmp and c_final */
    {
      GArray *tmp = c_tmp;

      c_tmp   = c_final;
      c_final = tmp;
    }
  }

  if (c_final != *coeffs)
  {
    /* If final coefficients are not in *coeffs, copy them */
    g_array_set_size (*coeffs, c_final->len);
    memcpy ((*coeffs)->data, c_final->data, sizeof (gdouble) * c_final->len);
  }
}

/**
 * ncm_spectral_chebT_to_gegenbauer_alpha1:
 * @c: (element-type gdouble): Chebyshev coefficients array
 * @g: (out callee-allocates) (transfer full) (element-type gdouble): Gegenbauer $C^{(1)}_n$ coefficients array
 *
 * Converts Chebyshev $T_n$ coefficients to Gegenbauer $C^{(1)}_n$ coefficients ($\alpha=1$).
 * Uses the relationship: $T_n = \frac{1}{2}(C^{(1)}_n + C^{(1)}_{n-2})$ for $n \geq 2$.
 *
 * If @g points to NULL, allocates a new GArray with same size as @c. If @g points to an
 * existing GArray, resizes it to match @c. Through bindings, @g always receives NULL.
 */
void
ncm_spectral_chebT_to_gegenbauer_alpha1 (GArray *c, GArray **g)
{
  const guint N = c->len;
  guint i;

  if (*g == NULL)
    *g = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), N);

  g_array_set_size (*g, N);

  if (N == 0)
    return;

  {
    const gdouble *c_data = (gdouble *) c->data;
    gdouble *g_data       = (gdouble *) (*g)->data;

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
 * @c: (element-type gdouble): Chebyshev coefficients array
 * @g: (out callee-allocates) (transfer full) (element-type gdouble): Gegenbauer $C^{(2)}_k$ coefficients array
 *
 * Converts Chebyshev $T_n$ coefficients to Gegenbauer $C^{(2)}_k$ coefficients ($\alpha=2$).
 *
 * Uses the projection formula:
 * $$g_k = \frac{1}{2} c_0 \delta_{k,0} + \frac{c_k}{2(k+1)} - \frac{(k+2) c_{k+2}}{(k+1)(k+3)} + \frac{c_{k+4}}{2(k+3)}$$
 * where $f(x) = \sum_n c_n T_n(x)$.
 *
 * If @g points to NULL, allocates a new GArray with same size as @c. If @g points to an
 * existing GArray, resizes it to match @c. Through bindings, @g always receives NULL.
 */
void
ncm_spectral_chebT_to_gegenbauer_alpha2 (GArray *c, GArray **g)
{
  const guint N = c->len;
  guint k;

  if (*g == NULL)
    *g = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), N);

  g_array_set_size (*g, N);

  if (N == 0)
    return;

  {
    const gdouble *c_data = (gdouble *) c->data;
    gdouble *g_data       = (gdouble *) (*g)->data;

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
 * @c: (element-type gdouble): Gegenbauer $C^{(1)}_n$ coefficients array
 * @t: point to evaluate in [-1, 1]
 *
 * Evaluates a Gegenbauer $C^{(1)}_n$ expansion at t using Clenshaw recurrence.
 * For $\alpha=1$, $C^{(1)}_n(t) = U_n(t)$ (Chebyshev polynomials of the second kind).
 * The variable t should be in the interval [-1, 1]. To evaluate at a point x in [a, b],
 * use ncm_spectral_gegenbauer_alpha1_eval_x() or first convert x to t using
 * ncm_spectral_x_to_t().
 *
 * Returns: the value of $\sum_{n=0}^{N-1} c_n C^{(1)}_n(t)$
 */
gdouble
ncm_spectral_gegenbauer_alpha1_eval (GArray *c, gdouble t)
{
  const guint N = c->len;

  if (N == 0)
    return 0.0;

  {
    const gdouble *c_data = (gdouble *) c->data;

    /* Endpoint handling: C_n^{(1)}(+/-1) = (n+1)*(+/-1)^n */
    if (fabs (t - 1.0) < 1e-15)
    {
      gdouble sum = 0.0;
      guint n;

      for (n = 0; n < N; n++)
        sum += c_data[n] * (gdouble) (n + 1);

      return sum;
    }

    if (fabs (t + 1.0) < 1e-15)
    {
      gdouble sum = 0.0;
      guint n;

      for (n = 0; n < N; n++)
        sum += c_data[n] * ((n & 1) ? -(gdouble) (n + 1) : (gdouble) (n + 1));

      return sum;
    }

    {
      /* Stable recurrence for interior t */
      gdouble Cnm1 = 1.0; /* U_0 */
      gdouble sum  = c_data[0] * Cnm1;

      if (N == 1)
        return sum;

      gdouble Cn = 2.0 * t; /* U_1 */

      sum += c_data[1] * Cn;

      for (guint n = 1; n < N - 1; n++)
      {
        gdouble Cnp1 = 2.0 * t * Cn - Cnm1; /* U_{n+1} */

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
 * @c: (element-type gdouble): Gegenbauer $C^{(2)}_n$ coefficients array
 * @t: point to evaluate in [-1, 1]
 *
 * Evaluates a Gegenbauer $C^{(2)}_n$ expansion at t using Clenshaw recurrence.
 * For $\alpha=2$, the recurrence relation is:
 * $(n+1) C^{(2)}_{n+1}(t) = 2(n+2)t C^{(2)}_n(t) - (n+3) C^{(2)}_{n-1}(t)$
 * The variable t should be in the interval [-1, 1]. To evaluate at a point x in [a, b],
 * use ncm_spectral_gegenbauer_alpha2_eval_x() or first convert x to t using
 * ncm_spectral_x_to_t().
 *
 * Returns: the value of $\sum_{n=0}^{N-1} c_n C^{(2)}_n(t)$
 */
gdouble
ncm_spectral_gegenbauer_alpha2_eval (GArray *c, gdouble t)
{
  const guint N = c->len;

  if (N == 0)
    return 0.0;

  {
    const gdouble *c_data = (gdouble *) c->data;

    /* Endpoint handling: C_n^{(2)}(+/-1) = binom(n+3,3)*(+/-1)^n = ((n+1)*(n+2)*(n+3)/6)*(+/-1)^n */
    if (fabs (t - 1.0) < 1e-15)
    {
      gdouble sum = 0.0;
      guint n;

      for (n = 0; n < N; n++)
        sum += c_data[n] * (gdouble) ((n + 1) * (n + 2) * (n + 3)) / 6.0;

      return sum;
    }

    if (fabs (t + 1.0) < 1e-15)
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
      /* Stable recurrence for interior t */
      gdouble Cnm1 = 1.0; /* C_0^{(2)} = 1 */
      gdouble sum  = c_data[0] * Cnm1;

      if (N == 1)
        return sum;

      gdouble Cn = 4.0 * t; /* C_1^{(2)} = 4t */

      sum += c_data[1] * Cn;

      for (guint n = 1; n < N - 1; n++)
      {
        /* (n+1) C_{n+1}^{(2)} = 2(n+2)t C_n^{(2)} - (n+3) C_{n-1}^{(2)} */
        gdouble Cnp1 = (2.0 * (gdouble) (n + 2) * t * Cn - (gdouble) (n + 3) * Cnm1) / (gdouble) (n + 1);

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
 * @a: (element-type gdouble): Chebyshev coefficients array
 * @t: point to evaluate in [-1, 1]
 *
 * Evaluates a Chebyshev expansion $f(t) = \sum_{k=0}^{N-1} a_k T_k(t)$ at t
 * using Clenshaw recurrence with Reinsch modification near endpoints.
 * The variable t should be in the interval [-1, 1]. To evaluate at a point x in [a, b],
 * use ncm_spectral_chebyshev_eval_x() or first convert x to t using ncm_spectral_x_to_t().
 *
 * Returns: the value of the Chebyshev expansion at t
 */
gdouble
ncm_spectral_chebyshev_eval (GArray *a, gdouble t)
{
  const guint N           = a->len;
  const gdouble *a_data   = (gdouble *) a->data;
  const gdouble threshold = 0.9;
  const gdouble eps       = 1e-15;

  if (N == 0)
    return 0.0;

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
 * @a: (element-type gdouble): Chebyshev coefficients array (a_j multiplies T_j)
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
ncm_spectral_chebyshev_deriv (GArray *a, gdouble t)
{
  const gint N = a->len;

  if (N <= 1)
    return 0.0;

  {
    const gdouble *a_data = (gdouble *) a->data;

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

/**
 * ncm_spectral_gegenbauer_alpha1_eval_x:
 * @c: (element-type gdouble): Gegenbauer $C^{(1)}_n$ coefficients array
 * @a: left endpoint of the interval
 * @b: right endpoint of the interval
 * @x: point to evaluate in [a, b]
 *
 * Evaluates a Gegenbauer $C^{(1)}_n$ expansion at a point x in [a, b].
 * This function converts x to t using $t = (2x - (a+b))/(b-a)$ and then
 * calls ncm_spectral_gegenbauer_alpha1_eval().
 *
 * Returns: the value of $\sum_{n=0}^{N-1} c_n C^{(1)}_n(t)$ where $t = (2x - (a+b))/(b-a)$
 */
gdouble
ncm_spectral_gegenbauer_alpha1_eval_x (GArray *c, gdouble a, gdouble b, gdouble x)
{
  const gdouble t = ncm_spectral_x_to_t (a, b, x);

  return ncm_spectral_gegenbauer_alpha1_eval (c, t);
}

/**
 * ncm_spectral_gegenbauer_alpha2_eval_x:
 * @c: (element-type gdouble): Gegenbauer $C^{(2)}_n$ coefficients array
 * @a: left endpoint of the interval
 * @b: right endpoint of the interval
 * @x: point to evaluate in [a, b]
 *
 * Evaluates a Gegenbauer $C^{(2)}_n$ expansion at a point x in [a, b].
 * This function converts x to t using $t = (2x - (a+b))/(b-a)$ and then
 * calls ncm_spectral_gegenbauer_alpha2_eval().
 *
 * Returns: the value of $\sum_{n=0}^{N-1} c_n C^{(2)}_n(t)$ where $t = (2x - (a+b))/(b-a)$
 */
gdouble
ncm_spectral_gegenbauer_alpha2_eval_x (GArray *c, gdouble a, gdouble b, gdouble x)
{
  const gdouble t = ncm_spectral_x_to_t (a, b, x);

  return ncm_spectral_gegenbauer_alpha2_eval (c, t);
}

/**
 * ncm_spectral_chebyshev_eval_x:
 * @a: (element-type gdouble): Chebyshev coefficients array
 * @a_v: left endpoint of the interval
 * @b: right endpoint of the interval
 * @x: point to evaluate in [a_v, b]
 *
 * Evaluates a Chebyshev expansion at a point x in [a_v, b].
 * This function converts x to t using $t = (2x - (a_v+b))/(b-a_v)$ and then
 * calls ncm_spectral_chebyshev_eval().
 *
 * Returns: the value of $\sum_{k=0}^{N-1} a_k T_k(t)$ where $t = (2x - (a_v+b))/(b-a_v)$
 */
gdouble
ncm_spectral_chebyshev_eval_x (GArray *a, gdouble a_v, gdouble b, gdouble x)
{
  const gdouble t = ncm_spectral_x_to_t (a_v, b, x);

  return ncm_spectral_chebyshev_eval (a, t);
}

/**
 * ncm_spectral_chebyshev_deriv_x:
 * @a: (element-type gdouble): Chebyshev coefficients array
 * @a_v: left endpoint of the interval
 * @b: right endpoint of the interval
 * @x: point to evaluate in [a_v, b]
 *
 * Evaluates the first derivative of a Chebyshev expansion at a point x in [a_v, b].
 * This function converts x to t using $t = (2x - (a_v+b))/(b-a_v)$, evaluates
 * the derivative with respect to t using ncm_spectral_chebyshev_deriv(), and then
 * applies the chain rule: $df/dx = (df/dt) \cdot (dt/dx) = (df/dt) \cdot 2/(b-a_v)$.
 *
 * Returns: the value of $df/dx$ at x
 */
gdouble
ncm_spectral_chebyshev_deriv_x (GArray *a, gdouble a_v, gdouble b, gdouble x)
{
  const gdouble t     = ncm_spectral_x_to_t (a_v, b, x);
  const gdouble df_dt = ncm_spectral_chebyshev_deriv (a, t);

  /* Chain rule: df/dx = (df/dt) * (dt/dx) = (df/dt) * 2/(b-a) */
  return df_dt * 2.0 / (b - a_v);
}

/**
 * ncm_spectral_get_proj_matrix:
 * @N: size of the matrix
 *
 * Returns the projection (identity) operator matrix that transforms Chebyshev $T_n$
 * coefficients to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Returns: (transfer full): the projection operator matrix
 */
NcmMatrix *
ncm_spectral_get_proj_matrix (guint N)
{
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  const glong bandwidth       = 9;
  gdouble * restrict row_data = g_new0 (gdouble, bandwidth);
  glong j, k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong cols_to_write = GSL_MIN (bandwidth, N - k);

    /* First entry is k, offset = 0 */
    memset (row_data, 0, sizeof (gdouble) * bandwidth);
    ncm_spectral_compute_proj_row (row_data, k, 0, 1.0);

    for (j = 0; j < cols_to_write; j++)
    {
      ncm_matrix_set (mat, k, k + j, row_data[j]);
    }
  }

  g_free (row_data);

  return mat;
}

/**
 * ncm_spectral_get_x_matrix:
 * @N: size of the matrix
 *
 * Returns the multiplication by $x$ operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $x \cdot f(x)$.
 *
 * Returns: (transfer full): the $x$ operator matrix
 */
NcmMatrix *
ncm_spectral_get_x_matrix (guint N)
{
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  const glong bandwidth       = 9;
  gdouble * restrict row_data = g_new0 (gdouble, bandwidth);
  glong j, k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong offset        = (k >= 1) ? 1 : 0;
    const glong cols_to_write = GSL_MIN (bandwidth, (glong) N + offset - k);

    memset (row_data, 0, sizeof (gdouble) * bandwidth);
    ncm_spectral_compute_x_row (row_data, k, offset, 1.0);

    for (j = 0; j < cols_to_write; j++)
    {
      ncm_matrix_set (mat, k, k - offset + j, row_data[j]);
    }
  }

  g_free (row_data);

  return mat;
}

/**
 * ncm_spectral_get_x2_matrix:
 * @N: size of the matrix
 *
 * Returns the multiplication by $x^2$ operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $x^2 \cdot f(x)$.
 *
 * Returns: (transfer full): the $x^2$ operator matrix
 */
NcmMatrix *
ncm_spectral_get_x2_matrix (guint N)
{
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  const glong bandwidth       = 9;
  gdouble * restrict row_data = g_new0 (gdouble, bandwidth);
  glong j, k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong offset        = (k >= 2) ? 2 : k;
    const glong cols_to_write = GSL_MIN (bandwidth, (glong) N + offset - k);

    memset (row_data, 0, sizeof (gdouble) * bandwidth);
    ncm_spectral_compute_x2_row (row_data, k, offset, 1.0);

    for (j = 0; j < cols_to_write; j++)
    {
      const glong col = k - offset + j;

      ncm_matrix_set (mat, k, col, row_data[j]);
    }
  }

  g_free (row_data);

  return mat;
}

/**
 * ncm_spectral_get_d_matrix:
 * @N: size of the matrix
 *
 * Returns the derivative operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $\frac{df}{dx}$.
 *
 * Returns: (transfer full): the derivative operator matrix
 */
NcmMatrix *
ncm_spectral_get_d_matrix (guint N)
{
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  const glong bandwidth       = 9;
  gdouble * restrict row_data = g_new0 (gdouble, bandwidth);
  glong j, k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong offset        = (k >= 1) ? 1 : k;
    const glong cols_to_write = GSL_MIN (bandwidth, (glong) N + offset - k);

    memset (row_data, 0, sizeof (gdouble) * bandwidth);
    ncm_spectral_compute_d_row (row_data, offset, 1.0);

    for (j = 0; j < cols_to_write; j++)
    {
      const glong col = k - offset + j;

      ncm_matrix_set (mat, k, col, row_data[j]);
    }
  }

  g_free (row_data);

  return mat;
}

/**
 * ncm_spectral_get_x_d_matrix:
 * @N: size of the matrix
 *
 * Returns the $x \cdot \frac{d}{dx}$ operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $x \cdot \frac{df}{dx}$.
 *
 * Returns: (transfer full): the $x \cdot d$ operator matrix
 */
NcmMatrix *
ncm_spectral_get_x_d_matrix (guint N)
{
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  const glong bandwidth       = 9;
  gdouble * restrict row_data = g_new0 (gdouble, bandwidth);
  glong j, k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong offset        = (k >= 1) ? 1 : k;
    const glong cols_to_write = GSL_MIN (bandwidth, (glong) N + offset - k);

    memset (row_data, 0, sizeof (gdouble) * bandwidth);
    ncm_spectral_compute_x_d_row (row_data, k, offset, 1.0);

    for (j = 0; j < cols_to_write; j++)
    {
      const glong col = k - offset + j;

      ncm_matrix_set (mat, k, col, row_data[j]);
    }
  }

  g_free (row_data);

  return mat;
}

/**
 * ncm_spectral_get_d2_matrix:
 * @N: size of the matrix
 *
 * Returns the second derivative operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $\frac{d^2f}{dx^2}$.
 *
 * Returns: (transfer full): the second derivative operator matrix
 */
NcmMatrix *
ncm_spectral_get_d2_matrix (guint N)
{
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  const glong bandwidth       = 9;
  gdouble * restrict row_data = g_new0 (gdouble, bandwidth);
  glong j, k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong offset        = (k >= 2) ? 2 : k;
    const glong cols_to_write = GSL_MIN (bandwidth, (glong) N + offset - k);

    memset (row_data, 0, sizeof (gdouble) * bandwidth);
    ncm_spectral_compute_d2_row (row_data, k, offset, 1.0);

    for (j = 0; j < cols_to_write; j++)
    {
      const glong col = k - offset + j;

      ncm_matrix_set (mat, k, col, row_data[j]);
    }
  }

  g_free (row_data);

  return mat;
}

/**
 * ncm_spectral_get_x_d2_matrix:
 * @N: size of the matrix
 *
 * Returns the $x \cdot \frac{d^2}{dx^2}$ operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $x \cdot \frac{d^2f}{dx^2}$.
 *
 * Returns: (transfer full): the $x \cdot d^2$ operator matrix
 */
NcmMatrix *
ncm_spectral_get_x_d2_matrix (guint N)
{
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  const glong bandwidth       = 9;
  gdouble * restrict row_data = g_new0 (gdouble, bandwidth);
  glong j, k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong offset        = (k >= 1) ? 1 : k;
    const glong cols_to_write = GSL_MIN (bandwidth, (glong) N + offset - k);

    memset (row_data, 0, sizeof (gdouble) * bandwidth);
    ncm_spectral_compute_x_d2_row (row_data, k, offset, 1.0);

    for (j = 0; j < cols_to_write; j++)
    {
      const glong col = k - offset + j;

      ncm_matrix_set (mat, k, col, row_data[j]);
    }
  }

  g_free (row_data);

  return mat;
}

/**
 * ncm_spectral_get_x2_d2_matrix:
 * @N: size of the matrix
 *
 * Returns the $x^2 \cdot \frac{d^2}{dx^2}$ operator matrix that transforms Chebyshev $T_n$
 * coefficients of $f(x)$ to Gegenbauer $C^{(2)}_k$ coefficients of $x^2 \cdot \frac{d^2f}{dx^2}$.
 *
 * Returns: (transfer full): the $x^2 \cdot d^2$ operator matrix
 */
NcmMatrix *
ncm_spectral_get_x2_d2_matrix (guint N)
{
  NcmMatrix *mat              = ncm_matrix_new (N, N);
  const glong bandwidth       = 9;
  gdouble * restrict row_data = g_new0 (gdouble, bandwidth);
  glong j, k;

  ncm_matrix_set_zero (mat);

  for (k = 0; k < N; k++)
  {
    const glong offset        = (k >= 2) ? 2 : k;
    const glong cols_to_write = GSL_MIN (bandwidth, (glong) N + offset - k);

    memset (row_data, 0, sizeof (gdouble) * bandwidth);
    ncm_spectral_compute_x2_d2_row (row_data, k, offset, 1.0);

    for (j = 0; j < cols_to_write; j++)
    {
      const glong col = k - offset + j;

      ncm_matrix_set (mat, k, col, row_data[j]);
    }
  }

  g_free (row_data);

  return mat;
}

/**
 * ncm_spectral_compute_proj_row:
 * @row_data: row structure to update
 * @k: row index (0-based)
 * @offset: row offset
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the projection operator that maps Chebyshev coefficients
 * to Gegenbauer $C^{(2)}_k$ coefficients (ultraspherical basis with $\lambda=2$).
 * This is the identity operator expressed in different bases.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: Gegenbauer $C^{(2)}_k(x)$ basis coefficient (row k)
 *
 * Mathematical formula:
 * $$
 * C^{(2)}_k = \frac{1}{2} c_0 \delta_{k,0} +
 * \frac{c_k}{2(k+1)} - \frac{(k+2) c_{k+2}}{(k+1)(k+3)} + \frac{c_{k+4}}{2(k+3)}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k: coeff * 1/(2*(k+1))
 * - column k+2: coeff * -(k+2)/((k+1)*(k+3))
 * - column k+4: coeff * 1/(2*(k+3))
 * - For k=0 only: additional value coeff * 1/2 at column 0
 */

/**
 * ncm_spectral_compute_x_row:
 * @row_data: row structure to update
 * @k: row index (0-based)
 * @offset: row offset
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the multiplication by x operator that maps Chebyshev
 * coefficients to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: Gegenbauer $C^{(2)}_k(x)$ basis coefficient for $x \cdot f(x)$ (row k)
 *
 * Mathematical formula:
 * $$
 * (x \cdot f)^{(2)}_k = \frac{\theta(k-1) c_{k-1}}{4(k+1)} - \frac{c_{k+1}}{4(k+3)} -
 * \frac{c_{k+3}}{4(k+1)} + \frac{c_{k+5}}{4(k+3)}
 * $$
 * plus special contributions: $\frac{c_1}{4}\delta_{k,0}$ and $\frac{c_0}{8}\delta_{k,1}$,
 * where $c_n$ are the input Chebyshev coefficients and $\theta$ is the Heaviside function.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - If k >= 1: column k-1: coeff * 1/(4*(k+1))
 * - column k+1: coeff * -1/(4*(k+3))
 * - column k+3: coeff * -1/(4*(k+1))
 * - column k+5: coeff * 1/(4*(k+3))
 * - For k=0 only: additional value coeff * 1/4 at column 1
 * - For k=1 only: additional value coeff * 1/8 at column 0
 */

/**
 * ncm_spectral_compute_x2_row:
 * @row_data: row structure to update
 * @k: row index (0-based)
 * @offset: row offset
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the multiplication by $x^2$ operator that maps Chebyshev
 * coefficients to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: Gegenbauer $C^{(2)}_k(x)$ basis coefficient for $x^2 \cdot f(x)$ (row k)
 *
 * Mathematical formula:
 * $$
 * (x^2 \cdot f)^{(2)}_k = \frac{\theta(k-2) c_{k-2}}{8(k+1)} + \frac{c_k}{4(k+1)(k+3)} -
 * \frac{(k+2) c_{k+2}}{4(k+1)(k+3)} - \frac{c_{k+4}}{4(k+1)(k+3)} + \frac{c_{k+6}}{8(k+3)}
 * $$
 * plus special contributions at k=0,1,2, where $c_n$ are the input Chebyshev
 * coefficients and $\theta$ is the Heaviside function.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - If k >= 2: column k-2: coeff * 1/(8*(k+1))
 * - column k: coeff * 1/(4*(k+1)*(k+3))
 * - column k+2: coeff * -(k+2)/(4*(k+1)*(k+3))
 * - column k+4: coeff * -1/(4*(k+1)*(k+3))
 * - column k+6: coeff * 1/(8*(k+3))
 * - For k=0 only: additional values coeff * 1/12 at column 0 and coeff * 1/8 at column 2
 * - For k=1 only: additional value coeff * 1/16 at column 1
 * - For k=2 only: additional value coeff * 1/24 at column 0
 */

/**
 * ncm_spectral_compute_d_row:
 * @row_data: row structure to update
 * @offset: row offset
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the first derivative operator $\frac{d}{dx}$ that maps Chebyshev
 * coefficients to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: $\langle C^{(2)}_k, f' \rangle$ - projection of $f'$ (row k)
 *
 * Mathematical formula:
 * $$
 * \langle C^{(2)}_k, f' \rangle = c_{k+1} - c_{k+3}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k+1: coeff * 1.0
 * - column k+3: coeff * (-1.0)
 */

/**
 * ncm_spectral_compute_x_d_row:
 * @row_data: row structure to update
 * @k: row index (0-based)
 * @offset: row offset
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the $x \cdot \frac{d}{dx}$ operator that maps Chebyshev coefficients
 * to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: $\langle C^{(2)}_k, x \cdot f' \rangle$ - projection of $x \cdot f'$ (row k)
 *
 * Mathematical formula:
 * $$
 * \langle C^{(2)}_k, x \cdot f' \rangle = \frac{k c_k}{2(k+1)} + \frac{(k+2) c_{k+2}}{(k+1)(k+3)} -
 * \frac{(k+4) c_{k+4}}{2(k+3)}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k: coeff * k/(2*(k+1))
 * - column k+2: coeff * (k+2)/((k+1)*(k+3))
 * - column k+4: coeff * -(k+4)/(2*(k+3))
 */

/**
 * ncm_spectral_compute_d2_row:
 * @row_data: row structure to update
 * @k: row index (0-based)
 * @offset: row offset
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the second derivative operator $\frac{d^2}{dx^2}$ that maps Chebyshev
 * coefficients to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: $\langle C^{(2)}_k, f'' \rangle$ - projection of $f''$ (row k)
 *
 * Mathematical formula:
 * $$
 * \langle C^{(2)}_k, f'' \rangle = 2(k+2) c_{k+2}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k+2: coeff * 2*(k+2) (single non-zero entry)
 */

/**
 * ncm_spectral_compute_x_d2_row:
 * @row_data: row structure to update
 * @k: row index (0-based)
 * @offset: row offset
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the $x \cdot \frac{d^2}{dx^2}$ operator that maps Chebyshev coefficients
 * to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns)
 * Output: $\langle C^{(2)}_k, x \cdot f'' \rangle$ - projection of $x \cdot f''$ (row k)
 *
 * Mathematical formula:
 * $$
 * \langle C^{(2)}_k, x \cdot f'' \rangle = k c_{k+1} + (k+4) c_{k+3}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k+1: coeff * k
 * - column k+3: coeff * (k+4)
 */

/**
 * ncm_spectral_compute_x2_d2_row:
 * @row_data: row structure to update
 * @k: row index (0-based)
 * @offset: row offset
 * @coeff: coefficient to multiply all elements
 *
 * Computes row k of the $x^2 \cdot \frac{d^2}{dx^2}$ operator that maps Chebyshev
 * coefficients to Gegenbauer $C^{(2)}_k$ coefficients.
 *
 * Input: Chebyshev $T_n(x)$ basis coefficients (columns) Output: $\langle C^{(2)}_k,
 * x^2 \cdot f'' \rangle$ - projection of $x^2 \cdot f''$ (row k)
 *
 * Mathematical formula:
 * $$
 * \langle C^{(2)}_k, x^2 \cdot f'' \rangle = \frac{k(k-1) c_k}{2(k+1)} +
 * \frac{(k+2)((k+2)^2-3) c_{k+2}}{(k+1)(k+3)} + \frac{(k+4)(k+5) c_{k+4}}{2(k+3)}
 * $$
 * where $c_n$ are the input Chebyshev coefficients.
 *
 * Adds to existing row data (for linear combinations of operators).
 *
 * Matrix entries for row k:
 * - column k: coeff * k*(k-1)/(2*(k+1))
 * - column k+2: coeff * (k+2)*((k+2)^2 - 3)/((k+1)*(k+3))
 * - column k+4: coeff * (k+4)*(k+5)/(2*(k+3))
 */

