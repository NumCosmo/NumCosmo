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
 * NcmSpectral:
 *
 * Chebyshev series of functions on an interval.
 *
 * Computes the Chebyshev coefficients of a function on $[a, b]$ from its values at the
 * Chebyshev-Lobatto nodes by a DCT-I, at a fixed order or adaptively on nested grids.
 * Evaluates, differentiates and integrates Chebyshev series, moves them to another
 * interval, and provides the banded ultraspherical (Gegenbauer) operators of the spectral
 * method for linear ODEs. The operator matrices and rows act on the Chebyshev variable,
 * written $x$ there, on $[-1, 1]$. See
 * <a href="../../theory/ncm/algebra/spectral.html">Spectral Methods</a> for the formulas.
 *
 * A coefficient array passed as a pointer to a pointer is reused when it is not %NULL
 * and allocated otherwise; through bindings a new one is always returned. An instance
 * holds the node and transform buffers of its expansions, so it must not be used by two
 * threads at once, and a batch expansion must not be called again from inside its own
 * callback.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_spectral.h"
#include "ncm/core/ncm_cfg.h"

#include <math.h>
#ifndef NUMCOSMO_GIR_SCAN
#include <fftw3.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* Largest max-order for which 2^max-order + 1 fits a guint */
#define NCM_SPECTRAL_MAX_ORDER_LIMIT (30)

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
  gdouble *f_vals_tmp;   /* Function values array, size: 2^max_order + 1 */
  gdouble *coeffs_work;  /* Coefficients work array, size: 2^max_order + 1 */
  GArray *coeffs;        /* Coefficients array, size: 2^max_order + 1 */
  GPtrArray *cos_arrays; /* Precomputed cosines for each k level */
  GPtrArray *fftw_plans; /* FFTW plans for each k level */

  /* Scratch for the interval rebase: three coefficient rows */
  gdouble *rebase_work;
  gsize rebase_work_len;

  /* Batch expansion of n_comp components, stored node-major, f_vals[j * n_comp + c] */
  guint batch_n_comp;
  gdouble *batch_f_vals;
  gdouble *batch_f_vals_tmp;
  gdouble *batch_coeffs_work;
  GPtrArray *batch_fftw_plans;

  /* Legacy fields (backward compatibility) */
  guint cheb_N_cached;     /* Cached N value */
  gdouble *cheb_f_vals;    /* Cached function values array */
  gdouble *cheb_c_vals;    /* Cached coefficient output array the plan is bound to */
  gdouble *cheb_cos_vals;  /* Cached cosine values at Chebyshev nodes */
  fftw_plan cheb_plan_r2r; /* Cached FFTW plan */
};

G_DEFINE_TYPE (NcmSpectral, ncm_spectral, G_TYPE_OBJECT)

static void
ncm_spectral_init (NcmSpectral *spectral)
{
  spectral->max_order   = 0;
  spectral->f_vals      = NULL;
  spectral->f_vals_tmp  = NULL;
  spectral->coeffs_work = NULL;
  spectral->coeffs      = g_array_new (FALSE, FALSE, sizeof (gdouble));
  spectral->cos_arrays  = NULL;
  spectral->fftw_plans  = NULL;

  spectral->batch_n_comp      = 0;
  spectral->batch_f_vals      = NULL;
  spectral->batch_f_vals_tmp  = NULL;
  spectral->batch_coeffs_work = NULL;
  spectral->batch_fftw_plans  = NULL;

  spectral->rebase_work     = NULL;
  spectral->rebase_work_len = 0;

  spectral->cheb_N_cached = 0;
  spectral->cheb_f_vals   = NULL;
  spectral->cheb_c_vals   = NULL;
  spectral->cheb_cos_vals = NULL;
  spectral->cheb_plan_r2r = NULL;
}

static void
ncm_spectral_finalize (GObject *object)
{
  NcmSpectral *spectral = NCM_SPECTRAL (object);

  /* Clean up adaptive refinement resources */
  g_clear_pointer (&spectral->fftw_plans, g_ptr_array_unref);
  g_clear_pointer (&spectral->cos_arrays, g_ptr_array_unref);
  g_clear_pointer (&spectral->f_vals, fftw_free);
  g_clear_pointer (&spectral->f_vals_tmp, fftw_free);
  g_clear_pointer (&spectral->coeffs_work, fftw_free);
  g_clear_pointer (&spectral->coeffs, g_array_unref);
  g_clear_pointer (&spectral->rebase_work, g_free);

  g_clear_pointer (&spectral->batch_fftw_plans, g_ptr_array_unref);
  g_clear_pointer (&spectral->batch_f_vals, fftw_free);
  g_clear_pointer (&spectral->batch_f_vals_tmp, fftw_free);
  g_clear_pointer (&spectral->batch_coeffs_work, fftw_free);

  g_clear_pointer (&spectral->cheb_plan_r2r, ncm_cfg_fftw_plan_destroy);
  g_clear_pointer (&spectral->cheb_f_vals, fftw_free);
  g_clear_pointer (&spectral->cheb_c_vals, fftw_free);
  g_clear_pointer (&spectral->cheb_cos_vals, g_free);

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
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
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
   * Highest refinement level $k$ of the adaptive expansions, at most $2^k + 1$ nodes, from 1
   * to 30. Buffers of that size are allocated when it is set. It bounds memory: the adaptive
   * expansions stop on their tolerance, and reaching this level without converging is an
   * error, except for the variants that report it.
   */
  g_object_class_install_property (object_class,
                                   PROP_MAX_ORDER,
                                   g_param_spec_uint ("max-order",
                                                      NULL,
                                                      "Maximum refinement order",
                                                      1, NCM_SPECTRAL_MAX_ORDER_LIMIT, 16,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_spectral_new:
 *
 * Creates a new #NcmSpectral with #NcmSpectral:max-order 16.
 *
 * Returns: (transfer full): a new #NcmSpectral.
 */
NcmSpectral *
ncm_spectral_new (void)
{
  return g_object_new (NCM_TYPE_SPECTRAL, NULL);
}

/**
 * ncm_spectral_new_with_max_order:
 * @max_order: the #NcmSpectral:max-order
 *
 * Returns: (transfer full): a new #NcmSpectral.
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
 * Returns: (transfer full): @spectral.
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
 * Decreases the reference count of @spectral by one.
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
 * If *@spectral is not %NULL, decreases its reference count by one and sets *@spectral
 * to %NULL.
 */
void
ncm_spectral_clear (NcmSpectral **spectral)
{
  g_clear_object (spectral);
}

/**
 * ncm_spectral_set_max_order:
 * @spectral: a #NcmSpectral
 * @max_order: the #NcmSpectral:max-order
 *
 * Sets #NcmSpectral:max-order, reallocating the buffers when it changes.
 */
void
ncm_spectral_set_max_order (NcmSpectral *spectral, guint max_order)
{
  g_return_if_fail (NCM_IS_SPECTRAL (spectral));
  g_assert_cmpuint (max_order, >=, 1);
  g_assert_cmpuint (max_order, <=, NCM_SPECTRAL_MAX_ORDER_LIMIT);

  if (spectral->max_order != max_order)
  {
    spectral->max_order = max_order;

    /* Clear cached plans and arrays as they depend on max_order */
    g_clear_pointer (&spectral->fftw_plans, g_ptr_array_unref);
    g_clear_pointer (&spectral->cos_arrays, g_ptr_array_unref);
    g_clear_pointer (&spectral->f_vals, fftw_free);
    g_clear_pointer (&spectral->f_vals_tmp, fftw_free);
    g_clear_pointer (&spectral->coeffs_work, fftw_free);

    /* Clear the coefficients array to prevent stale data */
    g_assert_nonnull (spectral->coeffs);
    g_array_set_size (spectral->coeffs, 0);

    {
      const guint N_max = (1 << spectral->max_order) + 1;

      spectral->f_vals      = fftw_malloc (sizeof (gdouble) * N_max);
      spectral->f_vals_tmp  = fftw_malloc (sizeof (gdouble) * N_max);
      spectral->coeffs_work = fftw_malloc (sizeof (gdouble) * N_max);
      spectral->fftw_plans  = g_ptr_array_new_with_free_func (ncm_cfg_fftw_plan_destroy);
      spectral->cos_arrays  = g_ptr_array_new_with_free_func (g_free);
    }
  }
}

/**
 * ncm_spectral_get_max_order:
 * @spectral: a #NcmSpectral
 *
 * Returns: the #NcmSpectral:max-order.
 */
guint
ncm_spectral_get_max_order (NcmSpectral *spectral)
{
  g_return_val_if_fail (NCM_IS_SPECTRAL (spectral), 0);

  return spectral->max_order;
}

static void _ncm_spectral_prepare_plan_for_k (NcmSpectral *spectral, guint k);
static void _ncm_spectral_normalize_coeffs (gdouble *coeffs_work, GArray *coeffs, guint N);

/* Buffers and plans for n_comp components; the plans depend on n_comp */
static void
_ncm_spectral_batch_prepare_buffers (NcmSpectral *spectral, guint n_comp)
{
  const guint N_max = (1 << spectral->max_order) + 1;

  if (spectral->batch_n_comp == n_comp)
    return;

  g_clear_pointer (&spectral->batch_fftw_plans, g_ptr_array_unref);
  g_clear_pointer (&spectral->batch_f_vals, fftw_free);
  g_clear_pointer (&spectral->batch_f_vals_tmp, fftw_free);
  g_clear_pointer (&spectral->batch_coeffs_work, fftw_free);

  spectral->batch_n_comp      = n_comp;
  spectral->batch_f_vals      = fftw_malloc (sizeof (gdouble) * N_max * n_comp);
  spectral->batch_f_vals_tmp  = fftw_malloc (sizeof (gdouble) * N_max * n_comp);
  spectral->batch_coeffs_work = fftw_malloc (sizeof (gdouble) * N_max * n_comp);
  spectral->batch_fftw_plans  = g_ptr_array_new_with_free_func (ncm_cfg_fftw_plan_destroy);
}

/* One DCT-I plan over the n_comp interleaved components, each a stride-n_comp vector */
static void
_ncm_spectral_batch_prepare_plan_for_k (NcmSpectral *spectral, guint k)
{
  const guint N      = (1 << k) + 1;
  const guint n_comp = spectral->batch_n_comp;

  if ((k < spectral->batch_fftw_plans->len) &&
      (g_ptr_array_index (spectral->batch_fftw_plans, k) != NULL))
    return;

  /* The cosine tables are shared with the scalar path, which owns them. */
  _ncm_spectral_prepare_plan_for_k (spectral, k);

  while (spectral->batch_fftw_plans->len <= k)
    g_ptr_array_add (spectral->batch_fftw_plans, NULL);

  {
    const fftw_r2r_kind kind[] = { FFTW_REDFT00 };
    const gint n[]             = { (gint) N };
    gboolean first;
    fftw_plan plan;

    memcpy (spectral->batch_f_vals_tmp, spectral->batch_f_vals,
            sizeof (gdouble) * N * n_comp);

    first = ncm_cfg_fftw_plan_begin ("ncm_spectral_batch_redft00_%u_%u", N, n_comp);

    plan = fftw_plan_many_r2r (1, n, (gint) n_comp,
                               spectral->batch_f_vals, NULL, (gint) n_comp, 1,
                               spectral->batch_coeffs_work, NULL, (gint) n_comp, 1,
                               kind, ncm_cfg_get_fftw_default_flag ());

    ncm_cfg_fftw_plan_end (first);

    memcpy (spectral->batch_f_vals, spectral->batch_f_vals_tmp,
            sizeof (gdouble) * N * n_comp);

    g_ptr_array_index (spectral->batch_fftw_plans, k) = plan;
  }
}

static void
_ncm_spectral_batch_evaluate_all_nodes (NcmSpectral *spectral, NcmSpectralFBatch F,
                                        NcmVector *y, gdouble a, gdouble b, guint k,
                                        gpointer user_data)
{
  const guint N           = (1 << k) + 1;
  const guint n_comp      = spectral->batch_n_comp;
  const gdouble mid       = 0.5 * (a + b);
  const gdouble half_h    = 0.5 * (b - a);
  const gdouble *cos_vals = g_ptr_array_index (spectral->cos_arrays, k);
  guint j, c;

  for (j = 0; j < N; j++)
  {
    const gdouble x = mid + half_h * cos_vals[j];

    F (user_data, x, y);

    for (c = 0; c < n_comp; c++)
      spectral->batch_f_vals[j * n_comp + c] = ncm_vector_get (y, c);
  }
}

/* Level k + 1 contains level k: old values move to even positions, odd ones are new */
static void
_ncm_spectral_batch_refine_to_k (NcmSpectral *spectral, NcmSpectralFBatch F,
                                 NcmVector *y, gdouble a, gdouble b,
                                 guint k_old, guint k_new, gpointer user_data)
{
  const guint N_old       = (1 << k_old) + 1;
  const guint N_new       = (1 << k_new) + 1;
  const guint n_comp      = spectral->batch_n_comp;
  const gdouble mid       = 0.5 * (a + b);
  const gdouble half_h    = 0.5 * (b - a);
  const gdouble *cos_vals = g_ptr_array_index (spectral->cos_arrays, k_new);
  gint j;
  guint jj, c;

  g_assert (k_new == k_old + 1);

  for (j = (gint) N_old - 1; j >= 0; j--)
    memmove (&spectral->batch_f_vals[2 * j * n_comp],
             &spectral->batch_f_vals[j * n_comp],
             sizeof (gdouble) * n_comp);

  for (jj = 1; jj < N_new; jj += 2)
  {
    const gdouble x = mid + half_h * cos_vals[jj];

    F (user_data, x, y);

    for (c = 0; c < n_comp; c++)
      spectral->batch_f_vals[jj * n_comp + c] = ncm_vector_get (y, c);
  }
}

static void
_ncm_spectral_batch_normalize_coeffs (NcmSpectral *spectral, NcmMatrix *coeffs, guint N)
{
  const guint n_comp     = spectral->batch_n_comp;
  const gdouble inv_2Nm1 = 1.0 / (2.0 * (N - 1.0));
  const gdouble inv_Nm1  = 2.0 * inv_2Nm1;
  guint i, c;

  for (c = 0; c < n_comp; c++)
  {
    ncm_matrix_set (coeffs, c, 0, spectral->batch_coeffs_work[c] * inv_2Nm1);
    ncm_matrix_set (coeffs, c, N - 1,
                    spectral->batch_coeffs_work[(N - 1) * n_comp + c] * inv_2Nm1);

    for (i = 1; i < N - 1; i++)
      ncm_matrix_set (coeffs, c, i, spectral->batch_coeffs_work[i * n_comp + c] * inv_Nm1);
  }
}

/* Every component must pass against its own norm, so a small one is still resolved */
static gboolean
_ncm_spectral_batch_check_convergence (NcmMatrix *c_2N, NcmMatrix *c_N, guint N,
                                       guint n_comp, gdouble reltol, gdouble abstol)
{
  guint i, c;

  for (c = 0; c < n_comp; c++)
  {
    gdouble norm2_diff = 0.0;
    gdouble norm2_2N   = 0.0;

    for (i = 0; i < N; i++)
    {
      const gdouble diff = ncm_matrix_get (c_2N, c, i) - ncm_matrix_get (c_N, c, i);

      norm2_diff += diff * diff;
      norm2_2N   += ncm_matrix_get (c_2N, c, i) * ncm_matrix_get (c_2N, c, i);
    }

    if (!(norm2_diff < MAX (reltol * reltol * norm2_2N, abstol * abstol) + 1.0e-100))
      return FALSE;
  }

  return TRUE;
}

/*
 * Whether the next doubling is predicted to fail. The envelope falls by
 * d = |c_{N-1}| / e_{N/2} over the top half of the spectrum, so the N - 1 modes the next
 * level adds carry an l2 mass of about |c_{N-1}| d sqrt(N), compared here with
 * max(reltol ||c||, abstol). The safety factor 10 abandons none of the panels accepted
 * on LSST-Y1 lensing and number counts, a Gaussian and two top-hat shells, at reltol
 * 1e-4 and 1e-6.
 */
#define NCM_SPECTRAL_ABANDON_SAFETY (10.0)

static gboolean
_ncm_spectral_batch_cannot_converge (NcmMatrix *c, guint N, guint n_comp,
                                     gdouble reltol, gdouble abstol)
{
  const guint half = N / 2;
  guint comp, i;

  for (comp = 0; comp < n_comp; comp++)
  {
    const gdouble c_end = fabs (ncm_matrix_get (c, comp, N - 1));
    gdouble env_half    = 0.0;
    gdouble norm2       = 0.0;
    gdouble d, tol_eff;

    for (i = 0; i < N; i++)
    {
      const gdouble a = fabs (ncm_matrix_get (c, comp, i));

      norm2 += a * a;

      if (i >= half)
        env_half = MAX (env_half, a);
    }

    if (env_half <= 0.0)
      continue;

    d       = c_end / env_half;
    tol_eff = MAX (reltol * sqrt (norm2), abstol);

    if (c_end * d * sqrt (N) > NCM_SPECTRAL_ABANDON_SAFETY * tol_eff)
      return TRUE;
  }

  return FALSE;
}

/**
 * ncm_spectral_compute_chebyshev_coeffs_batch_adaptive:
 * @spectral: a #NcmSpectral
 * @F: (scope call): vector-valued function to expand
 * @n_comp: number of components of @F
 * @a: interval lower bound
 * @b: interval upper bound
 * @k_min: starting refinement level
 * @reltol: relative tolerance
 * @abstol: absolute tolerance, in units of the coefficients
 * @coeffs: (out) (transfer full): an @n_comp by $N$ #NcmMatrix of coefficients
 * @user_data: user data for @F
 *
 * Computes the Chebyshev coefficients of the @n_comp components of $F$ on one shared
 * grid, as ncm_spectral_compute_chebyshev_coeffs_adaptive_full(). Each node is evaluated
 * once for all components. Every component must converge, so the one needing the highest
 * level sets it for all. Reaching #NcmSpectral:max-order without converging is an error.
 *
 * Returns: the level $k$ reached, with $N = 2^k + 1$.
 */
guint
ncm_spectral_compute_chebyshev_coeffs_batch_adaptive (NcmSpectral *spectral, NcmSpectralFBatch F,
                                                      guint n_comp, gdouble a, gdouble b,
                                                      guint k_min, gdouble reltol, gdouble abstol,
                                                      NcmMatrix **coeffs, gpointer user_data)
{
  return ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap (spectral, F, n_comp, a, b,
                                                                   k_min, spectral->max_order,
                                                                   reltol, abstol, TRUE,
                                                                   coeffs, user_data);
}

/**
 * ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap:
 * @spectral: a #NcmSpectral
 * @F: (scope call): vector-valued function to expand
 * @n_comp: number of components of @F
 * @a: interval lower bound
 * @b: interval upper bound
 * @k_min: starting refinement level
 * @k_cap: highest refinement level
 * @reltol: relative tolerance
 * @abstol: absolute tolerance, in units of the coefficients
 * @fatal: whether not converging is an error
 * @coeffs: (out) (transfer full): an @n_comp by $N$ #NcmMatrix of coefficients
 * @user_data: user data for @F
 *
 * As ncm_spectral_compute_chebyshev_coeffs_batch_adaptive(), stopping at the smaller of
 * @k_cap and #NcmSpectral:max-order. When @fatal is %FALSE, not converging leaves @coeffs
 * unchanged and returns 0, and the last doubling is skipped when the level below predicts
 * that it cannot converge. That suits a caller that splits its interval on failure.
 *
 * Returns: the level reached, or 0 when @fatal is %FALSE and the expansion did not converge.
 */
guint
ncm_spectral_compute_chebyshev_coeffs_batch_adaptive_cap (NcmSpectral *spectral, NcmSpectralFBatch F,
                                                          guint n_comp, gdouble a, gdouble b,
                                                          guint k_min, guint k_cap,
                                                          gdouble reltol, gdouble abstol,
                                                          gboolean fatal,
                                                          NcmMatrix **coeffs, gpointer user_data)
{
  guint k            = k_min;
  gboolean converged = FALSE;
  NcmMatrix *c_previous, *c_current;
  NcmVector *y;

  k_cap = MIN (k_cap, spectral->max_order);

  g_assert (k_min <= k_cap);
  g_assert_cmpuint (n_comp, >, 0);

  _ncm_spectral_batch_prepare_buffers (spectral, n_comp);

  y          = ncm_vector_new (n_comp);
  c_previous = ncm_matrix_new (n_comp, (1 << spectral->max_order) + 1);
  c_current  = ncm_matrix_new (n_comp, (1 << spectral->max_order) + 1);

  _ncm_spectral_batch_prepare_plan_for_k (spectral, k);
  _ncm_spectral_batch_evaluate_all_nodes (spectral, F, y, a, b, k, user_data);

  {
    const guint N = (1 << k) + 1;

    fftw_execute (g_ptr_array_index (spectral->batch_fftw_plans, k));
    _ncm_spectral_batch_normalize_coeffs (spectral, c_previous, N);
  }

  while (k < k_cap)
  {
    const guint N_prev = (1 << k) + 1;
    guint N;

    /* c_previous holds level k */
    if (!fatal && (k + 1 == k_cap) &&
        _ncm_spectral_batch_cannot_converge (c_previous, N_prev, n_comp,
                                             reltol, abstol))
      break;

    _ncm_spectral_batch_prepare_plan_for_k (spectral, k + 1);
    _ncm_spectral_batch_refine_to_k (spectral, F, y, a, b, k, k + 1, user_data);
    k++;
    N = (1 << k) + 1;

    fftw_execute (g_ptr_array_index (spectral->batch_fftw_plans, k));
    _ncm_spectral_batch_normalize_coeffs (spectral, c_current, N);

    /* Compared over the coefficients of the coarser level */
    if (_ncm_spectral_batch_check_convergence (c_current, c_previous, N_prev,
                                               n_comp, reltol, abstol))
    {
      converged = TRUE;
      break;
    }

    {
      NcmMatrix *tmp = c_previous;

      c_previous = c_current;
      c_current  = tmp;
    }
  }

  if (!converged && fatal)
    g_error ("ncm_spectral_compute_chebyshev_coeffs_batch_adaptive: reached the "
             "maximum order %u (N = %u) without converging to reltol %.3e, "
             "abstol %.3e. max-order is a memory guard, not a stopping rule: "
             "either the tolerances are past what the sampled function carries, "
             "or the interval needs splitting.",
             k_cap, (1 << k_cap) + 1, reltol, abstol);

  if (!converged)
  {
    ncm_matrix_free (c_previous);
    ncm_matrix_free (c_current);
    ncm_vector_free (y);

    return 0;
  }

  {
    const guint N = (1 << k) + 1;
    guint c, i;

    /* Copied, so the result does not keep the max-order buffer alive; a matrix of the
     * right shape is reused */
    if ((*coeffs != NULL) &&
        ((ncm_matrix_nrows (*coeffs) != n_comp) || (ncm_matrix_ncols (*coeffs) != N)))
      ncm_matrix_clear (coeffs);

    if (*coeffs == NULL)
      *coeffs = ncm_matrix_new (n_comp, N);

    for (c = 0; c < n_comp; c++)
      for (i = 0; i < N; i++)
        ncm_matrix_set (*coeffs, c, i, ncm_matrix_get (c_current, c, i));
  }

  ncm_matrix_free (c_previous);
  ncm_matrix_free (c_current);
  ncm_vector_free (y);

  return k;
}

/**
 * ncm_spectral_compute_chebyshev_coeffs:
 * @spectral: a #NcmSpectral
 * @F: (scope call): function to expand
 * @a: left endpoint of the interval
 * @b: right endpoint of the interval
 * @order: number of coefficients $N$
 * @coeffs: (out callee-allocates) (transfer full) (element-type gdouble): the coefficients
 * @user_data: user data for @F
 *
 * Computes the $N$ Chebyshev coefficients of $F$ on $[a, b]$ from its values at the $N$
 * Chebyshev-Lobatto nodes. Aborts if $N < 2$.
 */
void
ncm_spectral_compute_chebyshev_coeffs (NcmSpectral *spectral, NcmSpectralF F, gdouble a, gdouble b, guint order, GArray **coeffs, gpointer user_data)
{
  const gdouble mid    = 0.5 * (a + b);
  const gdouble half_h = 0.5 * (b - a);
  const guint N        = order;
  guint i;

  g_assert_cmpuint (N, >=, 2);

  if (*coeffs == NULL)
    *coeffs = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), N);

  g_array_set_size (*coeffs, N);

  /* Reallocate and replan if N has changed */
  if (spectral->cheb_N_cached != N)
  {
    /* Clean up old resources */
    g_clear_pointer (&spectral->cheb_plan_r2r, ncm_cfg_fftw_plan_destroy);

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

    g_clear_pointer (&spectral->cheb_c_vals, fftw_free);

    /* Allocate new resources */
    spectral->cheb_f_vals   = fftw_malloc (sizeof (gdouble) * N);
    spectral->cheb_c_vals   = fftw_malloc (sizeof (gdouble) * N);
    spectral->cheb_cos_vals = g_new (gdouble, N);

    /* Precompute cosine values at Chebyshev nodes */
    {
      const gdouble inv_Nm1 = 1.0 / (N - 1.0);
      const gdouble pi_Nm1  = M_PI * inv_Nm1;

      for (i = 0; i < N; i++)
        spectral->cheb_cos_vals[i] = cos (pi_Nm1 * i);
    }

    /* Planned and executed on owned buffers: FFTW's new-array execution requires the
     * alignment of the planned arrays, which a caller's array need not have */
    {
      const gboolean first = ncm_cfg_fftw_plan_begin ("ncm_spectral_redft00_%u", N);

      spectral->cheb_plan_r2r = fftw_plan_r2r_1d (N, spectral->cheb_f_vals, spectral->cheb_c_vals,
                                                  FFTW_REDFT00, ncm_cfg_get_fftw_default_flag ());
      ncm_cfg_fftw_plan_end (first);
    }

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

  fftw_execute (spectral->cheb_plan_r2r);
  _ncm_spectral_normalize_coeffs (spectral->cheb_c_vals, *coeffs, N);
}

static void
_ncm_spectral_prepare_plan_for_k (NcmSpectral *spectral, guint k)
{
  const guint N = (1 << k) + 1;
  guint j;

  /* Check if already prepared */
  if ((k < spectral->fftw_plans->len) &&
      (g_ptr_array_index (spectral->fftw_plans, k) != NULL))
    return;

  /* Ensure arrays are large enough */
  while (spectral->fftw_plans->len <= k)
  {
    g_ptr_array_add (spectral->fftw_plans, NULL);
    g_ptr_array_add (spectral->cos_arrays, NULL);
  }

  /* Chebyshev-Lobatto nodes cos(j pi / 2^k) */
  {
    gdouble *cos_vals        = g_new (gdouble, N);
    const gdouble pi_over_2k = M_PI / (1 << k);

    for (j = 0; j < N; j++)
      cos_vals[j] = cos (j * pi_over_2k);

    g_ptr_array_index (spectral->cos_arrays, k) = cos_vals;
  }

  /* Create out-of-place FFTW plan */
  {
    gboolean first;
    fftw_plan plan;

    memcpy (spectral->f_vals_tmp, spectral->f_vals, sizeof (gdouble) * N);

    first = ncm_cfg_fftw_plan_begin ("ncm_spectral_redft00_%u", N);

    plan = fftw_plan_r2r_1d (N,
                             spectral->f_vals,
                             spectral->coeffs_work,
                             FFTW_REDFT00,
                             ncm_cfg_get_fftw_default_flag ());

    ncm_cfg_fftw_plan_end (first);

    memcpy (spectral->f_vals, spectral->f_vals_tmp, sizeof (gdouble) * N);

    g_ptr_array_index (spectral->fftw_plans, k) = plan;
  }
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
  {
    spectral->f_vals[2 * j] = spectral->f_vals[j];
  }

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
_ncm_spectral_check_convergence (GArray *coeffs_2N, GArray *coeffs_N, gdouble tol, gdouble abstol)
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


  if (norm2_diff < MAX (tol * tol * norm2_2N, abstol * abstol) + 1.0e-100)
    return TRUE;

  return FALSE;
}

static guint
_ncm_spectral_compute_chebyshev_coeffs_adaptive_internal (NcmSpectral *spectral, NcmSpectralF F,
                                                          gdouble a, gdouble b, guint k_min,
                                                          gdouble tol, gdouble abstol,
                                                          GArray **coeffs, gpointer user_data,
                                                          gboolean require_convergence,
                                                          guint k_cap,
                                                          gboolean *converged_out)
{
  guint k            = k_min;
  gboolean converged = FALSE;
  GArray *c_previous, *c_current;

  k_cap = MIN (k_cap, spectral->max_order);
  g_assert (k_min <= k_cap);

  if (*coeffs == NULL)
    *coeffs = g_array_new (FALSE, FALSE, sizeof (gdouble));

  g_array_set_size (spectral->coeffs, 0);

  /* Initial evaluation at k_min */
  _ncm_spectral_prepare_plan_for_k (spectral, k);
  _ncm_spectral_evaluate_all_nodes (spectral, F, a, b, k, user_data);

  c_previous = spectral->coeffs;
  c_current  = *coeffs;

  /* Transform using N and store in coeffs_work */
  {
    const guint N  = (1 << k) + 1;
    fftw_plan plan = g_ptr_array_index (spectral->fftw_plans, k);

    /* Transform f_vals -> coeffs_work */
    fftw_execute (plan);

    g_array_set_size (c_previous, N);
    _ncm_spectral_normalize_coeffs (spectral->coeffs_work, c_previous, N);
  }

  while (k < k_cap)
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

      g_array_set_size (c_current, N);
      _ncm_spectral_normalize_coeffs (spectral->coeffs_work, c_current, N);
    }

    if (_ncm_spectral_check_convergence (c_current, c_previous, tol, abstol))
    {
      converged = TRUE;
      break;
    }

    /* Swap c_previous and c_current for next iteration */
    if (k < k_cap)
    {
      GArray *tmp = c_previous;

      c_previous = c_current;
      c_current  = tmp;
    }
  }

  /* max-order bounds memory; the tolerance is what ends the refinement */
  if (require_convergence && !converged)
    g_error ("_ncm_spectral_compute_chebyshev_coeffs_adaptive_internal: "
             "reached max-order %u (N = %u) without converging to tol = %.17g "
             "(abstol = %.17g). Raise max-order, relax tol, or supply an "
             "absolute scale.",
             spectral->max_order, (1 << k) + 1, tol, abstol);

  {
    /* Without a doubling, the only level computed is k_min, in c_previous */
    GArray *c_final = (k == k_min) ? c_previous : c_current;

    if (c_final != *coeffs)
    {
      g_array_set_size (*coeffs, c_final->len);
      memcpy ((*coeffs)->data, c_final->data, sizeof (gdouble) * c_final->len);
    }
  }

  if (converged_out != NULL)
    *converged_out = converged;


  return k;
}

/**
 * ncm_spectral_compute_chebyshev_coeffs_adaptive:
 * @spectral: a #NcmSpectral
 * @F: (scope call): function to expand
 * @a: left endpoint of the interval
 * @b: right endpoint of the interval
 * @k_min: starting refinement level
 * @tol: relative tolerance
 * @coeffs: (out callee-allocates) (transfer full) (element-type gdouble): the coefficients
 * @user_data: user data for @F
 *
 * Computes the Chebyshev coefficients of $F$ on $[a, b]$, doubling the nested
 * Chebyshev-Lobatto grid from $2^{k_\mathrm{min}} + 1$ nodes until, at level $k$,
 * the $\ell_2$ norm of the change of the first $2^{k-1} + 1$ coefficients is below @tol times the norm of the level-$k$ coefficients.
 * Each doubling evaluates $F$ only at the new nodes. Reaching #NcmSpectral:max-order
 * without converging is an error.
 *
 * Returns: the level reached, with $2^k + 1$ coefficients.
 */
guint
ncm_spectral_compute_chebyshev_coeffs_adaptive (NcmSpectral *spectral, NcmSpectralF F,
                                                gdouble a, gdouble b, guint k_min,
                                                gdouble tol, GArray **coeffs, gpointer user_data)
{
  return _ncm_spectral_compute_chebyshev_coeffs_adaptive_internal (spectral, F, a, b, k_min,
                                                                   tol, 0.0, coeffs, user_data,
                                                                   TRUE,
                                                                   spectral->max_order,
                                                                   NULL);
}

/**
 * ncm_spectral_compute_chebyshev_coeffs_adaptive_full:
 * @spectral: a #NcmSpectral
 * @F: (scope call): function to expand
 * @a: left endpoint of the interval
 * @b: right endpoint of the interval
 * @k_min: starting refinement level
 * @reltol: relative tolerance
 * @abstol: absolute tolerance, in units of the coefficients
 * @coeffs: (out callee-allocates) (transfer full) (element-type gdouble): the coefficients
 * @user_data: user data for @F
 *
 * As ncm_spectral_compute_chebyshev_coeffs_adaptive(), with the change compared to the
 * larger of @reltol times the norm and @abstol. An @abstol of 0.0 gives
 * ncm_spectral_compute_chebyshev_coeffs_adaptive(); a positive one stops the refinement of
 * a function known to be negligible.
 *
 * Returns: the level reached, with $2^k + 1$ coefficients.
 */
guint
ncm_spectral_compute_chebyshev_coeffs_adaptive_full (NcmSpectral *spectral, NcmSpectralF F,
                                                     gdouble a, gdouble b, guint k_min,
                                                     gdouble reltol, gdouble abstol,
                                                     GArray **coeffs, gpointer user_data)
{
  return _ncm_spectral_compute_chebyshev_coeffs_adaptive_internal (spectral, F, a, b, k_min,
                                                                   reltol, abstol, coeffs, user_data,
                                                                   TRUE,
                                                                   spectral->max_order,
                                                                   NULL);
}

/**
 * ncm_spectral_compute_chebyshev_coeffs_adaptive_try:
 * @spectral: a #NcmSpectral
 * @F: (scope call): function to expand
 * @a: interval lower bound
 * @b: interval upper bound
 * @k_min: starting refinement level
 * @k_cap: highest refinement level
 * @reltol: relative tolerance
 * @abstol: absolute tolerance, in units of the coefficients
 * @coeffs: (out callee-allocates) (transfer full) (element-type gdouble): the coefficients
 * @user_data: user data for @F
 * @converged: (out): whether the tolerance was met
 *
 * As ncm_spectral_compute_chebyshev_coeffs_adaptive_full(), stopping at the smaller of
 * @k_cap and #NcmSpectral:max-order and reporting the outcome through @converged instead
 * of failing. It uses the single-function buffers, so it may be called from inside the
 * callback of a batch expansion.
 *
 * Returns: the level reached.
 */
guint
ncm_spectral_compute_chebyshev_coeffs_adaptive_try (NcmSpectral *spectral, NcmSpectralF F,
                                                    gdouble a, gdouble b, guint k_min, guint k_cap,
                                                    gdouble reltol, gdouble abstol,
                                                    GArray **coeffs, gpointer user_data,
                                                    gboolean *converged)
{
  return _ncm_spectral_compute_chebyshev_coeffs_adaptive_internal (spectral, F, a, b, k_min,
                                                                   reltol, abstol, coeffs, user_data,
                                                                   FALSE,
                                                                   k_cap,
                                                                   converged);
}

/**
 * ncm_spectral_chebT_to_gegenbauer_alpha1:
 * @c: (element-type gdouble): Chebyshev coefficients
 * @g: (out callee-allocates) (transfer full) (element-type gdouble): the $C^{(1)}_n = U_n$ coefficients
 *
 * Converts the $T_n$ coefficients @c of a series to its $C^{(1)}_n = U_n$ coefficients, of the
 * same length.
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
 * @c: (element-type gdouble): Chebyshev coefficients
 * @g: (out callee-allocates) (transfer full) (element-type gdouble): the $C^{(2)}_k$ coefficients
 *
 * Converts the $T_n$ coefficients @c of a series to its $C^{(2)}_k$ coefficients, of the same
 * length.
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
 * ncm_spectral_chebT_deriv_to_gegenbauer_alpha2:
 * @c: (element-type gdouble): Chebyshev coefficients
 * @g: (out callee-allocates) (transfer full) (element-type gdouble): the $C^{(2)}_k$ coefficients of the derivative
 *
 * Computes the $C^{(2)}_k$ coefficients of $f'(t)$ for $f(t) = \sum_n c_n T_n(t)$,
 * $g_k = c_{k+1} - c_{k+3}$, from $T_n' = n\,U_{n-1}$ and
 * $U_m = (C^{(2)}_m - C^{(2)}_{m-2})/(m + 1)$. For $N$ coefficients in @c, @g has $N - 1$,
 * or one zero when $N \le 1$. The derivative is in $t$; on $[a, b]$ multiply by
 * $2/(b - a)$.
 */
void
ncm_spectral_chebT_deriv_to_gegenbauer_alpha2 (GArray *c, GArray **g)
{
  const guint N    = c->len;
  const guint Nout = (N > 1) ? N - 1 : 1;
  guint k;

  if (*g == NULL)
    *g = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), Nout);

  g_array_set_size (*g, Nout);

  {
    const gdouble *c_data = (gdouble *) c->data;
    gdouble *g_data       = (gdouble *) (*g)->data;

    memset (g_data, 0, Nout * sizeof (gdouble));

    for (k = 0; k + 1 < N; k++)
      g_data[k] = c_data[k + 1] - ((k + 3 < N) ? c_data[k + 3] : 0.0);
  }
}

/**
 * ncm_spectral_chebT_deriv2_to_gegenbauer_alpha2:
 * @c: (element-type gdouble): Chebyshev coefficients
 * @g: (out callee-allocates) (transfer full) (element-type gdouble): the $C^{(2)}_k$ coefficients of the second derivative
 *
 * Computes the $C^{(2)}_k$ coefficients of $f''(t)$ for $f(t) = \sum_n c_n T_n(t)$,
 * $g_k = 2(k + 2)\,c_{k+2}$, from $T_n'' = 2n\,C^{(2)}_{n-2}$. For $N$ coefficients in @c,
 * @g has $N - 2$, or one zero when $N \le 2$. The derivative is in $t$; on $[a, b]$
 * multiply by $(2/(b - a))^2$.
 */
void
ncm_spectral_chebT_deriv2_to_gegenbauer_alpha2 (GArray *c, GArray **g)
{
  const guint N    = c->len;
  const guint Nout = (N > 2) ? N - 2 : 1;
  guint k;

  if (*g == NULL)
    *g = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), Nout);

  g_array_set_size (*g, Nout);

  {
    const gdouble *c_data = (gdouble *) c->data;
    gdouble *g_data       = (gdouble *) (*g)->data;

    memset (g_data, 0, Nout * sizeof (gdouble));

    for (k = 0; k + 2 < N; k++)
      g_data[k] = 2.0 * (k + 2.0) * c_data[k + 2];
  }
}

/**
 * ncm_spectral_gegenbauer_alpha2_xmul:
 * @g: (element-type gdouble): $C^{(2)}_n$ coefficients
 * @alpha: linear coefficient of the factor
 * @beta: constant coefficient of the factor
 * @out: (out callee-allocates) (transfer full) (element-type gdouble): the $C^{(2)}_n$ coefficients of the product
 *
 * Multiplies $\sum_n g_n C^{(2)}_n(t)$ by $\alpha t + \beta$, using
 * $t\,C^{(2)}_n = [(n + 1)\,C^{(2)}_{n+1} + (n + 3)\,C^{(2)}_{n-1}]/(2(n + 2))$. For $N$
 * coefficients in @g, @out has $N + 1$. @out must not alias @g.
 */
void
ncm_spectral_gegenbauer_alpha2_xmul (GArray *g, gdouble alpha, gdouble beta, GArray **out)
{
  const guint N    = g->len;
  const guint Nout = N + 1;
  guint n;

  g_assert (*out != g);

  if (*out == NULL)
    *out = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), Nout);

  g_array_set_size (*out, Nout);

  {
    const gdouble *g_data = (gdouble *) g->data;
    gdouble *out_data     = (gdouble *) (*out)->data;

    memset (out_data, 0, Nout * sizeof (gdouble));

    for (n = 0; n < N; n++)
    {
      const gdouble nd  = (gdouble) n;
      const gdouble gn  = g_data[n];
      const gdouble den = 2.0 * (nd + 2.0);

      out_data[n + 1] += alpha * gn * (nd + 1.0) / den;
      out_data[n]     += beta * gn;

      if (n > 0)
        out_data[n - 1] += alpha * gn * (nd + 3.0) / den;
    }
  }
}

/**
 * ncm_spectral_chebyshev_rebase:
 * @spectral: a #NcmSpectral
 * @c: (element-type gdouble): Chebyshev coefficients on [@a_in, @b_in]
 * @len: number of leading coefficients of @c to use, 0 for all
 * @a_in: left endpoint of the interval of @c
 * @b_in: right endpoint of the interval of @c
 * @a_out: left endpoint of the target interval
 * @b_out: right endpoint of the target interval
 * @rebased: (out callee-allocates) (transfer full) (element-type gdouble): the
 *   coefficients on [@a_out, @b_out]
 *
 * Expresses the same polynomial as a Chebyshev series on [@a_out, @b_out]. The argument on
 * [@a_in, @b_in] is $s = \alpha t + \beta$ in terms of the one on [@a_out, @b_out], and each
 * $T_k(s)$ is expanded in the $T_j(t)$ by the Chebyshev recurrence, at a cost of $O(n^2)$.
 *
 * The target need not lie inside the source interval; outside it the result continues the
 * polynomial, and $T_k(s)$ grows as $(|s| + \sqrt{s^2 - 1})^k$ there. The returned
 * $\sum_j |b_j|$ bounds $|f|$ on the target interval and shows when that growth has
 * amplified roundoff. The scratch space belongs to @spectral.
 *
 * Returns: $\sum_j |b_j|$ over the rebased coefficients $b_j$, or infinity if one is not
 *   finite.
 */
gdouble
ncm_spectral_chebyshev_rebase (NcmSpectral *spectral, GArray *c, guint len,
                               gdouble a_in, gdouble b_in,
                               gdouble a_out, gdouble b_out,
                               GArray **rebased)
{
  const guint n       = (len == 0) ? c->len : len;
  const gdouble alpha = (b_out - a_out) / (b_in - a_in);
  const gdouble beta  = (b_out + a_out - b_in - a_in) / (b_in - a_in);
  const gdouble *a    = (const gdouble *) c->data;
  gdouble *previous, *current, *next, *b;
  gdouble norm = 0.0;
  guint degree, i;

  g_assert_cmpuint (len, <=, c->len);
  g_assert_cmpfloat (b_in, >, a_in);
  g_assert_cmpfloat (b_out, >, a_out);

  if (*rebased == NULL)
    *rebased = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n);

  g_array_set_size (*rebased, n);

  if (n == 0)
    return 0.0;

  if (spectral->rebase_work_len < 3 * n)
  {
    spectral->rebase_work     = g_realloc_n (spectral->rebase_work, 3 * n, sizeof (gdouble));
    spectral->rebase_work_len = 3 * n;
  }

  previous = spectral->rebase_work;
  current  = previous + n;
  next     = current + n;
  memset (spectral->rebase_work, 0, 3 * n * sizeof (gdouble));

  b = (gdouble *) (*rebased)->data;
  memset (b, 0, n * sizeof (gdouble));

  previous[0] = 1.0;
  b[0]        = a[0];

  if (n > 1)
  {
    current[0] = beta;
    current[1] = alpha;
    b[0]      += a[1] * beta;
    b[1]      += a[1] * alpha;
  }

  /* Recursively form T_degree(alpha t + beta) in the T_k(t) basis. */
  for (degree = 1; degree + 1 < n; degree++)
  {
    gdouble *tmp;

    memset (next, 0, n * sizeof (gdouble));

    for (i = 0; i <= degree; i++)
    {
      next[i] += 2.0 * beta * current[i] - previous[i];

      if (i == 0)
      {
        next[1] += 2.0 * alpha * current[0];
      }
      else
      {
        next[i - 1] += alpha * current[i];
        next[i + 1] += alpha * current[i];
      }
    }

    for (i = 0; i <= degree + 1; i++)
      b[i] += a[degree + 1] * next[i];

    tmp      = previous;
    previous = current;
    current  = next;
    next     = tmp;
  }

  for (i = 0; i < n; i++)
  {
    if (!isfinite (b[i]))
      return HUGE_VAL;

    norm += fabs (b[i]);
  }

  return norm;
}

/**
 * ncm_spectral_gegenbauer_alpha1_eval:
 * @c: (element-type gdouble): $C^{(1)}_n$ coefficients
 * @t: the point, in $[-1, 1]$
 *
 * Evaluates the series by the forward recurrence of $C^{(1)}_n = U_n$, with the closed
 * form $U_n(\pm 1) = (\pm 1)^n (n + 1)$ at the endpoints.
 *
 * Returns: $\sum_n c_n C^{(1)}_n(t)$.
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
 * @c: (element-type gdouble): $C^{(2)}_n$ coefficients
 * @t: the point, in $[-1, 1]$
 *
 * Evaluates the series by the forward recurrence of $C^{(2)}_n$, with the closed form
 * $C^{(2)}_n(\pm 1) = (\pm 1)^n \binom{n+3}{3}$ at the endpoints.
 *
 * Returns: $\sum_n c_n C^{(2)}_n(t)$.
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
 * @a: (element-type gdouble): Chebyshev coefficients
 * @t: the point, in $[-1, 1]$
 *
 * Evaluates the series by the Clenshaw recurrence for $|t| < 0.9$ and by Reinsch's
 * modification of it closer to the endpoints.
 *
 * Returns: $\sum_k a_k T_k(t)$.
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
 * @a: (element-type gdouble): Chebyshev coefficients
 * @t: the point, in $[-1, 1]$
 *
 * Evaluates the derivative in $t$ in one backward pass that builds the coefficients of the
 * derivative series and sums them by the Clenshaw recurrence.
 *
 * Returns: $\sum_k a_k T_k'(t)$.
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
 * @c: (element-type gdouble): $C^{(1)}_n$ coefficients
 * @a: left endpoint of the interval
 * @b: right endpoint of the interval
 * @x: the point, in $[a, b]$
 *
 * Same as ncm_spectral_gegenbauer_alpha1_eval() at $t$ given by ncm_spectral_x_to_t().
 *
 * Returns: $\sum_n c_n C^{(1)}_n(t)$.
 */
gdouble
ncm_spectral_gegenbauer_alpha1_eval_x (GArray *c, gdouble a, gdouble b, gdouble x)
{
  const gdouble t = ncm_spectral_x_to_t (a, b, x);

  return ncm_spectral_gegenbauer_alpha1_eval (c, t);
}

/**
 * ncm_spectral_gegenbauer_alpha2_eval_x:
 * @c: (element-type gdouble): $C^{(2)}_n$ coefficients
 * @a: left endpoint of the interval
 * @b: right endpoint of the interval
 * @x: the point, in $[a, b]$
 *
 * Same as ncm_spectral_gegenbauer_alpha2_eval() at $t$ given by ncm_spectral_x_to_t().
 *
 * Returns: $\sum_n c_n C^{(2)}_n(t)$.
 */
gdouble
ncm_spectral_gegenbauer_alpha2_eval_x (GArray *c, gdouble a, gdouble b, gdouble x)
{
  const gdouble t = ncm_spectral_x_to_t (a, b, x);

  return ncm_spectral_gegenbauer_alpha2_eval (c, t);
}

/**
 * ncm_spectral_chebyshev_eval_x:
 * @a: (element-type gdouble): Chebyshev coefficients
 * @a_v: left endpoint of the interval
 * @b: right endpoint of the interval
 * @x: the point, in [@a_v, @b]
 *
 * Same as ncm_spectral_chebyshev_eval() at $t$ given by ncm_spectral_x_to_t().
 *
 * Returns: $\sum_k a_k T_k(t)$.
 */
gdouble
ncm_spectral_chebyshev_eval_x (GArray *a, gdouble a_v, gdouble b, gdouble x)
{
  const gdouble t = ncm_spectral_x_to_t (a_v, b, x);

  return ncm_spectral_chebyshev_eval (a, t);
}

/**
 * ncm_spectral_chebyshev_deriv_x:
 * @a: (element-type gdouble): Chebyshev coefficients
 * @a_v: left endpoint of the interval
 * @b: right endpoint of the interval
 * @x: the point, in [@a_v, @b]
 *
 * Same as ncm_spectral_chebyshev_deriv() at $t$ given by ncm_spectral_x_to_t(), times
 * $\mathrm{d}t/\mathrm{d}x = 2/(b - a_v)$.
 *
 * Returns: $\mathrm{d}f/\mathrm{d}x$ at @x.
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
 * ncm_spectral_chebyshev_integrate:
 * @a: (element-type gdouble): Chebyshev coefficients
 * @a_v: left endpoint of the interval
 * @b: right endpoint of the interval
 *
 * Integrates the series over its interval, from $\int_{-1}^{1} T_k(t)\,\mathrm{d}t = 2/(1 - k^2)$
 * for even $k$ and zero for odd $k$. Applied to the coefficients of
 * ncm_spectral_compute_chebyshev_coeffs() or its adaptive variants, this is Clenshaw-Curtis
 * quadrature.
 *
 * Returns: $\int_{a_v}^{b} \sum_k a_k T_k(t(x))\,\mathrm{d}x$.
 */
gdouble
ncm_spectral_chebyshev_integrate (GArray *a, gdouble a_v, gdouble b)
{
  const gdouble *a_data = (gdouble *) a->data;
  gdouble sum           = 0.0;
  guint k;

  for (k = 0; k < a->len; k += 2)
    sum += a_data[k] / (1.0 - (gdouble) k * (gdouble) k);

  return (b - a_v) * sum;
}

/**
 * ncm_spectral_get_proj_matrix:
 * @N: size of the matrix
 *
 * Builds the $N \times N$ matrix of the identity from the $T_n$ coefficients of $f$ to the
 * $C^{(2)}_k$ coefficients of $f$, truncated to the first $N$ columns.
 *
 * Returns: (transfer full): the operator matrix.
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
 * Builds the $N \times N$ matrix of multiplication by $x$ from the $T_n$ coefficients of $f$ to the
 * $C^{(2)}_k$ coefficients of $x f$, truncated to the first $N$ columns.
 *
 * Returns: (transfer full): the operator matrix.
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
 * Builds the $N \times N$ matrix of multiplication by $x^2$ from the $T_n$ coefficients of $f$ to the
 * $C^{(2)}_k$ coefficients of $x^2 f$, truncated to the first $N$ columns.
 *
 * Returns: (transfer full): the operator matrix.
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
 * Builds the $N \times N$ matrix of the derivative from the $T_n$ coefficients of $f$ to the
 * $C^{(2)}_k$ coefficients of $f'$, truncated to the first $N$ columns.
 *
 * Returns: (transfer full): the operator matrix.
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
 * Builds the $N \times N$ matrix of $x\,\mathrm{d}/\mathrm{d}x$ from the $T_n$ coefficients of $f$ to the
 * $C^{(2)}_k$ coefficients of $x f'$, truncated to the first $N$ columns.
 *
 * Returns: (transfer full): the operator matrix.
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
 * Builds the $N \times N$ matrix of the second derivative from the $T_n$ coefficients of $f$ to the
 * $C^{(2)}_k$ coefficients of $f''$, truncated to the first $N$ columns.
 *
 * Returns: (transfer full): the operator matrix.
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
 * Builds the $N \times N$ matrix of $x\,\mathrm{d}^2/\mathrm{d}x^2$ from the $T_n$ coefficients of $f$ to the
 * $C^{(2)}_k$ coefficients of $x f''$, truncated to the first $N$ columns.
 *
 * Returns: (transfer full): the operator matrix.
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
 * Builds the $N \times N$ matrix of $x^2\,\mathrm{d}^2/\mathrm{d}x^2$ from the $T_n$ coefficients of $f$ to the
 * $C^{(2)}_k$ coefficients of $x^2 f''$, truncated to the first $N$ columns.
 *
 * Returns: (transfer full): the operator matrix.
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
 * @row_data: the row
 * @k: row index
 * @offset: position of column @k in @row_data
 * @coeff: factor of the row
 *
 * Adds @coeff times row $k$ of ncm_spectral_get_proj_matrix() to @row_data. Its nonzero
 * entries are
 *
 * - column $k$: $1/(2(k + 1))$, plus $1/2$ when $k = 0$;
 * - column $k + 2$: $-(k + 2)/((k + 1)(k + 3))$;
 * - column $k + 4$: $1/(2(k + 3))$.
 */

/**
 * ncm_spectral_compute_x_row:
 * @row_data: the row
 * @k: row index
 * @offset: position of column @k in @row_data
 * @coeff: factor of the row
 *
 * Adds @coeff times row $k$ of ncm_spectral_get_x_matrix() to @row_data. Its nonzero
 * entries are
 *
 * - column $k - 1$, for $k \ge 1$: $1/(4(k + 1))$, plus $1/8$ when $k = 1$;
 * - column $k + 1$: $-1/(4(k + 3))$, plus $1/4$ when $k = 0$;
 * - column $k + 3$: $-1/(4(k + 1))$;
 * - column $k + 5$: $1/(4(k + 3))$.
 */

/**
 * ncm_spectral_compute_x2_row:
 * @row_data: the row
 * @k: row index
 * @offset: position of column @k in @row_data
 * @coeff: factor of the row
 *
 * Adds @coeff times row $k$ of ncm_spectral_get_x2_matrix() to @row_data. Its nonzero
 * entries are
 *
 * - column $k - 2$, for $k \ge 2$: $1/(8(k + 1))$, plus $1/24$ when $k = 2$;
 * - column $k$: $1/(4(k + 1)(k + 3))$, plus $1/12$ when $k = 0$ and $1/16$ when $k = 1$;
 * - column $k + 2$: $-(k + 2)/(4(k + 1)(k + 3))$, plus $1/8$ when $k = 0$;
 * - column $k + 4$: $-1/(4(k + 1)(k + 3))$;
 * - column $k + 6$: $1/(8(k + 3))$.
 */

/**
 * ncm_spectral_compute_d_row:
 * @row_data: the row
 * @offset: position of column @k in @row_data
 * @coeff: factor of the row
 *
 * Adds @coeff times row $k$ of ncm_spectral_get_d_matrix() to @row_data. Its nonzero
 * entries are
 *
 * - column $k + 1$: $1$;
 * - column $k + 3$: $-1$.
 */

/**
 * ncm_spectral_compute_x_d_row:
 * @row_data: the row
 * @k: row index
 * @offset: position of column @k in @row_data
 * @coeff: factor of the row
 *
 * Adds @coeff times row $k$ of ncm_spectral_get_x_d_matrix() to @row_data. Its nonzero
 * entries are
 *
 * - column $k$: $k/(2(k + 1))$;
 * - column $k + 2$: $(k + 2)/((k + 1)(k + 3))$;
 * - column $k + 4$: $-(k + 4)/(2(k + 3))$.
 */

/**
 * ncm_spectral_compute_d2_row:
 * @row_data: the row
 * @k: row index
 * @offset: position of column @k in @row_data
 * @coeff: factor of the row
 *
 * Adds @coeff times row $k$ of ncm_spectral_get_d2_matrix() to @row_data. Its nonzero
 * entries are
 *
 * - column $k + 2$: $2(k + 2)$.
 */

/**
 * ncm_spectral_compute_x_d2_row:
 * @row_data: the row
 * @k: row index
 * @offset: position of column @k in @row_data
 * @coeff: factor of the row
 *
 * Adds @coeff times row $k$ of ncm_spectral_get_x_d2_matrix() to @row_data. Its nonzero
 * entries are
 *
 * - column $k + 1$: $k$;
 * - column $k + 3$: $k + 4$.
 */

/**
 * ncm_spectral_compute_x2_d2_row:
 * @row_data: the row
 * @k: row index
 * @offset: position of column @k in @row_data
 * @coeff: factor of the row
 *
 * Adds @coeff times row $k$ of ncm_spectral_get_x2_d2_matrix() to @row_data. Its nonzero
 * entries are
 *
 * - column $k$: $k(k - 1)/(2(k + 1))$;
 * - column $k + 2$: $(k + 2)((k + 2)^2 - 3)/((k + 1)(k + 3))$;
 * - column $k + 4$: $(k + 4)(k + 5)/(2(k + 3))$.
 */

