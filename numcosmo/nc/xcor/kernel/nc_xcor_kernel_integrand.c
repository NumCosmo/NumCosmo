/***************************************************************************
 *            nc_xcor_kernel_integrand.c
 *
 *  Tue July 14 12:00:00 2015
 *  Copyright  2015  Cyrille Doux
 *  <cdoux@apc.in2p3.fr>
 *  Thu October 08 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2015 Cyrille Doux <cdoux@apc.in2p3.fr>
 * Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 *
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


/*
 * NcXcorKernelIntegrand: the Chebyshev closure of W_l(k) of an ell block, as the k
 * integrals see it. A record of the panels, one row per multipole, with each row's
 * scale and support; the builder in nc_xcor_kernel_closure.c fills it.
 */

/**
 * NcXcorKernelIntegrand:
 *
 * The closure of $W_\ell(k)$ over a block of multipoles: a Chebyshev series on panels,
 * built by nc_xcor_kernel_get_eval_vectorized_full(). The panels hold
 * nc_xcor_kernel_integrand_get_n_terms() groups of one row per multipole, and multipole
 * $i$ of a block of $n$ is the sum of the panel rows $i + t n$ over the terms $t$. Panel
 * row $r$ is a function of $x = k / s_r$, with $s_r$ the scale of
 * nc_xcor_kernel_integrand_get_scale(), and is zero outside its support
 * nc_xcor_kernel_integrand_get_range_comp().
 *
 * Evaluation reads the panels and writes nothing, so one integrand may be evaluated by
 * several threads at once.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/algebra/ncm_spectral_panels.h"
#include "nc/xcor/nc_xcor_kernel.h"
#include "nc/xcor/kernel/nc_xcor_kernel_private.h"

G_DEFINE_BOXED_TYPE (NcXcorKernelIntegrand, nc_xcor_kernel_integrand, nc_xcor_kernel_integrand_ref, nc_xcor_kernel_integrand_unref)

/*
 * A closure of @len multipoles in @n_terms terms with empty panels and unset rows, for
 * the builder to fill: the panels, every panel row and the block's k range.
 */
NcXcorKernelIntegrand *
_nc_xcor_kernel_integrand_new (guint len, guint n_terms)
{
  NcXcorKernelIntegrand *integrand = g_new0 (NcXcorKernelIntegrand, 1);
  const guint n_rows               = len * n_terms;

  integrand->refcount = 1;
  integrand->len      = len;
  integrand->n_terms  = n_terms;
  integrand->sp       = ncm_spectral_panels_new (n_rows);
  integrand->rows     = g_array_sized_new (FALSE, TRUE, sizeof (NcXcorKernelIntegrandRow), n_rows);
  integrand->scales   = ncm_vector_new (n_rows);
  integrand->windows  = NULL;

  g_array_set_size (integrand->rows, n_rows);

  return integrand;
}

/*
 * Builds, for a block of several terms, each multipole's window as one row in k: the sum
 * of its term rows by ncm_spectral_panels_new_sum_rows() on their supports and scales,
 * trailing coefficients below @abstol dropped. The k integrals of nc_xcor_kquad.c read
 * these rows; nc_xcor_kernel_integrand_eval_comps() reads the term rows.
 */
void
_nc_xcor_kernel_integrand_sum_terms (NcXcorKernelIntegrand *integrand, NcmSpectral *spectral, gdouble abstol)
{
  const guint len     = integrand->len;
  const guint n_terms = integrand->n_terms;
  NcmSpectralPanels *sps[MAX_LIMBER_TERMS];
  guint rows[MAX_LIMBER_TERMS];
  gdouble scales[MAX_LIMBER_TERMS];
  gdouble x_min[MAX_LIMBER_TERMS];
  gdouble x_max[MAX_LIMBER_TERMS];
  guint i, t;

  if (n_terms > MAX_LIMBER_TERMS)
    g_error ("_nc_xcor_kernel_integrand_sum_terms: %u terms, at most %d fit.", n_terms, MAX_LIMBER_TERMS);

  g_clear_pointer (&integrand->windows, g_ptr_array_unref);
  integrand->windows = g_ptr_array_new_full (len, (GDestroyNotify) ncm_spectral_panels_free);

  for (i = 0; i < len; i++)
  {
    for (t = 0; t < n_terms; t++)
    {
      const guint r                       = i + t * len;
      const NcXcorKernelIntegrandRow *row = &g_array_index (integrand->rows, NcXcorKernelIntegrandRow, r);

      sps[t]    = integrand->sp;
      rows[t]   = r;
      scales[t] = ncm_vector_get (integrand->scales, r);
      x_min[t]  = row->x_min;
      x_max[t]  = row->x_max;
    }

    g_ptr_array_add (integrand->windows, ncm_spectral_panels_new_sum_rows (spectral, n_terms, sps, rows, scales, x_min, x_max, abstol));
  }
}

static gdouble _nc_xcor_kernel_integrand_cell_bound (NcmSpectralPanels *win, guint p);
static gdouble _nc_xcor_kernel_integrand_cell_peak (NcmSpectralPanels *win, guint p);

/*
 * The smallest over the multipoles of a block of several terms of the peak of |W_i(k)|, from
 * the summed windows: each cell's polynomial sampled at 2 N + 1 Chebyshev-Lobatto points, N
 * its coefficient count. The cell of the largest bound sum_j |c_j| is sampled first, and a
 * cell whose bound is below the largest value found is skipped. Multipoles whose window is
 * zero everywhere are left out, and a block whose windows all are gives 0.
 */
gdouble
_nc_xcor_kernel_integrand_window_peak_min (NcXcorKernelIntegrand *integrand)
{
  gdouble peak_min = G_MAXDOUBLE;
  guint i, p;

  g_assert_nonnull (integrand->windows);

  for (i = 0; i < integrand->len; i++)
  {
    NcmSpectralPanels *win = g_ptr_array_index (integrand->windows, i);
    const guint n_panels   = ncm_spectral_panels_get_n_panels (win);
    gdouble bound_max      = -1.0;
    gdouble peak           = 0.0;
    guint p_max            = 0;

    for (p = 0; p < n_panels; p++)
    {
      const gdouble bound = _nc_xcor_kernel_integrand_cell_bound (win, p);

      if (bound > bound_max)
      {
        bound_max = bound;
        p_max     = p;
      }
    }

    if (n_panels > 0)
      peak = _nc_xcor_kernel_integrand_cell_peak (win, p_max);

    for (p = 0; p < n_panels; p++)
      if ((p != p_max) && (_nc_xcor_kernel_integrand_cell_bound (win, p) > peak))
        peak = GSL_MAX (peak, _nc_xcor_kernel_integrand_cell_peak (win, p));

    if (peak > 0.0)
      peak_min = GSL_MIN (peak_min, peak);
  }

  return (peak_min < G_MAXDOUBLE) ? peak_min : 0.0;
}

/* sum_j |c_j| of cell p, a bound on |W| over the cell */
static gdouble
_nc_xcor_kernel_integrand_cell_bound (NcmSpectralPanels *win, guint p)
{
  gdouble a, b, bound = 0.0;
  guint N, j;
  const gdouble *c = ncm_spectral_panels_peek_panel (win, p, &a, &b, &N);

  for (j = 0; j < N; j++)
    bound += fabs (c[j]);

  return bound;
}

/* The largest |W| over 2 N + 1 Chebyshev-Lobatto points of cell p */
static gdouble
_nc_xcor_kernel_integrand_cell_peak (NcmSpectralPanels *win, guint p)
{
  gdouble a, b, peak = 0.0;
  guint N, j;

  ncm_spectral_panels_peek_panel (win, p, &a, &b, &N);

  for (j = 0; j < 2 * N + 1; j++)
  {
    const gdouble x = 0.5 * (a + b) - 0.5 * (b - a) * cos (M_PI * j / (2.0 * N));

    peak = GSL_MAX (peak, fabs (ncm_spectral_panels_eval_row_on (win, p, 0, x)));
  }

  return peak;
}

/**
 * nc_xcor_kernel_integrand_ref:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Increases the reference count of @integrand by one atomically.
 *
 * Returns: (transfer full): @integrand
 */
NcXcorKernelIntegrand *
nc_xcor_kernel_integrand_ref (NcXcorKernelIntegrand *integrand)
{
  g_atomic_int_inc (&integrand->refcount);

  return integrand;
}

/**
 * nc_xcor_kernel_integrand_unref:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Decreases the reference count of @integrand by one atomically, and frees it when
 * the count reaches zero.
 */
void
nc_xcor_kernel_integrand_unref (NcXcorKernelIntegrand *integrand)
{
  if (g_atomic_int_dec_and_test (&integrand->refcount))
  {
    ncm_spectral_panels_clear (&integrand->sp);
    g_array_unref (integrand->rows);
    ncm_vector_clear (&integrand->scales);
    g_clear_pointer (&integrand->windows, g_ptr_array_unref);
    g_free (integrand);
  }
}

/**
 * nc_xcor_kernel_integrand_clear:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * If *@integrand is not %NULL, decreases its reference count and sets the pointer to
 * %NULL.
 */
void
nc_xcor_kernel_integrand_clear (NcXcorKernelIntegrand **integrand)
{
  if (*integrand != NULL)
  {
    nc_xcor_kernel_integrand_unref (*integrand);
    *integrand = NULL;
  }
}

/**
 * nc_xcor_kernel_integrand_get_len:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Returns: the number of rows, one per multipole of the block
 */
guint
nc_xcor_kernel_integrand_get_len (NcXcorKernelIntegrand *integrand)
{
  return integrand->len;
}

/**
 * nc_xcor_kernel_integrand_get_n_terms:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Gets the number of terms: the panels hold this many groups of
 * nc_xcor_kernel_integrand_get_len() rows, and each multipole is the sum of its row in
 * every group.
 *
 * Returns: the number of terms
 */
guint
nc_xcor_kernel_integrand_get_n_terms (NcXcorKernelIntegrand *integrand)
{
  return integrand->n_terms;
}

/**
 * nc_xcor_kernel_integrand_peek_window:
 * @integrand: a #NcXcorKernelIntegrand
 * @i: a multipole of the block
 *
 * Peeks the window of multipole @i as one row in $k$: with several terms, the sum of its
 * term rows as a one-row panel set, trailing coefficients below a hundredth of the fit's
 * absolute floor dropped; the $k$ integrals read it. With one term the rows are the
 * windows and there is none.
 *
 * Returns: (transfer none) (nullable): the window of multipole @i as one row, or %NULL with one term
 */
NcmSpectralPanels *
nc_xcor_kernel_integrand_peek_window (NcXcorKernelIntegrand *integrand, guint i)
{
  if (integrand->windows == NULL)
    return NULL;

  g_assert_cmpuint (i, <, integrand->len);

  return g_ptr_array_index (integrand->windows, i);
}

/**
 * nc_xcor_kernel_integrand_get_range:
 * @integrand: a #NcXcorKernelIntegrand
 * @k_min: (out): minimum k value
 * @k_max: (out): maximum k value
 *
 * Gets the k range of the block, the union of the rows' ranges.
 */
void
nc_xcor_kernel_integrand_get_range (NcXcorKernelIntegrand *integrand, gdouble *k_min, gdouble *k_max)
{
  *k_min = integrand->k_min;
  *k_max = integrand->k_max;
}

/**
 * nc_xcor_kernel_integrand_get_range_comp:
 * @integrand: a #NcXcorKernelIntegrand
 * @i: panel row index
 * @k_min: (out): minimum k value
 * @k_max: (out): maximum k value
 *
 * Gets the k range panel row @i is supported on, which can be a part of the range
 * nc_xcor_kernel_integrand_get_range() reports for the block: under the Limber
 * approximation each multipole vanishes outside its own band.
 */
void
nc_xcor_kernel_integrand_get_range_comp (NcXcorKernelIntegrand *integrand, guint i, gdouble *k_min, gdouble *k_max)
{
  const NcXcorKernelIntegrandRow *row = &g_array_index (integrand->rows, NcXcorKernelIntegrandRow, i);

  *k_min = row->k_min;
  *k_max = row->k_max;
}

/**
 * nc_xcor_kernel_integrand_get_scale:
 * @integrand: a #NcXcorKernelIntegrand
 * @i: panel row index
 *
 * The panels are in the variable $x = k / s_i$, with $s_i$ the scale of panel row @i:
 * $s_i = 1$ for a closure built in $k$, and $s_i = \nu_i = \ell_i + 1/2$ for a Limber
 * closure built in $u = k / \nu$, whose panels are then common to the block.
 *
 * Returns: the scale of panel row @i
 */
gdouble
nc_xcor_kernel_integrand_get_scale (NcXcorKernelIntegrand *integrand, guint i)
{
  return ncm_vector_get (integrand->scales, i);
}

/**
 * nc_xcor_kernel_integrand_peek_scales:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Peeks the scales of every panel row, see nc_xcor_kernel_integrand_get_scale().
 *
 * Returns: (transfer none): the scales, one per panel row
 */
NcmVector *
nc_xcor_kernel_integrand_peek_scales (NcXcorKernelIntegrand *integrand)
{
  return integrand->scales;
}

/**
 * nc_xcor_kernel_integrand_peek_panels:
 * @integrand: a #NcXcorKernelIntegrand
 *
 * Peeks the panels, nc_xcor_kernel_integrand_get_n_terms() groups of one row per
 * multipole, row $i$ in the variable $x = k / s_i$ of
 * nc_xcor_kernel_integrand_get_scale().
 *
 * Returns: (transfer none): the panels
 */
NcmSpectralPanels *
nc_xcor_kernel_integrand_peek_panels (NcXcorKernelIntegrand *integrand)
{
  return integrand->sp;
}

/**
 * nc_xcor_kernel_integrand_eval: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @k: wavenumber
 * @W: (array) (out caller-allocates): array of length nc_xcor_kernel_integrand_get_len()
 *
 * Evaluates every row at wavenumber @k.
 */
void
nc_xcor_kernel_integrand_eval (NcXcorKernelIntegrand *integrand, gdouble k, gdouble *W)
{
  nc_xcor_kernel_integrand_eval_comps (integrand, k, 0, integrand->len, W);
}

/**
 * nc_xcor_kernel_integrand_eval_comps: (skip)
 * @integrand: a #NcXcorKernelIntegrand
 * @k: wavenumber
 * @offset: index of the first row to evaluate
 * @len: number of rows to evaluate
 * @W: (array) (out caller-allocates): full-length array to store results in
 *
 * Evaluates multipoles [@offset, @offset + @len) at wavenumber @k, the sum of their
 * terms, writing them at their own indices in @W and leaving the rest of @W untouched. A
 * panel row is zero outside its support.
 */
void
nc_xcor_kernel_integrand_eval_comps (NcXcorKernelIntegrand *integrand, gdouble k, guint offset, guint len, gdouble *W)
{
  guint i, t;

  /* In k, the support of every row is the panels' range. */
  if (!integrand->scaled)
  {
    gdouble x_min, x_max;

    ncm_spectral_panels_get_range (integrand->sp, &x_min, &x_max);

    for (i = offset; i < offset + len; i++)
    {
      W[i] = 0.0;

      if ((k >= x_min) && (k <= x_max))
        for (t = 0; t < integrand->n_terms; t++)
          W[i] += ncm_spectral_panels_eval_row (integrand->sp, i + t * integrand->len, k);
    }

    return;
  }

  /* With per-row scales each row lives in its own x = k / scale and is nonzero on its own
   * support [x_min, x_max], whose ends are panel edges. Whether k is inside is decided in
   * k, against that support times the scale, and x is then clamped to the support: k / scale
   * can round past an end that k reaches exactly. The upper end is read from inside the
   * support, on the panel that ends there, where the panel that starts there may already be
   * past the row's support. */
  for (i = offset; i < offset + len; i++)
  {
    W[i] = 0.0;

    for (t = 0; t < integrand->n_terms; t++)
    {
      const guint r                       = i + t * integrand->len;
      const NcXcorKernelIntegrandRow *row = &g_array_index (integrand->rows, NcXcorKernelIntegrandRow, r);

      if ((k >= row->k_min) && (k <= row->k_max))
      {
        const gdouble x = GSL_MIN (GSL_MAX (k / ncm_vector_fast_get (integrand->scales, r), row->x_min), row->x_max);
        guint p         = ncm_spectral_panels_find_panel (integrand->sp, x);

        if ((x == row->x_max) && (p > 0))
        {
          gdouble a, b;
          guint N;

          ncm_spectral_panels_peek_panel (integrand->sp, p, &a, &b, &N);

          if (a == x)
            p--;
        }

        W[i] += ncm_spectral_panels_eval_row_on (integrand->sp, p, r, x);
      }
    }
  }
}

/**
 * nc_xcor_kernel_integrand_eval_array:
 * @integrand: a #NcXcorKernelIntegrand
 * @k: wavenumber
 *
 * Evaluates every row at wavenumber @k, see nc_xcor_kernel_integrand_eval().
 *
 * Returns: (transfer full) (element-type gdouble): a #GArray of nc_xcor_kernel_integrand_get_len() values
 */
GArray *
nc_xcor_kernel_integrand_eval_array (NcXcorKernelIntegrand *integrand, gdouble k)
{
  GArray *arr = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), integrand->len);

  g_array_set_size (arr, integrand->len);
  nc_xcor_kernel_integrand_eval (integrand, k, &g_array_index (arr, gdouble, 0));

  return arr;
}

