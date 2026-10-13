/***************************************************************************
 *            test_nc_xcor_closure_check.c
 *
 *  Sun Oct 11 2026
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

/*
 * A Chebyshev closure checked against W_ell(k) computed directly.
 *
 * The direct value shares nothing with the closure but the components themselves:
 *
 *   W_ell(k) = sum_c prefactor_c(k, ell) int kernel_c(chi, k) j_ell^(d_c)(k chi) dchi,
 *
 * with its own Levin integrator over each component's chi limits, as the kernel's
 * non-Limber evaluator defines it, but without the evaluator, its truncation or any
 * fit. Points are taken inside every panel, between its Chebyshev-Lobatto nodes and
 * away from its edges, so each lies on a definite side of every cut, and on a
 * log-spaced grid over the closure's range, which a missed feature cannot avoid.
 *
 * The bound at a point is what the closure was asked for: the acceptance test of its
 * panel, max(reltol vscale, abstol) with vscale the largest |W| on the panel's nodes and
 * abstol peak-epsilon times the smallest multipole peak, plus adaptive-epsilon times
 * the block's norm, below which a component's tail is truncated to zero.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>

#include "test_nc_xcor_closure_check.h"

#define TEST_CLOSURE_CHECK_GRID 200

typedef struct _TestComponentArgs
{
  NcXcorKernelComponent *comp;
  NcHICosmo *cosmo;
} TestComponentArgs;

static gdouble
_test_component_kernel (gpointer user_data, gdouble chi, gdouble k)
{
  TestComponentArgs *args = user_data;

  return nc_xcor_kernel_component_eval_kernel (args->comp, args->cosmo, chi, k);
}

/* W_ell(k) of a non-Limber kernel from its components, with sbi set to the range [ell, ell] */
gdouble
test_nc_xcor_closure_direct (NcXcorKernel *xclk, NcHICosmo *cosmo, NcmSBesselIntegrator *sbi, gint ell, gdouble k)
{
  GPtrArray *comps = nc_xcor_kernel_get_component_list (xclk);
  NcmVector *out   = ncm_vector_new (1);
  gdouble total    = 0.0;
  guint c;

  for (c = 0; c < comps->len; c++)
  {
    NcXcorKernelComponent *comp = g_ptr_array_index (comps, c);
    TestComponentArgs args      = {comp, cosmo};
    gdouble chi_min, chi_max, k_min, k_max;

    nc_xcor_kernel_component_get_limits (comp, cosmo, &chi_min, &chi_max, &k_min, &k_max);
    ncm_sbessel_integrator_integrate_deriv (sbi, _test_component_kernel, chi_min, chi_max, k,
                                            nc_xcor_kernel_component_get_bessel_deriv (comp), out, &args);
    total += nc_xcor_kernel_component_eval_prefactor (comp, cosmo, k, ell) * ncm_vector_get (out, 0);
  }

  ncm_vector_free (out);
  g_ptr_array_unref (comps);

  return total;
}

static void _test_closure_points (NcmSpectralPanels *sp, GArray *points, GArray *panel_of);
static void _test_closure_vscale (NcXcorKernelIntegrand *closure, NcmSpectralPanels *sp, guint n_rows, NcmMatrix *vscale);

/* Compares a non-Limber Chebyshev closure of the multipoles lmin, lmin + 1, ... with the
 * direct values, and returns the worst deviation over its bound */
TestNcXcorClosureCheck
test_nc_xcor_closure_check (NcXcorKernel *xclk, NcHICosmo *cosmo, NcXcorKernelIntegrand *closure, gint lmin)
{
  NcmSpectralPanels *sp        = nc_xcor_kernel_integrand_peek_panels (closure);
  const guint n_rows           = ncm_spectral_panels_get_fdim (sp);
  const guint n_panels         = ncm_spectral_panels_get_n_panels (sp);
  GArray *points               = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *panel_of             = g_array_new (FALSE, FALSE, sizeof (guint));
  NcmMatrix *vscale            = ncm_matrix_new (n_panels, n_rows);
  NcmMatrix *got               = NULL;
  NcmMatrix *direct            = NULL;
  NcmVector *peaks             = ncm_vector_new (n_rows);
  TestNcXcorClosureCheck check = {0.0, lmin, 0.0, 0};
  gdouble abstol, truncation, peak_min = GSL_POSINF, norm2 = 0.0;
  guint p, r;

  _test_closure_points (sp, points, panel_of);
  _test_closure_vscale (closure, sp, n_rows, vscale);

  got    = ncm_matrix_new (points->len, n_rows);
  direct = ncm_matrix_new (points->len, n_rows);
  ncm_vector_set_zero (peaks);

  for (r = 0; r < n_rows; r++)
  {
    NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (lmin + r, lmin + r));

    for (p = 0; p < points->len; p++)
    {
      const gdouble d = test_nc_xcor_closure_direct (xclk, cosmo, sbi, lmin + r, g_array_index (points, gdouble, p));

      ncm_matrix_set (direct, p, r, d);
      ncm_vector_set (peaks, r, MAX (ncm_vector_get (peaks, r), fabs (d)));
    }

    ncm_sbessel_integrator_free (sbi);
  }

  for (p = 0; p < points->len; p++)
    nc_xcor_kernel_integrand_eval (closure, g_array_index (points, gdouble, p), ncm_matrix_ptr (got, p, 0));

  for (r = 0; r < n_rows; r++)
  {
    peak_min = MIN (peak_min, ncm_vector_get (peaks, r));
    norm2   += gsl_pow_2 (ncm_vector_get (peaks, r));
  }

  abstol     = nc_xcor_kernel_get_peak_epsilon (xclk) * peak_min;
  truncation = nc_xcor_kernel_get_adaptive_epsilon (xclk) * sqrt (norm2);

  for (p = 0; p < points->len; p++)
  {
    const guint i = g_array_index (panel_of, guint, p);

    for (r = 0; r < n_rows; r++)
    {
      const gdouble bound = MAX (nc_xcor_kernel_get_reltol (xclk) * ncm_matrix_get (vscale, i, r), abstol) + truncation;
      const gdouble ratio = fabs (ncm_matrix_get (got, p, r) - ncm_matrix_get (direct, p, r)) / bound;

      if (ratio > check.worst)
      {
        check.worst  = ratio;
        check.at_ell = lmin + r;
        check.at_k   = g_array_index (points, gdouble, p);
      }
    }
  }

  check.n_points = points->len;

  ncm_vector_free (peaks);
  ncm_matrix_free (direct);
  ncm_matrix_free (got);
  ncm_matrix_free (vscale);
  g_array_unref (panel_of);
  g_array_unref (points);

  return check;
}

/* Three points of every panel, halfway in angle between Chebyshev-Lobatto nodes, next to
 * each edge and in the middle, and a log-spaced grid over the range, with the panel of each */
static void
_test_closure_points (NcmSpectralPanels *sp, GArray *points, GArray *panel_of)
{
  const guint n_panels = ncm_spectral_panels_get_n_panels (sp);
  gdouble lo, hi;
  guint i, j;

  for (i = 0; i < n_panels; i++)
  {
    gdouble a, b;
    guint N, m, js[3];

    ncm_spectral_panels_peek_panel (sp, i, &a, &b, &N);

    m     = N - 1;
    js[0] = 0;
    js[1] = m / 2;
    js[2] = m - 1;

    for (j = 0; j < 3; j++)
    {
      const gdouble x = 0.5 * (a + b) + 0.5 * (b - a) * cos ((js[j] + 0.5) * M_PI / m);

      g_array_append_val (points, x);
      g_array_append_val (panel_of, i);
    }
  }

  ncm_spectral_panels_get_range (sp, &lo, &hi);

  for (j = 1; j <= TEST_CLOSURE_CHECK_GRID; j++)
  {
    const gdouble x = lo * pow (hi / lo, j / (TEST_CLOSURE_CHECK_GRID + 1.0));
    const guint ip  = ncm_spectral_panels_find_panel (sp, x);

    g_array_append_val (points, x);
    g_array_append_val (panel_of, ip);
  }
}

/* The largest |W| of each row on each panel's Chebyshev-Lobatto nodes */
static void
_test_closure_vscale (NcXcorKernelIntegrand *closure, NcmSpectralPanels *sp, guint n_rows, NcmMatrix *vscale)
{
  const guint n_panels = ncm_spectral_panels_get_n_panels (sp);
  GArray *W            = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), n_rows);
  guint i, j, r;

  g_array_set_size (W, n_rows);
  ncm_matrix_set_zero (vscale);

  for (i = 0; i < n_panels; i++)
  {
    gdouble a, b;
    guint N;

    ncm_spectral_panels_peek_panel (sp, i, &a, &b, &N);

    for (j = 0; j < N; j++)
    {
      const gdouble x = 0.5 * (a + b) + 0.5 * (b - a) * cos (j * M_PI / (N - 1));

      nc_xcor_kernel_integrand_eval (closure, x, &g_array_index (W, gdouble, 0));

      for (r = 0; r < n_rows; r++)
        ncm_matrix_set (vscale, i, r, MAX (ncm_matrix_get (vscale, i, r), fabs (g_array_index (W, gdouble, r))));
    }
  }

  g_array_unref (W);
}

