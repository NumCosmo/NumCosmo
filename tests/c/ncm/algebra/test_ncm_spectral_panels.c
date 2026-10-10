/***************************************************************************
 *            test_ncm_spectral_panels.c
 *
 *  Sat Oct 10 2026
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <string.h>
#include <glib.h>
#include <glib-object.h>

/* Three panels with a different number of coefficients each */
#define TEST_N_PANELS 3
static const gdouble test_edges[TEST_N_PANELS + 1] = {0.0, 1.0, 2.5, 4.0};
static const guint test_N[TEST_N_PANELS]           = {4, 9, 5};

/* Row 0 is x^3, exact on every panel since each has at least four coefficients. Row 1 is
 * i + 1 on panel i, so it jumps at every interior edge. */
static gdouble
_cubic (gpointer user_data, gdouble x)
{
  return x * x * x;
}

static gdouble
_step (gpointer user_data, gdouble x)
{
  return 1.0 + GPOINTER_TO_UINT (user_data);
}

static gdouble
_exact (guint r, guint i, gdouble x)
{
  return (r == 0) ? x * x * x : 1.0 + i;
}

/* The panel set above, the coefficients computed on the nodes of each panel */
static NcmSpectralPanels *
_panels_new (void)
{
  NcmSpectralPanels *sp = ncm_spectral_panels_new (2);
  NcmSpectral *spectral = ncm_spectral_new ();
  GArray *row           = NULL;
  guint i;

  for (i = 0; i < TEST_N_PANELS; i++)
  {
    const gdouble a   = test_edges[i];
    const gdouble b   = test_edges[i + 1];
    NcmMatrix *coeffs = ncm_matrix_new (2, test_N[i]);

    ncm_spectral_compute_chebyshev_coeffs (spectral, _cubic, a, b, test_N[i], &row, NULL);
    memcpy (ncm_matrix_ptr (coeffs, 0, 0), row->data, sizeof (gdouble) * test_N[i]);
    ncm_spectral_compute_chebyshev_coeffs (spectral, _step, a, b, test_N[i], &row, GUINT_TO_POINTER (i));
    memcpy (ncm_matrix_ptr (coeffs, 1, 0), row->data, sizeof (gdouble) * test_N[i]);

    ncm_spectral_panels_append (sp, a, b, coeffs);
    ncm_matrix_free (coeffs);
  }

  g_array_unref (row);
  ncm_spectral_free (spectral);

  return sp;
}

static void
test_ncm_spectral_panels_new (void)
{
  NcmSpectralPanels *sp = ncm_spectral_panels_new (3);

  g_assert_cmpuint (ncm_spectral_panels_get_fdim (sp), ==, 3);
  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (sp), ==, 0);

  ncm_spectral_panels_free (sp);
}

/* The edges, sizes and coefficients read back are the ones appended, from a matrix whose
 * rows are not contiguous too */
static void
test_ncm_spectral_panels_append (void)
{
  NcmSpectralPanels *sp = _panels_new ();
  NcmMatrix *wide       = ncm_matrix_new (2, 7);
  NcmMatrix *view       = ncm_matrix_get_submatrix (wide, 0, 0, 2, 3);
  gdouble a, b;
  guint i, N;

  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (sp), ==, TEST_N_PANELS);
  ncm_spectral_panels_get_range (sp, &a, &b);
  g_assert_cmpfloat (a, ==, test_edges[0]);
  g_assert_cmpfloat (b, ==, test_edges[TEST_N_PANELS]);

  for (i = 0; i < TEST_N_PANELS; i++)
  {
    const gdouble *c  = ncm_spectral_panels_peek_panel (sp, i, &a, &b, &N);
    NcmMatrix *coeffs = ncm_spectral_panels_get_panel (sp, i, &a, &b);
    guint r, n;

    g_assert_cmpfloat (a, ==, test_edges[i]);
    g_assert_cmpfloat (b, ==, test_edges[i + 1]);
    g_assert_cmpuint (N, ==, test_N[i]);
    g_assert_cmpuint (ncm_matrix_nrows (coeffs), ==, 2);
    g_assert_cmpuint (ncm_matrix_ncols (coeffs), ==, N);

    for (r = 0; r < 2; r++)
      for (n = 0; n < N; n++)
        g_assert_cmpfloat (ncm_matrix_get (coeffs, r, n), ==, c[r * N + n]);

    ncm_matrix_free (coeffs);
  }

  for (i = 0; i < 2 * 7; i++)
    ncm_matrix_data (wide)[i] = i;

  ncm_spectral_panels_append (sp, test_edges[TEST_N_PANELS], 5.0, view);

  {
    const gdouble *c = ncm_spectral_panels_peek_panel (sp, TEST_N_PANELS, &a, &b, &N);

    g_assert_cmpuint (N, ==, 3);
    g_assert_cmpfloat (c[0], ==, 0.0);
    g_assert_cmpfloat (c[2], ==, 2.0);
    g_assert_cmpfloat (c[3], ==, 7.0);
    g_assert_cmpfloat (c[5], ==, 9.0);
  }

  ncm_matrix_free (view);
  ncm_matrix_free (wide);
  ncm_spectral_panels_free (sp);
}

/* An interior edge belongs to the panel that starts there; outside the range, the end panels */
static void
test_ncm_spectral_panels_find_panel (void)
{
  NcmSpectralPanels *sp = _panels_new ();

  g_assert_cmpuint (ncm_spectral_panels_find_panel (sp, -1.0), ==, 0);
  g_assert_cmpuint (ncm_spectral_panels_find_panel (sp, 0.0), ==, 0);
  g_assert_cmpuint (ncm_spectral_panels_find_panel (sp, 0.5), ==, 0);
  g_assert_cmpuint (ncm_spectral_panels_find_panel (sp, 1.0), ==, 1);
  g_assert_cmpuint (ncm_spectral_panels_find_panel (sp, 2.0), ==, 1);
  g_assert_cmpuint (ncm_spectral_panels_find_panel (sp, 2.5), ==, 2);
  g_assert_cmpuint (ncm_spectral_panels_find_panel (sp, 4.0), ==, 2);
  g_assert_cmpuint (ncm_spectral_panels_find_panel (sp, 9.0), ==, 2);

  ncm_spectral_panels_free (sp);
}

/* Both rows match the exact functions inside the panels, and on an interior edge row 1
 * takes the value of the panel to its right */
static void
test_ncm_spectral_panels_eval (void)
{
  NcmSpectralPanels *sp = _panels_new ();
  NcmVector *fval       = ncm_vector_new (2);
  guint j;

  for (j = 0; j <= 400; j++)
  {
    const gdouble x = test_edges[TEST_N_PANELS] * j / 400.0;
    const guint i   = ncm_spectral_panels_find_panel (sp, x);
    guint r;

    ncm_spectral_panels_eval (sp, x, fval);

    for (r = 0; r < 2; r++)
    {
      const gdouble exact = _exact (r, i, x);

      g_assert_cmpfloat (fabs (ncm_vector_get (fval, r) - exact), <=, 1.0e-14 * MAX (1.0, fabs (exact)));
      g_assert_cmpfloat (ncm_spectral_panels_eval_row (sp, r, x), ==, ncm_vector_get (fval, r));
    }
  }

  g_assert_cmpfloat (fabs (ncm_spectral_panels_eval_row (sp, 1, 1.0) - 2.0), <=, 1.0e-14);
  g_assert_cmpfloat (fabs (ncm_spectral_panels_eval_row (sp, 1, 2.5) - 3.0), <=, 1.0e-14);

  ncm_vector_free (fval);
  ncm_spectral_panels_free (sp);
}

/* A subinterval of a panel gives the same polynomial; the whole panel gives its own
 * coefficients; an interval across an edge is refused */
static void
test_ncm_spectral_panels_restrict (void)
{
  NcmSpectralPanels *sp = _panels_new ();
  NcmSpectral *spectral = ncm_spectral_new ();
  NcmMatrix *coeffs     = NULL;
  NcmMatrix *kept;
  GArray *row = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const gdouble a = 1.2, b = 2.1;
  gdouble pa, pb;
  guint N, r, n, j;
  const gdouble *c;

  g_assert_true (ncm_spectral_panels_restrict (sp, spectral, a, b, &coeffs));
  g_assert_cmpuint (ncm_matrix_nrows (coeffs), ==, 2);
  g_assert_cmpuint (ncm_matrix_ncols (coeffs), ==, test_N[1]);

  g_array_set_size (row, test_N[1]);

  for (r = 0; r < 2; r++)
  {
    memcpy (row->data, ncm_matrix_ptr (coeffs, r, 0), sizeof (gdouble) * test_N[1]);

    for (j = 0; j <= 50; j++)
    {
      const gdouble x = a + (b - a) * j / 50.0;

      g_assert_cmpfloat (fabs (ncm_spectral_chebyshev_eval_x (row, a, b, x) - ncm_spectral_panels_eval_row (sp, r, x)), <=, 1.0e-13);
    }
  }

  /* The reused matrix has the shape of panel 1; panel 2 needs another one */
  kept = coeffs;
  g_assert_true (ncm_spectral_panels_restrict (sp, spectral, 1.0, 2.5, &coeffs));
  g_assert_true (coeffs == kept);

  c = ncm_spectral_panels_peek_panel (sp, 1, &pa, &pb, &N);

  for (r = 0; r < 2; r++)
    for (n = 0; n < N; n++)
      g_assert_cmpfloat (ncm_matrix_get (coeffs, r, n), ==, c[r * N + n]);

  g_assert_true (ncm_spectral_panels_restrict (sp, spectral, 3.0, 3.5, &coeffs));
  g_assert_cmpuint (ncm_matrix_ncols (coeffs), ==, test_N[2]);

  kept = coeffs;
  g_assert_false (ncm_spectral_panels_restrict (sp, spectral, 2.0, 3.0, &coeffs));
  g_assert_true (coeffs == kept);

  g_array_unref (row);
  ncm_matrix_free (coeffs);
  ncm_spectral_free (spectral);
  ncm_spectral_panels_free (sp);
}

/* After a reset the panel set holds new panels from any first edge */
static void
test_ncm_spectral_panels_reset (void)
{
  NcmSpectralPanels *sp = _panels_new ();
  NcmMatrix *coeffs     = ncm_matrix_new (2, 2);
  gdouble a, b;

  ncm_matrix_set_all (coeffs, 1.0);
  ncm_spectral_panels_reset (sp);
  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (sp), ==, 0);
  g_assert_cmpuint (ncm_spectral_panels_get_fdim (sp), ==, 2);

  ncm_spectral_panels_append (sp, -3.0, -2.0, coeffs);
  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (sp), ==, 1);
  ncm_spectral_panels_get_range (sp, &a, &b);
  g_assert_cmpfloat (a, ==, -3.0);
  g_assert_cmpfloat (b, ==, -2.0);
  g_assert_cmpfloat (ncm_spectral_panels_eval_row (sp, 0, -2.0), ==, 2.0);

  ncm_matrix_free (coeffs);
  ncm_spectral_panels_free (sp);
}

/* A copy through serialization has the same edges, sizes and coefficients, bit for bit */
static void
test_ncm_spectral_panels_serialize (void)
{
  NcmSpectralPanels *sp    = _panels_new ();
  NcmSpectralPanels *empty = ncm_spectral_panels_new (4);
  NcmSpectralPanels *dup   = NCM_SPECTRAL_PANELS (ncm_serialize_global_dup_obj (G_OBJECT (sp)));
  NcmSpectralPanels *edup  = NCM_SPECTRAL_PANELS (ncm_serialize_global_dup_obj (G_OBJECT (empty)));
  guint i;

  g_assert_cmpuint (ncm_spectral_panels_get_fdim (dup), ==, 2);
  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (dup), ==, TEST_N_PANELS);

  for (i = 0; i < TEST_N_PANELS; i++)
  {
    gdouble a, b, da, db;
    guint N, dN, n;
    const gdouble *c  = ncm_spectral_panels_peek_panel (sp, i, &a, &b, &N);
    const gdouble *dc = ncm_spectral_panels_peek_panel (dup, i, &da, &db, &dN);

    g_assert_cmpfloat (da, ==, a);
    g_assert_cmpfloat (db, ==, b);
    g_assert_cmpuint (dN, ==, N);

    for (n = 0; n < 2 * N; n++)
      g_assert_cmpfloat (dc[n], ==, c[n]);
  }

  g_assert_cmpuint (ncm_spectral_panels_get_fdim (edup), ==, 4);
  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (edup), ==, 0);

  ncm_spectral_panels_free (edup);
  ncm_spectral_panels_free (dup);
  ncm_spectral_panels_free (empty);
  ncm_spectral_panels_free (sp);
}

static void
test_ncm_spectral_panels_append_rows (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels/append/rows/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*expected fdim = 2*");
}

static void
test_ncm_spectral_panels_append_rows_subprocess (void)
{
  NcmSpectralPanels *sp = ncm_spectral_panels_new (2);
  NcmMatrix *coeffs     = ncm_matrix_new (3, 4);

  ncm_spectral_panels_append (sp, 0.0, 1.0, coeffs);
}

static void
test_ncm_spectral_panels_append_empty (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels/append/empty/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*is empty*");
}

static void
test_ncm_spectral_panels_append_empty_subprocess (void)
{
  NcmSpectralPanels *sp = ncm_spectral_panels_new (2);
  NcmMatrix *coeffs     = ncm_matrix_new (2, 4);

  ncm_spectral_panels_append (sp, 1.0, 1.0, coeffs);
}

static void
test_ncm_spectral_panels_append_gap (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels/append/gap/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*not at the last edge*");
}

static void
test_ncm_spectral_panels_append_gap_subprocess (void)
{
  NcmSpectralPanels *sp = _panels_new ();
  NcmMatrix *coeffs     = ncm_matrix_new (2, 4);

  ncm_spectral_panels_append (sp, 4.5, 5.0, coeffs);
}

static void
test_ncm_spectral_panels_properties_edges (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels/properties/edges/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*3 edges for 1 panels*");
}

static void
test_ncm_spectral_panels_properties_edges_subprocess (void)
{
  GVariant *panels = g_variant_ref_sink (g_variant_new_parsed ("([0.0, 1.0, 2.0], [[1.0, 2.0]])"));

  g_object_new (NCM_TYPE_SPECTRAL_PANELS, "fdim", 1, "panels", panels, NULL);
}

static void
test_ncm_spectral_panels_properties_rows (void)
{
  g_test_trap_subprocess ("/ncm/spectral_panels/properties/rows/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*not a multiple of fdim = 2*");
}

static void
test_ncm_spectral_panels_properties_rows_subprocess (void)
{
  GVariant *panels = g_variant_ref_sink (g_variant_new_parsed ("([0.0, 1.0], [[1.0, 2.0, 3.0]])"));

  g_object_new (NCM_TYPE_SPECTRAL_PANELS, "fdim", 2, "panels", panels, NULL);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/spectral_panels/new", &test_ncm_spectral_panels_new);
  g_test_add_func ("/ncm/spectral_panels/append", &test_ncm_spectral_panels_append);
  g_test_add_func ("/ncm/spectral_panels/find_panel", &test_ncm_spectral_panels_find_panel);
  g_test_add_func ("/ncm/spectral_panels/eval", &test_ncm_spectral_panels_eval);
  g_test_add_func ("/ncm/spectral_panels/restrict", &test_ncm_spectral_panels_restrict);
  g_test_add_func ("/ncm/spectral_panels/reset", &test_ncm_spectral_panels_reset);
  g_test_add_func ("/ncm/spectral_panels/serialize", &test_ncm_spectral_panels_serialize);
  g_test_add_func ("/ncm/spectral_panels/append/rows", &test_ncm_spectral_panels_append_rows);
  g_test_add_func ("/ncm/spectral_panels/append/rows/subprocess", &test_ncm_spectral_panels_append_rows_subprocess);
  g_test_add_func ("/ncm/spectral_panels/append/empty", &test_ncm_spectral_panels_append_empty);
  g_test_add_func ("/ncm/spectral_panels/append/empty/subprocess", &test_ncm_spectral_panels_append_empty_subprocess);
  g_test_add_func ("/ncm/spectral_panels/append/gap", &test_ncm_spectral_panels_append_gap);
  g_test_add_func ("/ncm/spectral_panels/append/gap/subprocess", &test_ncm_spectral_panels_append_gap_subprocess);
  g_test_add_func ("/ncm/spectral_panels/properties/edges", &test_ncm_spectral_panels_properties_edges);
  g_test_add_func ("/ncm/spectral_panels/properties/edges/subprocess", &test_ncm_spectral_panels_properties_edges_subprocess);
  g_test_add_func ("/ncm/spectral_panels/properties/rows", &test_ncm_spectral_panels_properties_rows);
  g_test_add_func ("/ncm/spectral_panels/properties/rows/subprocess", &test_ncm_spectral_panels_properties_rows_subprocess);

  g_test_run ();

  return 0;
}

