/***************************************************************************
 *            test_ncm_laurent_series.c
 *
 *  Tue Jul 8 2026
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

/* Independent tests for NcmLaurentSeries, deliberately with no reference to
 * the weak-lensing problem this was extracted from -- pure complex
 * Laurent-polynomial arithmetic and the Jacobi-Anger reduction to scaled
 * Bessel functions, cross-checked here against direct numerical
 * theta-integration (trustworthy for correctness) and against hand-computed
 * small cases. Also exercises the introspectable (NcmComplex-based) API
 * directly, confirming it agrees with the native (`_c`-suffixed) one. */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <complex.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

static void
test_add_scale_hand_computed (void)
{
  NcmLaurentSeries *a   = ncm_laurent_series_new_single (1, 2.0 + 1.0 * I);
  NcmLaurentSeries *b   = ncm_laurent_series_new_single (-1, 3.0 - 2.0 * I);
  NcmLaurentSeries *sum = ncm_laurent_series_add (a, b, 2.0);
  NcmLaurentSeries *scl = ncm_laurent_series_scale_c (a, I);

  /* a + 2*b: harmonic 1 -> (2+i), harmonic -1 -> 2*(3-2i) = (6-4i) */
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (sum, 1) - (2.0 + 1.0 * I)), <, 1.0e-13);
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (sum, -1) - (6.0 - 4.0 * I)), <, 1.0e-13);
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (sum, 0)), <, 1.0e-13);

  /* i*(2+i) = -1+2i */
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (scl, 1) - (-1.0 + 2.0 * I)), <, 1.0e-13);

  ncm_laurent_series_free (a);
  ncm_laurent_series_free (b);
  ncm_laurent_series_free (sum);
  ncm_laurent_series_free (scl);
}

static void
test_conv_hand_computed (void)
{
  /* (2 + 3w) * (1 - w^{-1}) = 2 - 2w^{-1} + 3w - 3 = -1 - 2w^{-1} + 3w */
  NcmLaurentSeries *a = ncm_laurent_series_new (0, 1);
  NcmLaurentSeries *b = ncm_laurent_series_new (-1, 0);
  NcmLaurentSeries *c;

  ncm_laurent_series_set_c (a, 0, 2.0);
  ncm_laurent_series_set_c (a, 1, 3.0);
  ncm_laurent_series_set_c (b, -1, -1.0);
  ncm_laurent_series_set_c (b, 0, 1.0);

  c = ncm_laurent_series_conv (a, b);

  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (c, -1) - (-2.0)), <, 1.0e-13);
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (c, 0) - (-1.0)), <, 1.0e-13);
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (c, 1) - (3.0)), <, 1.0e-13);

  ncm_laurent_series_free (a);
  ncm_laurent_series_free (b);
  ncm_laurent_series_free (c);
}

static void
test_conj_hand_computed (void)
{
  NcmLaurentSeries *a = ncm_laurent_series_new (-1, 2);
  NcmLaurentSeries *ab;

  ncm_laurent_series_set_c (a, -1, 1.0 + 2.0 * I);
  ncm_laurent_series_set_c (a, 0, 3.0);
  ncm_laurent_series_set_c (a, 1, -1.0 + I);
  ncm_laurent_series_set_c (a, 2, 0.5);

  ab = ncm_laurent_series_conj (a);

  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (ab, 1) - conj (1.0 + 2.0 * I)), <, 1.0e-13);
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (ab, 0) - conj (3.0)), <, 1.0e-13);
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (ab, -1) - conj (-1.0 + I)), <, 1.0e-13);
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (ab, -2) - conj (0.5)), <, 1.0e-13);

  ncm_laurent_series_free (a);
  ncm_laurent_series_free (ab);
}

static void
test_eval_hand_computed (void)
{
  /* a(w) = (1+2i) + 3i*w^{-1} - 0.5*w^2, evaluated at w=2-i. */
  NcmLaurentSeries *a    = ncm_laurent_series_new (-1, 2);
  const complex double w = 2.0 - I;
  complex double expected;

  ncm_laurent_series_set_c (a, -1, 3.0 * I);
  ncm_laurent_series_set_c (a, 0, 1.0 + 2.0 * I);
  ncm_laurent_series_set_c (a, 1, 0.0);
  ncm_laurent_series_set_c (a, 2, -0.5);

  expected = (1.0 + 2.0 * I) + 3.0 * I / w - 0.5 * w * w;

  g_assert_cmpfloat (cabs (ncm_laurent_series_eval_c (a, w) - expected), <, 1.0e-12);

  ncm_laurent_series_free (a);
}

/* Ground truth: direct numerical theta-integration of
 * cm(theta)*exp(z*(cos(theta-phi)-1)), no Jacobi-Anger identity used. */
struct _ThetaCtx
{
  const NcmLaurentSeries *cm;
  gdouble phi;
  gdouble z;
};

static gdouble
_theta_integrand (gdouble theta, gpointer p)
{
  struct _ThetaCtx *c = (struct _ThetaCtx *) p;
  complex double v    = 0.0;
  gint h;

  for (h = ncm_laurent_series_hmin (c->cm); h <= ncm_laurent_series_hmax (c->cm); h++)
    v += ncm_laurent_series_get_c (c->cm, h) * cexp (I * h * theta);

  return creal (v) * exp (c->z * (cos (theta - c->phi) - 1.0));
}

static gdouble
_direct_theta_integral (const NcmLaurentSeries *cm, gdouble phi, gdouble z)
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  gdouble result, abserr;
  struct _ThetaCtx ctx = {cm, phi, z};

  F.function = &_theta_integrand;
  F.params   = &ctx;

  gsl_integration_qag (&F, 0.0, 2.0 * M_PI, 1.0e-12, 1.0e-11, 1000, GSL_INTEG_GAUSS61, w, &result, &abserr);
  gsl_integration_workspace_free (w);

  return result;
}

static void
test_jacobi_anger_matches_direct_integration (void)
{
  GRand *rng = g_rand_new_with_seed (20260708);
  guint trial;

  for (trial = 0; trial < 20; trial++)
  {
    const gint maxh      = 4;
    NcmLaurentSeries *cm = ncm_laurent_series_new (-maxh, maxh);
    const gdouble phi    = g_rand_double_range (rng, -M_PI, M_PI);
    const gdouble z      = g_rand_double_range (rng, 0.1, 50.0);
    gdouble Ik[5];
    gdouble direct, reduced;
    gint k;

    for (k = -maxh; k <= maxh; k++)
      ncm_laurent_series_set_c (cm, k, g_rand_double_range (rng, -1.0, 1.0) + I * g_rand_double_range (rng, -1.0, 1.0));

    /* cm must be real-valued as a function of theta for the reduction's own
     * definition (2*Re(...)) to equal a real integral -- enforce Hermitian
     * symmetry c_{-h} = conj(c_h), matching how the physics code's own
     * series are always built (real observable quantities). */
    for (k = 1; k <= maxh; k++)
      ncm_laurent_series_set_c (cm, -k, conj (ncm_laurent_series_get_c (cm, k)));

    ncm_laurent_series_set_c (cm, 0, creal (ncm_laurent_series_get_c (cm, 0)));

    for (k = 0; k <= maxh; k++)
      Ik[k] = gsl_sf_bessel_In_scaled (k, z);

    direct  = _direct_theta_integral (cm, phi, z);
    reduced = ncm_laurent_series_jacobi_anger_reduce (cm, phi, Ik, maxh + 1);

    g_assert_cmpfloat (fabs (reduced - direct), <, 1.0e-8 * (1.0 + fabs (direct)));

    ncm_laurent_series_free (cm);
  }

  g_rand_free (rng);
}

static NcmLaurentSeries *
_track (GPtrArray *scratch, NcmLaurentSeries *x)
{
  g_ptr_array_add (scratch, x);

  return x;
}

/* Cross-language regression: chi_I(chi_L,g) Taylor coefficients, chi (TRACE)
 * convention, matching dev-notes/wl_shape_series_marginalization_derivation.py
 * verify_chi_closed_form_pieces() at rho=0.25, theta=0.7, g=0.09, N=12:
 * chi_I ~ 0.0131042280+0.1640678007j (see that function's own docstring).
 * Recursion (module docstring's own derivation): c_0 = chi_L;
 * c_k = n_k - d1*c_{k-1} - d2*c_{k-2} (k>=1), n_1=-2, n_2=chibar_L, n_k=0
 * otherwise, d1=-(chi_L+chibar_L), d2=1. */
static void
test_chi_taylor_matches_python_reference (void)
{
  const gdouble rho             = 0.25;
  const gdouble theta           = 0.7;
  const gdouble g               = 0.09;
  const gint N                  = 12;
  const complex double chi_L    = rho * cexp (I * theta);
  const complex double chibar_L = conj (chi_L);
  GPtrArray *scratch            = g_ptr_array_new_with_free_func ((GDestroyNotify) ncm_laurent_series_free);

  /* Laurent coefficients are theta-independent (chi_L=rho*w, so the plain
   * real rho is w^1's coefficient; theta-dependence comes entirely from the
   * harmonic exponential at evaluation time below) -- matches the Python
   * reference's own chi_I_taylor_chi(rho, N) signature, which never takes
   * theta at all. */
  NcmLaurentSeries *chi       = _track (scratch, ncm_laurent_series_new_single (1, rho));
  NcmLaurentSeries *chibar    = _track (scratch, ncm_laurent_series_new_single (-1, rho));
  NcmLaurentSeries *two_s     = _track (scratch, ncm_laurent_series_add (chi, chibar, 1.0));
  NcmLaurentSeries *d1        = _track (scratch, ncm_laurent_series_scale_c (two_s, -1.0));
  NcmLaurentSeries *n1        = _track (scratch, ncm_laurent_series_new_single (0, -2.0));
  NcmLaurentSeries **c        = g_new0 (NcmLaurentSeries *, N + 1);
  complex double chi_I_series = 0.0;
  complex double chi_I_exact;
  gint k;

  c[0] = chi;

  for (k = 1; k <= N; k++)
  {
    NcmLaurentSeries *base;
    NcmLaurentSeries *acc;

    if (k == 1)
      base = n1;
    else if (k == 2)
      base = chibar;
    else
      base = _track (scratch, ncm_laurent_series_new (0, 0));

    acc = _track (scratch, ncm_laurent_series_add (base, _track (scratch, ncm_laurent_series_conv (d1, c[k - 1])), -1.0));

    if (k >= 2)
      acc = _track (scratch, ncm_laurent_series_add (acc, c[k - 2], -1.0));

    c[k] = acc;
  }

  for (k = 0; k <= N; k++)
  {
    complex double term = 0.0;
    gint h;

    for (h = ncm_laurent_series_hmin (c[k]); h <= ncm_laurent_series_hmax (c[k]); h++)
      term += ncm_laurent_series_get_c (c[k], h) * cexp (I * h * theta);

    chi_I_series += term * pow (g, k);
  }

  chi_I_exact = (chi_L - 2.0 * g + g * g * chibar_L) / (1.0 - g * (chi_L + chibar_L) + g * g);

  g_assert_cmpfloat (cabs (chi_I_series - chi_I_exact), <, 1.0e-10);
  g_assert_cmpfloat (fabs (creal (chi_I_series) - 0.0131042280), <, 1.0e-6);
  g_assert_cmpfloat (fabs (cimag (chi_I_series) - 0.1640678007), <, 1.0e-6);

  g_free (c);
  g_ptr_array_unref (scratch);
}

static void
test_introspectable_api_matches_native (void)
{
  NcmComplex *va  = ncm_complex_new ();
  NcmComplex *vb  = ncm_complex_new ();
  NcmComplex *s   = ncm_complex_new ();
  NcmComplex *out = ncm_complex_new ();
  NcmLaurentSeries *a, *b, *sum, *conv, *cj, *scl;

  g_assert_true (ncm_laurent_series_get_type () != G_TYPE_INVALID);

  ncm_complex_set (va, 2.0, 1.0);
  ncm_complex_set (vb, 3.0, -2.0);
  ncm_complex_set (s, 0.0, 1.0);

  a = ncm_laurent_series_new (1, 1);
  ncm_laurent_series_set (a, 1, va);
  b = ncm_laurent_series_new (-1, -1);
  ncm_laurent_series_set (b, -1, vb);

  g_assert_cmpint (ncm_laurent_series_hmin (a), ==, 1);
  g_assert_cmpint (ncm_laurent_series_hmax (a), ==, 1);

  sum = ncm_laurent_series_add (a, b, 2.0);
  ncm_laurent_series_get (sum, 1, out);
  g_assert_cmpfloat (fabs (ncm_complex_Re (out) - 2.0), <, 1.0e-13);
  g_assert_cmpfloat (fabs (ncm_complex_Im (out) - 1.0), <, 1.0e-13);
  ncm_laurent_series_get (sum, -1, out);
  g_assert_cmpfloat (fabs (ncm_complex_Re (out) - 6.0), <, 1.0e-13);
  g_assert_cmpfloat (fabs (ncm_complex_Im (out) - (-4.0)), <, 1.0e-13);

  conv = ncm_laurent_series_conv (a, b);
  ncm_laurent_series_get (conv, 0, out);
  g_assert_cmpfloat (cabs (ncm_complex_c (out) - (2.0 + I) * (3.0 - 2.0 * I)), <, 1.0e-13);

  cj = ncm_laurent_series_conj (a);
  g_assert_cmpint (ncm_laurent_series_hmin (cj), ==, -1);
  ncm_laurent_series_get (cj, -1, out);
  g_assert_cmpfloat (cabs (ncm_complex_c (out) - conj (2.0 + I)), <, 1.0e-13);

  scl = ncm_laurent_series_scale (a, s);
  ncm_laurent_series_get (scl, 1, out);
  g_assert_cmpfloat (cabs (ncm_complex_c (out) - I * (2.0 + I)), <, 1.0e-13);

  ncm_laurent_series_eval (a, s, out);
  g_assert_cmpfloat (cabs (ncm_complex_c (out) - ncm_laurent_series_eval_c (a, ncm_complex_c (s))), <, 1.0e-13);

  ncm_complex_free (va);
  ncm_complex_free (vb);
  ncm_complex_free (s);
  ncm_complex_free (out);
  ncm_laurent_series_free (a);
  ncm_laurent_series_free (b);
  ncm_laurent_series_free (sum);
  ncm_laurent_series_free (conv);
  ncm_laurent_series_free (cj);
  ncm_laurent_series_free (scl);
}

static void
test_copy_is_independent (void)
{
  NcmLaurentSeries *a = ncm_laurent_series_new_single (2, 5.0 + 0.0 * I);
  NcmLaurentSeries *b = ncm_laurent_series_copy (a);

  ncm_laurent_series_set_c (a, 2, 7.0);

  g_assert_cmpfloat (fabs (creal (ncm_laurent_series_get_c (a, 2)) - 7.0), <, 1.0e-13);
  g_assert_cmpfloat (fabs (creal (ncm_laurent_series_get_c (b, 2)) - 5.0), <, 1.0e-13);

  ncm_laurent_series_free (a);
  ncm_laurent_series_free (b);
}

/* ncm_laurent_series_reset() is the grow-only reuse primitive
 * nc_galaxy_shape_factor_series_lensed.c's per-instance workspace pool
 * relies on: reusing a range within the existing allocation must not
 * reallocate (checked here by comparing the raw `c` pointer, direct field
 * access being explicitly sanctioned by the header's own doc comment on
 * the struct); growing beyond it must reallocate and preserve correctness;
 * either way every coefficient must come back zeroed. */
static void
test_reset_grows_only_when_needed (void)
{
  NcmLaurentSeries *a = ncm_laurent_series_new (0, 4); /* 5 slots */
  complex double *c0;

  ncm_laurent_series_set_c (a, 0, 1.0);
  ncm_laurent_series_set_c (a, 4, 2.0);
  c0 = a->c;

  /* Shrinks the logical range but stays within the 5-slot allocation --
   * must reuse the same buffer and zero it. */
  ncm_laurent_series_reset (a, 0, 2);
  g_assert_true (a->c == c0);
  g_assert_cmpint (a->c_cap, ==, 5);
  g_assert_cmpint (ncm_laurent_series_hmin (a), ==, 0);
  g_assert_cmpint (ncm_laurent_series_hmax (a), ==, 2);
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (a, 0)), <, 1.0e-13);

  ncm_laurent_series_set_c (a, 1, 3.0 + 4.0 * I);

  /* Still within the original 5-slot allocation -- no reallocation. */
  ncm_laurent_series_reset (a, -1, 3);
  g_assert_true (a->c == c0);
  g_assert_cmpint (a->c_cap, ==, 5);
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (a, 1)), <, 1.0e-13);

  ncm_laurent_series_set_c (a, -1, 5.0);

  /* Exceeds the 5-slot allocation -- must grow, preserving no stale data
   * (freshly zeroed) and remaining fully usable. */
  ncm_laurent_series_reset (a, -4, 4);
  g_assert_cmpint (a->c_cap, ==, 9);
  g_assert_cmpint (ncm_laurent_series_hmin (a), ==, -4);
  g_assert_cmpint (ncm_laurent_series_hmax (a), ==, 4);
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (a, -1)), <, 1.0e-13);

  ncm_laurent_series_set_c (a, -4, 9.0);
  g_assert_cmpfloat (fabs (creal (ncm_laurent_series_get_c (a, -4)) - 9.0), <, 1.0e-13);

  ncm_laurent_series_free (a);
}

/* Every `_into` variant must produce results bit-identical to its existing
 * "returns new" sibling, both when the output object already has enough
 * capacity (pure reuse) and when it must grow -- the two paths inside
 * ncm_laurent_series_reset(). */
static void
test_into_variants_match_returns_new (void)
{
  NcmLaurentSeries *a             = ncm_laurent_series_new_single (1, 2.0 + 1.0 * I);
  NcmLaurentSeries *b             = ncm_laurent_series_new_single (-1, 3.0 - 2.0 * I);
  NcmLaurentSeries *expected_sum  = ncm_laurent_series_add (a, b, 2.0);
  NcmLaurentSeries *expected_conv = ncm_laurent_series_conv (a, b);
  NcmLaurentSeries *expected_scl  = ncm_laurent_series_scale_c (a, I);
  NcmLaurentSeries *expected_conj = ncm_laurent_series_conj (a);
  gint h;

  /* Pure reuse: already big enough for every result below. */
  {
    NcmLaurentSeries *out_single = ncm_laurent_series_new (0, 4);
    NcmLaurentSeries *out_sum    = ncm_laurent_series_new (-4, 4);
    NcmLaurentSeries *out_conv   = ncm_laurent_series_new (-4, 4);
    NcmLaurentSeries *out_scl    = ncm_laurent_series_new (-4, 4);
    NcmLaurentSeries *out_conj   = ncm_laurent_series_new (-4, 4);

    ncm_laurent_series_set_single_into (out_single, 1, 2.0 + 1.0 * I);
    g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (out_single, 1) - (2.0 + 1.0 * I)), <, 1.0e-13);

    ncm_laurent_series_add_into (out_sum, a, b, 2.0);
    ncm_laurent_series_conv_into (out_conv, a, b);
    ncm_laurent_series_scale_c_into (out_scl, a, I);
    ncm_laurent_series_conj_into (out_conj, a);

    for (h = -4; h <= 4; h++)
    {
      g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (out_sum, h) - ncm_laurent_series_get_c (expected_sum, h)), <, 1.0e-13);
      g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (out_conv, h) - ncm_laurent_series_get_c (expected_conv, h)), <, 1.0e-13);
      g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (out_scl, h) - ncm_laurent_series_get_c (expected_scl, h)), <, 1.0e-13);
      g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (out_conj, h) - ncm_laurent_series_get_c (expected_conj, h)), <, 1.0e-13);
    }

    ncm_laurent_series_free (out_single);
    ncm_laurent_series_free (out_sum);
    ncm_laurent_series_free (out_conv);
    ncm_laurent_series_free (out_scl);
    ncm_laurent_series_free (out_conj);
  }

  /* Must-grow path: outputs start with zero capacity (a single [0,0] slot). */
  {
    NcmLaurentSeries *out_sum  = ncm_laurent_series_new (0, 0);
    NcmLaurentSeries *out_conv = ncm_laurent_series_new (0, 0);

    ncm_laurent_series_add_into (out_sum, a, b, 2.0);
    ncm_laurent_series_conv_into (out_conv, a, b);

    g_assert_cmpint (ncm_laurent_series_hmin (out_sum), ==, ncm_laurent_series_hmin (expected_sum));
    g_assert_cmpint (ncm_laurent_series_hmax (out_sum), ==, ncm_laurent_series_hmax (expected_sum));

    for (h = ncm_laurent_series_hmin (expected_sum); h <= ncm_laurent_series_hmax (expected_sum); h++)
      g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (out_sum, h) - ncm_laurent_series_get_c (expected_sum, h)), <, 1.0e-13);

    for (h = ncm_laurent_series_hmin (expected_conv); h <= ncm_laurent_series_hmax (expected_conv); h++)
      g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (out_conv, h) - ncm_laurent_series_get_c (expected_conv, h)), <, 1.0e-13);

    ncm_laurent_series_free (out_sum);
    ncm_laurent_series_free (out_conv);
  }

  ncm_laurent_series_free (a);
  ncm_laurent_series_free (b);
  ncm_laurent_series_free (expected_sum);
  ncm_laurent_series_free (expected_conv);
  ncm_laurent_series_free (expected_scl);
  ncm_laurent_series_free (expected_conj);
}

/* ncm_laurent_series_tps_new() must own order+1 independently-usable,
 * zero-initialized coefficients and report its own order back -- exactly
 * the self-describing-length property `NcmLaurentSeries **` plus a
 * separately-threaded `guint order` never had. */
static void
test_tps_new_owns_coefficients_and_reports_order (void)
{
  NcmLaurentSeriesTPS *tps = ncm_laurent_series_tps_new (2);
  NcmLaurentSeries *slot0, *slot2;

  g_assert_cmpuint (ncm_laurent_series_tps_order (tps), ==, 2);

  slot0 = ncm_laurent_series_tps_get (tps, 0);
  slot2 = ncm_laurent_series_tps_get (tps, 2);
  g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (slot0, 0)), <, 1.0e-13);

  ncm_laurent_series_set_single_into (slot0, 0, 5.0);
  ncm_laurent_series_set_single_into (slot2, 1, 7.0);
  g_assert_cmpfloat (fabs (creal (ncm_laurent_series_get_c (ncm_laurent_series_tps_get (tps, 0), 0)) - 5.0), <, 1.0e-13);
  g_assert_cmpfloat (fabs (creal (ncm_laurent_series_get_c (ncm_laurent_series_tps_get (tps, 2), 1)) - 7.0), <, 1.0e-13);

  ncm_laurent_series_tps_unref (tps);
}

/* The boxed "copy" is a ref-count increment (matching NcGalaxySDShapeData
 * and siblings), not a deep copy: a "copy" shares the same underlying
 * storage, so mutating through one reference must be visible through the
 * other. */
static void
test_tps_boxed_copy_is_refcount_not_deep_copy (void)
{
  NcmLaurentSeriesTPS *tps    = ncm_laurent_series_tps_new (1);
  NcmLaurentSeriesTPS *shared = (NcmLaurentSeriesTPS *) g_boxed_copy (NCM_TYPE_LAURENT_SERIES_TPS, tps);

  g_assert_true (shared == tps);

  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (tps, 0), 0, 3.0);
  g_assert_cmpfloat (fabs (creal (ncm_laurent_series_get_c (ncm_laurent_series_tps_get (shared, 0), 0)) - 3.0), <, 1.0e-13);

  g_boxed_free (NCM_TYPE_LAURENT_SERIES_TPS, shared);
  ncm_laurent_series_tps_unref (tps);
}

/* ncm_laurent_series_tps_conv() must be the truncated Cauchy product
 * (never extending to deg(a)+deg(b)) and must match plain scalar
 * convolution when every coefficient is a single-harmonic series. Also
 * exercises @out being reused as the destination of a second, independent
 * conv() call, since its private conv scratch must not leak stale state
 * across calls. */
static void
test_tps_conv_is_truncated_cauchy_product (void)
{
  const guint order        = 3;
  NcmLaurentSeriesTPS *a   = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *b   = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *out = ncm_laurent_series_tps_new (order);
  guint m;

  /* a(g) = 1 + 2g (higher terms zero), b(g) = 3 - g (higher terms zero),
   * every coefficient a plain real scalar (harmonic 0 only): the product
   * truncated at order 3 is 3 + 5g - 2g^2, zero above that. */
  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (a, 0), 0, 1.0);
  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (a, 1), 0, 2.0);

  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (b, 0), 0, 3.0);
  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (b, 1), 0, -1.0);

  ncm_laurent_series_tps_conv (out, a, b);

  {
    const gdouble expected[4] = {3.0, 5.0, -2.0, 0.0};

    for (m = 0; m <= order; m++)
      g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (ncm_laurent_series_tps_get (out, m), 0) - expected[m]), <, 1.0e-13);
  }

  /* Reuse @out for a second, unrelated conv(): a(g)=g, b(g)=g -> g^2. */
  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (a, 0), 0, 0.0);
  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (a, 1), 0, 1.0);
  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (b, 0), 0, 0.0);
  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (b, 1), 0, 1.0);

  ncm_laurent_series_tps_conv (out, a, b);

  {
    const gdouble expected2[4] = {0.0, 0.0, 1.0, 0.0};

    for (m = 0; m <= order; m++)
      g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (ncm_laurent_series_tps_get (out, m), 0) - expected2[m]), <, 1.0e-13);
  }

  ncm_laurent_series_tps_unref (a);
  ncm_laurent_series_tps_unref (b);
  ncm_laurent_series_tps_unref (out);
}

/* ncm_laurent_series_tps_conj/_add/_scale() are plain term-by-term
 * wrappers around the scalar `_into` primitives -- checked directly
 * against those primitives rather than re-deriving the math. */
static void
test_tps_conj_add_scale_match_scalar_ops (void)
{
  const guint order              = 1;
  NcmLaurentSeriesTPS *a         = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *b         = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *conj_out  = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *add_out   = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *scale_out = ncm_laurent_series_tps_new (order);
  guint m;

  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (a, 0), 1, 2.0 + 1.0 * I);
  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (a, 1), -1, 3.0 - 2.0 * I);
  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (b, 0), 1, -1.0 + 0.5 * I);
  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (b, 1), -1, 0.25);

  ncm_laurent_series_tps_conj (conj_out, a);
  ncm_laurent_series_tps_add (add_out, a, b, 2.0);
  ncm_laurent_series_tps_scale (scale_out, a, I);

  for (m = 0; m <= order; m++)
  {
    NcmLaurentSeries *expected_conj  = ncm_laurent_series_conj (ncm_laurent_series_tps_get (a, m));
    NcmLaurentSeries *expected_add   = ncm_laurent_series_add (ncm_laurent_series_tps_get (a, m), ncm_laurent_series_tps_get (b, m), 2.0);
    NcmLaurentSeries *expected_scale = ncm_laurent_series_scale_c (ncm_laurent_series_tps_get (a, m), I);
    gint h;

    for (h = -1; h <= 1; h++)
    {
      g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (ncm_laurent_series_tps_get (conj_out, m), h) - ncm_laurent_series_get_c (expected_conj, h)), <, 1.0e-13);
      g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (ncm_laurent_series_tps_get (add_out, m), h) - ncm_laurent_series_get_c (expected_add, h)), <, 1.0e-13);
      g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (ncm_laurent_series_tps_get (scale_out, m), h) - ncm_laurent_series_get_c (expected_scale, h)), <, 1.0e-13);
    }

    ncm_laurent_series_free (expected_conj);
    ncm_laurent_series_free (expected_add);
    ncm_laurent_series_free (expected_scale);
  }

  ncm_laurent_series_tps_unref (a);
  ncm_laurent_series_tps_unref (b);
  ncm_laurent_series_tps_unref (conj_out);
  ncm_laurent_series_tps_unref (add_out);
  ncm_laurent_series_tps_unref (scale_out);
}

static void
test_tps_eval_hand_computed (void)
{
  /* tps(w,g) = 1 + 2*w*g + (w+w^{-1})*g^2, evaluated at w=i, g=0.5. */
  NcmLaurentSeriesTPS *tps = ncm_laurent_series_tps_new (2);
  const complex double w   = I;
  const gdouble g          = 0.5;
  complex double expected;

  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (tps, 0), 0, 1.0);
  ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (tps, 1), 1, 2.0);
  {
    NcmLaurentSeries *c2 = ncm_laurent_series_tps_get (tps, 2);

    ncm_laurent_series_reset (c2, -1, 1);
    ncm_laurent_series_set_c (c2, -1, 1.0);
    ncm_laurent_series_set_c (c2, 1, 1.0);
  }

  expected = 1.0 + 2.0 * w * g + (w + 1.0 / w) * g * g;

  g_assert_cmpfloat (cabs (ncm_laurent_series_tps_eval_c (tps, w, g) - expected), <, 1.0e-12);

  ncm_laurent_series_tps_unref (tps);
}

/* Golden reference: sympy.series((2+0.5*g-0.3*g**2+0.1*g**3)**Rational(27,10), g, 0, 4),
 * i.e. a real (non-integer) exponent, exactly the case
 * ncm_laurent_series_tps_pow() exists for (NcGalaxyShapePopBeta's own
 * rho2^(alpha-1) composition). Every coefficient here is a plain real
 * scalar (harmonic 0 only), so this only checks the recursion itself, not
 * any Laurent-series bookkeeping (already covered by
 * test_tps_conv_is_truncated_cauchy_product() and friends). */
static void
test_tps_pow_matches_sympy_reference (void)
{
  const guint order         = 3;
  const gdouble a_coeffs[4] = {2.0, 0.5, -0.3, 0.1};
  const gdouble p           = 2.7;
  const gdouble expected[4] = {
    6.498019170849885, 4.386162940323673, -1.6996381393754232, -0.1868688169367064
  };
  NcmLaurentSeriesTPS *a   = ncm_laurent_series_tps_new (order);
  NcmLaurentSeriesTPS *out = ncm_laurent_series_tps_new (order);
  guint n;

  for (n = 0; n <= order; n++)
    ncm_laurent_series_set_single_into (ncm_laurent_series_tps_get (a, n), 0, a_coeffs[n]);

  ncm_laurent_series_tps_pow (out, a, p);

  for (n = 0; n <= order; n++)
    g_assert_cmpfloat (cabs (ncm_laurent_series_get_c (ncm_laurent_series_tps_get (out, n), 0) - expected[n]), <, 1.0e-10);

  ncm_laurent_series_tps_unref (a);
  ncm_laurent_series_tps_unref (out);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/laurent_series/add_scale_hand_computed", &test_add_scale_hand_computed);
  g_test_add_func ("/ncm/laurent_series/conv_hand_computed", &test_conv_hand_computed);
  g_test_add_func ("/ncm/laurent_series/conj_hand_computed", &test_conj_hand_computed);
  g_test_add_func ("/ncm/laurent_series/eval_hand_computed", &test_eval_hand_computed);
  g_test_add_func ("/ncm/laurent_series/jacobi_anger_matches_direct_integration", &test_jacobi_anger_matches_direct_integration);
  g_test_add_func ("/ncm/laurent_series/chi_taylor_matches_python_reference", &test_chi_taylor_matches_python_reference);
  g_test_add_func ("/ncm/laurent_series/introspectable_api_matches_native", &test_introspectable_api_matches_native);
  g_test_add_func ("/ncm/laurent_series/copy_is_independent", &test_copy_is_independent);
  g_test_add_func ("/ncm/laurent_series/reset_grows_only_when_needed", &test_reset_grows_only_when_needed);
  g_test_add_func ("/ncm/laurent_series/into_variants_match_returns_new", &test_into_variants_match_returns_new);
  g_test_add_func ("/ncm/laurent_series/tps_new_owns_coefficients_and_reports_order", &test_tps_new_owns_coefficients_and_reports_order);
  g_test_add_func ("/ncm/laurent_series/tps_boxed_copy_is_refcount_not_deep_copy", &test_tps_boxed_copy_is_refcount_not_deep_copy);
  g_test_add_func ("/ncm/laurent_series/tps_conv_is_truncated_cauchy_product", &test_tps_conv_is_truncated_cauchy_product);
  g_test_add_func ("/ncm/laurent_series/tps_conj_add_scale_match_scalar_ops", &test_tps_conj_add_scale_match_scalar_ops);
  g_test_add_func ("/ncm/laurent_series/tps_eval_hand_computed", &test_tps_eval_hand_computed);
  g_test_add_func ("/ncm/laurent_series/tps_pow_matches_sympy_reference", &test_tps_pow_matches_sympy_reference);

  return g_test_run ();
}

