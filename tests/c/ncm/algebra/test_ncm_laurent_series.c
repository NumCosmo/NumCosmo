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

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/laurent_series/add_scale_hand_computed", &test_add_scale_hand_computed);
  g_test_add_func ("/ncm/laurent_series/conv_hand_computed", &test_conv_hand_computed);
  g_test_add_func ("/ncm/laurent_series/conj_hand_computed", &test_conj_hand_computed);
  g_test_add_func ("/ncm/laurent_series/jacobi_anger_matches_direct_integration", &test_jacobi_anger_matches_direct_integration);
  g_test_add_func ("/ncm/laurent_series/chi_taylor_matches_python_reference", &test_chi_taylor_matches_python_reference);
  g_test_add_func ("/ncm/laurent_series/introspectable_api_matches_native", &test_introspectable_api_matches_native);
  g_test_add_func ("/ncm/laurent_series/copy_is_independent", &test_copy_is_independent);

  return g_test_run ();
}

