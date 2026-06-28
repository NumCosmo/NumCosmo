/***************************************************************************
 *            test_nc_wl_ellipticity.c
 *
 *  Tue Jun 24 2026
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

#include <complex.h>
#include <math.h>
#include <glib.h>
#include <glib-object.h>

typedef complex double (*ShearMap) (complex double g, complex double z);

typedef struct _TestShearConv
{
  ShearMap apply;
  ShearMap apply_inv;
  gdouble (*lndet_jac) (complex double g, complex double z_obs);
} TestShearConv;

/* The two conventions, addressed uniformly so the algebraic identities can be
 * checked with a single body. */
static const TestShearConv test_conv_trace = {
  nc_wl_ellipticity_apply_shear_trace_c,
  nc_wl_ellipticity_apply_shear_inv_trace_c,
  nc_wl_ellipticity_lndet_jac_trace_c,
};

static const TestShearConv test_conv_trace_det = {
  nc_wl_ellipticity_apply_shear_trace_det_c,
  nc_wl_ellipticity_apply_shear_inv_trace_det_c,
  nc_wl_ellipticity_lndet_jac_trace_det_c,
};

/* A spread of reduced shears that exercises both the |g| <= 1 and the |g| > 1
 * (strong-lensing) branches, including the trivial g = 0 case. */
static const complex double test_shears[] = {
  0.0 + 0.0 * I,
  0.10 + 0.05 * I,
  -0.30 + 0.20 * I,
  0.45 - 0.40 * I,
  1.50 + 0.30 * I,
  -1.20 - 0.60 * I,
};

/* A spread of observed ellipticities/distortions. */
static const complex double test_ellips[] = {
  0.00 + 0.00 * I,
  0.20 + 0.10 * I,
  -0.15 + 0.25 * I,
  0.30 - 0.22 * I,
};

/* Numerical log|det J| of @map (with @g held fixed) at the point @z0, by central
 * differences on the real and imaginary directions. */
static gdouble
test_numerical_lndet_jac (ShearMap map, complex double g, complex double z0)
{
  const gdouble h          = 1.0e-6;
  const complex double dfx = (map (g, z0 + h) - map (g, z0 - h)) / (2.0 * h);
  const complex double dfy = (map (g, z0 + I * h) - map (g, z0 - I * h)) / (2.0 * h);
  const gdouble J00        = creal (dfx);
  const gdouble J10        = cimag (dfx);
  const gdouble J01        = creal (dfy);
  const gdouble J11        = cimag (dfy);
  const gdouble det        = J00 * J11 - J01 * J10;

  return log (fabs (det));
}

static void
test_assert_cmplx_close (complex double a, complex double b, gdouble tol)
{
  g_assert_cmpfloat (cabs (a - b), <, tol);
}

/* apply_shear_inv must invert apply_shear for every shear/ellipticity pair, in
 * both branches. */
static void
test_nc_wl_ellipticity_round_trip (gconstpointer pdata)
{
  const TestShearConv *conv = pdata;
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (test_shears); i++)
  {
    const complex double g = test_shears[i];

    for (j = 0; j < G_N_ELEMENTS (test_ellips); j++)
    {
      const complex double z_obs = test_ellips[j];
      const complex double z_int = conv->apply_inv (g, z_obs);
      const complex double z_rec = conv->apply (g, z_int);

      test_assert_cmplx_close (z_rec, z_obs, 1.0e-9);
    }
  }
}

/* lndet_jac must equal the log|det| of the observed -> intrinsic map
 * (apply_shear_inv) evaluated at the observed point: this is exactly the
 * Jacobian factor the likelihood uses as a change of variables. Checked against
 * a finite-difference Jacobian, in both branches. */
static void
test_nc_wl_ellipticity_lndet_jac (gconstpointer pdata)
{
  const TestShearConv *conv = pdata;
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (test_shears); i++)
  {
    const complex double g = test_shears[i];

    for (j = 0; j < G_N_ELEMENTS (test_ellips); j++)
    {
      const complex double z_obs    = test_ellips[j];
      const gdouble lndet_ana       = conv->lndet_jac (g, z_obs);
      const gdouble lndet_num       = test_numerical_lndet_jac (conv->apply_inv, g, z_obs);

      g_assert_cmpfloat (fabs (lndet_ana - lndet_num), <, 1.0e-6);
    }
  }
}

/* The g = 0 limit must be the identity, with a unit Jacobian. */
static void
test_nc_wl_ellipticity_zero_shear (gconstpointer pdata)
{
  const TestShearConv *conv = pdata;
  const complex double g    = 0.0;
  guint j;

  for (j = 0; j < G_N_ELEMENTS (test_ellips); j++)
  {
    const complex double z = test_ellips[j];

    test_assert_cmplx_close (conv->apply (g, z), z, 1.0e-15);
    test_assert_cmplx_close (conv->apply_inv (g, z), z, 1.0e-15);
    g_assert_cmpfloat (fabs (conv->lndet_jac (g, z)), <, 1.0e-12);
  }
}

/* The introspectable NcmComplex wrappers must agree with the inline complex
 * double kernels, both branches. Exercises nc_wl_ellipticity.c. */
static void
test_nc_wl_ellipticity_introspectable_wrappers (void)
{
  /* one |g| <= 1 and one |g| > 1 shear */
  const complex double gs[] = { 0.2 + 0.1 * I, 1.4 - 0.3 * I };
  const complex double e    = 0.25 - 0.15 * I;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (gs); i++)
  {
    const complex double g = gs[i];
    NcmComplex *cg         = ncm_complex_new ();
    NcmComplex *ce         = ncm_complex_new ();
    NcmComplex *cout       = ncm_complex_new ();

    ncm_complex_set_c (cg, g);
    ncm_complex_set_c (ce, e);

    /* TRACE convention */
    nc_wl_ellipticity_apply_shear_trace (cg, ce, cout);
    test_assert_cmplx_close (ncm_complex_c (cout), nc_wl_ellipticity_apply_shear_trace_c (g, e), 1.0e-15);

    nc_wl_ellipticity_apply_shear_inv_trace (cg, ce, cout);
    test_assert_cmplx_close (ncm_complex_c (cout), nc_wl_ellipticity_apply_shear_inv_trace_c (g, e), 1.0e-15);

    g_assert_cmpfloat (fabs (nc_wl_ellipticity_lndet_jac_trace (cg, ce) -
                             nc_wl_ellipticity_lndet_jac_trace_c (g, e)), <, 1.0e-15);

    /* TRACE_DET convention */
    nc_wl_ellipticity_apply_shear_trace_det (cg, ce, cout);
    test_assert_cmplx_close (ncm_complex_c (cout), nc_wl_ellipticity_apply_shear_trace_det_c (g, e), 1.0e-15);

    nc_wl_ellipticity_apply_shear_inv_trace_det (cg, ce, cout);
    test_assert_cmplx_close (ncm_complex_c (cout), nc_wl_ellipticity_apply_shear_inv_trace_det_c (g, e), 1.0e-15);

    g_assert_cmpfloat (fabs (nc_wl_ellipticity_lndet_jac_trace_det (cg, ce) -
                             nc_wl_ellipticity_lndet_jac_trace_det_c (g, e)), <, 1.0e-15);

    ncm_complex_free (cg);
    ncm_complex_free (ce);
    ncm_complex_free (cout);
  }
}

/* The frame parity must be the identity for CELESTIAL, a complex conjugation for
 * CARTESIAN, an involution, and value-compatible with the legacy enum. */
static void
test_nc_wl_ellipticity_frame_parity (void)
{
  guint j;

  const gdouble phi = 0.37;

  /* values must match the legacy NcGalaxyWLObsCoord (CELESTIAL=0, EUCLIDEAN=1) */
  g_assert_cmpint (NC_WL_ELLIPTICITY_FRAME_CELESTIAL, ==, 0);
  g_assert_cmpint (NC_WL_ELLIPTICITY_FRAME_CARTESIAN, ==, 1);

  /* position angle: identity for celestial, reflected for cartesian */
  g_assert_cmpfloat (fabs (nc_wl_ellipticity_celestial_to_frame_angle (NC_WL_ELLIPTICITY_FRAME_CELESTIAL, phi) - phi), <, 1.0e-15);
  g_assert_cmpfloat (fabs (nc_wl_ellipticity_celestial_to_frame_angle (NC_WL_ELLIPTICITY_FRAME_CARTESIAN, phi) - (M_PI - phi)), <, 1.0e-15);

  for (j = 0; j < G_N_ELEMENTS (test_ellips); j++)
  {
    const complex double e = test_ellips[j];

    test_assert_cmplx_close (nc_wl_ellipticity_celestial_to_frame_c (NC_WL_ELLIPTICITY_FRAME_CELESTIAL, e), e, 1.0e-15);
    test_assert_cmplx_close (nc_wl_ellipticity_celestial_to_frame_c (NC_WL_ELLIPTICITY_FRAME_CARTESIAN, e), conj (e), 1.0e-15);

    /* involution: applying twice returns the original, in either frame */
    test_assert_cmplx_close (
      nc_wl_ellipticity_celestial_to_frame_c (NC_WL_ELLIPTICITY_FRAME_CARTESIAN,
                                              nc_wl_ellipticity_celestial_to_frame_c (NC_WL_ELLIPTICITY_FRAME_CARTESIAN, e)),
      e, 1.0e-15);
  }
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_data_func ("/nc/wl_ellipticity/trace/round_trip",
                        &test_conv_trace, &test_nc_wl_ellipticity_round_trip);
  g_test_add_data_func ("/nc/wl_ellipticity/trace/lndet_jac",
                        &test_conv_trace, &test_nc_wl_ellipticity_lndet_jac);
  g_test_add_data_func ("/nc/wl_ellipticity/trace/zero_shear",
                        &test_conv_trace, &test_nc_wl_ellipticity_zero_shear);

  g_test_add_data_func ("/nc/wl_ellipticity/trace_det/round_trip",
                        &test_conv_trace_det, &test_nc_wl_ellipticity_round_trip);
  g_test_add_data_func ("/nc/wl_ellipticity/trace_det/lndet_jac",
                        &test_conv_trace_det, &test_nc_wl_ellipticity_lndet_jac);
  g_test_add_data_func ("/nc/wl_ellipticity/trace_det/zero_shear",
                        &test_conv_trace_det, &test_nc_wl_ellipticity_zero_shear);

  g_test_add_func ("/nc/wl_ellipticity/introspectable_wrappers",
                   &test_nc_wl_ellipticity_introspectable_wrappers);
  g_test_add_func ("/nc/wl_ellipticity/frame_parity",
                   &test_nc_wl_ellipticity_frame_parity);

  return g_test_run ();
}
