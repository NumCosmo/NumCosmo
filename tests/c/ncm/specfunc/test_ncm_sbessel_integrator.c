/***************************************************************************
 *            test_ncm_sbessel_integrator.c
 *
 *  Sat September 26 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_sbessel_integrator.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * The NcmSBesselIntegrator contract, exercised through NcmSBesselIntegratorGL, and GL
 * against the Arb truth tables in the regime it is documented for.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

void test_ncm_sbessel_integrator_default_range (void);
void test_ncm_sbessel_integrator_ell_range (void);
void test_ncm_sbessel_integrator_integrate_ell (void);
void test_ncm_sbessel_integrator_gl_truth (void);
void test_ncm_sbessel_integrator_gl_far_bump (void);
void test_ncm_sbessel_integrator_gl_identity (void);
void test_ncm_sbessel_integrator_deriv_zero (void);
void test_ncm_sbessel_integrator_invalid_deriv_gl (void);
void test_ncm_sbessel_integrator_traps (void);
void test_ncm_sbessel_integrator_invalid_length (void);
void test_ncm_sbessel_integrator_invalid_ell (void);
void test_ncm_sbessel_integrator_invalid_deriv (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);

  /*
   * The GSL error handler stays off, the library default: deep in the evanescent tail GSL's
   * j_ell underflows to zero, the right value, and reports it as an error.
   */

  g_test_add_func ("/ncm/sbessel_integrator/default_range", test_ncm_sbessel_integrator_default_range);
  g_test_add_func ("/ncm/sbessel_integrator/ell_range", test_ncm_sbessel_integrator_ell_range);
  g_test_add_func ("/ncm/sbessel_integrator/integrate_ell", test_ncm_sbessel_integrator_integrate_ell);
  g_test_add_func ("/ncm/sbessel_integrator/gl/truth", test_ncm_sbessel_integrator_gl_truth);
  g_test_add_func ("/ncm/sbessel_integrator/gl/far_bump", test_ncm_sbessel_integrator_gl_far_bump);
  g_test_add_func ("/ncm/sbessel_integrator/gl/identity", test_ncm_sbessel_integrator_gl_identity);
  g_test_add_func ("/ncm/sbessel_integrator/deriv_zero", test_ncm_sbessel_integrator_deriv_zero);
  g_test_add_func ("/ncm/sbessel_integrator/traps", test_ncm_sbessel_integrator_traps);
  g_test_add_func ("/ncm/sbessel_integrator/invalid/length/subprocess", test_ncm_sbessel_integrator_invalid_length);
  g_test_add_func ("/ncm/sbessel_integrator/invalid/ell/subprocess", test_ncm_sbessel_integrator_invalid_ell);
  g_test_add_func ("/ncm/sbessel_integrator/invalid/deriv/subprocess", test_ncm_sbessel_integrator_invalid_deriv);
  g_test_add_func ("/ncm/sbessel_integrator/invalid/deriv_gl/subprocess", test_ncm_sbessel_integrator_invalid_deriv_gl);

  g_test_run ();
}

/* Without ell-range at construction the range is [0, 0] */
void
test_ncm_sbessel_integrator_default_range (void)
{
  NcmSBesselIntegrator *sbi = g_object_new (NCM_TYPE_SBESSEL_INTEGRATOR_GL, NULL);
  guint ell_min, ell_max;

  ncm_sbessel_integrator_get_ell_range (sbi, &ell_min, &ell_max);
  g_assert_cmpuint (ell_min, ==, 0);
  g_assert_cmpuint (ell_max, ==, 0);

  ncm_sbessel_integrator_free (sbi);
}

void
test_ncm_sbessel_integrator_ell_range (void)
{
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_gl_new (3, 17));
  NcmDTuple2 range          = NCM_DTUPLE2_STATIC_INIT (5.0, 9.0);
  NcmDTuple2 *got;
  guint ell_min, ell_max;

  ncm_sbessel_integrator_get_ell_range (sbi, &ell_min, &ell_max);
  g_assert_cmpuint (ell_min, ==, 3);
  g_assert_cmpuint (ell_max, ==, 17);

  g_object_set (sbi, "ell-range", &range, NULL);
  ncm_sbessel_integrator_get_ell_range (sbi, &ell_min, &ell_max);
  g_assert_cmpuint (ell_min, ==, 5);
  g_assert_cmpuint (ell_max, ==, 9);

  ncm_sbessel_integrator_set_ell_range (sbi, 2, 4);
  g_object_get (sbi, "ell-range", &got, NULL);
  g_assert_cmpfloat (got->elements[0], ==, 2.0);
  g_assert_cmpfloat (got->elements[1], ==, 4.0);

  ncm_dtuple2_free (got);
  ncm_sbessel_integrator_free (sbi);
}

/* The vector and the single-multipole calls agree, and the single one restores the range */
void
test_ncm_sbessel_integrator_integrate_ell (void)
{
  const guint ell_min       = 2;
  const guint ell_max       = 12;
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_gl_new (ell_min, ell_max));
  NcmVector *res            = ncm_vector_new (ell_max - ell_min + 1);
  guint ell, lo, hi;

  ncm_sbessel_integrator_integrate_gaussian (sbi, 0.5, 0.05, 0.1, 0.9, 7.0, res);

  for (ell = ell_min; ell <= ell_max; ell++)
  {
    const gdouble single = ncm_sbessel_integrator_integrate_gaussian_ell (sbi, 0.5, 0.05, 0.1, 0.9, 7.0, ell);

    g_assert_cmpfloat (single, ==, ncm_vector_get (res, ell - ell_min));
  }

  ncm_sbessel_integrator_get_ell_range (sbi, &lo, &hi);
  g_assert_cmpuint (lo, ==, ell_min);
  g_assert_cmpuint (hi, ==, ell_max);

  ncm_vector_free (res);
  ncm_sbessel_integrator_free (sbi);
}

typedef struct _TestTruth
{
  gboolean rational;
  guint ell;
  gdouble k;
  gdouble value;
} TestTruth;

/*
 * Entries of the Arb truth tables data/truth_tables/sbessel/gauss_jl_500.json.gz (Gaussian,
 * center 0.5, std 0.05, on [0.1, 0.9]) and rational_jl_500.json.gz (center 1.5, std 0.05,
 * on [0.01, 6.5]), with k b / nu between 0.05 and 1.45. Relative error measured at most
 * 2.4e-14.
 */
void
test_ncm_sbessel_integrator_gl_truth (void)
{
  const TestTruth truth[] = {
    {FALSE,   1,   0.5011872336272722,     0.010401531128943767},  /* gauss_jl_500[1][14]      */
    {FALSE,  10,    3.548133892335755,     3.930368033880559e-09}, /* gauss_jl_500[10][31]    */
    {FALSE,  10,   15.848931924611136,    0.0024113723097797055},  /* gauss_jl_500[10][44]     */
    {FALSE, 100,    5.623413251903491,   2.3565017876021603e-133}, /* gauss_jl_500[100][35]   */
    {FALSE, 100,   158.48931924611136,   2.8383605939301653e-05},  /* gauss_jl_500[100][64]    */
    {FALSE, 400,   125.89254117941672,   1.9474096904859334e-189}, /* gauss_jl_500[400][62]   */
    {FALSE, 400,    630.9573444801933,   2.6874336745760094e-06},  /* gauss_jl_500[400][76]    */
    {TRUE,    1,  0.19952623149688797,     0.013118394822050067},  /* rational_jl_500[1][6]    */
    {TRUE,   10,   0.5011872336272722,     5.632574931070647e-13}, /* rational_jl_500[10][14] */
    {TRUE,  100,    4.466835921509632,     5.402038523595476e-56}, /* rational_jl_500[100][33] */
    {TRUE,  100,   22.387211385683393,   2.4656192215801135e-12},  /* rational_jl_500[100][47] */
    {TRUE,  400,    89.12509381337455,     3.080838561205891e-13}, /* rational_jl_500[400][59] */
  };
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_gl_new (0, 0));
  guint i;

  for (i = 0; i < G_N_ELEMENTS (truth); i++)
  {
    const TestTruth *t = &truth[i];
    const gdouble val  = t->rational ?
                         ncm_sbessel_integrator_integrate_rational_ell (sbi, 1.5, 0.05, 0.01, 6.5, t->k, t->ell) :
                         ncm_sbessel_integrator_integrate_gaussian_ell (sbi, 0.5, 0.05, 0.1, 0.9, t->k, t->ell);

    g_assert_cmpfloat (fabs (val / t->value - 1.0), <, 2.0e-13);
  }

  ncm_sbessel_integrator_free (sbi);
}

/*
 * A bump far from the lower limit, small where the panels start: the old early exit
 * returned 4e-137 here. Levin, an independent method, agrees to 5e-13.
 */
void
test_ncm_sbessel_integrator_gl_far_bump (void)
{
  NcmSBesselIntegrator *gl = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_gl_new (10, 10));
  NcmSBesselIntegrator *lv = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (10, 10));
  const gdouble v_gl       = ncm_sbessel_integrator_integrate_gaussian_ell (gl, 500.0, 20.0, 0.0, 1000.0, 0.05, 10);
  const gdouble v_lv       = ncm_sbessel_integrator_integrate_gaussian_ell (lv, 500.0, 20.0, 0.0, 1000.0, 0.05, 10);

  g_assert_cmpfloat (fabs (v_gl / v_lv - 1.0), <, 1.0e-10);

  ncm_sbessel_integrator_free (gl);
  ncm_sbessel_integrator_free (lv);
}

typedef struct _TestPower
{
  gdouble b;
  gint ell;
} TestPower;

static gdouble
_test_power (gpointer user_data, gdouble chi, gdouble k)
{
  TestPower *p = (TestPower *) user_data;

  return gsl_pow_int (chi / p->b, p->ell + 2);
}

/*
 * d/dx [x^(l+2) j_(l+1)(x)] = x^(l+2) j_l(x), so with F = (chi / b)^(l+2) and k = 1 the
 * integral is j_(l+1)(b) - (a / b)^(l+2) j_(l+1)(a), exact, with j_(l+1) from
 * ncm_sf_sbessel(). Ranges up to 1.4 nu; relative error measured at most 5.4e-14.
 */
void
test_ncm_sbessel_integrator_gl_identity (void)
{
  const gint ell_a[]         = {5, 50, 200};
  const gdouble range_a[][2] = {
    {
      0.1, 0.9
    }, {
      0.1, 1.4
    }, {
      0.5, 1.2
    }
  };
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_gl_new (0, 0));
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (ell_a); i++)
  {
    const gdouble nu = ell_a[i] + 0.5;

    for (j = 0; j < G_N_ELEMENTS (range_a); j++)
    {
      const gdouble a     = range_a[j][0] * nu;
      const gdouble b     = range_a[j][1] * nu;
      TestPower p         = {b, ell_a[i]};
      const gdouble truth = ncm_sf_sbessel (ell_a[i] + 1, b) - gsl_pow_int (a / b, ell_a[i] + 2) * ncm_sf_sbessel (ell_a[i] + 1, a);
      const gdouble val   = ncm_sbessel_integrator_integrate_ell (sbi, &_test_power, a, b, 1.0, ell_a[i], &p);

      g_assert_cmpfloat (fabs (val / truth - 1.0), <, 5.0e-13);
    }
  }

  ncm_sbessel_integrator_free (sbi);
  ncm_mpsf_sbessel_free_cache ();
}

static gdouble
_test_gauss (gpointer user_data, gdouble chi, gdouble k)
{
  const gdouble z = (chi - 40.0) / 10.0;

  return exp (-0.5 * z * z);
}

/* Order zero dispatches to ncm_sbessel_integrator_integrate(), bit for bit */
void
test_ncm_sbessel_integrator_deriv_zero (void)
{
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_gl_new (0, 8));
  NcmVector *res_a          = ncm_vector_new (9);
  NcmVector *res_b          = ncm_vector_new (9);
  guint i;

  ncm_sbessel_integrator_integrate_deriv (sbi, &_test_gauss, 5.0, 80.0, 1.0, 0, res_a, NULL);
  ncm_sbessel_integrator_integrate (sbi, &_test_gauss, 5.0, 80.0, 1.0, res_b, NULL);

  for (i = 0; i < 9; i++)
    g_assert_cmpfloat (ncm_vector_get (res_a, i), ==, ncm_vector_get (res_b, i));

  ncm_vector_free (res_a);
  ncm_vector_free (res_b);
  ncm_sbessel_integrator_free (sbi);
}

void
test_ncm_sbessel_integrator_traps (void)
{
  g_test_trap_subprocess ("/ncm/sbessel_integrator/invalid/length/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*result has length*");

  g_test_trap_subprocess ("/ncm/sbessel_integrator/invalid/ell/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*negative multipole*");

  g_test_trap_subprocess ("/ncm/sbessel_integrator/invalid/deriv/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*not supported*");

  /* Only Levin implements the derivative weights */
  g_test_trap_subprocess ("/ncm/sbessel_integrator/invalid/deriv_gl/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*not implemented*");
}

void
test_ncm_sbessel_integrator_invalid_length (void)
{
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_gl_new (0, 4));
  NcmVector *res            = ncm_vector_new (4);

  ncm_sbessel_integrator_integrate_gaussian (sbi, 0.5, 0.05, 0.1, 0.9, 1.0, res);
}

void
test_ncm_sbessel_integrator_invalid_ell (void)
{
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_gl_new (0, 4));

  ncm_sbessel_integrator_integrate_gaussian_ell (sbi, 0.5, 0.05, 0.1, 0.9, 1.0, -1);
}

void
test_ncm_sbessel_integrator_invalid_deriv (void)
{
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_gl_new (0, 4));
  NcmVector *res            = ncm_vector_new (5);

  ncm_sbessel_integrator_integrate_deriv (sbi, NULL, 0.1, 0.9, 1.0, 3, res, NULL);
}

void
test_ncm_sbessel_integrator_invalid_deriv_gl (void)
{
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_gl_new (0, 3));
  NcmVector *res            = ncm_vector_new (4);

  ncm_sbessel_integrator_integrate_deriv (sbi, &_test_gauss, 1.0, 2.0, 1.0, 1, res, NULL);
}

