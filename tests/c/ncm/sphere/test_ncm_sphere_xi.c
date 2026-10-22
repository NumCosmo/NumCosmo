/***************************************************************************
 *            test_ncm_sphere_xi.c
 *
 *  Fri Oct 16 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_sphere_xi.c
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_LMAX 3000
#define TEST_N_SPEC 3

/* |K - K_ref| / ((2 ell + 1) / (4 pi)) on the certified rows; measured at most 5.8e-15 */
#define TEST_KERNEL_TOL 2.0e-14

/* |xi - sum_ell K C_ell| / sum_ell |K C_ell| between the BLAS products and a loop; at most
 * 7.3e-17 measured against an exact sum */
#define TEST_APPLY_TOL 1.0e-14

/* The Chebyshev test spectrum C_ell = 1 / ((ell + 1) (ell + 2)), ell = 2..300, as in the generator */
#define TEST_CHEB_LMIN 2
#define TEST_CHEB_LMAX 300

/* |c - c_ref| / max |c_ref| and |xi - xi_ref| / max |xi_ref| on the certified rows; measured at
 * most 2.2e-16 and 4.9e-16 */
#define TEST_CHEB_TOL 2.0e-15

typedef struct _TestNcmSphereXi
{
  NcmSphereXi *sxi;
  NcmVector *theta_edges;
} TestNcmSphereXi;

typedef struct _TestSphereXiRow
{
  guint bin;
  guint ell;
  gdouble K;
} TestSphereXiRow;

typedef struct _TestSphereXiChebRow
{
  guint k;
  gdouble c;
} TestSphereXiChebRow;

typedef struct _TestSphereXiThetaRow
{
  gdouble theta;
  gdouble xi;
} TestSphereXiThetaRow;

/*
 * Bin edges in radians and kernel entries K_{b ell} certified by
 * tests/tools/make_sphere_xi_cases.py (mpmath, 40 digits, rounded to double).
 */
static const gdouble test_edges[] = {
  /* BEGIN make_sphere_xi_cases.py edges */
  0.0,
  0.0002908882086657216,
  0.000727220521664304,
  0.0008726646259971647,
  0.002908882086657216,
  0.017453292519943295,
  0.07272205216643039,
  1.5,
  1.7,
  3.0,
  3.141592653589793,
  /* END make_sphere_xi_cases.py edges */
};

static const TestSphereXiRow test_rows[] = {
  /* BEGIN make_sphere_xi_cases.py rows */
  {0, 0, 0.07957747154594767},
  {0, 1, 0.23873240958770053},
  {0, 2, 0.3978873324790263},
  {0, 3, 0.5570422301196415},
  {0, 10, 1.671124958160788},
  {0, 100, 15.993363125861633},
  {0, 1000, 157.55455321093547},
  {0, 3000, 433.4911205648049},
  {1, 0, 0.07957747154594767},
  {1, 1, 0.23873237802431121},
  {1, 2, 0.39788717466209594},
  {1, 3, 0.5570417882323042},
  {1, 10, 1.6711128062894516},
  {1, 100, 15.982686409315699},
  {1, 1000, 147.28422122589203},
  {1, 3000, 208.43322150232683},
  {2, 0, 0.07957747154594767},
  {2, 1, 0.23873233762317397},
  {2, 2, 0.397886972656456},
  {2, 3, 0.5570412226167069},
  {2, 10, 1.6710972519478564},
  {2, 100, 15.969024591277325},
  {2, 1000, 134.55445657203532},
  {2, 3000, 0.31397024405232477},
  {3, 0, 0.07957747154594767},
  {3, 1, 0.23873186417266823},
  {3, 2, 0.3978846054077748},
  {3, 3, 0.5570345943365593},
  {3, 10, 1.6709149815020874},
  {3, 100, 15.809488717605772},
  {3, 1000, 30.40911501256944},
  {3, 3000, 15.475958130497435},
  {4, 0, 0.07957747154594767},
  {4, 1, 0.23871372957241008},
  {4, 2, 0.397793937149178},
  {4, 3, 0.5567807431305375},
  {4, 10, 1.663943013952967},
  {4, 100, 10.433857932895423},
  {4, 1000, -4.158833740445377},
  {4, 3000, 1.1670585909937032},
  {5, 0, 0.07957747154594767},
  {5, 1, 0.23839873975722445},
  {5, 2, 0.39622045782566945},
  {5, 3, 0.5523811671647694},
  {5, 10, 1.5456926669778683},
  {5, 100, -0.2527584477923056},
  {5, 1000, 0.2975725356104371},
  {5, 3000, -0.27088644632093706},
  {6, 0, 0.07957747154594767},
  {6, 1, 0.12749434396926823},
  {6, 2, 0.013980728902996781},
  {6, 3, -0.0744741212680555},
  {6, 10, 0.02416044017467019},
  {6, 100, -0.011150076712907658},
  {6, 1000, -0.004602275816306329},
  {6, 3000, 0.002740001469458611},
  {7, 0, 0.07957747154594767},
  {7, 1, -0.006936047138553807},
  {7, 2, -0.196458761616961},
  {7, 3, 0.02383909983562413},
  {7, 10, -0.32415464879640193},
  {7, 100, 0.07250919050065506},
  {7, 1000, 0.010931499738631426},
  {7, 3000, -0.021841000982825324},
  {8, 0, 0.07957747154594767},
  {8, 1, -0.13355132821136664},
  {8, 2, 0.02471698739219393},
  {8, 3, 0.07919580801509606},
  {8, 10, 0.029827691567201713},
  {8, 100, -0.010295562080515895},
  {8, 1000, 0.0012474640060778492},
  {8, 3000, 0.0025008220437200187},
  {9, 0, 0.07957747154594767},
  {9, 1, -0.23753785691230697},
  {9, 2, 0.3919344933362236},
  {9, 3, -0.5404576133675653},
  {9, 10, 1.2508917278799476},
  {9, 100, 0.3741199435454046},
  {9, 1000, 0.071056513287101},
  {9, 3000, 0.004416464066441027},
  /* END make_sphere_xi_cases.py rows */
};

/*
 * Chebyshev coefficients c_k and values xi(theta) of the Chebyshev test spectrum, certified by
 * tests/tools/make_sphere_xi_cases.py (mpmath, 40 digits, rounded to double).
 */
static const TestSphereXiChebRow test_cheb_rows[] = {
  /* BEGIN make_sphere_xi_cases.py cheb */
  {0, 0.01921884901620038},
  {1, 0.028479827313003938},
  {2, 0.04778511569190286},
  {3, 0.036675565163464124},
  {10, 0.013769659377249478},
  {100, 0.0012301208182954428},
  {299, 3.443238550627715e-05},
  {300, 3.42613637239612e-05},
  /* END make_sphere_xi_cases.py cheb */
};

static const TestSphereXiThetaRow test_xi_rows[] = {
  /* BEGIN make_sphere_xi_cases.py xi */
  {0.0, 0.6829263847987749},
  {0.0002908882086657216, 0.682775880073147},
  {0.01, 0.5479440497454116},
  {0.5, 0.014825360476154023},
  {1.5, -0.012505576595940605},
  {3.0, 0.017744473468416726},
  {3.141592653589793, 0.018359443529658556},
  /* END make_sphere_xi_cases.py xi */
};

void test_ncm_sphere_xi_new (TestNcmSphereXi *test, gconstpointer pdata);
void test_ncm_sphere_xi_free (TestNcmSphereXi *test, gconstpointer pdata);

void test_ncm_sphere_xi_kernel (TestNcmSphereXi *test, gconstpointer pdata);
void test_ncm_sphere_xi_apply (TestNcmSphereXi *test, gconstpointer pdata);
void test_ncm_sphere_xi_apply_vector (TestNcmSphereXi *test, gconstpointer pdata);
void test_ncm_sphere_xi_properties (TestNcmSphereXi *test, gconstpointer pdata);
void test_ncm_sphere_xi_serialize (TestNcmSphereXi *test, gconstpointer pdata);

void test_ncm_sphere_xi_chebyshev (void);
void test_ncm_sphere_xi_chebyshev_eval (void);
void test_ncm_sphere_xi_chebyshev_grid (void);
void test_ncm_sphere_xi_chebyshev_legendre (void);

void test_ncm_sphere_xi_traps (void);
void test_ncm_sphere_xi_invalid_lmax_zero (void);
void test_ncm_sphere_xi_invalid_edges_unset (void);
void test_ncm_sphere_xi_invalid_edges_short (void);
void test_ncm_sphere_xi_invalid_edges_order (void);
void test_ncm_sphere_xi_invalid_edges_range (void);
void test_ncm_sphere_xi_invalid_not_prepared (void);
void test_ncm_sphere_xi_invalid_stale (void);
void test_ncm_sphere_xi_invalid_peek (void);
void test_ncm_sphere_xi_invalid_ell_range (void);
void test_ncm_sphere_xi_invalid_xi_shape (void);
void test_ncm_sphere_xi_invalid_vector_ell_range (void);
void test_ncm_sphere_xi_invalid_vector_shape (void);
void test_ncm_sphere_xi_invalid_cheb_lmax_zero (void);
void test_ncm_sphere_xi_invalid_cheb_ell_range (void);
void test_ncm_sphere_xi_invalid_cheb_shape (void);
void test_ncm_sphere_xi_invalid_cheb_eval_shape (void);
void test_ncm_sphere_xi_invalid_cheb_grid_short (void);
void test_ncm_sphere_xi_invalid_cheb_grid_shape (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/sphere_xi/kernel", TestNcmSphereXi, NULL,
              &test_ncm_sphere_xi_new, &test_ncm_sphere_xi_kernel, &test_ncm_sphere_xi_free);
  g_test_add ("/ncm/sphere_xi/apply", TestNcmSphereXi, NULL,
              &test_ncm_sphere_xi_new, &test_ncm_sphere_xi_apply, &test_ncm_sphere_xi_free);
  g_test_add ("/ncm/sphere_xi/apply_vector", TestNcmSphereXi, NULL,
              &test_ncm_sphere_xi_new, &test_ncm_sphere_xi_apply_vector, &test_ncm_sphere_xi_free);
  g_test_add ("/ncm/sphere_xi/properties", TestNcmSphereXi, NULL,
              &test_ncm_sphere_xi_new, &test_ncm_sphere_xi_properties, &test_ncm_sphere_xi_free);
  g_test_add ("/ncm/sphere_xi/serialize", TestNcmSphereXi, NULL,
              &test_ncm_sphere_xi_new, &test_ncm_sphere_xi_serialize, &test_ncm_sphere_xi_free);

  g_test_add_func ("/ncm/sphere_xi/chebyshev", &test_ncm_sphere_xi_chebyshev);
  g_test_add_func ("/ncm/sphere_xi/chebyshev/eval", &test_ncm_sphere_xi_chebyshev_eval);
  g_test_add_func ("/ncm/sphere_xi/chebyshev/grid", &test_ncm_sphere_xi_chebyshev_grid);
  g_test_add_func ("/ncm/sphere_xi/chebyshev/legendre", &test_ncm_sphere_xi_chebyshev_legendre);

  g_test_add_func ("/ncm/sphere_xi/traps", &test_ncm_sphere_xi_traps);
  g_test_add_func ("/ncm/sphere_xi/invalid/lmax_zero/subprocess", &test_ncm_sphere_xi_invalid_lmax_zero);
  g_test_add_func ("/ncm/sphere_xi/invalid/edges_unset/subprocess", &test_ncm_sphere_xi_invalid_edges_unset);
  g_test_add_func ("/ncm/sphere_xi/invalid/edges_short/subprocess", &test_ncm_sphere_xi_invalid_edges_short);
  g_test_add_func ("/ncm/sphere_xi/invalid/edges_order/subprocess", &test_ncm_sphere_xi_invalid_edges_order);
  g_test_add_func ("/ncm/sphere_xi/invalid/edges_range/subprocess", &test_ncm_sphere_xi_invalid_edges_range);
  g_test_add_func ("/ncm/sphere_xi/invalid/not_prepared/subprocess", &test_ncm_sphere_xi_invalid_not_prepared);
  g_test_add_func ("/ncm/sphere_xi/invalid/stale/subprocess", &test_ncm_sphere_xi_invalid_stale);
  g_test_add_func ("/ncm/sphere_xi/invalid/peek/subprocess", &test_ncm_sphere_xi_invalid_peek);
  g_test_add_func ("/ncm/sphere_xi/invalid/ell_range/subprocess", &test_ncm_sphere_xi_invalid_ell_range);
  g_test_add_func ("/ncm/sphere_xi/invalid/xi_shape/subprocess", &test_ncm_sphere_xi_invalid_xi_shape);
  g_test_add_func ("/ncm/sphere_xi/invalid/vector_ell_range/subprocess", &test_ncm_sphere_xi_invalid_vector_ell_range);
  g_test_add_func ("/ncm/sphere_xi/invalid/vector_shape/subprocess", &test_ncm_sphere_xi_invalid_vector_shape);
  g_test_add_func ("/ncm/sphere_xi/invalid/cheb_lmax_zero/subprocess", &test_ncm_sphere_xi_invalid_cheb_lmax_zero);
  g_test_add_func ("/ncm/sphere_xi/invalid/cheb_ell_range/subprocess", &test_ncm_sphere_xi_invalid_cheb_ell_range);
  g_test_add_func ("/ncm/sphere_xi/invalid/cheb_shape/subprocess", &test_ncm_sphere_xi_invalid_cheb_shape);
  g_test_add_func ("/ncm/sphere_xi/invalid/cheb_eval_shape/subprocess", &test_ncm_sphere_xi_invalid_cheb_eval_shape);
  g_test_add_func ("/ncm/sphere_xi/invalid/cheb_grid_short/subprocess", &test_ncm_sphere_xi_invalid_cheb_grid_short);
  g_test_add_func ("/ncm/sphere_xi/invalid/cheb_grid_shape/subprocess", &test_ncm_sphere_xi_invalid_cheb_grid_shape);

  g_test_run ();

  return 0;
}

/* A prepared NcmSphereXi on two bins, for the abort tests */
static NcmSphereXi *
_test_two_bins (void)
{
  NcmVector *edges = ncm_vector_new_data_static ((gdouble *) test_edges, 3, 1);
  NcmSphereXi *sxi = ncm_sphere_xi_new (edges, 10);

  ncm_sphere_xi_prepare (sxi);
  ncm_vector_free (edges);

  return sxi;
}

/* Random spectra C_ell in [-1, 1] for ell = lmin, ..., lmax, one per row, from a fixed seed */
static NcmMatrix *
_test_spectra (const guint lmin, const guint lmax)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, 20261016);
  NcmMatrix *C_ells = ncm_matrix_new (TEST_N_SPEC, lmax - lmin + 1);
  guint r, j;

  for (r = 0; r < TEST_N_SPEC; r++)
    for (j = 0; j < lmax - lmin + 1; j++)
      ncm_matrix_set (C_ells, r, j, ncm_rng_uniform_gen (rng, -1.0, 1.0));

  ncm_rng_free (rng);

  return C_ells;
}

/*
 * Checks xi[b xi_stride] against sum_ell K_{b ell} C_ell summed in a loop, relative to
 * sum_ell |K C_ell|, with C_ell = C[j C_stride] for ell = lmin + j, j < n_ell
 */
static void
_test_check_sum (NcmMatrix *K, const gdouble *C, const guint C_stride, const guint n_ell, const guint lmin, const gdouble *xi, const guint xi_stride)
{
  guint b, j;

  for (b = 0; b < ncm_matrix_nrows (K); b++)
  {
    gdouble sum = 0.0, abs_sum = 0.0;

    for (j = 0; j < n_ell; j++)
    {
      const gdouble t = ncm_matrix_get (K, b, lmin + j) * C[j * C_stride];

      sum     += t;
      abs_sum += fabs (t);
    }

    g_assert_cmpfloat (fabs (xi[b * xi_stride] - sum), <=, TEST_APPLY_TOL * abs_sum);
  }
}

/* The Chebyshev coefficients of the Chebyshev test spectrum, one row */
static NcmMatrix *
_test_cheb_coeffs (NcmSphereXi *sxi)
{
  const guint n_ell = TEST_CHEB_LMAX - TEST_CHEB_LMIN + 1;
  NcmMatrix *C_ells = ncm_matrix_new (1, n_ell);
  NcmMatrix *c      = ncm_matrix_new (1, TEST_CHEB_LMAX + 1);
  guint j;

  for (j = 0; j < n_ell; j++)
  {
    const gdouble ell = TEST_CHEB_LMIN + j;

    ncm_matrix_set (C_ells, 0, j, 1.0 / ((ell + 1.0) * (ell + 2.0)));
  }

  ncm_sphere_xi_to_chebyshev (sxi, C_ells, TEST_CHEB_LMIN, c);
  ncm_matrix_free (C_ells);

  return c;
}

/* The largest |xi| among the certified values */
static gdouble
_test_xi_scale (void)
{
  gdouble scale = 0.0;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_xi_rows); i++)
    scale = MAX (scale, fabs (test_xi_rows[i].xi));

  return scale;
}

void
test_ncm_sphere_xi_new (TestNcmSphereXi *test, gconstpointer pdata)
{
  test->theta_edges = ncm_vector_new_data_static ((gdouble *) test_edges, G_N_ELEMENTS (test_edges), 1);
  test->sxi         = ncm_sphere_xi_new (test->theta_edges, TEST_LMAX);

  ncm_sphere_xi_prepare (test->sxi);
}

void
test_ncm_sphere_xi_free (TestNcmSphereXi *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_sphere_xi_free, test->sxi);
  ncm_vector_free (test->theta_edges);
}

/* The kernel against the certified entries */
void
test_ncm_sphere_xi_kernel (TestNcmSphereXi *test, gconstpointer pdata)
{
  NcmMatrix *K = ncm_sphere_xi_peek_kernel (test->sxi);
  guint i;

  g_assert_cmpuint (ncm_matrix_nrows (K), ==, G_N_ELEMENTS (test_edges) - 1);
  g_assert_cmpuint (ncm_matrix_ncols (K), ==, TEST_LMAX + 1);

  for (i = 0; i < G_N_ELEMENTS (test_rows); i++)
  {
    const TestSphereXiRow *row = &test_rows[i];
    const gdouble scale        = (2.0 * row->ell + 1.0) / (4.0 * M_PI);

    g_assert_cmpfloat (fabs (ncm_matrix_get (K, row->bin, row->ell) - row->K), <=, TEST_KERNEL_TOL * scale);
  }
}

/* ncm_sphere_xi_apply() against the loop sum, for spectra from ell = 2 and for a sub-range */
void
test_ncm_sphere_xi_apply (TestNcmSphereXi *test, gconstpointer pdata)
{
  NcmMatrix *K            = ncm_sphere_xi_peek_kernel (test->sxi);
  const guint n_bins      = ncm_sphere_xi_get_n_bins (test->sxi);
  const guint ranges[][2] = {
    {
      2, TEST_LMAX
    }, {
      0, TEST_LMAX
    }, {
      5, 100
    }
  };
  guint i, r;

  for (i = 0; i < G_N_ELEMENTS (ranges); i++)
  {
    NcmMatrix *C_ells = _test_spectra (ranges[i][0], ranges[i][1]);
    NcmMatrix *xi     = ncm_matrix_new (TEST_N_SPEC, n_bins);

    ncm_sphere_xi_apply (test->sxi, C_ells, ranges[i][0], xi);

    for (r = 0; r < TEST_N_SPEC; r++)
      _test_check_sum (K, ncm_matrix_ptr (C_ells, r, 0), 1, ncm_matrix_ncols (C_ells), ranges[i][0], ncm_matrix_ptr (xi, r, 0), 1);

    ncm_matrix_free (C_ells);
    ncm_matrix_free (xi);
  }
}

/* ncm_sphere_xi_apply_vector() against the loop sum, on a strided spectrum (a column) */
void
test_ncm_sphere_xi_apply_vector (TestNcmSphereXi *test, gconstpointer pdata)
{
  NcmMatrix *K       = ncm_sphere_xi_peek_kernel (test->sxi);
  const guint n_bins = ncm_sphere_xi_get_n_bins (test->sxi);
  const guint lmin   = 2;
  const guint n_ell  = TEST_LMAX - lmin + 1;
  NcmMatrix *C_ells  = _test_spectra (lmin, TEST_LMAX);
  NcmMatrix *C_ellsT = ncm_matrix_new (n_ell, TEST_N_SPEC);
  NcmVector *C_ell;
  NcmVector *xi = ncm_vector_new (n_bins);

  ncm_matrix_transpose_memcpy (C_ellsT, C_ells);
  C_ell = ncm_matrix_get_col (C_ellsT, 1);

  ncm_sphere_xi_apply_vector (test->sxi, C_ell, lmin, xi);
  _test_check_sum (K, ncm_matrix_ptr (C_ellsT, 0, 1), ncm_matrix_tda (C_ellsT), n_ell, lmin, ncm_vector_data (xi), 1);

  ncm_vector_free (C_ell);
  ncm_matrix_free (C_ells);
  ncm_matrix_free (C_ellsT);
  ncm_vector_free (xi);
}

/* Getters, the copy of the edges, and the kernel shape after changing lmax */
void
test_ncm_sphere_xi_properties (TestNcmSphereXi *test, gconstpointer pdata)
{
  NcmVector *edges = ncm_sphere_xi_peek_theta_edges (test->sxi);
  guint lmax;

  g_assert_cmpuint (ncm_sphere_xi_get_lmax (test->sxi), ==, TEST_LMAX);
  g_assert_cmpuint (ncm_sphere_xi_get_n_bins (test->sxi), ==, G_N_ELEMENTS (test_edges) - 1);
  g_assert_true (edges != test->theta_edges);
  ncm_assert_cmpdouble (ncm_vector_get (edges, 1), ==, test_edges[1]);

  g_object_get (test->sxi, "lmax", &lmax, NULL);
  g_assert_cmpuint (lmax, ==, TEST_LMAX);

  g_object_set (test->sxi, "lmax", 10, NULL);
  ncm_sphere_xi_prepare (test->sxi);
  g_assert_cmpuint (ncm_matrix_ncols (ncm_sphere_xi_peek_kernel (test->sxi)), ==, 11);

  ncm_sphere_xi_set_theta_edges (test->sxi, NULL);
  g_assert_cmpuint (ncm_sphere_xi_get_n_bins (test->sxi), ==, 0);
}

/* A serialized copy prepares the same kernel */
void
test_ncm_sphere_xi_serialize (TestNcmSphereXi *test, gconstpointer pdata)
{
  NcmSphereXi *dup = NCM_SPHERE_XI (ncm_serialize_global_dup_obj (G_OBJECT (test->sxi)));
  NcmMatrix *K     = ncm_sphere_xi_peek_kernel (test->sxi);
  NcmMatrix *K_dup;
  guint b, ell;

  ncm_sphere_xi_prepare (dup);
  K_dup = ncm_sphere_xi_peek_kernel (dup);

  for (b = 0; b < ncm_matrix_nrows (K); b++)
    for (ell = 0; ell <= TEST_LMAX; ell++)
      ncm_assert_cmpdouble (ncm_matrix_get (K_dup, b, ell), ==, ncm_matrix_get (K, b, ell));

  NCM_TEST_FREE (ncm_sphere_xi_free, dup);
}

/* ncm_sphere_xi_to_chebyshev() against the certified coefficients */
void
test_ncm_sphere_xi_chebyshev (void)
{
  NcmSphereXi *sxi = ncm_sphere_xi_new (NULL, TEST_CHEB_LMAX);
  NcmMatrix *c     = _test_cheb_coeffs (sxi);
  gdouble scale    = 0.0;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_cheb_rows); i++)
    scale = MAX (scale, fabs (test_cheb_rows[i].c));

  for (i = 0; i < G_N_ELEMENTS (test_cheb_rows); i++)
    g_assert_cmpfloat (fabs (ncm_matrix_get (c, 0, test_cheb_rows[i].k) - test_cheb_rows[i].c), <=, TEST_CHEB_TOL * scale);

  ncm_matrix_free (c);
  NCM_TEST_FREE (ncm_sphere_xi_free, sxi);
}

/* ncm_sphere_xi_eval_chebyshev() against the certified values, including theta = 0 and pi */
void
test_ncm_sphere_xi_chebyshev_eval (void)
{
  NcmSphereXi *sxi    = ncm_sphere_xi_new (NULL, TEST_CHEB_LMAX);
  NcmMatrix *c        = _test_cheb_coeffs (sxi);
  NcmVector *theta    = ncm_vector_new (G_N_ELEMENTS (test_xi_rows));
  NcmMatrix *xi       = ncm_matrix_new (1, G_N_ELEMENTS (test_xi_rows));
  const gdouble scale = _test_xi_scale ();
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_xi_rows); i++)
    ncm_vector_set (theta, i, test_xi_rows[i].theta);

  ncm_sphere_xi_eval_chebyshev (sxi, c, theta, xi);

  for (i = 0; i < G_N_ELEMENTS (test_xi_rows); i++)
    g_assert_cmpfloat (fabs (ncm_matrix_get (xi, 0, i) - test_xi_rows[i].xi), <=, TEST_CHEB_TOL * scale);

  ncm_vector_free (theta);
  ncm_matrix_free (xi);
  ncm_matrix_free (c);
  NCM_TEST_FREE (ncm_sphere_xi_free, sxi);
}

/*
 * ncm_sphere_xi_eval_chebyshev_grid(): the ends against the certified xi(0) and xi(pi), every
 * grid angle against ncm_sphere_xi_eval_chebyshev(), and a second call reusing the plan
 */
void
test_ncm_sphere_xi_chebyshev_grid (void)
{
  const guint n       = TEST_CHEB_LMAX + 1;
  NcmSphereXi *sxi    = ncm_sphere_xi_new (NULL, TEST_CHEB_LMAX);
  NcmMatrix *c        = _test_cheb_coeffs (sxi);
  NcmMatrix *grid     = ncm_matrix_new (1, n);
  NcmMatrix *xi       = ncm_matrix_new (1, n);
  NcmVector *theta    = ncm_vector_new (n);
  const gdouble scale = _test_xi_scale ();
  guint j, pass;

  for (j = 0; j < n; j++)
    ncm_vector_set (theta, j, j * M_PI / (n - 1));

  ncm_sphere_xi_eval_chebyshev (sxi, c, theta, xi);

  for (pass = 0; pass < 2; pass++)
  {
    ncm_sphere_xi_eval_chebyshev_grid (sxi, c, grid);

    g_assert_cmpfloat (fabs (ncm_matrix_get (grid, 0, 0) - test_xi_rows[0].xi), <=, TEST_CHEB_TOL * scale);
    g_assert_cmpfloat (fabs (ncm_matrix_get (grid, 0, n - 1) - test_xi_rows[G_N_ELEMENTS (test_xi_rows) - 1].xi), <=, TEST_CHEB_TOL * scale);

    for (j = 0; j < n; j++)
      g_assert_cmpfloat (fabs (ncm_matrix_get (grid, 0, j) - ncm_matrix_get (xi, 0, j)), <=, TEST_CHEB_TOL * scale);
  }

  ncm_vector_free (theta);
  ncm_matrix_free (grid);
  ncm_matrix_free (xi);
  ncm_matrix_free (c);
  NCM_TEST_FREE (ncm_sphere_xi_free, sxi);
}

/*
 * A single multipole, C_ell = 4 pi / (2 ell + 1) at ell = n, gives the Chebyshev coefficients
 * of P_n: P_2 = (T_0 + 3 T_2) / 4 and P_3 = (3 T_1 + 5 T_3) / 8.
 */
void
test_ncm_sphere_xi_chebyshev_legendre (void)
{
  const gdouble P2[] = {0.25, 0.0, 0.75};
  const gdouble P3[] = {0.0, 0.375, 0.0, 0.625};
  const gdouble *P[] = {P2, P3};
  guint n;

  for (n = 2; n <= 3; n++)
  {
    NcmSphereXi *sxi  = ncm_sphere_xi_new (NULL, n);
    NcmMatrix *C_ells = ncm_matrix_new (1, 1);
    NcmMatrix *c      = ncm_matrix_new (1, n + 1);
    guint k;

    ncm_matrix_set (C_ells, 0, 0, 4.0 * M_PI / (2.0 * n + 1.0));
    ncm_sphere_xi_to_chebyshev (sxi, C_ells, n, c);

    for (k = 0; k <= n; k++)
      g_assert_cmpfloat (fabs (ncm_matrix_get (c, 0, k) - P[n - 2][k]), <=, 1.0e-15);

    ncm_matrix_free (C_ells);
    ncm_matrix_free (c);
    NCM_TEST_FREE (ncm_sphere_xi_free, sxi);
  }
}

void
test_ncm_sphere_xi_traps (void)
{
  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/lmax_zero/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*lmax is zero*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/edges_unset/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*theta-edges is unset*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/edges_short/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*theta-edges has 1 entries, at least two are needed*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/edges_order/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*not strictly increasing at entry 1*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/edges_range/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*outside [0, pi]*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/not_prepared/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_xi_apply: the kernel is not prepared*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/stale/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_xi_apply_vector: the kernel is not prepared*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/peek/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_xi_peek_kernel: the kernel is not prepared*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/ell_range/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the spectra reach ell = 11, above lmax = 10*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/xi_shape/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*xi is 3 x 3, the result is 3 x 2*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/vector_ell_range/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the spectrum reaches ell = 11, above lmax = 10*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/vector_shape/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*xi has 3 entries, the result has 2*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/cheb_lmax_zero/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_xi_to_chebyshev: lmax is zero*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/cheb_ell_range/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_xi_to_chebyshev: the spectra reach ell = 11, above lmax = 10*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/cheb_shape/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_xi_to_chebyshev: c is 3 x 10, the result is 3 x 11*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/cheb_eval_shape/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_xi_eval_chebyshev: xi is 3 x 4, the result is 3 x 5*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/cheb_grid_short/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*c has 1 columns, at least two are needed*");

  g_test_trap_subprocess ("/ncm/sphere_xi/invalid/cheb_grid_shape/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_sphere_xi_eval_chebyshev_grid: xi is 3 x 10, the result is 3 x 11*");
}

void
test_ncm_sphere_xi_invalid_lmax_zero (void)
{
  NcmVector *edges = ncm_vector_new_data_static ((gdouble *) test_edges, 3, 1);
  NcmSphereXi *sxi = ncm_sphere_xi_new (edges, 0);

  ncm_sphere_xi_prepare (sxi);
}

void
test_ncm_sphere_xi_invalid_edges_unset (void)
{
  NcmSphereXi *sxi = ncm_sphere_xi_new (NULL, 10);

  ncm_sphere_xi_prepare (sxi);
}

void
test_ncm_sphere_xi_invalid_edges_short (void)
{
  NcmVector *edges = ncm_vector_new_data_static ((gdouble *) test_edges, 1, 1);
  NcmSphereXi *sxi = ncm_sphere_xi_new (edges, 10);

  ncm_sphere_xi_prepare (sxi);
}

void
test_ncm_sphere_xi_invalid_edges_order (void)
{
  const gdouble bad[] = {0.0, 0.2, 0.2, 0.3};
  NcmVector *edges    = ncm_vector_new_data_static ((gdouble *) bad, G_N_ELEMENTS (bad), 1);
  NcmSphereXi *sxi    = ncm_sphere_xi_new (edges, 10);

  ncm_sphere_xi_prepare (sxi);
}

void
test_ncm_sphere_xi_invalid_edges_range (void)
{
  const gdouble bad[] = {0.1, 0.2, 3.2};
  NcmVector *edges    = ncm_vector_new_data_static ((gdouble *) bad, G_N_ELEMENTS (bad), 1);
  NcmSphereXi *sxi    = ncm_sphere_xi_new (edges, 10);

  ncm_sphere_xi_prepare (sxi);
}

void
test_ncm_sphere_xi_invalid_not_prepared (void)
{
  NcmVector *edges  = ncm_vector_new_data_static ((gdouble *) test_edges, 3, 1);
  NcmSphereXi *sxi  = ncm_sphere_xi_new (edges, 10);
  NcmMatrix *C_ells = ncm_matrix_new (TEST_N_SPEC, 11);
  NcmMatrix *xi     = ncm_matrix_new (TEST_N_SPEC, 2);

  ncm_sphere_xi_apply (sxi, C_ells, 0, xi);
}

/* Changing lmax after the preparation leaves the kernel unprepared */
void
test_ncm_sphere_xi_invalid_stale (void)
{
  NcmSphereXi *sxi = _test_two_bins ();
  NcmVector *C_ell = ncm_vector_new (11);
  NcmVector *xi    = ncm_vector_new (2);

  ncm_sphere_xi_set_lmax (sxi, 20);
  ncm_sphere_xi_apply_vector (sxi, C_ell, 0, xi);
}

void
test_ncm_sphere_xi_invalid_peek (void)
{
  NcmSphereXi *sxi = _test_two_bins ();

  ncm_sphere_xi_set_theta_edges (sxi, NULL);
  ncm_sphere_xi_peek_kernel (sxi);
}

void
test_ncm_sphere_xi_invalid_ell_range (void)
{
  NcmSphereXi *sxi  = _test_two_bins ();
  NcmMatrix *C_ells = ncm_matrix_new (TEST_N_SPEC, 10);
  NcmMatrix *xi     = ncm_matrix_new (TEST_N_SPEC, 2);

  ncm_sphere_xi_apply (sxi, C_ells, 2, xi);
}

void
test_ncm_sphere_xi_invalid_xi_shape (void)
{
  NcmSphereXi *sxi  = _test_two_bins ();
  NcmMatrix *C_ells = ncm_matrix_new (TEST_N_SPEC, 11);
  NcmMatrix *xi     = ncm_matrix_new (TEST_N_SPEC, 3);

  ncm_sphere_xi_apply (sxi, C_ells, 0, xi);
}

void
test_ncm_sphere_xi_invalid_vector_ell_range (void)
{
  NcmSphereXi *sxi = _test_two_bins ();
  NcmVector *C_ell = ncm_vector_new (10);
  NcmVector *xi    = ncm_vector_new (2);

  ncm_sphere_xi_apply_vector (sxi, C_ell, 2, xi);
}

void
test_ncm_sphere_xi_invalid_vector_shape (void)
{
  NcmSphereXi *sxi = _test_two_bins ();
  NcmVector *C_ell = ncm_vector_new (11);
  NcmVector *xi    = ncm_vector_new (3);

  ncm_sphere_xi_apply_vector (sxi, C_ell, 0, xi);
}

void
test_ncm_sphere_xi_invalid_cheb_lmax_zero (void)
{
  NcmSphereXi *sxi  = ncm_sphere_xi_new (NULL, 0);
  NcmMatrix *C_ells = ncm_matrix_new (TEST_N_SPEC, 1);
  NcmMatrix *c      = ncm_matrix_new (TEST_N_SPEC, 1);

  ncm_sphere_xi_to_chebyshev (sxi, C_ells, 0, c);
}

void
test_ncm_sphere_xi_invalid_cheb_ell_range (void)
{
  NcmSphereXi *sxi  = ncm_sphere_xi_new (NULL, 10);
  NcmMatrix *C_ells = ncm_matrix_new (TEST_N_SPEC, 10);
  NcmMatrix *c      = ncm_matrix_new (TEST_N_SPEC, 11);

  ncm_sphere_xi_to_chebyshev (sxi, C_ells, 2, c);
}

void
test_ncm_sphere_xi_invalid_cheb_shape (void)
{
  NcmSphereXi *sxi  = ncm_sphere_xi_new (NULL, 10);
  NcmMatrix *C_ells = ncm_matrix_new (TEST_N_SPEC, 11);
  NcmMatrix *c      = ncm_matrix_new (TEST_N_SPEC, 10);

  ncm_sphere_xi_to_chebyshev (sxi, C_ells, 0, c);
}

void
test_ncm_sphere_xi_invalid_cheb_eval_shape (void)
{
  NcmSphereXi *sxi = ncm_sphere_xi_new (NULL, 10);
  NcmMatrix *c     = ncm_matrix_new (TEST_N_SPEC, 11);
  NcmVector *theta = ncm_vector_new (5);
  NcmMatrix *xi    = ncm_matrix_new (TEST_N_SPEC, 4);

  ncm_matrix_set_zero (c);
  ncm_vector_set_zero (theta);
  ncm_sphere_xi_eval_chebyshev (sxi, c, theta, xi);
}

void
test_ncm_sphere_xi_invalid_cheb_grid_short (void)
{
  NcmSphereXi *sxi = ncm_sphere_xi_new (NULL, 10);
  NcmMatrix *c     = ncm_matrix_new (TEST_N_SPEC, 1);
  NcmMatrix *xi    = ncm_matrix_new (TEST_N_SPEC, 1);

  ncm_sphere_xi_eval_chebyshev_grid (sxi, c, xi);
}

void
test_ncm_sphere_xi_invalid_cheb_grid_shape (void)
{
  NcmSphereXi *sxi = ncm_sphere_xi_new (NULL, 10);
  NcmMatrix *c     = ncm_matrix_new (TEST_N_SPEC, 11);
  NcmMatrix *xi    = ncm_matrix_new (TEST_N_SPEC, 10);

  ncm_sphere_xi_eval_chebyshev_grid (sxi, c, xi);
}

