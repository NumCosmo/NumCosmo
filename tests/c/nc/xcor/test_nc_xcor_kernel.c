/***************************************************************************
 *            test_nc_xcor_kernel.c
 *
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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
 * The #NcXcorKernel surface every concrete kernel has to implement: redshift range,
 * Limber evaluation, the component list and its limits, the k range, noise. Each kernel
 * is built from the cheapest input that is still its own -- a Gaussian dn/dz spline
 * rather than a photo-z binning, an analytic power spectrum rather than CLASS -- because
 * what is under test is the interface, not the astrophysics. Accuracy of the windows is
 * checked in tests/python/nc/xcor against certified values.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_DNDZ_MEAN 0.6
#define TEST_DNDZ_SIGMA 0.15
#define TEST_DNDZ_KNOTS 128

/* Most kernels live in a photo-z bin; the CMB ones reach last scattering, so they get a
 * distance object that covers it. */
#define TEST_ZMAX 6.0
#define TEST_ZMAX_LSS 1200.0
#define TEST_LMAX 64

typedef struct _TestNcXcorKernel
{
  NcHICosmo *cosmo;
  NcDistance *dist;
  NcmPowspec *ps;
  NcRecomb *recomb;
  NcXcorKernel *xclk;
} TestNcXcorKernel;

/**
 * TestKernel:
 * @name: path component naming the kernel
 * @build: builds it on the fixture's distance, power spectrum and recombination
 * @zf: how far the fixture's distance object has to reach for this kernel
 * @has_noise: whether add_noise() is expected to change anything
 */
typedef struct _TestKernel
{
  const gchar *name;

  NcXcorKernel *(*build) (TestNcXcorKernel *test);

  gdouble zf;
  gboolean has_noise;
} TestKernel;

/* A normalized Gaussian dn/dz on [0, 2]: a photo-z bin's shape without a photo-z model. */
static NcmSpline *
_dndz_new (void)
{
  NcmVector *z = ncm_vector_new (TEST_DNDZ_KNOTS);
  NcmVector *n = ncm_vector_new (TEST_DNDZ_KNOTS);
  NcmSpline *s;
  guint i;

  for (i = 0; i < TEST_DNDZ_KNOTS; i++)
  {
    const gdouble zi = 2.0 * i / (gdouble) (TEST_DNDZ_KNOTS - 1);
    const gdouble t  = (zi - TEST_DNDZ_MEAN) / TEST_DNDZ_SIGMA;

    ncm_vector_set (z, i, zi);
    ncm_vector_set (n, i, exp (-0.5 * t * t));
  }

  s = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (z, n, TRUE));

  ncm_vector_free (z);
  ncm_vector_free (n);

  return s;
}

/* Flat noise spectrum, long enough for the multipoles the noise check touches. */
static NcmVector *
_nl_new (void)
{
  NcmVector *nl = ncm_vector_new (TEST_LMAX + 1);

  ncm_vector_set_all (nl, 1.0e-8);

  return nl;
}

static NcXcorKernel *
_build_gal (TestNcXcorKernel *test)
{
  NcmSpline *dndz    = _dndz_new ();
  NcXcorKernelGal *g = nc_xcor_kernel_gal_new (test->dist, test->ps, 1, 1.234, dndz, FALSE);

  ncm_model_orig_vparam_set (NCM_MODEL (g), NC_XCOR_KERNEL_GAL_BIAS, 0, 1.5);
  ncm_spline_free (dndz);

  return NC_XCOR_KERNEL (g);
}

static NcXcorKernel *
_build_gal_magbias (TestNcXcorKernel *test)
{
  NcmSpline *dndz    = _dndz_new ();
  NcXcorKernelGal *g = nc_xcor_kernel_gal_new (test->dist, test->ps, 1, 1.234, dndz, TRUE);

  ncm_model_orig_vparam_set (NCM_MODEL (g), NC_XCOR_KERNEL_GAL_BIAS, 0, 1.5);
  ncm_spline_free (dndz);

  return NC_XCOR_KERNEL (g);
}

static NcXcorKernel *
_build_weak_lensing (TestNcXcorKernel *test)
{
  NcmSpline *dndz = _dndz_new ();
  NcXcorKernel *k = NC_XCOR_KERNEL (nc_xcor_kernel_weak_lensing_new (test->dist, test->ps, dndz, 3.0, 0.3));

  ncm_spline_free (dndz);

  return k;
}

static NcXcorKernel *
_build_cmb_lensing (TestNcXcorKernel *test)
{
  NcmVector *nl   = _nl_new ();
  NcXcorKernel *k = NC_XCOR_KERNEL (nc_xcor_kernel_cmb_lensing_new (test->dist, test->ps, test->recomb, nl));

  ncm_vector_free (nl);

  return k;
}

static NcXcorKernel *
_build_cmb_isw (TestNcXcorKernel *test)
{
  NcmVector *nl   = _nl_new ();
  NcXcorKernel *k = NC_XCOR_KERNEL (nc_xcor_kernel_cmb_isw_new (test->dist, test->ps, test->recomb, nl));

  ncm_vector_free (nl);

  return k;
}

static NcXcorKernel *
_build_tsz (TestNcXcorKernel *test)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_tsz_new (test->dist, test->ps, 2.0));
}

static NcXcorKernel *
_build_cluster_tophat (TestNcXcorKernel *test)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_cluster_tophat_new (test->dist, test->ps, 0.2, 0.8));
}

static NcXcorKernel *
_build_table (TestNcXcorKernel *test)
{
  NcmVector *chi = ncm_vector_new (64);
  NcmVector *W   = ncm_vector_new (64);
  NcXcorKernel *k;
  guint i;

  for (i = 0; i < 64; i++)
  {
    const gdouble chi_i = 1000.0 + 20.0 * i;
    const gdouble t     = (chi_i - 1630.0) / 300.0;

    ncm_vector_set (chi, i, chi_i);
    ncm_vector_set (W, i, exp (-0.5 * t * t));
  }

  k = NC_XCOR_KERNEL (nc_xcor_kernel_table_new (test->dist, test->ps, chi, W));

  ncm_vector_free (chi);
  ncm_vector_free (W);

  return k;
}

static const TestKernel test_kernels[] = {
  {"gal",            &_build_gal,            TEST_ZMAX,     TRUE },
  {"gal_magbias",    &_build_gal_magbias,    TEST_ZMAX,     TRUE },
  {"weak_lensing",   &_build_weak_lensing,   TEST_ZMAX,     TRUE },
  {"cmb_lensing",    &_build_cmb_lensing,    TEST_ZMAX_LSS, TRUE },
  {"cmb_isw",        &_build_cmb_isw,        TEST_ZMAX_LSS, TRUE },
  {"tsz",            &_build_tsz,            TEST_ZMAX,     FALSE},
  {"cluster_tophat", &_build_cluster_tophat, TEST_ZMAX,     FALSE},
  {"table",          &_build_table,          TEST_ZMAX,     FALSE},
};

static void
test_nc_xcor_kernel_new (TestNcXcorKernel *test, gconstpointer pdata)
{
  const TestKernel *tk = pdata;
  NcHICosmo *cosmo     = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_H0,      70.0);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_C, 0.255);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_X, 0.7);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_OMEGA_B, 0.045);
  ncm_model_orig_param_set (NCM_MODEL (cosmo), NC_HICOSMO_DE_XCDM_W, -1.0);
  nc_hicosmo_de_omega_x2omega_k (NC_HICOSMO_DE (cosmo), NULL);
  ncm_model_param_set_by_name (NCM_MODEL (cosmo), "Omegak", 0.0, NULL);

  test->cosmo = cosmo;
  test->dist  = nc_distance_new (tk->zf);
  test->ps    = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                       NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  test->recomb = NULL;

  nc_distance_prepare (test->dist, cosmo);
  ncm_powspec_prepare (test->ps, NCM_MODEL (cosmo));

  /* Recombination is only what the CMB kernels need, and solving it is the most
   * expensive thing in this fixture -- about as much as everything else together. Built
   * on demand, it is paid by two kernels instead of all eight. */
  if (tk->zf == TEST_ZMAX_LSS)
  {
    test->recomb = NC_RECOMB (nc_recomb_seager_new ());
    nc_recomb_prepare (test->recomb, cosmo);
  }

  test->xclk = tk->build (test);

  g_assert_true (NC_IS_XCOR_KERNEL (test->xclk));

  nc_xcor_kernel_prepare (test->xclk, cosmo);
}

static void
test_nc_xcor_kernel_free (TestNcXcorKernel *test, gconstpointer pdata)
{
  nc_recomb_clear (&test->recomb);
  ncm_powspec_free (test->ps);
  nc_distance_free (test->dist);
  nc_hicosmo_free (test->cosmo);

  NCM_TEST_FREE (nc_xcor_kernel_free, test->xclk);
}

static void
test_nc_xcor_kernel_basic (TestNcXcorKernel *test, gconstpointer pdata)
{
  NcXcorKernel *xclk2 = nc_xcor_kernel_ref (test->xclk);

  g_assert_true (nc_xcor_kernel_peek_dist (test->xclk) == test->dist);
  g_assert_true (nc_xcor_kernel_peek_powspec (test->xclk) == test->ps);

  g_assert_cmpuint (nc_xcor_kernel_obs_len (test->xclk), >, 0);
  g_assert_cmpuint (nc_xcor_kernel_obs_params_len (test->xclk), >=, 0);

  nc_xcor_kernel_clear (&xclk2);
  g_assert_true (xclk2 == NULL);
}

/* Every knob the adaptive machinery reads is settable, and reads back what was set. The
 * values are not the library's defaults precisely so a setter that silently ignores its
 * argument fails here. */
static void
test_nc_xcor_kernel_knobs (TestNcXcorKernel *test, gconstpointer pdata)
{
  NcXcorKernel *xclk = test->xclk;

  nc_xcor_kernel_set_lmax (xclk, TEST_LMAX);
  nc_xcor_kernel_set_adaptive_epsilon (xclk, 1.0e-5);
  nc_xcor_kernel_set_adaptive_boundary_tries (xclk, 3);
  nc_xcor_kernel_set_reltol (xclk, 1.0e-3);
  nc_xcor_kernel_set_scaled_abstol (xclk, 1.0e-4);
  nc_xcor_kernel_set_max_border_expansions (xclk, 1);
  nc_xcor_kernel_set_max_iter (xclk, 5);
  nc_xcor_kernel_set_expansion_factor (xclk, 0.5);
  nc_xcor_kernel_set_panel_order_cap (xclk, 12);
  nc_xcor_kernel_set_track_fit_residual (xclk, TRUE);

  g_assert_cmpuint (nc_xcor_kernel_get_lmax (xclk), ==, TEST_LMAX);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_get_adaptive_epsilon (xclk), ==, 1.0e-5, 1.0e-15, 0.0);
  g_assert_cmpuint (nc_xcor_kernel_get_adaptive_boundary_tries (xclk), ==, 3);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_get_reltol (xclk), ==, 1.0e-3, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_get_scaled_abstol (xclk), ==, 1.0e-4, 1.0e-15, 0.0);
  g_assert_cmpuint (nc_xcor_kernel_get_max_border_expansions (xclk), ==, 1);
  g_assert_cmpuint (nc_xcor_kernel_get_max_iter (xclk), ==, 5);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_get_expansion_factor (xclk), ==, 0.5, 1.0e-15, 0.0);
  g_assert_cmpuint (nc_xcor_kernel_get_panel_order_cap (xclk), ==, 12);
  g_assert_true (nc_xcor_kernel_get_track_fit_residual (xclk));
}

static void
test_nc_xcor_kernel_z_range (TestNcXcorKernel *test, gconstpointer pdata)
{
  const TestKernel *tk = pdata;
  gdouble zmin, zmax, zmid;

  nc_xcor_kernel_get_z_range (test->xclk, &zmin, &zmax, &zmid);

  g_assert_true (gsl_finite (zmin));
  g_assert_true (gsl_finite (zmax));
  g_assert_true (gsl_finite (zmid));
  g_assert_cmpfloat (zmin, >=, 0.0);
  g_assert_cmpfloat (zmin, <=, zmid);
  g_assert_cmpfloat (zmid, <=, zmax);
  g_assert_cmpfloat (zmax, <=, tk->zf);

  /* A CMB kernel integrates from here to last scattering; one that stopped inside the
   * photo-z range would be silently truncating most of its own line of sight. */
  if (tk->zf == TEST_ZMAX_LSS)
    g_assert_cmpfloat (zmax, >, 100.0);
}

static void
test_nc_xcor_kernel_limber (TestNcXcorKernel *test, gconstpointer pdata)
{
  const gint l = 10;
  gdouble zmin, zmax, zmid;
  guint j;

  nc_xcor_kernel_get_z_range (test->xclk, &zmin, &zmax, &zmid);

  g_assert_true (gsl_finite (nc_xcor_kernel_eval_limber_z_prefactor (test->xclk, test->cosmo, l)));

  for (j = 1; j < 8; j++)
  {
    const gdouble z = zmin + (zmax - zmin) * j / 8.0;
    const gdouble W = nc_xcor_kernel_eval_limber_z_full (test->xclk, test->cosmo, z, test->dist, l);
    NcXcorKinetic xck;

    g_assert_true (gsl_finite (W));

    /* The _full form looks the kinetic quantities up and folds in the prefactor; the
     * plain form takes them and does not. Composing them must reproduce it exactly --
     * and xi_z is in Hubble-radius units, which is the easy thing to get wrong. */
    xck.xi_z = nc_distance_comoving (test->dist, test->cosmo, z);
    xck.E_z  = nc_hicosmo_E (test->cosmo, z);

    ncm_assert_cmpdouble_e (nc_xcor_kernel_eval_limber_z (test->xclk, test->cosmo, z, &xck, l) *
                            nc_xcor_kernel_eval_limber_z_prefactor (test->xclk, test->cosmo, l),
                            ==, W, 1.0e-14, 0.0);
  }
}

static void
test_nc_xcor_kernel_k_range (TestNcXcorKernel *test, gconstpointer pdata)
{
  const gint ls[] = { 2, 10, 50 };
  guint j;

  for (j = 0; j < G_N_ELEMENTS (ls); j++)
  {
    gdouble kmin, kmax;

    nc_xcor_kernel_get_k_range (test->xclk, test->cosmo, ls[j], &kmin, &kmax);

    g_assert_true (gsl_finite (kmin));
    g_assert_true (gsl_finite (kmax));
    g_assert_cmpfloat (kmin, >, 0.0);
    g_assert_cmpfloat (kmin, <, kmax);
  }
}

static void
test_nc_xcor_kernel_components (TestNcXcorKernel *test, gconstpointer pdata)
{
  GPtrArray *comps = nc_xcor_kernel_get_component_list (test->xclk);
  guint i;

  g_assert_nonnull (comps);
  g_assert_cmpuint (comps->len, >, 0);

  for (i = 0; i < comps->len; i++)
  {
    NcXcorKernelComponent *comp = g_ptr_array_index (comps, i);
    gdouble xi_min, xi_max, k_min, k_max;

    g_assert_true (NC_IS_XCOR_KERNEL_COMPONENT (comp));

    nc_xcor_kernel_component_set_epsilon (comp, 1.0e-6);
    nc_xcor_kernel_component_set_ny (comp, 8);
    nc_xcor_kernel_component_set_max_iter (comp, 5);
    nc_xcor_kernel_component_set_tol (comp, 1.0e-4);

    ncm_assert_cmpdouble_e (nc_xcor_kernel_component_get_epsilon (comp), ==, 1.0e-6, 1.0e-15, 0.0);
    g_assert_cmpuint (nc_xcor_kernel_component_get_ny (comp), ==, 8);
    g_assert_cmpuint (nc_xcor_kernel_component_get_max_iter (comp), ==, 5);
    ncm_assert_cmpdouble_e (nc_xcor_kernel_component_get_tol (comp), ==, 1.0e-4, 1.0e-15, 0.0);

    /* The same four through the property interface, which serialization uses. */
    {
      NcXcorKernelComponent *comp2 = nc_xcor_kernel_component_ref (comp);
      gdouble epsilon, tol;
      guint ny, max_iter;

      g_object_get (comp, "epsilon", &epsilon, "ny", &ny,
                    "max-iter", &max_iter, "tol", &tol, NULL);

      ncm_assert_cmpdouble_e (epsilon, ==, 1.0e-6, 1.0e-15, 0.0);
      g_assert_cmpuint (ny, ==, 8);
      g_assert_cmpuint (max_iter, ==, 5);
      ncm_assert_cmpdouble_e (tol, ==, 1.0e-4, 1.0e-15, 0.0);

      nc_xcor_kernel_component_clear (&comp2);
      g_assert_true (comp2 == NULL);
    }

    nc_xcor_kernel_component_prepare (comp, test->cosmo);
    nc_xcor_kernel_component_get_limits (comp, test->cosmo, &xi_min, &xi_max, &k_min, &k_max);

    g_assert_true (gsl_finite (xi_min));
    g_assert_true (gsl_finite (xi_max));
    g_assert_cmpfloat (xi_min, <, xi_max);
    g_assert_cmpfloat (k_min, >, 0.0);
    g_assert_cmpfloat (k_min, <, k_max);

    g_assert_true (gsl_finite (nc_xcor_kernel_component_eval_kernel (comp, test->cosmo,
                                                                     0.5 * (xi_min + xi_max),
                                                                     sqrt (k_min * k_max))));
    g_assert_true (gsl_finite (nc_xcor_kernel_component_eval_prefactor (comp, test->cosmo,
                                                                        sqrt (k_min * k_max), 10)));

    /* The k-limit searches: bounded by max_iter above, so they take their branches
     * without being run to convergence. */
    g_assert_true (gsl_finite (nc_xcor_kernel_component_eval_k_max (comp, 10.0)));
    g_assert_true (gsl_finite (nc_xcor_kernel_component_eval_KL_max (comp, 10.0)));
    g_assert_true (gsl_finite (nc_xcor_kernel_component_eval_k_epsilon (comp, 10.0)));
  }

  g_ptr_array_unref (comps);
}

static void
test_nc_xcor_kernel_noise (TestNcXcorKernel *test, gconstpointer pdata)
{
  const TestKernel *tk = pdata;
  const guint len      = 16;
  NcmVector *v1        = ncm_vector_new (len);
  NcmVector *v2        = ncm_vector_new (len);
  guint i;

  ncm_vector_set_all (v1, 1.0e-9);
  ncm_vector_set_all (v2, 1.0e-9);

  nc_xcor_kernel_add_noise (test->xclk, v1, v2, 2);

  for (i = 0; i < len; i++)
  {
    g_assert_true (gsl_finite (ncm_vector_get (v2, i)));
    g_assert_cmpfloat (ncm_vector_get (v2, i), >=, ncm_vector_get (v1, i));
  }

  /* A kernel that declares noise must actually add some, or the auto-spectrum it
   * contributes to is silently noiseless. */
  if (tk->has_noise)
    g_assert_cmpfloat (ncm_vector_get (v2, len - 1), >, ncm_vector_get (v1, len - 1));

  ncm_vector_free (v1);
  ncm_vector_free (v2);
}

static void
test_nc_xcor_kernel_serialize (TestNcXcorKernel *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcXcorKernel *dup = NC_XCOR_KERNEL (ncm_serialize_dup_obj (ser, G_OBJECT (test->xclk)));
  gdouble zmin, zmax, zmid;
  guint j;

  g_assert_true (G_OBJECT_TYPE (dup) == G_OBJECT_TYPE (test->xclk));
  g_assert_cmpuint (nc_xcor_kernel_obs_len (dup), ==, nc_xcor_kernel_obs_len (test->xclk));

  nc_xcor_kernel_prepare (dup, test->cosmo);
  nc_xcor_kernel_get_z_range (test->xclk, &zmin, &zmax, &zmid);

  for (j = 1; j < 8; j++)
  {
    const gdouble z = zmin + (zmax - zmin) * j / 8.0;

    ncm_assert_cmpdouble_e (nc_xcor_kernel_eval_limber_z_full (dup, test->cosmo, z, test->dist, 10),
                            ==,
                            nc_xcor_kernel_eval_limber_z_full (test->xclk, test->cosmo, z, test->dist, 10),
                            1.0e-9, 0.0);
  }

  nc_xcor_kernel_free (dup);
  ncm_serialize_free (ser);
}

/* The _new() helpers leave the integrator unset, which pins the kernel to the Limber
 * tier -- nc_xcor_kernel_set_l_limber() refuses to move it off. The property form is how
 * a caller opts in, and it is the only path to the non-Limber knobs. */
static void
test_nc_xcor_kernel_integrator (void)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps   = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                            NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, 8));
  NcmSpline *dndz           = _dndz_new ();
  NcXcorKernel *xclk;

  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));

  xclk = NC_XCOR_KERNEL (g_object_new (NC_TYPE_XCOR_KERNEL_GAL,
                                       "dist", dist,
                                       "powspec", ps,
                                       "dndz", dndz,
                                       "integrator", sbi,
                                       NULL));

  g_assert_true (nc_xcor_kernel_peek_integrator (xclk) == sbi);

  nc_xcor_kernel_set_l_limber (xclk, 7);
  g_assert_cmpint (nc_xcor_kernel_get_l_limber (xclk), ==, 7);

  nc_xcor_kernel_free (xclk);
  ncm_spline_free (dndz);
  ncm_sbessel_integrator_free (sbi);
  ncm_powspec_free (ps);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
}

/* The galaxy kernel caches its bias so a likelihood stepping only b(z) can skip
 * rebuilding the window. The cache is what set_bias_old()/get_bias() expose, and getting
 * it wrong means either a stale window or no speedup at all. */
static void
test_nc_xcor_kernel_gal_bias (void)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps   = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                            NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcmSpline *dndz    = _dndz_new ();
  NcXcorKernelGal *g = nc_xcor_kernel_gal_new (dist, ps, 1, 1.234, dndz, FALSE);
  gdouble bias, bias_old, noise_bias_old;

  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));

  ncm_model_orig_vparam_set (NCM_MODEL (g), NC_XCOR_KERNEL_GAL_BIAS, 0, 1.5);
  nc_xcor_kernel_prepare (NC_XCOR_KERNEL (g), cosmo);

  nc_xcor_kernel_gal_set_fast_update (g, TRUE);
  g_assert_true (nc_xcor_kernel_gal_get_fast_update (g));
  nc_xcor_kernel_gal_set_fast_update (g, FALSE);
  g_assert_false (nc_xcor_kernel_gal_get_fast_update (g));

  nc_xcor_kernel_gal_set_bias_old (g, 1.5, 1.234);
  nc_xcor_kernel_gal_get_bias (g, &bias, &bias_old, &noise_bias_old);

  ncm_assert_cmpdouble_e (bias, ==, 1.5, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (bias_old, ==, 1.5, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (noise_bias_old, ==, 1.234, 1.0e-15, 0.0);

  nc_xcor_kernel_free (NC_XCOR_KERNEL (g));
  ncm_spline_free (dndz);
  ncm_powspec_free (ps);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
}

/* nc_xcor_kernel_table_new() takes the density kind, cubic order and normalization as
 * given; _new_full() is where a caller chooses them, and the shear kind is a different
 * radial factor rather than a different scaling. */
static void
test_nc_xcor_kernel_table_full (void)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps   = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                            NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (0, 8));
  NcmVector *chi            = ncm_vector_new (64);
  NcmVector *W              = ncm_vector_new (64);
  NcXcorKernelTable *xckt;
  guint i;

  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));

  for (i = 0; i < 64; i++)
  {
    const gdouble chi_i = 1000.0 + 20.0 * i;
    const gdouble t     = (chi_i - 1630.0) / 300.0;

    ncm_vector_set (chi, i, chi_i);
    ncm_vector_set (W, i, exp (-0.5 * t * t));
  }

  xckt = nc_xcor_kernel_table_new_full (dist, ps, chi, W,
                                        NC_XCOR_KERNEL_TABLE_KIND_SHEAR, 6, FALSE, sbi);

  g_assert_cmpuint (nc_xcor_kernel_table_get_kind (xckt), ==, NC_XCOR_KERNEL_TABLE_KIND_SHEAR);
  g_assert_cmpuint (nc_xcor_kernel_table_get_order (xckt), ==, 6);
  g_assert_false (nc_xcor_kernel_table_get_normalize (xckt));
  g_assert_true (gsl_finite (nc_xcor_kernel_table_get_norm (xckt)));
  g_assert_nonnull (nc_xcor_kernel_table_peek_spline (xckt));
  g_assert_nonnull (nc_xcor_kernel_table_peek_knots (xckt));

  nc_xcor_kernel_prepare (NC_XCOR_KERNEL (xckt), cosmo);
  g_assert_true (gsl_finite (nc_xcor_kernel_radial_eval_W (NC_XCOR_KERNEL_RADIAL (xckt), 1600.0)));
  g_assert_true (gsl_finite (nc_xcor_kernel_eval_limber_z_full (NC_XCOR_KERNEL (xckt), cosmo, 0.5, dist, 10)));

  nc_xcor_kernel_free (NC_XCOR_KERNEL (xckt));
  ncm_vector_free (chi);
  ncm_vector_free (W);
  ncm_sbessel_integrator_free (sbi);
  ncm_powspec_free (ps);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
}

/* An optional scale-dependent factor multiplying the radial integrand. Unset it is the
 * identity, which is why a radial kernel never has to know whether one is attached. */
static void
test_nc_xcor_kernel_radial_kdep (void)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps   = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                            NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcXcorKernelRadialKDepGrowth *kdepg = nc_xcor_kernel_radial_kdep_growth_new (0.3, 0.05, 1500.0);
  NcXcorKernelRadialKDep *kdep        = NC_XCOR_KERNEL_RADIAL_KDEP (kdepg);
  NcXcorKernel *plain;
  NcXcorKernel *withk;
  gdouble amplitude, k_transition, chi_ref;

  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));

  ncm_assert_cmpdouble_e (nc_xcor_kernel_radial_kdep_growth_get_amplitude (kdepg), ==, 0.3, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_radial_kdep_growth_get_k_transition (kdepg), ==, 0.05, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_radial_kdep_growth_get_chi_ref (kdepg), ==, 1500.0, 1.0e-15, 0.0);

  g_object_get (kdepg, "amplitude", &amplitude, "k-transition", &k_transition,
                "chi-ref", &chi_ref, NULL);
  ncm_assert_cmpdouble_e (amplitude, ==, 0.3, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (k_transition, ==, 0.05, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (chi_ref, ==, 1500.0, 1.0e-15, 0.0);

  g_assert_true (gsl_finite (nc_xcor_kernel_radial_kdep_eval (kdep, 1500.0, 0.05)));

  /* At the reference distance the growth factor is normalized away, whatever k. */
  ncm_assert_cmpdouble_e (nc_xcor_kernel_radial_kdep_eval (kdep, 1500.0, 0.5), ==, 1.0, 1.0e-14, 0.0);

  plain = NC_XCOR_KERNEL (nc_xcor_kernel_analytic_gauss_new (dist, ps, 1500.0, 300.0, 4.0));
  withk = NC_XCOR_KERNEL (g_object_new (NC_TYPE_XCOR_KERNEL_ANALYTIC_GAUSS,
                                        "dist", dist, "powspec", ps,
                                        "chi-mean", 1500.0, "chi-sigma", 300.0, "n-sigma", 4.0,
                                        "scale-dependence", kdep,
                                        NULL));

  g_assert_null (nc_xcor_kernel_radial_peek_kdep (NC_XCOR_KERNEL_RADIAL (plain)));
  g_assert_true (nc_xcor_kernel_radial_peek_kdep (NC_XCOR_KERNEL_RADIAL (withk)) == kdep);

  nc_xcor_kernel_prepare (plain, cosmo);
  nc_xcor_kernel_prepare (withk, cosmo);

  /* The factor multiplies the k-dependent integrand, not the window itself. */
  ncm_assert_cmpdouble_e (nc_xcor_kernel_radial_eval_W (NC_XCOR_KERNEL_RADIAL (withk), 1600.0), ==,
                          nc_xcor_kernel_radial_eval_W (NC_XCOR_KERNEL_RADIAL (plain), 1600.0),
                          1.0e-14, 0.0);

  g_assert_true (gsl_finite (nc_xcor_kernel_radial_eval_kernel_factor (NC_XCOR_KERNEL_RADIAL (withk),
                                                                       cosmo, 1600.0, 0.5)));

  nc_xcor_kernel_free (plain);
  nc_xcor_kernel_free (withk);
  nc_xcor_kernel_radial_kdep_free (kdep);
  ncm_powspec_free (ps);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
}

typedef struct _TestCheck
{
  const gchar *name;

  void (*func) (TestNcXcorKernel *test, gconstpointer pdata);
} TestCheck;

static const TestCheck test_checks[] = {
  {"basic",      &test_nc_xcor_kernel_basic},
  {"knobs",      &test_nc_xcor_kernel_knobs},
  {"z_range",    &test_nc_xcor_kernel_z_range},
  {"limber",     &test_nc_xcor_kernel_limber},
  {"k_range",    &test_nc_xcor_kernel_k_range},
  {"components", &test_nc_xcor_kernel_components},
  {"noise",      &test_nc_xcor_kernel_noise},
  {"serialize",  &test_nc_xcor_kernel_serialize},
};

gint
main (gint argc, gchar *argv[])
{
  guint i, j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < G_N_ELEMENTS (test_kernels); i++)
  {
    for (j = 0; j < G_N_ELEMENTS (test_checks); j++)
    {
      gchar *path = g_strdup_printf ("/nc/xcor/kernel/%s/%s", test_kernels[i].name, test_checks[j].name);

      g_test_add (path, TestNcXcorKernel, &test_kernels[i],
                  &test_nc_xcor_kernel_new,
                  test_checks[j].func,
                  &test_nc_xcor_kernel_free);

      g_free (path);
    }
  }

  g_test_add_func ("/nc/xcor/kernel/integrator", &test_nc_xcor_kernel_integrator);
  g_test_add_func ("/nc/xcor/kernel/gal/bias", &test_nc_xcor_kernel_gal_bias);
  g_test_add_func ("/nc/xcor/kernel/table/new_full", &test_nc_xcor_kernel_table_full);
  g_test_add_func ("/nc/xcor/kernel/radial/kdep", &test_nc_xcor_kernel_radial_kdep);

  g_test_run ();
}

