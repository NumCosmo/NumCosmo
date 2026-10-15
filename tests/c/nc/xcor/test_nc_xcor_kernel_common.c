/***************************************************************************
 *            test_nc_xcor_kernel_common.c
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
#include "test_nc_xcor_kernel_common.h"

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_integration.h>

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
 * @own_binary: whether the kernel's checks run in a test binary of their own
 *
 * The CMB kernels' checks cost about 30 s each in the coverage build, nine tenths of the
 * file, so each runs in its own binary and the rest share test_nc_xcor_kernel.c.
 */
typedef struct _TestKernel
{
  const gchar *name;

  NcXcorKernel *(*build) (TestNcXcorKernel *test);

  gdouble zf;
  gboolean has_noise;
  gboolean own_binary;
} TestKernel;

/* The crosses the backends are compared on, by the name of the fixture's kernel */
typedef struct _TestCross
{
  const gchar *name;

  NcXcorKernel *(*build_partner) (TestNcXcorKernel *test);
} TestCross;

typedef NcXcorKernel *(*TestDndzBuild) (NcDistance *dist, NcmPowspec *ps, NcmSpline *dndz);

typedef struct _TestCheck
{
  const gchar *name;

  void (*func) (TestNcXcorKernel *test, gconstpointer pdata);
} TestCheck;

static void test_nc_xcor_kernel_new (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_free (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_basic (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_knobs (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_z_range (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_limber (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_k_range (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_components (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_noise (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_serialize (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_limber_backends (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_limber_backends_cross (TestNcXcorKernel *test, gconstpointer pdata);
static void test_nc_xcor_kernel_dndz_untouched (void);
static void test_nc_xcor_kernel_gal_magbias_null (void);
static void test_nc_xcor_kernel_integrator (void);
static void test_nc_xcor_kernel_cluster_limber_range_ends (void);
static void test_nc_xcor_kernel_gal_limber_capped_powspec (void);
static void test_nc_xcor_kernel_gal_rsd_limber_steps (void);
static void test_nc_xcor_kernel_gal_rsd_limber_z_st (void);
static void test_nc_xcor_kernel_gal_rsd_limber_z_trap (void);
static void test_nc_xcor_kernel_gal_bias (void);
static void test_nc_xcor_kernel_cluster_tophat_volume_normalized (void);
static void test_nc_xcor_kernel_table_full (void);
static void test_nc_xcor_kernel_radial_kdep (void);
static NcXcorKernel *_build_gal (TestNcXcorKernel *test);
static NcXcorKernel *_build_gal_magbias (TestNcXcorKernel *test);
static NcXcorKernel *_build_weak_lensing (TestNcXcorKernel *test);
static NcXcorKernel *_build_cmb_lensing (TestNcXcorKernel *test);
static NcXcorKernel *_build_cmb_isw (TestNcXcorKernel *test);
static NcXcorKernel *_build_tsz (TestNcXcorKernel *test);
static NcXcorKernel *_build_cluster_tophat (TestNcXcorKernel *test);
static NcXcorKernel *_build_table (TestNcXcorKernel *test);

static const TestKernel test_kernels[] = {
  {"gal",            &_build_gal,            TEST_ZMAX,     TRUE,  FALSE},
  {"gal_magbias",    &_build_gal_magbias,    TEST_ZMAX,     TRUE,  FALSE},
  {"weak_lensing",   &_build_weak_lensing,   TEST_ZMAX,     TRUE,  FALSE},
  {"cmb_lensing",    &_build_cmb_lensing,    TEST_ZMAX_LSS, TRUE,  TRUE },
  {"cmb_isw",        &_build_cmb_isw,        TEST_ZMAX_LSS, TRUE,  TRUE },
  {"tsz",            &_build_tsz,            TEST_ZMAX,     FALSE, FALSE},
  {"cluster_tophat", &_build_cluster_tophat, TEST_ZMAX,     FALSE, FALSE},
  {"table",          &_build_table,          TEST_ZMAX,     FALSE, FALSE},
};

static const TestCross test_crosses[] = {
  {"gal",         &_build_weak_lensing},
  {"cmb_lensing", &_build_cmb_isw     },
};

static const TestCheck test_checks[] = {
  {"basic",      &test_nc_xcor_kernel_basic},
  {"knobs",      &test_nc_xcor_kernel_knobs},
  {"z_range",    &test_nc_xcor_kernel_z_range},
  {"limber",     &test_nc_xcor_kernel_limber},
  {"k_range",    &test_nc_xcor_kernel_k_range},
  {"components", &test_nc_xcor_kernel_components},
  {"noise",      &test_nc_xcor_kernel_noise},
  {"serialize",  &test_nc_xcor_kernel_serialize},
  {"limber_backends", &test_nc_xcor_kernel_limber_backends},
};

/**
 * test_nc_xcor_kernel_main:
 * @argc: the argument count of main ()
 * @argv: the arguments of main ()
 * @kernel: (nullable): the kernel this binary runs, one marked own_binary in test_kernels
 *
 * Registers and runs the checks of @kernel, or with %NULL those of every kernel not marked
 * own_binary together with the kernel-specific tests.
 *
 * Returns: the result of g_test_run ()
 */
gint
test_nc_xcor_kernel_main (gint argc, gchar *argv[], const gchar *kernel)
{
  gboolean found = (kernel == NULL);
  guint i, j;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  for (i = 0; i < G_N_ELEMENTS (test_kernels); i++)
  {
    if (kernel == NULL)
    {
      if (test_kernels[i].own_binary)
        continue;
    }
    else if (g_strcmp0 (kernel, test_kernels[i].name) != 0)
    {
      continue;
    }
    else if (!test_kernels[i].own_binary)
    {
      g_error ("test_nc_xcor_kernel_main: kernel `%s' is not marked own_binary.", kernel);
    }
    else
    {
      found = TRUE;
    }

    for (j = 0; j < G_N_ELEMENTS (test_checks); j++)
    {
      gchar *path = g_strdup_printf ("/nc/xcor/kernel/%s/%s", test_kernels[i].name, test_checks[j].name);

      g_test_add (path, TestNcXcorKernel, &test_kernels[i],
                  &test_nc_xcor_kernel_new,
                  test_checks[j].func,
                  &test_nc_xcor_kernel_free);

      g_free (path);
    }

    for (j = 0; j < G_N_ELEMENTS (test_crosses); j++)
    {
      if (g_strcmp0 (test_kernels[i].name, test_crosses[j].name) == 0)
      {
        gchar *path = g_strdup_printf ("/nc/xcor/kernel/%s/limber_backends_cross", test_kernels[i].name);

        g_test_add (path, TestNcXcorKernel, &test_kernels[i],
                    &test_nc_xcor_kernel_new,
                    &test_nc_xcor_kernel_limber_backends_cross,
                    &test_nc_xcor_kernel_free);

        g_free (path);
      }
    }
  }

  if (!found)
    g_error ("test_nc_xcor_kernel_main: no kernel `%s' in the table.", kernel);

  if (kernel != NULL)
    return g_test_run ();

  g_test_add_func ("/nc/xcor/kernel/integrator", &test_nc_xcor_kernel_integrator);
  g_test_add_func ("/nc/xcor/kernel/dndz_untouched", &test_nc_xcor_kernel_dndz_untouched);
  g_test_add_func ("/nc/xcor/kernel/gal/magbias_null", &test_nc_xcor_kernel_gal_magbias_null);
  g_test_add_func ("/nc/xcor/kernel/gal/bias", &test_nc_xcor_kernel_gal_bias);
  g_test_add_func ("/nc/xcor/kernel/cluster_tophat/volume_normalized", &test_nc_xcor_kernel_cluster_tophat_volume_normalized);
  g_test_add_func ("/nc/xcor/kernel/gal/rsd_limber_steps", &test_nc_xcor_kernel_gal_rsd_limber_steps);
  g_test_add_func ("/nc/xcor/kernel/gal/limber_capped_powspec", &test_nc_xcor_kernel_gal_limber_capped_powspec);
  g_test_add_func ("/nc/xcor/kernel/cluster/limber_range_ends", &test_nc_xcor_kernel_cluster_limber_range_ends);
  g_test_add_func ("/nc/xcor/kernel/gal/rsd_limber_z/subprocess", &test_nc_xcor_kernel_gal_rsd_limber_z_st);
  g_test_add_func ("/nc/xcor/kernel/gal/rsd_limber_z/trap", &test_nc_xcor_kernel_gal_rsd_limber_z_trap);
  g_test_add_func ("/nc/xcor/kernel/table/new_full", &test_nc_xcor_kernel_table_full);
  g_test_add_func ("/nc/xcor/kernel/radial/kdep", &test_nc_xcor_kernel_radial_kdep);

  return g_test_run ();
}

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

/* The Limber window of @xclk at eight redshifts from the lower end of its range and the
 * multipoles in @ells, z outermost */
static GArray *
_test_limber_values (TestNcXcorKernel *test, NcXcorKernel *xclk, const gint *ells, guint n_ell)
{
  GArray *values = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble zmin, zmax, zmid;
  guint i, j;

  nc_xcor_kernel_get_z_range (test->xclk, &zmin, &zmax, &zmid);

  for (j = 0; j < 8; j++)
  {
    const gdouble z = zmin + (zmax - zmin) * j / 8.0;

    for (i = 0; i < n_ell; i++)
    {
      const gdouble W = nc_xcor_kernel_eval_limber_z_full (xclk, test->cosmo, z, test->dist, ells[i]);

      g_array_append_val (values, W);
    }
  }

  return values;
}

static void
_test_assert_same_values (GArray *a, GArray *b)
{
  guint i;

  g_assert_cmpuint (a->len, ==, b->len);

  for (i = 0; i < a->len; i++)
    g_assert_cmpfloat (g_array_index (b, gdouble, i), ==, g_array_index (a, gdouble, i));
}

/* The two Limber backends integrate the same z integrand with different quadratures, so
 * on ell 2 to 100 they agree to 1e-5 with GSL asked for 1e-7. An auto spectrum (k1 == k2)
 * integrates W^2 P and is positive at every multipole. */
static void
_test_limber_backends_pair (TestNcXcorKernel *test, NcXcorKernel *k1, NcXcorKernel *k2)
{
  const guint lmin  = 2;
  const guint lmax  = 100;
  NcmVector *cl_gsl = ncm_vector_new (lmax - lmin + 1);
  NcmVector *cl_cub = ncm_vector_new (lmax - lmin + 1);
  NcXcor *xc_gsl    = nc_xcor_new (test->dist, test->ps, NC_XCOR_METHOD_LIMBER_Z_GSL);
  NcXcor *xc_cub    = nc_xcor_new (test->dist, test->ps, NC_XCOR_METHOD_LIMBER_Z_CUBATURE);
  guint i;

  nc_xcor_set_reltol (xc_gsl, 1.0e-7);
  nc_xcor_prepare (xc_gsl, test->cosmo);
  nc_xcor_prepare (xc_cub, test->cosmo);

  nc_xcor_compute (xc_gsl, k1, k2, test->cosmo, lmin, lmax, cl_gsl);
  nc_xcor_compute (xc_cub, k1, k2, test->cosmo, lmin, lmax, cl_cub);

  /* A spectrum that is zero everywhere would agree trivially. */
  {
    gdouble absmin, absmax;

    ncm_vector_get_absminmax (cl_gsl, &absmin, &absmax);
    g_assert_cmpfloat (absmax, >, 0.0);
  }

  for (i = 0; i < lmax - lmin + 1; i++)
  {
    ncm_assert_cmpdouble_e (ncm_vector_get (cl_cub, i), ==, ncm_vector_get (cl_gsl, i), 1.0e-5, 1.0e-50);

    if (k1 == k2)
    {
      g_assert_cmpfloat (ncm_vector_get (cl_gsl, i), >, 0.0);
      g_assert_cmpfloat (ncm_vector_get (cl_cub, i), >, 0.0);
    }
  }

  nc_xcor_free (xc_gsl);
  nc_xcor_free (xc_cub);
  ncm_vector_free (cl_gsl);
  ncm_vector_free (cl_cub);
}

/* A normalized Gaussian dN/dz like _dndz_new (), centred at @mean */
static NcmSpline *
_dndz_at_new (gdouble mean)
{
  NcmVector *z = ncm_vector_new (TEST_DNDZ_KNOTS);
  NcmVector *n = ncm_vector_new (TEST_DNDZ_KNOTS);
  NcmSpline *s;
  guint i;

  for (i = 0; i < TEST_DNDZ_KNOTS; i++)
  {
    const gdouble zi = 2.0 * i / (gdouble) (TEST_DNDZ_KNOTS - 1);
    const gdouble t  = (zi - mean) / TEST_DNDZ_SIGMA;

    ncm_vector_set (z, i, zi);
    ncm_vector_set (n, i, exp (-0.5 * t * t));
  }

  s = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (z, n, TRUE));

  ncm_vector_free (z);
  ncm_vector_free (n);

  return s;
}

static NcXcorKernel *
_dndz_build_weak_lensing (NcDistance *dist, NcmPowspec *ps, NcmSpline *dndz)
{
  return NC_XCOR_KERNEL (nc_xcor_kernel_weak_lensing_new (dist, ps, dndz, 3.0, 0.3));
}

static NcXcorKernel *
_dndz_build_gal (NcDistance *dist, NcmPowspec *ps, NcmSpline *dndz)
{
  NcXcorKernelGal *g = nc_xcor_kernel_gal_new (dist, ps, 1, 1.234, dndz, TRUE);

  ncm_model_orig_vparam_set (NCM_MODEL (g), NC_XCOR_KERNEL_GAL_BIAS, 0, 1.5);

  return NC_XCOR_KERNEL (g);
}

/* The Limber window of @xclk at six redshifts and two multipoles */
static void
_dndz_values (NcXcorKernel *xclk, NcHICosmo *cosmo, NcDistance *dist, gdouble *out)
{
  const gint ells[2] = {10, 200};
  guint j, i;

  for (j = 0; j < 6; j++)
    for (i = 0; i < 2; i++)
      out[2 * j + i] = nc_xcor_kernel_eval_limber_z_full (xclk, cosmo, 0.1 + 0.3 * j, dist, ells[i]);
}

/* A Gaussian dN/dz truncated on [0.4, 0.9], nonzero at both ends as a photometric bin is. */
static NcmSpline *
_dndz_truncated_new (void)
{
  NcmVector *z = ncm_vector_new (TEST_DNDZ_KNOTS);
  NcmVector *n = ncm_vector_new (TEST_DNDZ_KNOTS);
  NcmSpline *s;
  guint i;

  for (i = 0; i < TEST_DNDZ_KNOTS; i++)
  {
    const gdouble zi = 0.4 + 0.5 * i / (gdouble) (TEST_DNDZ_KNOTS - 1);
    const gdouble t  = (zi - TEST_DNDZ_MEAN) / TEST_DNDZ_SIGMA;

    ncm_vector_set (z, i, zi);
    ncm_vector_set (n, i, exp (-0.5 * t * t));
  }

  s = NCM_SPLINE (ncm_spline_cubic_notaknot_new_full (z, n, TRUE));

  ncm_vector_free (z);
  ncm_vector_free (n);

  return s;
}

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
  NcXcorKernel *xclk                = test->xclk;
  NcmSpectralPanelsBuilder *builder = ncm_spectral_panels_builder_new ();

  nc_xcor_kernel_set_lmax (xclk, TEST_LMAX);
  nc_xcor_kernel_set_adaptive_epsilon (xclk, 1.0e-5);
  nc_xcor_kernel_set_adaptive_boundary_tries (xclk, 3);
  nc_xcor_kernel_set_reltol (xclk, 1.0e-3);
  nc_xcor_kernel_set_peak_epsilon (xclk, 1.0e-4);
  nc_xcor_kernel_set_max_border_expansions (xclk, 1);
  nc_xcor_kernel_set_expansion_factor (xclk, 0.5);
  nc_xcor_kernel_set_panels_builder (xclk, builder);

  g_assert_cmpuint (nc_xcor_kernel_get_lmax (xclk), ==, TEST_LMAX);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_get_adaptive_epsilon (xclk), ==, 1.0e-5, 1.0e-15, 0.0);
  g_assert_cmpuint (nc_xcor_kernel_get_adaptive_boundary_tries (xclk), ==, 3);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_get_reltol (xclk), ==, 1.0e-3, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_get_peak_epsilon (xclk), ==, 1.0e-4, 1.0e-15, 0.0);
  g_assert_cmpuint (nc_xcor_kernel_get_max_border_expansions (xclk), ==, 1);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_get_expansion_factor (xclk), ==, 0.5, 1.0e-15, 0.0);
  g_assert_true (nc_xcor_kernel_peek_panels_builder (xclk) == builder);

  /* The GObject properties are the same knobs: serialization reads and writes them, so a
   * property and its accessor have to agree in both directions. */
  {
    gdouble reltol, peak_epsilon, expansion_factor;
    guint max_border_expansions;

    g_object_get (xclk,
                  "reltol", &reltol,
                  "peak-epsilon", &peak_epsilon,
                  "max-border-expansions", &max_border_expansions,
                  "expansion-factor", &expansion_factor,
                  NULL);

    ncm_assert_cmpdouble_e (reltol, ==, 1.0e-3, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (peak_epsilon, ==, 1.0e-4, 1.0e-15, 0.0);
    g_assert_cmpuint (max_border_expansions, ==, 1);
    ncm_assert_cmpdouble_e (expansion_factor, ==, 0.5, 1.0e-15, 0.0);

    g_object_set (xclk,
                  "reltol", 2.0e-3,
                  "peak-epsilon", 2.0e-4,
                  "max-border-expansions", 2,
                  "expansion-factor", 0.25,
                  NULL);

    ncm_assert_cmpdouble_e (nc_xcor_kernel_get_reltol (xclk), ==, 2.0e-3, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (nc_xcor_kernel_get_peak_epsilon (xclk), ==, 2.0e-4, 1.0e-15, 0.0);
    g_assert_cmpuint (nc_xcor_kernel_get_max_border_expansions (xclk), ==, 2);
    ncm_assert_cmpdouble_e (nc_xcor_kernel_get_expansion_factor (xclk), ==, 0.25, 1.0e-15, 0.0);
  }

  /* Below its useful floor peak-epsilon is accepted and warned about, never clamped: it
   * enters the C_l integrand squared, so 1e-8 there is 1e-16. At the floor there is no
   * warning, and an unexpected one would fail this test. */
  g_test_expect_message (G_LOG_DOMAIN, G_LOG_LEVEL_WARNING, "*below the useful floor*");
  nc_xcor_kernel_set_peak_epsilon (xclk, 1.0e-8);
  g_test_assert_expected_messages ();
  ncm_assert_cmpdouble_e (nc_xcor_kernel_get_peak_epsilon (xclk), ==, 1.0e-8, 1.0e-15, 0.0);

  nc_xcor_kernel_set_peak_epsilon (xclk, NC_XCOR_KERNEL_MIN_USEFUL_PEAK_EPSILON);
  ncm_assert_cmpdouble_e (nc_xcor_kernel_get_peak_epsilon (xclk), ==, NC_XCOR_KERNEL_MIN_USEFUL_PEAK_EPSILON, 1.0e-15, 0.0);

  ncm_spectral_panels_builder_free (builder);
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

  /* The lower end of the range is a point the integrators may evaluate. For the lensing
   * kernels it is z = 0, where chi and the transverse distance vanish and the window is
   * its limit 0, not 0 * inf or 0 / 0. */
  {
    const gdouble W_min = nc_xcor_kernel_eval_limber_z_full (test->xclk, test->cosmo, zmin, test->dist, l);

    g_assert_true (gsl_finite (W_min));

    if (NC_IS_XCOR_KERNEL_WEAK_LENSING (test->xclk) || NC_IS_XCOR_KERNEL_CMB_LENSING (test->xclk))
    {
      g_assert_cmpfloat (zmin, ==, 0.0);
      g_assert_cmpfloat (W_min, ==, 0.0);
    }
  }

  for (j = 1; j < 8; j++)
  {
    const gdouble z = zmin + (zmax - zmin) * j / 8.0;
    const gdouble W = nc_xcor_kernel_eval_limber_z_full (test->xclk, test->cosmo, z, test->dist, l);
    NcXcorKinetic xck;

    g_assert_true (gsl_finite (W));

    /* The _full form looks the kinetic quantities up and folds in the prefactor; the
     * plain form takes them and does not. Composing them must reproduce it exactly --
     * and chi_z is in Hubble-radius units, which is the easy thing to get wrong. */
    xck.chi_z = nc_distance_comoving (test->dist, test->cosmo, z);
    xck.E_z   = nc_hicosmo_E (test->cosmo, z);
    xck.z     = z;

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

    /* The defaults are usable before anything is set */
    g_assert_true (gsl_finite (nc_xcor_kernel_component_get_epsilon (comp)));
    g_assert_cmpfloat (nc_xcor_kernel_component_get_epsilon (comp), >, 0.0);
    g_assert_cmpuint (nc_xcor_kernel_component_get_ny (comp), >, 0);
    g_assert_cmpuint (nc_xcor_kernel_component_get_max_iter (comp), >, 0);
    g_assert_true (gsl_finite (nc_xcor_kernel_component_get_tol (comp)));
    g_assert_cmpfloat (nc_xcor_kernel_component_get_tol (comp), >, 0.0);

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
    g_assert_true (gsl_finite (k_max));
    g_assert_cmpfloat (xi_min, >, 0.0);
    g_assert_cmpfloat (xi_min, <, xi_max);
    g_assert_cmpfloat (k_min, >, 0.0);
    g_assert_cmpfloat (k_min, <, k_max);

    g_assert_true (gsl_finite (nc_xcor_kernel_component_eval_kernel (comp, test->cosmo,
                                                                     0.5 * (xi_min + xi_max),
                                                                     sqrt (k_min * k_max))));
    g_assert_true (gsl_finite (nc_xcor_kernel_component_eval_prefactor (comp, test->cosmo,
                                                                        sqrt (k_min * k_max), 10)));
    g_assert_true (gsl_finite (nc_xcor_kernel_component_eval_prefactor (comp, test->cosmo,
                                                                        sqrt (k_min * k_max), 100)));

    /* The kernel factorizes into the point factor and the wave number factor. */
    {
      const guint n_chi = 5;
      const guint n_k   = 5;
      guint a, b;

      for (a = 0; a < n_chi; a++)
      {
        const gdouble chi = xi_min + (xi_max - xi_min) * (a + 0.5) / n_chi;
        NcXcorKinetic xck;

        xck.chi_z = chi;
        xck.z     = nc_distance_inv_comoving (test->dist, test->cosmo, chi);
        xck.E_z   = nc_hicosmo_E (test->cosmo, xck.z);

        for (b = 0; b < n_k; b++)
        {
          const gdouble k      = k_min * exp (log (k_max / k_min) * (b + 0.5) / n_k);
          const gdouble kernel = nc_xcor_kernel_component_eval_kernel (comp, test->cosmo, chi, k);
          const gdouble split  = nc_xcor_kernel_component_eval_window (comp, test->cosmo, &xck) *
                                 nc_xcor_kernel_component_eval_kfactor (comp, test->cosmo, &xck, k);

          ncm_assert_cmpdouble_e (split, ==, kernel, 1.0e-13, 0.0);
        }

        /* The vector form reproduces the scalar factor at every wave number. */
        {
          NcmVector *kv  = ncm_vector_new (n_k);
          NcmVector *out = ncm_vector_new (n_k);

          for (b = 0; b < n_k; b++)
            ncm_vector_set (kv, b, k_min * exp (log (k_max / k_min) * (b + 0.5) / n_k));

          nc_xcor_kernel_component_eval_kfactor_vec (comp, test->cosmo, &xck, kv, out);

          for (b = 0; b < n_k; b++)
            ncm_assert_cmpdouble_e (ncm_vector_get (out, b), ==,
                                    nc_xcor_kernel_component_eval_kfactor (comp, test->cosmo, &xck, ncm_vector_get (kv, b)),
                                    1.0e-13, 0.0);

          ncm_vector_free (kv);
          ncm_vector_free (out);
        }
      }
    }

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

/* A serialized copy is the same kernel, to the last bit: type, observables, knobs,
 * range, and the Limber window at every redshift and multipole checked; and preparing
 * either of them again changes nothing. */
static void
test_nc_xcor_kernel_serialize (TestNcXcorKernel *test, gconstpointer pdata)
{
  const gint ells[] = {2, 10, 50, 100, 500};
  const guint n_ell = G_N_ELEMENTS (ells);
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_KEEP_SHARED);
  NcXcorKernel *dup = NC_XCOR_KERNEL (ncm_serialize_dup_obj (ser, G_OBJECT (test->xclk)));
  GArray *orig      = _test_limber_values (test, test->xclk, ells, n_ell);
  gdouble zmin, zmax, zmid, dmin, dmax, dmid;
  GArray *first, *again;

  g_assert_true (dup != test->xclk);
  g_assert_true (G_OBJECT_TYPE (dup) == G_OBJECT_TYPE (test->xclk));
  g_assert_cmpuint (nc_xcor_kernel_obs_len (dup), ==, nc_xcor_kernel_obs_len (test->xclk));
  g_assert_cmpuint (nc_xcor_kernel_obs_params_len (dup), ==, nc_xcor_kernel_obs_params_len (test->xclk));
  g_assert_cmpfloat (nc_xcor_kernel_get_reltol (dup), ==, nc_xcor_kernel_get_reltol (test->xclk));
  g_assert_cmpfloat (nc_xcor_kernel_get_peak_epsilon (dup), ==, nc_xcor_kernel_get_peak_epsilon (test->xclk));
  g_assert_cmpuint (nc_xcor_kernel_get_max_border_expansions (dup), ==, nc_xcor_kernel_get_max_border_expansions (test->xclk));
  g_assert_cmpfloat (nc_xcor_kernel_get_expansion_factor (dup), ==, nc_xcor_kernel_get_expansion_factor (test->xclk));

  nc_xcor_kernel_prepare (dup, test->cosmo);

  nc_xcor_kernel_get_z_range (test->xclk, &zmin, &zmax, &zmid);
  nc_xcor_kernel_get_z_range (dup, &dmin, &dmax, &dmid);
  g_assert_cmpfloat (dmin, ==, zmin);
  g_assert_cmpfloat (dmax, ==, zmax);

  first = _test_limber_values (test, dup, ells, n_ell);
  _test_assert_same_values (orig, first);

  nc_xcor_kernel_prepare (dup, test->cosmo);
  again = _test_limber_values (test, dup, ells, n_ell);
  _test_assert_same_values (first, again);
  g_array_unref (again);

  nc_xcor_kernel_prepare (test->xclk, test->cosmo);
  again = _test_limber_values (test, test->xclk, ells, n_ell);
  _test_assert_same_values (orig, again);

  g_array_unref (again);
  g_array_unref (first);
  g_array_unref (orig);
  nc_xcor_kernel_free (dup);
  ncm_serialize_free (ser);
}

static void
test_nc_xcor_kernel_limber_backends (TestNcXcorKernel *test, gconstpointer pdata)
{
  _test_limber_backends_pair (test, test->xclk, test->xclk);
}

static void
test_nc_xcor_kernel_limber_backends_cross (TestNcXcorKernel *test, gconstpointer pdata)
{
  const TestKernel *tk = pdata;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_crosses); i++)
  {
    if (g_strcmp0 (tk->name, test_crosses[i].name) == 0)
    {
      NcXcorKernel *partner = test_crosses[i].build_partner (test);

      nc_xcor_kernel_prepare (partner, test->cosmo);
      _test_limber_backends_pair (test, test->xclk, partner);
      nc_xcor_kernel_free (partner);

      return;
    }
  }

  g_assert_not_reached ();
}

/* The kernels that normalize a dN/dz do it on a copy: the caller's spline comes out of
 * construction unchanged, two kernels built on one spline are the same kernel, and a
 * spline set after construction gives the kernel built with it. */
static void
test_nc_xcor_kernel_dndz_untouched (void)
{
  const TestDndzBuild builds[2] = {&_dndz_build_weak_lensing, &_dndz_build_gal};
  NcHICosmo *cosmo              = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist              = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps                = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                                         NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  guint b;

  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));

  for (b = 0; b < 2; b++)
  {
    NcmSpline *dndz      = _dndz_at_new (TEST_DNDZ_MEAN);
    NcmSpline *dndz2     = _dndz_at_new (0.8);
    NcmVector *y_before  = ncm_vector_dup (ncm_spline_peek_yv (dndz));
    NcXcorKernel *first  = builds[b](dist, ps, dndz);
    NcXcorKernel *second = builds[b](dist, ps, dndz);
    NcXcorKernel *moved  = builds[b](dist, ps, dndz2);
    gdouble v1[12], v2[12];
    guint i;

    for (i = 0; i < ncm_vector_len (y_before); i++)
      g_assert_cmpfloat (ncm_vector_get (ncm_spline_peek_yv (dndz), i), ==, ncm_vector_get (y_before, i));

    nc_xcor_kernel_prepare (first, cosmo);
    nc_xcor_kernel_prepare (second, cosmo);
    _dndz_values (first, cosmo, dist, v1);
    _dndz_values (second, cosmo, dist, v2);

    for (i = 0; i < 12; i++)
      g_assert_cmpfloat (v2[i], ==, v1[i]);

    /* Set after construction: normalized as at construction, and the kernel knows it
     * has to prepare again. */
    g_object_set (first, "dndz", dndz2, NULL);
    nc_xcor_kernel_prepare_if_needed (first, cosmo);
    nc_xcor_kernel_prepare (moved, cosmo);
    _dndz_values (first, cosmo, dist, v1);
    _dndz_values (moved, cosmo, dist, v2);

    for (i = 0; i < 12; i++)
      g_assert_cmpfloat (v1[i], ==, v2[i]);

    nc_xcor_kernel_free (first);
    nc_xcor_kernel_free (second);
    nc_xcor_kernel_free (moved);
    ncm_vector_free (y_before);
    ncm_spline_free (dndz);
    ncm_spline_free (dndz2);
  }

  ncm_powspec_free (ps);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
}

/* At s = 0.4 the magnification weight 5s - 2 is zero (exactly, in doubles): a galaxy
 * kernel with magnification on is then the kernel without it, though the magnification
 * component is still built and evaluated. */
static void
test_nc_xcor_kernel_gal_magbias_null (void)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps   = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                            NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcmSpline *dndz      = _dndz_new ();
  NcXcorKernelGal *on  = nc_xcor_kernel_gal_new (dist, ps, 1, 1.234, dndz, TRUE);
  NcXcorKernelGal *off = nc_xcor_kernel_gal_new (dist, ps, 1, 1.234, dndz, FALSE);
  GPtrArray *comps;
  guint i, j;

  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));

  ncm_model_orig_vparam_set (NCM_MODEL (on), NC_XCOR_KERNEL_GAL_BIAS, 0, 1.5);
  ncm_model_orig_vparam_set (NCM_MODEL (off), NC_XCOR_KERNEL_GAL_BIAS, 0, 1.5);
  ncm_model_orig_param_set (NCM_MODEL (on), NC_XCOR_KERNEL_GAL_MAG_BIAS, 0.4);

  nc_xcor_kernel_prepare (NC_XCOR_KERNEL (on), cosmo);
  nc_xcor_kernel_prepare (NC_XCOR_KERNEL (off), cosmo);

  comps = nc_xcor_kernel_get_component_list (NC_XCOR_KERNEL (on));
  g_assert_cmpuint (comps->len, ==, 2);
  g_ptr_array_unref (comps);

  for (j = 0; j < 6; j++)
  {
    const gdouble z    = 0.1 + 0.3 * j;
    const gint ells[3] = {2, 10, 200};

    for (i = 0; i < 3; i++)
      g_assert_cmpfloat (nc_xcor_kernel_eval_limber_z_full (NC_XCOR_KERNEL (on), cosmo, z, dist, ells[i]), ==,
                         nc_xcor_kernel_eval_limber_z_full (NC_XCOR_KERNEL (off), cosmo, z, dist, ells[i]));
  }

  /* And the term is there at any other s, or the equality above would prove nothing. */
  {
    gboolean differs = FALSE;

    ncm_model_orig_param_set (NCM_MODEL (on), NC_XCOR_KERNEL_GAL_MAG_BIAS, 0.6);
    nc_xcor_kernel_prepare_if_needed (NC_XCOR_KERNEL (on), cosmo);

    for (j = 0; j < 6; j++)
    {
      const gdouble z = 0.1 + 0.3 * j;

      differs = differs || (nc_xcor_kernel_eval_limber_z_full (NC_XCOR_KERNEL (on), cosmo, z, dist, 10) !=
                            nc_xcor_kernel_eval_limber_z_full (NC_XCOR_KERNEL (off), cosmo, z, dist, 10));
    }

    g_assert_true (differs);
  }

  nc_xcor_kernel_free (NC_XCOR_KERNEL (on));
  nc_xcor_kernel_free (NC_XCOR_KERNEL (off));
  ncm_spline_free (dndz);
  ncm_powspec_free (ps);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
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

/* A closure in u read at exactly the ends of a multipole's k range gives its value from
 * inside the range, as just inside it does, however k / nu rounds: the cluster top-hat is
 * far from zero at its band edges, where it steps to zero */
static void
test_nc_xcor_kernel_cluster_limber_range_ends (void)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps   = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                            NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcXcorKernel *xclk;
  NcXcorKernelIntegrand *closure;
  const gint lmin = 10;
  const gint lmax = 17;
  gdouble W_end[TEST_LMAX], W_in[TEST_LMAX];
  gint l;

  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));

  xclk = NC_XCOR_KERNEL (nc_xcor_kernel_cluster_tophat_new (dist, ps, 0.3, 0.6));
  nc_xcor_kernel_set_l_limber (xclk, 0);
  nc_xcor_kernel_prepare (xclk, cosmo);

  closure = nc_xcor_kernel_get_eval_vectorized (xclk, cosmo, lmin, lmax);

  for (l = lmin; l <= lmax; l++)
  {
    const guint r = l - lmin;
    gdouble k_lo, k_hi;

    nc_xcor_kernel_integrand_get_range_comp (closure, r, &k_lo, &k_hi);

    nc_xcor_kernel_integrand_eval (closure, k_hi, W_end);
    nc_xcor_kernel_integrand_eval (closure, k_hi * (1.0 - 1.0e-12), W_in);
    g_assert_cmpfloat (fabs (W_in[r]), >, 0.0);
    g_assert_cmpfloat (fabs (W_end[r] - W_in[r]), <=, 1.0e-9 * fabs (W_in[r]));

    nc_xcor_kernel_integrand_eval (closure, k_lo, W_end);
    nc_xcor_kernel_integrand_eval (closure, k_lo * (1.0 + 1.0e-12), W_in);
    g_assert_cmpfloat (fabs (W_in[r]), >, 0.0);
    g_assert_cmpfloat (fabs (W_end[r] - W_in[r]), <=, 1.0e-9 * fabs (W_in[r]));
  }

  nc_xcor_kernel_integrand_unref (closure);
  nc_xcor_kernel_free (xclk);
  ncm_powspec_free (ps);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
}

/* Under Limber in u each multipole's support ends at k_max / nu of the power spectrum,
 * one end per multipole. With the power spectrum capped where the window is still large,
 * every multipole must reach the cap: between (nu / nu_max) k_max and k_max it matches
 * the Limber formula in z within what the closure was asked for. */
static void
test_nc_xcor_kernel_gal_limber_capped_powspec (void)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps   = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                            NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcmSpline *dndz    = _dndz_new ();
  const gint lmin    = 100;
  const gint lmax    = 131;
  const gdouble nu_0 = lmin + 0.5;
  const gdouble nu_n = lmax + 0.5;
  NcXcorKernelGal *g;
  NcXcorKernelIntegrand *closure;
  gdouble RH, k_cap, peak = 0.0, bound;
  gdouble W[TEST_LMAX];
  gint l;
  guint j;

  nc_distance_prepare (dist, cosmo);
  RH = nc_hicosmo_RH_Mpc (cosmo);

  /* The cap at chi = chi(z = 0.6) for the middle multipole, near the peak of dN/dz */
  k_cap = 0.5 * (nu_0 + nu_n) / nc_distance_comoving (dist, cosmo, 0.6);
  ncm_powspec_set_kmax (ps, k_cap / RH);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));

  g = g_object_new (NC_TYPE_XCOR_KERNEL_GAL,
                    "dist", dist,
                    "powspec", ps,
                    "bparam-length", (gsize) 1,
                    "nbarm1", 0.0,
                    "dndz", dndz,
                    "domagbias", FALSE,
                    "dorsd", FALSE,
                    NULL);
  ncm_model_orig_vparam_set (NCM_MODEL (g), NC_XCOR_KERNEL_GAL_BIAS, 0, 1.5);
  nc_xcor_kernel_set_l_limber (NC_XCOR_KERNEL (g), 0);
  nc_xcor_kernel_prepare (NC_XCOR_KERNEL (g), cosmo);

  closure = nc_xcor_kernel_get_eval_vectorized (NC_XCOR_KERNEL (g), cosmo, lmin, lmax);

  /* The block's peak, from the Limber formula between the far edge of dN/dz, z = 2, and the cap */
  for (l = lmin; l <= lmax; l++)
  {
    const gdouble nu    = l + 0.5;
    const gdouble k_far = nu / nc_distance_comoving (dist, cosmo, 2.0);

    for (j = 0; j <= 200; j++)
    {
      const gdouble k  = k_far + (0.999 * k_cap - k_far) * j / 200.0;
      const gdouble z  = nc_distance_inv_comoving (dist, cosmo, nu / k);
      const gdouble Wz = nc_xcor_kernel_eval_limber_z_full (NC_XCOR_KERNEL (g), cosmo, z, dist, l) *
                         sqrt (ncm_powspec_eval (ps, NCM_MODEL (cosmo), z, k / RH)) * sqrt (M_PI / (2.0 * nu)) / k;

      peak = MAX (peak, fabs (Wz));
    }
  }

  bound = MAX (nc_xcor_kernel_get_reltol (NC_XCOR_KERNEL (g)), nc_xcor_kernel_get_peak_epsilon (NC_XCOR_KERNEL (g))) * peak;

  for (l = lmin; l < lmax; l++)
  {
    const gdouble nu    = l + 0.5;
    const gdouble k_low = k_cap * nu / nu_n;

    for (j = 1; j <= 3; j++)
    {
      const gdouble k  = k_low + (k_cap - k_low) * j / 4.0;
      const gdouble z  = nc_distance_inv_comoving (dist, cosmo, nu / k);
      const gdouble Wz = nc_xcor_kernel_eval_limber_z_full (NC_XCOR_KERNEL (g), cosmo, z, dist, l) *
                         sqrt (ncm_powspec_eval (ps, NCM_MODEL (cosmo), z, k / RH)) * sqrt (M_PI / (2.0 * nu)) / k;

      nc_xcor_kernel_integrand_eval (closure, k, W);
      g_assert_cmpfloat (fabs (W[l - lmin] - Wz), <=, bound);
    }
  }

  nc_xcor_kernel_integrand_unref (closure);
  nc_xcor_kernel_free (NC_XCOR_KERNEL (g));
  ncm_spline_free (dndz);
  ncm_powspec_free (ps);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
}

/* The j_m peak term of a component under Limber, sqrt(pi / (2 (m + 1/2))) K(chi, k) / k
 * at chi = (m + 1/2) / k, K the window factor times the wave number factor, zero outside
 * the component's chi range */
static gdouble
_limber_peak (NcXcorKernelComponent *comp, NcHICosmo *cosmo, NcDistance *dist, gdouble m_half,
              gdouble k, gdouble chi_min, gdouble chi_max)
{
  const gdouble chi = m_half / k;
  NcXcorKinetic xck;

  if ((chi < chi_min) || (chi > chi_max))
    return 0.0;

  xck.chi_z = chi;
  xck.z     = nc_distance_inv_comoving (dist, cosmo, chi);
  xck.E_z   = nc_hicosmo_E (cosmo, xck.z);

  return sqrt (M_PI / (2.0 * m_half)) * nc_xcor_kernel_component_eval_window (comp, cosmo, &xck) *
         nc_xcor_kernel_component_eval_kfactor (comp, cosmo, &xck, k) / k;
}

/* The coefficient of j_{l+t} in the Bessel weight of order deriv, from the order
 * recurrence of DLMF 10.51.1: j_l' = (l j_{l-1} - (l + 1) j_{l+1}) / (2 l + 1), applied
 * twice for j_l''. Zero for a negative order. */
static gdouble
_limber_order_coef (guint deriv, gint l, gint t)
{
  const gdouble dl = l;

  if (l + t < 0)
    return 0.0;

  if (deriv == 0)
    return (t == 0) ? 1.0 : 0.0;

  if (deriv == 1)
    return (t == -1) ? dl / (2.0 * dl + 1.0) : ((t == 1) ? -(dl + 1.0) / (2.0 * dl + 1.0) : 0.0);

  switch (t)
  {
    case -2:
      return dl * (dl - 1.0) / ((2.0 * dl - 1.0) * (2.0 * dl + 1.0));

    case 0:
      return -(2.0 * dl * dl + 2.0 * dl - 1.0) / ((2.0 * dl - 1.0) * (2.0 * dl + 3.0));

    case 2:
      return (dl + 1.0) * (dl + 2.0) / ((2.0 * dl + 1.0) * (2.0 * dl + 3.0));

    default:
      return 0.0;
  }
}

/* W_l(k) under Limber from the components: each Bessel weight is its combination of the
 * orders l + t, and each order its peak term P_{l+t} */
static gdouble
_limber_direct (GPtrArray *comps, NcHICosmo *cosmo, NcDistance *dist, gint l, gdouble k)
{
  gdouble W = 0.0;
  guint i;

  for (i = 0; i < comps->len; i++)
  {
    NcXcorKernelComponent *comp = g_ptr_array_index (comps, i);
    const guint deriv           = nc_xcor_kernel_component_get_bessel_deriv (comp);
    gdouble chi_min, chi_max, k_min, k_max;
    gdouble val = 0.0;
    gint t;

    nc_xcor_kernel_component_get_limits (comp, cosmo, &chi_min, &chi_max, &k_min, &k_max);

    if ((k < k_min) || (k > k_max))
      continue;

    for (t = -2; t <= 2; t++)
    {
      const gdouble coef = _limber_order_coef (deriv, l, t);

      if (coef != 0.0)
        val += coef * _limber_peak (comp, cosmo, dist, l + t + 0.5, k, chi_min, chi_max);
    }

    W += nc_xcor_kernel_component_eval_prefactor (comp, cosmo, k, l) * val;
  }

  return W;
}

/* Under Limber an RSD component's j_l'' weight is the combination of the orders l - 2, l
 * and l + 2, each a closure term of its own sampled in x = k / (nu + t), starting at
 * x = 1 / chi_max for every multipole. With a truncated dN/dz the window is nonzero there, a
 * step that is a panel edge of every term; the closure must build, hold three terms, and
 * their sum must be the Limber W_l(k) of the components across the band, its ends and steps
 * included, within the closure tolerance of each multipole's peak. */
static void
test_nc_xcor_kernel_gal_rsd_limber_steps (void)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps   = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                            NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcmSpline *dndz    = _dndz_truncated_new ();
  NcXcorKernelGal *g = g_object_new (NC_TYPE_XCOR_KERNEL_GAL,
                                     "dist", dist,
                                     "powspec", ps,
                                     "bparam-length", (gsize) 1,
                                     "nbarm1", 0.0,
                                     "dndz", dndz,
                                     "domagbias", FALSE,
                                     "dorsd", TRUE,
                                     NULL);
  const gint lmin = 2;
  const gint lmax = 33;
  const gint ls[] = { 2, 17, 33 };
  const guint n_k = 400;
  gdouble W[33 - 2 + 1];
  NcXcorKernelIntegrand *xclki;
  GPtrArray *comps;
  gdouble chi_lo = G_MAXDOUBLE, chi_hi = 0.0, tol;
  guint i, j, n;

  ncm_model_orig_vparam_set (NCM_MODEL (g), NC_XCOR_KERNEL_GAL_BIAS, 0, 1.5);
  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));
  nc_xcor_kernel_set_l_limber (NC_XCOR_KERNEL (g), 0);
  nc_xcor_kernel_prepare (NC_XCOR_KERNEL (g), cosmo);

  xclki = nc_xcor_kernel_get_eval_vectorized (NC_XCOR_KERNEL (g), cosmo, lmin, lmax);
  comps = nc_xcor_kernel_get_component_list (NC_XCOR_KERNEL (g));
  tol   = GSL_MAX (nc_xcor_kernel_get_reltol (NC_XCOR_KERNEL (g)), nc_xcor_kernel_get_peak_epsilon (NC_XCOR_KERNEL (g)));

  g_assert_cmpuint (ncm_spectral_panels_get_n_panels (nc_xcor_kernel_integrand_peek_panels (xclki)), >, 0);
  g_assert_cmpuint (nc_xcor_kernel_integrand_get_n_terms (xclki), ==, 3);

  for (i = 0; i < comps->len; i++)
  {
    gdouble chi_min, chi_max, k_min, k_max;

    nc_xcor_kernel_component_get_limits (g_ptr_array_index (comps, i), cosmo, &chi_min, &chi_max, &k_min, &k_max);
    chi_lo = GSL_MIN (chi_lo, chi_min);
    chi_hi = GSL_MAX (chi_hi, chi_max);
  }

  for (n = 0; n < G_N_ELEMENTS (ls); n++)
  {
    const gint l       = ls[n];
    const gdouble nu   = l + 0.5;
    const gdouble k_lo = 0.9 * (nu - 2.0) / chi_hi;
    const gdouble k_hi = 1.1 * (nu + 2.0) / chi_lo;
    gdouble peak       = 0.0, dev = 0.0;

    for (j = 0; j < n_k; j++)
    {
      const gdouble k = k_lo * pow (k_hi / k_lo, j / (n_k - 1.0));

      peak = GSL_MAX (peak, fabs (_limber_direct (comps, cosmo, dist, l, k)));
    }

    for (j = 0; j < n_k; j++)
    {
      const gdouble k = k_lo * pow (k_hi / k_lo, j / (n_k - 1.0));

      nc_xcor_kernel_integrand_eval (xclki, k, W);
      dev = GSL_MAX (dev, fabs (W[l - lmin] - _limber_direct (comps, cosmo, dist, l, k)));
    }

    g_assert_cmpfloat (peak, >, 0.0);
    g_assert_cmpfloat (dev, <=, tol * peak);
  }

  g_ptr_array_unref (comps);
  nc_xcor_kernel_integrand_unref (xclki);
  nc_xcor_kernel_free (NC_XCOR_KERNEL (g));
  ncm_spline_free (dndz);
  ncm_powspec_free (ps);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
}

/* RSD kernels have no representation in the redshift-space Limber tier, and asking
 * for one must abort rather than silently drop the term. */
static void
test_nc_xcor_kernel_gal_rsd_limber_z_st (void)
{
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps   = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                            NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcmSpline *dndz    = _dndz_new ();
  NcXcorKernelGal *g = g_object_new (NC_TYPE_XCOR_KERNEL_GAL,
                                     "dist", dist,
                                     "powspec", ps,
                                     "bparam-length", (gsize) 1,
                                     "nbarm1", 1.234,
                                     "dndz", dndz,
                                     "domagbias", FALSE,
                                     "dorsd", TRUE,
                                     NULL);

  ncm_model_orig_vparam_set (NCM_MODEL (g), NC_XCOR_KERNEL_GAL_BIAS, 0, 1.5);
  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));
  nc_xcor_kernel_prepare (NC_XCOR_KERNEL (g), cosmo);

  nc_xcor_kernel_eval_limber_z_full (NC_XCOR_KERNEL (g), cosmo, 0.5, dist, 10);

  g_assert_not_reached ();
}

static void
test_nc_xcor_kernel_gal_rsd_limber_z_trap (void)
{
  g_test_trap_subprocess ("/nc/xcor/kernel/gal/rsd_limber_z/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*redshift-space distortions are not supported*");
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

/* The cluster top-hat's radial window is normalized in comoving volume alone: its
 * Limber-z window chi_t^2 integrated over dchi = dz / E(z), times the prefactor
 * 1 / (V_upper - V_lower), is one, so that C_0 / 4 pi is the bin's super-sample
 * covariance. An extra 1 / (z_upper - z_lower) would rescale every C_l. GL(32) is exact
 * here to rounding, the integrand being smooth inside the bin. */
static void
test_nc_xcor_kernel_cluster_tophat_volume_normalized (void)
{
  const gdouble z_lower = 0.1;
  const gdouble z_upper = 0.2;
  NcHICosmo *cosmo      = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist      = nc_distance_new (TEST_ZMAX);
  NcmPowspec *ps        = NCM_POWSPEC (ncm_powspec_analytic_new (NCM_POWSPEC_ANALYTIC_SHAPE_BBKS,
                                                                 NCM_POWSPEC_ANALYTIC_GROWTH_LCDM));
  NcXcorKernel *xclk                = NC_XCOR_KERNEL (nc_xcor_kernel_cluster_tophat_new (dist, ps, z_lower, z_upper));
  gsl_integration_glfixed_table *gl = gsl_integration_glfixed_table_alloc (32);
  NcXcorKinetic xck                 = {0.0, 0.0, 0.0};
  gdouble integral                  = 0.0;
  guint i;

  nc_distance_prepare (dist, cosmo);
  ncm_powspec_prepare (ps, NCM_MODEL (cosmo));
  nc_xcor_kernel_prepare (xclk, cosmo);

  for (i = 0; i < gl->n; i++)
  {
    gdouble z, w;

    gsl_integration_glfixed_point (z_lower, z_upper, i, &z, &w, gl);
    integral += w * nc_xcor_kernel_eval_limber_z (xclk, cosmo, z, &xck, 0) / nc_hicosmo_E (cosmo, z);
  }

  ncm_assert_cmpdouble_e (integral * nc_xcor_kernel_eval_limber_z_prefactor (xclk, cosmo, 0), ==, 1.0, 1.0e-10, 0.0);

  gsl_integration_glfixed_table_free (gl);
  nc_xcor_kernel_free (xclk);
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

  g_assert_true (gsl_finite (nc_xcor_kernel_radial_eval_kernel_factor (NC_XCOR_KERNEL_RADIAL (withk), 0,
                                                                       cosmo, 1600.0, 0.5)));

  nc_xcor_kernel_free (plain);
  nc_xcor_kernel_free (withk);
  nc_xcor_kernel_radial_kdep_free (kdep);
  ncm_powspec_free (ps);
  nc_distance_free (dist);
  nc_hicosmo_free (cosmo);
}

