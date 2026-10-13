/***************************************************************************
 *            test_nc_xcor_window_levin.c
 *
 *  Tue Oct 13 2026
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
 * The radial integral against certified Arb values, through the Levin integrator itself.
 *
 * The sparse table (data/truth_tables/xcor/xcor_window_ilk.bin, written by
 * tests/tools/make_xcor_window_bin.py from xcor_window_ilk*.json.gz) holds
 *
 *   I_ell^(d)(k) = int W(chi) j_ell^(d)(k chi) dchi,  chi in Mpc, k in 1/Mpc,
 *
 * for the analytic windows, 11 points in x/nu per multipole up to ell = 1000 (d = 2 to
 * 500), with the kernels that define them. Here each component's window is integrated
 * over its support by the integrator directly, not through a closure, so this is the
 * guard to run before and after touching the sbessel integrators. The window is a hard
 * zero outside its support and the k chi round trip puts nodes a hair beyond the edge, a
 * cliff the Chebyshev fit cannot resolve, so the window is clamped to the support, as the
 * library's own path does.
 *
 * Every point is held to rtol |I| + atol peak, the peak being the multipole's largest
 * certified |I|. Measured at d = 0 over the table: worst 6.9e-10 relative within three
 * orders of the peak and 5.7e-11 of the peak, so rtol 1e-8 and atol 1e-11. At d = 1, 2:
 * worst 1.7e-10 at d = 1, 2.1e-9 relative and 2.8e-9 of the peak at d = 2 (j_ell'' is
 * reached from neighbouring orders whose leading terms cancel), so rtol and atol 1e-8.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define TEST_SPARSE_FILE "truth_tables/xcor/xcor_window_ilk.bin"
#define TEST_RTOL 1.0e-8
#define TEST_ATOL_FRAC 1.0e-11
#define TEST_DERIV_RTOL 1.0e-8
#define TEST_DERIV_ATOL_FRAC 1.0e-8
#define TEST_TAU_MIN_OSC 200.0
#define TEST_BATCH 8

/* How the integral is taken */
typedef enum _TestPath
{
  TEST_PATH_ELL,         /* integrate_ell, one multipole */
  TEST_PATH_DERIV,       /* integrate_deriv over a block, its first multipole */
  TEST_PATH_TAU,         /* integrate_ell with the tau constraint */
  TEST_PATH_TAU_BATCHED, /* integrate over a block of TEST_BATCH, with the tau constraint */
} TestPath;

/* One case: the window, the path, the Bessel-derivative order, the block size and, for
 * the tau paths, whether the constraint is forced everywhere instead of by the rule */
typedef struct _TestCase
{
  gchar *shape;
  TestPath path;
  guint d;
  guint n_block;
  gboolean forced;
} TestCase;

/* The window of one component, clamped to its support */
typedef struct _TestWindow
{
  NcXcorKernelRadial *kernel;
  guint comp;
  gdouble chi_min;
  gdouble chi_max;
} TestWindow;

static NcmObjDictStr *test_table = NULL;

static gint
_test_cmp_int (gconstpointer a, gconstpointer b)
{
  return *(const gint *) a - *(const gint *) b;
}

static gint
_test_cmp_str (gconstpointer a, gconstpointer b)
{
  return strcmp (*(const gchar * const *) a, *(const gchar * const *) b);
}

static void
_test_case_free (gpointer data)
{
  TestCase *tc = data;

  g_free (tc->shape);
  g_free (tc);
}

static gdouble
_test_window (gpointer user_data, gdouble chi, gdouble k)
{
  TestWindow *w = user_data;

  return nc_xcor_kernel_radial_eval_W_comp (w->kernel, w->comp, GSL_MIN (GSL_MAX (chi, w->chi_min), w->chi_max));
}

/* The multipoles the table certifies for @shape, ascending */
static GArray *
_test_ells (const gchar *shape)
{
  GStrv keys    = ncm_obj_dict_str_keys (test_table);
  gchar *prefix = g_strdup_printf ("table/%s|", shape);
  GArray *ells  = g_array_new (FALSE, FALSE, sizeof (gint));
  guint i;

  for (i = 0; keys[i] != NULL; i++)
  {
    if (g_str_has_prefix (keys[i], prefix))
    {
      const gint ell = atoi (keys[i] + strlen (prefix));

      g_array_append_val (ells, ell);
    }
  }

  g_array_sort (ells, _test_cmp_int);

  g_free (prefix);
  g_free (keys); /* transfer container: the strings are the dictionary's */

  return ells;
}

/* The integrator of @tc for multipole @ell, with its constraint set */
static NcmSBesselIntegrator *
_test_integrator (const TestCase *tc, gint ell)
{
  const guint n_block              = ((tc->path == TEST_PATH_DERIV) || (tc->path == TEST_PATH_TAU_BATCHED)) ? tc->n_block : 1;
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new (ell, ell + n_block - 1);

  if ((tc->path == TEST_PATH_TAU) || (tc->path == TEST_PATH_TAU_BATCHED))
  {
    if (tc->forced)
    {
      /* The rule reapplies its own decision to every operator it configures, so it has
       * to be off for the solver-wide flag to reach the shallow panels. */
      ncm_sbessel_integrator_levin_set_tau_constraint_min_osc (sbilv, 0.0);
      ncm_sbessel_integrator_levin_set_tau_constraint (sbilv, TRUE);
    }
    else
    {
      ncm_sbessel_integrator_levin_set_tau_constraint_min_osc (sbilv, TEST_TAU_MIN_OSC);
    }
  }

  return NCM_SBESSEL_INTEGRATOR (sbilv);
}

/* I_ell^(d)(k) summed over the components, through @tc's path */
static gdouble
_test_integral (const TestCase *tc, NcmSBesselIntegrator *sbi, NcXcorKernelRadial *kernel, gint ell, gdouble k, NcmVector *block)
{
  const guint n_comps = nc_xcor_kernel_radial_get_n_comps (kernel);
  gdouble total       = 0.0;
  guint comp;

  for (comp = 0; comp < n_comps; comp++)
  {
    TestWindow w = {kernel, comp, 0.0, 0.0};

    nc_xcor_kernel_radial_get_comp_support (kernel, comp, &w.chi_min, &w.chi_max);

    switch (tc->path)
    {
      case TEST_PATH_ELL:
      case TEST_PATH_TAU:
        total += ncm_sbessel_integrator_integrate_ell (sbi, &_test_window, w.chi_min, w.chi_max, k, ell, &w);
        break;
      case TEST_PATH_DERIV:
        ncm_sbessel_integrator_integrate_deriv (sbi, &_test_window, w.chi_min, w.chi_max, k, tc->d, block, &w);
        total += ncm_vector_get (block, 0);
        break;
      case TEST_PATH_TAU_BATCHED:
        ncm_sbessel_integrator_integrate (sbi, &_test_window, w.chi_min, w.chi_max, k, block, &w);
        total += ncm_vector_get (block, 0);
        break;
      default:
        g_assert_not_reached ();
        break;
    }
  }

  return total;
}

/*
 * One case: every multipole the table certifies at the order, every point compared.
 */
static void
test_nc_xcor_window_levin (gconstpointer pdata)
{
  const TestCase *tc         = pdata;
  const gdouble rtol         = (tc->d == 0) ? TEST_RTOL : TEST_DERIV_RTOL;
  const gdouble atol_frac    = (tc->d == 0) ? TEST_ATOL_FRAC : TEST_DERIV_ATOL_FRAC;
  gchar *kkey                = g_strdup_printf ("kernel/%s/d0", tc->shape);
  NcmSerialize *ser          = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcXcorKernel *proto        = NC_XCOR_KERNEL (ncm_obj_dict_str_peek (test_table, kkey));
  NcXcorKernel *xclk         = NC_XCOR_KERNEL (ncm_serialize_dup_obj (ser, G_OBJECT (proto)));
  NcXcorKernelRadial *kernel = NC_XCOR_KERNEL_RADIAL (xclk);
  NcHICosmo *cosmo           = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  GArray *ells               = _test_ells (tc->shape);
  NcmVector *block           = ncm_vector_new (GSL_MAX (tc->n_block, 1));
  guint fallbacks            = 0;
  guint i;

  nc_xcor_kernel_set_l_limber (xclk, -1);
  nc_xcor_kernel_prepare (xclk, cosmo);

  for (i = 0; i < ells->len; i++)
  {
    const gint ell   = g_array_index (ells, gint, i);
    gchar *key       = g_strdup_printf ("table/%s|%d", tc->shape, ell);
    NcmMatrix *table = NCM_MATRIX (ncm_obj_dict_str_peek (test_table, key));
    const guint n    = ncm_matrix_nrows (table);
    gdouble peak     = 0.0;
    gdouble worst    = 0.0;
    NcmSBesselIntegrator *sbi;
    guint j;

    if (!gsl_finite (ncm_matrix_get (table, 0, 1 + tc->d)))
    {
      g_free (key);
      continue;
    }

    for (j = 0; j < n; j++)
      peak = GSL_MAX (peak, fabs (ncm_matrix_get (table, j, 1 + tc->d)));

    sbi = _test_integrator (tc, ell);

    for (j = 0; j < n; j++)
    {
      const gdouble k       = ncm_matrix_get (table, j, 0);
      const gdouble ref     = ncm_matrix_get (table, j, 1 + tc->d);
      const gdouble got     = _test_integral (tc, sbi, kernel, ell, k, block);
      const gdouble allowed = rtol * fabs (ref) + atol_frac * peak;

      worst = GSL_MAX (worst, fabs (got - ref) / allowed);
    }

    if (worst > 1.0)
      g_test_message ("%s d = %u ell = %d: deviation %.3f times rtol %.0e plus atol %.0e of peak",
                      tc->shape, tc->d, ell, worst, rtol, atol_frac);

    g_assert_cmpfloat (worst, <=, 1.0);

    if ((tc->path == TEST_PATH_TAU) || (tc->path == TEST_PATH_TAU_BATCHED))
      fallbacks += ncm_sbessel_integrator_levin_get_n_constraint_fallbacks (NCM_SBESSEL_INTEGRATOR_LEVIN (sbi));

    ncm_sbessel_integrator_free (sbi);
    g_free (key);
  }

  /* The table samples below the turning point for every shape; without the guard those
   * entries were wrong by 1e10 and more. */
  if (tc->forced)
    g_assert_cmpuint (fallbacks, >, 0);

  ncm_vector_free (block);
  g_array_unref (ells);
  nc_hicosmo_free (cosmo);
  nc_xcor_kernel_free (xclk);
  ncm_serialize_free (ser);
  g_free (kkey);
}

/* Registers @tc under @path, taking ownership */
static void
_test_add (const gchar *path, const gchar *shape, TestPath tpath, guint d, guint n_block, gboolean forced)
{
  TestCase *tc = g_new0 (TestCase, 1);

  tc->shape   = g_strdup (shape);
  tc->path    = tpath;
  tc->d       = d;
  tc->n_block = n_block;
  tc->forced  = forced;

  g_test_add_data_func_full (path, tc, &test_nc_xcor_window_levin, &_test_case_free);
}

gint
main (gint argc, gchar *argv[])
{
  NcmSerialize *ser;
  GStrv all_keys;
  GPtrArray *shapes;
  gchar *filename;
  guint i;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  filename   = ncm_cfg_get_data_filename (TEST_SPARSE_FILE, TRUE);
  ser        = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  test_table = ncm_serialize_dict_str_from_binfile (ser, filename);
  all_keys   = ncm_obj_dict_str_keys (test_table);
  shapes     = g_ptr_array_new ();

  for (i = 0; all_keys[i] != NULL; i++)
    if (g_str_has_prefix (all_keys[i], "kernel/") && g_str_has_suffix (all_keys[i], "/d0"))
      g_ptr_array_add (shapes, all_keys[i]);

  g_ptr_array_sort (shapes, _test_cmp_str);

  for (i = 0; i < shapes->len; i++)
  {
    gchar **parts     = g_strsplit (g_ptr_array_index (shapes, i), "/", 3);
    const gchar *name = parts[1];
    gchar *path;
    guint d, n;

    path = g_strdup_printf ("/nc/xcor/window_levin/d0/%s", name);
    _test_add (path, name, TEST_PATH_ELL, 0, 1, FALSE);
    g_free (path);

    for (d = 1; d <= 2; d++)
    {
      for (n = 1; n <= TEST_BATCH; n += TEST_BATCH - 1)
      {
        path = g_strdup_printf ("/nc/xcor/window_levin/deriv/d%u/block%u/%s", d, n, name);
        _test_add (path, name, TEST_PATH_DERIV, d, n, FALSE);
        g_free (path);
      }
    }

    path = g_strdup_printf ("/nc/xcor/window_levin/tau/rule/%s", name);
    _test_add (path, name, TEST_PATH_TAU, 0, 1, FALSE);
    g_free (path);

    path = g_strdup_printf ("/nc/xcor/window_levin/tau/forced/%s", name);
    _test_add (path, name, TEST_PATH_TAU, 0, 1, TRUE);
    g_free (path);

    path = g_strdup_printf ("/nc/xcor/window_levin/tau_batched/rule/%s", name);
    _test_add (path, name, TEST_PATH_TAU_BATCHED, 0, TEST_BATCH, FALSE);
    g_free (path);

    path = g_strdup_printf ("/nc/xcor/window_levin/tau_batched/forced/%s", name);
    _test_add (path, name, TEST_PATH_TAU_BATCHED, 0, TEST_BATCH, TRUE);
    g_free (path);

    g_strfreev (parts);
  }

  g_ptr_array_unref (shapes);
  g_free (all_keys); /* transfer container */
  g_free (filename);

  g_test_run ();

  ncm_obj_dict_str_unref (test_table);
  ncm_serialize_free (ser);
}

