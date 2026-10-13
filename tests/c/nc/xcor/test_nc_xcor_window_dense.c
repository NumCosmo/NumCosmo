/***************************************************************************
 *            test_nc_xcor_window_dense.c
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
 * The Chebyshev closures against the certified windows.
 *
 * Two tables hold I_ell^(d)(k) = int W(chi) j_ell^(d)(k chi) dchi of the analytic windows,
 * with the kernels that define them (tests/tools/make_xcor_window_bin.py): the dense one,
 * certified in Arb to a relative radius of 1e-25 on 50 k per decade over the closure's
 * domain, for d = 0, 1, 2 and ell 2..10, 50, 100, 200 (a data-file release file); and the
 * sparse one, 11 points in x/nu per multipole up to ell = 1000 (data/truth_tables/xcor),
 * of which only the multipoles the dense one lacks are checked. A closure of an analytic
 * kernel holds sqrt(P(k)) I_ell(k), so each certified point is compared with it directly.
 *
 * A multipole's scale is its peak. The dense table's points find it; the sparse table's
 * x/nu grid finds it for d = 0 but misses the peak of the derivative orders by 10 to 100
 * times, so there the scale is the closure's own largest |W| over a dense scan of its
 * domain. The deviations are still taken against the certified values.
 *
 * Per multipole, as a fraction of its peak: inside the closure's domain the deviation is
 * held to a multiple of #NcXcorKernel:peak-epsilon, measured at 1.9 at the defaults and
 * 1.2 at the tight setting over the 16 shapes and three orders (the excess at the
 * defaults is the windows' hard-edge oscillation in k, left unresolved at the size of
 * the tolerance); beyond it, the certified values the domain drops are held below
 * #NcXcorKernel:adaptive-epsilon times the block's norm, the truncation that sets the
 * domain.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define TEST_DENSE_FILE "xcor_window_ilk_dense.bin"
#define TEST_SPARSE_FILE "truth_tables/xcor/xcor_window_ilk.bin"
#define TEST_PEAK_SCAN 4096

/* One setting the closures are built at, with the deviation it is held to */
typedef struct _TestSetting
{
  const gchar *name;
  gdouble reltol;           /* <= 0: the kernel's default */
  gdouble peak_epsilon;     /* <= 0: the kernel's default */
  gdouble adaptive_epsilon; /* <= 0: the kernel's default */
  gdouble factor;           /* the deviation bound, in units of peak-epsilon */
} TestSetting;

static const TestSetting test_settings[] = {
  {"default", 0.0,    0.0,    0.0,    3.0},
  {"tight",   1.0e-8, 1.0e-7, 1.0e-9, 2.0},
};

/* One (shape, order) case: the kernel's key in its table, the table whose multipoles it
 * leaves out (NULL: none), and whether the scale is the closure's own peak */
typedef struct _TestCase
{
  gchar *shape;
  guint d;
  NcmObjDictStr *table;
  NcmObjDictStr *exclude;
  gboolean closure_peak;
} TestCase;

static NcmObjDictStr *test_dense  = NULL;
static NcmObjDictStr *test_sparse = NULL;

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

/* The multipoles @table certifies for @shape, ascending */
static GArray *
_test_ells (NcmObjDictStr *table, const gchar *shape)
{
  GStrv keys    = ncm_obj_dict_str_keys (table);
  gchar *prefix = g_strdup_printf ("table/%s|", shape);
  GArray *ells  = g_array_new (FALSE, FALSE, sizeof (gint));
  guint i;

  for (i = 0; keys[i] != NULL; i++)
  {
    const gchar *key = keys[i];

    if (g_str_has_prefix (key, prefix))
    {
      const gint ell = atoi (key + strlen (prefix));

      g_array_append_val (ells, ell);
    }
  }

  g_array_sort (ells, _test_cmp_int);

  g_free (prefix);
  g_free (keys); /* transfer container: the strings are the dictionary's */

  return ells;
}

/*
 * Builds the closure of @kernel over [lmin, lmax] at @ts and compares each multipole
 * with its certified values.
 */
static void
_test_block (NcXcorKernel *proto, NcHICosmo *cosmo, const TestCase *tc, const TestSetting *ts, gint lmin, gint lmax)
{
  NcmSerialize *ser         = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcXcorKernel *xclk        = NC_XCOR_KERNEL (ncm_serialize_dup_obj (ser, G_OBJECT (proto)));
  NcmSBesselIntegrator *sbi = NCM_SBESSEL_INTEGRATOR (ncm_sbessel_integrator_levin_new (lmin, lmax));
  const gdouble RH          = nc_hicosmo_RH_Mpc (cosmo);
  GArray *W                 = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *peak              = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *deviation         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *dropped           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble norm2             = 0.0;
  NcXcorKernelIntegrand *xclki;
  NcmPowspec *ps;
  gdouble k_lo, k_hi, peak_epsilon, adaptive_epsilon;
  gint ell;

  g_array_set_size (W, lmax - lmin + 1);
  g_array_set_size (peak, lmax - lmin + 1);
  g_array_set_size (deviation, lmax - lmin + 1);
  g_array_set_size (dropped, lmax - lmin + 1);

  nc_xcor_kernel_set_l_limber (xclk, -1);

  if (ts->reltol > 0.0)
    nc_xcor_kernel_set_reltol (xclk, ts->reltol);

  /* The tight setting is below the floor nc_xcor_kernel_set_peak_epsilon() warns at: on
   * the C_l integrand it is past what the outer integral carries, but on W_l(k) it is
   * what this test measures. */
  if (ts->peak_epsilon > 0.0)
  {
    if (ts->peak_epsilon < NC_XCOR_KERNEL_MIN_USEFUL_PEAK_EPSILON)
      g_test_expect_message (G_LOG_DOMAIN, G_LOG_LEVEL_WARNING, "*below the useful floor*");

    nc_xcor_kernel_set_peak_epsilon (xclk, ts->peak_epsilon);
    g_test_assert_expected_messages ();
  }

  if (ts->adaptive_epsilon > 0.0)
    nc_xcor_kernel_set_adaptive_epsilon (xclk, ts->adaptive_epsilon);

  peak_epsilon     = nc_xcor_kernel_get_peak_epsilon (xclk);
  adaptive_epsilon = nc_xcor_kernel_get_adaptive_epsilon (xclk);

  nc_xcor_kernel_prepare (xclk, cosmo);
  ps    = nc_xcor_kernel_peek_powspec (xclk);
  xclki = nc_xcor_kernel_get_eval_vectorized_full (xclk, cosmo, lmin, lmax, sbi);
  nc_xcor_kernel_integrand_get_range (xclki, &k_lo, &k_hi);

  /* Each multipole's deviation inside the domain, its certified values beyond it, and
   * its peak; the domain is truncated against the block's norm, so the dropped values
   * are held to that. */
  for (ell = lmin; ell <= lmax; ell++)
  {
    gchar *key       = g_strdup_printf ("table/%s|%d", tc->shape, ell);
    NcmMatrix *table = NCM_MATRIX (ncm_obj_dict_str_peek (tc->table, key));
    const guint n    = ncm_matrix_nrows (table);
    const guint i    = ell - lmin;
    guint j;

    g_array_index (peak, gdouble, i)      = 0.0;
    g_array_index (deviation, gdouble, i) = 0.0;
    g_array_index (dropped, gdouble, i)   = 0.0;

    for (j = 0; j < n; j++)
    {
      const gdouble k   = ncm_matrix_get (table, j, 0);
      const gdouble ref = ncm_matrix_get (table, j, 1 + tc->d) * sqrt (ncm_powspec_eval (ps, NCM_MODEL (cosmo), 0.0, k));

      g_array_index (peak, gdouble, i) = MAX (g_array_index (peak, gdouble, i), fabs (ref));

      if ((k * RH >= k_lo) && (k * RH <= k_hi))
      {
        nc_xcor_kernel_integrand_eval_comps (xclki, k * RH, i, 1, &g_array_index (W, gdouble, 0));
        g_array_index (deviation, gdouble, i) = MAX (g_array_index (deviation, gdouble, i), fabs (g_array_index (W, gdouble, i) - ref));
      }
      else
      {
        g_array_index (dropped, gdouble, i) = MAX (g_array_index (dropped, gdouble, i), fabs (ref));
      }
    }

    g_free (key);
  }

  if (tc->closure_peak)
  {
    for (ell = lmin; ell <= lmax; ell++)
    {
      const guint i = ell - lmin;
      guint j;

      for (j = 0; j < TEST_PEAK_SCAN; j++)
      {
        const gdouble k = k_lo * pow (k_hi / k_lo, (j + 0.5) / TEST_PEAK_SCAN);

        nc_xcor_kernel_integrand_eval_comps (xclki, k, i, 1, &g_array_index (W, gdouble, 0));
        g_array_index (peak, gdouble, i) = MAX (g_array_index (peak, gdouble, i), fabs (g_array_index (W, gdouble, i)));
      }
    }
  }

  for (ell = lmin; ell <= lmax; ell++)
    norm2 += gsl_pow_2 (g_array_index (peak, gdouble, ell - lmin));

  for (ell = lmin; ell <= lmax; ell++)
  {
    const guint i          = ell - lmin;
    const gdouble peak_i   = g_array_index (peak, gdouble, i);
    const gdouble dev_i    = g_array_index (deviation, gdouble, i);
    const gdouble drop_i   = g_array_index (dropped, gdouble, i);
    const gdouble dev_max  = ts->factor * peak_epsilon * peak_i;
    const gdouble drop_max = adaptive_epsilon * sqrt (norm2);

    if ((dev_i > dev_max) || (drop_i > drop_max))
      g_test_message ("%s d = %u %s ell = %d (block %d..%d): deviation %.3e of peak (bound %.1f peak-epsilon), "
                      "dropped %.3e of the block's norm (bound adaptive-epsilon)",
                      tc->shape, tc->d, ts->name, ell, lmin, lmax, dev_i / peak_i, ts->factor, drop_i / sqrt (norm2));

    g_assert_cmpfloat (dev_i, <=, dev_max);
    g_assert_cmpfloat (drop_i, <=, drop_max);
  }

  g_array_unref (peak);
  g_array_unref (deviation);
  g_array_unref (dropped);
  g_array_unref (W);
  nc_xcor_kernel_integrand_unref (xclki);
  ncm_sbessel_integrator_free (sbi);
  nc_xcor_kernel_free (xclk);
  ncm_serialize_free (ser);
}

/* Whether @tc checks @ell: present in its table for its order, absent from the excluded one */
static gboolean
_test_checks_ell (const TestCase *tc, gint ell)
{
  gchar *key       = g_strdup_printf ("table/%s|%d", tc->shape, ell);
  NcmMatrix *table = NCM_MATRIX (ncm_obj_dict_str_peek (tc->table, key));
  gboolean checks  = gsl_finite (ncm_matrix_get (table, 0, 1 + tc->d)) &&
                     ((tc->exclude == NULL) || (ncm_obj_dict_str_peek (tc->exclude, key) == NULL));

  g_free (key);

  return checks;
}

/*
 * One (shape, order): the closures at every setting, the multipoles 2..9 as one block,
 * as the solver builds them, and every other multipole alone.
 */
static void
test_nc_xcor_window (gconstpointer pdata)
{
  const TestCase *tc = pdata;
  gchar *kkey        = g_strdup_printf ("kernel/%s/d%u", tc->shape, tc->d);
  NcXcorKernel *xclk = NC_XCOR_KERNEL (ncm_obj_dict_str_peek (tc->table, kkey));
  NcHICosmo *cosmo   = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  GArray *ells       = _test_ells (tc->table, tc->shape);
  guint s, i;

  g_assert_nonnull (xclk);

  for (s = 0; s < G_N_ELEMENTS (test_settings); s++)
  {
    guint n_block = 0;

    for (i = 0; i < ells->len; i++)
      if ((g_array_index (ells, gint, i) >= 2) && (g_array_index (ells, gint, i) <= 9) &&
          _test_checks_ell (tc, g_array_index (ells, gint, i)))
        n_block++;

    if (n_block == 8)
      _test_block (xclk, cosmo, tc, &test_settings[s], 2, 9);

    for (i = 0; i < ells->len; i++)
    {
      const gint ell = g_array_index (ells, gint, i);

      if (((n_block != 8) || (ell > 9)) && _test_checks_ell (tc, ell))
        _test_block (xclk, cosmo, tc, &test_settings[s], ell, ell);
    }
  }

  g_array_unref (ells);
  nc_hicosmo_free (cosmo);
  g_free (kkey);
}

/* Registers a test per kernel of @table under @prefix, leaving out @exclude's multipoles */
static void
_test_register (NcmObjDictStr *table, NcmObjDictStr *exclude, gboolean closure_peak, const gchar *prefix)
{
  GStrv all_keys  = ncm_obj_dict_str_keys (table);
  GPtrArray *keys = g_ptr_array_new ();
  guint i;

  for (i = 0; all_keys[i] != NULL; i++)
    g_ptr_array_add (keys, all_keys[i]);

  g_ptr_array_sort (keys, _test_cmp_str);

  for (i = 0; i < keys->len; i++)
  {
    const gchar *key = g_ptr_array_index (keys, i);

    if (g_str_has_prefix (key, "kernel/"))
    {
      gchar **parts = g_strsplit (key, "/", 3);
      TestCase *tc  = g_new0 (TestCase, 1);
      gchar *path;

      tc->shape        = g_strdup (parts[1]);
      tc->d            = atoi (parts[2] + 1);
      tc->table        = table;
      tc->exclude      = exclude;
      tc->closure_peak = closure_peak;
      path             = g_strdup_printf ("%s/d%u/%s", prefix, tc->d, tc->shape);

      g_test_add_data_func_full (path, tc, &test_nc_xcor_window, &_test_case_free);

      g_free (path);
      g_strfreev (parts);
    }
  }

  g_ptr_array_unref (keys);
  g_free (all_keys); /* transfer container */
}

gint
main (gint argc, gchar *argv[])
{
  NcmSerialize *ser;
  gchar *dense_file;
  gchar *sparse_file;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  dense_file  = nc_data_download_release_file (TEST_DENSE_FILE);
  sparse_file = ncm_cfg_get_data_filename (TEST_SPARSE_FILE, TRUE);
  ser         = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  test_dense  = ncm_serialize_dict_str_from_binfile (ser, dense_file);
  ncm_serialize_reset (ser, TRUE);
  test_sparse = ncm_serialize_dict_str_from_binfile (ser, sparse_file);

  _test_register (test_dense, NULL, FALSE, "/nc/xcor/window_dense");
  _test_register (test_sparse, test_dense, TRUE, "/nc/xcor/window_sparse");

  g_free (dense_file);
  g_free (sparse_file);

  g_test_run ();

  ncm_obj_dict_str_unref (test_dense);
  ncm_obj_dict_str_unref (test_sparse);
  ncm_serialize_free (ser);
}

