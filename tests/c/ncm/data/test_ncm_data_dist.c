/***************************************************************************
 *            test_ncm_data_dist.c
 *
 *  Tue Sep 29 03:00:00 2026
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
#include <glib.h>
#include <glib-object.h>

/*
 * Toy distributions with the rate lambda = mu_0 of a NcmModelMVND: an exponential,
 * -2 ln p (x) = -2 ln lambda + 2 lambda x, and in two dimensions the product of two.
 * inv_pdf returns the uniform numbers themselves, so a resample lies in [0, 1].
 */

static gdouble
_test_rate (NcmMSet *mset)
{
  return ncm_model_orig_vparam_get (ncm_mset_peek (mset, ncm_model_mvnd_id ()), NCM_MODEL_MVND_MEAN, 0);
}

#define TEST_TYPE_DIST1D (test_dist1d_get_type ())
G_DECLARE_FINAL_TYPE (TestDist1d, test_dist1d, TEST, DIST1D, NcmDataDist1d)

struct _TestDist1d
{
  NcmDataDist1d parent_instance;
};

G_DEFINE_TYPE (TestDist1d, test_dist1d, NCM_TYPE_DATA_DIST1D)

static void
test_dist1d_init (TestDist1d *td)
{
}

static gdouble
_test_dist1d_m2lnL_val (NcmDataDist1d *dist1d, NcmMSet *mset, gdouble x)
{
  const gdouble lambda = _test_rate (mset);

  return -2.0 * log (lambda) + 2.0 * lambda * x;
}

static gdouble
_test_dist1d_inv_pdf (NcmDataDist1d *dist1d, NcmMSet *mset, gdouble u)
{
  return u;
}

static void
test_dist1d_class_init (TestDist1dClass *klass)
{
  NCM_DATA_DIST1D_CLASS (klass)->dist1d_m2lnL_val = &_test_dist1d_m2lnL_val;
  NCM_DATA_DIST1D_CLASS (klass)->inv_pdf          = &_test_dist1d_inv_pdf;
}

#define TEST_TYPE_DIST2D (test_dist2d_get_type ())
G_DECLARE_FINAL_TYPE (TestDist2d, test_dist2d, TEST, DIST2D, NcmDataDist2d)

struct _TestDist2d
{
  NcmDataDist2d parent_instance;
};

G_DEFINE_TYPE (TestDist2d, test_dist2d, NCM_TYPE_DATA_DIST2D)

static void
test_dist2d_init (TestDist2d *td)
{
}

static gdouble
_test_dist2d_m2lnL_val (NcmDataDist2d *dist2d, NcmMSet *mset, gdouble x, gdouble y)
{
  const gdouble lambda = _test_rate (mset);

  return -4.0 * log (lambda) + 2.0 * lambda * (x + y);
}

static void
_test_dist2d_inv_pdf (NcmDataDist2d *dist2d, NcmMSet *mset, gdouble u, gdouble v, gdouble *x, gdouble *y)
{
  *x = u;
  *y = v;
}

static void
test_dist2d_class_init (TestDist2dClass *klass)
{
  NCM_DATA_DIST2D_CLASS (klass)->dist2d_m2lnL_val = &_test_dist2d_m2lnL_val;
  NCM_DATA_DIST2D_CLASS (klass)->inv_pdf          = &_test_dist2d_inv_pdf;
}

/* Subclasses that implement nothing. */

#define TEST_TYPE_DIST1D_EMPTY (test_dist1d_empty_get_type ())
G_DECLARE_FINAL_TYPE (TestDist1dEmpty, test_dist1d_empty, TEST, DIST1D_EMPTY, NcmDataDist1d)

struct _TestDist1dEmpty
{
  NcmDataDist1d parent_instance;
};

G_DEFINE_TYPE (TestDist1dEmpty, test_dist1d_empty, NCM_TYPE_DATA_DIST1D)

static void
test_dist1d_empty_init (TestDist1dEmpty *td)
{
}

static void
test_dist1d_empty_class_init (TestDist1dEmptyClass *klass)
{
}

#define TEST_TYPE_DIST2D_EMPTY (test_dist2d_empty_get_type ())
G_DECLARE_FINAL_TYPE (TestDist2dEmpty, test_dist2d_empty, TEST, DIST2D_EMPTY, NcmDataDist2d)

struct _TestDist2dEmpty
{
  NcmDataDist2d parent_instance;
};

G_DEFINE_TYPE (TestDist2dEmpty, test_dist2d_empty, NCM_TYPE_DATA_DIST2D)

static void
test_dist2d_empty_init (TestDist2dEmpty *td)
{
}

static void
test_dist2d_empty_class_init (TestDist2dEmptyClass *klass)
{
}

void test_ncm_data_dist1d (void);
void test_ncm_data_dist2d (void);
void test_ncm_data_dist_errors (void);
void test_ncm_data_dist_errors_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/data_dist/1d", &test_ncm_data_dist1d);
  g_test_add_func ("/ncm/data_dist/2d", &test_ncm_data_dist2d);
  g_test_add_func ("/ncm/data_dist/errors", &test_ncm_data_dist_errors);
  g_test_add_func ("/ncm/data_dist/errors/subprocess", &test_ncm_data_dist_errors_subprocess);

  g_test_run ();
}

static NcmMSet *
_test_mset_new (void)
{
  NcmModelMVND *mvnd = ncm_model_mvnd_new (1);
  NcmMSet *mset      = ncm_mset_new (mvnd, NULL, NULL);

  ncm_model_orig_vparam_set (NCM_MODEL (mvnd), NCM_MODEL_MVND_MEAN, 0, 1.0);
  ncm_model_mvnd_free (mvnd);

  return mset;
}

void
test_ncm_data_dist1d (void)
{
  NcmMSet *mset = _test_mset_new ();
  NcmRNG *rng   = ncm_rng_seeded_new (NULL, 1);
  NcmData *data = g_object_new (TEST_TYPE_DIST1D, "n-points", 2, NULL);
  NcmVector *x  = ncm_data_dist1d_get_data (NCM_DATA_DIST1D (data));
  gdouble m2lnL;
  guint i;

  /* With lambda = 1, -2 ln L = 2 (x_1 + x_2). */
  ncm_vector_set (x, 0, 1.0);
  ncm_vector_set (x, 1, 2.0);
  ncm_data_set_init (data, TRUE);
  ncm_data_m2lnL_val (data, mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 6.0);

  /* A bootstrap of one draw sums one term. */
  ncm_data_bootstrap_create (data);
  ncm_bootstrap_set_bsize (ncm_data_peek_bootstrap (data), 1);
  ncm_data_bootstrap_resample (data, rng);
  ncm_data_m2lnL_val (data, mset, &m2lnL);
  g_assert_true ((m2lnL == 2.0) || (m2lnL == 4.0));
  ncm_data_bootstrap_remove (data);

  /* Resampling goes through inv_pdf. */
  ncm_data_resample (data, mset, rng);

  for (i = 0; i < 2; i++)
  {
    g_assert_cmpfloat (ncm_vector_get (x, i), >=, 0.0);
    g_assert_cmpfloat (ncm_vector_get (x, i), <=, 1.0);
  }

  ncm_vector_free (x);
  ncm_data_free (data);
  ncm_rng_free (rng);
  ncm_mset_free (mset);
}

void
test_ncm_data_dist2d (void)
{
  NcmMSet *mset = _test_mset_new ();
  NcmRNG *rng   = ncm_rng_seeded_new (NULL, 1);
  NcmData *data = g_object_new (TEST_TYPE_DIST2D, "n-points", 2, NULL);
  NcmMatrix *m  = ncm_data_dist2d_get_data (NCM_DATA_DIST2D (data));
  gdouble m2lnL, x, y;
  guint i;

  ncm_matrix_set (m, 0, 0, 1.0);
  ncm_matrix_set (m, 0, 1, 2.0);
  ncm_matrix_set (m, 1, 0, 0.5);
  ncm_matrix_set (m, 1, 1, 0.5);
  ncm_data_set_init (data, TRUE);
  ncm_data_m2lnL_val (data, mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 8.0);

  ncm_data_dist2d_inv_pdf (NCM_DATA_DIST2D (data), mset, 0.25, 0.75, &x, &y);
  g_assert_cmpfloat (x, ==, 0.25);
  g_assert_cmpfloat (y, ==, 0.75);

  ncm_data_resample (data, mset, rng);

  for (i = 0; i < 2; i++)
  {
    g_assert_cmpfloat (ncm_matrix_get (m, i, 0), <=, 1.0);
    g_assert_cmpfloat (ncm_matrix_get (m, i, 1), <=, 1.0);
  }

  ncm_matrix_free (m);
  ncm_data_free (data);
  ncm_rng_free (rng);
  ncm_mset_free (mset);
}

void
test_ncm_data_dist_errors (void)
{
  const gchar *cases[][2] = {
    {"1d_m2lnL",   "*`TestDist1dEmpty' does not implement dist1d_m2lnL_val*"},
    {"1d_inv_pdf", "*`TestDist1dEmpty' does not implement inv_pdf, so it cannot be resampled*"},
    {"2d_m2lnL",   "*`TestDist2dEmpty' does not implement dist2d_m2lnL_val*"},
    {"2d_inv_pdf", "*`TestDist2dEmpty' does not implement inv_pdf, so it cannot be resampled*"},
    {"no_points",  "*ncm_data_dist1d_get_data: data `TestDist1d' has no points*"},
    {"no_realization", "*data `TestDist1d': the bootstrap has no realization*"},
  };
  guint i;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    g_setenv ("TEST_NCM_DATA_DIST_ERROR", cases[i][0], TRUE);
    g_test_trap_subprocess ("/ncm/data_dist/errors/subprocess", 0, 0);
    g_test_trap_assert_failed ();
    g_test_trap_assert_stderr (cases[i][1]);
  }

  g_unsetenv ("TEST_NCM_DATA_DIST_ERROR");
}

void
test_ncm_data_dist_errors_subprocess (void)
{
  const gchar *which = g_getenv ("TEST_NCM_DATA_DIST_ERROR");
  NcmMSet *mset      = _test_mset_new ();
  NcmRNG *rng        = ncm_rng_seeded_new (NULL, 1);
  gdouble m2lnL;

  if (g_str_has_prefix (which, "1d"))
  {
    NcmData *data = g_object_new (TEST_TYPE_DIST1D_EMPTY, "n-points", 2, NULL);

    ncm_data_set_init (data, TRUE);

    if (g_str_has_suffix (which, "m2lnL"))
      ncm_data_m2lnL_val (data, mset, &m2lnL);
    else
      ncm_data_resample (data, mset, rng);
  }
  else if (g_str_has_prefix (which, "2d"))
  {
    NcmData *data = g_object_new (TEST_TYPE_DIST2D_EMPTY, "n-points", 2, NULL);

    ncm_data_set_init (data, TRUE);

    if (g_str_has_suffix (which, "m2lnL"))
      ncm_data_m2lnL_val (data, mset, &m2lnL);
    else
      ncm_data_resample (data, mset, rng);
  }
  else if (g_str_equal (which, "no_realization"))
  {
    NcmData *data = g_object_new (TEST_TYPE_DIST1D, "n-points", 2, NULL);

    ncm_data_set_init (data, TRUE);
    ncm_data_bootstrap_create (data);
    ncm_data_m2lnL_val (data, mset, &m2lnL);
  }
  else
  {
    NcmDataDist1d *dist1d = g_object_new (TEST_TYPE_DIST1D, NULL);

    ncm_data_dist1d_get_data (dist1d);
  }
}

