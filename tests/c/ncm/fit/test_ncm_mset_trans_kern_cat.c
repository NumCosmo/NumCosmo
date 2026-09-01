/***************************************************************************
 *            test_ncm_mset_trans_kern_cat.c
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
 * #NcmMSetTransKernCat's CHOOSE sampler, which draws a starting point by picking a row
 * of a catalog rather than from a fitted density. It needs no #NcmStatsDist, so a
 * catalog built by hand is the whole input -- nothing has to be sampled to test it.
 *
 * Two properties of the sampler decide what the fixture has to look like. It draws
 * without replacement, keyed on each row's m2lnL, so the rows need distinct m2lnL values
 * and no more draws than rows may be asked for. And it rejects any row outside the
 * mset's bounds, so the rows are laid inside them by construction rather than drawn and
 * hoped over.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_TKC_DIM 3
#define TEST_TKC_NROWS 64

typedef struct _TestNcmMSetTransKernCat
{
  NcmMSet *mset;
  NcmMSetCatalog *mcat;
  NcmMSetTransKernCat *tcat;
  NcmRNG *rng;
} TestNcmMSetTransKernCat;

static void
test_ncm_mset_trans_kern_cat_new (TestNcmMSetTransKernCat *test, gconstpointer pdata)
{
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (TEST_TKC_DIM);
  NcmMSet *mset            = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmRNG *rng              = ncm_rng_seeded_new (NULL, 24681357);
  NcmMSetCatalog *mcat;
  guint i, j;

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  mcat = ncm_mset_catalog_new (mset, 1, 1, FALSE, "m2lnL", "-2\\ln(L)", NULL);
  ncm_mset_catalog_set_m2lnp_var (mcat, 0);
  ncm_mset_catalog_set_run_type (mcat, "trans-kern-cat");

  {
    NcmVector *x = ncm_vector_new (ncm_mset_fparams_len (mset));

    for (i = 0; i < TEST_TKC_NROWS; i++)
    {
      /* Distinct, so the without-replacement bookkeeping has something to key on. */
      gdouble ax[1] = { 10.0 + i };

      for (j = 0; j < ncm_mset_fparams_len (mset); j++)
      {
        const gdouble lb = ncm_mset_fparam_get_lower_bound (mset, j);
        const gdouble ub = ncm_mset_fparam_get_upper_bound (mset, j);
        const gdouble u  = ncm_rng_uniform_gen (rng, 0.2, 0.8);

        ncm_vector_set (x, j, lb + (ub - lb) * u);
      }

      g_assert_true (ncm_mset_fparam_valid_bounds (mset, x));
      ncm_mset_catalog_add_from_vector_array (mcat, x, ax);
    }

    ncm_vector_clear (&x);
  }

  test->mset = mset;
  test->mcat = mcat;
  test->rng  = rng;
  test->tcat = ncm_mset_trans_kern_cat_new (mcat, NULL);

  ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (test->tcat), mset);

  /* The default is RBF interpolation, which would need a NcmStatsDist; CHOOSE is what
   * these check and it has to be asked for. */
  ncm_mset_trans_kern_cat_set_sampling (test->tcat, NCM_MSET_TRANS_KERN_CAT_SAMPLING_CHOOSE);

  ncm_model_mvnd_clear (&model_mvnd);
}

static void
test_ncm_mset_trans_kern_cat_free (TestNcmMSetTransKernCat *test, gconstpointer pdata)
{
  ncm_mset_catalog_clear (&test->mcat);
  ncm_mset_clear (&test->mset);
  ncm_rng_free (test->rng);

  NCM_TEST_FREE (ncm_mset_trans_kern_free, NCM_MSET_TRANS_KERN (test->tcat));
}

/* Is @thetastar one of the catalog's rows? The sampler picks, it does not interpolate,
 * so every point it returns has to be a row and not something between them. */
static gboolean
_test_tkc_is_catalog_row (TestNcmMSetTransKernCat *test, NcmVector *thetastar)
{
  const guint nadd = ncm_mset_catalog_nadd_vals (test->mcat);
  const guint len  = ncm_vector_len (thetastar);
  guint i, j;

  for (i = 0; i < ncm_mset_catalog_len (test->mcat); i++)
  {
    NcmVector *row = ncm_mset_catalog_peek_row (test->mcat, i);
    gboolean same  = TRUE;

    for (j = 0; j < len; j++)
    {
      if (ncm_vector_get (row, nadd + j) != ncm_vector_get (thetastar, j))
      {
        same = FALSE;
        break;
      }
    }

    if (same)
      return TRUE;
  }

  return FALSE;
}

static void
test_ncm_mset_trans_kern_cat_choose (TestNcmMSetTransKernCat *test, gconstpointer pdata)
{
  NcmVector *theta     = ncm_vector_new (ncm_mset_fparams_len (test->mset));
  NcmVector *thetastar = ncm_vector_new (ncm_mset_fparams_len (test->mset));
  const guint ndraws   = 16;
  guint i;

  g_assert_cmpuint (ncm_mset_trans_kern_cat_get_sampling (test->tcat), ==,
                    NCM_MSET_TRANS_KERN_CAT_SAMPLING_CHOOSE);

  ncm_mset_fparams_get_vector (test->mset, theta);

  for (i = 0; i < ndraws; i++)
  {
    ncm_vector_set_all (thetastar, GSL_NAN);
    ncm_mset_trans_kern_generate (NCM_MSET_TRANS_KERN (test->tcat), theta, thetastar, test->rng);

    g_assert_true (_test_tkc_is_catalog_row (test, thetastar));
    g_assert_true (ncm_mset_fparam_valid_bounds (test->mset, thetastar));
  }

  ncm_vector_free (theta);
  ncm_vector_free (thetastar);
}

/* The percentile cut restricts the draw to the best rows: with the cut on, every point
 * has to come from below the m2lnL threshold, which by construction is the low end of
 * the values written above. */
static void
test_ncm_mset_trans_kern_cat_choose_cut (TestNcmMSetTransKernCat *test, gconstpointer pdata)
{
  NcmVector *theta     = ncm_vector_new (ncm_mset_fparams_len (test->mset));
  NcmVector *thetastar = ncm_vector_new (ncm_mset_fparams_len (test->mset));
  const gdouble pct    = 0.5;
  guint nth            = 0;
  gdouble cut;
  guint i;

  g_object_set (test->tcat, "choose-cut", TRUE, "choose-percentile", pct, NULL);

  cut = ncm_mset_catalog_get_nth_m2lnL_percentile (test->mcat, pct, &nth);

  g_assert_true (gsl_finite (cut));
  g_assert_cmpuint (nth, <, TEST_TKC_NROWS);

  ncm_mset_fparams_get_vector (test->mset, theta);

  for (i = 0; i < nth; i++)
  {
    ncm_mset_trans_kern_generate (NCM_MSET_TRANS_KERN (test->tcat), theta, thetastar, test->rng);
    g_assert_true (_test_tkc_is_catalog_row (test, thetastar));
  }

  ncm_vector_free (theta);
  ncm_vector_free (thetastar);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/mset/trans_kern/cat/choose", TestNcmMSetTransKernCat, NULL,
              &test_ncm_mset_trans_kern_cat_new,
              &test_ncm_mset_trans_kern_cat_choose,
              &test_ncm_mset_trans_kern_cat_free);

  g_test_add ("/ncm/mset/trans_kern/cat/choose/cut", TestNcmMSetTransKernCat, NULL,
              &test_ncm_mset_trans_kern_cat_new,
              &test_ncm_mset_trans_kern_cat_choose_cut,
              &test_ncm_mset_trans_kern_cat_free);

  g_test_run ();
}

