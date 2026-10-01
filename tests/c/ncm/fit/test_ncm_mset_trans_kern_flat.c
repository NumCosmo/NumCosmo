/***************************************************************************
 *            test_ncm_mset_trans_kern_flat.c
 *
 *  Thu Oct 01 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

void test_ncm_mset_trans_kern_flat_generate (void);

int
main (int argc, char *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/mset_trans_kern_flat/generate", test_ncm_mset_trans_kern_flat_generate);

  g_test_run ();
}

void
test_ncm_mset_trans_kern_flat_generate (void)
{
  /* Two parameters with bounds [-1, 3] and [1/2, 1]: every draw is inside the box, the
   * draws differ, and the density is the inverse of the box area, 1/2, wherever it is
   * evaluated. */
  const gdouble lb[2]          = {-1.0, 0.5};
  const gdouble ub[2]          = {3.0, 1.0};
  NcmModelMVND *model_mvnd     = ncm_model_mvnd_new (2);
  NcmMSet *mset                = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  NcmMSetTransKernFlat *tkernf = ncm_mset_trans_kern_flat_new ();
  NcmMSetTransKern *tkern      = NCM_MSET_TRANS_KERN (tkernf);
  NcmRNG *rng                  = ncm_rng_seeded_new (NULL, 20261001);
  NcmVector *theta             = ncm_vector_new (2);
  NcmVector *thetastar         = ncm_vector_new (2);
  NcmVector *first             = ncm_vector_new (2);
  gboolean differ              = FALSE;
  guint i, k;

  for (k = 0; k < 2; k++)
  {
    ncm_model_param_set_lower_bound (NCM_MODEL (model_mvnd), k, lb[k]);
    ncm_model_param_set_upper_bound (NCM_MODEL (model_mvnd), k, ub[k]);
    ncm_model_orig_param_set (NCM_MODEL (model_mvnd), k, 0.5 * (lb[k] + ub[k]));
  }

  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);
  ncm_mset_trans_kern_set_mset (tkern, mset);
  ncm_mset_fparams_get_vector (mset, theta);

  g_assert_cmpstr (ncm_mset_trans_kern_get_name (tkern), ==, "Multivariate Flat Sampler");

  for (i = 0; i < 1000; i++)
  {
    ncm_mset_trans_kern_generate (tkern, theta, thetastar, rng);
    g_assert_true (ncm_mset_fparam_valid_bounds (mset, thetastar));

    if (i == 0)
      ncm_vector_memcpy (first, thetastar);
    else
      differ = differ || (ncm_vector_get (thetastar, 0) != ncm_vector_get (first, 0));

    g_assert_cmpfloat (ncm_mset_trans_kern_pdf (tkern, theta, thetastar), ==, 0.5);
  }

  g_assert_true (differ);

  ncm_vector_free (theta);
  ncm_vector_free (thetastar);
  ncm_vector_free (first);
  ncm_rng_free (rng);
  ncm_mset_trans_kern_free (tkern);
  ncm_mset_free (mset);
  ncm_model_mvnd_free (model_mvnd);
}

