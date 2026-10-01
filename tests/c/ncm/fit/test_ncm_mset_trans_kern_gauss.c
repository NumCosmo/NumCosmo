/***************************************************************************
 *            test_ncm_mset_trans_kern_gauss.c
 *
 *  Wed Jul 22 2026
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

void test_ncm_mset_trans_kern_gauss_generate_unbounded (void);
void test_ncm_mset_trans_kern_gauss_exhaustion (void);

int
main (int argc, char *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/mset_trans_kern_gauss/generate_unbounded", test_ncm_mset_trans_kern_gauss_generate_unbounded);
  g_test_add_func ("/ncm/mset_trans_kern_gauss/exhaustion", test_ncm_mset_trans_kern_gauss_exhaustion);

  g_test_run ();
}

/* A one-parameter model with bounds [-0.01, 0.01] and a proposal of variance 1e10. */
static void
_test_ncm_mset_trans_kern_gauss_tight (NcmMSet **mset, NcmMSetTransKernGauss **tkerng)
{
  NcmModelMVND *model_mvnd = ncm_model_mvnd_new (1);
  NcmMatrix *cov           = ncm_matrix_new (1, 1);

  *mset   = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
  *tkerng = ncm_mset_trans_kern_gauss_new (0);

  ncm_mset_param_set_all_ftype (*mset, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_lower_bound (NCM_MODEL (model_mvnd), 0, -0.01);
  ncm_model_param_set_upper_bound (NCM_MODEL (model_mvnd), 0,  0.01);
  ncm_model_orig_param_set (NCM_MODEL (model_mvnd), 0, 0.0);
  ncm_mset_prepare_fparam_map (*mset);
  ncm_matrix_set (cov, 0, 0, 1.0e10);

  ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (*tkerng), *mset);
  ncm_mset_trans_kern_gauss_set_cov (*tkerng, cov);

  ncm_matrix_free (cov);
  ncm_model_mvnd_free (model_mvnd);
}

void
test_ncm_mset_trans_kern_gauss_generate_unbounded (void)
{
  /* generate() draws once, so the proposal stays symmetric; a draw outside the bounds is
   * returned as it is, for the sampler to reject. */
  NcmMSet *mset                 = NULL;
  NcmMSetTransKernGauss *tkerng = NULL;
  NcmRNG *rng                   = ncm_rng_seeded_new (NULL, 0);
  NcmVector *theta              = ncm_vector_new (1);
  NcmVector *thetastar          = ncm_vector_new (1);

  _test_ncm_mset_trans_kern_gauss_tight (&mset, &tkerng);
  ncm_vector_set (theta, 0, 0.0);

  ncm_mset_trans_kern_generate (NCM_MSET_TRANS_KERN (tkerng), theta, thetastar, rng);
  g_assert_false (ncm_mset_fparam_valid_bounds (mset, thetastar));

  ncm_vector_free (theta);
  ncm_vector_free (thetastar);
  ncm_rng_free (rng);
  ncm_mset_trans_kern_free (NCM_MSET_TRANS_KERN (tkerng));
  ncm_mset_free (mset);
}

/*
 * ncm_mset_trans_kern_prior_sample() draws again outside the bounds; with tight bounds
 * and a huge proposal covariance it gives up after 1000 draws. The g_error() abort is
 * exercised via g_test_trap_subprocess().
 */
void
test_ncm_mset_trans_kern_gauss_exhaustion (void)
{
  if (g_test_subprocess ())
  {
    /* g_test_init() escalates WARNING/CRITICAL to fatal by default, so
     * without this reset the g_warning() below would itself abort the
     * subprocess before the g_error() ever runs, leaving it uncovered. */
    g_log_set_always_fatal (G_LOG_LEVEL_ERROR);

    NcmMSet *mset                 = NULL;
    NcmMSetTransKernGauss *tkerng = NULL;
    NcmRNG *rng                   = ncm_rng_seeded_new (NULL, 0);
    NcmVector *thetastar          = ncm_vector_new (1);

    _test_ncm_mset_trans_kern_gauss_tight (&mset, &tkerng);
    ncm_mset_trans_kern_set_prior_from_mset (NCM_MSET_TRANS_KERN (tkerng));
    ncm_mset_trans_kern_prior_sample (NCM_MSET_TRANS_KERN (tkerng), thetastar, rng);

    return; /* LCOV_EXCL_LINE */
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*is out of bounds*");
  g_test_trap_assert_stderr ("*failed to draw a sample within the bounds after 1000 draws*");
}

