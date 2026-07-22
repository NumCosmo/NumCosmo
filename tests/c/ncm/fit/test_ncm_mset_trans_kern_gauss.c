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

void test_ncm_mset_trans_kern_gauss_max_iter_property (void);
void test_ncm_mset_trans_kern_gauss_exhaustion (void);

int
main (int argc, char *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/mset_trans_kern_gauss/max_iter_property", test_ncm_mset_trans_kern_gauss_max_iter_property);
  g_test_add_func ("/ncm/mset_trans_kern_gauss/exhaustion", test_ncm_mset_trans_kern_gauss_exhaustion);

  g_test_run ();
}

void
test_ncm_mset_trans_kern_gauss_max_iter_property (void)
{
  NcmMSetTransKernGauss *tkerng = ncm_mset_trans_kern_gauss_new (1);
  guint max_iter                = 0;

  /* Default (see the "max-iter" g_param_spec_uint in class_init). */
  g_object_get (tkerng, "max-iter", &max_iter, NULL);
  g_assert_cmpuint (max_iter, ==, 1000);

  g_object_set (tkerng, "max-iter", 7u, NULL);
  g_object_get (tkerng, "max-iter", &max_iter, NULL);
  g_assert_cmpuint (max_iter, ==, 7);

  ncm_mset_trans_kern_free (NCM_MSET_TRANS_KERN (tkerng));
}

/*
 * _ncm_mset_trans_kern_gauss_generate() redraws until the proposal lands
 * within the free parameters' bounds or max-iter rounds are exhausted, in
 * which case it g_warning()s each still-out-of-bounds parameter and
 * g_error()s -- a fatal abort (see g_test_trap_subprocess() below, the
 * established idiom in this test suite for exercising a g_error() path
 * without crashing the outer test process). Forced here with max-iter=1,
 * bounds tightened to a width the huge proposal covariance essentially
 * never lands inside, and a starting theta safely centered in-bounds (only
 * the *proposal* needs to escape).
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

    NcmModelMVND *model_mvnd      = ncm_model_mvnd_new (1);
    NcmMSet *mset                 = ncm_mset_new (NCM_MODEL (model_mvnd), NULL, NULL);
    NcmMSetTransKernGauss *tkerng = ncm_mset_trans_kern_gauss_new (0);
    NcmRNG *rng                   = ncm_rng_seeded_new (NULL, 0);
    NcmMatrix *cov                = ncm_matrix_new (1, 1);
    NcmVector *theta              = ncm_vector_new (1);
    NcmVector *thetastar          = ncm_vector_new (1);

    ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
    ncm_model_param_set_lower_bound (NCM_MODEL (model_mvnd), 0, -0.01);
    ncm_model_param_set_upper_bound (NCM_MODEL (model_mvnd), 0,  0.01);
    ncm_mset_prepare_fparam_map (mset);
    ncm_vector_set (theta, 0, 0.0);
    ncm_matrix_set (cov, 0, 0, 1.0e10);

    g_object_set (tkerng, "max-iter", 1u, NULL);
    ncm_mset_trans_kern_set_mset (NCM_MSET_TRANS_KERN (tkerng), mset);
    ncm_mset_trans_kern_gauss_set_cov (tkerng, cov);

    ncm_mset_trans_kern_generate (NCM_MSET_TRANS_KERN (tkerng), theta, thetastar, rng);

    return; /* LCOV_EXCL_LINE */
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*is out of bounds*");
  g_test_trap_assert_stderr ("*failed to generate a valid sample after*");
}

