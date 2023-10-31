/***************************************************************************
 *            test_ncm_fit.c
 *
 *  Sun February 04 16:02:57 2018
 *  Copyright  2018  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2018 <vitenti@uel.br>
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

typedef struct _TestNcmFit
{
  NcmFit *fit;
  NcmRNG *rng;
  NcmDataGaussCovMVND *data_mvnd;
  guint ntests;
} TestNcmFit;

#define TESTS_NCM_DECL(lib, algo) \
        void test_ncm_fit_ ## lib ## _ ## algo ## _new (TestNcmFit * test, gconstpointer pdata); \
        void test_ncm_fit_ ## lib ## _ ## algo ## _new_empty (TestNcmFit * test, gconstpointer pdata); \
        void test_ncm_fit_ ## lib ## _ ## algo ## _traps (TestNcmFit * test, gconstpointer pdata);

#define TESTS_NCM_ADD(lib, algo) \
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_run, \
                    &test_ncm_fit_free); \
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run/empty", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new_empty, \
                    &test_ncm_fit_run_empty, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/traps", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _traps, \
                    &test_ncm_fit_free);

#define TESTS_NCM_ADD_INVALID(lib, algo) \
        g_test_add ("/ncm/fit/" #lib "/" #algo "/invalid/run/subprocess", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_invalid_run, \
                    &test_ncm_fit_free);

#define TESTS_NCM_NEW(lib, algo, lib_enum, algo_str, max_dim, max_iter) \
        void \
        test_ncm_fit_ ## lib ## _ ## algo ## _new (TestNcmFit * test, gconstpointer pdata) \
        { \
          const gint dim                 = g_test_rand_int_range (1, max_dim); \
          NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ()); \
          NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (dim, 1.0e-2, 1.0e0, 50.0, -1.0, 1.0, rng); \
          NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (dim); \
          NcmDataset *dset               = ncm_dataset_new_list (data_mvnd, NULL); \
          NcmLikelihood *lh              = ncm_likelihood_new (dset); \
          NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL); \
          NcmFit *fit; \
 \
          ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE); \
   \
          fit = ncm_fit_new (lib_enum, algo_str, lh, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL); \
          ncm_fit_set_maxiter (fit, max_iter); \
 \
          test->data_mvnd = ncm_data_gauss_cov_mvnd_ref (data_mvnd); \
          test->fit       = ncm_fit_ref (fit); \
          test->rng       = rng; \
   \
          g_assert_true (NCM_IS_FIT (fit)); \
 \
          ncm_data_gauss_cov_mvnd_clear (&data_mvnd); \
          ncm_model_mvnd_clear (&model_mvnd); \
          ncm_dataset_clear (&dset); \
          ncm_likelihood_clear (&lh); \
          ncm_mset_clear (&mset); \
          ncm_fit_clear (&fit); \
        } \
        void \
        test_ncm_fit_ ## lib ## _ ## algo ## _new_empty (TestNcmFit * test, gconstpointer pdata) \
        { \
          const gint dim                 = g_test_rand_int_range (1, max_dim); \
          NcmRNG *rng                    = ncm_rng_seeded_new (NULL, g_test_rand_int ()); \
          NcmDataGaussCovMVND *data_mvnd = ncm_data_gauss_cov_mvnd_new_full (dim, 1.0e-2, 1.0e0, 50.0, -1.0, 1.0, rng); \
          NcmModelMVND *model_mvnd       = ncm_model_mvnd_new (dim); \
          NcmDataset *dset               = ncm_dataset_new_list (data_mvnd, NULL); \
          NcmLikelihood *lh              = ncm_likelihood_new (dset); \
          NcmMSet *mset                  = ncm_mset_new (NCM_MODEL (model_mvnd), NULL); \
          NcmFit *fit; \
 \
          ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FIXED); \
   \
          fit = ncm_fit_new (lib_enum, algo_str, lh, mset, NCM_FIT_GRAD_NUMDIFF_CENTRAL); \
          ncm_fit_set_maxiter (fit, max_iter); \
 \
          test->data_mvnd = ncm_data_gauss_cov_mvnd_ref (data_mvnd); \
          test->fit       = ncm_fit_ref (fit); \
          test->rng       = rng; \
   \
          g_assert_true (NCM_IS_FIT (fit)); \
 \
          ncm_data_gauss_cov_mvnd_clear (&data_mvnd); \
          ncm_model_mvnd_clear (&model_mvnd); \
          ncm_dataset_clear (&dset); \
          ncm_likelihood_clear (&lh); \
          ncm_mset_clear (&mset); \
          ncm_fit_clear (&fit); \
        }

#if GLIB_CHECK_VERSION (2, 38, 0)
#define TESTS_NCM_TRAPS(lib, algo) \
        void \
        test_ncm_fit_ ## lib ## _ ## algo ## _traps (TestNcmFit * test, gconstpointer pdata) \
        { \
          g_test_trap_subprocess ("/ncm/fit/" #lib "/" #algo "/invalid/run/subprocess", 0, 0); \
          g_test_trap_assert_failed (); \
        }
#else
#define TESTS_NCM_TRAPS(lib, algo) \
        void \
        test_ncm_fit_ ## lib ## _ ## algo ## _traps (TestNcmFit * test, gconstpointer pdata) \
        { \
        }
#endif

#ifdef NUMCOSMO_HAVE_NLOPT
TESTS_NCM_DECL (nlopt, neldermead)
TESTS_NCM_DECL (nlopt, slsqp)
#endif /* NUMCOSMO_HAVE_NLOPT */

TESTS_NCM_DECL (gsl, ls)

TESTS_NCM_DECL (gsl, mm_conjugate_fr)
TESTS_NCM_DECL (gsl, mm_conjugate_pr)
TESTS_NCM_DECL (gsl, mm_vector_bfgs)
TESTS_NCM_DECL (gsl, mm_vector_bfgs2)
TESTS_NCM_DECL (gsl, mm_steepest_descent)

TESTS_NCM_DECL (gsl, nmsimplex)
TESTS_NCM_DECL (gsl, nmsimplex2)
TESTS_NCM_DECL (gsl, nmsimplex2rand)

void test_ncm_fit_free (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run_empty (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_invalid_run (TestNcmFit *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);

  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

#ifdef NUMCOSMO_HAVE_NLOPT
  TESTS_NCM_ADD (nlopt, neldermead)
  TESTS_NCM_ADD (nlopt, slsqp)
#endif /* NUMCOSMO_HAVE_NLOPT */

  TESTS_NCM_ADD (gsl, ls)

  TESTS_NCM_ADD (gsl, mm_conjugate_fr)
  TESTS_NCM_ADD (gsl, mm_conjugate_pr)
  TESTS_NCM_ADD (gsl, mm_vector_bfgs)
  TESTS_NCM_ADD (gsl, mm_vector_bfgs2)
  TESTS_NCM_ADD (gsl, mm_steepest_descent)

  TESTS_NCM_ADD (gsl, nmsimplex)
  TESTS_NCM_ADD (gsl, nmsimplex2)
  TESTS_NCM_ADD (gsl, nmsimplex2rand)

#if GLIB_CHECK_VERSION (2, 38, 0)
#ifdef NUMCOSMO_HAVE_NLOPT
  TESTS_NCM_ADD_INVALID (nlopt, neldermead)
  TESTS_NCM_ADD_INVALID (nlopt, slsqp)
#endif /* NUMCOSMO_HAVE_NLOPT */

  TESTS_NCM_ADD_INVALID (gsl, ls)

  TESTS_NCM_ADD_INVALID (gsl, mm_conjugate_fr)
  TESTS_NCM_ADD_INVALID (gsl, mm_conjugate_pr)
  TESTS_NCM_ADD_INVALID (gsl, mm_vector_bfgs)
  TESTS_NCM_ADD_INVALID (gsl, mm_vector_bfgs2)
  TESTS_NCM_ADD_INVALID (gsl, mm_steepest_descent)

  TESTS_NCM_ADD_INVALID (gsl, nmsimplex)
  TESTS_NCM_ADD_INVALID (gsl, nmsimplex2)
  TESTS_NCM_ADD_INVALID (gsl, nmsimplex2rand)
#endif
  g_test_run ();
}

#ifdef NUMCOSMO_HAVE_NLOPT
TESTS_NCM_NEW (nlopt, neldermead, NCM_FIT_TYPE_NLOPT, "ln-neldermead", 20, NCM_FIT_DEFAULT_MAXITER)
TESTS_NCM_NEW (nlopt, slsqp,      NCM_FIT_TYPE_NLOPT, "ld-slsqp",      20, NCM_FIT_DEFAULT_MAXITER)
#endif /* NUMCOSMO_HAVE_NLOPT */

TESTS_NCM_NEW (gsl, ls, NCM_FIT_TYPE_GSL_LS, NULL, 20, 10000000)

TESTS_NCM_NEW (gsl, mm_conjugate_fr,     NCM_FIT_TYPE_GSL_MM, "conjugate-fr",     20, 10000000)
TESTS_NCM_NEW (gsl, mm_conjugate_pr,     NCM_FIT_TYPE_GSL_MM, "conjugate-pr",     20, 10000000)
TESTS_NCM_NEW (gsl, mm_vector_bfgs,      NCM_FIT_TYPE_GSL_MM, "vector-bfgs",      20, 10000000)
TESTS_NCM_NEW (gsl, mm_vector_bfgs2,     NCM_FIT_TYPE_GSL_MM, "vector-bfgs2",     20, 10000000)
TESTS_NCM_NEW (gsl, mm_steepest_descent, NCM_FIT_TYPE_GSL_MM, "steepest-descent", 20, 10000000)

TESTS_NCM_NEW (gsl, nmsimplex, NCM_FIT_TYPE_GSL_MMS,      "nmsimplex",     20, 10000000)
TESTS_NCM_NEW (gsl, nmsimplex2, NCM_FIT_TYPE_GSL_MMS,     "nmsimplex2",     5, 10000000)
TESTS_NCM_NEW (gsl, nmsimplex2rand, NCM_FIT_TYPE_GSL_MMS, "nmsimplex2rand", 5, 10000000)

void
test_ncm_fit_free (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit = test->fit;

  NCM_TEST_FREE (ncm_fit_free, fit);
  NCM_TEST_FREE (ncm_data_free, NCM_DATA (test->data_mvnd));
  NCM_TEST_FREE (ncm_rng_free, test->rng);
}

void
test_ncm_fit_run (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit = test->fit;

  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);
  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

  {
    NcmMSet *mset   = ncm_fit_peek_mset (fit);
    NcmModel *model = NCM_MODEL (ncm_mset_peek (mset, ncm_model_mvnd_id ()));
    NcmVector *ym   = ncm_model_orig_vparam_get_vector (model, NCM_MODEL_MVND_MEAN);
    NcmVector *y    = ncm_data_gauss_cov_mvnd_peek_mean (test->data_mvnd);
    gint i;

    for (i = 0; i < ncm_vector_len (y); i++)
    {
      ncm_assert_cmpdouble_e (ncm_vector_get (y, i), ==, ncm_vector_get (ym, i), 5.0e-2, 5.0e-2);
      /* printf ("[%4d] % 22.15g % 22.15g %e\n", i, ncm_vector_get (y, i), ncm_vector_get (ym, i), fabs (ncm_vector_get (y, i) / ncm_vector_get (ym, i) - 1.0)); */
    }
  }
}

void
test_ncm_fit_run_empty (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit = test->fit;

  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);
  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);
}

#ifdef NUMCOSMO_HAVE_NLOPT
TESTS_NCM_TRAPS (nlopt, neldermead)
TESTS_NCM_TRAPS (nlopt, slsqp)
#endif /* NUMCOSMO_HAVE_NLOPT */

TESTS_NCM_TRAPS (gsl, ls)

TESTS_NCM_TRAPS (gsl, mm_conjugate_fr)
TESTS_NCM_TRAPS (gsl, mm_conjugate_pr)
TESTS_NCM_TRAPS (gsl, mm_vector_bfgs)
TESTS_NCM_TRAPS (gsl, mm_vector_bfgs2)
TESTS_NCM_TRAPS (gsl, mm_steepest_descent)

TESTS_NCM_TRAPS (gsl, nmsimplex)
TESTS_NCM_TRAPS (gsl, nmsimplex2)
TESTS_NCM_TRAPS (gsl, nmsimplex2rand)

void
test_ncm_fit_invalid_run (TestNcmFit *test, gconstpointer pdata)
{
  /*NcmFit *fit = test->fit;*/
  g_assert_not_reached ();
}

