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
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run/simple", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_run_simple, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run/full", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_run_full, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/set_get", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_set_get, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/params/set_get", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_params_set_get, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run/grad/forward", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_run_grad_forward, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run/grad/accurate", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_run_grad_accurate, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run/grad/wrong/type", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_run_grad_wrong_type, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run/empty", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new_empty, \
                    &test_ncm_fit_run_empty, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run/empty/restart", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new_empty, \
                    &test_ncm_fit_run_restart, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run/restart", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_run_restart, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run/restart/save", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_run_restart_save, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/run/restart/save/file", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_run_restart_save_file, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/serialize", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_serialize, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/copy_new", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_copy_new, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/sub_fit/wrong/fit", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_sub_fit_wrong_fit, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/sub_fit/wrong/mset", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_sub_fit_wrong_mset, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/sub_fit/wrong/param", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_sub_fit_wrong_param, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/sub_fit/run", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_sub_fit_run, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/constraints/equality", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_equality_constraints, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/constraints/inequality", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_inequality_constraints, \
                    &test_ncm_fit_free); \
\
        g_test_add ("/ncm/fit/" #lib "/" #algo "/serialize/constraints", TestNcmFit, NULL, \
                    &test_ncm_fit_ ## lib ## _ ## algo ## _new, \
                    &test_ncm_fit_serialize_constraints, \
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

#define TESTS_NCM_NEW(lib, algo, lib_enum, fit_type, algo_str, max_dim, max_iter) \
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
          g_assert_true (NCM_IS_ ## fit_type (fit)); \
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

TESTS_NCM_DECL (levmar, der)
TESTS_NCM_DECL (levmar, dif)
TESTS_NCM_DECL (levmar, bc_dif)
TESTS_NCM_DECL (levmar, bc_der)

void test_ncm_fit_free (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_set_get (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_params_set_get (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run_simple (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run_full (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run_grad_forward (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run_grad_accurate (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run_grad_wrong_type (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run_empty (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run_restart (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run_restart_save (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_run_restart_save_file (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_serialize (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_copy_new (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_sub_fit_wrong_fit (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_sub_fit_wrong_mset (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_sub_fit_wrong_param (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_sub_fit_run (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_equality_constraints (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_inequality_constraints (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_serialize_constraints (TestNcmFit *test, gconstpointer pdata);
void test_ncm_fit_invalid_run (TestNcmFit *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);

  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* g_test_set_nonfatal_assertions (); */

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

  TESTS_NCM_ADD (levmar, der)
  TESTS_NCM_ADD (levmar, dif)
  TESTS_NCM_ADD (levmar, bc_der)
  TESTS_NCM_ADD (levmar, bc_dif)

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

  TESTS_NCM_ADD_INVALID (levmar, der)
  TESTS_NCM_ADD_INVALID (levmar, dif)
  TESTS_NCM_ADD_INVALID (levmar, bc_der)
  TESTS_NCM_ADD_INVALID (levmar, bc_dif)
#endif
  g_test_run ();
}

#ifdef NUMCOSMO_HAVE_NLOPT
TESTS_NCM_NEW (nlopt, neldermead, NCM_FIT_TYPE_NLOPT, FIT_NLOPT, "ln-neldermead", 20, NCM_FIT_DEFAULT_MAXITER)
TESTS_NCM_NEW (nlopt, slsqp,      NCM_FIT_TYPE_NLOPT, FIT_NLOPT, "ld-slsqp",      20, NCM_FIT_DEFAULT_MAXITER)
#endif /* NUMCOSMO_HAVE_NLOPT */

TESTS_NCM_NEW (gsl, ls, NCM_FIT_TYPE_GSL_LS, FIT_GSL_LS, NULL, 20, 10000000)

TESTS_NCM_NEW (gsl, mm_conjugate_fr,     NCM_FIT_TYPE_GSL_MM, FIT_GSL_MM, "conjugate-fr",     20, 10000000)
TESTS_NCM_NEW (gsl, mm_conjugate_pr,     NCM_FIT_TYPE_GSL_MM, FIT_GSL_MM, "conjugate-pr",     20, 10000000)
TESTS_NCM_NEW (gsl, mm_vector_bfgs,      NCM_FIT_TYPE_GSL_MM, FIT_GSL_MM, "vector-bfgs",      20, 10000000)
TESTS_NCM_NEW (gsl, mm_vector_bfgs2,     NCM_FIT_TYPE_GSL_MM, FIT_GSL_MM, "vector-bfgs2",     20, 10000000)
TESTS_NCM_NEW (gsl, mm_steepest_descent, NCM_FIT_TYPE_GSL_MM, FIT_GSL_MM, "steepest-descent", 20, 10000000)

TESTS_NCM_NEW (gsl, nmsimplex,      NCM_FIT_TYPE_GSL_MMS, FIT_GSL_MMS, "nmsimplex",     20, 10000000)
TESTS_NCM_NEW (gsl, nmsimplex2,     NCM_FIT_TYPE_GSL_MMS, FIT_GSL_MMS, "nmsimplex2",     5, 10000000)
TESTS_NCM_NEW (gsl, nmsimplex2rand, NCM_FIT_TYPE_GSL_MMS, FIT_GSL_MMS, "nmsimplex2rand", 5, 10000000)

TESTS_NCM_NEW (levmar, der,    NCM_FIT_TYPE_LEVMAR, FIT_LEVMAR, "der",    20, 10000000)
TESTS_NCM_NEW (levmar, dif,    NCM_FIT_TYPE_LEVMAR, FIT_LEVMAR, "dif",    20, 10000000)
TESTS_NCM_NEW (levmar, bc_der, NCM_FIT_TYPE_LEVMAR, FIT_LEVMAR, "bc-der", 20, 10000000)
TESTS_NCM_NEW (levmar, bc_dif, NCM_FIT_TYPE_LEVMAR, FIT_LEVMAR, "bc-dif", 20, 10000000)

void
test_ncm_fit_free (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit = test->fit;

  NCM_TEST_FREE (ncm_fit_free, fit);
  NCM_TEST_FREE (ncm_data_free, NCM_DATA (test->data_mvnd));
  NCM_TEST_FREE (ncm_rng_free, test->rng);
}

void
test_ncm_fit_set_get (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit                 = test->fit;
  const guint maxiter         = g_test_rand_int_range (1, 1000);
  const gdouble m2lnL_abstol  = g_test_rand_double_range (1.0e-6, 1.0e-3);
  const gdouble m2lnL_reltol  = g_test_rand_double_range (1.0e-6, 1.0e-3);
  const gdouble params_reltol = g_test_rand_double_range (1.0e-6, 1.0e-3);

  ncm_fit_set_maxiter (fit, maxiter);
  g_assert_true (ncm_fit_get_maxiter (fit) == maxiter);

  ncm_fit_set_m2lnL_abstol (fit, m2lnL_abstol);
  g_assert_true (ncm_fit_get_m2lnL_abstol (fit) == m2lnL_abstol);

  ncm_fit_set_m2lnL_reltol (fit, m2lnL_reltol);
  g_assert_true (ncm_fit_get_m2lnL_reltol (fit) == m2lnL_reltol);

  ncm_fit_set_params_reltol (fit, params_reltol);
  g_assert_true (ncm_fit_get_params_reltol (fit) == params_reltol);

  ncm_fit_set_messages (fit, NCM_FIT_RUN_MSGS_NONE);
  g_assert_true (ncm_fit_get_messages (fit) == NCM_FIT_RUN_MSGS_NONE);

  ncm_fit_set_messages (fit, NCM_FIT_RUN_MSGS_SIMPLE);
  g_assert_true (ncm_fit_get_messages (fit) == NCM_FIT_RUN_MSGS_SIMPLE);

  ncm_fit_set_messages (fit, NCM_FIT_RUN_MSGS_FULL);
  g_assert_true (ncm_fit_get_messages (fit) == NCM_FIT_RUN_MSGS_FULL);

  ncm_fit_set_grad_type (fit, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  g_assert_true (ncm_fit_get_grad_type (fit) == NCM_FIT_GRAD_NUMDIFF_FORWARD);

  ncm_fit_set_grad_type (fit, NCM_FIT_GRAD_NUMDIFF_CENTRAL);
  g_assert_true (ncm_fit_get_grad_type (fit) == NCM_FIT_GRAD_NUMDIFF_CENTRAL);

  ncm_fit_set_grad_type (fit, NCM_FIT_GRAD_NUMDIFF_ACCURATE);
  g_assert_true (ncm_fit_get_grad_type (fit) == NCM_FIT_GRAD_NUMDIFF_ACCURATE);

  if (NCM_IS_FIT_LEVMAR (fit) || NCM_IS_FIT_GSL_LS (fit))
    g_assert_true (ncm_fit_is_least_squares (fit));
  else
    g_assert_false (ncm_fit_is_least_squares (fit));

  g_assert_true (NCM_IS_MSET (ncm_fit_peek_mset (fit)));
  g_assert_true (NCM_IS_FIT_STATE (ncm_fit_peek_state (fit)));
  g_assert_true (NCM_IS_LIKELIHOOD (ncm_fit_peek_likelihood (fit)));
  g_assert_true (NCM_IS_DIFF (ncm_fit_peek_diff (fit)));
}

void
test_ncm_fit_params_set_get (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit             = test->fit;
  NcmMSet *mset           = ncm_fit_peek_mset (fit);
  const gdouble x0        = g_test_rand_double_range (-1.0, 1.0);
  const guint fparams_len = ncm_mset_fparam_len (mset);

  ncm_fit_params_set (fit, 0, x0);
  g_assert_true (ncm_mset_fparam_get (mset, 0) == x0);

  /* Testing setting vector */
  {
    NcmVector *x_vec = ncm_vector_new (fparams_len);
    gint i;

    for (i = 0; i < fparams_len; i++)
      ncm_vector_set (x_vec, i, g_test_rand_double_range (-1.0, 1.0));

    ncm_fit_params_set_vector (fit, x_vec);

    for (i = 0; i < fparams_len; i++)
      g_assert_true (ncm_mset_fparam_get (mset, i) == ncm_vector_get (x_vec, i));
  }

  /* Testing setting vector with offset */
  {
    guint offset     = g_test_rand_int_range (1, 10);
    NcmVector *x_vec = ncm_vector_new (fparams_len + offset);
    gint i;

    for (i = 0; i < fparams_len + offset; i++)
      ncm_vector_set (x_vec, i, g_test_rand_double_range (-1.0, 1.0));

    ncm_fit_params_set_vector_offset (fit, x_vec, offset);

    for (i = 0; i < fparams_len; i++)
      g_assert_true (ncm_mset_fparam_get (mset, i) == ncm_vector_get (x_vec, i + offset));
  }

  /* Testing setting array */
  {
    GArray *x_array = g_array_new (FALSE, FALSE, sizeof (gdouble));
    gint i;

    g_array_set_size (x_array, fparams_len);

    for (i = 0; i < fparams_len; i++)
      g_array_index (x_array, gdouble, i) = g_test_rand_double_range (-1.0, 1.0);

    ncm_fit_params_set_array (fit, (gdouble *) x_array->data);

    for (i = 0; i < fparams_len; i++)
      g_assert_true (ncm_mset_fparam_get (mset, i) == g_array_index (x_array, gdouble, i));
  }

  /* Testing setting gsl_vector */
  {
    NcmVector *x_vec = ncm_vector_new (fparams_len);
    gint i;

    for (i = 0; i < fparams_len; i++)
      ncm_vector_set (x_vec, i, g_test_rand_double_range (-1.0, 1.0));

    ncm_fit_params_set_gsl_vector (fit, ncm_vector_gsl (x_vec));

    for (i = 0; i < fparams_len; i++)
      g_assert_true (ncm_mset_fparam_get (mset, i) == ncm_vector_get (x_vec, i));
  }
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

      /*
       *  printf ("[%4d] % 22.15g % 22.15g %e\n", i,
       *       ncm_vector_get (y, i),
       *       ncm_vector_get (ym, i),
       *       fabs (ncm_vector_get (y, i) / ncm_vector_get (ym, i) - 1.0));
       */
    }
  }
}

void
test_ncm_fit_run_simple (TestNcmFit *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    NcmFit *fit = test->fit;

    ncm_fit_run (fit, NCM_FIT_RUN_MSGS_SIMPLE);
    ncm_fit_run (fit, NCM_FIT_RUN_MSGS_SIMPLE);

    {
      NcmMSet *mset   = ncm_fit_peek_mset (fit);
      NcmModel *model = NCM_MODEL (ncm_mset_peek (mset, ncm_model_mvnd_id ()));
      NcmVector *ym   = ncm_model_orig_vparam_get_vector (model, NCM_MODEL_MVND_MEAN);
      NcmVector *y    = ncm_data_gauss_cov_mvnd_peek_mean (test->data_mvnd);
      gint i;

      for (i = 0; i < ncm_vector_len (y); i++)
      {
        ncm_assert_cmpdouble_e (ncm_vector_get (y, i), ==, ncm_vector_get (ym, i), 5.0e-2, 5.0e-2);

        /*
         *  printf ("[%4d] % 22.15g % 22.15g %e\n", i,
         *       ncm_vector_get (y, i),
         *       ncm_vector_get (ym, i),
         *       fabs (ncm_vector_get (y, i) / ncm_vector_get (ym, i) - 1.0));
         */
      }
    }

    return;
  }

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, G_TEST_SUBPROCESS_INHERIT_STDOUT);
  g_test_trap_assert_passed ();
  g_test_trap_assert_stdout ("*Minimum found with precision*");
}

void
test_ncm_fit_run_full (TestNcmFit *test, gconstpointer pdata)
{
  if (g_test_subprocess ())
  {
    NcmFit *fit = test->fit;

    ncm_fit_run (fit, NCM_FIT_RUN_MSGS_FULL);
    ncm_fit_run (fit, NCM_FIT_RUN_MSGS_FULL);

    {
      NcmMSet *mset   = ncm_fit_peek_mset (fit);
      NcmModel *model = NCM_MODEL (ncm_mset_peek (mset, ncm_model_mvnd_id ()));
      NcmVector *ym   = ncm_model_orig_vparam_get_vector (model, NCM_MODEL_MVND_MEAN);
      NcmVector *y    = ncm_data_gauss_cov_mvnd_peek_mean (test->data_mvnd);
      gint i;

      for (i = 0; i < ncm_vector_len (y); i++)
      {
        ncm_assert_cmpdouble_e (ncm_vector_get (y, i), ==, ncm_vector_get (ym, i), 5.0e-2, 5.0e-2);

        /*
         *  printf ("[%4d] % 22.15g % 22.15g %e\n", i,
         *       ncm_vector_get (y, i),
         *       ncm_vector_get (ym, i),
         *       fabs (ncm_vector_get (y, i) / ncm_vector_get (ym, i) - 1.0));
         */
      }
    }

    return;
  }

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, G_TEST_SUBPROCESS_INHERIT_STDOUT);
  g_test_trap_assert_passed ();
  g_test_trap_assert_stdout ("*Minimum found with precision*");
}

void
test_ncm_fit_run_grad_forward (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit = test->fit;

  ncm_fit_set_grad_type (fit, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  g_assert_true (ncm_fit_get_grad_type (fit) == NCM_FIT_GRAD_NUMDIFF_FORWARD);

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

      /*
       *  printf ("[%4d] % 22.15g % 22.15g %e\n", i,
       *       ncm_vector_get (y, i),
       *       ncm_vector_get (ym, i),
       *       fabs (ncm_vector_get (y, i) / ncm_vector_get (ym, i) - 1.0));
       */
    }
  }
}

void
test_ncm_fit_run_grad_accurate (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit = test->fit;

  ncm_fit_set_grad_type (fit, NCM_FIT_GRAD_NUMDIFF_ACCURATE);
  g_assert_true (ncm_fit_get_grad_type (fit) == NCM_FIT_GRAD_NUMDIFF_ACCURATE);

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

      /*
       *  printf ("[%4d] % 22.15g % 22.15g %e\n", i,
       *       ncm_vector_get (y, i),
       *       ncm_vector_get (ym, i),
       *       fabs (ncm_vector_get (y, i) / ncm_vector_get (ym, i) - 1.0));
       */
    }
  }
}

void
test_ncm_fit_run_grad_wrong_type (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit = test->fit;

  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_fit_set_grad_type (fit, 1000);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fit_run_empty (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit = test->fit;

  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);
  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);
}

void
test_ncm_fit_run_restart (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit = test->fit;

  ncm_fit_run_restart (fit, NCM_FIT_RUN_MSGS_NONE, 1.0e-3, 0.0, NULL, NULL);
}

void
test_ncm_fit_run_restart_save (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit       = test->fit;
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmMSet *mset     = ncm_fit_peek_mset (fit);
  NcmMSet *mset_dup = ncm_mset_dup (mset, ser);

  ncm_fit_run_restart (fit, NCM_FIT_RUN_MSGS_NONE, 1.0e-3, 0.0, mset_dup, NULL);

  g_assert_true (ncm_mset_cmp (mset, mset_dup, TRUE));

  {
    const guint fparams_len = ncm_mset_fparams_len (mset);
    NcmVector *x            = ncm_vector_new (fparams_len);
    NcmVector *x_dup        = ncm_vector_new (fparams_len);

    g_assert_true (fparams_len == ncm_mset_fparams_len (mset_dup));

    ncm_mset_fparams_get_vector (mset, x);
    ncm_mset_fparams_get_vector (mset_dup, x_dup);

    g_assert_true (ncm_vector_cmp2 (x, x_dup, 1.0e-3, 0.0) == 0);
  }

  ncm_mset_free (mset_dup);
  ncm_serialize_free (ser);
}

void
test_ncm_fit_run_restart_save_file (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit       = test->fit;
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmMSet *mset     = ncm_fit_peek_mset (fit);
  NcmMSet *mset_dup;

  ncm_fit_run_restart (fit, NCM_FIT_RUN_MSGS_NONE, 1.0e-3, 0.0, NULL, "tmp_test_ncm_fit_run_restart_save_file.mset");

  mset_dup = ncm_mset_load ("tmp_test_ncm_fit_run_restart_save_file.mset", ser);

  g_assert_true (ncm_mset_cmp (mset, mset_dup, TRUE));

  {
    const guint fparams_len = ncm_mset_fparams_len (mset);
    NcmVector *x            = ncm_vector_new (fparams_len);
    NcmVector *x_dup        = ncm_vector_new (fparams_len);

    g_assert_true (fparams_len == ncm_mset_fparams_len (mset_dup));

    ncm_mset_fparams_get_vector (mset, x);
    ncm_mset_fparams_get_vector (mset_dup, x_dup);

    g_assert_true (ncm_vector_cmp2 (x, x_dup, 1.0e-3, 0.0) == 0);
  }

  ncm_mset_free (mset_dup);
  ncm_serialize_free (ser);
}

void
test_ncm_fit_serialize (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit       = test->fit;
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmFit *fit_dup   = NCM_FIT (ncm_serialize_dup_obj (ser, G_OBJECT (fit)));

  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);
  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

  ncm_fit_run (fit_dup, NCM_FIT_RUN_MSGS_NONE);
  ncm_fit_run (fit_dup, NCM_FIT_RUN_MSGS_NONE);

  {
    NcmMSet *mset       = ncm_fit_peek_mset (fit);
    NcmMSet *mset_dup   = ncm_fit_peek_mset (fit_dup);
    NcmModel *model     = NCM_MODEL (ncm_mset_peek (mset, ncm_model_mvnd_id ()));
    NcmModel *model_dup = NCM_MODEL (ncm_mset_peek (mset_dup, ncm_model_mvnd_id ()));
    NcmVector *y        = ncm_model_orig_vparam_get_vector (model, NCM_MODEL_MVND_MEAN);
    NcmVector *y_dup    = ncm_model_orig_vparam_get_vector (model_dup, NCM_MODEL_MVND_MEAN);
    gint i;

    g_assert_true (NCM_IS_FIT (fit_dup));
    g_assert_true (mset != mset_dup);
    g_assert_true (model != model_dup);
    g_assert_true (y != y_dup);

    for (i = 0; i < ncm_vector_len (y); i++)
    {
      ncm_assert_cmpdouble_e (ncm_vector_get (y, i), ==, ncm_vector_get (y_dup, i), 5.0e-2, 5.0e-2);

      /*
       *  printf ("[%4d] % 22.15g % 22.15g %e\n", i,
       *       ncm_vector_get (y, i),
       *       ncm_vector_get (y_dup, i),
       *       fabs (ncm_vector_get (y, i) / ncm_vector_get (y_dup, i) - 1.0));
       */
    }
  }

  ncm_serialize_free (ser);
}

void
test_ncm_fit_copy_new (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit     = test->fit;
  NcmFit *fit_dup = ncm_fit_copy_new (fit,
                                      ncm_fit_peek_likelihood (fit),
                                      ncm_fit_peek_mset (fit),
                                      ncm_fit_get_grad_type (fit));

  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);
  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

  ncm_fit_run (fit_dup, NCM_FIT_RUN_MSGS_NONE);
  ncm_fit_run (fit_dup, NCM_FIT_RUN_MSGS_NONE);

  {
    NcmMSet *mset       = ncm_fit_peek_mset (fit);
    NcmMSet *mset_dup   = ncm_fit_peek_mset (fit_dup);
    NcmModel *model     = NCM_MODEL (ncm_mset_peek (mset, ncm_model_mvnd_id ()));
    NcmModel *model_dup = NCM_MODEL (ncm_mset_peek (mset_dup, ncm_model_mvnd_id ()));
    NcmVector *y        = ncm_model_orig_vparam_get_vector (model, NCM_MODEL_MVND_MEAN);
    NcmVector *y_dup    = ncm_model_orig_vparam_get_vector (model_dup, NCM_MODEL_MVND_MEAN);
    gint i;

    g_assert_true (NCM_IS_FIT (fit_dup));
    g_assert_true (mset == mset_dup);
    g_assert_true (model == model_dup);

    for (i = 0; i < ncm_vector_len (y); i++)
    {
      ncm_assert_cmpdouble_e (ncm_vector_get (y, i), ==, ncm_vector_get (y_dup, i), 5.0e-2, 5.0e-2);

      /*
       *  printf ("[%4d] % 22.15g % 22.15g %e\n", i,
       *       ncm_vector_get (y, i),
       *       ncm_vector_get (y_dup, i),
       *       fabs (ncm_vector_get (y, i) / ncm_vector_get (y_dup, i) - 1.0));
       */
    }
  }

  ncm_fit_free (fit_dup);
}

void
test_ncm_fit_sub_fit_wrong_fit (TestNcmFit *test, gconstpointer pdata)
{
  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_fit_set_sub_fit (test->fit, test->fit);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fit_sub_fit_wrong_mset (TestNcmFit *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmFit *fit_dup   = NCM_FIT (ncm_serialize_dup_obj (ser, G_OBJECT (test->fit)));

  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_mset_remove (ncm_fit_peek_mset (fit_dup), ncm_model_mvnd_id ());
    ncm_fit_set_sub_fit (fit_dup, test->fit);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();

  ncm_serialize_free (ser);
  ncm_fit_free (fit_dup);
}

void
test_ncm_fit_sub_fit_wrong_param (TestNcmFit *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmFit *fit_dup   = NCM_FIT (ncm_serialize_dup_obj (ser, G_OBJECT (test->fit)));

  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_fit_set_sub_fit (fit_dup, test->fit);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();

  ncm_serialize_free (ser);
  ncm_fit_free (fit_dup);
}

void
test_ncm_fit_sub_fit_run (TestNcmFit *test, gconstpointer pdata)
{
  NcmMSet *mset     = ncm_fit_peek_mset (test->fit);
  NcmMSet *mset_dup = ncm_mset_shallow_copy (mset);
  NcmFit *fit_dup   = ncm_fit_copy_new (test->fit,
                                        ncm_fit_peek_likelihood (test->fit),
                                        mset_dup,
                                        ncm_fit_get_grad_type (test->fit));

  ncm_mset_param_set_ftype (mset, ncm_model_mvnd_id (), 0, NCM_PARAM_TYPE_FIXED);

  ncm_mset_param_set_all_ftype (mset_dup, NCM_PARAM_TYPE_FIXED);
  ncm_mset_param_set_ftype (mset_dup, ncm_model_mvnd_id (), 0, NCM_PARAM_TYPE_FREE);

  ncm_fit_set_sub_fit (test->fit, fit_dup);
  ncm_fit_run (test->fit, NCM_FIT_RUN_MSGS_NONE);

  g_assert_true (ncm_fit_get_sub_fit (test->fit) == fit_dup);

  {
    NcmMSet *mset       = ncm_fit_peek_mset (test->fit);
    NcmMSet *mset_dup   = ncm_fit_peek_mset (fit_dup);
    NcmModel *model     = NCM_MODEL (ncm_mset_peek (mset, ncm_model_mvnd_id ()));
    NcmModel *model_dup = NCM_MODEL (ncm_mset_peek (mset_dup, ncm_model_mvnd_id ()));
    NcmVector *y        = ncm_model_orig_vparam_get_vector (model, NCM_MODEL_MVND_MEAN);
    NcmVector *y_dup    = ncm_model_orig_vparam_get_vector (model_dup, NCM_MODEL_MVND_MEAN);
    gint i;

    g_assert_true (NCM_IS_FIT (fit_dup));
    g_assert_false (mset == mset_dup);
    g_assert_true (model == model_dup);

    for (i = 0; i < ncm_vector_len (y); i++)
    {
      ncm_assert_cmpdouble_e (ncm_vector_get (y, i), ==, ncm_vector_get (y_dup, i), 5.0e-2, 5.0e-2);

      /*
       *  printf ("[%4d] % 22.15g % 22.15g %e\n", i,
       *       ncm_vector_get (y, i),
       *       ncm_vector_get (y_dup, i),
       *       fabs (ncm_vector_get (y, i) / ncm_vector_get (y_dup, i) - 1.0));
       */
    }
  }

  ncm_fit_free (fit_dup);
  ncm_fit_free (fit_dup);
}

void
test_ncm_fit_equality_constraints (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit       = test->fit;
  NcmMSetFunc *func = NCM_MSET_FUNC (ncm_prior_gauss_param_new (ncm_model_mvnd_id (), 0, 1.0, 1.0));

  ncm_fit_add_equality_constraint (fit, func, 1.0e-5);

#ifdef NUMCOSMO_HAVE_NLOPT

  if (NCM_IS_FIT_NLOPT (fit))
  {
    NcmFitNloptAlgorithm algorithm;
    g_object_get (fit, "algorithm", &algorithm, NULL);

    if (algorithm != NCM_FIT_NLOPT_LD_SLSQP)
    {
      g_test_skip ("Only NLOpt:slsqp supports equality constraints");

      return;
    }

    ncm_fit_set_grad_type (fit, NCM_FIT_GRAD_NUMDIFF_CENTRAL);
  }
  else
  {
    g_test_skip ("Only NLOpt:slsqp supports equality constraints");

    return;
  }

  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

  {
    NcmMSet *mset   = ncm_fit_peek_mset (fit);
    NcmModel *model = NCM_MODEL (ncm_mset_peek (mset, ncm_model_mvnd_id ()));
    NcmVector *ym   = ncm_model_orig_vparam_get_vector (model, NCM_MODEL_MVND_MEAN);

    ncm_assert_cmpdouble_e (ncm_vector_get (ym, 0), ==, 1.0, 5.0e-2, 5.0e-2);
  }

  g_assert_true (ncm_fit_equality_constraints_len (fit) == 1);
  {
    NcmMSetFunc *func2 = NULL;
    gdouble tol;
    ncm_fit_get_equality_constraint (fit, 0, &func2, &tol);

    g_assert_true (func == func2);
    g_assert_true (tol == 1.0e-5);
  }
  ncm_fit_remove_equality_constraints (fit);
  g_assert_true (ncm_fit_equality_constraints_len (fit) == 0);

#endif /* NUMCOSMO_HAVE_NLOPT */
  ncm_mset_func_free (func);
}

void
test_ncm_fit_inequality_constraints (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit       = test->fit;
  NcmMSetFunc *func = NCM_MSET_FUNC (ncm_prior_gauss_param_new (ncm_model_mvnd_id (), 0, 1.0, 1.0));

  ncm_fit_add_inequality_constraint (fit, func, 1.0e-5);

#ifdef NUMCOSMO_HAVE_NLOPT

  if (NCM_IS_FIT_NLOPT (fit))
  {
    NcmFitNloptAlgorithm algorithm;
    g_object_get (fit, "algorithm", &algorithm, NULL);

    if (algorithm != NCM_FIT_NLOPT_LD_SLSQP)
    {
      g_test_skip ("Only NLOpt:slsqp supports equality constraints");

      return;
    }

    ncm_fit_set_grad_type (fit, NCM_FIT_GRAD_NUMDIFF_CENTRAL);
  }
  else
  {
    g_test_skip ("Only NLOpt:slsqp supports equality constraints");

    return;
  }

  ncm_fit_run (fit, NCM_FIT_RUN_MSGS_NONE);

  {
    NcmMSet *mset   = ncm_fit_peek_mset (fit);
    NcmModel *model = NCM_MODEL (ncm_mset_peek (mset, ncm_model_mvnd_id ()));
    NcmVector *ym   = ncm_model_orig_vparam_get_vector (model, NCM_MODEL_MVND_MEAN);

    ncm_assert_cmpdouble_e (ncm_vector_get (ym, 0), ==, 1.0, 5.0e-2, 5.0e-2);
  }

  g_assert_true (ncm_fit_inequality_constraints_len (fit) == 1);
  {
    NcmMSetFunc *func2 = NULL;
    gdouble tol;
    ncm_fit_get_inequality_constraint (fit, 0, &func2, &tol);

    g_assert_true (func == func2);
    g_assert_true (tol == 1.0e-5);
  }
  ncm_fit_remove_inequality_constraints (fit);
  g_assert_true (ncm_fit_inequality_constraints_len (fit) == 0);

#endif /* NUMCOSMO_HAVE_NLOPT */
  ncm_mset_func_free (func);
}

void
test_ncm_fit_serialize_constraints (TestNcmFit *test, gconstpointer pdata)
{
  NcmFit *fit        = test->fit;
  NcmSerialize *ser  = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmMSetFunc *func0 = NCM_MSET_FUNC (ncm_prior_gauss_param_new (ncm_model_mvnd_id (), 0, 1.0, 1.0));
  NcmMSetFunc *func1 = NCM_MSET_FUNC (ncm_prior_gauss_param_new (ncm_model_mvnd_id (), 1, 1.1, 1.1));
  NcmFit *fit_dup;

  ncm_fit_add_inequality_constraint (fit, func0, 1.0e-5);
  ncm_fit_add_equality_constraint (fit, func1, 1.0e-5);

  fit_dup = ncm_fit_dup (fit, ser);

  g_assert_true (ncm_fit_inequality_constraints_len (fit_dup) == 1);
  g_assert_true (ncm_fit_equality_constraints_len (fit_dup) == 1);

  {
    NcmMSetFunc *func0_dup = NULL;
    NcmMSetFunc *func1_dup = NULL;
    gdouble mu, sigma;
    guint mid, pid;
    gdouble tol0_dup;
    gdouble tol1_dup;

    ncm_fit_get_inequality_constraint (fit_dup, 0, &func0_dup, &tol0_dup);
    ncm_fit_get_equality_constraint (fit_dup, 0, &func1_dup, &tol1_dup);

    g_assert_true (NCM_IS_MSET_FUNC (func0_dup));
    g_assert_true (NCM_IS_MSET_FUNC (func1_dup));

    g_assert_true (NCM_IS_PRIOR_GAUSS (func0_dup));
    g_assert_true (NCM_IS_PRIOR_GAUSS (func1_dup));

    g_assert_true (NCM_IS_PRIOR_GAUSS_PARAM (func0_dup));
    g_assert_true (NCM_IS_PRIOR_GAUSS_PARAM (func1_dup));

    g_object_get (func0_dup,
                  "mu", &mu,
                  "sigma", &sigma,
                  "mid", &mid,
                  "pid", &pid,
                  NULL);

    g_assert_true (mu == 1.0);
    g_assert_true (sigma == 1.0);
    g_assert_true (mid == ncm_model_mvnd_id ());
    g_assert_true (pid == 0);

    g_object_get (func1_dup,
                  "mu", &mu,
                  "sigma", &sigma,
                  "mid", &mid,
                  "pid", &pid,
                  NULL);

    g_assert_true (mu == 1.1);
    g_assert_true (sigma == 1.1);
    g_assert_true (mid == ncm_model_mvnd_id ());
    g_assert_true (pid == 1);

    g_assert_true (tol0_dup == 1.0e-5);
    g_assert_true (tol1_dup == 1.0e-5);
  }

  g_object_set (fit_dup,
                "equality-constraints", NULL,
                "equality-constraints-tot", NULL,
                "inequality-constraints", NULL,
                "inequality-constraints-tot", NULL,
                NULL);

  g_assert_true (ncm_fit_inequality_constraints_len (fit_dup) == 0);
  g_assert_true (ncm_fit_equality_constraints_len (fit_dup) == 0);
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

TESTS_NCM_TRAPS (levmar, der)
TESTS_NCM_TRAPS (levmar, dif)
TESTS_NCM_TRAPS (levmar, bc_der)
TESTS_NCM_TRAPS (levmar, bc_dif)

void
test_ncm_fit_invalid_run (TestNcmFit *test, gconstpointer pdata)
{
  /*NcmFit *fit = test->fit;*/
  g_assert_not_reached ();
}

