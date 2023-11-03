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
void test_ncm_fit_run (TestNcmFit *test, gconstpointer pdata);
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
void test_ncm_fit_invalid_run (TestNcmFit *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);

  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

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

  ncm_mset_param_set_all_ftype (mset_dup, NCM_PARAM_TYPE_FIXED);
  ncm_mset_param_set_ftype (mset_dup, ncm_model_mvnd_id (), 0, NCM_PARAM_TYPE_FREE);
  ncm_mset_param_set_ftype (mset, ncm_model_mvnd_id (), 0, NCM_PARAM_TYPE_FIXED);

  ncm_fit_set_sub_fit (test->fit, fit_dup);

  ncm_fit_run (test->fit, NCM_FIT_RUN_MSGS_NONE);

  ncm_fit_free (fit_dup);
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

