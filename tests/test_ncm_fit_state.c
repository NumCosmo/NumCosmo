/***************************************************************************
 *            test_ncm_fit_state.c
 *
 *  Thu November 29 15:27:17 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) 2012 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

typedef struct _TestNcmFitState
{
  NcmFitState *fit_state;
} TestNcmFitState;

void test_ncm_fit_state_new (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_free (TestNcmFitState *test, gconstpointer pdata);

void test_ncm_fit_state_sanity (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_serialize (TestNcmFitState *test, gconstpointer pdata);

void test_ncm_fit_state_set_ls (TestNcmFitState *test, gconstpointer pdata);

void test_ncm_fit_state_fparam_len (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_data_len (TestNcmFitState *test, gconstpointer pdata);

void test_ncm_fit_state_dof (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_niter (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_func_eval (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_grad_eval (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_m2lnL_prec (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_m2lnL_curval (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_params_prec (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_elapsed_time (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_has_covar (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_is_best_fit (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_is_least_squares (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_fparams (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_hessian (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_covar (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_f (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_J (TestNcmFitState *test, gconstpointer pdata);

void test_ncm_fit_state_set_ls_wrong (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_set_m2lnL_prec_wrong_lower_bound (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_set_m2lnL_prec_wrong_upper_bound (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_set_params_prec_wrong_lower_bound (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_set_params_prec_wrong_upper_bound (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_set_elapsed_time_wrong_lower_bound (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_peek_f_wrong (TestNcmFitState *test, gconstpointer pdata);
void test_ncm_fit_state_peek_J_wrong (TestNcmFitState *test, gconstpointer pdata);

int
main (int argc, char *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/fit_state/sanity", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_sanity,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/serialize", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_serialize,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/set_ls", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_set_ls,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/fparam_len", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_fparam_len,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/data_len", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_data_len,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/dof", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_dof,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/niter", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_niter,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/func_eval", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_func_eval,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/grad_eval", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_grad_eval,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/m2lnL_prec", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_m2lnL_prec,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/m2lnL_curval", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_m2lnL_curval,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/params_prec", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_params_prec,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/elapsed_time", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_elapsed_time,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/has_covar", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_has_covar,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/is_best_fit", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_is_best_fit,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/is_least_squares", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_is_least_squares,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/fparams", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_fparams,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/hessian", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_hessian,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/covar", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_covar,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/f", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_f,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/J", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_J,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/set_ls/wrong", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_set_ls_wrong,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/set_m2lnL_prec/wrong/lower_bound", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_set_m2lnL_prec_wrong_lower_bound,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/set_m2lnL_prec/wrong/upper_bound", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_set_m2lnL_prec_wrong_upper_bound,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/set_params_prec/wrong/lower_bound", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_set_params_prec_wrong_lower_bound,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/set_params_prec/wrong/upper_bound", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_set_params_prec_wrong_upper_bound,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/set_elapsed_time/wrong/lower_bound", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_set_elapsed_time_wrong_lower_bound,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/peek_f/wrong", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_peek_f_wrong,
              &test_ncm_fit_state_free);

  g_test_add ("/ncm/fit_state/peek_J/wrong", TestNcmFitState, NULL,
              &test_ncm_fit_state_new,
              &test_ncm_fit_state_peek_J_wrong,
              &test_ncm_fit_state_free);

  return g_test_run ();
}

void
test_ncm_fit_state_new (TestNcmFitState *test, gconstpointer pdata)
{
  NcmFitState *fit_state = ncm_fit_state_new (g_test_rand_int_range (1, 100),
                                              g_test_rand_int_range (1, 100),
                                              g_test_rand_int_range (1, 100),
                                              g_test_rand_int_range (0, 1));

  g_assert_true (NCM_IS_FIT_STATE (fit_state));

  test->fit_state = fit_state;
}

void
test_ncm_fit_state_free (TestNcmFitState *test, gconstpointer pdata)
{
  NcmFitState *fit_state = test->fit_state;

  NCM_TEST_FREE (ncm_fit_state_free, fit_state);
}

void
test_ncm_fit_state_sanity (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  g_assert_true (ncm_fit_state_ref (test->fit_state) == test->fit_state);
  ncm_fit_state_free (test->fit_state);

  {
    NcmFitState *fit_state = ncm_fit_state_ref (test->fit_state);

    ncm_fit_state_clear (&fit_state);

    g_assert_true (fit_state == NULL);
  }
}

void
test_ncm_fit_state_serialize (TestNcmFitState *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);

  NcmFitState *dup_fit_state = NCM_FIT_STATE (ncm_serialize_dup_obj (ser, G_OBJECT (test->fit_state)));

  g_assert_true (NCM_IS_FIT_STATE (dup_fit_state));
  g_assert_true (ncm_fit_state_get_data_len (dup_fit_state) == ncm_fit_state_get_data_len (test->fit_state));
  g_assert_true (ncm_fit_state_get_fparam_len (dup_fit_state) == ncm_fit_state_get_fparam_len (test->fit_state));
  g_assert_true (ncm_fit_state_get_dof (dup_fit_state) == ncm_fit_state_get_dof (test->fit_state));
  g_assert_true (ncm_fit_state_is_least_squares (dup_fit_state) == ncm_fit_state_is_least_squares (test->fit_state));
  g_assert_true (ncm_fit_state_get_niter (dup_fit_state) == ncm_fit_state_get_niter (test->fit_state));
  g_assert_true (ncm_fit_state_get_func_eval (dup_fit_state) == ncm_fit_state_get_func_eval (test->fit_state));
  g_assert_true (ncm_fit_state_get_grad_eval (dup_fit_state) == ncm_fit_state_get_grad_eval (test->fit_state));
  g_assert_true (ncm_fit_state_is_best_fit (dup_fit_state) == ncm_fit_state_is_best_fit (test->fit_state));

  ncm_fit_state_free (dup_fit_state);
  ncm_serialize_free (ser);
}

void
test_ncm_fit_state_set_ls (TestNcmFitState *test, gconstpointer pdata)
{
  const guint fparam_len = g_test_rand_int_range (1, 100);
  const guint data_len   = g_test_rand_int_range (1, 100);

  NcmVector *f = ncm_vector_new (data_len);
  NcmMatrix *J = ncm_matrix_new (data_len, fparam_len);

  g_assert_true (f != NULL);

  ncm_fit_state_set_fparam_len (test->fit_state, fparam_len);
  ncm_fit_state_set_data_len (test->fit_state, data_len);
  ncm_fit_state_set_is_least_squares (test->fit_state, TRUE);

  ncm_fit_state_set_ls (test->fit_state, f, J);

  {
    NcmVector *f_ = ncm_fit_state_peek_f (test->fit_state);
    NcmMatrix *J_ = ncm_fit_state_peek_J (test->fit_state);

    g_assert_true (f_ != NULL);
    g_assert_true (J_ != NULL);

    g_assert_true (NCM_IS_VECTOR (f_));
    g_assert_true (NCM_IS_MATRIX (J_));

    g_assert_cmpuint (ncm_vector_len (f_), ==, data_len);
    g_assert_cmpuint (ncm_matrix_nrows (J_), ==, data_len);
    g_assert_cmpuint (ncm_matrix_ncols (J_), ==, fparam_len);

    g_assert_cmpuint (ncm_vector_len (f_), ==, ncm_fit_state_get_data_len (test->fit_state));
    g_assert_cmpuint (ncm_matrix_nrows (J_), ==, ncm_fit_state_get_data_len (test->fit_state));
    g_assert_cmpuint (ncm_matrix_ncols (J_), ==, ncm_fit_state_get_fparam_len (test->fit_state));

    g_assert_cmpuint (ncm_vector_len (f_), ==, ncm_vector_len (f));
    g_assert_cmpuint (ncm_matrix_nrows (J_), ==, ncm_matrix_nrows (J));
    g_assert_cmpuint (ncm_matrix_ncols (J_), ==, ncm_matrix_ncols (J));

    g_assert_cmpuint (ncm_vector_len (f_), ==, ncm_vector_len (f));
    g_assert_cmpuint (ncm_matrix_nrows (J_), ==, ncm_matrix_nrows (J));
    g_assert_cmpuint (ncm_matrix_ncols (J_), ==, ncm_matrix_ncols (J));

    g_assert_cmpfloat (ncm_vector_get (f_, 0), ==, ncm_vector_get (f, 0));
    g_assert_cmpfloat (ncm_matrix_get (J_, 0, 0), ==, ncm_matrix_get (J, 0, 0));
  }
}

void
test_ncm_fit_state_fparam_len (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    /* Set the fparam_len to a random integer between 1 and 100 */
    guint expected_fparam_len = g_test_rand_int_range (1, 100);
    guint actual_fparam_len;

    ncm_fit_state_set_fparam_len (test->fit_state, expected_fparam_len);
    /* Check if the value was set correctly */
    actual_fparam_len = ncm_fit_state_get_fparam_len (test->fit_state);

    g_assert_cmpuint (actual_fparam_len, ==, expected_fparam_len);
  }
}

void
test_ncm_fit_state_data_len (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    /* Set the data_len to a random integer between 1 and 100 */
    guint expected_data_len = g_test_rand_int_range (1, 100);
    guint actual_data_len;

    ncm_fit_state_set_data_len (test->fit_state, expected_data_len);
    /* Check if the value was set correctly */
    actual_data_len = ncm_fit_state_get_data_len (test->fit_state);

    g_assert_cmpuint (actual_data_len, ==, expected_data_len);
  }
}

void
test_ncm_fit_state_dof (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    /* Set the dof to a random integer between 1 and 100 */
    gint expected_dof = g_test_rand_int_range (1, 100);
    gint actual_dof;

    ncm_fit_state_set_dof (test->fit_state, expected_dof);
    /* Check if the value was set correctly */
    actual_dof = ncm_fit_state_get_dof (test->fit_state);

    g_assert_cmpint (actual_dof, ==, expected_dof);
  }
}

void
test_ncm_fit_state_niter (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    /* Set the niter to a random integer between 1 and 100 */
    guint expected_niter = g_test_rand_int_range (1, 100);
    guint actual_niter;

    ncm_fit_state_set_niter (test->fit_state, expected_niter);
    /* Check if the value was set correctly */
    actual_niter = ncm_fit_state_get_niter (test->fit_state);

    g_assert_cmpuint (actual_niter, ==, expected_niter);

    /* Add 2 to the niter */
    expected_niter += 2;
    ncm_fit_state_add_iter (test->fit_state, 2);
    /* Check if the value was set correctly */
    actual_niter = ncm_fit_state_get_niter (test->fit_state);

    g_assert_cmpuint (actual_niter, ==, expected_niter);
  }
}

void
test_ncm_fit_state_func_eval (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    /* Set the func_eval to a random integer between 1 and 100 */
    guint expected_func_eval = g_test_rand_int_range (1, 100);
    guint actual_func_eval;

    ncm_fit_state_set_func_eval (test->fit_state, expected_func_eval);
    /* Check if the value was set correctly */
    actual_func_eval = ncm_fit_state_get_func_eval (test->fit_state);

    g_assert_cmpuint (actual_func_eval, ==, expected_func_eval);

    /* Add 2 to the func_eval */
    expected_func_eval += 2;
    ncm_fit_state_add_func_eval (test->fit_state, 2);
    /* Check if the value was set correctly */
    actual_func_eval = ncm_fit_state_get_func_eval (test->fit_state);

    g_assert_cmpuint (actual_func_eval, ==, expected_func_eval);
  }
}

void
test_ncm_fit_state_grad_eval (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    /* Set the grad_eval to a random integer between 1 and 100 */
    guint expected_grad_eval = g_test_rand_int_range (1, 100);
    guint actual_grad_eval;

    ncm_fit_state_set_grad_eval (test->fit_state, expected_grad_eval);
    /* Check if the value was set correctly */
    actual_grad_eval = ncm_fit_state_get_grad_eval (test->fit_state);

    g_assert_cmpuint (actual_grad_eval, ==, expected_grad_eval);

    /* Add 2 to the grad_eval */
    expected_grad_eval += 2;
    ncm_fit_state_add_grad_eval (test->fit_state, 2);
    /* Check if the value was set correctly */
    actual_grad_eval = ncm_fit_state_get_grad_eval (test->fit_state);

    g_assert_cmpuint (actual_grad_eval, ==, expected_grad_eval);
  }
}

void
test_ncm_fit_state_m2lnL_prec (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    /* Set the m2lnL_prec to a random double between 1.0 and 100.0 */
    gdouble expected_m2lnL_prec = g_test_rand_double_range (0.0, 1.0);
    gdouble actual_m2lnL_prec;

    ncm_fit_state_set_m2lnL_prec (test->fit_state, expected_m2lnL_prec);
    /* Check if the value was set correctly */
    actual_m2lnL_prec = ncm_fit_state_get_m2lnL_prec (test->fit_state);

    g_assert_cmpfloat (actual_m2lnL_prec, ==, expected_m2lnL_prec);
  }
}

void
test_ncm_fit_state_m2lnL_curval (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    /* Set the m2lnL_curval to a random double between -100.0 and 100.0 */
    gdouble expected_m2lnL_curval = g_test_rand_double_range (-100.0, +100.0);
    gdouble actual_m2lnL_curval;

    ncm_fit_state_set_m2lnL_curval (test->fit_state, expected_m2lnL_curval);
    /* Check if the value was set correctly */
    actual_m2lnL_curval = ncm_fit_state_get_m2lnL_curval (test->fit_state);

    g_assert_cmpfloat (actual_m2lnL_curval, ==, expected_m2lnL_curval);
  }
}

void
test_ncm_fit_state_params_prec (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    /* Set the params_prec to random double between 1.0 and 100.0 */
    gdouble expected_params_prec = g_test_rand_double_range (0.0, 1.0);

    gdouble actual_params_prec;

    ncm_fit_state_set_params_prec (test->fit_state, expected_params_prec);
    /* Check if the value was set correctly */
    actual_params_prec = ncm_fit_state_get_params_prec (test->fit_state);

    g_assert_cmpfloat (actual_params_prec, ==, expected_params_prec);
  }
}

void
test_ncm_fit_state_elapsed_time (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    /* Set the elapsed_time to random double between 1.0 and 100.0 */
    gdouble expected_elapsed_time = g_test_rand_double_range (1.0, 100.0);
    gdouble actual_elapsed_time;

    ncm_fit_state_set_elapsed_time (test->fit_state, expected_elapsed_time);
    /* Check if the value was set correctly */
    actual_elapsed_time = ncm_fit_state_get_elapsed_time (test->fit_state);

    g_assert_cmpfloat (actual_elapsed_time, ==, expected_elapsed_time);
  }
}

void
test_ncm_fit_state_has_covar (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    gboolean expected_has_covar = g_test_rand_int_range (0, 1);
    gboolean actual_has_covar;

    ncm_fit_state_set_has_covar (test->fit_state, expected_has_covar);
    /* Check if the value was set correctly */
    actual_has_covar = ncm_fit_state_has_covar (test->fit_state);

    g_assert_cmpint (actual_has_covar, ==, expected_has_covar);
  }
}

void
test_ncm_fit_state_is_best_fit (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    gboolean expected_is_best_fit = g_test_rand_int_range (0, 1);
    gboolean actual_is_best_fit;

    ncm_fit_state_set_is_best_fit (test->fit_state, expected_is_best_fit);
    /* Check if the value was set correctly */
    actual_is_best_fit = ncm_fit_state_is_best_fit (test->fit_state);

    g_assert_cmpint (actual_is_best_fit, ==, expected_is_best_fit);
  }
}

void
test_ncm_fit_state_is_least_squares (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));

  {
    gboolean expected_is_least_squares = g_test_rand_int_range (0, 1);
    gboolean actual_is_least_squares;

    ncm_fit_state_set_is_least_squares (test->fit_state, expected_is_least_squares);
    /* Check if the value was set correctly */
    actual_is_least_squares = ncm_fit_state_is_least_squares (test->fit_state);

    g_assert_cmpint (actual_is_least_squares, ==, expected_is_least_squares);
  }
}

void
test_ncm_fit_state_fparams (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));
  {
    guint fparam_len = g_test_rand_int_range (1, 100);

    ncm_fit_state_set_fparam_len (test->fit_state, fparam_len);

    {
      NcmVector *fparam = ncm_fit_state_peek_fparams (test->fit_state);

      g_assert_true (fparam != NULL);
      g_assert_true (NCM_IS_VECTOR (fparam));
      g_assert_cmpuint (ncm_vector_len (fparam), ==, fparam_len);
      g_assert_cmpuint (ncm_vector_len (fparam), ==, ncm_fit_state_get_fparam_len (test->fit_state));
    }
  }
}

void
test_ncm_fit_state_hessian (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));
  {
    guint fparam_len = g_test_rand_int_range (1, 100);

    ncm_fit_state_set_fparam_len (test->fit_state, fparam_len);
    {
      NcmMatrix *hessian = ncm_fit_state_peek_hessian (test->fit_state);

      g_assert_true (hessian != NULL);
      g_assert_true (NCM_IS_MATRIX (hessian));
      g_assert_cmpuint (ncm_matrix_nrows (hessian), ==, fparam_len);
      g_assert_cmpuint (ncm_matrix_ncols (hessian), ==, fparam_len);
      g_assert_cmpuint (ncm_matrix_nrows (hessian), ==, ncm_fit_state_get_fparam_len (test->fit_state));
      g_assert_cmpuint (ncm_matrix_ncols (hessian), ==, ncm_fit_state_get_fparam_len (test->fit_state));
    }
  }
  ncm_fit_state_set_fparam_len (test->fit_state, 0);
  {
    NcmMatrix *hessian = ncm_fit_state_peek_hessian (test->fit_state);

    g_assert_true (hessian == NULL);
  }
}

void
test_ncm_fit_state_covar (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));
  {
    guint fparam_len = g_test_rand_int_range (1, 100);

    ncm_fit_state_set_fparam_len (test->fit_state, fparam_len);
    {
      NcmMatrix *covar = ncm_matrix_new (fparam_len, fparam_len);

      g_assert_true (covar != NULL);
      g_assert_true (NCM_IS_MATRIX (covar));
      g_assert_cmpuint (ncm_matrix_nrows (covar), ==, fparam_len);
      g_assert_cmpuint (ncm_matrix_ncols (covar), ==, fparam_len);
      g_assert_cmpuint (ncm_matrix_nrows (covar), ==, ncm_fit_state_get_fparam_len (test->fit_state));
      g_assert_cmpuint (ncm_matrix_ncols (covar), ==, ncm_fit_state_get_fparam_len (test->fit_state));
    }
  }
  ncm_fit_state_set_fparam_len (test->fit_state, 0);
  {
    NcmMatrix *covar = ncm_fit_state_peek_covar (test->fit_state);

    g_assert_true (covar == NULL);
  }
}

void
test_ncm_fit_state_f (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));
  {
    guint fparam_len = g_test_rand_int_range (1, 100);
    guint data_len   = g_test_rand_int_range (1, 100);

    ncm_fit_state_set_fparam_len (test->fit_state, fparam_len);
    ncm_fit_state_set_data_len (test->fit_state, data_len);

    ncm_fit_state_set_is_least_squares (test->fit_state, TRUE);
    {
      NcmVector *f = ncm_fit_state_peek_f (test->fit_state);

      g_assert_true (f != NULL);
      g_assert_true (NCM_IS_VECTOR (f));
      g_assert_cmpuint (ncm_vector_len (f), ==, data_len);
      g_assert_cmpuint (ncm_vector_len (f), ==, ncm_fit_state_get_data_len (test->fit_state));
    }

    fparam_len = g_test_rand_int_range (1, 100);
    ncm_fit_state_set_fparam_len (test->fit_state, fparam_len);
    {
      NcmVector *f = ncm_fit_state_peek_f (test->fit_state);

      g_assert_true (f != NULL);
      g_assert_true (NCM_IS_VECTOR (f));
      g_assert_cmpuint (ncm_vector_len (f), ==, data_len);
      g_assert_cmpuint (ncm_vector_len (f), ==, ncm_fit_state_get_data_len (test->fit_state));
    }

    data_len = g_test_rand_int_range (1, 100);
    ncm_fit_state_set_data_len (test->fit_state, data_len);
    {
      NcmVector *f = ncm_fit_state_peek_f (test->fit_state);

      g_assert_true (f != NULL);
      g_assert_true (NCM_IS_VECTOR (f));
      g_assert_cmpuint (ncm_vector_len (f), ==, data_len);
      g_assert_cmpuint (ncm_vector_len (f), ==, ncm_fit_state_get_data_len (test->fit_state));
    }

    ncm_fit_state_set_fparam_len (test->fit_state, 0);
    {
      NcmVector *f = ncm_fit_state_peek_f (test->fit_state);

      g_assert_true (f != NULL);
      g_assert_true (NCM_IS_VECTOR (f));
      g_assert_cmpuint (ncm_vector_len (f), ==, data_len);
      g_assert_cmpuint (ncm_vector_len (f), ==, ncm_fit_state_get_data_len (test->fit_state));
    }

    ncm_fit_state_set_data_len (test->fit_state, 0);
    {
      NcmVector *f = ncm_fit_state_peek_f (test->fit_state);

      g_assert_true (f == NULL);
    }
  }
}

void
test_ncm_fit_state_J (TestNcmFitState *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_FIT_STATE (test->fit_state));
  {
    guint fparam_len = g_test_rand_int_range (1, 100);
    guint data_len   = g_test_rand_int_range (1, 100);

    ncm_fit_state_set_fparam_len (test->fit_state, fparam_len);
    ncm_fit_state_set_data_len (test->fit_state, data_len);

    ncm_fit_state_set_is_least_squares (test->fit_state, TRUE);
    {
      NcmMatrix *J = ncm_fit_state_peek_J (test->fit_state);

      g_assert_true (J != NULL);
      g_assert_true (NCM_IS_MATRIX (J));
      g_assert_cmpuint (ncm_matrix_nrows (J), ==, data_len);
      g_assert_cmpuint (ncm_matrix_ncols (J), ==, fparam_len);
      g_assert_cmpuint (ncm_matrix_nrows (J), ==, ncm_fit_state_get_data_len (test->fit_state));
      g_assert_cmpuint (ncm_matrix_ncols (J), ==, ncm_fit_state_get_fparam_len (test->fit_state));
    }

    ncm_fit_state_set_fparam_len (test->fit_state, 0);
    {
      NcmMatrix *J = ncm_fit_state_peek_J (test->fit_state);

      g_assert_true (J == NULL);
    }

    ncm_fit_state_set_fparam_len (test->fit_state, fparam_len);
    ncm_fit_state_set_data_len (test->fit_state, 0);
    {
      NcmMatrix *J = ncm_fit_state_peek_J (test->fit_state);

      g_assert_true (J == NULL);
    }
  }
}

void
test_ncm_fit_state_set_ls_wrong (TestNcmFitState *test, gconstpointer pdata)
{
  ncm_fit_state_set_is_least_squares (test->fit_state, FALSE);

  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_fit_state_set_ls (test->fit_state, NULL, NULL);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fit_state_set_m2lnL_prec_wrong_lower_bound (TestNcmFitState *test, gconstpointer pdata)
{
  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_fit_state_set_m2lnL_prec (test->fit_state, -1.0);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fit_state_set_m2lnL_prec_wrong_upper_bound (TestNcmFitState *test, gconstpointer pdata)
{
  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_fit_state_set_m2lnL_prec (test->fit_state, 1.1);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fit_state_set_params_prec_wrong_lower_bound (TestNcmFitState *test, gconstpointer pdata)
{
  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_fit_state_set_params_prec (test->fit_state, -1.0);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fit_state_set_params_prec_wrong_upper_bound (TestNcmFitState *test, gconstpointer pdata)
{
  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_fit_state_set_params_prec (test->fit_state, 1.1);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fit_state_set_elapsed_time_wrong_lower_bound (TestNcmFitState *test, gconstpointer pdata)
{
  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_fit_state_set_elapsed_time (test->fit_state, -1.0);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fit_state_peek_f_wrong (TestNcmFitState *test, gconstpointer pdata)
{
  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_fit_state_set_is_least_squares (test->fit_state, FALSE);
    ncm_fit_state_peek_f (test->fit_state);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_fit_state_peek_J_wrong (TestNcmFitState *test, gconstpointer pdata)
{
  /* LCOV_EXCL_START */
  if (g_test_subprocess ())
  {
    ncm_fit_state_set_is_least_squares (test->fit_state, FALSE);
    ncm_fit_state_peek_J (test->fit_state);

    return;
  }

  /* LCOV_EXCL_STOP */

  /* Reruns this same test in a subprocess */
  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
}

