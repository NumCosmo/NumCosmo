/***************************************************************************
 *            test_ncm_mset_func.c
 *
 *  Mon Sep 28 10:00:00 2026
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
 * A minimal concrete NcmMSetFunc: it returns the sum of its arguments.
 */

#define TEST_TYPE_MSET_FUNC_SUM (test_mset_func_sum_get_type ())
G_DECLARE_FINAL_TYPE (TestMSetFuncSum, test_mset_func_sum, TEST, MSET_FUNC_SUM, NcmMSetFunc)

struct _TestMSetFuncSum
{
  NcmMSetFunc parent_instance;
};

G_DEFINE_TYPE (TestMSetFuncSum, test_mset_func_sum, NCM_TYPE_MSET_FUNC)

static void
test_mset_func_sum_init (TestMSetFuncSum *fs)
{
}

static void
_test_mset_func_sum_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  const guint nvar = ncm_mset_func_get_nvar (func);
  guint i;

  res[0] = 0.0;

  for (i = 0; i < nvar; i++)
    res[0] += x[i];
}

static void
test_mset_func_sum_class_init (TestMSetFuncSumClass *klass)
{
  NcmMSetFuncClass *func_class = NCM_MSET_FUNC_CLASS (klass);

  func_class->eval = &_test_mset_func_sum_eval;
}

static NcmMSetFunc *
test_mset_func_sum_new (const guint nvar)
{
  NcmMSetFunc *func = g_object_new (TEST_TYPE_MSET_FUNC_SUM, NULL);

  ncm_mset_func_set_meta (func, "f", "f", "Test", "Sum of the arguments", nvar, 1);

  return func;
}

/*
 * A concrete NcmMSetFunc of the free parameters and its arguments:
 * f = prod_j x_j sum_i (i + 1) p_i^2.
 */

#define TEST_TYPE_MSET_FUNC_FPARAMS (test_mset_func_fparams_get_type ())
G_DECLARE_FINAL_TYPE (TestMSetFuncFParams, test_mset_func_fparams, TEST, MSET_FUNC_FPARAMS, NcmMSetFunc)

struct _TestMSetFuncFParams
{
  NcmMSetFunc parent_instance;
};

G_DEFINE_TYPE (TestMSetFuncFParams, test_mset_func_fparams, NCM_TYPE_MSET_FUNC)

static void
test_mset_func_fparams_init (TestMSetFuncFParams *fp)
{
}

static void
_test_mset_func_fparams_eval (NcmMSetFunc *func, NcmMSet *mset, const gdouble *x, gdouble *res)
{
  const guint fparam_len = ncm_mset_fparam_len (mset);
  const guint nvar       = ncm_mset_func_get_nvar (func);
  guint i;

  res[0] = 0.0;

  for (i = 0; i < fparam_len; i++)
  {
    const gdouble p_i = ncm_mset_fparam_get (mset, i);

    res[0] += (i + 1.0) * p_i * p_i;
  }

  for (i = 0; i < nvar; i++)
    res[0] *= x[i];
}

static void
test_mset_func_fparams_class_init (TestMSetFuncFParamsClass *klass)
{
  NcmMSetFuncClass *func_class = NCM_MSET_FUNC_CLASS (klass);

  func_class->eval = &_test_mset_func_fparams_eval;
}

/*
 * A concrete NcmMSetFunc1: its values are sum_j x_j + k for k = 0, ..., dim - 1.
 */

#define TEST_TYPE_MSET_FUNC1_SUM (test_mset_func1_sum_get_type ())
G_DECLARE_FINAL_TYPE (TestMSetFunc1Sum, test_mset_func1_sum, TEST, MSET_FUNC1_SUM, NcmMSetFunc1)

struct _TestMSetFunc1Sum
{
  NcmMSetFunc1 parent_instance;
  guint nret;
};

G_DEFINE_TYPE (TestMSetFunc1Sum, test_mset_func1_sum, NCM_TYPE_MSET_FUNC1)

static void
test_mset_func1_sum_init (TestMSetFunc1Sum *f1s)
{
  f1s->nret = 0;
}

static GArray *
_test_mset_func1_sum_eval1 (NcmMSetFunc1 *f1, NcmMSet *mset, GArray *x)
{
  TestMSetFunc1Sum *f1s = TEST_MSET_FUNC1_SUM (f1);
  GArray *res           = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), f1s->nret);
  gdouble sum           = 0.0;
  guint i;

  for (i = 0; i < x->len; i++)
    sum += g_array_index (x, gdouble, i);

  for (i = 0; i < f1s->nret; i++)
  {
    const gdouble v = sum + i;

    g_array_append_val (res, v);
  }

  return res;
}

static void
test_mset_func1_sum_class_init (TestMSetFunc1SumClass *klass)
{
  NcmMSetFunc1Class *func1_class = NCM_MSET_FUNC1_CLASS (klass);

  func1_class->eval1 = &_test_mset_func1_sum_eval1;
}

static NcmMSetFunc *
test_mset_func1_sum_new (const guint nvar, const guint dim, const guint nret)
{
  TestMSetFunc1Sum *f1s = g_object_new (TEST_TYPE_MSET_FUNC1_SUM, NULL);
  NcmMSetFunc *func     = NCM_MSET_FUNC (f1s);

  f1s->nret = nret;
  ncm_mset_func_set_meta (func, "f", "f", "Test", "Shifted sums of the arguments", nvar, dim);

  return func;
}

void test_ncm_mset_func_unames (void);
void test_ncm_mset_func_numdiff_fparams (void);
void test_ncm_mset_func_numdiff_fparams_eval_x (void);
void test_ncm_mset_func_eval_args (void);
void test_ncm_mset_func_eval_no_args (void);
void test_ncm_mset_func_eval_no_args_subprocess (void);
void test_ncm_mset_func_unames_set_meta (void);
void test_ncm_mset_func_eval_array (void);
void test_ncm_mset_func_eval_array_bad_len (void);
void test_ncm_mset_func_eval_array_bad_len_subprocess (void);
void test_ncm_mset_func1_eval (void);
void test_ncm_mset_func1_bad_dim (void);
void test_ncm_mset_func1_bad_dim_subprocess (void);
void test_ncm_mset_func_eval_x_prop (void);
void test_ncm_mset_func_eval_x_bad_len (void);
void test_ncm_mset_func_eval_x_bad_len_subprocess (void);
void test_ncm_mset_func_not_scalar (void);
void test_ncm_mset_func_not_scalar_eval0_subprocess (void);
void test_ncm_mset_func_not_scalar_eval1_subprocess (void);
void test_ncm_mset_func_not_scalar_eval_vector_subprocess (void);
void test_ncm_mset_func_eval_vector_bad_len_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/mset_func/unames", &test_ncm_mset_func_unames);
  g_test_add_func ("/ncm/mset_func/unames/set_meta", &test_ncm_mset_func_unames_set_meta);
  g_test_add_func ("/ncm/mset_func/numdiff_fparams", &test_ncm_mset_func_numdiff_fparams);
  g_test_add_func ("/ncm/mset_func/numdiff_fparams/eval_x", &test_ncm_mset_func_numdiff_fparams_eval_x);
  g_test_add_func ("/ncm/mset_func/eval/args", &test_ncm_mset_func_eval_args);
  g_test_add_func ("/ncm/mset_func/eval/no_args", &test_ncm_mset_func_eval_no_args);
  g_test_add_func ("/ncm/mset_func/eval/no_args/subprocess", &test_ncm_mset_func_eval_no_args_subprocess);
  g_test_add_func ("/ncm/mset_func/eval_array", &test_ncm_mset_func_eval_array);
  g_test_add_func ("/ncm/mset_func/eval_array/bad_len", &test_ncm_mset_func_eval_array_bad_len);
  g_test_add_func ("/ncm/mset_func/eval_array/bad_len/subprocess", &test_ncm_mset_func_eval_array_bad_len_subprocess);
  g_test_add_func ("/ncm/mset_func1/eval", &test_ncm_mset_func1_eval);
  g_test_add_func ("/ncm/mset_func1/bad_dim", &test_ncm_mset_func1_bad_dim);
  g_test_add_func ("/ncm/mset_func1/bad_dim/subprocess", &test_ncm_mset_func1_bad_dim_subprocess);
  g_test_add_func ("/ncm/mset_func/eval_x/prop", &test_ncm_mset_func_eval_x_prop);
  g_test_add_func ("/ncm/mset_func/eval_x/bad_len", &test_ncm_mset_func_eval_x_bad_len);
  g_test_add_func ("/ncm/mset_func/eval_x/bad_len/subprocess", &test_ncm_mset_func_eval_x_bad_len_subprocess);
  g_test_add_func ("/ncm/mset_func/not_scalar", &test_ncm_mset_func_not_scalar);
  g_test_add_func ("/ncm/mset_func/not_scalar/eval0/subprocess", &test_ncm_mset_func_not_scalar_eval0_subprocess);
  g_test_add_func ("/ncm/mset_func/not_scalar/eval1/subprocess", &test_ncm_mset_func_not_scalar_eval1_subprocess);
  g_test_add_func ("/ncm/mset_func/not_scalar/eval_vector/subprocess", &test_ncm_mset_func_not_scalar_eval_vector_subprocess);
  g_test_add_func ("/ncm/mset_func/eval_vector/bad_len/subprocess", &test_ncm_mset_func_eval_vector_bad_len_subprocess);

  g_test_run ();
}

typedef struct _TestUNameCase
{
  guint len;
  gdouble x[2];
  const gchar *uname;
  const gchar *usymbol;
} TestUNameCase;

void
test_ncm_mset_func_unames (void)
{
  const TestUNameCase cases[] = {
    {1, {1.0,         0.0},    "f_1",                   "f(1)"                  },
    {1, {-1.0,        0.0},    "f_m1",                  "f(-1)"                 },
    {1, {1.0e-5,      0.0},    "f_1em05",               "f(1e-05)"              },
    {1, {1.0e5,       0.0},    "f_100000",              "f(100000)"             },
    {1, {1.0e20,      0.0},    "f_1e20",                "f(1e+20)"              },
    {1, {0.3,         0.0},    "f_0p3",                 "f(0.3)"                },
    {1, {0.1 + 0.2,   0.0},    "f_0p30000000000000004", "f(0.30000000000000004)"},
    {1, {1.52,        0.0},    "f_1p52",                "f(1.52)"               },
    {2, {1.0,         2.0},    "f_1_2",                 "f(1,2)"                },
    {2, {12.0,        0.0},    "f_12_0",                "f(12,0)"               },
    {2, {1.5,         -2.0e-5}, "f_1p5_m2em05",         "f(1.5,-2e-05)"         },
  };
  guint i;

  for (i = 0; i < G_N_ELEMENTS (cases); i++)
  {
    NcmMSetFunc *func = test_mset_func_sum_new (cases[i].len);

    g_assert_cmpstr (ncm_mset_func_peek_uname (func), ==, "f");
    g_assert_cmpstr (ncm_mset_func_peek_usymbol (func), ==, "f");

    ncm_mset_func_set_eval_x (func, (gdouble *) cases[i].x, cases[i].len);

    g_assert_cmpstr (ncm_mset_func_peek_uname (func), ==, cases[i].uname);
    g_assert_cmpstr (ncm_mset_func_peek_usymbol (func), ==, cases[i].usymbol);

    ncm_mset_func_free (func);
  }

  /* The shortest form reads back to the same double. */
  {
    const gdouble x   = 0.1 + 0.2;
    NcmMSetFunc *func = test_mset_func_sum_new (1);
    const gchar *sym;

    ncm_mset_func_set_eval_x (func, (gdouble *) &x, 1);
    sym = ncm_mset_func_peek_usymbol (func);

    g_assert_cmpfloat (g_ascii_strtod (sym + 2, NULL), ==, x);

    ncm_mset_func_free (func);
  }
}

void
test_ncm_mset_func_unames_set_meta (void)
{
  NcmMSetFunc *func = test_mset_func_sum_new (1);
  gdouble x         = 2.5;

  ncm_mset_func_set_eval_x (func, &x, 1);
  g_assert_cmpstr (ncm_mset_func_peek_uname (func), ==, "f_2p5");

  ncm_mset_func_set_meta (func, "g", "g", "Test", "Sum of the arguments", 1, 1);

  g_assert_cmpstr (ncm_mset_func_peek_uname (func), ==, "g_2p5");
  g_assert_cmpstr (ncm_mset_func_peek_usymbol (func), ==, "g(2.5)");

  ncm_mset_func_free (func);
}

void
test_ncm_mset_func_numdiff_fparams (void)
{
  NcmModelRosenbrock *mrb = ncm_model_rosenbrock_new ();
  NcmMSet *mset           = ncm_mset_new (mrb, NULL, NULL);
  NcmMSetFunc *func       = g_object_new (TEST_TYPE_MSET_FUNC_FPARAMS, NULL);
  const gdouble p[2]      = {1.5, -0.5};
  NcmVector *grad         = NULL;
  NcmVector *grad_data    = NULL;
  NcmVector *grad_in;
  guint i;

  ncm_mset_func_set_meta (func, "f", "f", "Test", "Weighted sum of squares", 0, 1);

  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X1, p[0]);
  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X2, p[1]);
  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);
  g_assert_cmpuint (ncm_mset_fparam_len (mset), ==, 2);

  /* A NULL *out allocates a new vector. */
  ncm_mset_func_numdiff_fparams (func, mset, NULL, &grad);
  g_assert_nonnull (grad);
  g_assert_cmpuint (ncm_vector_len (grad), ==, 2);

  for (i = 0; i < 2; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (grad, i), ==, 2.0 * (i + 1.0) * p[i], 1.0e-12, 0.0);

  /* The free parameters are restored. */
  for (i = 0; i < 2; i++)
    g_assert_cmpfloat (ncm_mset_fparam_get (mset, i), ==, p[i]);

  /* A non-NULL *out is overwritten in place. */
  grad_data = ncm_vector_new (2);
  grad_in   = grad_data;
  ncm_vector_set_all (grad_data, -7.0);

  ncm_mset_func_numdiff_fparams (func, mset, NULL, &grad_data);
  g_assert_true (grad_data == grad_in);

  for (i = 0; i < 2; i++)
    g_assert_cmpfloat (ncm_vector_get (grad_data, i), ==, ncm_vector_get (grad, i));

  ncm_vector_free (grad);
  ncm_vector_free (grad_data);
  ncm_mset_func_free (func);
  ncm_mset_free (mset);
  ncm_model_rosenbrock_free (mrb);
}

void
test_ncm_mset_func_numdiff_fparams_eval_x (void)
{
  NcmModelRosenbrock *mrb = ncm_model_rosenbrock_new ();
  NcmMSet *mset           = ncm_mset_new (mrb, NULL, NULL);
  NcmMSetFunc *func       = g_object_new (TEST_TYPE_MSET_FUNC_FPARAMS, NULL);
  const gdouble p[2]      = {1.5, -0.5};
  gdouble eval_x          = 3.0;
  gdouble x               = 5.0;
  NcmVector *grad         = NULL;
  guint i;

  ncm_mset_func_set_meta (func, "f", "f", "Test", "Scaled weighted sum of squares", 1, 1);
  ncm_mset_func_set_eval_x (func, &eval_x, 1);

  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X1, p[0]);
  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X2, p[1]);
  ncm_mset_param_set_all_ftype (mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  /* A NULL x differentiates at the evaluation point. */
  ncm_mset_func_numdiff_fparams (func, mset, NULL, &grad);

  for (i = 0; i < 2; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (grad, i), ==, eval_x * 2.0 * (i + 1.0) * p[i], 1.0e-12, 0.0);

  /* An explicit x wins over the evaluation point. */
  ncm_mset_func_numdiff_fparams (func, mset, &x, &grad);

  for (i = 0; i < 2; i++)
    ncm_assert_cmpdouble_e (ncm_vector_get (grad, i), ==, x * 2.0 * (i + 1.0) * p[i], 1.0e-12, 0.0);

  ncm_vector_free (grad);
  ncm_mset_func_free (func);
  ncm_mset_free (mset);
  ncm_model_rosenbrock_free (mrb);
}

void
test_ncm_mset_func_eval_args (void)
{
  NcmMSet *mset       = ncm_mset_empty_new ();
  NcmMSetFunc *func2  = test_mset_func_sum_new (2);
  NcmMSetFunc *func1  = test_mset_func_sum_new (1);
  NcmMSetFunc *func0  = test_mset_func_sum_new (0);
  gdouble eval_x2[2]  = {1.0, 2.0};
  gdouble x2[2]       = {5.0, 6.0};
  gdouble eval_x1     = 2.0;
  gdouble x_v_data[3] = {1.0, 4.0, 8.0};
  NcmVector *x_v      = ncm_vector_new_data_static (x_v_data, 3, 1);
  NcmVector *res_v    = ncm_vector_new (3);
  gdouble res;
  guint i;

  /* Without an evaluation point an explicit x is used. */
  ncm_mset_func_eval (func2, mset, x2, &res);
  g_assert_cmpfloat (res, ==, 11.0);
  g_assert_cmpfloat (ncm_mset_func_eval_nvar (func2, mset, x2), ==, 11.0);

  /* With one, a NULL x uses it and an explicit x wins. */
  ncm_mset_func_set_eval_x (func2, eval_x2, 2);

  ncm_mset_func_eval (func2, mset, NULL, &res);
  g_assert_cmpfloat (res, ==, 3.0);
  g_assert_cmpfloat (ncm_mset_func_eval_nvar (func2, mset, NULL), ==, 3.0);
  g_assert_cmpfloat (ncm_mset_func_eval0 (func2, mset), ==, 3.0);

  ncm_mset_func_eval (func2, mset, x2, &res);
  g_assert_cmpfloat (res, ==, 11.0);
  g_assert_cmpfloat (ncm_mset_func_eval_nvar (func2, mset, x2), ==, 11.0);

  /* eval1 and eval_vector always use their arguments. */
  ncm_mset_func_set_eval_x (func1, &eval_x1, 1);
  g_assert_cmpfloat (ncm_mset_func_eval0 (func1, mset), ==, 2.0);
  g_assert_cmpfloat (ncm_mset_func_eval1 (func1, mset, 7.0), ==, 7.0);

  ncm_mset_func_eval_vector (func1, mset, x_v, res_v);

  for (i = 0; i < 3; i++)
    g_assert_cmpfloat (ncm_vector_get (res_v, i), ==, x_v_data[i]);

  /* A function without variables takes no arguments. */
  g_assert_cmpfloat (ncm_mset_func_eval0 (func0, mset), ==, 0.0);
  g_assert_cmpfloat (ncm_mset_func_eval_nvar (func0, mset, NULL), ==, 0.0);

  ncm_vector_free (x_v);
  ncm_vector_free (res_v);
  ncm_mset_func_free (func0);
  ncm_mset_func_free (func1);
  ncm_mset_func_free (func2);
  ncm_mset_free (mset);
}

void
test_ncm_mset_func_eval_no_args (void)
{
  g_test_trap_subprocess ("/ncm/mset_func/eval/no_args/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*function `f' takes 1 variable(s), but it was called without arguments and no evaluation point is set*");
}

void
test_ncm_mset_func_eval_no_args_subprocess (void)
{
  NcmMSet *mset     = ncm_mset_empty_new ();
  NcmMSetFunc *func = test_mset_func_sum_new (1);

  ncm_mset_func_eval0 (func, mset);

  ncm_mset_func_free (func);
  ncm_mset_free (mset);
}

void
test_ncm_mset_func_eval_array (void)
{
  NcmMSet *mset      = ncm_mset_empty_new ();
  NcmMSetFunc *func2 = test_mset_func_sum_new (2);
  NcmMSetFunc *func0 = test_mset_func_sum_new (0);
  gdouble eval_x2[2] = {1.0, 2.0};
  gdouble x2_data[2] = {5.0, 6.0};
  GArray *x2         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *res;

  g_array_append_vals (x2, x2_data, 2);

  res = ncm_mset_func_eval_array (func2, mset, x2);
  g_assert_cmpuint (res->len, ==, 1);
  g_assert_cmpfloat (g_array_index (res, gdouble, 0), ==, 11.0);
  g_array_unref (res);

  ncm_mset_func_set_eval_x (func2, eval_x2, 2);

  res = ncm_mset_func_eval_array (func2, mset, NULL);
  g_assert_cmpfloat (g_array_index (res, gdouble, 0), ==, 3.0);
  g_array_unref (res);

  res = ncm_mset_func_eval_array (func2, mset, x2);
  g_assert_cmpfloat (g_array_index (res, gdouble, 0), ==, 11.0);
  g_array_unref (res);

  res = ncm_mset_func_eval_array (func0, mset, NULL);
  g_assert_cmpuint (res->len, ==, 1);
  g_assert_cmpfloat (g_array_index (res, gdouble, 0), ==, 0.0);
  g_array_unref (res);

  g_array_unref (x2);
  ncm_mset_func_free (func0);
  ncm_mset_func_free (func2);
  ncm_mset_free (mset);
}

void
test_ncm_mset_func_eval_array_bad_len (void)
{
  g_test_trap_subprocess ("/ncm/mset_func/eval_array/bad_len/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*function `f' takes 2 variable(s), but 1 argument(s) were given*");
}

void
test_ncm_mset_func_eval_array_bad_len_subprocess (void)
{
  NcmMSet *mset     = ncm_mset_empty_new ();
  NcmMSetFunc *func = test_mset_func_sum_new (2);
  GArray *x         = g_array_new (FALSE, FALSE, sizeof (gdouble));
  const gdouble x0  = 1.0;

  g_array_append_val (x, x0);
  g_array_unref (ncm_mset_func_eval_array (func, mset, x));

  g_array_unref (x);
  ncm_mset_func_free (func);
  ncm_mset_free (mset);
}

void
test_ncm_mset_func1_eval (void)
{
  NcmMSet *mset      = ncm_mset_empty_new ();
  NcmMSetFunc *func  = test_mset_func1_sum_new (2, 3, 3);
  NcmMSetFunc *func1 = test_mset_func1_sum_new (1, 1, 1);
  gdouble x[2]       = {1.0, 2.0};
  gdouble eval_x[2]  = {10.0, 20.0};
  gdouble res[3];
  GArray *res_a;
  guint i;

  g_assert_true (NCM_IS_MSET_FUNC1 (func));

  ncm_mset_func_eval (func, mset, x, res);

  for (i = 0; i < 3; i++)
    g_assert_cmpfloat (res[i], ==, 3.0 + i);

  ncm_mset_func_set_eval_x (func, eval_x, 2);
  res_a = ncm_mset_func_eval_array (func, mset, NULL);
  g_assert_cmpuint (res_a->len, ==, 3);

  for (i = 0; i < 3; i++)
    g_assert_cmpfloat (g_array_index (res_a, gdouble, i), ==, 30.0 + i);

  g_array_unref (res_a);

  g_assert_cmpfloat (ncm_mset_func_eval1 (func1, mset, 6.0), ==, 6.0);

  ncm_mset_func_free (func1);
  ncm_mset_func_free (func);
  ncm_mset_free (mset);
}

void
test_ncm_mset_func1_bad_dim (void)
{
  g_test_trap_subprocess ("/ncm/mset_func1/bad_dim/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*function `f' has dimension 2, but eval1 returned 1 value(s)*");
}

void
test_ncm_mset_func1_bad_dim_subprocess (void)
{
  NcmMSet *mset     = ncm_mset_empty_new ();
  NcmMSetFunc *func = test_mset_func1_sum_new (0, 2, 1);

  g_array_unref (ncm_mset_func_eval_array (func, mset, NULL));

  ncm_mset_func_free (func);
  ncm_mset_free (mset);
}

void
test_ncm_mset_func_eval_x_prop (void)
{
  NcmMSet *mset     = ncm_mset_empty_new ();
  NcmMSetFunc *func = test_mset_func_sum_new (2);
  NcmMatrix *m      = ncm_matrix_new (2, 2);
  NcmVector *col;
  NcmVector *eval_x;

  /* A strided column: its components are 1 and 2, the row holds 1 and 9. */
  ncm_matrix_set (m, 0, 0, 1.0);
  ncm_matrix_set (m, 1, 0, 2.0);
  ncm_matrix_set (m, 0, 1, 9.0);
  ncm_matrix_set (m, 1, 1, 9.0);
  col = ncm_matrix_get_col (m, 0);
  g_assert_cmpuint (ncm_vector_stride (col), ==, 2);

  g_object_set (func, "eval-x", col, NULL);

  g_assert_true (ncm_mset_func_is_const (func));
  g_assert_cmpfloat (ncm_mset_func_eval0 (func, mset), ==, 3.0);
  g_assert_cmpstr (ncm_mset_func_peek_uname (func), ==, "f_1_2");

  /* The function holds a copy. */
  ncm_vector_set (col, 0, 5.0);
  g_assert_cmpfloat (ncm_mset_func_eval0 (func, mset), ==, 3.0);
  g_assert_cmpstr (ncm_mset_func_peek_uname (func), ==, "f_1_2");

  g_object_get (func, "eval-x", &eval_x, NULL);
  g_assert_true (eval_x != col);
  g_assert_cmpuint (ncm_vector_len (eval_x), ==, 2);
  g_assert_cmpfloat (ncm_vector_get (eval_x, 0), ==, 1.0);
  g_assert_cmpfloat (ncm_vector_get (eval_x, 1), ==, 2.0);
  ncm_vector_free (eval_x);

  /* NULL clears the evaluation point. */
  g_object_set (func, "eval-x", NULL, NULL);
  g_assert_false (ncm_mset_func_is_const (func));
  g_assert_cmpstr (ncm_mset_func_peek_uname (func), ==, "f");

  ncm_vector_free (col);
  ncm_matrix_free (m);
  ncm_mset_func_free (func);
  ncm_mset_free (mset);
}

void
test_ncm_mset_func_eval_x_bad_len (void)
{
  g_test_trap_subprocess ("/ncm/mset_func/eval_x/bad_len/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*function `f' takes 2 variable(s), but the evaluation point has 3*");
}

void
test_ncm_mset_func_eval_x_bad_len_subprocess (void)
{
  NcmMSetFunc *func = test_mset_func_sum_new (2);
  NcmVector *x      = ncm_vector_new (3);

  ncm_vector_set_all (x, 1.0);
  g_object_set (func, "eval-x", x, NULL);

  ncm_vector_free (x);
  ncm_mset_func_free (func);
}

void
test_ncm_mset_func_not_scalar (void)
{
  g_test_trap_subprocess ("/ncm/mset_func/not_scalar/eval0/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*function `f' has dimension 2, but only scalar functions return a single value*");

  g_test_trap_subprocess ("/ncm/mset_func/not_scalar/eval1/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*function `f' takes 1 variable(s) and has dimension 2*");

  g_test_trap_subprocess ("/ncm/mset_func/not_scalar/eval_vector/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*function `f' takes 1 variable(s) and has dimension 2*");

  g_test_trap_subprocess ("/ncm/mset_func/eval_vector/bad_len/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*3 argument(s) but room for 2 value(s)*");
}

void
test_ncm_mset_func_not_scalar_eval0_subprocess (void)
{
  NcmMSet *mset     = ncm_mset_empty_new ();
  NcmMSetFunc *func = test_mset_func1_sum_new (0, 2, 2);

  ncm_mset_func_eval0 (func, mset);

  ncm_mset_func_free (func);
  ncm_mset_free (mset);
}

void
test_ncm_mset_func_not_scalar_eval1_subprocess (void)
{
  NcmMSet *mset     = ncm_mset_empty_new ();
  NcmMSetFunc *func = test_mset_func1_sum_new (1, 2, 2);

  ncm_mset_func_eval1 (func, mset, 1.0);

  ncm_mset_func_free (func);
  ncm_mset_free (mset);
}

void
test_ncm_mset_func_not_scalar_eval_vector_subprocess (void)
{
  NcmMSet *mset     = ncm_mset_empty_new ();
  NcmMSetFunc *func = test_mset_func1_sum_new (1, 2, 2);
  NcmVector *x_v    = ncm_vector_new (2);
  NcmVector *res_v  = ncm_vector_new (2);

  ncm_vector_set_all (x_v, 1.0);
  ncm_mset_func_eval_vector (func, mset, x_v, res_v);

  ncm_vector_free (x_v);
  ncm_vector_free (res_v);
  ncm_mset_func_free (func);
  ncm_mset_free (mset);
}

void
test_ncm_mset_func_eval_vector_bad_len_subprocess (void)
{
  NcmMSet *mset     = ncm_mset_empty_new ();
  NcmMSetFunc *func = test_mset_func_sum_new (1);
  NcmVector *x_v    = ncm_vector_new (3);
  NcmVector *res_v  = ncm_vector_new (2);

  ncm_vector_set_all (x_v, 1.0);
  ncm_mset_func_eval_vector (func, mset, x_v, res_v);

  ncm_vector_free (x_v);
  ncm_vector_free (res_v);
  ncm_mset_func_free (func);
  ncm_mset_free (mset);
}

