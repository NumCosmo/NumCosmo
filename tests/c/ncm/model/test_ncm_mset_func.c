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
 * A concrete NcmMSetFunc of the free parameters: f = sum_i (i + 1) p_i^2.
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
  guint i;

  res[0] = 0.0;

  for (i = 0; i < fparam_len; i++)
  {
    const gdouble p_i = ncm_mset_fparam_get (mset, i);

    res[0] += (i + 1.0) * p_i * p_i;
  }
}

static void
test_mset_func_fparams_class_init (TestMSetFuncFParamsClass *klass)
{
  NcmMSetFuncClass *func_class = NCM_MSET_FUNC_CLASS (klass);

  func_class->eval = &_test_mset_func_fparams_eval;
}

void test_ncm_mset_func_unames (void);
void test_ncm_mset_func_numdiff_fparams (void);
void test_ncm_mset_func_unames_set_meta (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/mset_func/unames", &test_ncm_mset_func_unames);
  g_test_add_func ("/ncm/mset_func/unames/set_meta", &test_ncm_mset_func_unames_set_meta);
  g_test_add_func ("/ncm/mset_func/numdiff_fparams", &test_ncm_mset_func_numdiff_fparams);

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

