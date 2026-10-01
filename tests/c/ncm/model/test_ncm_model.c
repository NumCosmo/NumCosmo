/***************************************************************************
 *            test_ncm_model.c
 *
 *  Mon May 7 16:39:21 2012
 *  Copyright  2012  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima & Sandro Dias Pinto Vitenti 2012 <pennalima@gmail.com>, <vitenti@uel.br>
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
#include "ncm_model_test.h"

typedef struct _TestNcmModel
{
  GType type;
  guint sparam_len, vparam_len;
  gchar *name, *nick;
  NcmModelTest *tm;
  NcmReparam *reparam;
} TestNcmModel;

void test_ncm_model_new (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_new_modify (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_child_new (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_child_child_new (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_reparam_new (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_free (TestNcmModel *test, gconstpointer pdata);

void test_ncm_model_test_new (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_test_length (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_test_defval (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_test_name_symbol (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_test_setget (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_test_setget_prop (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_test_setget_vector (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_test_setget_model (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_test_name_index (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_test_dup (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_test_impl (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_param_names (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_test_finite (TestNcmModel *test, gconstpointer pdata);
void test_ncm_model_svparams_len (TestNcmModel *test, gconstpointer pdata);

void test_ncm_model_reparam_current (void);
void test_ncm_model_reparam_is_equal (void);
void test_ncm_model_reparam_remove (void);
void test_ncm_model_reparam_remove_subprocess (void);
void test_ncm_model_set_default_desc (void);
void test_ncm_model_id_by_type_error (void);
void test_ncm_model_param_desc (void);
void test_ncm_model_renamed_submodel_param (void);
void test_ncm_model_vparam_set_vector_len (void);
void test_ncm_model_vparam_set_vector_len_subprocess (void);

#define TEST_NCM_MODEL_NTYPES 5

gint
main (gint argc, gchar *argv[])
{
  gint i;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  gpointer ccc[TEST_NCM_MODEL_NTYPES][3] = {
    {"model",            &test_ncm_model_new,             &test_ncm_model_free},
    {"model/modified",   &test_ncm_model_new_modify,      &test_ncm_model_free},
    {"model/child",      &test_ncm_model_child_new,       &test_ncm_model_free},
    {"model/grandchild", &test_ncm_model_child_child_new, &test_ncm_model_free},
    {"model/reparam",    &test_ncm_model_reparam_new,     &test_ncm_model_free},
  };

  for (i = 0; i < TEST_NCM_MODEL_NTYPES; i++)
  {
    gchar *d;

    d = g_strdup_printf ("/ncm/%s/new", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_new, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/length", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_length, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/defval", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_defval, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/name_symbol", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_name_symbol, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/setget", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_setget, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/setget/prop", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_setget_prop, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/setget/vector", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_setget_vector, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/setget/model", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_setget_model, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/name_index", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_name_index, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/dup", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_dup, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/impl", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_impl, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/param_names", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_param_names, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/finite", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_test_finite, ccc[i][2]);
    g_free (d);

    d = g_strdup_printf ("/ncm/%s/svparams_len", (gchar *) ccc[i][0]);
    g_test_add (d, TestNcmModel, NULL, ccc[i][1], &test_ncm_model_svparams_len, ccc[i][2]);
    g_free (d);
  }

  g_test_add_func ("/ncm/model/reparam/current", &test_ncm_model_reparam_current);
  g_test_add_func ("/ncm/model/reparam/is_equal", &test_ncm_model_reparam_is_equal);
  g_test_add_func ("/ncm/model/reparam/remove", &test_ncm_model_reparam_remove);
  g_test_add_func ("/ncm/model/reparam/remove/subprocess", &test_ncm_model_reparam_remove_subprocess);
  g_test_add_func ("/ncm/model/set_default_desc", &test_ncm_model_set_default_desc);
  g_test_add_func ("/ncm/model/id_by_type_error", &test_ncm_model_id_by_type_error);
  g_test_add_func ("/ncm/model/param_desc", &test_ncm_model_param_desc);
  g_test_add_func ("/ncm/model/renamed_submodel_param", &test_ncm_model_renamed_submodel_param);
  g_test_add_func ("/ncm/model/vparam_set_vector_len", &test_ncm_model_vparam_set_vector_len);
  g_test_add_func ("/ncm/model/vparam_set_vector_len/subprocess", &test_ncm_model_vparam_set_vector_len_subprocess);

  g_test_run ();
}

void
test_ncm_model_new (TestNcmModel *test, gconstpointer pdata)
{
  test->type       = NCM_TYPE_MODEL_TEST;
  test->tm         = g_object_new (test->type, NULL);
  test->sparam_len = SPARAM_LEN1;
  test->vparam_len = VPARAM_LEN1;
  test->name       = name_tot[0];
  test->nick       = nick_tot[0];
  test->reparam    = NULL;

  g_assert_true (test->type != 0);
}

void
test_ncm_model_new_modify (TestNcmModel *test, gconstpointer pdata)
{
  test->type       = NCM_TYPE_MODEL_TEST;
  test->tm         = g_object_new (test->type, NULL);
  test->sparam_len = SPARAM_LEN1;
  test->vparam_len = VPARAM_LEN1;
  test->name       = name_tot[0];
  test->nick       = nick_tot[0];
  test->reparam    = NULL;

  {
    NcmModel *model = NCM_MODEL (test->tm);
    guint model_len = ncm_model_len (model);
    guint i;

    for (i = 0; i < model_len; i++)
    {
      const gdouble curval   = ncm_model_param_get (model, i);
      const gdouble s_lb     = g_test_rand_double_range (curval - fabs (curval) * 1.5, curval + fabs (curval) * 1.5);
      const gdouble s_ub     = g_test_rand_double_range (s_lb, curval + fabs (curval) * 1.5);
      const gdouble s_scale  = fabs (((s_ub + s_lb) != 0 ? (s_ub + s_lb) : 1.0) * pow (10.0, -g_test_rand_double_range (1.0,  2.0)));
      const gdouble s_abstol = s_scale * pow (10.0, -g_test_rand_double_range (1.0, 11.0));

      ncm_model_param_set_upper_bound (model, i, GSL_POSINF);
      ncm_model_param_set_lower_bound (model, i, s_lb);
      ncm_model_param_set_upper_bound (model, i, s_ub);
      ncm_model_param_set_scale (model, i, s_scale);
      ncm_model_param_set_abstol (model, i, s_abstol);
    }
  }

  g_assert_true (test->type != 0);
}

void
test_ncm_model_child_new (TestNcmModel *test, gconstpointer pdata)
{
  test->type       = NCM_TYPE_MODEL_TEST_CHILD;
  test->tm         = g_object_new (test->type, NULL);
  test->sparam_len = SPARAM_LEN1 + SPARAM_LEN2;
  test->vparam_len = VPARAM_LEN1 + VPARAM_LEN2;
  test->name       = name_tot[1];
  test->nick       = nick_tot[1];
  test->reparam    = NULL;

  g_assert_true (test->type != 0);
}

void
test_ncm_model_child_child_new (TestNcmModel *test, gconstpointer pdata)
{
  test->type       = NCM_TYPE_MODEL_TEST_CHILD_CHILD;
  test->tm         = g_object_new (test->type, NULL);
  test->sparam_len = SPARAM_LEN1 + SPARAM_LEN2 + SPARAM_LEN3;
  test->vparam_len = VPARAM_LEN1 + VPARAM_LEN2 + VPARAM_LEN3;
  test->name       = name_tot[2];
  test->nick       = nick_tot[2];
  test->reparam    = NULL;

  g_assert_true (test->type != 0);
}

NcmReparam *
_test_ncm_model_create_reparam (TestNcmModel *test)
{
  const guint size     = ncm_model_len (NCM_MODEL (test->tm));
  NcmMatrix *T         = ncm_matrix_new (size, size);
  NcmVector *v         = ncm_vector_new (size);
  NcmBootstrap *bstrap = ncm_bootstrap_sized_new (size);
  NcmRNG *rng          = ncm_rng_seeded_new (NULL, g_test_rand_int ());
  guint cdesc_n        = g_test_rand_int_range (size / 2, size);
  NcmReparamLinear *relin;
  guint i;

  if (cdesc_n == 0)
    cdesc_n = 1;

  ncm_bootstrap_remix (bstrap, rng);

  ncm_matrix_set_zero (T);
  ncm_vector_set_zero (v);

  for (i = 0; i < size; i++)
  {
    ncm_matrix_set (T, i, ncm_bootstrap_get (bstrap, i), g_test_rand_double_range (1.0, 10.0));
    ncm_vector_set (v, i, g_test_rand_double ());
  }

  relin = ncm_reparam_linear_new (size, T, v);
  ncm_reparam_linear_set_compat_type (relin, NCM_TYPE_MODEL_TEST);

  for (i = 0; i < cdesc_n; i++)
  {
    gchar *new_param        = g_strdup_printf ("new_param_%u", i);
    gchar *new_param_symbol = g_strdup_printf ("NP_%u", i);

    ncm_reparam_set_param_desc_full (NCM_REPARAM (relin),
                                     i,
                                     new_param,
                                     new_param_symbol,
                                     -10.0,
                                     10.0,
                                     1.0,
                                     0.0,
                                     1.0,
                                     NCM_PARAM_TYPE_FIXED);
    g_free (new_param);
    g_free (new_param_symbol);
  }

  ncm_model_set_reparam (NCM_MODEL (test->tm), NCM_REPARAM (relin), NULL);
  ncm_reparam_free (NCM_REPARAM (relin));

  ncm_vector_free (v);
  ncm_matrix_free (T);
  ncm_rng_free (rng);
  ncm_bootstrap_free (bstrap);

  return NCM_REPARAM (relin);
}

void
test_ncm_model_reparam_new (TestNcmModel *test, gconstpointer pdata)
{
  test->type       = NCM_TYPE_MODEL_TEST;
  test->tm         = g_object_new (test->type, NULL);
  test->sparam_len = SPARAM_LEN1;
  test->vparam_len = VPARAM_LEN1;
  test->name       = name_tot[0];
  test->nick       = nick_tot[0];
  test->reparam    = _test_ncm_model_create_reparam (test);

  g_assert_true (test->type != 0);
}

void
test_ncm_model_free (TestNcmModel *test, gconstpointer pdata)
{
  NCM_TEST_FREE (g_object_unref, test->tm);
}

void
test_ncm_model_test_new (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model;

  g_assert_true (NCM_IS_MODEL (tm));
  g_assert_true (NCM_IS_MODEL_TEST (tm));
  model = NCM_MODEL (tm);

  g_assert_true (ncm_model_name (model) != test->name);
  g_assert_true (ncm_model_nick (model) != test->nick);

  g_assert_cmpstr (ncm_model_name (model), ==, test->name);
  g_assert_cmpstr (ncm_model_nick (model), ==, test->nick);

  g_assert_true (ncm_model_peek_reparam (model) == test->reparam);
}

void
test_ncm_model_test_length (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model  = NCM_MODEL (tm);
  guint model_len  = 0;
  guint i;

  model_len = test->sparam_len;

  for (i = 0; i < test->vparam_len; i++)
  {
    g_assert_cmpint (ncm_model_vparam_len (model, i), ==, v_len_tot[i]);
    model_len += v_len_tot[i];
  }

  g_assert_cmpint (ncm_model_len (model), ==, model_len);
  g_assert_cmpint (ncm_model_sparam_len (model), ==, test->sparam_len);
  g_assert_cmpint (ncm_model_vparam_array_len (model), ==, test->vparam_len);
}

void
test_ncm_model_test_defval (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model  = NCM_MODEL (tm);
  guint i;

  for (i = 0; i < test->sparam_len; i++)
  {
    gdouble p = ncm_model_orig_param_get (model, i);

    ncm_assert_cmpdouble (p, ==, s_defval_tot[i]);
  }

  for (i = 0; i < test->vparam_len; i++)
  {
    guint j;

    for (j = 0; j < v_len_tot[i]; j++)
    {
      gdouble p = ncm_model_orig_vparam_get (model, i, j);

      ncm_assert_cmpdouble (p, ==, v_defval_tot[i]);
    }
  }
}

void
test_ncm_model_test_name_symbol (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model  = NCM_MODEL (tm);
  guint i;

  for (i = 0; i < test->sparam_len; i++)
  {
    const gchar *pname   = ncm_model_orig_param_name (model, i);
    const gchar *psymbol = ncm_model_orig_param_symbol (model, i);

    g_assert_true (pname != s_name_tot[i]);
    g_assert_true (psymbol != s_symbol_tot[i]);
    g_assert_cmpstr (pname, ==, s_name_tot[i]);
    g_assert_cmpstr (psymbol, ==, s_symbol_tot[i]);

    if (ncm_model_peek_reparam (model) != NULL)
    {
      NcmReparam *reparam       = ncm_model_peek_reparam (model);
      NcmSParam *reparam_sparam = ncm_reparam_get_param_desc (reparam, i);

      if (reparam_sparam != NULL)
      {
        const gchar *pname_new   = ncm_model_param_name (model, i);
        const gchar *psymbol_new = ncm_model_param_symbol (model, i);

        g_assert_cmpstr (pname_new, !=, s_name_tot[i]);
        g_assert_cmpstr (psymbol_new, !=, s_symbol_tot[i]);

        ncm_sparam_free (reparam_sparam);
      }
    }
  }

  for (i = 0; i < test->vparam_len; i++)
  {
    guint j;

    for (j = 0; j < v_len_tot[i]; j++)
    {
      gchar *vp_name       = g_strdup_printf ("%s_%u", v_name_tot[i], j);
      gchar *vp_symbol     = g_strdup_printf ("{%s}_%u", v_symbol_tot[i], j);
      const gchar *pname   = ncm_model_orig_param_name (model, ncm_model_vparam_index (model, i, j));
      const gchar *psymbol = ncm_model_orig_param_symbol (model, ncm_model_vparam_index (model, i, j));

      g_assert_true (pname != vp_name);
      g_assert_true (psymbol != vp_symbol);
      g_assert_cmpstr (pname, ==, vp_name);
      g_assert_cmpstr (psymbol, ==, vp_symbol);

      if (ncm_model_peek_reparam (model) != NULL)
      {
        guint pindex              = ncm_model_vparam_index (model, i, j);
        NcmReparam *reparam       = ncm_model_peek_reparam (model);
        NcmSParam *reparam_sparam = ncm_reparam_get_param_desc (reparam, pindex);

        if (reparam_sparam != NULL)
        {
          const gchar *pname_new   = ncm_model_param_name (model, pindex);
          const gchar *psymbol_new = ncm_model_param_symbol (model, pindex);

          g_assert_cmpstr (pname_new, !=, vp_name);
          g_assert_cmpstr (psymbol_new, !=, vp_symbol);

          ncm_sparam_free (reparam_sparam);
        }
      }

      g_free (vp_name);
      g_free (vp_symbol);
    }
  }
}

void
test_ncm_model_test_finite (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model  = NCM_MODEL (tm);
  guint model_len  = ncm_model_len (model);
  guint i;

  for (i = 0; i < model_len; i++)
  {
    g_assert_true (ncm_model_param_finite (model, i));
    g_assert_true (ncm_model_params_finite (model));

    ncm_model_param_set (model, i, GSL_NAN);
    g_assert_true (!ncm_model_param_finite (model, i));
    g_assert_true (!ncm_model_params_finite (model));
    ncm_model_param_set_default (model, i);
    g_assert_true (ncm_model_param_finite (model, i));
    g_assert_true (ncm_model_params_finite (model));

    ncm_model_param_set (model, i, GSL_POSINF);
    g_assert_true (!ncm_model_param_finite (model, i));
    g_assert_true (!ncm_model_params_finite (model));
    ncm_model_param_set_default (model, i);
    g_assert_true (ncm_model_param_finite (model, i));
    g_assert_true (ncm_model_params_finite (model));

    ncm_model_param_set (model, i, GSL_NEGINF);
    g_assert_true (!ncm_model_param_finite (model, i));
    g_assert_true (!ncm_model_params_finite (model));
    ncm_model_param_set_default (model, i);
    g_assert_true (ncm_model_param_finite (model, i));
    g_assert_true (ncm_model_params_finite (model));
  }
}

void
test_ncm_model_test_setget (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model  = NCM_MODEL (tm);
  guint model_len  = ncm_model_len (model);
  NcmVector *tmp   = ncm_vector_new (model_len);
  guint i;

  /* ncm_model_param_get / ncm_model_param_set */

  for (i = 0; i < model_len; i++)
  {
    const gdouble lb = ncm_model_param_get_lower_bound (model, i);
    const gdouble ub = ncm_model_param_get_upper_bound (model, i);
    gdouble val;

    while ((val = g_test_rand_double_range (lb, ub)))
    {
      if (val != ncm_model_param_get (model, i))
        break;
    }

    ncm_model_param_set (model, i, val);
    ncm_assert_cmpdouble (ncm_model_param_get (model, i), ==, val);

    ncm_model_param_set_default (model, i);
    ncm_assert_cmpdouble (ncm_model_param_get (model, i), !=, val);

    ncm_model_param_set (model, i, val);
    ncm_vector_set (tmp, i, val);
  }

  /* ncm_model_params_save_as_default */
  ncm_model_params_save_as_default (model);

  for (i = 0; i < model_len; i++)
  {
    ncm_assert_cmpdouble (ncm_model_param_get (model, i), ==, ncm_vector_get (tmp, i));
  }

  for (i = 0; i < test->sparam_len; i++)
  {
    ncm_model_param_set (model, i, s_defval_tot[i]);
  }

  for (i = 0; i < test->vparam_len; i++)
  {
    guint j;

    for (j = 0; j < v_len_tot[i]; j++)
    {
      ncm_model_param_set (model, ncm_model_vparam_index (model, i, j), v_defval_tot[i]);
    }
  }

  ncm_model_params_set_default (model);

  for (i = 0; i < model_len; i++)
  {
    ncm_assert_cmpdouble (ncm_model_param_get (model, i), ==, ncm_vector_get (tmp, i));
  }

  for (i = 0; i < test->sparam_len; i++)
  {
    ncm_model_param_set (model, i, s_defval_tot[i]);
  }

  for (i = 0; i < test->vparam_len; i++)
  {
    guint j;

    for (j = 0; j < v_len_tot[i]; j++)
    {
      ncm_model_param_set (model, ncm_model_vparam_index (model, i, j), v_defval_tot[i]);
    }
  }

  ncm_model_params_save_as_default (model);

  /* ncm_model_params_set_all_data */
  ncm_model_params_set_all_data (model, ncm_vector_data (tmp));

  for (i = 0; i < model_len; i++)
  {
    ncm_assert_cmpdouble (ncm_model_param_get (model, i), ==, ncm_vector_get (tmp, i));
  }

  ncm_vector_free (tmp);
}

void
test_ncm_model_test_setget_prop (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model  = NCM_MODEL (tm);
  guint i;

  /* ncm_model_param_get / ncm_model_param_set */

  for (i = 0; i < test->sparam_len; i++)
  {
    gdouble lb = ncm_model_orig_param_get_lower_bound (model, i);
    gdouble ub = ncm_model_orig_param_get_upper_bound (model, i);
    gdouble val, val_out;

    while ((val = g_test_rand_double_range (lb, ub)))
    {
      if (val != ncm_model_orig_param_get (model, i))
        break;
    }

    g_object_set (model, s_name_tot[i], val, NULL);
    g_object_get (model, s_name_tot[i], &val_out, NULL);
    ncm_assert_cmpdouble (val_out, ==, val);
  }

  for (i = 0; i < test->vparam_len; i++)
  {
    NcmVector *tmp     = ncm_vector_new (v_len_tot[i]);
    NcmVector *tmp_out = NULL;
    guint j;

    for (j = 0; j < v_len_tot[i]; j++)
    {
      guint n    = ncm_model_vparam_index (model, i, j);
      gdouble lb = ncm_model_orig_param_get_lower_bound (model, n);
      gdouble ub = ncm_model_orig_param_get_upper_bound (model, n);
      gdouble val;

      while ((val = g_test_rand_double_range (lb, ub)))
      {
        if (val != ncm_model_orig_param_get (model, n))
          break;
      }

      ncm_vector_set (tmp, j, val);
    }

    tmp_out = NULL;
    g_object_set (model, v_name_tot[i], tmp, NULL);
    g_object_get (model, v_name_tot[i], &tmp_out, NULL);

    for (j = 0; j < v_len_tot[i]; j++)
    {
      guint n = ncm_model_vparam_index (model, i, j);

      ncm_assert_cmpdouble (ncm_vector_get (tmp, j), ==, ncm_model_orig_param_get (model, n));
      ncm_assert_cmpdouble (ncm_vector_get (tmp, j), ==, ncm_vector_get (tmp_out, j));
    }

    NCM_TEST_FREE (ncm_vector_free, tmp);
    NCM_TEST_FREE (ncm_vector_free, tmp_out);
  }
}

void
test_ncm_model_test_setget_vector (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model  = NCM_MODEL (tm);
  guint model_len  = ncm_model_len (model);

  NCM_TEST_FAIL (G_STMT_START {
    NcmVector *tmp2 = ncm_vector_new (model_len + 1);
    ncm_model_params_set_vector (model, tmp2);
    ncm_vector_free (tmp2);
  } G_STMT_END);
}

void
test_ncm_model_test_setget_model (TestNcmModel *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_global ();
  NcmModelTest *tm1 = test->tm;
  NcmModel *model1  = NCM_MODEL (tm1);
  NcmModel *model2  = ncm_model_dup (model1, ser);
  NcmModelTest *tm3 = g_object_new (test->type,
                                    "VPBase0-length", ncm_model_vparam_len (model1, 0) + 1,
                                    "VPBase1-length", ncm_model_vparam_len (model1, 1) + 1,
                                    NULL);
  NcmModel *model3 = NCM_MODEL (tm3);
  guint model_len  = ncm_model_len (model1);
  guint i;

  ncm_serialize_free (ser);

  g_assert_true (ncm_model_is_equal (model1, model2));
  g_assert_true (!ncm_model_is_equal (model1, model3));

  for (i = 0; i < model_len; i++)
  {
    gdouble lb = ncm_model_param_get_lower_bound (model2, i);
    gdouble ub = ncm_model_param_get_upper_bound (model2, i);
    gdouble val;

    while ((val = g_test_rand_double_range (lb, ub)))
    {
      if (val != ncm_model_param_get (model2, i))
        break;
    }

    ncm_model_param_set (model2, i, val);
  }

  ncm_model_params_set_model (model1, model2);

  for (i = 0; i < model_len; i++)
  {
    ncm_assert_cmpdouble (ncm_model_param_get (model1, i), ==, ncm_model_param_get (model2, i));
  }

  g_assert_true (!ncm_model_is_equal (model1, model3));
  g_assert_true (!ncm_model_is_equal (model3, model1));

  NCM_TEST_FREE (ncm_model_free, model2);
  NCM_TEST_FREE (ncm_model_free, model3);
}

void
test_ncm_model_test_name_index (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model  = NCM_MODEL (tm);
  guint i;

  for (i = 0; i < test->sparam_len; i++)
  {
    guint n;
    gboolean found          = ncm_model_orig_param_index_from_name (model, s_name_tot[i], &n);
    const gchar *s_name_n   = ncm_model_orig_param_name (model, n);
    const gchar *s_symbol_n = ncm_model_orig_param_symbol (model, n);

    g_assert_true (found);
    g_assert_cmpuint (n, ==, i);
    g_assert_cmpstr (s_name_n, ==, s_name_tot[n]);
    g_assert_cmpstr (s_symbol_n, ==, s_symbol_tot[n]);
  }

  for (i = 0; i < test->vparam_len; i++)
  {
    guint j;

    for (j = 0; j < v_len_tot[i]; j++)
    {
      gchar *v_name_ij   = g_strdup_printf ("%s_%u", v_name_tot[i], j);
      gchar *v_symbol_ij = g_strdup_printf ("{%s}_%u", v_symbol_tot[i], j);
      guint n;
      gboolean found          = ncm_model_orig_param_index_from_name (model, v_name_ij, &n);
      const gchar *v_name_n   = ncm_model_orig_param_name (model, n);
      const gchar *v_symbol_n = ncm_model_orig_param_symbol (model, n);

      g_assert_true (found);
      g_assert_cmpstr (v_name_n, ==, v_name_ij);
      g_assert_cmpstr (v_symbol_n, ==, v_symbol_ij);

      g_free (v_name_ij);
      g_free (v_symbol_ij);
    }
  }
}

void
test_ncm_model_test_dup (TestNcmModel *test, gconstpointer pdata)
{
  NcmSerialize *ser   = ncm_serialize_global ();
  NcmModelTest *tm    = test->tm;
  NcmModel *model     = NCM_MODEL (tm);
  NcmModel *model_dup = ncm_model_dup (model, ser);
  guint model_len     = ncm_model_len (model);
  guint i;

  ncm_serialize_free (ser);
  g_assert_true (ncm_model_is_equal (model, model_dup));

  for (i = 0; i < model_len; i++)
  {
    ncm_assert_cmpdouble (ncm_model_param_get (model, i),             ==, ncm_model_param_get (model_dup, i));
    ncm_assert_cmpdouble (ncm_model_param_get_scale (model, i),       ==, ncm_model_param_get_scale (model_dup, i));
    ncm_assert_cmpdouble (ncm_model_param_get_lower_bound (model, i), ==, ncm_model_param_get_lower_bound (model_dup, i));
    ncm_assert_cmpdouble (ncm_model_param_get_upper_bound (model, i), ==, ncm_model_param_get_upper_bound (model_dup, i));
    ncm_assert_cmpdouble (ncm_model_param_get_abstol (model, i),      ==, ncm_model_param_get_abstol (model_dup, i));

    ncm_assert_cmpdouble (ncm_model_orig_param_get (model, i),             ==, ncm_model_orig_param_get (model_dup, i));
    ncm_assert_cmpdouble (ncm_model_orig_param_get_scale (model, i),       ==, ncm_model_orig_param_get_scale (model_dup, i));
    ncm_assert_cmpdouble (ncm_model_orig_param_get_lower_bound (model, i), ==, ncm_model_orig_param_get_lower_bound (model_dup, i));
    ncm_assert_cmpdouble (ncm_model_orig_param_get_upper_bound (model, i), ==, ncm_model_orig_param_get_upper_bound (model_dup, i));
    ncm_assert_cmpdouble (ncm_model_orig_param_get_abstol (model, i),      ==, ncm_model_orig_param_get_abstol (model_dup, i));
  }

  ncm_model_free (model_dup);
}

void
test_ncm_model_test_impl (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model  = NCM_MODEL (tm);

  g_assert_true (ncm_model_check_impl_flag (model, 1 << 0));
  g_assert_true (ncm_model_check_impl_flag (model, 1 << 1));
  g_assert_true (ncm_model_check_impl_flag (model, 1 << 2));
  g_assert_false (ncm_model_check_impl_flag (model, 1 << 3));
  g_assert_false (ncm_model_check_impl_flag (model, 1 << 4));
  g_assert_false (ncm_model_check_impl_flag (model, 1 << 5));

  g_assert_true (ncm_model_check_impl_opt (model, 0));
  g_assert_true (ncm_model_check_impl_opt (model, 1));
  g_assert_true (ncm_model_check_impl_opt (model, 2));
  g_assert_false (ncm_model_check_impl_opt (model, 3));
  g_assert_false (ncm_model_check_impl_opt (model, 4));
  g_assert_false (ncm_model_check_impl_opt (model, 5));

  g_assert_true (ncm_model_check_impl_opts (model, 0, 1, 2, -1));
  g_assert_true (ncm_model_check_impl_opts (model, 2, 1, 2, -1));
  g_assert_false (ncm_model_check_impl_opts (model, 2, 1, 4, -1));

  {
    gint64 flags = 0;

    g_object_get (model, "implementation", &flags, NULL);
    g_assert_true (flags & (1 << 0));
    g_assert_true (flags & (1 << 1));
    g_assert_true (flags & (1 << 2));
    g_assert_false (flags & (1 << 3));
    g_assert_false (flags & (1 << 4));
    g_assert_false (flags & (1 << 5));
  }
}

void
test_ncm_model_param_names (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model  = NCM_MODEL (tm);
  GPtrArray *names = ncm_model_param_names (model);
  guint nnames     = test->sparam_len;
  guint i;

  for (i = 0; i < test->vparam_len; i++)
  {
    nnames += v_len_tot[i];
  }

  g_assert_cmpuint (names->len, ==, nnames);

  {
    guint name_index = 0;
    guint i;

    for (i = 0; i < test->sparam_len; i++)
    {
      const gchar *name = g_ptr_array_index (names, i);

      g_assert_cmpstr (name, ==, ncm_model_param_name (model, i));
    }

    name_index += test->sparam_len;

    for (i = 0; i < test->vparam_len; i++)
    {
      guint j;

      for (j = 0; j < v_len_tot[i]; j++)
      {
        const gchar *name   = g_ptr_array_index (names, name_index + j);
        const gchar *v_name = ncm_model_param_name (model, name_index + j);

        g_assert_cmpstr (name, ==, v_name);
      }

      name_index += v_len_tot[i];
    }
  }

  g_ptr_array_unref (names);
}

void
test_ncm_model_svparams_len (TestNcmModel *test, gconstpointer pdata)
{
  NcmModelTest *tm = test->tm;
  NcmModel *model  = NCM_MODEL (tm);
  guint sparam_len, vparam_len;

  g_object_get (model, "scalar-params-len", &sparam_len, NULL);
  g_object_get (model, "vector-params-len", &vparam_len, NULL);

  g_assert_cmpuint (sparam_len, ==, test->sparam_len);
  g_assert_cmpuint (vparam_len, ==, test->vparam_len);
}

/* A NcmModelMVND of dimension 2 under p_n = T p with T = [[1, 0.5], [0, 1]]. */
static NcmReparam *
_test_ncm_model_reparam_linear (void)
{
  NcmMatrix *T = ncm_matrix_new (2, 2);
  NcmVector *v = ncm_vector_new (2);
  NcmReparam *reparam;

  ncm_matrix_set_identity (T);
  ncm_matrix_set (T, 0, 1, 0.5);
  ncm_vector_set_zero (v);

  reparam = NCM_REPARAM (ncm_reparam_linear_new (2, T, v));

  ncm_matrix_free (T);
  ncm_vector_free (v);

  return reparam;
}

/* The param functions, the finite checks included, read the new parameters. */
void
test_ncm_model_reparam_current (void)
{
  NcmModel *model     = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmReparam *reparam = _test_ncm_model_reparam_linear ();

  ncm_model_orig_param_set (model, 0, 1.0);
  ncm_model_orig_param_set (model, 1, 2.0);
  ncm_model_set_reparam (model, reparam, NULL);
  g_assert_true (ncm_model_peek_reparam (model) == reparam);

  ncm_assert_cmpdouble_e (ncm_model_param_get (model, 0), ==, 2.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (ncm_model_param_get (model, 1), ==, 2.0, 1.0e-15, 0.0);

  ncm_model_param_set (model, 0, 3.0);
  ncm_assert_cmpdouble_e (ncm_model_orig_param_get (model, 0), ==, 2.0, 1.0e-15, 0.0);

  /* A non-finite new parameter set without update leaves the original ones finite. */
  ncm_model_param_set0 (model, 1, GSL_NAN);
  g_assert_false (ncm_model_param_finite (model, 1));
  g_assert_true (ncm_model_param_finite (model, 0));
  g_assert_false (ncm_model_params_finite (model));
  g_assert_true (gsl_finite (ncm_model_orig_param_get (model, 1)));

  ncm_reparam_free (reparam);
  ncm_model_free (model);
}

/* is_equal is symmetric: a model with a reparametrization differs from one without. */
void
test_ncm_model_reparam_is_equal (void)
{
  NcmModel *a         = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmModel *b         = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmModel *c         = NCM_MODEL (ncm_model_mvnd_new (3));
  NcmReparam *reparam = _test_ncm_model_reparam_linear ();

  g_assert_true (ncm_model_is_equal (a, b));
  g_assert_false (ncm_model_is_equal (a, c));

  ncm_model_set_reparam (b, reparam, NULL);
  g_assert_false (ncm_model_is_equal (a, b));
  g_assert_false (ncm_model_is_equal (b, a));

  ncm_model_set_reparam (a, reparam, NULL);
  g_assert_true (ncm_model_is_equal (a, b));

  ncm_reparam_free (reparam);
  ncm_model_free (a);
  ncm_model_free (b);
  ncm_model_free (c);
}

/* NULL does nothing on a model without a reparametrization; removing one aborts. */
void
test_ncm_model_reparam_remove (void)
{
  NcmModel *model = NCM_MODEL (ncm_model_mvnd_new (2));

  ncm_model_set_reparam (model, NULL, NULL);
  g_assert_null (ncm_model_peek_reparam (model));
  ncm_model_free (model);

  g_test_trap_subprocess ("/ncm/model/reparam/remove/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*cannot be removed*");
}

void
test_ncm_model_reparam_remove_subprocess (void)
{
  NcmModel *model     = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmReparam *reparam = _test_ncm_model_reparam_linear ();

  ncm_model_set_reparam (model, reparam, NULL);
  ncm_model_set_reparam (model, NULL, NULL);
}

/* Setting a parameter to its default does not mark its description as modified. */
void
test_ncm_model_set_default_desc (void)
{
  NcmModel *model = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmObjDictInt *modified;

  ncm_model_orig_param_set (model, 0, 0.7);
  ncm_model_param_set_default (model, 0);
  ncm_assert_cmpdouble_e (ncm_model_orig_param_get (model, 0), ==, ncm_model_param_get (model, 0), 0.0, 0.0);

  g_object_get (model, "sparam-array", &modified, NULL);
  g_assert_cmpuint (ncm_obj_dict_int_len (modified), ==, 0);
  ncm_obj_dict_int_unref (modified);

  ncm_model_params_save_as_default (model);
  g_object_get (model, "sparam-array", &modified, NULL);
  g_assert_cmpuint (ncm_obj_dict_int_len (modified), ==, 2);
  ncm_obj_dict_int_unref (modified);

  ncm_model_free (model);
}

/* A type that is not a model is an error, with -1 as the id. */
void
test_ncm_model_id_by_type_error (void)
{
  GError *error = NULL;

  g_assert_cmpint (ncm_model_id_by_type (G_TYPE_OBJECT, &error), ==, -1);
  g_assert_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_INVALID_TYPE);
  g_clear_error (&error);
}

static GValue *
_test_ncm_model_gvalue_double (const gdouble x)
{
  GValue *value = g_new0 (GValue, 1);

  g_value_init (value, G_TYPE_DOUBLE);
  g_value_set_double (value, x);

  return value;
}

static void
_test_ncm_model_gvalue_free (gpointer data)
{
  g_value_unset (data);
  g_free (data);
}

/* get_desc reports a parameter and set_desc changes it; unknown keys are one error. */
void
test_ncm_model_param_desc (void)
{
  NcmModel *model = NCM_MODEL (ncm_model_mvnd_new (2));
  GError *error   = NULL;
  GHashTable *desc, *set;

  ncm_model_param_set (model, 1, 0.25);
  desc = ncm_model_param_get_desc (model, "mu_1", &error);
  g_assert_no_error (error);
  g_assert_cmpstr (g_value_get_string (g_hash_table_lookup (desc, "name")), ==, "mu_1");
  g_assert_cmpfloat (g_value_get_double (g_hash_table_lookup (desc, "value")), ==, 0.25);
  g_assert_cmpfloat (g_value_get_double (g_hash_table_lookup (desc, "upper-bound")), ==, ncm_model_param_get_upper_bound (model, 1));
  g_assert_false (g_value_get_boolean (g_hash_table_lookup (desc, "fit")));
  g_hash_table_unref (desc);

  set = g_hash_table_new_full (g_str_hash, g_str_equal, NULL, _test_ncm_model_gvalue_free);
  g_hash_table_insert (set, "scale", _test_ncm_model_gvalue_double (0.125));
  g_hash_table_insert (set, "value", _test_ncm_model_gvalue_double (0.5));
  ncm_model_param_set_desc (model, "mu_1", set, &error);
  g_assert_no_error (error);
  g_assert_cmpfloat (ncm_model_param_get_scale (model, 1), ==, 0.125);
  g_assert_cmpfloat (ncm_model_param_get (model, 1), ==, 0.5);

  g_hash_table_insert (set, "bogus", _test_ncm_model_gvalue_double (1.0));
  g_hash_table_insert (set, "name", _test_ncm_model_gvalue_double (1.0));
  ncm_model_param_set_desc (model, "mu_1", set, &error);
  g_assert_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_INVALID_KEY);
  g_clear_error (&error);

  g_hash_table_unref (set);
  ncm_model_free (model);
}

/* An original name renamed by a submodel's reparametrization is an error, qualified or
 * not, never an abort. */
void
test_ncm_model_renamed_submodel_param (void)
{
  NcHIReion *reion     = NC_HIREION (nc_hireion_camb_new ());
  NcHICosmoDEXcdm *cde = nc_hicosmo_de_xcdm_new_full (reion, NULL, NULL);
  NcmModel *cosmo      = NCM_MODEL (cde);
  NcmReparam *tau      = NCM_REPARAM (nc_hireion_camb_reparam_tau_new (ncm_model_len (NCM_MODEL (reion))));
  GError *error        = NULL;

  ncm_model_set_reparam (NCM_MODEL (reion), tau, NULL);

  ncm_model_param_get_by_name (cosmo, "reion:z_re", &error);
  g_assert_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_CHANGED);
  g_clear_error (&error);

  ncm_model_param_get_by_name (cosmo, "z_re", &error);
  g_assert_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_CHANGED);
  g_clear_error (&error);

  g_assert_true (gsl_finite (ncm_model_param_get_by_name (cosmo, "tau_reion", &error)));
  g_assert_no_error (error);

  ncm_reparam_free (tau);
  ncm_model_free (cosmo);
  nc_hireion_free (reion);
}

static void
_test_ncm_model_vparam_set_vector_len (void)
{
  NcmModel *model = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmVector *v    = ncm_vector_new (3);

  ncm_vector_set_zero (v);
  ncm_model_orig_vparam_set_vector (model, 0, v);
}

void
test_ncm_model_vparam_set_vector_len (void)
{
  g_test_trap_subprocess ("/ncm/model/vparam_set_vector_len/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the vector has 3 elements but the vector parameter 0 has 2*");
}

void
test_ncm_model_vparam_set_vector_len_subprocess (void)
{
  _test_ncm_model_vparam_set_vector_len ();
}

