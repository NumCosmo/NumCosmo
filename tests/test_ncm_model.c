/***************************************************************************
 *            test_ncm_model.c
 *
 *  Mon May 7 16:39:21 2012
 *  Copyright  2012  Mariana Penna Lima & Sandro Dias Pinto Vitenti
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima & Sandro Dias Pinto Vitenti 2012 <pennalima@gmail.com>, <sandro@isoftware.com.br>
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

void test_ncm_model_test_type (void);
void test_ncm_model_test_child_type (void);
void test_ncm_model_test_child_child_type (void);
void test_ncm_model_test_new (void);
void test_ncm_model_test_length (void);
void test_ncm_model_test_defval (void);
void test_ncm_model_test_name_symbol (void);
void test_ncm_model_test_setget (void);
void test_ncm_model_test_setget_prop (void);
void test_ncm_model_test_setget_vector (void);
void test_ncm_model_test_setget_model (void);
void test_ncm_model_test_name_index (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/numcosmo/ncm_model/type", &test_ncm_model_test_type);
  g_test_add_func ("/numcosmo/ncm_model/new", &test_ncm_model_test_new);
  g_test_add_func ("/numcosmo/ncm_model/length", &test_ncm_model_test_length);
  g_test_add_func ("/numcosmo/ncm_model/defval", &test_ncm_model_test_defval);
  g_test_add_func ("/numcosmo/ncm_model/name_symbol", &test_ncm_model_test_name_symbol);
  g_test_add_func ("/numcosmo/ncm_model/setget", &test_ncm_model_test_setget);
  g_test_add_func ("/numcosmo/ncm_model/setget/prop", &test_ncm_model_test_setget_prop);
  g_test_add_func ("/numcosmo/ncm_model/setget/vector", &test_ncm_model_test_setget_vector);
  g_test_add_func ("/numcosmo/ncm_model/setget/model", &test_ncm_model_test_setget_model);
  g_test_add_func ("/numcosmo/ncm_model/name_index", &test_ncm_model_test_name_index);

  g_test_add_func ("/numcosmo/ncm_model/child/type", &test_ncm_model_test_child_type);
  g_test_add_func ("/numcosmo/ncm_model/child/new", &test_ncm_model_test_new);
  g_test_add_func ("/numcosmo/ncm_model/child/length", &test_ncm_model_test_length);
  g_test_add_func ("/numcosmo/ncm_model/child/defval", &test_ncm_model_test_defval);
  g_test_add_func ("/numcosmo/ncm_model/child/name_symbol", &test_ncm_model_test_name_symbol);
  g_test_add_func ("/numcosmo/ncm_model/child/setget", &test_ncm_model_test_setget);
  g_test_add_func ("/numcosmo/ncm_model/child/setget/prop", &test_ncm_model_test_setget_prop);
  g_test_add_func ("/numcosmo/ncm_model/child/setget/vector", &test_ncm_model_test_setget_vector);
  g_test_add_func ("/numcosmo/ncm_model/child/setget/model", &test_ncm_model_test_setget_model);
  g_test_add_func ("/numcosmo/ncm_model/child/name_index", &test_ncm_model_test_name_index);

  g_test_add_func ("/numcosmo/ncm_model/grandchild/type", &test_ncm_model_test_child_child_type);
  g_test_add_func ("/numcosmo/ncm_model/grandchild/new", &test_ncm_model_test_new);
  g_test_add_func ("/numcosmo/ncm_model/grandchild/length", &test_ncm_model_test_length);
  g_test_add_func ("/numcosmo/ncm_model/grandchild/defval", &test_ncm_model_test_defval);
  g_test_add_func ("/numcosmo/ncm_model/grandchild/name_symbol", &test_ncm_model_test_name_symbol);
  g_test_add_func ("/numcosmo/ncm_model/grandchild/setget", &test_ncm_model_test_setget);
  g_test_add_func ("/numcosmo/ncm_model/grandchild/setget/prop", &test_ncm_model_test_setget_prop);
  g_test_add_func ("/numcosmo/ncm_model/grandchild/setget/vector", &test_ncm_model_test_setget_vector);
  g_test_add_func ("/numcosmo/ncm_model/grandchild/setget/model", &test_ncm_model_test_setget_model);
  g_test_add_func ("/numcosmo/ncm_model/grandchild/name_index", &test_ncm_model_test_name_index);

  g_test_run ();
}

static GType test_type;

static guint sparam_len, vparam_len;
gchar *name, *nick;

void
test_ncm_model_test_type (void)
{
  test_type = NCM_TYPE_MODEL_TEST;
  {
    NcmModelTest *tm = g_object_new (test_type, NULL);
    ncm_model_free (NCM_MODEL (tm));
  }
  sparam_len = SPARAM_LEN1;
  vparam_len = VPARAM_LEN1;
  name = name_tot[0];
  nick = nick_tot[0];

  g_assert (test_type != 0);
}

void
test_ncm_model_test_child_type (void)
{
  test_type = NCM_TYPE_MODEL_TEST_CHILD;
  {
    NcmModelTest *tm = g_object_new (test_type, NULL);
    ncm_model_free (NCM_MODEL (tm));
  }
  sparam_len += SPARAM_LEN2;
  vparam_len += VPARAM_LEN2;
  name = name_tot[1];
  nick = nick_tot[1];

  g_assert (test_type != 0);
}

void
test_ncm_model_test_child_child_type (void)
{
  test_type = NCM_TYPE_MODEL_TEST_CHILD_CHILD;
  {
    NcmModelTest *tm = g_object_new (test_type, NULL);
    ncm_model_free (NCM_MODEL (tm));
  }
  sparam_len += SPARAM_LEN3;
  vparam_len += VPARAM_LEN3;
  name = name_tot[2];
  nick = nick_tot[2];

  g_assert (test_type != 0);
}

void
test_ncm_model_test_new (void)
{
  NcmModelTest *tm = g_object_new (test_type, NULL);
  NcmModel *model;
  g_assert (NCM_IS_MODEL (tm));
  g_assert (NCM_IS_MODEL_TEST (tm));
  model = NCM_MODEL (tm);

  g_assert (ncm_model_name (model) != name);
  g_assert (ncm_model_nick (model) != nick);

  g_assert_cmpstr (ncm_model_name (model), ==, name);
  g_assert_cmpstr (ncm_model_nick (model), ==, nick);
  g_assert (ncm_model_impl (model) == 0);
  g_assert (ncm_model_peek_reparam (model) == NULL);

  g_object_unref (tm);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    g_object_unref (tm);
    exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_model_test_length (void)
{
  NcmModelTest *tm = g_object_new (test_type, NULL);
  NcmModel *model = NCM_MODEL (tm);
  guint model_len = 0;
  guint i;

  model_len = sparam_len;
  for (i = 0; i < vparam_len; i++)
  {
    g_assert_cmpint (ncm_model_vparam_len (model, i), ==, v_len_tot[i]);
    model_len += v_len_tot[i];
  }

  g_assert_cmpint (ncm_model_len (model), ==, model_len);
  g_assert_cmpint (ncm_model_sparam_len (model), ==, sparam_len);
  g_assert_cmpint (ncm_model_vparam_array_len (model), ==, vparam_len);

  g_object_unref (tm);
}

void
test_ncm_model_test_defval (void)
{
  NcmModelTest *tm = g_object_new (test_type, NULL);
  NcmModel *model = NCM_MODEL (tm);
  guint i;

  for (i = 0; i < sparam_len; i++)
  {
    gdouble p = ncm_model_orig_param_get (model, i);
    ncm_assert_cmpdouble (p, ==, s_defval_tot[i]);
  }

  for (i = 0; i < vparam_len; i++)
  {
    guint j;
    for (j = 0; j < v_len_tot[i]; j++)
    {
      gdouble p = ncm_model_orig_vparam_get (model, i, j);
      ncm_assert_cmpdouble (p, ==, v_defval_tot[i]);
    }
  }

  g_object_unref (tm);
}

void
test_ncm_model_test_name_symbol (void)
{
  NcmModelTest *tm = g_object_new (test_type, NULL);
  NcmModel *model = NCM_MODEL (tm);
  guint i;

  for (i = 0; i < sparam_len; i++)
  {
    const gchar *pname = ncm_model_param_name (model, i);
    const gchar *psymbol = ncm_model_param_symbol (model, i);
    g_assert (pname != s_name_tot[i]);
    g_assert (psymbol != s_symbol_tot[i]);
    g_assert_cmpstr (pname, ==, s_name_tot[i]);
    g_assert_cmpstr (psymbol, ==, s_symbol_tot[i]);
  }

  for (i = 0; i < vparam_len; i++)
  {
    guint j;
    for (j = 0; j < v_len_tot[i]; j++)
    {
      gchar *vp_name = g_strdup_printf ("%s_%u", v_name_tot[i], j);
      gchar *vp_symbol = g_strdup_printf ("%s_%u", v_symbol_tot[i], j);
      const gchar *pname = ncm_model_param_name (model, ncm_model_vparam_index (model, i, j));
      const gchar *psymbol = ncm_model_param_symbol (model, ncm_model_vparam_index (model, i, j));

      g_assert (pname != vp_name);
      g_assert (psymbol != vp_symbol);
      g_assert_cmpstr (pname, ==, vp_name);
      g_assert_cmpstr (psymbol, ==, vp_symbol);

      g_free (vp_name);
      g_free (vp_symbol);
    }
  }

  g_object_unref (tm);
}

void
test_ncm_model_test_finite (void)
{
  NcmModelTest *tm = g_object_new (test_type, NULL);
  NcmModel *model = NCM_MODEL (tm);
  guint model_len = ncm_model_len (model);
  guint i;

  for (i = 0; i < model_len; i++)
  {
    g_assert (ncm_model_param_finite (model, i));
    g_assert (ncm_model_params_finite (model));

    ncm_model_param_set (model, i, GSL_NAN);
    g_assert (!ncm_model_param_finite (model, i));
    g_assert (!ncm_model_params_finite (model));
    ncm_model_param_set_default (model, i);
    g_assert (ncm_model_param_finite (model, i));
    g_assert (ncm_model_params_finite (model));

    ncm_model_param_set (model, i, GSL_POSINF);
    g_assert (!ncm_model_param_finite (model, i));
    g_assert (!ncm_model_params_finite (model));
    ncm_model_param_set_default (model, i);
    g_assert (ncm_model_param_finite (model, i));
    g_assert (ncm_model_params_finite (model));

    ncm_model_param_set (model, i, GSL_NEGINF);
    g_assert (!ncm_model_param_finite (model, i));
    g_assert (!ncm_model_params_finite (model));
    ncm_model_param_set_default (model, i);
    g_assert (ncm_model_param_finite (model, i));
    g_assert (ncm_model_params_finite (model));
  }

  g_object_unref (tm);
}

void
test_ncm_model_test_setget (void)
{
  NcmModelTest *tm = g_object_new (test_type, NULL);
  NcmModel *model = NCM_MODEL (tm);
  guint model_len = ncm_model_len (model);
  NcmVector *tmp = ncm_vector_new (model_len);
  guint i;

  /* ncm_model_param_get / ncm_model_param_set */

  for (i = 0; i < model_len; i++)
  {
    gdouble lb = ncm_model_param_get_lower_bound (model, i);
    gdouble ub = ncm_model_param_get_upper_bound (model, i);
    gdouble val;
    while ((val = g_test_rand_double_range (lb, ub)))
      if (val != ncm_model_param_get (model, i))
      break;

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

  for (i = 0; i < sparam_len; i++)
  {
    ncm_model_param_set (model, i, s_defval_tot[i]);
  }
  for (i = 0; i < vparam_len; i++)
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

  for (i = 0; i < sparam_len; i++)
  {
    ncm_model_param_set (model, i, s_defval_tot[i]);
  }
  for (i = 0; i < vparam_len; i++)
  {
    guint j;
    for (j = 0; j < v_len_tot[i]; j++)
    {
      ncm_model_param_set (model, ncm_model_vparam_index (model, i, j), v_defval_tot[i]);
    }
  }
  ncm_model_params_save_as_default (model);

  /* ncm_model_params_set_all_data */
  ncm_model_params_set_all_data (model, NCM_VECTOR_DATA (tmp));
  for (i = 0; i < model_len; i++)
  {
    ncm_assert_cmpdouble (ncm_model_param_get (model, i), ==, ncm_vector_get (tmp, i));
  }

  ncm_vector_free (tmp);
  g_object_unref (tm);
}

void
test_ncm_model_test_setget_prop (void)
{
  NcmModelTest *tm = g_object_new (test_type, NULL);
  NcmModel *model = NCM_MODEL (tm);
  guint i;

  /* ncm_model_param_get / ncm_model_param_set */

  for (i = 0; i < sparam_len; i++)
  {
    gdouble lb = ncm_model_param_get_lower_bound (model, i);
    gdouble ub = ncm_model_param_get_upper_bound (model, i);
    gdouble val, val_out;
    while ((val = g_test_rand_double_range (lb, ub)))
    {
      if (val != ncm_model_param_get (model, i))
        break;
    }
    
    g_object_set (model, s_name_tot[i], val, NULL);
    g_object_get (model, s_name_tot[i], &val_out, NULL);
    ncm_assert_cmpdouble (val_out, ==, val);
  }

  for (i = 0; i < vparam_len; i++)
  {
    NcmVector *tmp = ncm_vector_new (v_len_tot[i]);
    NcmVector *tmp_out = NULL;
    guint j;
    for (j = 0; j < v_len_tot[i]; j++)
    {
      guint n = ncm_model_vparam_index (model, i, j);
      gdouble lb = ncm_model_param_get_lower_bound (model, n);
      gdouble ub = ncm_model_param_get_upper_bound (model, n);
      gdouble val;

      while ((val = g_test_rand_double_range (lb, ub)))
      {
        if (val != ncm_model_param_get (model, n))
          break;
      }
      ncm_vector_set (tmp, j, val);
    }

    {
      GVariant *var = ncm_vector_get_variant (tmp);
      GVariant *var_out = NULL;
      g_object_set (model, v_name_tot[i], var, NULL);
      g_object_get (model, v_name_tot[i], &var_out, NULL);

      tmp_out = ncm_vector_new_variant (var_out);
      g_variant_unref (var);
      g_variant_unref (var_out);
    }

    for (j = 0; j < v_len_tot[i]; j++)
    {
      guint n = ncm_model_vparam_index (model, i, j);
      ncm_assert_cmpdouble (ncm_vector_get (tmp, j), ==, ncm_model_param_get (model, n));
      ncm_assert_cmpdouble (ncm_vector_get (tmp, j), ==, ncm_vector_get (tmp_out, j));
    }
    ncm_vector_free (tmp);
    if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
    {
      ncm_vector_free (tmp);
      exit (0);
    }
    g_test_trap_assert_failed ();
    ncm_vector_free (tmp_out);
    if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
    {
      ncm_vector_free (tmp_out);
      exit (0);
    }
    g_test_trap_assert_failed ();
  }

  g_object_unref (tm);
}

void
test_ncm_model_test_setget_vector (void)
{
  NcmModelTest *tm = g_object_new (test_type, NULL);
  NcmModel *model = NCM_MODEL (tm);
  guint model_len = ncm_model_len (model);

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    NcmVector *tmp2 = ncm_vector_new (model_len + 1);
    ncm_model_params_set_vector (model, tmp2);
    ncm_vector_free (tmp2);
    exit (0);
  }
  g_test_trap_assert_failed ();

  g_object_unref (tm);
}

void
test_ncm_model_test_setget_model (void)
{
  NcmModelTest *tm1 = g_object_new (test_type, NULL);
  NcmModel *model1 = NCM_MODEL (tm1);
  NcmModel *model2 = ncm_model_copy (model1);
  NcmModelTest *tm3 = g_object_new (test_type,
                                    "VPBase0-length", ncm_model_vparam_len (model1, 0) + 1,
                                    "VPBase1-length", ncm_model_vparam_len (model1, 1) + 1,
                                    NULL);
  NcmModel *model3 = NCM_MODEL (tm3);
  guint model_len = ncm_model_len (model1);
  guint i;

  g_assert (ncm_model_is_equal (model1, model2));
  g_assert (!ncm_model_is_equal (model1, model3));

  for (i = 0; i < model_len; i++)
  {
    gdouble lb = ncm_model_param_get_lower_bound (model2, i);
    gdouble ub = ncm_model_param_get_upper_bound (model2, i);
    gdouble val;
    while ((val = g_test_rand_double_range (lb, ub)))
      if (val != ncm_model_param_get (model2, i))
      break;

    ncm_model_param_set (model2, i, val);
  }

  ncm_model_params_set_model (model1, model2);

  for (i = 0; i < model_len; i++)
  {
    ncm_assert_cmpdouble (ncm_model_param_get (model1, i), ==, ncm_model_param_get (model2, i));
  }

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_model_params_set_model (model1, model3);
    exit (0);
  }
  g_test_trap_assert_failed ();

  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_model_params_set_model (model3, model1);
    exit (0);
  }
  g_test_trap_assert_failed ();

  ncm_model_free (model1);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_model_free (model1);
    exit (0);
  }
  g_test_trap_assert_failed ();

  ncm_model_free (model2);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_model_free (model2);
    exit (0);
  }
  g_test_trap_assert_failed ();

  ncm_model_free (model3);
  if (g_test_trap_fork (0, G_TEST_TRAP_SILENCE_STDOUT | G_TEST_TRAP_SILENCE_STDERR))
  {
    ncm_model_free (model3);
    exit (0);
  }
  g_test_trap_assert_failed ();
}

void
test_ncm_model_test_name_index (void)
{
  NcmModelTest *tm = g_object_new (test_type, NULL);
  NcmModel *model = NCM_MODEL (tm);
  guint i;

  for (i = 0; i < sparam_len; i++)
  {
    guint n = ncm_model_param_index_from_name (model, s_name_tot[i]);
    const gchar *s_name_n = ncm_model_param_name (model, n);
    const gchar *s_symbol_n = ncm_model_param_symbol (model, n);

    g_assert_cmpuint (n, ==, i);
    g_assert_cmpstr (s_name_n, ==, s_name_tot[n]);
    g_assert_cmpstr (s_symbol_n, ==, s_symbol_tot[n]);
  }

  for (i = 0; i < vparam_len; i++)
  {
    guint j;

    for (j = 0; j < v_len_tot[i]; j++)
    {
      gchar *v_name_ij = g_strdup_printf ("%s_%u", v_name_tot[i], j);
      gchar *v_symbol_ij = g_strdup_printf ("%s_%u", v_symbol_tot[i], j);
      guint n = ncm_model_param_index_from_name (model, v_name_ij);
      const gchar *v_name_n = ncm_model_param_name (model, n);
      const gchar *v_symbol_n = ncm_model_param_symbol (model, n);

      g_assert_cmpstr (v_name_n, ==, v_name_ij);
      g_assert_cmpstr (v_symbol_n, ==, v_symbol_ij);

      g_free (v_name_ij);
      g_free (v_symbol_ij);
    }
  }

  g_object_unref (tm);
}
