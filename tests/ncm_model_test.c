/***************************************************************************
 *            ncm_model_test.c
 *
 *  Fri May 18 16:00:21 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

G_DEFINE_TYPE (NcmModelTest, ncm_model_test, NCM_TYPE_MODEL);
G_DEFINE_TYPE (NcmModelTestChild, ncm_model_test_child, NCM_TYPE_MODEL_TEST);
G_DEFINE_TYPE (NcmModelTestChildChild, ncm_model_test_child_child, NCM_TYPE_MODEL_TEST_CHILD);

enum {
  PROP_0,
  PROP_A,
  PROP_SIZE,
};

static void
ncm_model_test_init (NcmModelTest *tm)
{
  tm->A = 0.0;
}

static void
ncm_model_test_child_init (NcmModelTestChild *tmc)
{
  tmc->B = 0.0;
}

static void
ncm_model_test_child_child_init (NcmModelTestChildChild *tmcc)
{
  tmcc->C = 0.0;
}

static void
ncm_model_test_finalize (GObject *object)
{
  G_OBJECT_CLASS (ncm_model_test_parent_class)->finalize (object);
}

static void
ncm_model_test_child_finalize (GObject *object)
{
  G_OBJECT_CLASS (ncm_model_test_child_parent_class)->finalize (object);
}

static void
ncm_model_test_child_child_finalize (GObject *object)
{
  G_OBJECT_CLASS (ncm_model_test_child_child_parent_class)->finalize (object);
}


static void
_ncm_model_test_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmModelTest *tm = NCM_MODEL_TEST (object);
  g_return_if_fail (NCM_IS_MODEL_TEST (object));

  switch (prop_id)
  {
    case PROP_A:
      g_value_set_double (value, tm->A);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_model_test_child_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmModelTestChild *tmc = NCM_MODEL_TEST_CHILD (object);
  g_return_if_fail (NCM_IS_MODEL_TEST_CHILD (object));

  switch (prop_id)
  {
    case PROP_A:
      g_value_set_double (value, tmc->B);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_model_test_child_child_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmModelTestChildChild *tmcc = NCM_MODEL_TEST_CHILD_CHILD (object);
  g_return_if_fail (NCM_IS_MODEL_TEST_CHILD_CHILD (object));

  switch (prop_id)
  {
    case PROP_A:
      g_value_set_double (value, tmcc->C);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_model_test_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmModelTest *tm = NCM_MODEL_TEST (object);
  g_return_if_fail (NCM_IS_MODEL_TEST (object));

  switch (prop_id)
  {
    case PROP_A:
      tm->A = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_model_test_child_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmModelTestChild *tmc = NCM_MODEL_TEST_CHILD (object);
  g_return_if_fail (NCM_IS_MODEL_TEST_CHILD (object));

  switch (prop_id)
  {
    case PROP_A:
      tmc->B = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_ncm_model_test_child_child_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmModelTestChildChild *tmcc = NCM_MODEL_TEST_CHILD_CHILD (object);
  g_return_if_fail (NCM_IS_MODEL_TEST_CHILD_CHILD (object));

  switch (prop_id)
  {
    case PROP_A:
      tmcc->C = g_value_get_double (value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

gchar *name_tot[3] = {"Test Model", "Test Model Child", "Test Model Grandchild"};
gchar *nick_tot[3] = {"TModel",     "TModelChild",      "TModelChildChild"};
gchar  *s_symbol_tot[SPARAM_LEN_TOT],        *v_symbol_tot[VPARAM_LEN_TOT];
gchar  *s_name_tot[SPARAM_LEN_TOT],          *v_name_tot[VPARAM_LEN_TOT];
gdouble s_lb_tot[SPARAM_LEN_TOT],             v_lb_tot[VPARAM_LEN_TOT];
gdouble s_ub_tot[SPARAM_LEN_TOT],             v_ub_tot[VPARAM_LEN_TOT];
gdouble s_scale_tot[SPARAM_LEN_TOT],          v_scale_tot[VPARAM_LEN_TOT];
gdouble s_abstol_tot[SPARAM_LEN_TOT],         v_abstol_tot[VPARAM_LEN_TOT];
gdouble s_defval_tot[SPARAM_LEN_TOT],         v_defval_tot[VPARAM_LEN_TOT];
NcmParamType s_ftype_tot[SPARAM_LEN_TOT], v_ftype_tot[VPARAM_LEN_TOT];
guint v_len_tot[VPARAM_LEN_TOT];

static gchar *name, *nick, **s_symbol, **v_symbol, **s_name, **v_name;
static gdouble *s_lb, *v_lb, *s_ub, *v_ub, *s_scale, *v_scale;
static gdouble *s_abstol, *v_abstol, *s_defval, *v_defval;
static NcmParamType *s_ftype, *v_ftype;
static guint *v_len, ci_sparam_len, ci_vparam_len;
static gchar *name_ext;

static void
ncm_model_test_base_class_init (NcmModelClass* model_class)
{
  guint i;

  for (i = 0; i < ci_sparam_len; i++)
  {
    const gdouble base_number = pow (10.0, g_test_rand_double_range (-10.0, 10.0));
    s_symbol[i] = g_strdup_printf ("SP%s_%u", name_ext, i);
    s_name[i]   = g_strdup_printf ("SP%s%u", name_ext, i);
    s_lb[i]     = g_test_rand_double_range (-base_number, base_number);
    s_ub[i]     = g_test_rand_double_range (     s_lb[i], base_number);
    s_scale[i]  = fabs (((s_ub[i] + s_lb[i]) != 0 ? (s_ub[i] + s_lb[i]) : 1.0) * pow (10.0, -g_test_rand_double_range (1.0,  2.0)));
    s_abstol[i] = s_scale[i] * pow (10.0, -g_test_rand_double_range (1.0, 11.0));
    s_defval[i] = g_test_rand_double_range ( s_lb[i], s_ub[i]);
    s_ftype[i]  = g_test_rand_bit () ? NCM_PARAM_TYPE_FREE  : NCM_PARAM_TYPE_FIXED;
  }

  for (i = 0; i < ci_vparam_len; i++)
  {
    const gdouble base_number = pow (10.0, g_test_rand_double_range (-10.0, 10.0));
    v_len[i]    = g_test_rand_int_range (1, 43);
    v_symbol[i] = g_strdup_printf ("VP%s_%u", name_ext, i);
    v_name[i]   = g_strdup_printf ("VP%s%u", name_ext, i);
    v_lb[i]     = g_test_rand_double_range (-base_number, base_number);
    v_ub[i]     = g_test_rand_double_range (     v_lb[i], base_number);
    v_scale[i]  = fabs (((v_ub[i] + v_lb[i]) != 0 ? (v_ub[i] + v_lb[i]) : 1.0)  * pow (10.0, -g_test_rand_double_range (1.0,  2.0)));
    v_abstol[i] = v_scale[i] * pow (10.0, -g_test_rand_double_range (1.0, 11.0));
    v_defval[i] = g_test_rand_double_range ( v_lb[i], v_ub[i]);
    v_ftype[i]  = g_test_rand_bit () ? NCM_PARAM_TYPE_FREE  : NCM_PARAM_TYPE_FIXED;
  }

  ncm_model_class_add_params (model_class, ci_sparam_len, ci_vparam_len, PROP_SIZE);

  g_assert (model_class->sparam_len == model_class->parent_sparam_len + ci_sparam_len);
  g_assert (model_class->vparam_len == model_class->parent_vparam_len + ci_vparam_len);

  ncm_model_class_set_name_nick (model_class, name, nick);

  g_assert (model_class->name != name);
  g_assert (model_class->nick != nick);
  g_assert_cmpstr (model_class->name, ==, name);
  g_assert_cmpstr (model_class->nick, ==, nick);

  NCM_TEST_FAIL (ncm_model_class_check_params_info (model_class));   

  for (i = 0; i < ci_sparam_len; i++)
  {
    ncm_model_class_set_sparam (model_class, model_class->parent_sparam_len + i, s_symbol[i], s_name[i],
                                s_lb[i], s_ub[i], s_scale[i], s_abstol[i],
                                s_defval[i], s_ftype[i]);
  }

  for (i = 0; i < ci_vparam_len; i++)
  {
    ncm_model_class_set_vparam (model_class, model_class->parent_vparam_len + i, v_len[i], v_symbol[i], v_name[i],
                                v_lb[i], v_ub[i], v_scale[i], v_abstol[i],
                                v_defval[i], v_ftype[i]);
  }

  ncm_model_class_check_params_info (model_class);

  for (i = 0; i < ci_sparam_len; i++)
  {
    NcmSParam *sp = g_ptr_array_index (model_class->sparam, model_class->parent_sparam_len + i);
    g_assert (sp->name   != s_name[i]);
    g_assert (sp->symbol != s_symbol[i]);
    g_assert_cmpstr (sp->name,   ==, s_name[i]);
    g_assert_cmpstr (sp->symbol, ==, s_symbol[i]);

    ncm_assert_cmpdouble (ncm_sparam_get_lower_bound (sp),        ==, s_lb[i]);
    ncm_assert_cmpdouble (ncm_sparam_get_upper_bound (sp),        ==, s_ub[i]);
    ncm_assert_cmpdouble (ncm_sparam_get_scale (sp),              ==, s_scale[i]);
    ncm_assert_cmpdouble (ncm_sparam_get_absolute_tolerance (sp), ==, s_abstol[i]);
    ncm_assert_cmpdouble (ncm_sparam_get_default_value (sp),      ==, s_defval[i]);
    g_assert_cmpint      (ncm_sparam_get_fit_type (sp),           ==, s_ftype[i]);
  }

  for (i = 0; i < ci_vparam_len; i++)
  {
    NcmVParam *vp = g_ptr_array_index (model_class->vparam, model_class->parent_vparam_len + i);
    g_assert (vp->default_sparam->name   != v_name[i]);
    g_assert (vp->default_sparam->symbol != v_symbol[i]);
    g_assert_cmpstr (vp->default_sparam->name,   ==, v_name[i]);
    g_assert_cmpstr (vp->default_sparam->symbol, ==, v_symbol[i]);

    ncm_assert_cmpdouble (ncm_vparam_get_lower_bound (vp, 0),        ==, v_lb[i]);
    ncm_assert_cmpdouble (ncm_vparam_get_upper_bound (vp, 0),        ==, v_ub[i]);
    ncm_assert_cmpdouble (ncm_vparam_get_scale (vp, 0),              ==, v_scale[i]);
    ncm_assert_cmpdouble (ncm_vparam_get_absolute_tolerance (vp, 0), ==, v_abstol[i]);
    ncm_assert_cmpdouble (ncm_vparam_get_default_value (vp, 0),      ==, v_defval[i]);
    g_assert_cmpint      (ncm_vparam_get_fit_type (vp, 0),           ==, v_ftype[i]);
  }
}

static void
ncm_model_test_class_init (NcmModelTestClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_ncm_model_test_set_property;
  model_class->get_property = &_ncm_model_test_get_property;

  ncm_model_class_add_params (model_class, 0, 0, 2);

  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("A",
                                                        NULL,
                                                        "Test A",
                                                        0.0, 100.0, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  object_class->finalize = ncm_model_test_finalize;

  name     = name_tot[0];
  nick     = nick_tot[0];
  s_symbol = &s_symbol_tot[0];
  v_symbol = &v_symbol_tot[0];
  s_name   = &s_name_tot[0];
  v_name   = &v_name_tot[0];
  s_lb     = &s_lb_tot[0];
  v_lb     = &v_lb_tot[0];
  s_ub     = &s_ub_tot[0];
  v_ub     = &v_ub_tot[0];
  s_scale  = &s_scale_tot[0];
  v_scale  = &v_scale_tot[0];
  s_abstol = &s_abstol_tot[0];
  v_abstol = &v_abstol_tot[0];
  s_defval = &s_defval_tot[0];
  v_defval = &v_defval_tot[0];
  s_ftype  = &s_ftype_tot[0];
  v_ftype  = &v_ftype_tot[0];
  v_len    = &v_len_tot[0];

  ci_sparam_len = SPARAM_LEN1;
  ci_vparam_len = VPARAM_LEN1;
  name_ext = "Base";

  ncm_model_test_base_class_init (model_class);
}

static void
ncm_model_test_child_class_init (NcmModelTestChildClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  model_class->set_property = &_ncm_model_test_child_set_property;
  model_class->get_property = &_ncm_model_test_child_get_property;

  ncm_model_class_add_params (model_class, 0, 0, 2);
  
  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("B",
                                                        NULL,
                                                        "Test B",
                                                        0.0, 100.0, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  object_class->finalize = ncm_model_test_child_finalize;

  name     = name_tot[1];
  nick     = nick_tot[1];
  s_symbol = &s_symbol_tot[SPARAM_LEN1];
  v_symbol = &v_symbol_tot[VPARAM_LEN1];
  s_name   = &s_name_tot[SPARAM_LEN1];
  v_name   = &v_name_tot[VPARAM_LEN1];
  s_lb     = &s_lb_tot[SPARAM_LEN1];
  v_lb     = &v_lb_tot[VPARAM_LEN1];
  s_ub     = &s_ub_tot[SPARAM_LEN1];
  v_ub     = &v_ub_tot[VPARAM_LEN1];
  s_scale  = &s_scale_tot[SPARAM_LEN1];
  v_scale  = &v_scale_tot[VPARAM_LEN1];
  s_abstol = &s_abstol_tot[SPARAM_LEN1];
  v_abstol = &v_abstol_tot[VPARAM_LEN1];
  s_defval = &s_defval_tot[SPARAM_LEN1];
  v_defval = &v_defval_tot[VPARAM_LEN1];
  s_ftype  = &s_ftype_tot[SPARAM_LEN1];
  v_ftype  = &v_ftype_tot[VPARAM_LEN1];
  v_len    = &v_len_tot[VPARAM_LEN1];

  ci_sparam_len = SPARAM_LEN2;
  ci_vparam_len = VPARAM_LEN2;
  name_ext = "Child";

  ncm_model_test_base_class_init (model_class);
}

static void
ncm_model_test_child_child_class_init (NcmModelTestChildChildClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);
  NcmModelClass* model_class = NCM_MODEL_CLASS (klass);

  ncm_model_class_add_params (model_class, 0, 0, 2);
  
  model_class->set_property = &_ncm_model_test_child_child_set_property;
  model_class->get_property = &_ncm_model_test_child_child_get_property;

  g_object_class_install_property (object_class,
                                   PROP_A,
                                   g_param_spec_double ("C",
                                                        NULL,
                                                        "Test C",
                                                        0.0, 100.0, 1.0,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  object_class->finalize = ncm_model_test_child_child_finalize;

  name     = name_tot[2];
  nick     = nick_tot[2];
  s_symbol = &s_symbol_tot[SPARAM_LEN1 + SPARAM_LEN2];
  v_symbol = &v_symbol_tot[VPARAM_LEN1 + VPARAM_LEN2];
  s_name   = &s_name_tot[SPARAM_LEN1 + SPARAM_LEN2];
  v_name   = &v_name_tot[VPARAM_LEN1 + VPARAM_LEN2];
  s_lb     = &s_lb_tot[SPARAM_LEN1 + SPARAM_LEN2];
  v_lb     = &v_lb_tot[VPARAM_LEN1 + VPARAM_LEN2];
  s_ub     = &s_ub_tot[SPARAM_LEN1 + SPARAM_LEN2];
  v_ub     = &v_ub_tot[VPARAM_LEN1 + VPARAM_LEN2];
  s_scale  = &s_scale_tot[SPARAM_LEN1 + SPARAM_LEN2];
  v_scale  = &v_scale_tot[VPARAM_LEN1 + VPARAM_LEN2];
  s_abstol = &s_abstol_tot[SPARAM_LEN1 + SPARAM_LEN2];
  v_abstol = &v_abstol_tot[VPARAM_LEN1 + VPARAM_LEN2];
  s_defval = &s_defval_tot[SPARAM_LEN1 + SPARAM_LEN2];
  v_defval = &v_defval_tot[VPARAM_LEN1 + VPARAM_LEN2];
  s_ftype  = &s_ftype_tot[SPARAM_LEN1 + SPARAM_LEN2];
  v_ftype  = &v_ftype_tot[VPARAM_LEN1 + VPARAM_LEN2];
  v_len    = &v_len_tot[VPARAM_LEN1 + VPARAM_LEN2];

  ci_sparam_len = SPARAM_LEN3;
  ci_vparam_len = VPARAM_LEN3;
  name_ext = "Grandchild";

  ncm_model_test_base_class_init (model_class);
}
