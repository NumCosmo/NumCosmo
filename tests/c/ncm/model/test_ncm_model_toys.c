/***************************************************************************
 *            test_ncm_model_toys.c
 *
 *  Mon Sep 28 20:00:00 2026
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

#include <glib.h>
#include <glib-object.h>

void test_ncm_model_funnel (void);
void test_ncm_model_rosenbrock (void);
void test_ncm_model_mvnd (void);
void test_ncm_model_mvnd_errors (void);
void test_ncm_model_mvnd_dim_mismatch_subprocess (void);
void test_ncm_model_mvnd_mean_len_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/model/funnel", &test_ncm_model_funnel);
  g_test_add_func ("/ncm/model/rosenbrock", &test_ncm_model_rosenbrock);
  g_test_add_func ("/ncm/model/mvnd", &test_ncm_model_mvnd);
  g_test_add_func ("/ncm/model/mvnd/errors", &test_ncm_model_mvnd_errors);
  g_test_add_func ("/ncm/model/mvnd/errors/dim_mismatch/subprocess", &test_ncm_model_mvnd_dim_mismatch_subprocess);
  g_test_add_func ("/ncm/model/mvnd/errors/mean_len/subprocess", &test_ncm_model_mvnd_mean_len_subprocess);

  g_test_run ();
}

void
test_ncm_model_funnel (void)
{
  NcmModelFunnel *mfu = ncm_model_funnel_new (4);
  NcmModel *model     = NCM_MODEL (mfu);
  NcmModelFunnel *mfu_ref;

  g_assert_cmpstr (ncm_model_name (model), ==, "Funnel distribution");
  g_assert_cmpstr (ncm_model_nick (model), ==, "Funnel");
  g_assert_cmpstr (ncm_mset_get_ns_by_id (ncm_model_funnel_id ()), ==, "NcmModelFunnel");

  g_assert_cmpuint (ncm_model_len (model), ==, 5);
  g_assert_cmpuint (ncm_model_vparam_len (model, NCM_MODEL_FUNNEL_X), ==, 4);
  g_assert_cmpstr (ncm_model_param_name (model, NCM_MODEL_FUNNEL_NU), ==, "nu");
  g_assert_cmpstr (ncm_model_param_name (model, ncm_model_vparam_index (model, NCM_MODEL_FUNNEL_X, 3)), ==, "x_3");

  mfu_ref = ncm_model_funnel_ref (mfu);
  g_assert_true (mfu_ref == mfu);
  ncm_model_funnel_free (mfu_ref);

  ncm_model_funnel_clear (&mfu);
  g_assert_null (mfu);
}

void
test_ncm_model_rosenbrock (void)
{
  NcmModelRosenbrock *mrb = ncm_model_rosenbrock_new ();
  NcmModel *model         = NCM_MODEL (mrb);
  NcmModelRosenbrock *mrb_ref;

  g_assert_cmpstr (ncm_model_name (model), ==, "Rosenbrock distribution");
  g_assert_cmpstr (ncm_model_nick (model), ==, "Rosenbrock");
  g_assert_cmpstr (ncm_mset_get_ns_by_id (ncm_model_rosenbrock_id ()), ==, "NcmModelRosenbrock");

  g_assert_cmpuint (ncm_model_len (model), ==, 2);
  g_assert_cmpstr (ncm_model_param_name (model, NCM_MODEL_ROSENBROCK_X1), ==, "x1");
  g_assert_cmpstr (ncm_model_param_name (model, NCM_MODEL_ROSENBROCK_X2), ==, "x2");

  mrb_ref = ncm_model_rosenbrock_ref (mrb);
  g_assert_true (mrb_ref == mrb);
  ncm_model_rosenbrock_free (mrb_ref);

  ncm_model_rosenbrock_clear (&mrb);
  g_assert_null (mrb);
}

void
test_ncm_model_mvnd (void)
{
  NcmModelMVND *mvnd = ncm_model_mvnd_new (3);
  NcmModel *model    = NCM_MODEL (mvnd);
  NcmVector *y       = ncm_vector_new (3);
  NcmModelMVND *mvnd_ref;
  guint dim, i;

  g_assert_cmpstr (ncm_model_name (model), ==, "Multivariate normal mean");
  g_assert_cmpstr (ncm_model_nick (model), ==, "MVND");
  g_assert_cmpstr (ncm_mset_get_ns_by_id (ncm_model_mvnd_id ()), ==, "NcmModelMVND");

  g_object_get (mvnd, "dim", &dim, NULL);
  g_assert_cmpuint (dim, ==, 3);
  g_assert_cmpuint (ncm_model_vparam_len (model, NCM_MODEL_MVND_MEAN), ==, 3);

  for (i = 0; i < 3; i++)
    ncm_model_orig_vparam_set (model, NCM_MODEL_MVND_MEAN, i, 0.5 * (i + 1.0));

  ncm_model_mvnd_mean (mvnd, y);

  for (i = 0; i < 3; i++)
    g_assert_cmpfloat (ncm_vector_get (y, i), ==, 0.5 * (i + 1.0));

  /* The mean follows later changes of the parameters. */
  ncm_model_orig_vparam_set (model, NCM_MODEL_MVND_MEAN, 1, -2.0);
  ncm_model_mvnd_mean (mvnd, y);
  g_assert_cmpfloat (ncm_vector_get (y, 1), ==, -2.0);

  mvnd_ref = ncm_model_mvnd_ref (mvnd);
  g_assert_true (mvnd_ref == mvnd);
  ncm_model_mvnd_free (mvnd_ref);

  ncm_vector_free (y);
  ncm_model_mvnd_clear (&mvnd);
  g_assert_null (mvnd);
}

void
test_ncm_model_mvnd_errors (void)
{
  g_test_trap_subprocess ("/ncm/model/mvnd/errors/dim_mismatch/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*dimension 3 differs from the mean length 1; set both `dim' and `mu-length', "
                             "or use ncm_model_mvnd_new()*");

  g_test_trap_subprocess ("/ncm/model/mvnd/errors/mean_len/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the mean has 3 components, but the vector has 2*");
}

void
test_ncm_model_mvnd_dim_mismatch_subprocess (void)
{
  g_object_unref (g_object_new (NCM_TYPE_MODEL_MVND, "dim", 3, NULL));
}

void
test_ncm_model_mvnd_mean_len_subprocess (void)
{
  NcmModelMVND *mvnd = ncm_model_mvnd_new (3);
  NcmVector *y       = ncm_vector_new (2);

  ncm_model_mvnd_mean (mvnd, y);

  ncm_vector_free (y);
  ncm_model_mvnd_free (mvnd);
}

