/***************************************************************************
 *            test_ncm_data_toys.c
 *
 *  Tue Sep 29 04:00:00 2026
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

void test_ncm_data_rosenbrock (void);
void test_ncm_data_funnel (void);
void test_ncm_data_gaussmix2d (void);
void test_ncm_data_gaussmix2d_no_model (void);
void test_ncm_data_gaussmix2d_no_model_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/ncm/data_toys/rosenbrock", &test_ncm_data_rosenbrock);
  g_test_add_func ("/ncm/data_toys/funnel", &test_ncm_data_funnel);
  g_test_add_func ("/ncm/data_toys/gaussmix2d", &test_ncm_data_gaussmix2d);
  g_test_add_func ("/ncm/data_toys/gaussmix2d/no_model", &test_ncm_data_gaussmix2d_no_model);
  g_test_add_func ("/ncm/data_toys/gaussmix2d/no_model/subprocess", &test_ncm_data_gaussmix2d_no_model_subprocess);

  g_test_run ();
}

void
test_ncm_data_rosenbrock (void)
{
  NcmModelRosenbrock *mrb = ncm_model_rosenbrock_new ();
  NcmMSet *mset           = ncm_mset_new (mrb, NULL, NULL);
  NcmData *data           = NCM_DATA (ncm_data_rosenbrock_new ());
  gdouble m2lnL;

  g_assert_true (ncm_data_is_init (data));
  g_assert_cmpuint (ncm_data_get_length (data), ==, 10);
  g_assert_cmpuint (ncm_data_get_dof (data), ==, 10);

  /* The minimum 0 at (1, 1); (1 - x_1)^2 / 10 = 0.1 at the origin. */
  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X1, 1.0);
  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X2, 1.0);
  ncm_data_m2lnL_val (data, mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 0.0);

  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X1, 0.0);
  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X2, 0.0);
  ncm_data_m2lnL_val (data, mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 0.1);

  ncm_data_free (data);
  ncm_mset_free (mset);
  ncm_model_rosenbrock_free (mrb);
}

void
test_ncm_data_funnel (void)
{
  NcmModelFunnel *mfu = ncm_model_funnel_new (4);
  NcmMSet *mset       = ncm_mset_new (mfu, NULL, NULL);
  NcmData *data       = NCM_DATA (ncm_data_funnel_new ());
  NcmModel *model     = NCM_MODEL (mfu);
  gdouble m2lnL;
  guint i;

  g_assert_true (ncm_data_is_init (data));

  /* At nu = 0 the x_i have unit variance: -2 ln L = sum x_i^2. */
  ncm_model_param_set (model, NCM_MODEL_FUNNEL_NU, 0.0);

  for (i = 0; i < 4; i++)
    ncm_model_orig_vparam_set (model, NCM_MODEL_FUNNEL_X, i, 0.0);

  ncm_data_m2lnL_val (data, mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 0.0);

  for (i = 0; i < 4; i++)
    ncm_model_orig_vparam_set (model, NCM_MODEL_FUNNEL_X, i, 1.0);

  ncm_data_m2lnL_val (data, mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 4.0);

  /* At x = 0 only n nu + (nu / 3)^2 remains: 4 * 3 + 1 at nu = 3. */
  for (i = 0; i < 4; i++)
    ncm_model_orig_vparam_set (model, NCM_MODEL_FUNNEL_X, i, 0.0);

  ncm_model_param_set (model, NCM_MODEL_FUNNEL_NU, 3.0);
  ncm_data_m2lnL_val (data, mset, &m2lnL);
  g_assert_cmpfloat (m2lnL, ==, 13.0);

  ncm_data_free (data);
  ncm_mset_free (mset);
  ncm_model_funnel_free (mfu);
}

void
test_ncm_data_gaussmix2d (void)
{
  NcmModelRosenbrock *mrb = ncm_model_rosenbrock_new ();
  NcmMSet *mset           = ncm_mset_new (mrb, NULL, NULL);
  NcmData *data           = NCM_DATA (ncm_data_gaussmix2d_new ());
  gdouble m2lnL_a, m2lnL_b;

  g_assert_true (ncm_data_is_init (data));

  /* Each component mean is a local maximum of the density: moving away increases
   * -2 ln L. */
  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X1, 1.5);
  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X2, 0.0);
  ncm_data_m2lnL_val (data, mset, &m2lnL_a);

  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X1, 1.6);
  ncm_data_m2lnL_val (data, mset, &m2lnL_b);
  g_assert_cmpfloat (m2lnL_b, >, m2lnL_a);

  /* Between the modes the density is lower than at either mean. */
  ncm_model_param_set (NCM_MODEL (mrb), NCM_MODEL_ROSENBROCK_X1, 0.0);
  ncm_data_m2lnL_val (data, mset, &m2lnL_b);
  g_assert_cmpfloat (m2lnL_b, >, m2lnL_a);
  g_assert_true (gsl_finite (m2lnL_b));

  ncm_data_free (data);
  ncm_mset_free (mset);
  ncm_model_rosenbrock_free (mrb);
}

void
test_ncm_data_gaussmix2d_no_model (void)
{
  g_test_trap_subprocess ("/ncm/data_toys/gaussmix2d/no_model/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the model set needs a NcmModelRosenbrock*");
}

void
test_ncm_data_gaussmix2d_no_model_subprocess (void)
{
  NcmMSet *mset = ncm_mset_empty_new ();
  NcmData *data = NCM_DATA (ncm_data_gaussmix2d_new ());

  ncm_data_prepare (data, mset);
}

