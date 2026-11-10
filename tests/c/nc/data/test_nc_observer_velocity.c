/***************************************************************************
 *            test_nc_observer_velocity.c
 *
 *  Sun Oct 18 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_observer_velocity.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

typedef struct _TestNcObserverVelocity
{
  NcObserverVelocity *ov;
} TestNcObserverVelocity;

void test_nc_observer_velocity_new (TestNcObserverVelocity *test, gconstpointer pdata);
void test_nc_observer_velocity_free (TestNcObserverVelocity *test, gconstpointer pdata);

void test_nc_observer_velocity_params (TestNcObserverVelocity *test, gconstpointer pdata);
void test_nc_observer_velocity_doppler (TestNcObserverVelocity *test, gconstpointer pdata);
void test_nc_observer_velocity_serialize (TestNcObserverVelocity *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/nc/observer_velocity/params", TestNcObserverVelocity, NULL,
              &test_nc_observer_velocity_new,
              &test_nc_observer_velocity_params,
              &test_nc_observer_velocity_free);

  g_test_add ("/nc/observer_velocity/doppler", TestNcObserverVelocity, NULL,
              &test_nc_observer_velocity_new,
              &test_nc_observer_velocity_doppler,
              &test_nc_observer_velocity_free);

  g_test_add ("/nc/observer_velocity/serialize", TestNcObserverVelocity, NULL,
              &test_nc_observer_velocity_new,
              &test_nc_observer_velocity_serialize,
              &test_nc_observer_velocity_free);

  g_test_run ();
}

void
test_nc_observer_velocity_new (TestNcObserverVelocity *test, gconstpointer pdata)
{
  test->ov = nc_observer_velocity_new ();

  g_assert_true (NC_IS_OBSERVER_VELOCITY (test->ov));
}

void
test_nc_observer_velocity_free (TestNcObserverVelocity *test, gconstpointer pdata)
{
  NCM_TEST_FREE (nc_observer_velocity_free, test->ov);
}

void
test_nc_observer_velocity_params (TestNcObserverVelocity *test, gconstpointer pdata)
{
  NcmModel *model = NCM_MODEL (test->ov);
  guint i;

  g_assert_cmpuint (ncm_model_len (model), ==, NC_OBSERVER_VELOCITY_SPARAM_LEN);
  g_assert_cmpuint (ncm_model_vparam_array_len (model), ==, 0);

  for (i = 0; i < NC_OBSERVER_VELOCITY_SPARAM_LEN; i++)
  {
    ncm_assert_cmpdouble (ncm_model_orig_param_get (model, i), ==, NC_OBSERVER_VELOCITY_DEFAULT_V);
    ncm_assert_cmpdouble (ncm_model_param_get_lower_bound (model, i), ==, -NC_OBSERVER_VELOCITY_MAX_V);
    ncm_assert_cmpdouble (ncm_model_param_get_upper_bound (model, i), ==, NC_OBSERVER_VELOCITY_MAX_V);
    g_assert_cmpint (ncm_model_param_get_ftype (model, i), ==, NCM_PARAM_TYPE_FIXED);
  }

  /* At rest: unit Doppler factor and no change of redshift in any direction */
  ncm_assert_cmpdouble (nc_observer_velocity_speed (test->ov), ==, 0.0);
  ncm_assert_cmpdouble (nc_observer_velocity_doppler (test->ov, 0.3, -0.4, 0.866), ==, 1.0);
  ncm_assert_cmpdouble_e (nc_observer_velocity_z_rest (test->ov, 0.05, 0.3, -0.4, 0.866), ==, 0.05, 1.0e-15, 0.0);
}

void
test_nc_observer_velocity_doppler (TestNcObserverVelocity *test, gconstpointer pdata)
{
  NcmModel *model      = NCM_MODEL (test->ov);
  const gdouble v      = 600.0;
  const gdouble beta   = v / (1.0e-3 * ncm_c_c ());
  const gdouble gamma1 = sqrt (1.0 - beta * beta);
  const gdouble z      = 0.03;

  ncm_model_orig_param_set (model, NC_OBSERVER_VELOCITY_VX, 0.0);
  ncm_model_orig_param_set (model, NC_OBSERVER_VELOCITY_VY, v);
  ncm_model_orig_param_set (model, NC_OBSERVER_VELOCITY_VZ, 0.0);

  ncm_assert_cmpdouble_e (nc_observer_velocity_speed (test->ov), ==, v, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_observer_velocity_beta_n (test->ov, 0.0, 1.0, 0.0), ==, beta, 1.0e-15, 0.0);

  /* Towards the motion: blueshift, 1 + z_rest larger than 1 + z_obs by sqrt((1 + beta) / (1 - beta)) */
  ncm_assert_cmpdouble_e (nc_observer_velocity_doppler (test->ov, 0.0, 1.0, 0.0), ==, sqrt ((1.0 + beta) / (1.0 - beta)), 1.0e-15, 0.0);
  /* Away from it: the inverse; across it: the transverse factor */
  ncm_assert_cmpdouble_e (nc_observer_velocity_doppler (test->ov, 0.0, -1.0, 0.0), ==, sqrt ((1.0 - beta) / (1.0 + beta)), 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_observer_velocity_doppler (test->ov, 1.0, 0.0, 0.0), ==, gamma1, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_observer_velocity_doppler (test->ov, 0.0, 0.0, 1.0), ==, gamma1, 1.0e-15, 0.0);

  /* The rest-frame redshift and its first order in beta */
  ncm_assert_cmpdouble_e (nc_observer_velocity_z_rest (test->ov, z, 0.0, 1.0, 0.0), ==,
                          (1.0 + z) * sqrt ((1.0 + beta) / (1.0 - beta)) - 1.0, 1.0e-15, 0.0);
  ncm_assert_cmpdouble_e (nc_observer_velocity_z_rest (test->ov, z, 0.0, 1.0, 0.0), ==,
                          z + (1.0 + z) * beta, 0.0, 2.0 * beta * beta);
  ncm_assert_cmpdouble_e (nc_observer_velocity_z_rest (test->ov, z, 0.0, -1.0, 0.0), ==,
                          z - (1.0 + z) * beta, 0.0, 2.0 * beta * beta);

  /* A general direction: beta . n from the components */
  {
    const gdouble n[3] = {0.48, 0.6, 0.64};
    const gdouble bn   = beta * n[1];

    ncm_assert_cmpdouble_e (nc_observer_velocity_beta_n (test->ov, n[0], n[1], n[2]), ==, bn, 1.0e-15, 0.0);
    ncm_assert_cmpdouble_e (nc_observer_velocity_doppler (test->ov, n[0], n[1], n[2]), ==, gamma1 / (1.0 - bn), 1.0e-15, 0.0);
  }
}

void
test_nc_observer_velocity_serialize (TestNcObserverVelocity *test, gconstpointer pdata)
{
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  NcmModel *model   = NCM_MODEL (test->ov);
  NcObserverVelocity *dup;

  ncm_model_orig_param_set (model, NC_OBSERVER_VELOCITY_VX, -120.0);
  ncm_model_orig_param_set (model, NC_OBSERVER_VELOCITY_VY, 30.0);
  ncm_model_orig_param_set (model, NC_OBSERVER_VELOCITY_VZ, 250.0);

  dup = NC_OBSERVER_VELOCITY (ncm_serialize_dup_obj (ser, G_OBJECT (test->ov)));

  ncm_assert_cmpdouble (ncm_model_orig_param_get (NCM_MODEL (dup), NC_OBSERVER_VELOCITY_VX), ==, -120.0);
  ncm_assert_cmpdouble (ncm_model_orig_param_get (NCM_MODEL (dup), NC_OBSERVER_VELOCITY_VY), ==, 30.0);
  ncm_assert_cmpdouble (ncm_model_orig_param_get (NCM_MODEL (dup), NC_OBSERVER_VELOCITY_VZ), ==, 250.0);
  ncm_assert_cmpdouble (nc_observer_velocity_speed (dup), ==, nc_observer_velocity_speed (test->ov));

  nc_observer_velocity_free (dup);
  ncm_serialize_free (ser);
}

