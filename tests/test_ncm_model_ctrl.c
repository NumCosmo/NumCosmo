/***************************************************************************
 *            test_ncm_model_ctrl.c
 *
 *  Tue February 16 14:02:12 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2016 <sandro@isoftware.com.br>
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

typedef struct _TestNcmModelCtrl
{
  NcmModelCtrl *ctrl;
  NcmModel *model;
  NcmModel *submodel1;
  NcmModel *submodel2;
} TestNcmModelCtrl;

void test_ncm_model_ctrl_new (TestNcmModelCtrl *test, gconstpointer pdata);
void test_ncm_model_ctrl_free (TestNcmModelCtrl *test, gconstpointer pdata);

void test_ncm_model_ctrl_model_update (TestNcmModelCtrl *test, gconstpointer pdata);
void test_ncm_model_ctrl_update (TestNcmModelCtrl *test, gconstpointer pdata);
void test_ncm_model_ctrl_submodel_update (TestNcmModelCtrl *test, gconstpointer pdata);

void test_ncm_model_ctrl_traps (TestNcmModelCtrl *test, gconstpointer pdata);
void test_ncm_model_ctrl_invalid_submodel_last_update (TestNcmModelCtrl *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();
  
  g_test_add ("/numcosmo/ncm_model_ctrl/model_update", TestNcmModelCtrl, NULL, 
              &test_ncm_model_ctrl_new, 
              &test_ncm_model_ctrl_model_update, 
              &test_ncm_model_ctrl_free);

  g_test_add ("/numcosmo/ncm_model_ctrl/update", TestNcmModelCtrl, NULL, 
              &test_ncm_model_ctrl_new, 
              &test_ncm_model_ctrl_update, 
              &test_ncm_model_ctrl_free);

  g_test_add ("/numcosmo/ncm_model_ctrl/submodel_update", TestNcmModelCtrl, NULL, 
              &test_ncm_model_ctrl_new, 
              &test_ncm_model_ctrl_submodel_update, 
              &test_ncm_model_ctrl_free);

  g_test_add ("/numcosmo/ncm_model_ctrl/traps", TestNcmModelCtrl, NULL,
              &test_ncm_model_ctrl_new,
              &test_ncm_model_ctrl_traps,
              &test_ncm_model_ctrl_free);
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 38))
  g_test_add ("/numcosmo/ncm_model_ctrl/invalid/submodel_last_update/subprocess", TestNcmModelCtrl, NULL, 
              &test_ncm_model_ctrl_new, 
              &test_ncm_model_ctrl_invalid_submodel_last_update, 
              &test_ncm_model_ctrl_free);
#endif
  g_test_run ();
}

void
test_ncm_model_ctrl_new (TestNcmModelCtrl *test, gconstpointer pdata)
{
  test->ctrl      = ncm_model_ctrl_new (NULL);
  test->model     = NCM_MODEL (nc_hicosmo_lcdm_new ());
  test->submodel1 = NCM_MODEL (nc_hiprim_power_law_new ());
  test->submodel2 = NCM_MODEL (nc_hireion_camb_new ());

  g_assert (test->ctrl != NULL);
  g_assert (NCM_IS_MODEL_CTRL (test->ctrl));

  g_assert (test->model != NULL);
  g_assert (NCM_IS_MODEL (test->model));

  g_assert (test->submodel1 != NULL);
  g_assert (NCM_IS_MODEL (test->submodel1));

  g_assert (test->submodel2 != NULL);
  g_assert (NCM_IS_MODEL (test->submodel2));

}

void _set_destroyed (gpointer b) { gboolean *destroyed = b; *destroyed = TRUE; }

void
test_ncm_model_ctrl_free (TestNcmModelCtrl *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_model_ctrl_free, test->ctrl);
  NCM_TEST_FREE (ncm_model_free, test->model);
  NCM_TEST_FREE (ncm_model_free, test->submodel1);
  NCM_TEST_FREE (ncm_model_free, test->submodel2);
}

void
test_ncm_model_ctrl_model_update (TestNcmModelCtrl *test, gconstpointer pdata)
{
  g_assert (ncm_model_ctrl_set_model (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_set_model (test->ctrl, test->model));

  ncm_model_add_submodel (test->model, test->submodel1);
  g_assert (ncm_model_ctrl_set_model (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_set_model (test->ctrl, test->model));

  ncm_model_add_submodel (test->model, test->submodel2);
  g_assert (ncm_model_ctrl_set_model (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_set_model (test->ctrl, test->model));

  {
    NcmModel *submodel2 = NCM_MODEL (nc_hireion_camb_new ());

    ncm_model_add_submodel (test->model, submodel2);
    g_assert (ncm_model_ctrl_set_model (test->ctrl, test->model));
    g_assert (!ncm_model_ctrl_set_model (test->ctrl, test->model));

    ncm_model_add_submodel (test->model, test->submodel2);
    g_assert (ncm_model_ctrl_set_model (test->ctrl, test->model));
    g_assert (!ncm_model_ctrl_set_model (test->ctrl, test->model));

    NCM_TEST_FREE (ncm_model_free, submodel2);
  }
  
}

void
test_ncm_model_ctrl_update (TestNcmModelCtrl *test, gconstpointer pdata)
{
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_update (test->ctrl, test->model));
  
  ncm_model_orig_param_set (test->model, 0, 
                            ncm_model_orig_param_get (test->model, 0) * 0.999);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));

  g_assert (ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (ncm_model_ctrl_model_last_update (test->ctrl));

  g_assert (!ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_model_last_update (test->ctrl));

  /* Now testing with submodels added */
  ncm_model_ctrl_force_update (test->ctrl);
  test_ncm_model_ctrl_model_update (test, pdata);

  g_assert (!ncm_model_ctrl_update (test->ctrl, test->model));
  
  ncm_model_orig_param_set (test->model, 0, 
                            ncm_model_orig_param_get (test->model, 0) * 0.999);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));

  g_assert (ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (ncm_model_ctrl_model_last_update (test->ctrl));

  g_assert (!ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_model_last_update (test->ctrl));  
}

void
test_ncm_model_ctrl_submodel_update (TestNcmModelCtrl *test, gconstpointer pdata)
{
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (!ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_model_last_update (test->ctrl));

  ncm_model_orig_param_set (test->model, 0, 
                            ncm_model_orig_param_get (test->model, 0) * 0.999);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (ncm_model_ctrl_model_last_update (test->ctrl));

  g_assert (!ncm_model_ctrl_model_has_submodel (test->ctrl, nc_hiprim_id ()));
  g_assert (!ncm_model_ctrl_model_has_submodel (test->ctrl, nc_hireion_id ()));

  ncm_model_add_submodel (test->model, test->submodel1);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (ncm_model_ctrl_model_has_submodel (test->ctrl, nc_hiprim_id ()));
  g_assert (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));

  g_assert (!ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));

  ncm_model_orig_param_set (test->model, 0, 
                            ncm_model_orig_param_get (test->model, 0) * 0.999);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));

  ncm_model_orig_param_set (test->submodel1, 0, 
                            ncm_model_orig_param_get (test->submodel1, 0) * 0.999);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));

  ncm_model_add_submodel (test->model, test->submodel2);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (ncm_model_ctrl_model_has_submodel (test->ctrl, nc_hireion_id ()));
  g_assert (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  ncm_model_orig_param_set (test->model, 0, 
                            ncm_model_orig_param_get (test->model, 0) * 0.999);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  ncm_model_orig_param_set (test->submodel1, 0, 
                            ncm_model_orig_param_get (test->submodel1, 0) * 0.999);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  ncm_model_orig_param_set (test->submodel2, 0, 
                            ncm_model_orig_param_get (test->submodel2, 0) * 0.999);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  ncm_model_orig_param_set (test->submodel1, 0, 
                            ncm_model_orig_param_get (test->submodel1, 0) * 0.999);
  ncm_model_orig_param_set (test->submodel2, 0, 
                            ncm_model_orig_param_get (test->submodel2, 0) * 0.999);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (!ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  ncm_model_orig_param_set (test->model, 0, 
                            ncm_model_orig_param_get (test->model, 0) * 0.999);
  ncm_model_orig_param_set (test->submodel1, 0, 
                            ncm_model_orig_param_get (test->submodel1, 0) * 0.999);
  ncm_model_orig_param_set (test->submodel2, 0, 
                            ncm_model_orig_param_get (test->submodel2, 0) * 0.999);
  g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert (ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  {
    NcmModel *submodel2 = NCM_MODEL (nc_hireion_camb_new ());

    ncm_model_add_submodel (test->model, submodel2);

    g_assert (ncm_model_ctrl_update (test->ctrl, test->model));
    g_assert (!ncm_model_ctrl_model_last_update (test->ctrl));
    g_assert (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
    g_assert (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

    ncm_model_add_submodel (test->model, test->submodel2);
    g_assert (ncm_model_ctrl_set_model (test->ctrl, test->model));
    g_assert (!ncm_model_ctrl_set_model (test->ctrl, test->model));

    NCM_TEST_FREE (ncm_model_free, submodel2);
  }
}

void
test_ncm_model_ctrl_traps (TestNcmModelCtrl *test, gconstpointer pdata)
{
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 38))
  g_test_trap_subprocess ("/numcosmo/ncm_model_ctrl/invalid/submodel_last_update/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_model_ctrl_invalid_submodel_last_update (TestNcmModelCtrl *test, gconstpointer pdata)
{
  ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ());
}
