/***************************************************************************
 *            test_ncm_model_ctrl.c
 *
 *  Tue February 16 14:02:12 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2016 <vitenti@uel.br>
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

void test_ncm_model_peek_host (TestNcmModelCtrl *test, gconstpointer pdata);
void test_ncm_model_peek_host_lifecycle (void);
void test_ncm_model_peek_host_cross_host_traps (void);
void test_ncm_model_peek_host_cross_host_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/model_ctrl/model_update", TestNcmModelCtrl, NULL,
              &test_ncm_model_ctrl_new,
              &test_ncm_model_ctrl_model_update,
              &test_ncm_model_ctrl_free);

  g_test_add ("/ncm/model_ctrl/update", TestNcmModelCtrl, NULL,
              &test_ncm_model_ctrl_new,
              &test_ncm_model_ctrl_update,
              &test_ncm_model_ctrl_free);

  g_test_add ("/ncm/model_ctrl/submodel_update", TestNcmModelCtrl, NULL,
              &test_ncm_model_ctrl_new,
              &test_ncm_model_ctrl_submodel_update,
              &test_ncm_model_ctrl_free);

  g_test_add ("/ncm/model_ctrl/traps", TestNcmModelCtrl, NULL,
              &test_ncm_model_ctrl_new,
              &test_ncm_model_ctrl_traps,
              &test_ncm_model_ctrl_free);

  g_test_add ("/ncm/model_ctrl/invalid/submodel_last_update/subprocess", TestNcmModelCtrl, NULL,
              &test_ncm_model_ctrl_new,
              &test_ncm_model_ctrl_invalid_submodel_last_update,
              &test_ncm_model_ctrl_free);

  g_test_add ("/ncm/model/peek_host", TestNcmModelCtrl, NULL,
              &test_ncm_model_ctrl_new,
              &test_ncm_model_peek_host,
              &test_ncm_model_ctrl_free);

  g_test_add_func ("/ncm/model/peek_host/lifecycle", &test_ncm_model_peek_host_lifecycle);
  g_test_add_func ("/ncm/model/peek_host/cross_host/traps", &test_ncm_model_peek_host_cross_host_traps);
  g_test_add_func ("/ncm/model/peek_host/cross_host/subprocess", &test_ncm_model_peek_host_cross_host_subprocess);

  g_test_run ();
}

void
test_ncm_model_ctrl_new (TestNcmModelCtrl *test, gconstpointer pdata)
{
  test->ctrl      = ncm_model_ctrl_new (NULL);
  test->submodel1 = NCM_MODEL (nc_hiprim_power_law_new ());
  test->submodel2 = NCM_MODEL (nc_hireion_camb_new ());

  /* NcHICosmo declares construction-only typed slots for prim/reion, so
   * submodels must be attached at construction time -- not via
   * ncm_model_add_submodel() afterward. */
  test->model = NCM_MODEL (nc_hicosmo_lcdm_new_full (NC_HIREION (test->submodel2),
                                                     NC_HIPRIM (test->submodel1), NULL));

  g_assert_true (test->ctrl != NULL);
  g_assert_true (NCM_IS_MODEL_CTRL (test->ctrl));

  g_assert_true (test->model != NULL);
  g_assert_true (NCM_IS_MODEL (test->model));

  g_assert_true (test->submodel1 != NULL);
  g_assert_true (NCM_IS_MODEL (test->submodel1));

  g_assert_true (test->submodel2 != NULL);
  g_assert_true (NCM_IS_MODEL (test->submodel2));
}

void
_set_destroyed (gpointer b)
{
  gboolean *destroyed = b;

  *destroyed = TRUE;
}

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
  /* Submodels are attached at construction now (see test_ncm_model_ctrl_new),
   * so the very first set_model() call already needs to pick up both prim
   * and reion in one shot -- there is no legal way to attach or replace a
   * submodel on an already-constructed NcHICosmo anymore, so what remains
   * testable here is: (1) first-attach detection covering the host and its
   * pre-existing submodels together, (2) idempotency, and (3) detecting a
   * genuinely different host object (e.g. swapping between two complete
   * cosmologies, which is still a legitimate operation). */
  g_assert_true (ncm_model_ctrl_set_model (test->ctrl, test->model));
  g_assert_true (!ncm_model_ctrl_set_model (test->ctrl, test->model));

  {
    NcmModel *other_model = NCM_MODEL (nc_hicosmo_lcdm_new ());

    g_assert_true (ncm_model_ctrl_set_model (test->ctrl, other_model));
    g_assert_true (!ncm_model_ctrl_set_model (test->ctrl, other_model));

    NCM_TEST_FREE (ncm_model_free, other_model);
  }

  /* Restore tracking to the fixture's own model -- callers of this helper
   * (e.g. test_ncm_model_ctrl_update()) assume ctrl is left watching
   * test->model afterward. */
  g_assert_true (ncm_model_ctrl_set_model (test->ctrl, test->model));
}

void
test_ncm_model_ctrl_update (TestNcmModelCtrl *test, gconstpointer pdata)
{
  g_assert_true (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert_true (!ncm_model_ctrl_update (test->ctrl, test->model));

  ncm_model_orig_param_set (test->model, 0,
                            ncm_model_orig_param_get (test->model, 0) * 0.999);
  g_assert_true (ncm_model_ctrl_update (test->ctrl, test->model));

  g_assert_true (ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert_true (ncm_model_ctrl_model_last_update (test->ctrl));

  g_assert_true (!ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert_true (!ncm_model_ctrl_model_last_update (test->ctrl));

  /* Now testing with submodels added */
  ncm_model_ctrl_force_update (test->ctrl);
  test_ncm_model_ctrl_model_update (test, pdata);

  g_assert_true (!ncm_model_ctrl_update (test->ctrl, test->model));

  ncm_model_orig_param_set (test->model, 0,
                            ncm_model_orig_param_get (test->model, 0) * 0.999);
  g_assert_true (ncm_model_ctrl_update (test->ctrl, test->model));

  g_assert_true (ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert_true (ncm_model_ctrl_model_last_update (test->ctrl));

  g_assert_true (!ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert_true (!ncm_model_ctrl_model_last_update (test->ctrl));
}

void
test_ncm_model_ctrl_submodel_update (TestNcmModelCtrl *test, gconstpointer pdata)
{
  /* Submodels are attached at construction now (see test_ncm_model_ctrl_new),
   * so the very first update() call already picks up the host together with
   * both prim and reion -- there is no legal way to attach a submodel to an
   * already-constructed NcHICosmo anymore. What remains fully testable:
   * has_submodel() reflecting the construction-time attachment, and each
   * submodel's own param changes being detected independently of the host's
   * and of each other.
   *
   * Note: on this very first call, ctrl_model != model triggers
   * ncm_model_ctrl_update()'s internal ncm_model_ctrl_set_model() call,
   * which itself creates and fully syncs the per-submodel sub-controllers
   * (pre-existing behaviour of ncm_model_ctrl_update(), not specific to
   * construction-time attachment) -- so submodel_last_update() is FALSE
   * right after this first call, even though has_submodel() is already
   * TRUE. */
  g_assert_true (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert_true (ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert_true (ncm_model_ctrl_model_has_submodel (test->ctrl, nc_hiprim_id ()));
  g_assert_true (ncm_model_ctrl_model_has_submodel (test->ctrl, nc_hireion_id ()));
  g_assert_true (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert_true (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  g_assert_true (!ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert_true (!ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert_true (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert_true (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  ncm_model_orig_param_set (test->model, 0,
                            ncm_model_orig_param_get (test->model, 0) * 0.999);
  g_assert_true (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert_true (ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert_true (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert_true (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  ncm_model_orig_param_set (test->submodel1, 0,
                            ncm_model_orig_param_get (test->submodel1, 0) * 0.999);
  g_assert_true (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert_true (!ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert_true (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert_true (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  ncm_model_orig_param_set (test->submodel2, 0,
                            ncm_model_orig_param_get (test->submodel2, 0) * 0.999);
  g_assert_true (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert_true (!ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert_true (!ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert_true (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  ncm_model_orig_param_set (test->submodel1, 0,
                            ncm_model_orig_param_get (test->submodel1, 0) * 0.999);
  ncm_model_orig_param_set (test->submodel2, 0,
                            ncm_model_orig_param_get (test->submodel2, 0) * 0.999);
  g_assert_true (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert_true (!ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert_true (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert_true (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));

  ncm_model_orig_param_set (test->model, 0,
                            ncm_model_orig_param_get (test->model, 0) * 0.999);
  ncm_model_orig_param_set (test->submodel1, 0,
                            ncm_model_orig_param_get (test->submodel1, 0) * 0.999);
  ncm_model_orig_param_set (test->submodel2, 0,
                            ncm_model_orig_param_get (test->submodel2, 0) * 0.999);
  g_assert_true (ncm_model_ctrl_update (test->ctrl, test->model));
  g_assert_true (ncm_model_ctrl_model_last_update (test->ctrl));
  g_assert_true (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hiprim_id ()));
  g_assert_true (ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ()));
}

void
test_ncm_model_ctrl_traps (TestNcmModelCtrl *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/model_ctrl/invalid/submodel_last_update/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_model_ctrl_invalid_submodel_last_update (TestNcmModelCtrl *test, gconstpointer pdata)
{
  ncm_model_ctrl_submodel_last_update (test->ctrl, nc_hireion_id ());
}

void
test_ncm_model_peek_host (TestNcmModelCtrl *test, gconstpointer pdata)
{
  /* Both submodels were attached to test->model at construction time (see
   * test_ncm_model_ctrl_new()) -- their host backpointer must reflect that,
   * while the host itself (not anyone's submodel) has no host of its own. */
  g_assert_true (ncm_model_peek_host (test->submodel1) == test->model);
  g_assert_true (ncm_model_peek_host (test->submodel2) == test->model);
  g_assert_true (ncm_model_peek_host (test->model) == NULL);
}

void
test_ncm_model_peek_host_lifecycle (void)
{
  NcHIReion *reion     = NC_HIREION (nc_hireion_camb_new ());
  NcHICosmoLCDM *cosmo = nc_hicosmo_lcdm_new_full (reion, NULL, NULL);

  g_assert_true (ncm_model_peek_host (NCM_MODEL (reion)) == NCM_MODEL (cosmo));

  /* reion keeps its own external ref (from nc_hireion_camb_new() above), so
   * it survives cosmo's disposal -- its host backpointer is a weak ref and
   * must correctly go stale (peek_host returns NULL), not dangle. */
  nc_hicosmo_free (NC_HICOSMO (cosmo));

  g_assert_true (ncm_model_peek_host (NCM_MODEL (reion)) == NULL);

  nc_hireion_free (reion);
}

void
test_ncm_model_peek_host_cross_host_subprocess (void)
{
  NcHIReion *reion   = NC_HIREION (nc_hireion_camb_new ());
  NcHICosmo *cosmo_a = NC_HICOSMO (nc_hicosmo_lcdm_new_full (reion, NULL, NULL));
  NcHICosmo *cosmo_b;

  /* cosmo_a is still alive (held above), so reion's host backpointer is
   * live and points elsewhere -- attaching the same slotted-type submodel
   * instance to a second host must be rejected, not silently reassigned. */
  cosmo_b = NC_HICOSMO (nc_hicosmo_lcdm_new_full (reion, NULL, NULL));

  nc_hicosmo_free (cosmo_a);
  nc_hicosmo_free (cosmo_b);
  nc_hireion_free (reion);
}

void
test_ncm_model_peek_host_cross_host_traps (void)
{
  g_test_trap_subprocess ("/ncm/model/peek_host/cross_host/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

