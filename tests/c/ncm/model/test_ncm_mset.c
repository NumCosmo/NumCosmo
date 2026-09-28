/***************************************************************************
 *            test_ncm_mset.c
 *
 *  Wed May 13 15:19:36 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2015 <vitenti@uel.br>
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
#include <numcosmo/ncm/model/ncm_reparam_linear.h>

#include <math.h>
#include <string.h>
#include <glib.h>
#include <glib-object.h>
#include <glib/gstdio.h>

typedef struct _TestNcmMSet
{
  NcmMSet *mset;
  GPtrArray *ma;
  GArray *ma_destroyed;
} TestNcmMSet;

void test_ncm_mset_new (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_free (TestNcmMSet *test, gconstpointer pdata);

void test_ncm_mset_setpeek (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_setpospeek (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_pushpeek (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_fparams (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_dup (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_shallow_copy (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_saveload (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_saveload_submodel_reparam (void);
void test_ncm_mset_unslotted_submodel_attach (void);
void test_ncm_mset_unslotted_submodel_attach_subprocess (void);
void test_ncm_mset_load_unslotted_submodel_group (void);
void test_ncm_mset_two_level_submodel_slots (void);
void test_ncm_mset_peek_by_name_unknown (void);
void test_ncm_mset_remove_host (void);
void test_ncm_mset_split_full_name (void);
void test_ncm_mset_set_fmap_updates_models (void);
void test_ncm_mset_set_fmap_invalid (void);
void test_ncm_mset_stack_bound (void);
void test_ncm_mset_id_lookups (void);
void test_ncm_mset_ns_by_negative_id (void);
void test_ncm_mset_set_fmap_update_models (void);
void test_ncm_mset_max_model_nick (void);
void test_ncm_mset_params_pretty_print (void);
void test_ncm_mset_fparam_lookups (void);
void test_ncm_mset_saveload_fmap (void);
void test_ncm_mset_load_twice_same_ser (void);
void test_ncm_mset_getsetitem (void);
void test_ncm_mset_fparam_get_fpi_range (void);
void test_ncm_mset_fparam_get_fpi_range_subprocess (void);
void test_ncm_mset_ns_by_negative_id_subprocess (void);

void test_ncm_mset_traps (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_invalid_get (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_invalid_stack (TestNcmMSet *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/mset/setpeek", TestNcmMSet, NULL,
              &test_ncm_mset_new,
              &test_ncm_mset_setpeek,
              &test_ncm_mset_free);

  g_test_add ("/ncm/mset/setpospeek", TestNcmMSet, NULL,
              &test_ncm_mset_new,
              &test_ncm_mset_setpospeek,
              &test_ncm_mset_free);

  g_test_add ("/ncm/mset/pushpeek", TestNcmMSet, NULL,
              &test_ncm_mset_new,
              &test_ncm_mset_pushpeek,
              &test_ncm_mset_free);

  g_test_add ("/ncm/mset/fparams", TestNcmMSet, NULL,
              &test_ncm_mset_new,
              &test_ncm_mset_fparams,
              &test_ncm_mset_free);

  g_test_add ("/ncm/mset/dup", TestNcmMSet, NULL,
              &test_ncm_mset_new,
              &test_ncm_mset_dup,
              &test_ncm_mset_free);

  g_test_add ("/ncm/mset/shallow_copy", TestNcmMSet, NULL,
              &test_ncm_mset_new,
              &test_ncm_mset_shallow_copy,
              &test_ncm_mset_free);

  g_test_add ("/ncm/mset/saveload", TestNcmMSet, NULL,
              &test_ncm_mset_new,
              &test_ncm_mset_saveload,
              &test_ncm_mset_free);

  g_test_add_func ("/ncm/mset/saveload/submodel_reparam", &test_ncm_mset_saveload_submodel_reparam);
  g_test_add_func ("/ncm/mset/submodel/unslotted_attach", &test_ncm_mset_unslotted_submodel_attach);
  g_test_add_func ("/ncm/mset/submodel/unslotted_attach/subprocess", &test_ncm_mset_unslotted_submodel_attach_subprocess);
  g_test_add_func ("/ncm/mset/load/unslotted_submodel_group", &test_ncm_mset_load_unslotted_submodel_group);
  g_test_add_func ("/ncm/mset/submodel/two_level_slots", &test_ncm_mset_two_level_submodel_slots);
  g_test_add_func ("/ncm/mset/peek_by_name/unknown", &test_ncm_mset_peek_by_name_unknown);
  g_test_add_func ("/ncm/mset/remove/host", &test_ncm_mset_remove_host);
  g_test_add_func ("/ncm/mset/split_full_name", &test_ncm_mset_split_full_name);
  g_test_add_func ("/ncm/mset/set_fmap/updates_models", &test_ncm_mset_set_fmap_updates_models);
  g_test_add_func ("/ncm/mset/set_fmap/invalid", &test_ncm_mset_set_fmap_invalid);
  g_test_add_func ("/ncm/mset/stack_bound", &test_ncm_mset_stack_bound);
  g_test_add_func ("/ncm/mset/id_lookups", &test_ncm_mset_id_lookups);
  g_test_add_func ("/ncm/mset/ns_by_negative_id", &test_ncm_mset_ns_by_negative_id);
  g_test_add_func ("/ncm/mset/set_fmap/update_models", &test_ncm_mset_set_fmap_update_models);
  g_test_add_func ("/ncm/mset/max_model_nick", &test_ncm_mset_max_model_nick);
  g_test_add_func ("/ncm/mset/params_pretty_print", &test_ncm_mset_params_pretty_print);
  g_test_add_func ("/ncm/mset/fparam/lookups", &test_ncm_mset_fparam_lookups);
  g_test_add_func ("/ncm/mset/saveload/fmap", &test_ncm_mset_saveload_fmap);
  g_test_add_func ("/ncm/mset/load/twice_same_ser", &test_ncm_mset_load_twice_same_ser);
  g_test_add_func ("/ncm/mset/getsetitem", &test_ncm_mset_getsetitem);
  g_test_add_func ("/ncm/mset/fparam/get_fpi_range", &test_ncm_mset_fparam_get_fpi_range);
  g_test_add_func ("/ncm/mset/fparam/get_fpi_range/subprocess", &test_ncm_mset_fparam_get_fpi_range_subprocess);
  g_test_add_func ("/ncm/mset/ns_by_negative_id/subprocess", &test_ncm_mset_ns_by_negative_id_subprocess);

  g_test_add ("/ncm/mset/traps", TestNcmMSet, NULL,
              &test_ncm_mset_new,
              &test_ncm_mset_traps,
              &test_ncm_mset_free);

  g_test_add ("/ncm/mset/invalid/get/subprocess", TestNcmMSet, NULL,
              &test_ncm_mset_new,
              &test_ncm_mset_invalid_get,
              &test_ncm_mset_free);
  g_test_add ("/ncm/mset/invalid/stack/subprocess", TestNcmMSet, NULL,
              &test_ncm_mset_new,
              &test_ncm_mset_invalid_stack,
              &test_ncm_mset_free);

  g_test_run ();
}

void
test_ncm_mset_new (TestNcmMSet *test, gconstpointer pdata)
{
  test->mset         = ncm_mset_empty_new ();
  test->ma           = g_ptr_array_new ();
  test->ma_destroyed = g_array_new (FALSE, TRUE, sizeof (gboolean));

  g_assert_true (test->mset != NULL);
  g_assert_true (NCM_IS_MSET (test->mset));

  g_assert_cmpuint (ncm_mset_total_len (test->mset), ==, 0);

  {
    NcHICosmoLCDM *cosmo = nc_hicosmo_lcdm_new ();
    gboolean f           = FALSE;

    ncm_mset_set (test->mset, NCM_MODEL (cosmo), NULL);
    ncm_mset_set (test->mset, NCM_MODEL (cosmo), NULL);

    g_ptr_array_add (test->ma, cosmo);
    g_array_append_val (test->ma_destroyed, f);

    nc_hicosmo_free (NC_HICOSMO (cosmo));
  }

  {
    NcDistance *dist    = nc_distance_new (5.0);
    NcSNIADistCov *snia = nc_snia_dist_cov_new (dist, 4);
    gboolean f          = FALSE;

    ncm_mset_set (test->mset, NCM_MODEL (snia), NULL);

    g_ptr_array_add (test->ma, snia);
    g_array_append_val (test->ma_destroyed, f);

    nc_snia_dist_cov_free (snia);
    nc_distance_free (dist);
  }
}

void
_set_destroyed (gpointer b)
{
  gboolean *destroyed = b;

  *destroyed = TRUE;
}

void
test_ncm_mset_free (TestNcmMSet *test, gconstpointer pdata)
{
  NcmMSet *mset      = test->mset;
  gboolean destroyed = FALSE;
  guint i;

  g_object_set_data_full (G_OBJECT (mset), "test-destroy", &destroyed, _set_destroyed);

  g_assert_cmpuint (test->ma_destroyed->len, ==, test->ma->len);

  for (i = 0; i < test->ma_destroyed->len; i++)
  {
    NcmModel *model = g_ptr_array_index (test->ma, i);

    g_array_index (test->ma_destroyed, gboolean, i) = FALSE;
    g_object_set_data_full (G_OBJECT (model), "test-destroy", &g_array_index (test->ma_destroyed, gboolean, i), _set_destroyed);
  }

  ncm_mset_free (mset);
  g_assert_true (destroyed);

  for (i = 0; i < test->ma_destroyed->len; i++)
  {
    g_assert_true (g_array_index (test->ma_destroyed, gboolean, i));
  }

  g_ptr_array_unref (test->ma);
  g_array_unref (test->ma_destroyed);
}

void
test_ncm_mset_setpeek (TestNcmMSet *test, gconstpointer pdata)
{
  NcClusterMass *mass = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassLnnormal"));
  gboolean f          = FALSE;

  g_assert_true (mass != NULL);
  g_assert_true (NC_IS_CLUSTER_MASS (mass));

  ncm_mset_set (test->mset, NCM_MODEL (mass), NULL);

  g_ptr_array_add (test->ma, mass);
  g_array_append_val (test->ma_destroyed, f);

  g_assert_true (ncm_mset_peek (test->mset, nc_cluster_mass_id ()) == NCM_MODEL (mass));

  nc_cluster_mass_free (mass);
}

void
test_ncm_mset_setpospeek (TestNcmMSet *test, gconstpointer pdata)
{
  NcClusterMass *mass = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassLnnormal"));
  gboolean f          = FALSE;

  g_assert_true (mass != NULL);
  g_assert_true (NC_IS_CLUSTER_MASS (mass));

  ncm_mset_set_pos (test->mset, NCM_MODEL (mass), 5, NULL);

  g_ptr_array_add (test->ma, mass);
  g_array_append_val (test->ma_destroyed, f);

  g_assert_true (ncm_mset_peek_pos (test->mset, nc_cluster_mass_id (), 5) == NCM_MODEL (mass));

  nc_cluster_mass_free (mass);
}

void
test_ncm_mset_pushpeek (TestNcmMSet *test, gconstpointer pdata)
{
  NcClusterMass *mass = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassLnnormal"));
  gboolean f          = FALSE;

  g_assert_true (mass != NULL);
  g_assert_true (NC_IS_CLUSTER_MASS (mass));

  ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
  ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
  g_ptr_array_add (test->ma, mass);
  g_array_append_val (test->ma_destroyed, f);

  g_assert_true (ncm_mset_peek (test->mset, nc_cluster_mass_id ()) == NCM_MODEL (mass));
  g_assert_true (ncm_mset_peek_pos (test->mset, nc_cluster_mass_id (), 1) == NCM_MODEL (mass));

  nc_cluster_mass_free (mass);
}

void
test_ncm_mset_fparams (TestNcmMSet *test, gconstpointer pdata)
{
  NcClusterMass *mass   = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassLnnormal"));
  NcClusterMass *benson = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassBenson"));
  gboolean f            = FALSE;

  g_assert_true (mass != NULL);
  g_assert_true (NC_IS_CLUSTER_MASS (mass));
  g_assert_true (benson != NULL);
  g_assert_true (NC_IS_CLUSTER_MASS (benson));

  ncm_mset_set_pos (test->mset, NCM_MODEL (mass), 10, NULL);

  ncm_mset_push (test->mset, NCM_MODEL (benson), NULL);
  ncm_mset_push (test->mset, NCM_MODEL (benson), NULL);
  ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
  ncm_mset_push (test->mset, NCM_MODEL (benson), NULL);

  ncm_mset_set_pos (test->mset, NCM_MODEL (mass), 1, NULL);

  g_ptr_array_add (test->ma, mass);
  g_array_append_val (test->ma_destroyed, f);

  g_ptr_array_add (test->ma, benson);
  g_array_append_val (test->ma_destroyed, f);

  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);

  g_assert_cmpuint (ncm_mset_total_len (test->mset), ==, ncm_mset_fparam_len (test->mset));

  ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
  g_assert_cmpuint (ncm_mset_total_len (test->mset), ==, ncm_mset_fparam_len (test->mset));

  ncm_mset_param_set_ftype (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 1), 0, NCM_PARAM_TYPE_FIXED);

  g_assert_cmpint (ncm_mset_param_get_ftype (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 1), 0, NULL), ==, NCM_PARAM_TYPE_FIXED);
  g_assert_cmpint (ncm_mset_param_get_ftype (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 2), 0, NULL), ==, NCM_PARAM_TYPE_FIXED);
  g_assert_cmpint (ncm_mset_param_get_ftype (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 4), 0, NULL), ==, NCM_PARAM_TYPE_FIXED);
  g_assert_cmpint (ncm_mset_param_get_ftype (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 10), 0, NULL), ==, NCM_PARAM_TYPE_FIXED);

  g_assert_cmpuint (ncm_mset_total_len (test->mset), ==, ncm_mset_fparam_len (test->mset) + 1);

  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FIXED);
  g_assert_cmpuint (ncm_mset_fparam_len (test->mset), ==, 0);

  ncm_mset_param_set_ftype (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 1), 0, NCM_PARAM_TYPE_FREE);
  g_assert_cmpuint (ncm_mset_fparam_len (test->mset), ==, 1);

  {
    const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi (test->mset, 0);

    g_assert_cmpuint (pi->mid, ==, NCM_MSET_MID (nc_cluster_mass_id (), 1));
    g_assert_cmpuint (pi->pid, ==, 0);
  }

  ncm_mset_fparam_set (test->mset, 0, 123.505);
  ncm_assert_cmpdouble (ncm_mset_fparam_get (test->mset, 0), ==, 123.505);
  ncm_assert_cmpdouble (ncm_mset_param_get (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 1), 0), ==, 123.505);
  ncm_assert_cmpdouble (ncm_mset_param_get (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 2), 0), ==, 123.505);
  ncm_assert_cmpdouble (ncm_mset_param_get (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 4), 0), ==, 123.505);
  ncm_assert_cmpdouble (ncm_mset_param_get (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 10), 0), ==, 123.505);

  nc_cluster_mass_free (mass);
  nc_cluster_mass_free (benson);
}

void
test_ncm_mset_dup (TestNcmMSet *test, gconstpointer pdata)
{
  NcClusterMass *mass = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassLnnormal"));
  gboolean f          = FALSE;

  g_assert_true (mass != NULL);
  g_assert_true (NC_IS_CLUSTER_MASS (mass));

  ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
  ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
  g_ptr_array_add (test->ma, mass);
  g_array_append_val (test->ma_destroyed, f);

  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);
  {
    guint i;

    for (i = 0; i < ncm_mset_fparam_len (test->mset); i++)
    {
      const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi (test->mset, i);
      gboolean free           = g_test_rand_double_range (0.0, 1.0) < 0.5;

      ncm_mset_param_set_ftype (test->mset, pi->mid, pi->pid, free ? NCM_PARAM_TYPE_FREE : NCM_PARAM_TYPE_FIXED);
    }
  }
  ncm_mset_prepare_fparam_map (test->mset);

  {
    NcClusterMass *benson = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassBenson"));
    NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    NcmMSet *mset_dup     = ncm_mset_dup (test->mset, ser);

    g_assert_true (benson != NULL);
    g_assert_true (NC_IS_CLUSTER_MASS (benson));

    g_assert_true (ncm_mset_cmp (test->mset, mset_dup, FALSE));
    g_assert_true (ncm_mset_cmp (test->mset, mset_dup, TRUE));

    {
      guint i;

      g_assert_cmpuint (ncm_mset_nmodels (test->mset), ==, ncm_mset_nmodels (mset_dup));

      for (i = 0; i < ncm_mset_nmodels (test->mset); i++)
      {
        NcmModel *model0 = ncm_mset_peek_array_pos (test->mset, i);
        NcmModel *model1 = ncm_mset_peek_array_pos (mset_dup, i);
        guint pid;

        for (pid = 0; pid < ncm_model_len (model0); pid++)
        {
          ncm_assert_cmpdouble (ncm_model_param_get (model0, pid), ==, ncm_model_param_get (model1, pid));
        }
      }

      g_assert_cmpuint (ncm_mset_fparams_len (test->mset), ==, ncm_mset_fparams_len (mset_dup));

      for (i = 0; i < ncm_mset_fparams_len (test->mset); i++)
      {
        g_assert_cmpstr (ncm_mset_fparam_full_name (test->mset, i), ==, ncm_mset_fparam_full_name (mset_dup, i));
        ncm_assert_cmpdouble (ncm_mset_fparam_get (test->mset, i), ==, ncm_mset_fparam_get (mset_dup, i));
      }
    }

    ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
    g_assert_true (!ncm_mset_cmp (test->mset, mset_dup, FALSE));
    g_assert_true (!ncm_mset_cmp (test->mset, mset_dup, TRUE));

    /*g_ptr_array_add (test->ma, benson);*/
    /*g_array_append_val (test->ma_destroyed, f);*/

    ncm_mset_push (mset_dup, NCM_MODEL (benson), NULL);
    g_assert_true (ncm_mset_cmp (test->mset, mset_dup, FALSE));
    g_assert_true (!ncm_mset_cmp (test->mset, mset_dup, TRUE));

    {
      gboolean destroyed = FALSE;

      g_object_set_data_full (G_OBJECT (ser), "test-destroy", &destroyed, _set_destroyed);
      ncm_serialize_clear (&ser);
      g_assert_true (destroyed);
    }

    {
      gboolean destroyed = FALSE;

      g_object_set_data_full (G_OBJECT (mset_dup), "test-destroy", &destroyed, _set_destroyed);
      ncm_mset_clear (&mset_dup);
      g_assert_true (destroyed);
    }

    {
      gboolean destroyed = FALSE;

      g_object_set_data_full (G_OBJECT (benson), "test-destroy", &destroyed, _set_destroyed);
      nc_cluster_mass_free (benson);
      g_assert_true (destroyed);
    }
  }

  nc_cluster_mass_free (mass);
}

void
test_ncm_mset_shallow_copy (TestNcmMSet *test, gconstpointer pdata)
{
  NcClusterMass *mass = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassLnnormal"));
  gboolean f          = FALSE;

  g_assert_true (mass != NULL);
  g_assert_true (NC_IS_CLUSTER_MASS (mass));

  ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
  ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
  g_ptr_array_add (test->ma, mass);
  g_array_append_val (test->ma_destroyed, f);

  {
    NcmMSet *mset_dup = ncm_mset_shallow_copy (test->mset, NULL);

    g_assert_true (ncm_mset_cmp (test->mset, mset_dup, FALSE));
    g_assert_true (ncm_mset_cmp (test->mset, mset_dup, TRUE));

    g_assert_true (ncm_mset_is_subset (test->mset, mset_dup));

    {
      guint i;

      g_assert_cmpuint (ncm_mset_nmodels (test->mset), ==, ncm_mset_nmodels (mset_dup));

      for (i = 0; i < ncm_mset_nmodels (test->mset); i++)
      {
        NcmModel *model0 = ncm_mset_peek_array_pos (test->mset, i);
        NcmModel *model1 = ncm_mset_peek_array_pos (mset_dup, i);

        g_assert_true (model0 == model1);
      }
    }

    while (ncm_mset_nmodels (mset_dup) > 0)
    {
      ncm_mset_remove (mset_dup, ncm_mset_get_mid_array_pos (mset_dup, 0));
      g_assert_true (ncm_mset_is_subset (test->mset, mset_dup));
    }

    {
      gboolean destroyed = FALSE;

      g_object_set_data_full (G_OBJECT (mset_dup), "test-destroy", &destroyed, _set_destroyed);
      ncm_mset_clear (&mset_dup);
      g_assert_true (destroyed);
    }
  }

  nc_cluster_mass_free (mass);
}

void
test_ncm_mset_saveload (TestNcmMSet *test, gconstpointer pdata)
{
  NcClusterMass *benson = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassBenson"));
  NcClusterMass *mass   = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassLnnormal"));
  gboolean f            = FALSE;

  g_assert_true (mass != NULL);
  g_assert_true (NC_IS_CLUSTER_MASS (mass));
  g_assert_true (benson != NULL);
  g_assert_true (NC_IS_CLUSTER_MASS (benson));

  ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
  ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
  ncm_mset_push (test->mset, NCM_MODEL (benson), NULL);

  g_ptr_array_add (test->ma, mass);
  g_array_append_val (test->ma_destroyed, f);

  g_ptr_array_add (test->ma, benson);
  g_array_append_val (test->ma_destroyed, f);

  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    gchar *tmp_dir    = g_dir_make_tmp ("test_ncm_mset_saved_XXXXXX", NULL);
    gchar *filename   = g_strdup_printf ("%s/test_ncm_mset_saved.mset", tmp_dir);
    NcmMSet *mset_dup = NULL;

    ncm_mset_save (test->mset, ser, filename, TRUE, NULL);
    mset_dup = ncm_mset_load (filename, ser, NULL);

    g_assert_true (ncm_mset_cmp (test->mset, mset_dup, FALSE));
    g_assert_true (ncm_mset_cmp (test->mset, mset_dup, TRUE));

    {
      guint i;

      g_assert_cmpuint (ncm_mset_nmodels (test->mset), ==, ncm_mset_nmodels (mset_dup));

      for (i = 0; i < ncm_mset_nmodels (test->mset); i++)
      {
        NcmModel *model0 = ncm_mset_peek_array_pos (test->mset, i);
        NcmModel *model1 = ncm_mset_peek_array_pos (mset_dup, i);
        guint pid;

        for (pid = 0; pid < ncm_model_len (model0); pid++)
        {
          ncm_assert_cmpdouble (ncm_model_param_get (model0, pid), ==, ncm_model_param_get (model1, pid));
        }
      }
    }

    ncm_mset_push (test->mset, NCM_MODEL (mass), NULL);
    g_assert_true (!ncm_mset_cmp (test->mset, mset_dup, FALSE));
    g_assert_true (!ncm_mset_cmp (test->mset, mset_dup, TRUE));

    /*g_ptr_array_add (test->ma, benson);*/
    /*g_array_append_val (test->ma_destroyed, f);*/

    ncm_mset_push (mset_dup, NCM_MODEL (benson), NULL);
    g_assert_true (ncm_mset_cmp (test->mset, mset_dup, FALSE));
    g_assert_true (!ncm_mset_cmp (test->mset, mset_dup, TRUE));

    {
      gboolean destroyed = FALSE;

      g_object_set_data_full (G_OBJECT (mset_dup), "test-destroy", &destroyed, _set_destroyed);
      ncm_mset_clear (&mset_dup);
      g_assert_true (destroyed);
    }

    {
      gboolean destroyed = FALSE;

      g_object_set_data_full (G_OBJECT (ser), "test-destroy", &destroyed, _set_destroyed);
      ncm_serialize_clear (&ser);
      g_assert_true (destroyed);
    }

    g_unlink (filename);
    g_rmdir (tmp_dir);

    g_free (filename);
    g_free (tmp_dir);
  }

  nc_cluster_mass_free (benson);
  nc_cluster_mass_free (mass);
}

/* A two-level slotted host family: the child class extends its parent's
 * declared submodel slots, exercising the parent-slot copy in
 * ncm_model_class_add_submodels(). */
#define TEST_TYPE_HOST_SLOTTED (test_host_slotted_get_type ())
G_DECLARE_DERIVABLE_TYPE (TestHostSlotted, test_host_slotted, TEST, HOST_SLOTTED, NcmModel)

struct _TestHostSlottedClass
{
  NcmModelClass parent_class;
};

G_DEFINE_TYPE (TestHostSlotted, test_host_slotted, NCM_TYPE_MODEL)
NCM_MSET_MODEL_REGISTER_ID (test_host_slotted, TEST_TYPE_HOST_SLOTTED)

#define TEST_TYPE_SUB_P (test_sub_p_get_type ())
G_DECLARE_FINAL_TYPE (TestSubP, test_sub_p, TEST, SUB_P, NcmModel)

struct _TestSubP
{
  NcmModel parent_instance;
};

G_DEFINE_TYPE (TestSubP, test_sub_p, NCM_TYPE_MODEL)
NCM_MSET_MODEL_REGISTER_ID (test_sub_p, TEST_TYPE_SUB_P)

#define TEST_TYPE_SUB_Q (test_sub_q_get_type ())
G_DECLARE_FINAL_TYPE (TestSubQ, test_sub_q, TEST, SUB_Q, NcmModel)

struct _TestSubQ
{
  NcmModel parent_instance;
};

G_DEFINE_TYPE (TestSubQ, test_sub_q, NCM_TYPE_MODEL)
NCM_MSET_MODEL_REGISTER_ID (test_sub_q, TEST_TYPE_SUB_Q)

static void
test_sub_p_init (TestSubP *sub)
{
}

static void
test_sub_p_class_init (TestSubPClass *klass)
{
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  ncm_model_class_set_name_nick (model_class, "Test slotted submodel P", "TestSubP");
  ncm_model_class_add_params (model_class, 1, 0, 1);
  ncm_model_class_set_sparam (model_class, 0, "x_0^{(p)}", "x0", -10.0, 10.0, 1.0, 0.0, 1.0, NCM_PARAM_TYPE_FIXED);

  ncm_mset_model_register_id (model_class, "TestSubP", "Test slotted submodel P", NULL, FALSE,
                              test_host_slotted_id ());

  ncm_model_class_check_params_info (model_class);
}

static void
test_sub_q_init (TestSubQ *sub)
{
}

static void
test_sub_q_class_init (TestSubQClass *klass)
{
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  ncm_model_class_set_name_nick (model_class, "Test slotted submodel Q", "TestSubQ");
  ncm_model_class_add_params (model_class, 1, 0, 1);
  ncm_model_class_set_sparam (model_class, 0, "x_0^{(q)}", "x0", -10.0, 10.0, 1.0, 0.0, 1.0, NCM_PARAM_TYPE_FIXED);

  ncm_mset_model_register_id (model_class, "TestSubQ", "Test slotted submodel Q", NULL, FALSE,
                              test_host_slotted_id ());

  ncm_model_class_check_params_info (model_class);
}

static void
test_host_slotted_init (TestHostSlotted *host)
{
}

static void
test_host_slotted_class_init (TestHostSlottedClass *klass)
{
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  ncm_model_class_set_name_nick (model_class, "Test slotted host", "TestHostSlotted");
  ncm_model_class_add_params (model_class, 0, 0, 1);

  ncm_model_class_add_submodels (model_class, 1);
  ncm_model_class_set_submodel (model_class, 0, "sub-p", "sub-p", TEST_TYPE_SUB_P);

  ncm_mset_model_register_id (model_class, "TestHostSlotted", "Test slotted host", NULL, FALSE,
                              NCM_MSET_MODEL_MAIN);

  ncm_model_class_check_params_info (model_class);
}

#define TEST_TYPE_HOST_SLOTTED_CHILD (test_host_slotted_child_get_type ())
G_DECLARE_FINAL_TYPE (TestHostSlottedChild, test_host_slotted_child, TEST, HOST_SLOTTED_CHILD, TestHostSlotted)

struct _TestHostSlottedChild
{
  TestHostSlotted parent_instance;
};

G_DEFINE_TYPE (TestHostSlottedChild, test_host_slotted_child, TEST_TYPE_HOST_SLOTTED)

static void
test_host_slotted_child_init (TestHostSlottedChild *host)
{
}

static void
test_host_slotted_child_class_init (TestHostSlottedChildClass *klass)
{
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  ncm_model_class_set_name_nick (model_class, "Test slotted host child", "TestHostSlottedChild");
  ncm_model_class_add_params (model_class, 0, 0, 1);

  /* Extends the parent's one declared slot with a second one. */
  ncm_model_class_add_submodels (model_class, 2);
  ncm_model_class_set_submodel (model_class, 1, "sub-q", "sub-q", TEST_TYPE_SUB_Q);

  ncm_model_class_check_params_info (model_class);
}

/*
 * A minimal submodel type with NcHICosmo as its main model but NO typed slot
 * declared on NcHICosmo -- the only in-tree way to exercise ncm_mset_load()'s
 * pass 3, which attaches pass-1 submodels that pass-2 slot injection did not
 * consume.
 */
#define TEST_TYPE_UNSLOTTED_SUB (test_unslotted_sub_get_type ())
G_DECLARE_FINAL_TYPE (TestUnslottedSub, test_unslotted_sub, TEST, UNSLOTTED_SUB, NcmModel)

struct _TestUnslottedSub
{
  NcmModel parent_instance;
};

G_DEFINE_TYPE (TestUnslottedSub, test_unslotted_sub, NCM_TYPE_MODEL)
NCM_MSET_MODEL_REGISTER_ID (test_unslotted_sub, TEST_TYPE_UNSLOTTED_SUB)

static void
test_unslotted_sub_init (TestUnslottedSub *sub)
{
}

static void
test_unslotted_sub_class_init (TestUnslottedSubClass *klass)
{
  NcmModelClass *model_class = NCM_MODEL_CLASS (klass);

  ncm_model_class_set_name_nick (model_class, "Unslotted test submodel", "TestUnslottedSub");
  ncm_model_class_add_params (model_class, 1, 0, 1);
  ncm_model_class_set_sparam (model_class, 0, "x_0", "x0", -10.0, 10.0, 1.0, 0.0, 1.0, NCM_PARAM_TYPE_FIXED);

  ncm_mset_model_register_id (model_class,
                              "TestUnslottedSub",
                              "Unslotted test submodel",
                              NULL,
                              FALSE,
                              nc_hicosmo_id ());

  ncm_model_class_check_params_info (model_class);
}

void
test_ncm_mset_unslotted_submodel_attach_subprocess (void)
{
  NcHICosmo *cosmo      = NC_HICOSMO (nc_hicosmo_lcdm_new_full (NULL, NULL, NULL));
  TestUnslottedSub *sub = g_object_new (TEST_TYPE_UNSLOTTED_SUB, NULL);

  /* Submodels are construction-fixed for every type: a type without a
   * declared slot has no legal attachment path at all. */
  ncm_model_add_submodel (NCM_MODEL (cosmo), NCM_MODEL (sub));

  ncm_model_free (NCM_MODEL (sub));
  nc_hicosmo_free (cosmo);
}

void
test_ncm_mset_unslotted_submodel_attach (void)
{
  g_test_trap_subprocess ("/ncm/mset/submodel/unslotted_attach/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*construction-fixed*");
}

void
test_ncm_mset_load_unslotted_submodel_group (void)
{
  /* A file carrying a submodel group whose type has no slot on its main
   * model must fail loudly, never be silently dropped. */
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_lcdm_new_full (NULL, NULL, NULL));
  GError *error    = NULL;

  /* Class init registers the model id; it must exist before ncm_mset_load
   * meets the group. */
  GType sub_class = (GType) g_type_class_ref (TEST_TYPE_UNSLOTTED_SUB);

  {
    NcmMSet *mset     = ncm_mset_empty_new ();
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    gchar *tmp_dir    = g_dir_make_tmp ("test_ncm_mset_unslotted_XXXXXX", NULL);
    gchar *filename   = g_build_filename (tmp_dir, "unslotted.mset", NULL);

    ncm_mset_set (mset, NCM_MODEL (cosmo), &error);
    g_assert_no_error (error);
    ncm_mset_save (mset, ser, filename, FALSE, &error);
    g_assert_no_error (error);

    {
      gchar *contents = NULL;
      gchar *extended = NULL;

      g_assert_true (g_file_get_contents (filename, &contents, NULL, NULL));
      extended = g_strconcat (contents,
                              "\n[TestUnslottedSub]\n"
                              "TestUnslottedSub=TestUnslottedSub\n",
                              NULL);
      g_assert_true (g_file_set_contents (filename, extended, -1, NULL));
      g_free (contents);
      g_free (extended);
    }

    {
      NcmSerialize *ser2 = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
      NcmMSet *mset2     = ncm_mset_load (filename, ser2, &error);

      g_assert_true (mset2 == NULL);
      g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_KEY_FILE_INVALID);
      g_assert_nonnull (strstr (error->message, "was not consumed by a construction-time"));
      g_clear_error (&error);
      ncm_serialize_free (ser2);
    }

    g_unlink (filename);
    g_rmdir (tmp_dir);
    g_free (filename);
    g_free (tmp_dir);
    ncm_serialize_free (ser);
    ncm_mset_free (mset);
  }

  g_type_class_unref ((gpointer) sub_class);
  nc_hicosmo_free (cosmo);
}

void
test_ncm_mset_two_level_submodel_slots (void)
{
  TestSubP *sub_p             = g_object_new (TEST_TYPE_SUB_P, NULL);
  TestSubQ *sub_q             = g_object_new (TEST_TYPE_SUB_Q, NULL);
  TestHostSlottedChild *child = g_object_new (TEST_TYPE_HOST_SLOTTED_CHILD,
                                              "sub-p", sub_p,
                                              "sub-q", sub_q,
                                              NULL);
  GError *error = NULL;

  /* Both the inherited and the child-added slot deliver at construction. */
  g_assert_true (ncm_model_peek_submodel_by_mid (NCM_MODEL (child), test_sub_p_id ()) == NCM_MODEL (sub_p));
  g_assert_true (ncm_model_peek_submodel_by_mid (NCM_MODEL (child), test_sub_q_id ()) == NCM_MODEL (sub_q));

  /* Qualified access works through both slots. */
  {
    NcmModel *target = NULL;
    guint i          = 99;

    g_assert_true (ncm_model_param_index_from_name_full (NCM_MODEL (child), "sub-p:x0", &target, &i, &error));
    g_assert_no_error (error);
    g_assert_true (target == NCM_MODEL (sub_p));

    g_assert_true (ncm_model_param_index_from_name_full (NCM_MODEL (child), "sub-q:x0", &target, &i, &error));
    g_assert_no_error (error);
    g_assert_true (target == NCM_MODEL (sub_q));

    /* Both attached submodels carry "x0": the bare name is ambiguous. */
    g_assert_false (ncm_model_param_index_from_name_full (NCM_MODEL (child), "x0", &target, &i, &error));
    g_assert_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_NAME_AMBIGUOUS);
    g_clear_error (&error);

    /* Known slot, unknown parameter: the resolver's own error. */
    g_assert_false (ncm_model_param_index_from_name_full (NCM_MODEL (child), "sub-p:nope", &target, &i, &error));
    g_assert_error (error, NCM_MODEL_ERROR, NCM_MODEL_ERROR_PARAM_NAME_NOT_FOUND);
    g_clear_error (&error);
  }

  ncm_model_free (NCM_MODEL (sub_p));
  ncm_model_free (NCM_MODEL (sub_q));
  ncm_model_free (NCM_MODEL (child));
}

void
test_ncm_mset_saveload_submodel_reparam (void)
{
  /*
   * Regression test for a real, pre-existing bug found while testing the
   * reparam backpointer migration: every NcmModel carries a generic
   * "submodel-array" boxed property (independent of any typed slot) that
   * ncm_serialize_to_variant() walks and recursively, fully serializes --
   * so when a host (cosmo) is serialized, its submodels (e.g. reion) get
   * fully walked and registered as already-seen instances in @ser's
   * persistent, cross-call name table (populated by
   * NCM_SERIALIZE_OPT_AUTOSAVE_SER, which ncm_mset_save()'s own @ser
   * always sets). Since a submodel is *also* independently serialized as
   * its own top-level NcmMSet entry, ncm_mset_save() already unregistered
   * the submodel itself from @ser right after the host's own walk, to
   * force a full re-serialization instead of a bare name-only reference
   * -- but didn't recurse into objects NESTED inside the submodel's own
   * serialization (e.g. an attached NcmReparam, and that reparam's own
   * "T"/"v" properties), so anything nested came out as an empty,
   * unusable reference after a save/load round-trip. Fixed in
   * ncm_mset_save() via a recursive unregistration helper
   * (_ncm_mset_serialize_unset_recursive()) that walks every
   * GObject-valued readable property reachable from the submodel, not
   * just the submodel itself.
   *
   * Uses a host-independent NcmReparamLinear here (not
   * NcHIReionCambReparamTau) specifically to isolate this bug from a
   * separate, known limitation: ncm_mset_load()'s two-pass reorder builds
   * every submodel (including setting its "reparam" property) in pass 1,
   * before pass 2 attaches it to its host in pass 2 -- so a reparam whose
   * own old2new/new2old needs ncm_model_peek_host() (like the tau reparam)
   * cannot currently be pre-attached before a save/load round-trip; it
   * must be (re-)attached via the public API (e.g. z_to_tau()) after
   * loading, not expected to survive serialization pre-attached.
   */
  NcHIReion *reion = NC_HIREION (nc_hireion_camb_new ());
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_lcdm_new_full (reion, NULL, NULL));
  GError *error    = NULL;

  {
    NcmMatrix *T = ncm_matrix_new (2, 2);
    NcmVector *v = ncm_vector_new (2);
    NcmReparam *rp;

    ncm_matrix_set_identity (T);
    ncm_vector_set_zero (v);

    rp = NCM_REPARAM (ncm_reparam_linear_new (2, T, v));
    ncm_model_set_reparam (NCM_MODEL (reion), rp, &error);
    g_assert_no_error (error);
    ncm_reparam_free (rp);

    ncm_matrix_free (T);
    ncm_vector_free (v);
  }

  ncm_model_param_set (NCM_MODEL (reion), NC_HIREION_CAMB_HII_HEII_Z, 9.5);

  {
    NcmMSet *mset     = ncm_mset_empty_new ();
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    gchar *tmp_dir    = g_dir_make_tmp ("test_ncm_mset_saved_reparam_XXXXXX", NULL);
    gchar *filename   = g_strdup_printf ("%s/test_ncm_mset_saved_reparam.mset", tmp_dir);
    NcmMSet *mset_dup = NULL;
    NcmModel *cosmo_dup;
    NcmModel *reion_dup;
    NcmReparam *reparam_dup;

    ncm_mset_set (mset, NCM_MODEL (cosmo), &error);
    g_assert_no_error (error);

    ncm_mset_save (mset, ser, filename, TRUE, NULL);
    mset_dup = ncm_mset_load (filename, ser, &error);
    g_assert_no_error (error);
    g_assert_true (mset_dup != NULL);

    g_assert_true (ncm_mset_cmp (mset, mset_dup, TRUE));

    cosmo_dup = ncm_mset_peek (mset_dup, nc_hicosmo_id ());
    g_assert_true (cosmo_dup != NULL);

    reion_dup = ncm_model_peek_submodel_by_mid (cosmo_dup, nc_hireion_id ());
    g_assert_true (reion_dup != NULL);

    /* The loaded reion must have its host resolvable again -- the whole
     * point of the backpointer -- and its reparam must have survived as a
     * full, working definition, not an empty reference. */
    g_assert_true (ncm_model_peek_host (reion_dup) == cosmo_dup);

    reparam_dup = ncm_model_peek_reparam (reion_dup);
    g_assert_true (NCM_IS_REPARAM_LINEAR (reparam_dup));
    g_assert_cmpuint (ncm_model_len (reion_dup), ==, ncm_model_len (NCM_MODEL (reion)));

    ncm_assert_cmpdouble (ncm_model_param_get (reion_dup, NC_HIREION_CAMB_HII_HEII_Z), ==,
                          ncm_model_param_get (NCM_MODEL (reion), NC_HIREION_CAMB_HII_HEII_Z));

    ncm_mset_free (mset_dup);
    ncm_serialize_free (ser);

    g_unlink (filename);
    g_rmdir (tmp_dir);

    g_free (filename);
    g_free (tmp_dir);

    ncm_mset_free (mset);
  }

  nc_hicosmo_free (cosmo);
  nc_hireion_free (reion);
}

void
test_ncm_mset_traps (TestNcmMSet *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/mset/invalid/get/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/mset/invalid/stack/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_mset_invalid_get (TestNcmMSet *test, gconstpointer pdata)
{
  g_assert_true (ncm_mset_get (test->mset, 34 * NCM_MSET_MAX_STACKSIZE + 5) != NULL);
}

void
test_ncm_mset_invalid_stack (TestNcmMSet *test, gconstpointer pdata)
{
  NcHICosmoLCDM *cosmo = nc_hicosmo_lcdm_new ();

  ncm_mset_push (test->mset, NCM_MODEL (cosmo), NULL);

  nc_hicosmo_free (NC_HICOSMO (cosmo));
}

/* An unregistered namespace gives no model, with or without a stack position, even when
 * namespace id -1 plus that position is the id of a model in the set. */
void
test_ncm_mset_peek_by_name_unknown (void)
{
  NcmModel *mvnd = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmMSet *mset  = ncm_mset_new (mvnd, NULL, NULL);
  const gint id  = ncm_model_mvnd_id ();
  gchar *name    = g_strdup_printf ("NoSuchModel:%d", id + 1);
  GError *error  = NULL;

  g_assert_true (ncm_mset_peek_by_name (mset, "NcmModelMVND", NULL) == mvnd);
  g_assert_null (ncm_mset_peek_by_name (mset, "NoSuchModel", NULL));
  g_assert_null (ncm_mset_peek_by_name (mset, name, NULL));

  g_assert_null (ncm_mset_fetch_by_name (mset, name, &error));
  g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_NOT_SET);
  g_clear_error (&error);

  g_free (name);
  ncm_mset_free (mset);
  ncm_model_free (mvnd);
}

/* Removing a host removes its submodels; the set then matches its serialized copy. */
void
test_ncm_mset_remove_host (void)
{
  NcHIReion *reion     = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim       = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcmModel *cosmo      = NCM_MODEL (nc_hicosmo_de_xcdm_new_full (reion, prim, NULL));
  NcmModel *mvnd       = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmMSet *mset        = ncm_mset_new (cosmo, NULL, mvnd, NULL);
  NcmSerialize *ser    = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  const guint mvnd_len = ncm_model_len (mvnd);
  NcmMSet *dup;

  g_assert_nonnull (ncm_mset_peek (mset, nc_hiprim_id ()));
  g_assert_cmpuint (ncm_mset_total_len (mset), >, mvnd_len);

  ncm_mset_remove (mset, nc_hicosmo_id ());

  g_assert_null (ncm_mset_peek (mset, nc_hicosmo_id ()));
  g_assert_null (ncm_mset_peek (mset, nc_hiprim_id ()));
  g_assert_null (ncm_mset_peek (mset, nc_hireion_id ()));
  g_assert_cmpuint (ncm_mset_nmodels (mset), ==, 1);
  g_assert_cmpuint (ncm_mset_total_len (mset), ==, mvnd_len);

  dup = ncm_mset_dup (mset, ser);
  g_assert_cmpuint (ncm_mset_nmodels (dup), ==, 1);
  g_assert_cmpuint (ncm_mset_total_len (dup), ==, mvnd_len);

  ncm_mset_free (dup);
  ncm_serialize_free (ser);
  ncm_mset_free (mset);
  ncm_model_free (mvnd);
  ncm_model_free (cosmo);
  nc_hiprim_free (prim);
  nc_hireion_free (reion);
}

/* Full names split with and without a stack position; a stack position out of range is
 * an error that leaves the outputs unset. */
void
test_ncm_mset_split_full_name (void)
{
  gchar *ns       = NULL;
  gchar *pname    = NULL;
  guint stackpos  = 99;
  GError *error   = NULL;
  gchar *too_high = g_strdup_printf ("NcHICosmo:%d:H0", NCM_MSET_MAX_STACKSIZE);

  g_assert_true (ncm_mset_split_full_name ("NcHICosmo:H0", &ns, &stackpos, &pname, &error));
  g_assert_no_error (error);
  g_assert_cmpstr (ns, ==, "NcHICosmo");
  g_assert_cmpstr (pname, ==, "H0");
  g_assert_cmpuint (stackpos, ==, 0);
  g_clear_pointer (&ns, g_free);
  g_clear_pointer (&pname, g_free);

  g_assert_true (ncm_mset_split_full_name ("NcHICosmo:2:Omegac", &ns, &stackpos, &pname, &error));
  g_assert_cmpuint (stackpos, ==, 2);
  g_assert_cmpstr (pname, ==, "Omegac");
  g_clear_pointer (&ns, g_free);
  g_clear_pointer (&pname, g_free);

  g_assert_false (ncm_mset_split_full_name ("not a name", &ns, &stackpos, &pname, &error));
  g_assert_no_error (error);

  g_assert_false (ncm_mset_split_full_name (too_high, &ns, &stackpos, &pname, &error));
  g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_FULLNAME_INVALID);
  g_assert_null (ns);
  g_assert_null (pname);
  g_clear_error (&error);

  g_free (too_high);
}

/* After set_fmap moves the free parameter to another model, setting the free
 * parameters updates that model (its pkey changes). */
void
test_ncm_mset_set_fmap_updates_models (void)
{
  NcmModel *a         = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmModel *b         = NCM_MODEL (ncm_model_rosenbrock_new ());
  NcmMSet *mset       = ncm_mset_new (a, NULL, b, NULL);
  NcmVector *x        = ncm_vector_new (1);
  const gchar *fmap[] = {"NcmModelRosenbrock:x1", NULL};
  gchar **got;
  guint64 pkey;

  ncm_model_param_set_ftype (a, 0, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);
  g_assert_cmpuint (ncm_mset_fparams_len (mset), ==, 1);

  ncm_mset_set_fmap (mset, fmap, FALSE, NULL);
  got = ncm_mset_get_fmap (mset);
  g_assert_cmpstr (got[0], ==, "NcmModelRosenbrock:x1");
  g_assert_null (got[1]);
  g_strfreev (got);

  pkey = ncm_model_state_get_pkey (b);
  ncm_vector_set (x, 0, 1.25);
  ncm_mset_fparams_set_vector (mset, x);
  g_assert_cmpfloat (ncm_model_param_get (b, 0), ==, 1.25);
  g_assert_cmpuint (ncm_model_state_get_pkey (b), !=, pkey);

  ncm_vector_free (x);
  ncm_mset_free (mset);
  ncm_model_free (a);
  ncm_model_free (b);
}

/* An unknown or repeated name is an error that leaves the previous map in place. */
void
test_ncm_mset_set_fmap_invalid (void)
{
  NcmModel *a          = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmMSet *mset        = ncm_mset_new (a, NULL, NULL);
  const gchar *good[]  = {"NcmModelMVND:mu_0", NULL};
  const gchar *bad[]   = {"NcmModelMVND:mu_1", "NcmModelMVND:nope", NULL};
  const gchar *twice[] = {"NcmModelMVND:mu_1", "NcmModelMVND:mu_1", NULL};
  GError *error        = NULL;

  ncm_mset_set_fmap (mset, good, FALSE, &error);
  g_assert_no_error (error);

  ncm_mset_set_fmap (mset, bad, FALSE, &error);
  g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_FULLNAME_NOT_FOUND);
  g_clear_error (&error);

  ncm_mset_set_fmap (mset, twice, FALSE, &error);
  g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_FULLNAME_INVALID);
  g_clear_error (&error);

  g_assert_true (ncm_mset_fparam_map_valid (mset));
  g_assert_cmpuint (ncm_mset_fparams_len (mset), ==, 1);
  g_assert_cmpstr (ncm_mset_fparam_full_name (mset, 0), ==, "NcmModelMVND:mu_0");
  g_assert_cmpint (ncm_mset_fparam_get_fpi (mset, ncm_model_mvnd_id (), 0), ==, 0);
  g_assert_cmpint (ncm_mset_fparam_get_fpi (mset, ncm_model_mvnd_id (), 1), ==, -1);

  ncm_mset_free (mset);
  ncm_model_free (a);
}

/* A stack position of NCM_MSET_MAX_STACKSIZE is rejected (it is the next class's
 * id); a non-stackable model at position 1 is rejected without keeping a reference. */
void
test_ncm_mset_stack_bound (void)
{
  NcDistance *dist   = nc_distance_new (3.0);
  NcHaloPosition *hp = nc_halo_position_new (dist);
  NcmModel *mvnd     = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmMSet *mset      = ncm_mset_empty_new ();
  GError *error      = NULL;

  ncm_mset_set_pos (mset, NCM_MODEL (hp), NCM_MSET_MAX_STACKSIZE - 1, &error);
  g_assert_no_error (error);

  ncm_mset_set_pos (mset, NCM_MODEL (hp), NCM_MSET_MAX_STACKSIZE, &error);
  g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_INVALID_ID);
  g_clear_error (&error);
  g_assert_null (ncm_mset_peek (mset, ncm_model_id (NCM_MODEL (hp)) + NCM_MSET_MAX_STACKSIZE));

  ncm_mset_set_pos (mset, mvnd, 1, &error);
  g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_NOT_STACKABLE);
  g_clear_error (&error);
  g_assert_cmpuint (G_OBJECT (mvnd)->ref_count, ==, 1);

  ncm_mset_free (mset);
  ncm_model_free (mvnd);
  nc_halo_position_free (hp);
  nc_distance_free (dist);
}

/* get_id_by_type finds a subclass of a registered class; ns and type by id agree. */
void
test_ncm_mset_id_lookups (void)
{
  g_assert_cmpint (ncm_mset_get_id_by_type (NC_TYPE_HICOSMO_DE_XCDM), ==, nc_hicosmo_id ());
  g_assert_cmpint (ncm_mset_get_id_by_type (NC_TYPE_HICOSMO), ==, nc_hicosmo_id ());
  g_assert_cmpint (ncm_mset_get_id_by_ns ("NcHICosmo"), ==, nc_hicosmo_id ());
  g_assert_cmpstr (ncm_mset_get_ns_by_id (nc_hicosmo_id () + 1), ==, "NcHICosmo");
  g_assert_true (ncm_mset_get_type_by_id (nc_hicosmo_id ()) == NC_TYPE_HICOSMO);
}

void
test_ncm_mset_ns_by_negative_id (void)
{
  g_test_trap_subprocess ("/ncm/mset/ns_by_negative_id/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_ncm_mset_ns_by_negative_id_subprocess (void)
{
  ncm_mset_get_ns_by_id (-1);
}

/* With update_models, the fit types follow the map and the map keeps its order. */
void
test_ncm_mset_set_fmap_update_models (void)
{
  NcmModel *a         = NCM_MODEL (ncm_model_mvnd_new (3));
  NcmMSet *mset       = ncm_mset_new (a, NULL, NULL);
  const gchar *fmap[] = {"NcmModelMVND:mu_2", "NcmModelMVND:mu_1", NULL};

  ncm_model_param_set_ftype (a, 0, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);

  ncm_mset_set_fmap (mset, fmap, TRUE, NULL);

  g_assert_cmpuint (ncm_mset_fparams_len (mset), ==, 2);
  g_assert_cmpstr (ncm_mset_fparam_full_name (mset, 0), ==, "NcmModelMVND:mu_2");
  g_assert_cmpstr (ncm_mset_fparam_full_name (mset, 1), ==, "NcmModelMVND:mu_1");
  g_assert_cmpint (ncm_model_param_get_ftype (a, 0), ==, NCM_PARAM_TYPE_FIXED);
  g_assert_cmpint (ncm_model_param_get_ftype (a, 1), ==, NCM_PARAM_TYPE_FREE);
  g_assert_cmpint (ncm_model_param_get_ftype (a, 2), ==, NCM_PARAM_TYPE_FREE);

  ncm_mset_free (mset);
  ncm_model_free (a);
}

/* Models without parameters count: NcBBNParthenope, a submodel of the cosmology. */
void
test_ncm_mset_max_model_nick (void)
{
  NcmModel *cosmo = NCM_MODEL (nc_hicosmo_de_xcdm_new ());
  NcmMSet *mset   = ncm_mset_new (cosmo, NULL, NULL);
  guint longest   = 0;
  guint i;

  for (i = 0; i < ncm_mset_nmodels (mset); i++)
    longest = GSL_MAX (longest, strlen (ncm_model_nick (ncm_mset_peek_array_pos (mset, i))));

  g_assert_nonnull (ncm_mset_peek (mset, nc_bbn_id ()));
  g_assert_cmpuint (ncm_model_len (ncm_mset_peek (mset, nc_bbn_id ())), ==, 0);
  g_assert_cmpuint (ncm_mset_max_model_nick (mset), ==, longest);

  ncm_mset_free (mset);
  ncm_model_free (cosmo);
}

/* The line after the header starts with the model nick. */
void
test_ncm_mset_params_pretty_print (void)
{
  NcmModel *a   = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmMSet *mset = ncm_mset_new (a, NULL, NULL);
  FILE *out     = tmpfile ();
  gchar line[256];

  ncm_mset_params_pretty_print (mset, out, "a header");
  rewind (out);

  g_assert_nonnull (fgets (line, sizeof (line), out));
  g_assert_cmpstr (line, ==, "# a header\n");
  g_assert_nonnull (fgets (line, sizeof (line), out));
  g_assert_true (g_str_has_prefix (line, ncm_model_nick (a)));

  fclose (out);
  ncm_mset_free (mset);
  ncm_model_free (a);
}

/* Free parameters by index, by name and by full name, with two stacked models sharing
 * parameter names. */
void
test_ncm_mset_fparam_lookups (void)
{
  NcDistance *dist    = nc_distance_new (3.0);
  NcHaloPosition *hp0 = nc_halo_position_new (dist);
  NcHaloPosition *hp1 = nc_halo_position_new (dist);
  NcmMSet *mset       = ncm_mset_empty_new ();
  const NcmModelID id = ncm_model_id (NCM_MODEL (hp0));
  GError *error       = NULL;
  const NcmMSetPIndex *pi;
  NcmMSetPIndex *pif;
  guint ra;

  ncm_mset_push (mset, NCM_MODEL (hp0), NULL);
  ncm_mset_push (mset, NCM_MODEL (hp1), NULL);
  g_assert_true (ncm_model_orig_param_index_from_name (NCM_MODEL (hp0), "ra", &ra));

  ncm_model_param_set_ftype (NCM_MODEL (hp0), ra, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (NCM_MODEL (hp1), ra, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (mset);
  g_assert_cmpuint (ncm_mset_fparams_len (mset), ==, 2);
  g_assert_cmpuint (ncm_mset_fparam_len (mset), ==, 2);

  ncm_mset_fparam_set (mset, 1, 12.5);
  g_assert_cmpfloat (ncm_mset_fparam_get (mset, 1), ==, 12.5);
  g_assert_cmpfloat (ncm_model_param_get (NCM_MODEL (hp1), ra), ==, 12.5);
  g_assert_cmpstr (ncm_mset_fparam_full_name (mset, 1), ==, "NcHaloPosition:01:ra");

  g_assert_null (ncm_mset_fparam_get_pi_by_name (mset, "ra", &error));
  g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_PARAM_NAME_AMBIGUOUS);
  g_clear_error (&error);

  pi = ncm_mset_fparam_get_pi_by_name (mset, "NcHaloPosition:01:ra", &error);
  g_assert_no_error (error);
  g_assert_cmpint (pi->mid, ==, id + 1);
  g_assert_null (ncm_mset_fparam_get_pi_by_name (mset, "nope", &error));
  g_assert_no_error (error);

  pif = ncm_mset_param_get_by_full_name (mset, "NcHaloPosition:1:ra", &error);
  g_assert_no_error (error);
  g_assert_cmpint (pif->mid, ==, id + 1);
  g_assert_cmpuint (pif->pid, ==, ra);
  ncm_mset_pindex_free (pif);

  pif = ncm_mset_param_get_by_full_name (mset, "NcHaloPosition:0", &error);
  g_assert_no_error (error);
  g_assert_cmpint (pif->mid, ==, id);
  g_assert_cmpuint (pif->pid, ==, 0);
  ncm_mset_pindex_free (pif);

  g_assert_null (ncm_mset_param_get_by_full_name (mset, "NcHaloPosition:5:ra", &error));
  g_assert_no_error (error);
  g_assert_null (ncm_mset_param_get_by_full_name (mset, "NcHaloPosition:nope", &error));
  g_assert_no_error (error);
  g_assert_null (ncm_mset_param_get_by_full_name (mset, "NoSuchModel:ra", &error));
  g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_NAMESPACE_NOT_FOUND);
  g_clear_error (&error);

  ncm_mset_free (mset);
  nc_halo_position_free (hp0);
  nc_halo_position_free (hp1);
  nc_distance_free (dist);
}

void
test_ncm_mset_fparam_get_fpi_range (void)
{
  g_test_trap_subprocess ("/ncm/mset/fparam/get_fpi_range/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*parameter 2 is out of range*");
}

void
test_ncm_mset_fparam_get_fpi_range_subprocess (void)
{
  NcmModel *a   = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmMSet *mset = ncm_mset_new (a, NULL, NULL);

  ncm_mset_prepare_fparam_map (mset);
  ncm_mset_fparam_get_fpi (mset, ncm_model_mvnd_id (), 2);
}

/* A map set in its own order survives save and load; a file without the fmap key, the
 * format before it, gets the map prepared from the fit types. */
void
test_ncm_mset_saveload_fmap (void)
{
  NcmModel *a         = NCM_MODEL (ncm_model_mvnd_new (3));
  NcmMSet *mset       = ncm_mset_new (a, NULL, NULL);
  NcmSerialize *ser   = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  const gchar *fmap[] = {"NcmModelMVND:mu_2", "NcmModelMVND:mu_0", NULL};
  gchar *filename     = NULL;
  gint fd             = g_file_open_tmp ("test_ncm_mset_fmap_XXXXXX.mset", &filename, NULL);
  NcmMSet *loaded;

  g_assert_cmpint (fd, >=, 0);
  g_close (fd, NULL);

  ncm_mset_set_fmap (mset, fmap, TRUE, NULL);

  {
    NcmSerialize *ser_save = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);

    ncm_mset_save (mset, ser_save, filename, TRUE, NULL);
    ncm_serialize_free (ser_save);
  }

  loaded = ncm_mset_load (filename, ser, NULL);
  g_assert_true (ncm_mset_fparam_map_valid (loaded));
  g_assert_cmpuint (ncm_mset_fparams_len (loaded), ==, 2);
  g_assert_cmpstr (ncm_mset_fparam_full_name (loaded, 0), ==, "NcmModelMVND:mu_2");
  g_assert_cmpstr (ncm_mset_fparam_full_name (loaded, 1), ==, "NcmModelMVND:mu_0");
  ncm_mset_free (loaded);

  {
    GKeyFile *kf = g_key_file_new ();

    g_assert_true (g_key_file_load_from_file (kf, filename, G_KEY_FILE_NONE, NULL));
    g_assert_true (g_key_file_remove_key (kf, "NcmMSet", "fmap", NULL));
    g_assert_true (g_key_file_save_to_file (kf, filename, NULL));
    g_key_file_free (kf);
  }

  loaded = ncm_mset_load (filename, ser, NULL);
  g_assert_true (ncm_mset_fparam_map_valid (loaded));
  g_assert_cmpuint (ncm_mset_fparams_len (loaded), ==, 2);
  g_assert_cmpstr (ncm_mset_fparam_full_name (loaded, 0), ==, "NcmModelMVND:mu_0");
  g_assert_cmpstr (ncm_mset_fparam_full_name (loaded, 1), ==, "NcmModelMVND:mu_2");
  ncm_mset_free (loaded);

  g_unlink (filename);
  g_free (filename);
  ncm_serialize_free (ser);
  ncm_mset_free (mset);
  ncm_model_free (a);
}

/* Loading leaves no submodel names in the serializer, so a second load with the same
 * serializer works (one without NCM_SERIALIZE_OPT_AUTOSAVE_SER, which would keep the
 * names of the models it builds). */
void
test_ncm_mset_load_twice_same_ser (void)
{
  NcHIReion *reion  = NC_HIREION (nc_hireion_camb_new ());
  NcHIPrim *prim    = NC_HIPRIM (nc_hiprim_power_law_new ());
  NcmModel *cosmo   = NCM_MODEL (nc_hicosmo_de_xcdm_new_full (reion, prim, NULL));
  NcmMSet *mset     = ncm_mset_new (cosmo, NULL, NULL);
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  gchar *filename   = NULL;
  gint fd           = g_file_open_tmp ("test_ncm_mset_twice_XXXXXX.mset", &filename, NULL);
  guint i;

  g_assert_cmpint (fd, >=, 0);
  g_close (fd, NULL);

  {
    NcmSerialize *ser_save = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);

    ncm_mset_save (mset, ser_save, filename, FALSE, NULL);
    ncm_serialize_free (ser_save);
  }

  for (i = 0; i < 2; i++)
  {
    GError *error   = NULL;
    NcmMSet *loaded = ncm_mset_load (filename, ser, &error);

    g_assert_no_error (error);
    g_assert_cmpuint (ncm_mset_nmodels (loaded), ==, ncm_mset_nmodels (mset));
    g_assert_nonnull (ncm_mset_peek (loaded, nc_hiprim_id ()));
    g_assert_false (ncm_serialize_contain_name (ser, "NcHIPrim"));
    ncm_mset_free (loaded);
  }

  g_unlink (filename);
  g_free (filename);
  ncm_serialize_free (ser);
  ncm_mset_free (mset);
  ncm_model_free (cosmo);
  nc_hiprim_free (prim);
  nc_hireion_free (reion);
}

static void
_test_ncm_mset_value_int (GValue *v, const gint i)
{
  g_value_init (v, G_TYPE_INT);
  g_value_set_int (v, i);
}

static void
_test_ncm_mset_value_str (GValue *v, const gchar *str)
{
  g_value_init (v, G_TYPE_STRING);
  g_value_set_string (v, str);
}

/* __getitem__ and __setitem__ by integer id, stack position included, and by name;
 * invalid stack positions are errors on both. */
void
test_ncm_mset_getsetitem (void)
{
  NcDistance *dist    = nc_distance_new (3.0);
  NcHaloPosition *h0  = nc_halo_position_new (dist);
  NcHaloPosition *h1  = nc_halo_position_new (dist);
  NcHaloPosition *h2  = nc_halo_position_new (dist);
  NcmModel *mvnd      = NCM_MODEL (ncm_model_mvnd_new (2));
  NcmMSet *mset       = ncm_mset_empty_new ();
  const NcmModelID id = ncm_model_id (NCM_MODEL (h0));
  const gchar *bad[]  = {"NcHaloPosition:", "NcHaloPosition:+1", "NcHaloPosition:1x", "NcHaloPosition:1000"};
  GError *error       = NULL;
  GValue v            = G_VALUE_INIT;
  guint i;

  _test_ncm_mset_value_int (&v, id);
  ncm_mset___setitem__ (mset, &v, NCM_MODEL (h0), &error);
  g_assert_no_error (error);
  g_value_unset (&v);

  _test_ncm_mset_value_int (&v, id + 2);
  ncm_mset___setitem__ (mset, &v, NCM_MODEL (h2), &error);
  g_assert_no_error (error);
  g_assert_true (ncm_mset___getitem__ (mset, &v, &error) == NCM_MODEL (h2));
  g_value_unset (&v);

  _test_ncm_mset_value_str (&v, "NcHaloPosition:01");
  ncm_mset___setitem__ (mset, &v, NCM_MODEL (h1), &error);
  g_assert_no_error (error);
  g_assert_true (ncm_mset___getitem__ (mset, &v, &error) == NCM_MODEL (h1));
  g_value_unset (&v);
  g_assert_true (ncm_mset_peek_pos (mset, id, 1) == NCM_MODEL (h1));
  g_assert_true (ncm_mset_peek_by_name (mset, "NcHaloPosition", NULL) == NCM_MODEL (h0));

  for (i = 0; i < G_N_ELEMENTS (bad); i++)
  {
    _test_ncm_mset_value_str (&v, bad[i]);

    ncm_mset___setitem__ (mset, &v, NCM_MODEL (h1), &error);
    g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_NAMESPACE_INVALID);
    g_clear_error (&error);

    g_assert_null (ncm_mset___getitem__ (mset, &v, &error));
    g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_NAMESPACE_INVALID);
    g_clear_error (&error);

    g_value_unset (&v);
  }

  /* The positions set before are untouched by the rejected names. */
  g_assert_true (ncm_mset_peek_pos (mset, id, 0) == NCM_MODEL (h0));
  g_assert_true (ncm_mset_peek_pos (mset, id, 1) == NCM_MODEL (h1));
  g_assert_true (ncm_mset_peek_pos (mset, id, 2) == NCM_MODEL (h2));

  _test_ncm_mset_value_int (&v, id + 1);
  ncm_mset___setitem__ (mset, &v, mvnd, &error);
  g_assert_error (error, NCM_MSET_ERROR, NCM_MSET_ERROR_MODEL_ID_MISMATCH);
  g_clear_error (&error);
  g_value_unset (&v);

  ncm_mset_free (mset);
  ncm_model_free (mvnd);
  nc_halo_position_free (h0);
  nc_halo_position_free (h1);
  nc_halo_position_free (h2);
  nc_distance_free (dist);
}

