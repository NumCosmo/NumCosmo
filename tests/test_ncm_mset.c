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

#include <math.h>
#include <glib.h>
#include <glib-object.h>

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
void test_ncm_mset_fparams_validate_all (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_dup (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_shallow_copy (TestNcmMSet *test, gconstpointer pdata);
void test_ncm_mset_saveload (TestNcmMSet *test, gconstpointer pdata);

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

  g_test_add ("/ncm/mset/fparams/validate_all", TestNcmMSet, NULL,
              &test_ncm_mset_new,
              &test_ncm_mset_fparams_validate_all,
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

    ncm_mset_set (test->mset, NCM_MODEL (cosmo));
    ncm_mset_set (test->mset, NCM_MODEL (cosmo));

    g_ptr_array_add (test->ma, cosmo);
    g_array_append_val (test->ma_destroyed, f);

    nc_hicosmo_free (NC_HICOSMO (cosmo));
  }

  {
    NcDistance *dist    = nc_distance_new (5.0);
    NcSNIADistCov *snia = nc_snia_dist_cov_new (dist, 4);
    gboolean f          = FALSE;

    ncm_mset_set (test->mset, NCM_MODEL (snia));

    g_ptr_array_add (test->ma, snia);
    g_array_append_val (test->ma_destroyed, f);

    nc_snia_dist_cov_free (snia);
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

  g_assert (mass != NULL);
  g_assert (NC_IS_CLUSTER_MASS (mass));

  ncm_mset_set (test->mset, NCM_MODEL (mass));

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

  g_assert (mass != NULL);
  g_assert (NC_IS_CLUSTER_MASS (mass));

  ncm_mset_set_pos (test->mset, NCM_MODEL (mass), 5);

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

  g_assert (mass != NULL);
  g_assert (NC_IS_CLUSTER_MASS (mass));

  ncm_mset_push (test->mset, NCM_MODEL (mass));
  ncm_mset_push (test->mset, NCM_MODEL (mass));
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

  g_assert (mass != NULL);
  g_assert (NC_IS_CLUSTER_MASS (mass));
  g_assert (benson != NULL);
  g_assert (NC_IS_CLUSTER_MASS (benson));

  ncm_mset_set_pos (test->mset, NCM_MODEL (mass), 10);

  ncm_mset_push (test->mset, NCM_MODEL (benson));
  ncm_mset_push (test->mset, NCM_MODEL (benson));
  ncm_mset_push (test->mset, NCM_MODEL (mass));
  ncm_mset_push (test->mset, NCM_MODEL (benson));

  ncm_mset_set_pos (test->mset, NCM_MODEL (mass), 1);

  g_ptr_array_add (test->ma, mass);
  g_array_append_val (test->ma_destroyed, f);

  g_ptr_array_add (test->ma, benson);
  g_array_append_val (test->ma_destroyed, f);

  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);

  g_assert_cmpuint (ncm_mset_total_len (test->mset), ==, ncm_mset_fparam_len (test->mset));

  ncm_mset_push (test->mset, NCM_MODEL (mass));
  g_assert_cmpuint (ncm_mset_total_len (test->mset), ==, ncm_mset_fparam_len (test->mset));

  ncm_mset_param_set_ftype (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 1), 0, NCM_PARAM_TYPE_FIXED);

  g_assert_cmpint (ncm_mset_param_get_ftype (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 1), 0), ==, NCM_PARAM_TYPE_FIXED);
  g_assert_cmpint (ncm_mset_param_get_ftype (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 2), 0), ==, NCM_PARAM_TYPE_FIXED);
  g_assert_cmpint (ncm_mset_param_get_ftype (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 4), 0), ==, NCM_PARAM_TYPE_FIXED);
  g_assert_cmpint (ncm_mset_param_get_ftype (test->mset, NCM_MSET_MID (nc_cluster_mass_id (), 10), 0), ==, NCM_PARAM_TYPE_FIXED);

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
test_ncm_mset_fparams_validate_all (TestNcmMSet *test, gconstpointer pdata)
{
  test_ncm_mset_fparams (test, pdata);

  ncm_mset_fparam_set (test->mset, 0, 1.0);
  ncm_mset_param_set_all_ftype (test->mset, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);

  {
    const gint ntests     = 100000;
    const gint fparam_len = ncm_mset_fparam_len (test->mset);
    NcmVector *theta      = ncm_vector_new (fparam_len);
    gint i;

    ncm_mset_fparams_get_vector (test->mset, theta);

    for (i = 0; i < ntests; i++)
    {
      g_assert_true (ncm_mset_fparam_validate_all (test->mset, theta));
    }

    ncm_vector_free (theta);
  }
}

void
test_ncm_mset_dup (TestNcmMSet *test, gconstpointer pdata)
{
  NcClusterMass *mass = NC_CLUSTER_MASS (ncm_serialize_global_from_string ("NcClusterMassLnnormal"));
  gboolean f          = FALSE;

  g_assert (mass != NULL);
  g_assert (NC_IS_CLUSTER_MASS (mass));

  ncm_mset_push (test->mset, NCM_MODEL (mass));
  ncm_mset_push (test->mset, NCM_MODEL (mass));
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

    g_assert (benson != NULL);
    g_assert (NC_IS_CLUSTER_MASS (benson));

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

    ncm_mset_push (test->mset, NCM_MODEL (mass));
    g_assert_true (!ncm_mset_cmp (test->mset, mset_dup, FALSE));
    g_assert_true (!ncm_mset_cmp (test->mset, mset_dup, TRUE));

    /*g_ptr_array_add (test->ma, benson);*/
    /*g_array_append_val (test->ma_destroyed, f);*/

    ncm_mset_push (mset_dup, NCM_MODEL (benson));
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

  g_assert (mass != NULL);
  g_assert (NC_IS_CLUSTER_MASS (mass));

  ncm_mset_push (test->mset, NCM_MODEL (mass));
  ncm_mset_push (test->mset, NCM_MODEL (mass));
  g_ptr_array_add (test->ma, mass);
  g_array_append_val (test->ma_destroyed, f);

  {
    NcmMSet *mset_dup = ncm_mset_shallow_copy (test->mset);

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

  g_assert (mass != NULL);
  g_assert (NC_IS_CLUSTER_MASS (mass));
  g_assert (benson != NULL);
  g_assert (NC_IS_CLUSTER_MASS (benson));

  ncm_mset_push (test->mset, NCM_MODEL (mass));
  ncm_mset_push (test->mset, NCM_MODEL (mass));
  ncm_mset_push (test->mset, NCM_MODEL (benson));

  g_ptr_array_add (test->ma, mass);
  g_array_append_val (test->ma_destroyed, f);

  g_ptr_array_add (test->ma, benson);
  g_array_append_val (test->ma_destroyed, f);

  {
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);

    ncm_mset_save (test->mset, ser, "test_ncm_mset_saved.mset", TRUE);

    NcmMSet *mset_dup = ncm_mset_load ("test_ncm_mset_saved.mset", ser);

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

    ncm_mset_push (test->mset, NCM_MODEL (mass));
    g_assert_true (!ncm_mset_cmp (test->mset, mset_dup, FALSE));
    g_assert_true (!ncm_mset_cmp (test->mset, mset_dup, TRUE));

    /*g_ptr_array_add (test->ma, benson);*/
    /*g_array_append_val (test->ma_destroyed, f);*/

    ncm_mset_push (mset_dup, NCM_MODEL (benson));
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
  }

  nc_cluster_mass_free (benson);
  nc_cluster_mass_free (mass);
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

  ncm_mset_push (test->mset, NCM_MODEL (cosmo));

  nc_hicosmo_free (NC_HICOSMO (cosmo));
}

