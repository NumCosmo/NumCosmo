/***************************************************************************
 *            ncm_obj_array.c
 *
 *  Fri March 23 14:26:42 2012
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

#include <math.h>
#include <glib.h>
#include <glib-object.h>

typedef struct _TestNcmObjArray
{
  NcmObjArray *oa;
   guint ntests;
} TestNcmObjArray;

void test_ncm_obj_array_new (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_free (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_add (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_saveload (TestNcmObjArray *test, gconstpointer pdata);

void test_ncm_obj_array_traps (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_invalid_add (TestNcmObjArray *test, gconstpointer pdata);
void test_ncm_obj_array_invalid_set (TestNcmObjArray *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/obj_array/add", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_add,
              &test_ncm_obj_array_free);

  g_test_add ("/ncm/obj_array/saveload", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_saveload,
              &test_ncm_obj_array_free);

  g_test_add ("/ncm/obj_array/traps", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_traps,
              &test_ncm_obj_array_free);
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 38))
  g_test_add ("/ncm/obj_array/invalid/add/subprocess", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_invalid_add,
              &test_ncm_obj_array_free);

  g_test_add ("/ncm/obj_array/invalid/set/subprocess", TestNcmObjArray, NULL,
              &test_ncm_obj_array_new,
              &test_ncm_obj_array_invalid_set,
              &test_ncm_obj_array_free);
#endif
  g_test_run ();
}

void
test_ncm_obj_array_new (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa;

  test->ntests = 1000;
  oa = test->oa = ncm_obj_array_new ();

  g_assert (oa != NULL);
}

void
test_ncm_obj_array_free (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa = test->oa;
  ncm_obj_array_unref (oa);
}

void
test_ncm_obj_array_add (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa = test->oa;

  g_assert (oa != NULL);

  {
    NcmVector *v     = ncm_vector_new (10);
    NcmMatrix *m     = ncm_matrix_new (5, 8);
    NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_lcdm_new ());
    NcDistance *dist = nc_distance_new (1.0);
    NcDataBaoA *baoa = nc_data_bao_a_new_from_id (dist, 0);

    ncm_vector_set_all (v, 1.2);
    ncm_matrix_set_all (m, 2.0);
    
    ncm_obj_array_add (oa, G_OBJECT (v));
    ncm_obj_array_add (oa, G_OBJECT (m));
    ncm_obj_array_add (oa, G_OBJECT (cosmo));

    ncm_obj_array_set (oa, 0, G_OBJECT (dist));
    ncm_obj_array_add (oa, G_OBJECT (dist));

    ncm_obj_array_add (oa, G_OBJECT (baoa));
    
    ncm_vector_free (v);
    ncm_matrix_free (m);
    nc_hicosmo_free (cosmo);
    nc_distance_free (dist);
  }
}

void
test_ncm_obj_array_saveload (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa   = test->oa;
  NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);

  test_ncm_obj_array_add (test, pdata);

  ncm_obj_array_save (oa, ser, "test_ncm_obj_array_saved.oa", TRUE);

  {
    NcmObjArray *oa_load = ncm_obj_array_load ("test_ncm_obj_array_saved.oa", ser);

    ncm_obj_array_unref (oa_load);
  }

  ncm_serialize_clear (&ser);
}

void
test_ncm_obj_array_traps (TestNcmObjArray *test, gconstpointer pdata)
{
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 38))
  g_test_trap_subprocess ("/ncm/obj_array/invalid/add/subprocess", 0, 0);
  g_test_trap_assert_failed ();

  g_test_trap_subprocess ("/ncm/obj_array/invalid/set/subprocess", 0, 0);
  g_test_trap_assert_failed ();
#endif
}

void
test_ncm_obj_array_invalid_add (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa = test->oa;
  ncm_obj_array_add (oa, NULL);
}

void
test_ncm_obj_array_invalid_set (TestNcmObjArray *test, gconstpointer pdata)
{
  NcmObjArray *oa = test->oa;
  ncm_obj_array_add (oa, NULL);

  test_ncm_obj_array_add (test, pdata);

  {
    NcmVector *v     = ncm_vector_new (10);

    ncm_obj_array_set (oa, 10, G_OBJECT (v));
    
    ncm_vector_free (v);
  }
}

