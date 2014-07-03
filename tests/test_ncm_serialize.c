/***************************************************************************
 *            test_ncm_object_serialization.c
 *
 *  Wed June 20 21:56:09 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti & Mariana Penna Lima
 *  <sandro@isoftware.com.br> <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

void test_ncm_serialize_global_from_string_plain (void);
void test_ncm_serialize_global_from_string_params (void);
void test_ncm_serialize_global_from_string_nest (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init ();
#if !((GLIB_MAJOR_VERSION == 2) && (GLIB_MINOR_VERSION < 30))
  g_test_add_func ("/numcosmo/ncm_serialize/from_string/plain", &test_ncm_serialize_global_from_string_plain);
  g_test_add_func ("/numcosmo/ncm_serialize/from_string/params", &test_ncm_serialize_global_from_string_params);
  g_test_add_func ("/numcosmo/ncm_serialize/from_string/nest", &test_ncm_serialize_global_from_string_nest);
#endif /* has glib >= 2.30 */
  g_test_run ();
}

void
test_ncm_serialize_global_from_string_plain ()
{
  GObject *obj = ncm_serialize_global_from_string ("NcHICosmoLCDM");
  NcHICosmo *hic = NC_HICOSMO (obj);
  gchar *obj_ser = ncm_serialize_global_to_string (obj, TRUE);
  GObject *obj_new = ncm_serialize_global_from_string (obj_ser);
  gchar *obj_new_ser = ncm_serialize_global_to_string (obj_new, TRUE);

  g_assert (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);

  obj_ser = ncm_serialize_global_to_string (obj, FALSE);
  obj_new_ser = ncm_serialize_global_to_string (obj_new, FALSE);
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);

  NCM_TEST_FREE (nc_hicosmo_free, hic);
}

void
test_ncm_serialize_global_from_string_params ()
{
  NcmModel *m = NCM_MODEL (ncm_serialize_global_from_string ("NcHICosmoLCDM{'H0':<12.3>,'Omegac':<0.2>}"));
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_H0), ==, 12.3);
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_OMEGA_C), ==, 0.2);

  NCM_TEST_FREE (ncm_model_free, m);

  m = NCM_MODEL (ncm_serialize_global_from_string ("{'NcHICosmoLCDM', {'H0':<12.3>,'Omegac':<0.2>}}"));
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_H0), ==, 12.3);
  ncm_assert_cmpdouble (ncm_model_param_get (m, NC_HICOSMO_DE_OMEGA_C), ==, 0.2);

  NCM_TEST_FREE (ncm_model_free, m);
}

void
test_ncm_serialize_global_from_string_nest ()
{
  GObject *obj = ncm_serialize_global_from_string ("NcMatterVar{'strategy':<2>,'window':<{'NcWindowTophat',@a{sv} {}}>,'transfer':<{'NcTransferFuncEH',@a{sv} {}}>}");
  NcMatterVar *mv = NC_MATTER_VAR (obj);
  gchar *obj_ser = ncm_serialize_global_to_string (obj, TRUE);
  GObject *obj_new = ncm_serialize_global_from_string (obj_ser);
  gchar *obj_new_ser = ncm_serialize_global_to_string (obj_new, TRUE);

  g_assert (G_OBJECT_TYPE (obj) == G_OBJECT_TYPE (obj_new));
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);

  obj_ser = ncm_serialize_global_to_string (obj, FALSE);
  obj_new_ser = ncm_serialize_global_to_string (obj_new, FALSE);
  g_assert_cmpstr (obj_ser, ==, obj_new_ser);
  g_free (obj_ser);
  g_free (obj_new_ser);
  g_object_unref (obj_new);

  g_assert (NC_IS_WINDOW_TOPHAT (mv->wp));
  g_assert (NC_IS_TRANSFER_FUNC_EH (mv->tf));

  NCM_TEST_FREE (nc_matter_var_free, mv);
}

