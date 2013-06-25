/***************************************************************************
 *            nc_data_snia.c
 *
 *  Mon December 10 00:20:48 2012
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

/**
 * SECTION:nc_data_snia
 * @title: SN Ia Data
 * @short_description: Helper function for obtaining Supernovae Ia data
 * 
 * FIXME
 * 
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc_data_snia.h"

#include "nc_data_dist_mu.h"
#include "nc_data_snia_cov.h"
#include "nc_snia_dist_cov.h"
#include "math/ncm_cfg.h"

#include <gio/gio.h>

static const gchar *_nc_data_snia_cats[] = {
  "snls3_conley_2011_sys_stat.fits",
  "snls3_conley_2011_stat_only.fits",
};

/**
 * nc_data_snia_get_catalog:
 * @id: FIXME
 * 
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
gchar * 
nc_data_snia_get_catalog (gchar *id)
{
  const GEnumValue *snia_id = 
    ncm_cfg_get_enum_by_id_name_nick (NC_TYPE_DATA_SNIA_ID, id);
  if (snia_id == NULL)
    g_error ("nc_data_snia_get_catalog: Cannot find id ``%s'' in catalogs.", id);
  return nc_data_snia_get_catalog_by_id (snia_id->value);
}

/**
 * nc_data_snia_get_catalog_by_id:
 * @id: FIXME
 * 
 * FIXME
 * 
 * Returns: (transfer full): FIXME
 */
gchar * 
nc_data_snia_get_catalog_by_id (NcDataSNIAId id)
{
  g_assert (id <= NC_DATA_SNIA_COV_END && id >= NC_DATA_SNIA_COV_START);
  {
    gint i = id - NC_DATA_SNIA_COV_START;
    gchar *full_filename = nc_data_snia_get_fits (_nc_data_snia_cats[i], FALSE);
    return full_filename;
  }
}

void 
nc_data_snia_load_cat (NcSNIADistCov *dcov, NcDataSNIAId id)
{
  gchar *full_filename = nc_data_snia_get_catalog_by_id (id);
#ifdef NUMCOSMO_HAVE_CFITSIO
  nc_snia_dist_cov_load (dcov, full_filename);
#else
  g_error ("nc_data_snia_load_cat: catalog load not available, recompile "PACKAGE_NAME" with cfitsio support.");
#endif /* NUMCOSMO_HAVE_CFITSIO */
  g_free (full_filename);
}

void 
_nc_data_snia_copy_prog (goffset current_num_bytes, goffset total_num_bytes, gpointer user_data)
{
  gint *old_prog = (gint *) user_data;
  gint prog = (100 * current_num_bytes) / total_num_bytes;
  if (prog > *old_prog)
  {
    ncm_message ("# % 3d%%\r", prog);
    *old_prog = prog;
  }
}

gchar *
nc_data_snia_get_fits (const gchar *filename, gboolean check_size)
{
  gchar *full_filename = ncm_cfg_get_fullpath (filename);
  gchar *url_str = g_strdup_printf ("http://download.savannah.gnu.org/releases/numcosmo/%s", filename);
  GFile *local  = g_file_new_for_path (full_filename);
  GFile *remote = g_file_new_for_uri (url_str);
  GError *error = NULL;
  gint prog = 0;
  gboolean download = FALSE;

  if (g_file_test (full_filename, G_FILE_TEST_EXISTS))
  {
    if (check_size)
    {
      GFileInfo *local_info, *remote_info;
      local_info = g_file_query_info (local, G_FILE_ATTRIBUTE_STANDARD_SIZE, 
                                      G_FILE_QUERY_INFO_NONE, 
                                      NULL, &error);
      if (local_info == NULL)
        g_error ("nc_data_snia_get_fits: cannot get info for %s: %s.", full_filename, error->message);

      remote_info = g_file_query_info (remote, G_FILE_ATTRIBUTE_STANDARD_SIZE, 
                                       G_FILE_QUERY_INFO_NONE, 
                                       NULL, &error);
      if (remote_info == NULL)
        g_error ("nc_data_snia_get_fits: cannot get info for %s: %s."
                 " To use this catalog, download the file from the url and copy "
                 "to ~/.numcosmo directory.", url_str, error->message);

      if (g_file_info_get_size (local_info) != g_file_info_get_size (remote_info))
        download = TRUE;
    }
  }
  else
    download = TRUE;

  if (download)
  {
    ncm_message ("# Downloading file [%s]...\n", url_str);
    if (!g_file_copy (remote, local, G_FILE_COPY_OVERWRITE, NULL, 
                      &_nc_data_snia_copy_prog, &prog, &error))
      g_error ("nc_data_snia_get_fits: cannot get fits file from %s: %s."
               " To use this catalog, download the file from the url and copy "
               "to ~/.numcosmo directory.", 
               url_str, error->message);
  }
  g_free (url_str);
  
  g_object_unref (local);
  g_object_unref (remote);
  
  return full_filename;
}
