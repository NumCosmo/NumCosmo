/***************************************************************************
 *            nc_data_download.c
 *
 *  Sun August 30 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
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

/*
 * One place where downloading a data file is done correctly, because the
 * destination is shared by every NumCosmo process on the machine and a test
 * suite under pytest-xdist runs several at once. See nc_data_download_priv.h.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "nc/data/nc_data_download_priv.h"
#include "ncm/core/ncm_cfg.h"

#include <glib/gstdio.h>
#include <errno.h>
#include <unistd.h>

gboolean
_nc_data_download_lock (const gchar *path, gint max_wait_s, gchar **lockdir)
{
  gint waited = 0;
  gchar *dir  = g_path_get_dirname (path);

  /* The lock sits beside the target, whose directory may not exist yet on a
   * fresh install -- mkdir would then fail with ENOENT and the caller would
   * silently skip the download it was asked for. */
  if (g_mkdir_with_parents (dir, 0755) != 0)
  {
    g_free (dir);
    *lockdir = NULL;

    return FALSE;
  }

  g_free (dir);

  *lockdir = g_strdup_printf ("%s.lock", path);

  while (g_mkdir (*lockdir, 0755) != 0)
  {
    if (errno != EEXIST)
    {
      g_clear_pointer (lockdir, g_free);

      return FALSE;
    }

    /* Whoever holds it may be finishing the very download we want. */
    if (g_file_test (path, G_FILE_TEST_EXISTS))
    {
      g_clear_pointer (lockdir, g_free);

      return FALSE;
    }

    if (waited >= max_wait_s)
      break;  /* stale lock from a killed run: take it over */

    g_usleep (G_USEC_PER_SEC);
    waited++;
  }

  if (g_file_test (path, G_FILE_TEST_EXISTS))
  {
    _nc_data_download_unlock (*lockdir);
    *lockdir = NULL;

    return FALSE;
  }

  return TRUE;
}

void
_nc_data_download_unlock (gchar *lockdir)
{
  if (lockdir != NULL)
  {
    g_rmdir (lockdir);
    g_free (lockdir);
  }
}

void
_nc_data_download_file (const gchar *url, const gchar *dest, const gchar *what)
{
  gchar *tmp    = g_strdup_printf ("%s.%d.part", dest, (gint) getpid ());
  GError *error = NULL;
  gint status   = 0;

  ncm_message ("# Downloading %s from [%s]...\n", what, url);

  {
    gchar *cmd[] = { "wget", "--tries=3", "--timeout=30", "-O", tmp, (gchar *) url, NULL };

    /* g_spawn_sync() returning TRUE only says the child started; without the
     * exit status a 404, a DNS failure or a full disk all looked like success,
     * and the first symptom was a corrupt file read much later. */
    if (!g_spawn_sync (NULL, cmd, NULL,
                       G_SPAWN_SEARCH_PATH | G_SPAWN_STDOUT_TO_DEV_NULL | G_SPAWN_STDERR_TO_DEV_NULL,
                       NULL, NULL, NULL, NULL, &status, &error))
    {
      g_unlink (tmp);
      g_error ("_nc_data_download_file: cannot run wget for %s. Error: %s. "
               "Please download %s by hand and place it at %s.",
               what, error->message, url, dest);
    }

    if (status != 0)
    {
      g_unlink (tmp);
      g_error ("_nc_data_download_file: wget failed (status %d) fetching %s. "
               "Please download %s by hand and place it at %s.",
               status, what, url, dest);
    }
  }

  /* Rename last: until this point no other process can see a partial file. */
  if (g_rename (tmp, dest) != 0)
  {
    g_unlink (tmp);
    g_error ("_nc_data_download_file: cannot move %s into place at %s.", what, dest);
  }

  g_free (tmp);
}

