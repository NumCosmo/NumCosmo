/***************************************************************************
 *            nc_data_download_priv.h
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

#ifndef _NC_DATA_DOWNLOAD_PRIV_H_
#define _NC_DATA_DOWNLOAD_PRIV_H_

#include <glib.h>

G_BEGIN_DECLS

/*
 * Fetching a data file is not a private act: the destination is a directory
 * shared by every NumCosmo process on the machine, and a test suite under
 * pytest-xdist runs several at once. Downloading straight to the final name
 * lets one process read what another is still writing -- a truncated FITS, or
 * a tar unpacked from an archive that is not finished -- so these serialize on
 * a lock, fetch to a per-process temporary, and rename into place only once
 * the transfer has been checked.
 */

/*
 * Serializes on @lockpath and reports whether there is work to do. Returns
 * TRUE when the caller owns the lock and must pass @lockdir to
 * _nc_data_download_unlock(); FALSE when @readypath already exists, which is
 * the common case once one of N workers has finished.
 *
 * The two are separate because "done" is not always "the destination exists".
 * An extracted tree is only complete when a marker inside it says so -- a
 * half-unpacked one has its first files and not its last -- and that marker
 * cannot double as the lock, since the tree it lives in gets replaced.
 *
 * A lock *directory* rather than a lock file: mkdir() is atomic everywhere
 * this runs, NFS included. A lock left by a killed process is taken over
 * after @max_wait_s rather than waited on forever.
 */
gboolean _nc_data_download_lock (const gchar *lockpath, const gchar *readypath, gint max_wait_s, gchar **lockdir);

void _nc_data_download_unlock (gchar *lockdir);

/*
 * Fetches @url to @dest. **Call with the lock held.**
 *
 * The transfer goes to a per-process temporary beside @dest and is renamed
 * over it only once wget has reported success, so no other process can observe
 * a partial file and an interrupted run leaves nothing that looks like good
 * data. @what names the data in the error message; failure is fatal, and says
 * the URL and directory so the file can be fetched by hand.
 */
void _nc_data_download_file (const gchar *url, const gchar *dest, const gchar *what);

G_END_DECLS

#endif /* _NC_DATA_DOWNLOAD_PRIV_H_ */

