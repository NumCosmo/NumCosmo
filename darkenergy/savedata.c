/***************************************************************************
 *            savedata.c
 *
 *  Tue Apr 26 23:13:55 2011
 *  Copyright  2011  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include "de_options.h"
#include "data_cluster.h"

FILE *
nc_de_open_dataout_file (NcHICosmo *model, gchar *prefix, gchar **filename)
{
  FILE *f;
  time_t now = time (NULL);
  GDate *date = g_date_new ();
  struct tm *time_ptr = localtime (&now);
  gchar date_string[500];
  g_date_set_time_t (date, now);
  g_assert (*filename == NULL);

  g_date_strftime(date_string, 500, "%G_%m_%d", date);

  *filename = g_strdup_printf ("%s_%s_%s_%02d_%02d.dat", prefix, ncm_model_name (NCM_MODEL (model)), date_string, time_ptr->tm_hour, time_ptr->tm_min);

  f = fopen (*filename, "w");

  fprintf (f, "# NumCosmo Version -- "NUMCOSMO_VERSION"\n");

  return f;
}
