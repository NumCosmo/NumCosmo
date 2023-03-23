/***************************************************************************
 *            mcat_join.c
 *
 *  Thu October 27 10:31:13 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * mcat_join.c
 *
 * Copyright (C) 2016 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <gsl/gsl_cdf.h>

gint
main (gint argc, gchar *argv[])
{
  gchar *cat_filename = NULL;
  gchar *out           = NULL;
  
  GError *error = NULL;
  GOptionContext *context;
  GOptionEntry entries[] =
  {
    { "catalog",        'c', 0, G_OPTION_ARG_STRING, &cat_filename,   "Input catalog filename.", NULL },
    { "out",            'o', 0, G_OPTION_ARG_STRING, &out,            "Output catalog.", NULL },
    { NULL }
  };

  ncm_cfg_init_full_ptr (&argc, &argv);
  
  context = g_option_context_new ("- trim out of bounds points.");
  g_option_context_set_summary (context, "catalog trim oob");
  g_option_context_set_description (context, "CJ Description <FIXME>");

  g_option_context_add_main_entries (context, entries, NULL);

  if (!g_option_context_parse (context, &argc, &argv, &error))
  {
    g_print ("option parsing failed: %s\n", error->message);
    exit (1);
  }

  if ((cat_filename == NULL) || (out == NULL))
  {
    g_print ("No input catalog and/or output file, use --catalog/-c and --out/-o.\n");
    exit (1);
  }
  else
  {    
    NcmMSetCatalog *mcat = ncm_mset_catalog_new_from_file_ro (cat_filename, 0);
    guint ndel = ncm_mset_catalog_trim_oob (mcat, out);

    g_message ("# Trimming `%s', %u points deleted, saving to `%s'.\n", cat_filename, ndel, out);
  }
  
  return 0;
}

