/***************************************************************************
 *            mcat_analyze.c
 *
 *  Mon March 16 11:04:23 2015
 *  Copyright  2015  Sandro Dias Pinto Vitenti
 *  <sandro@iaoftware.com.br>
 ****************************************************************************/
/*
 * mcat_analyze.c
 *
 * Copyright (C) 2015 - Sandro Dias Pinto Vitenti
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


gint
main (gint argc, gchar *argv[])
{
  gchar **models = NULL;
  gchar *outfile = NULL;
  gboolean save_comments = FALSE;
  gboolean overwrite = FALSE;

  GError *error = NULL;
  GOptionContext *context;
  GOptionEntry entries[] =
  {
    { "model",     'm', 0, G_OPTION_ARG_STRING_ARRAY, &models,        "Model (also accepts a serialized version), repeat for several models.", NULL },
    { "out",       'o', 0, G_OPTION_ARG_FILENAME,     &outfile,       "Filename where the mset should be written", NULL},
    { "comments",  'c', 0, G_OPTION_ARG_NONE,         &save_comments, "Whether comments must be saved in the .mset file.", NULL},
    { "overwrite", 'w', 0, G_OPTION_ARG_NONE,         &overwrite,     "Whether it should overwrite an already existing .mset file.", NULL},
    { NULL }
  };

  ncm_cfg_init ();

  context = g_option_context_new ("- generate a NcmMSet catalog.");
  g_option_context_set_summary (context, "MSet generator");
  g_option_context_set_description (context, "Generate an .mset file containing all chosen models");

  g_option_context_add_main_entries (context, entries, NULL);

  if (!g_option_context_parse (context, &argc, &argv, &error))
  {
    g_print ("option parsing failed: %s\n", error->message);
    exit (1);
  }

  if (models == NULL)
  {
    g_print ("No model chosen, use --model/-m to choose at least one model.\n");
    exit (1);
  }
  else
  {
    NcmMSet *mset = ncm_mset_empty_new ();
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    guint nmodels = g_strv_length (models);
    guint i;
    for (i = 0; i < nmodels; i++)
    {
      NcmModel *model = NCM_MODEL (ncm_serialize_from_string (ser, models[i]));
      g_assert (NCM_IS_MODEL (model));
      ncm_mset_set (mset, model);
    }

    outfile = outfile != NULL ? outfile : "models.mset";
    if (g_file_test (outfile, G_FILE_TEST_EXISTS) && !overwrite)
    {
      g_print ("File `%s' already exists, choose --overwrite/-w to overwrite.\n", outfile);
      exit (1);
    }

    ncm_mset_save (mset, ser, outfile, save_comments);

    ncm_mset_clear (&mset);
    ncm_serialize_clear (&ser);
  }

  return 0;
}
