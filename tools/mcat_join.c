/***************************************************************************
 *            mcat_join.c
 *
 *  Thu October 27 10:31:13 2016
 *  Copyright  2016  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
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
  gchar **cat_filename = NULL;
  gchar **burnins      = NULL;
  gchar *out           = NULL;
  
  GError *error = NULL;
  GOptionContext *context;
  GOptionEntry entries[] =
  {
    { "catalog",        'c', 0, G_OPTION_ARG_STRING_ARRAY, &cat_filename,   "Input catalog filename.", NULL },
    { "burnin",         'b', 0, G_OPTION_ARG_STRING_ARRAY, &burnins,        "Burnin for the input catalogs.", NULL },
    { "out",            'o', 0, G_OPTION_ARG_STRING,       &out,            "Output catalog.", NULL },
    { NULL }
  };

  ncm_cfg_init_full_ptr (&argc, &argv);
  
  context = g_option_context_new ("- join different compatible catalogs in a single one.");
  g_option_context_set_summary (context, "catalog join");
  g_option_context_set_description (context, "CJ Description <FIXME>");

  g_option_context_add_main_entries (context, entries, NULL);

  if (!g_option_context_parse (context, &argc, &argv, &error))
  {
    g_print ("option parsing failed: %s\n", error->message);
    exit (1);
  }

  if (cat_filename == NULL)
  {
    g_print ("No input catalogs, use --catalog/-c.\n");
    exit (1);
  }
  else
  {
    const guint nmcats     = g_strv_length (cat_filename);
    const guint nburnins   = (burnins != NULL) ? g_strv_length (burnins) : 0;

    if (nburnins > nmcats)
      g_warning ("mcat_join: more burnins than catalogs nburnins %u > nmcats %u!", nburnins, nmcats);
    
    if (nmcats == 1)
    {
      g_print ("A single input catalog was passed, nothing to do!.\n");
    }
    else
    {
      GPtrArray *mcats_array   = g_ptr_array_new ();
      GArray *max_time_array   = g_array_new (TRUE, TRUE, sizeof (guint));
      NcmMSetCatalog *mcat_out = ncm_mset_catalog_new_from_file_ro (cat_filename[0], 0);
      const guint nchains      = ncm_mset_catalog_nchains (mcat_out);
      NcmMSet *mset            = ncm_mset_catalog_get_mset (mcat_out);
      guint mmax_time          = 0;
      guint i, t;
      
      g_ptr_array_set_free_func (mcats_array, (GDestroyNotify)ncm_mset_catalog_free);

      ncm_mset_catalog_set_file (mcat_out, NULL);
      ncm_mset_catalog_reset (mcat_out);
      
      for (i = 0; i < nmcats; i++)
      {
        const glong burnin   = (i < nburnins) ? atol (burnins[i]) : 0;
        NcmMSetCatalog *mcat = ncm_mset_catalog_new_from_file_ro (cat_filename[i], burnin);
        NcmMSet *mset_i      = ncm_mset_catalog_get_mset (mcat);
        guint max_time_i     = ncm_mset_catalog_max_time (mcat);

        g_assert_cmpuint (ncm_mset_catalog_nchains (mcat), ==, nchains);
        g_assert_cmpuint (ncm_mset_cmp_all (mset, mset_i), ==, 0);

        if (burnin % nchains != 0)
          g_error ("Burnin[%u] not multiple of nchains %u.", i, nchains);

        mmax_time = GSL_MAX (mmax_time, max_time_i);
        
        ncm_mset_free (mset_i);
        g_ptr_array_add (mcats_array, mcat);
        g_array_append_val (max_time_array, max_time_i);
      }

      for (t = 0; t < mmax_time; t++)
      {
        const guint offset = t * nchains;
        for (i = 0; i < nmcats; i++)
        {
          if (t < g_array_index (max_time_array, guint, i))
          {
            NcmMSetCatalog *mcat_i = g_ptr_array_index (mcats_array, i);
            guint wi;
            for (wi = 0; wi < nchains; wi++)
            {
              NcmVector *row = ncm_mset_catalog_peek_row (mcat_i, offset + wi);
              ncm_mset_catalog_add_from_vector (mcat_out, row);
            }
          }
        }
      }

      if (out == NULL)
      {
        gchar *out_default = g_strdup_printf ("joined_mcat.fits");
        guint i = 0;
        
        while (g_file_test (out_default, G_FILE_TEST_EXISTS))
        {
          g_free (out_default);
          out_default = g_strdup_printf ("joined_mcat_%d.fits", i++);
        }
        
        ncm_mset_catalog_set_file (mcat_out, out_default);
        g_free (out_default);
      }
      else
      {
        if (g_file_test (out, G_FILE_TEST_EXISTS))
        {
          g_warning ("out file already exists, renaming.");
          gchar *out_base = ncm_util_basename_fits (out);
          gchar *out_ren  = g_strdup_printf ("%s_ren.fits", out_base);
          guint i = 0;
          
          while (g_file_test (out_ren, G_FILE_TEST_EXISTS))
          {
            g_free (out_ren);
            out_ren = g_strdup_printf ("%s_ren%d.fits", out_base, i++);
          }

          ncm_mset_catalog_set_file (mcat_out, out_ren);
          
          g_free (out_base);
          g_free (out_ren);
        }
        else
          ncm_mset_catalog_set_file (mcat_out, out);
      }

      g_ptr_array_unref (mcats_array);
      g_array_unref (max_time_array);
      ncm_mset_free (mset);
    }
  }

  return 0;
}

