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
    { NULL, 0, 0, 0, NULL, NULL, NULL }
  };

  ncm_cfg_init_full_ptr (&argc, &argv);

  context = g_option_context_new ("- join different compatible catalogs in a single one.");
  g_option_context_set_summary (context, "catalog join");
  g_option_context_set_description (context, "Joins different compatible catalogs into a single one.");

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
    const guint nmcats   = g_strv_length (cat_filename);
    const guint nburnins = (burnins != NULL) ? g_strv_length (burnins) : 0;

    if (nburnins > nmcats)
      g_warning ("mcat_join: more burnins than catalogs nburnins %u > nmcats %u!", nburnins, nmcats);

    if (nmcats == 1)
    {
      g_print ("A single input catalog was passed, nothing to do!.\n");
    }
    else
    {
      NcmMSetCatalogStorageFits *first_fits = ncm_mset_catalog_storage_fits_new (cat_filename[0], TRUE, NULL);
      NcmMSetCatalogStorage *first          = NCM_MSET_CATALOG_STORAGE (first_fits);
      NcmMSetCatalogDesc *desc              = ncm_mset_catalog_storage_peek_desc (first, NULL);
      GPtrArray *mcats_array                = g_ptr_array_new ();
      GArray *max_time_array                = g_array_new (TRUE, TRUE, sizeof (guint));
      gchar *out_name                       = NULL;
      NcmMSetCatalogStorageFits *out_fits;
      NcmMSetCatalog *mcat_out;
      guint nchains;
      guint mmax_time = 0;
      guint i, t;

      if (desc == NULL)
        g_error ("mcat_join: file `%s' holds no catalog.", cat_filename[0]);

      nchains = ncm_mset_catalog_desc_nchains (desc);
      g_ptr_array_set_free_func (mcats_array, (GDestroyNotify) ncm_mset_catalog_free);

      /* Every input must describe the same catalog, with the parameters in the same order. */
      for (i = 0; i < nmcats; i++)
      {
        const glong burnin                      = (i < nburnins) ? atol (burnins[i]) : 0;
        NcmMSetCatalogStorageFits *storage_fits = ncm_mset_catalog_storage_fits_new (cat_filename[i], TRUE, NULL);
        NcmMSetCatalogStorage *storage          = NCM_MSET_CATALOG_STORAGE (storage_fits);
        NcmMSetCatalogDesc *desc_i              = ncm_mset_catalog_storage_peek_desc (storage, NULL);
        gchar *reason                           = NULL;
        NcmMSetCatalog *mcat;
        guint max_time_i;

        if (desc_i == NULL)
          g_error ("mcat_join: file `%s' holds no catalog.", cat_filename[i]);

        if (!ncm_mset_catalog_desc_is_compatible (desc, desc_i, &reason))
          g_error ("mcat_join: catalog `%s' is not compatible with `%s': %s.", cat_filename[i], cat_filename[0], reason);

        if (!ncm_mset_catalog_desc_has_same_order (desc, desc_i))
          g_error ("mcat_join: catalog `%s' does not have the parameters of `%s' in the same order.", cat_filename[i], cat_filename[0]);

        if (burnin % nchains != 0)
          g_error ("Burnin[%u] not multiple of nchains %u.", i, nchains);

        mcat       = ncm_mset_catalog_new_from_storage_full (storage, burnin, 1);
        max_time_i = ncm_mset_catalog_max_time (mcat);
        mmax_time  = GSL_MAX (mmax_time, max_time_i);

        g_ptr_array_add (mcats_array, mcat);
        g_array_append_val (max_time_array, max_time_i);
        ncm_mset_catalog_storage_fits_free (storage_fits);
        g_free (reason);
      }

      if (out == NULL)
      {
        guint k = 0;

        out_name = g_strdup_printf ("joined_mcat.fits");

        while (g_file_test (out_name, G_FILE_TEST_EXISTS))
        {
          g_free (out_name);
          out_name = g_strdup_printf ("joined_mcat_%d.fits", k++);
        }
      }
      else if (g_file_test (out, G_FILE_TEST_EXISTS))
      {
        gchar *out_base = ncm_util_basename_fits (out);
        guint k         = 0;

        g_warning ("out file already exists, renaming.");
        out_name = g_strdup_printf ("%s_ren.fits", out_base);

        while (g_file_test (out_name, G_FILE_TEST_EXISTS))
        {
          g_free (out_name);
          out_name = g_strdup_printf ("%s_ren%d.fits", out_base, k++);
        }

        g_free (out_base);
      }
      else
      {
        out_name = g_strdup (out);
      }

      /* The joined catalog is written at the end, in one recording. */
      out_fits = ncm_mset_catalog_storage_fits_new (out_name, FALSE, NULL);
      g_object_set (out_fits, "smode", NCM_MSET_CATALOG_SYNC_DISABLE, NULL);
      ncm_mset_catalog_storage_set_desc (NCM_MSET_CATALOG_STORAGE (out_fits), desc, NULL);
      mcat_out = ncm_mset_catalog_new_from_storage (NCM_MSET_CATALOG_STORAGE (out_fits));
      ncm_mset_catalog_set_sampler (mcat_out, "mcat_join");

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

      ncm_mset_catalog_clear (&mcat_out);
      ncm_mset_catalog_storage_fits_free (out_fits);
      ncm_mset_catalog_storage_fits_free (first_fits);
      g_ptr_array_unref (mcats_array);
      g_array_unref (max_time_array);
      g_free (out_name);
    }
  }

  return 0;
}

