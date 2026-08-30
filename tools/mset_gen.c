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

/*
 * Extracts the leading type-name token of a --model string (e.g.
 * "NcHIReionCamb" out of "NcHIReionCamb{'z-re':<13.0>}") and checks
 * whether that type is a submodel, without constructing anything.
 */
static gboolean
_mset_gen_str_is_submodel (const gchar *model_str)
{
  gchar *brace         = strchr (model_str, '{');
  gchar *type_name     = brace != NULL ? g_strndup (model_str, brace - model_str) : g_strdup (model_str);
  GType gtype          = g_type_from_name (g_strstrip (type_name));
  gboolean is_submodel = FALSE;

  if ((gtype != 0) && g_type_is_a (gtype, NCM_TYPE_MODEL))
  {
    NcmModelClass *model_class = NCM_MODEL_CLASS (g_type_class_ref (gtype));

    is_submodel = model_class->is_submodel;
    g_type_class_unref (model_class);
  }

  g_free (type_name);

  return is_submodel;
}

/*
 * For every construction-only, #NcmModel-typed property declared on
 * @host_gtype's class (a submodel slot installed via
 * ncm_model_class_set_submodel()), appends a reference to any matching
 * submodel already registered in @ser (by its own namespace, see the
 * pass-1 loop in main()) to @obj_ser, letting it be attached at
 * construction time instead of via the no-longer-legal
 * ncm_model_add_submodel() post-hoc path.
 */
static void
_mset_gen_inject_submodel_props (NcmSerialize *ser, GType host_gtype, GString *obj_ser, gboolean *needs_comma)
{
  NcmModelClass *host_class = NCM_MODEL_CLASS (g_type_class_ref (host_gtype));
  guint n_props             = 0;
  GParamSpec **props        = g_object_class_list_properties (G_OBJECT_CLASS (host_class), &n_props);
  guint i;

  for (i = 0; i < n_props; i++)
  {
    GParamSpec *pspec = props[i];

    if (G_IS_PARAM_SPEC_OBJECT (pspec) &&
        (pspec->flags & G_PARAM_CONSTRUCT_ONLY) &&
        g_type_is_a (G_PARAM_SPEC_VALUE_TYPE (pspec), NCM_TYPE_MODEL))
    {
      NcmModelID slot_mid = ncm_mset_get_id_by_type (G_PARAM_SPEC_VALUE_TYPE (pspec));

      if (slot_mid >= 0)
      {
        const gchar *slot_ns = ncm_mset_get_ns_by_id (slot_mid);
        gpointer submodel    = ncm_serialize_peek_by_name (ser, slot_ns);

        if (submodel != NULL)
        {
          if (*needs_comma)
            g_string_append (obj_ser, ", ");

          g_string_append_printf (obj_ser, "\'%s\':<(\'%s[%s]\', @a{sv} {})>",
                                  pspec->name, G_OBJECT_TYPE_NAME (submodel), slot_ns);
          *needs_comma = TRUE;
        }
      }
    }
  }

  g_free (props);
  g_type_class_unref (host_class);
}

/*
 * Splices any injected submodel-slot properties into @model_str's own
 * property dict (creating one if it has none), returning a new string
 * ready for ncm_serialize_from_string().
 */
static gchar *
_mset_gen_str_with_submodels (NcmSerialize *ser, const gchar *model_str)
{
  gchar *brace     = strchr (model_str, '{');
  gchar *type_name = brace != NULL ? g_strndup (model_str, brace - model_str) : g_strdup (model_str);
  GType gtype      = g_type_from_name (g_strstrip (type_name));
  GString *full    = g_string_new (model_str);
  gboolean needs_comma;

  if (brace != NULL)
  {
    /* Already has an explicit (possibly empty) property dict -- drop the
     * trailing '}' and only prefix a comma if it wasn't empty. */
    g_string_truncate (full, full->len - 1);
    needs_comma = full->str[full->len - 1] != '{';

    if (needs_comma)
      g_string_append (full, ", ");
  }
  else
  {
    g_string_append (full, "{");
    needs_comma = FALSE;
  }

  if (gtype != 0)
    _mset_gen_inject_submodel_props (ser, gtype, full, &needs_comma);

  g_string_append (full, "}");
  g_free (type_name);

  return g_string_free (full, FALSE);
}

gint
main (gint argc, gchar *argv[])
{
  gchar **models         = NULL;
  gchar *outfile         = NULL;
  gboolean save_comments = FALSE;
  gboolean overwrite     = FALSE;

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

  ncm_cfg_init_full_ptr (&argc, &argv);

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
    NcmMSet *mset     = ncm_mset_empty_new ();
    NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
    guint nmodels     = g_strv_length (models);
    guint i;

    /* Pass 1: build every submodel arg first and register it in @ser
     * under its own namespace, so pass 2's main-model args can reference
     * it as a construction-time property -- ncm_model_add_submodel() can
     * no longer attach a submodel to an already-constructed host
     * declaring a typed slot for it. */
    for (i = 0; i < nmodels; i++)
    {
      if (_mset_gen_str_is_submodel (models[i]))
      {
        NcmModel *model = NCM_MODEL (ncm_serialize_from_string (ser, models[i]));
        const gchar *ns;

        g_assert (NCM_IS_MODEL (model));
        ns = ncm_mset_get_ns_by_id (ncm_model_id (model));
        ncm_serialize_set (ser, model, ns, FALSE);
        ncm_model_free (model);
      }
    }

    /* Pass 2: build every main-model arg, with any pass-1-registered
     * submodel matching one of its declared typed slots spliced into its
     * own construction string. */
    for (i = 0; i < nmodels; i++)
    {
      NcmModel *model;
      gchar *full_str;

      if (_mset_gen_str_is_submodel (models[i]))
        continue;

      full_str = _mset_gen_str_with_submodels (ser, models[i]);
      model    = NCM_MODEL (ncm_serialize_from_string (ser, full_str));
      g_free (full_str);

      g_assert (NCM_IS_MODEL (model));
      g_assert (!ncm_model_is_submodel (model));

      ncm_mset_set (mset, model, NULL);
      ncm_model_free (model);
    }

    outfile = outfile != NULL ? outfile : "models.mset";

    if (g_file_test (outfile, G_FILE_TEST_EXISTS) && !overwrite)
    {
      g_print ("File `%s' already exists, choose --overwrite/-w to overwrite.\n", outfile);
      exit (1);
    }

    ncm_mset_save (mset, ser, outfile, save_comments, NULL);

    ncm_mset_clear (&mset);
    ncm_serialize_clear (&ser);
  }

  return 0;
}

