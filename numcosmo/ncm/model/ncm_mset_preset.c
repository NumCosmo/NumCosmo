/***************************************************************************
 *            ncm_mset_preset.c
 *
 *  Thu Nov 19 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_preset.c
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

/**
 * NcmMSetPreset:
 *
 * Parameter values, scales, bounds, absolute tolerances and fit flags taken from a
 * publication, applied to the parameters of an #NcmMSet.
 *
 * A preset describes one publication and one data combination: #NcmMSetPreset:desc
 * gives the reference and the table the numbers come from,
 * #NcmMSetPreset:data-combination the data used and #NcmMSetPreset:statistic the
 * statistic the values are (for example "best fit" or "mean").
 *
 * The entries are six dictionaries, one per parameter attribute, each keyed on the
 * full parameter names of ncm_mset_param_get_by_full_name() ("NcHICosmo:H0"):
 * #NcmMSetPreset:value, #NcmMSetPreset:scale, #NcmMSetPreset:lower-bound,
 * #NcmMSetPreset:upper-bound, #NcmMSetPreset:abstol (doubles) and #NcmMSetPreset:fit
 * (booleans). A parameter may appear in any subset of them; the attributes it does not
 * appear in are left unchanged.
 *
 * Parameter names depend on the reparametrization in force. #NcmMSetPreset:reparam maps
 * a model name (as in ncm_mset_peek_by_name()) to the type name of the
 * #NcmReparam that model must have. A model that has parameters in the preset and no
 * entry in #NcmMSetPreset:reparam must have no reparametrization.
 *
 * The data directory holds the distributed presets in its #NCM_MSET_PRESET_DATA_DIR
 * subdirectory, one YAML file per preset, written in the form of
 * ncm_serialize_to_yaml(). ncm_mset_preset_new_from_name() loads one by its file name
 * without the ".yaml" suffix and ncm_mset_preset_list_names() lists them.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/model/ncm_mset_preset.h"
#include "ncm/model/ncm_reparam.h"
#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <string.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* *INDENT-OFF* */
G_DEFINE_QUARK (ncm-mset-preset-error, ncm_mset_preset_error)
/* *INDENT-ON* */

typedef enum _NcmMSetPresetAttr
{
  NCM_MSET_PRESET_ATTR_VALUE = 0,
  NCM_MSET_PRESET_ATTR_SCALE,
  NCM_MSET_PRESET_ATTR_LOWER_BOUND,
  NCM_MSET_PRESET_ATTR_UPPER_BOUND,
  NCM_MSET_PRESET_ATTR_ABSTOL,
  NCM_MSET_PRESET_ATTR_FIT,
  NCM_MSET_PRESET_ATTR_LEN,
} NcmMSetPresetAttr;

struct _NcmMSetPreset
{
  GObject parent_instance;
  gchar *desc;
  gchar *data_combination;
  gchar *statistic;
  NcmVarDict *reparam;
  NcmVarDict *attr[NCM_MSET_PRESET_ATTR_LEN];
};

/*
 * One parameter touched by a preset: where it lives, its attributes after the preset
 * is applied and which of them the preset sets.
 */
typedef struct _NcmMSetPresetEntry
{
  NcmModel *model;
  guint pid;
  gdouble dval[NCM_MSET_PRESET_ATTR_FIT];
  gboolean fit;
  gboolean set[NCM_MSET_PRESET_ATTR_LEN];
} NcmMSetPresetEntry;

enum
{
  PROP_0,
  PROP_DESC,
  PROP_DATA_COMBINATION,
  PROP_STATISTIC,
  PROP_REPARAM,
  PROP_VALUE,
  PROP_SCALE,
  PROP_LOWER_BOUND,
  PROP_UPPER_BOUND,
  PROP_ABSTOL,
  PROP_FIT,
};

static const gchar *_ncm_mset_preset_attr_names[NCM_MSET_PRESET_ATTR_LEN] = {
  "value", "scale", "lower-bound", "upper-bound", "abstol", "fit",
};

G_DEFINE_TYPE (NcmMSetPreset, ncm_mset_preset, G_TYPE_OBJECT)

static void
ncm_mset_preset_init (NcmMSetPreset *preset)
{
  guint a;

  preset->desc             = NULL;
  preset->data_combination = NULL;
  preset->statistic        = NULL;
  preset->reparam          = NULL;

  for (a = 0; a < NCM_MSET_PRESET_ATTR_LEN; a++)
    preset->attr[a] = NULL;
}

/* Replaces *@dest by a new dictionary holding the entries of @src (none when NULL). */
static void
_ncm_mset_preset_copy_var_dict (NcmVarDict **dest, NcmVarDict *src)
{
  NcmVarDict *copy = ncm_var_dict_new ();

  if (src != NULL)
  {
    GHashTableIter iter;
    gpointer key, value;

    g_hash_table_iter_init (&iter, (GHashTable *) src);

    while (g_hash_table_iter_next (&iter, &key, &value))
      g_hash_table_insert ((GHashTable *) copy, g_strdup (key), g_variant_ref (value));
  }

  ncm_var_dict_clear (dest);
  *dest = copy;
}

static void
_ncm_mset_preset_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetPreset *preset = NCM_MSET_PRESET (object);

  g_return_if_fail (NCM_IS_MSET_PRESET (object));

  switch (prop_id)
  {
    case PROP_DESC:
      g_free (preset->desc);
      preset->desc = g_value_dup_string (value);
      break;
    case PROP_DATA_COMBINATION:
      g_free (preset->data_combination);
      preset->data_combination = g_value_dup_string (value);
      break;
    case PROP_STATISTIC:
      g_free (preset->statistic);
      preset->statistic = g_value_dup_string (value);
      break;
    case PROP_REPARAM:
      _ncm_mset_preset_copy_var_dict (&preset->reparam, g_value_get_boxed (value));
      break;
    case PROP_VALUE:
    case PROP_SCALE:
    case PROP_LOWER_BOUND:
    case PROP_UPPER_BOUND:
    case PROP_ABSTOL:
    case PROP_FIT:
      _ncm_mset_preset_copy_var_dict (&preset->attr[prop_id - PROP_VALUE], g_value_get_boxed (value));
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_preset_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetPreset *preset = NCM_MSET_PRESET (object);

  g_return_if_fail (NCM_IS_MSET_PRESET (object));

  switch (prop_id)
  {
    case PROP_DESC:
      g_value_set_string (value, preset->desc);
      break;
    case PROP_DATA_COMBINATION:
      g_value_set_string (value, preset->data_combination);
      break;
    case PROP_STATISTIC:
      g_value_set_string (value, preset->statistic);
      break;
    case PROP_REPARAM:
      g_value_set_boxed (value, preset->reparam);
      break;
    case PROP_VALUE:
    case PROP_SCALE:
    case PROP_LOWER_BOUND:
    case PROP_UPPER_BOUND:
    case PROP_ABSTOL:
    case PROP_FIT:
      g_value_set_boxed (value, preset->attr[prop_id - PROP_VALUE]);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_preset_finalize (GObject *object)
{
  NcmMSetPreset *preset = NCM_MSET_PRESET (object);
  guint a;

  g_clear_pointer (&preset->desc, g_free);
  g_clear_pointer (&preset->data_combination, g_free);
  g_clear_pointer (&preset->statistic, g_free);
  ncm_var_dict_clear (&preset->reparam);

  for (a = 0; a < NCM_MSET_PRESET_ATTR_LEN; a++)
    ncm_var_dict_clear (&preset->attr[a]);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_preset_parent_class)->finalize (object);
}

static void
ncm_mset_preset_class_init (NcmMSetPresetClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_mset_preset_set_property;
  object_class->get_property = &_ncm_mset_preset_get_property;
  object_class->finalize     = &_ncm_mset_preset_finalize;

  /**
   * NcmMSetPreset:desc:
   *
   * The publication and the table the entries come from.
   */
  g_object_class_install_property (object_class,
                                   PROP_DESC,
                                   g_param_spec_string ("desc",
                                                        NULL,
                                                        "Reference of the preset",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetPreset:data-combination:
   *
   * The data combination the entries were obtained from.
   */
  g_object_class_install_property (object_class,
                                   PROP_DATA_COMBINATION,
                                   g_param_spec_string ("data-combination",
                                                        NULL,
                                                        "Data combination of the preset",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetPreset:statistic:
   *
   * The statistic of the posterior given in #NcmMSetPreset:value, for example "best fit"
   * or "mean".
   */
  g_object_class_install_property (object_class,
                                   PROP_STATISTIC,
                                   g_param_spec_string ("statistic",
                                                        NULL,
                                                        "Statistic of the values",
                                                        NULL,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetPreset:reparam:
   *
   * Model name to the type name of the #NcmReparam the model must have.
   */
  g_object_class_install_property (object_class,
                                   PROP_REPARAM,
                                   g_param_spec_boxed ("reparam",
                                                       NULL,
                                                       "Required reparametrization of each model",
                                                       NCM_TYPE_VAR_DICT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetPreset:value:
   *
   * Full parameter name to the parameter value.
   */
  g_object_class_install_property (object_class,
                                   PROP_VALUE,
                                   g_param_spec_boxed ("value",
                                                       NULL,
                                                       "Parameter values",
                                                       NCM_TYPE_VAR_DICT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetPreset:scale:
   *
   * Full parameter name to the parameter scale, #NcmSParam:scale.
   */
  g_object_class_install_property (object_class,
                                   PROP_SCALE,
                                   g_param_spec_boxed ("scale",
                                                       NULL,
                                                       "Parameter scales",
                                                       NCM_TYPE_VAR_DICT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetPreset:lower-bound:
   *
   * Full parameter name to the parameter lower bound, #NcmSParam:lower-bound.
   */
  g_object_class_install_property (object_class,
                                   PROP_LOWER_BOUND,
                                   g_param_spec_boxed ("lower-bound",
                                                       NULL,
                                                       "Parameter lower bounds",
                                                       NCM_TYPE_VAR_DICT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetPreset:upper-bound:
   *
   * Full parameter name to the parameter upper bound, #NcmSParam:upper-bound.
   */
  g_object_class_install_property (object_class,
                                   PROP_UPPER_BOUND,
                                   g_param_spec_boxed ("upper-bound",
                                                       NULL,
                                                       "Parameter upper bounds",
                                                       NCM_TYPE_VAR_DICT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetPreset:abstol:
   *
   * Full parameter name to the parameter absolute tolerance,
   * #NcmSParam:absolute-tolerance.
   */
  g_object_class_install_property (object_class,
                                   PROP_ABSTOL,
                                   g_param_spec_boxed ("abstol",
                                                       NULL,
                                                       "Parameter absolute tolerances",
                                                       NCM_TYPE_VAR_DICT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetPreset:fit:
   *
   * Full parameter name to whether the parameter is free (%TRUE) or fixed (%FALSE).
   */
  g_object_class_install_property (object_class,
                                   PROP_FIT,
                                   g_param_spec_boxed ("fit",
                                                       NULL,
                                                       "Parameter fit flags",
                                                       NCM_TYPE_VAR_DICT,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/**
 * ncm_mset_preset_new:
 * @desc: (nullable): the publication and table the entries come from
 *
 * Creates a preset with no entries.
 *
 * Returns: (transfer full): a new #NcmMSetPreset.
 */
NcmMSetPreset *
ncm_mset_preset_new (const gchar *desc)
{
  NcmMSetPreset *preset = g_object_new (NCM_TYPE_MSET_PRESET,
                                        "desc", desc,
                                        NULL);

  return preset;
}

/**
 * ncm_mset_preset_new_from_file:
 * @filename: path of a YAML file
 * @error: a #GError
 *
 * Loads the preset serialized in @filename by ncm_serialize_to_yaml_file(). A missing
 * file sets %NCM_MSET_PRESET_ERROR_NOT_FOUND and a file holding another object
 * %NCM_MSET_PRESET_ERROR_INVALID_FILE; malformed YAML aborts in the deserializer.
 *
 * Returns: (transfer full) (nullable): the preset, or %NULL on error.
 */
NcmMSetPreset *
ncm_mset_preset_new_from_file (const gchar *filename, GError **error)
{
  NcmSerialize *ser;
  GObject *obj;

  if (!g_file_test (filename, G_FILE_TEST_IS_REGULAR))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_NOT_FOUND,
                                "ncm_mset_preset_new_from_file: file `%s' not found.", filename);

    return NULL;
  }

  ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
  obj = ncm_serialize_from_yaml_file (ser, filename);
  ncm_serialize_free (ser);

  if (!NCM_IS_MSET_PRESET (obj))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_INVALID_FILE,
                                "ncm_mset_preset_new_from_file: `%s' holds a `%s', not a NcmMSetPreset.",
                                filename, G_OBJECT_TYPE_NAME (obj));
    g_object_unref (obj);

    return NULL;
  }

  return NCM_MSET_PRESET (obj);
}

/**
 * ncm_mset_preset_new_from_name:
 * @name: a preset name, as listed by ncm_mset_preset_list_names()
 * @error: a #GError
 *
 * Loads the preset @name from the #NCM_MSET_PRESET_DATA_DIR subdirectory of the data
 * directory, searched as in ncm_cfg_get_data_filename(). A name with a directory
 * separator, or with no preset file, sets %NCM_MSET_PRESET_ERROR_NOT_FOUND.
 *
 * Returns: (transfer full) (nullable): the preset, or %NULL on error.
 */
NcmMSetPreset *
ncm_mset_preset_new_from_name (const gchar *name, GError **error)
{
  NcmMSetPreset *preset;
  gchar *basename;
  gchar *relname;
  gchar *filename;

  if ((name[0] == '\0') || (strchr (name, '/') != NULL) || (strchr (name, G_DIR_SEPARATOR) != NULL))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_NOT_FOUND,
                                "ncm_mset_preset_new_from_name: invalid preset name `%s'.", name);

    return NULL;
  }

  basename = g_strconcat (name, ".yaml", NULL);
  relname  = g_build_filename (NCM_MSET_PRESET_DATA_DIR, basename, NULL);
  filename = ncm_cfg_get_data_filename (relname, FALSE);

  g_free (basename);
  g_free (relname);

  if (filename == NULL)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_NOT_FOUND,
                                "ncm_mset_preset_new_from_name: preset `%s' not found.", name);

    return NULL;
  }

  preset = ncm_mset_preset_new_from_file (filename, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, g_free (filename), NULL, "ncm_mset_preset_new_from_name: ");

  g_free (filename);

  return preset;
}

/**
 * ncm_mset_preset_ref:
 * @preset: a #NcmMSetPreset
 *
 * Increases the reference count of @preset by one.
 *
 * Returns: (transfer full): @preset.
 */
NcmMSetPreset *
ncm_mset_preset_ref (NcmMSetPreset *preset)
{
  return g_object_ref (preset);
}

/**
 * ncm_mset_preset_free:
 * @preset: a #NcmMSetPreset
 *
 * Decreases the reference count of @preset by one.
 *
 */
void
ncm_mset_preset_free (NcmMSetPreset *preset)
{
  g_object_unref (preset);
}

/**
 * ncm_mset_preset_clear:
 * @preset: a #NcmMSetPreset
 *
 * If *@preset is not %NULL, decreases its reference count by one and sets *@preset to
 * %NULL.
 *
 */
void
ncm_mset_preset_clear (NcmMSetPreset **preset)
{
  g_clear_object (preset);
}

static gint
_ncm_mset_preset_cmp_str (gconstpointer a, gconstpointer b)
{
  return g_strcmp0 (*(const gchar * const *) a, *(const gchar * const *) b);
}

/**
 * ncm_mset_preset_list_names:
 *
 * Lists the presets in the #NCM_MSET_PRESET_DATA_DIR subdirectory of the first data
 * directory, searched as in ncm_cfg_get_data_filename(), that has one: the names of
 * its ".yaml" files without the suffix, sorted. The list is empty when no data
 * directory has the subdirectory.
 *
 * Returns: (transfer full) (array zero-terminated=1): the preset names.
 */
GStrv
ncm_mset_preset_list_names (void)
{
  GPtrArray *names = g_ptr_array_new ();
  gchar *dirname   = ncm_cfg_get_data_filename (NCM_MSET_PRESET_DATA_DIR, FALSE);

  if (dirname != NULL)
  {
    GDir *dir = g_dir_open (dirname, 0, NULL);

    if (dir != NULL)
    {
      const gchar *entry;

      while ((entry = g_dir_read_name (dir)) != NULL)
      {
        if (g_str_has_suffix (entry, ".yaml"))
          g_ptr_array_add (names, g_strndup (entry, strlen (entry) - strlen (".yaml")));
      }

      g_dir_close (dir);
    }

    g_free (dirname);
  }

  g_ptr_array_sort (names, &_ncm_mset_preset_cmp_str);
  g_ptr_array_add (names, NULL);

  return (GStrv) g_ptr_array_free (names, FALSE);
}

/**
 * ncm_mset_preset_peek_desc:
 * @preset: a #NcmMSetPreset
 *
 * Returns: (transfer none) (nullable): #NcmMSetPreset:desc.
 */
const gchar *
ncm_mset_preset_peek_desc (NcmMSetPreset *preset)
{
  return preset->desc;
}

/**
 * ncm_mset_preset_peek_data_combination:
 * @preset: a #NcmMSetPreset
 *
 * Returns: (transfer none) (nullable): #NcmMSetPreset:data-combination.
 */
const gchar *
ncm_mset_preset_peek_data_combination (NcmMSetPreset *preset)
{
  return preset->data_combination;
}

/**
 * ncm_mset_preset_peek_statistic:
 * @preset: a #NcmMSetPreset
 *
 * Returns: (transfer none) (nullable): #NcmMSetPreset:statistic.
 */
const gchar *
ncm_mset_preset_peek_statistic (NcmMSetPreset *preset)
{
  return preset->statistic;
}

/**
 * ncm_mset_preset_set_reparam:
 * @preset: a #NcmMSetPreset
 * @model_name: a model name, as in ncm_mset_peek_by_name()
 * @reparam_type: the type name of the #NcmReparam the model must have
 *
 * Sets the reparametrization required for @model_name.
 *
 */
void
ncm_mset_preset_set_reparam (NcmMSetPreset *preset, const gchar *model_name, const gchar *reparam_type)
{
  ncm_var_dict_set_string (preset->reparam, model_name, reparam_type);
}

/**
 * ncm_mset_preset_set_value:
 * @preset: a #NcmMSetPreset
 * @fullname: a full parameter name
 * @value: the parameter value
 *
 * Sets the value of @fullname in #NcmMSetPreset:value.
 *
 */
void
ncm_mset_preset_set_value (NcmMSetPreset *preset, const gchar *fullname, gdouble value)
{
  ncm_var_dict_set_double (preset->attr[NCM_MSET_PRESET_ATTR_VALUE], fullname, value);
}

/**
 * ncm_mset_preset_set_scale:
 * @preset: a #NcmMSetPreset
 * @fullname: a full parameter name
 * @scale: the parameter scale
 *
 * Sets the scale of @fullname in #NcmMSetPreset:scale.
 *
 */
void
ncm_mset_preset_set_scale (NcmMSetPreset *preset, const gchar *fullname, gdouble scale)
{
  ncm_var_dict_set_double (preset->attr[NCM_MSET_PRESET_ATTR_SCALE], fullname, scale);
}

/**
 * ncm_mset_preset_set_lower_bound:
 * @preset: a #NcmMSetPreset
 * @fullname: a full parameter name
 * @lower_bound: the parameter lower bound
 *
 * Sets the lower bound of @fullname in #NcmMSetPreset:lower-bound.
 *
 */
void
ncm_mset_preset_set_lower_bound (NcmMSetPreset *preset, const gchar *fullname, gdouble lower_bound)
{
  ncm_var_dict_set_double (preset->attr[NCM_MSET_PRESET_ATTR_LOWER_BOUND], fullname, lower_bound);
}

/**
 * ncm_mset_preset_set_upper_bound:
 * @preset: a #NcmMSetPreset
 * @fullname: a full parameter name
 * @upper_bound: the parameter upper bound
 *
 * Sets the upper bound of @fullname in #NcmMSetPreset:upper-bound.
 *
 */
void
ncm_mset_preset_set_upper_bound (NcmMSetPreset *preset, const gchar *fullname, gdouble upper_bound)
{
  ncm_var_dict_set_double (preset->attr[NCM_MSET_PRESET_ATTR_UPPER_BOUND], fullname, upper_bound);
}

/**
 * ncm_mset_preset_set_abstol:
 * @preset: a #NcmMSetPreset
 * @fullname: a full parameter name
 * @abstol: the parameter absolute tolerance
 *
 * Sets the absolute tolerance of @fullname in #NcmMSetPreset:abstol.
 *
 */
void
ncm_mset_preset_set_abstol (NcmMSetPreset *preset, const gchar *fullname, gdouble abstol)
{
  ncm_var_dict_set_double (preset->attr[NCM_MSET_PRESET_ATTR_ABSTOL], fullname, abstol);
}

/**
 * ncm_mset_preset_set_fit:
 * @preset: a #NcmMSetPreset
 * @fullname: a full parameter name
 * @fit: whether the parameter is free
 *
 * Sets the fit flag of @fullname in #NcmMSetPreset:fit.
 *
 */
void
ncm_mset_preset_set_fit (NcmMSetPreset *preset, const gchar *fullname, gboolean fit)
{
  ncm_var_dict_set_boolean (preset->attr[NCM_MSET_PRESET_ATTR_FIT], fullname, fit);
}

/**
 * ncm_mset_preset_remove:
 * @preset: a #NcmMSetPreset
 * @fullname: a full parameter name
 *
 * Removes @fullname from every attribute dictionary of @preset, so that
 * ncm_mset_preset_apply() leaves that parameter unchanged.
 *
 * Returns: whether @preset had an entry for @fullname.
 */
gboolean
ncm_mset_preset_remove (NcmMSetPreset *preset, const gchar *fullname)
{
  gboolean found = FALSE;
  guint a;

  for (a = 0; a < NCM_MSET_PRESET_ATTR_LEN; a++)
    found = g_hash_table_remove ((GHashTable *) preset->attr[a], fullname) || found;

  return found;
}

static void _ncm_mset_preset_check_named_reparam (NcmMSetPreset *preset, NcmMSet *mset, GHashTable *checked, GError **error);
static void _ncm_mset_preset_collect (NcmMSetPreset *preset, NcmMSet *mset, GArray *entries, GError **error);
static void _ncm_mset_preset_check_unnamed_reparam (GArray *entries, GHashTable *checked, GError **error);
static void _ncm_mset_preset_check_entries (GArray *entries, GError **error);

/**
 * ncm_mset_preset_apply:
 * @preset: a #NcmMSetPreset
 * @mset: a #NcmMSet
 * @error: a #GError
 *
 * Sets, for every parameter in the preset, the attributes the preset gives. The
 * whole preset is checked before any parameter changes, so on error @mset is
 * unchanged. The checks are, in order: every model in #NcmMSetPreset:reparam is in
 * @mset (%NCM_MSET_PRESET_ERROR_PARAM_NOT_FOUND) and has the required
 * reparametrization (%NCM_MSET_PRESET_ERROR_REPARAM_MISMATCH); every full name is a
 * parameter of @mset (%NCM_MSET_PRESET_ERROR_PARAM_NOT_FOUND) and every entry has its
 * type (%NCM_MSET_PRESET_ERROR_INVALID_ENTRY); every other model with parameters in the
 * preset has no reparametrization (%NCM_MSET_PRESET_ERROR_REPARAM_MISMATCH); and, with
 * the preset applied, each
 * parameter has lower bound < upper bound, its value within the bounds, a positive
 * scale and a non-negative absolute tolerance (%NCM_MSET_PRESET_ERROR_INVALID_ENTRY).
 * As with ncm_model_param_set_ftype(), a changed fit flag takes effect in @mset at the
 * next ncm_mset_prepare_fparam_map().
 *
 */
void
ncm_mset_preset_apply (NcmMSetPreset *preset, NcmMSet *mset, GError **error)
{
  GArray *entries     = g_array_new (FALSE, FALSE, sizeof (NcmMSetPresetEntry));
  GHashTable *checked = g_hash_table_new (g_direct_hash, g_direct_equal);
  guint i;

  _ncm_mset_preset_check_named_reparam (preset, mset, checked, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, { g_array_unref (entries);
                                      g_hash_table_unref (checked);
                             }, , "ncm_mset_preset_apply: ");

  _ncm_mset_preset_collect (preset, mset, entries, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, { g_array_unref (entries);
                                      g_hash_table_unref (checked);
                             }, , "ncm_mset_preset_apply: ");

  _ncm_mset_preset_check_unnamed_reparam (entries, checked, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, { g_array_unref (entries);
                                      g_hash_table_unref (checked);
                             }, , "ncm_mset_preset_apply: ");

  _ncm_mset_preset_check_entries (entries, error);
  NCM_UTIL_ON_ERROR_FORWARD (error, { g_array_unref (entries);
                                      g_hash_table_unref (checked);
                             }, , "ncm_mset_preset_apply: ");

  for (i = 0; i < entries->len; i++)
  {
    const NcmMSetPresetEntry *e = &g_array_index (entries, NcmMSetPresetEntry, i);

    if (e->set[NCM_MSET_PRESET_ATTR_LOWER_BOUND])
      ncm_model_param_set_lower_bound (e->model, e->pid, e->dval[NCM_MSET_PRESET_ATTR_LOWER_BOUND]);

    if (e->set[NCM_MSET_PRESET_ATTR_UPPER_BOUND])
      ncm_model_param_set_upper_bound (e->model, e->pid, e->dval[NCM_MSET_PRESET_ATTR_UPPER_BOUND]);

    if (e->set[NCM_MSET_PRESET_ATTR_SCALE])
      ncm_model_param_set_scale (e->model, e->pid, e->dval[NCM_MSET_PRESET_ATTR_SCALE]);

    if (e->set[NCM_MSET_PRESET_ATTR_ABSTOL])
      ncm_model_param_set_abstol (e->model, e->pid, e->dval[NCM_MSET_PRESET_ATTR_ABSTOL]);

    if (e->set[NCM_MSET_PRESET_ATTR_FIT])
      ncm_model_param_set_ftype (e->model, e->pid, e->fit ? NCM_PARAM_TYPE_FREE : NCM_PARAM_TYPE_FIXED);

    if (e->set[NCM_MSET_PRESET_ATTR_VALUE])
      ncm_model_param_set (e->model, e->pid, e->dval[NCM_MSET_PRESET_ATTR_VALUE]);
  }

  g_array_unref (entries);
  g_hash_table_unref (checked);
}

/*
 * Checks that every model named in #NcmMSetPreset:reparam is in @mset and has the
 * required reparametrization, and adds it to @checked.
 */
static void
_ncm_mset_preset_check_named_reparam (NcmMSetPreset *preset, NcmMSet *mset, GHashTable *checked, GError **error)
{
  GHashTableIter iter;
  gpointer key, value;

  g_hash_table_iter_init (&iter, (GHashTable *) preset->reparam);

  while (g_hash_table_iter_next (&iter, &key, &value))
  {
    const gchar *model_name = key;
    GVariant *var           = value;
    NcmModel *model         = ncm_mset_peek_by_name (mset, model_name, error);
    NcmReparam *reparam;
    const gchar *required;
    const gchar *actual;

    NCM_UTIL_ON_ERROR_RETURN (error, , );

    if (model == NULL)
    {
      ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_PARAM_NOT_FOUND,
                                  "model `%s' not found in the NcmMSet.", model_name);

      return;
    }

    if (!g_variant_is_of_type (var, G_VARIANT_TYPE_STRING))
    {
      ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_INVALID_ENTRY,
                                  "reparam of `%s' must be a string.", model_name);

      return;
    }

    reparam  = ncm_model_peek_reparam (model);
    required = g_variant_get_string (var, NULL);
    actual   = (reparam != NULL) ? G_OBJECT_TYPE_NAME (reparam) : "none";

    if ((reparam == NULL) || (g_strcmp0 (required, actual) != 0))
    {
      ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_REPARAM_MISMATCH,
                                  "model `%s' has reparametrization `%s', the preset requires `%s'.",
                                  model_name, actual, required);

      return;
    }

    g_hash_table_add (checked, model);
  }
}

/*
 * Fills @entries with one entry per full name in the attribute dictionaries: the
 * model and index of the parameter, its current attributes overwritten by those the
 * preset gives.
 */
static void
_ncm_mset_preset_collect (NcmMSetPreset *preset, NcmMSet *mset, GArray *entries, GError **error)
{
  GHashTable *index = g_hash_table_new (g_str_hash, g_str_equal);
  guint a;

  for (a = 0; a < NCM_MSET_PRESET_ATTR_LEN; a++)
  {
    GHashTableIter iter;
    gpointer key, value;

    g_hash_table_iter_init (&iter, (GHashTable *) preset->attr[a]);

    while (g_hash_table_iter_next (&iter, &key, &value))
    {
      const gchar *fullname = key;
      GVariant *var         = value;
      gpointer pos;
      NcmMSetPresetEntry *e;

      if (!g_hash_table_lookup_extended (index, fullname, NULL, &pos))
      {
        NcmMSetPIndex *pi = ncm_mset_param_get_by_full_name (mset, fullname, error);
        NcmMSetPresetEntry new_entry;
        guint b;

        NCM_UTIL_ON_ERROR_RETURN (error, g_hash_table_unref (index), );

        if (pi == NULL)
        {
          ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_PARAM_NOT_FOUND,
                                      "parameter `%s' not found in the NcmMSet.", fullname);
          g_hash_table_unref (index);

          return;
        }

        new_entry.model                                  = ncm_mset_peek (mset, pi->mid);
        new_entry.pid                                    = pi->pid;
        new_entry.dval[NCM_MSET_PRESET_ATTR_VALUE]       = ncm_model_param_get (new_entry.model, pi->pid);
        new_entry.dval[NCM_MSET_PRESET_ATTR_SCALE]       = ncm_model_param_get_scale (new_entry.model, pi->pid);
        new_entry.dval[NCM_MSET_PRESET_ATTR_LOWER_BOUND] = ncm_model_param_get_lower_bound (new_entry.model, pi->pid);
        new_entry.dval[NCM_MSET_PRESET_ATTR_UPPER_BOUND] = ncm_model_param_get_upper_bound (new_entry.model, pi->pid);
        new_entry.dval[NCM_MSET_PRESET_ATTR_ABSTOL]      = ncm_model_param_get_abstol (new_entry.model, pi->pid);
        new_entry.fit                                    = (ncm_model_param_get_ftype (new_entry.model, pi->pid) == NCM_PARAM_TYPE_FREE);

        for (b = 0; b < NCM_MSET_PRESET_ATTR_LEN; b++)
          new_entry.set[b] = FALSE;

        ncm_mset_pindex_free (pi);

        pos = GUINT_TO_POINTER (entries->len);
        g_array_append_val (entries, new_entry);
        g_hash_table_insert (index, (gpointer) fullname, pos);
      }

      e = &g_array_index (entries, NcmMSetPresetEntry, GPOINTER_TO_UINT (pos));

      if (a == NCM_MSET_PRESET_ATTR_FIT)
      {
        if (!g_variant_is_of_type (var, G_VARIANT_TYPE_BOOLEAN))
        {
          ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_INVALID_ENTRY,
                                      "fit of `%s' must be a boolean.", fullname);
          g_hash_table_unref (index);

          return;
        }

        e->fit = g_variant_get_boolean (var);
      }
      else
      {
        if (!g_variant_is_of_type (var, G_VARIANT_TYPE_DOUBLE))
        {
          ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_INVALID_ENTRY,
                                      "%s of `%s' must be a double.", _ncm_mset_preset_attr_names[a], fullname);
          g_hash_table_unref (index);

          return;
        }

        e->dval[a] = g_variant_get_double (var);
      }

      e->set[a] = TRUE;
    }
  }

  g_hash_table_unref (index);
}

/* Checks that every model with entries and not in @checked has no reparametrization. */
static void
_ncm_mset_preset_check_unnamed_reparam (GArray *entries, GHashTable *checked, GError **error)
{
  guint i;

  for (i = 0; i < entries->len; i++)
  {
    const NcmMSetPresetEntry *e = &g_array_index (entries, NcmMSetPresetEntry, i);
    NcmReparam *reparam         = ncm_model_peek_reparam (e->model);

    if ((reparam != NULL) && !g_hash_table_contains (checked, e->model))
    {
      ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_REPARAM_MISMATCH,
                                  "model `%s' has reparametrization `%s', the preset requires none.",
                                  G_OBJECT_TYPE_NAME (e->model), G_OBJECT_TYPE_NAME (reparam));

      return;
    }
  }
}

/*
 * Checks that each entry, with the preset applied, has lower bound < upper bound, its
 * value within the bounds, a positive scale and a non-negative absolute tolerance.
 */
static void
_ncm_mset_preset_check_entries (GArray *entries, GError **error)
{
  guint i;

  for (i = 0; i < entries->len; i++)
  {
    const NcmMSetPresetEntry *e = &g_array_index (entries, NcmMSetPresetEntry, i);
    const gdouble value         = e->dval[NCM_MSET_PRESET_ATTR_VALUE];
    const gdouble lb            = e->dval[NCM_MSET_PRESET_ATTR_LOWER_BOUND];
    const gdouble ub            = e->dval[NCM_MSET_PRESET_ATTR_UPPER_BOUND];
    const gdouble scale         = e->dval[NCM_MSET_PRESET_ATTR_SCALE];
    const gdouble abstol        = e->dval[NCM_MSET_PRESET_ATTR_ABSTOL];
    const gchar *pname          = ncm_model_param_name (e->model, e->pid);
    const gchar *mname          = G_OBJECT_TYPE_NAME (e->model);

    if (!(lb < ub))
    {
      ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_INVALID_ENTRY,
                                  "`%s' of `%s' has lower bound % 22.15g not below upper bound % 22.15g.",
                                  pname, mname, lb, ub);

      return;
    }

    if (!((value >= lb) && (value <= ub)))
    {
      ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_INVALID_ENTRY,
                                  "`%s' of `%s' has value % 22.15g outside [% 22.15g, % 22.15g].",
                                  pname, mname, value, lb, ub);

      return;
    }

    if (!(scale > 0.0))
    {
      ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_INVALID_ENTRY,
                                  "`%s' of `%s' has non-positive scale % 22.15g.",
                                  pname, mname, scale);

      return;
    }

    if (!(abstol >= 0.0))
    {
      ncm_util_set_or_call_error (error, NCM_MSET_PRESET_ERROR, NCM_MSET_PRESET_ERROR_INVALID_ENTRY,
                                  "`%s' of `%s' has negative absolute tolerance % 22.15g.",
                                  pname, mname, abstol);

      return;
    }
  }
}

