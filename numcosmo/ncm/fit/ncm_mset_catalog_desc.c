/***************************************************************************
 *            ncm_mset_catalog_desc.c
 *
 *  Thu Oct 29 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * ncm_mset_catalog_desc.c
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
 * NcmMSetCatalogDesc:
 *
 * Immutable description of the columns and metadata of an #NcmMSetCatalog.
 *
 * The description includes the number of chains, whether rows are weighted,
 * the additional value representing $-2\ln(p)$ and the run type.
 *
 * Each row contains the additional values followed by the free parameters
 * of the #NcmMSet, ordered according to its free-parameter map. The additional
 * values are those provided at construction, followed by the weight column
 * when the catalog is weighted.
 *
 * All properties are set at construction and cannot be modified afterward.
 * The description stores a copy of the #NcmMSet, ensuring that subsequent
 * changes to the original do not affect it.
 *
 * The properties #NcmMSetCatalogDesc:nadd-val-names and
 * #NcmMSetCatalogDesc:nadd-val-symbols exclude the weight column, whereas
 * ncm_mset_catalog_desc_nadd_vals() and the column accessors include it.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/fit/ncm_mset_catalog_desc.h"
#include "ncm/fit/ncm_mset_catalog.h"
#include "ncm/core/ncm_serialize.h"
#include "ncm/core/ncm_util.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <stdlib.h>
#include <string.h>
#endif /* NUMCOSMO_GIR_SCAN */

/* *INDENT-OFF* */
G_DEFINE_QUARK (ncm-mset-catalog-desc-error, ncm_mset_catalog_desc_error)
/* *INDENT-ON* */

struct _NcmMSetCatalogDesc
{
  GObject parent_instance;
  NcmMSet *mset;
  guint nchains;
  gboolean weighted;
  gint m2lnp_var;
  gchar *run_type;
  /* The user's additional values, then the weight column when weighted. */
  GPtrArray *add_val_names;
  GPtrArray *add_val_symbols;
};

enum
{
  PROP_0,
  PROP_MSET,
  PROP_NCHAINS,
  PROP_WEIGHTED,
  PROP_ADD_VAL_NAMES,
  PROP_ADD_VAL_SYMBS,
  PROP_M2LNP_VAR,
  PROP_RUN_TYPE_STR,
};

G_DEFINE_TYPE (NcmMSetCatalogDesc, ncm_mset_catalog_desc, G_TYPE_OBJECT)

static void
ncm_mset_catalog_desc_init (NcmMSetCatalogDesc *desc)
{
  desc->mset            = NULL;
  desc->nchains         = 0;
  desc->weighted        = FALSE;
  desc->m2lnp_var       = 0;
  desc->run_type        = NULL;
  desc->add_val_names   = g_ptr_array_new_with_free_func (g_free);
  desc->add_val_symbols = g_ptr_array_new_with_free_func (g_free);
}

static void
_ncm_mset_catalog_desc_set_strv (GPtrArray *array, const gchar * const *strv)
{
  g_ptr_array_set_size (array, 0);

  if (strv != NULL)
  {
    guint i;

    for (i = 0; strv[i] != NULL; i++)
      g_ptr_array_add (array, g_strdup (strv[i]));
  }
}

static gchar **
_ncm_mset_catalog_desc_dup_user_strv (NcmMSetCatalogDesc *desc, GPtrArray *array)
{
  const guint len = array->len - (desc->weighted ? 1 : 0);
  gchar **strv    = g_new0 (gchar *, len + 1);
  guint i;

  for (i = 0; i < len; i++)
    strv[i] = g_strdup (g_ptr_array_index (array, i));

  return strv;
}

static void
_ncm_mset_catalog_desc_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalogDesc *desc = NCM_MSET_CATALOG_DESC (object);

  g_return_if_fail (NCM_IS_MSET_CATALOG_DESC (object));

  switch (prop_id)
  {
    case PROP_MSET:
    {
      NcmMSet *mset = g_value_get_object (value);

      if (mset != NULL)
      {
        NcmSerialize *ser = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);

        desc->mset = ncm_mset_dup (mset, ser);
        ncm_serialize_free (ser);
      }

      break;
    }
    case PROP_NCHAINS:
      desc->nchains = g_value_get_uint (value);
      break;
    case PROP_WEIGHTED:
      desc->weighted = g_value_get_boolean (value);
      break;
    case PROP_ADD_VAL_NAMES:
      _ncm_mset_catalog_desc_set_strv (desc->add_val_names, g_value_get_boxed (value));
      break;
    case PROP_ADD_VAL_SYMBS:
      _ncm_mset_catalog_desc_set_strv (desc->add_val_symbols, g_value_get_boxed (value));
      break;
    case PROP_M2LNP_VAR:
      desc->m2lnp_var = g_value_get_int (value);
      break;
    case PROP_RUN_TYPE_STR:
      desc->run_type = g_value_dup_string (value);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static void
_ncm_mset_catalog_desc_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcmMSetCatalogDesc *desc = NCM_MSET_CATALOG_DESC (object);

  g_return_if_fail (NCM_IS_MSET_CATALOG_DESC (object));

  switch (prop_id)
  {
    case PROP_MSET:
      g_value_set_object (value, desc->mset);
      break;
    case PROP_NCHAINS:
      g_value_set_uint (value, desc->nchains);
      break;
    case PROP_WEIGHTED:
      g_value_set_boolean (value, desc->weighted);
      break;
    case PROP_ADD_VAL_NAMES:
      g_value_take_boxed (value, _ncm_mset_catalog_desc_dup_user_strv (desc, desc->add_val_names));
      break;
    case PROP_ADD_VAL_SYMBS:
      g_value_take_boxed (value, _ncm_mset_catalog_desc_dup_user_strv (desc, desc->add_val_symbols));
      break;
    case PROP_M2LNP_VAR:
      g_value_set_int (value, desc->m2lnp_var);
      break;
    case PROP_RUN_TYPE_STR:
      g_value_set_string (value, desc->run_type);
      break;
    default:                                                      /* LCOV_EXCL_LINE */
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec); /* LCOV_EXCL_LINE */
      break;                                                      /* LCOV_EXCL_LINE */
  }
}

static gboolean _ncm_mset_catalog_desc_check_input (NcmMSet *mset, guint names_len, guint symbols_len, gint m2lnp_var, GError **error);

static void
_ncm_mset_catalog_desc_constructed (GObject *object)
{
  /* Chain up : start */
  G_OBJECT_CLASS (ncm_mset_catalog_desc_parent_class)->constructed (object);
  {
    NcmMSetCatalogDesc *desc = NCM_MSET_CATALOG_DESC (object);

    /* Built from properties there is no error location: an invalid input aborts. */
    _ncm_mset_catalog_desc_check_input (desc->mset,
                                        desc->add_val_names->len,
                                        desc->add_val_symbols->len,
                                        desc->m2lnp_var,
                                        NULL);

    if (desc->weighted)
    {
      g_ptr_array_add (desc->add_val_names, g_strdup (NCM_MSET_CATALOG_DESC_WEIGHT_COLNAME));
      g_ptr_array_add (desc->add_val_symbols, g_strdup (NCM_MSET_CATALOG_DESC_WEIGHT_SYMBOL));
    }
  }
}

static void
_ncm_mset_catalog_desc_dispose (GObject *object)
{
  NcmMSetCatalogDesc *desc = NCM_MSET_CATALOG_DESC (object);

  ncm_mset_clear (&desc->mset);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_desc_parent_class)->dispose (object);
}

static void
_ncm_mset_catalog_desc_finalize (GObject *object)
{
  NcmMSetCatalogDesc *desc = NCM_MSET_CATALOG_DESC (object);

  g_clear_pointer (&desc->run_type, g_free);
  g_clear_pointer (&desc->add_val_names, g_ptr_array_unref);
  g_clear_pointer (&desc->add_val_symbols, g_ptr_array_unref);

  /* Chain up : end */
  G_OBJECT_CLASS (ncm_mset_catalog_desc_parent_class)->finalize (object);
}

static void
ncm_mset_catalog_desc_class_init (NcmMSetCatalogDescClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_ncm_mset_catalog_desc_set_property;
  object_class->get_property = &_ncm_mset_catalog_desc_get_property;
  object_class->constructed  = &_ncm_mset_catalog_desc_constructed;
  object_class->dispose      = &_ncm_mset_catalog_desc_dispose;
  object_class->finalize     = &_ncm_mset_catalog_desc_finalize;

  /**
   * NcmMSetCatalogDesc:mset:
   *
   * The #NcmMSet whose free parameters are the parameter columns. The description holds a
   * copy, made when the property is set; its free-parameter map must be prepared.
   */
  g_object_class_install_property (object_class,
                                   PROP_MSET,
                                   g_param_spec_object ("mset",
                                                        NULL,
                                                        "NcmMSet object",
                                                        NCM_TYPE_MSET,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_NCHAINS,
                                   g_param_spec_uint ("nchains",
                                                      NULL,
                                                      "Number of different chains in the catalog",
                                                      1, G_MAXUINT32, 1,
                                                      G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_WEIGHTED,
                                   g_param_spec_boolean ("weighted",
                                                         NULL,
                                                         "Catalog with weighted rows",
                                                         FALSE,
                                                         G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogDesc:nadd-val-names:
   *
   * The names of the user's additional values, without the weight column.
   */
  g_object_class_install_property (object_class,
                                   PROP_ADD_VAL_NAMES,
                                   g_param_spec_boxed ("nadd-val-names",
                                                       NULL,
                                                       "Additional value names",
                                                       G_TYPE_STRV,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogDesc:nadd-val-symbols:
   *
   * The symbols of the user's additional values, without the weight column.
   */
  g_object_class_install_property (object_class,
                                   PROP_ADD_VAL_SYMBS,
                                   g_param_spec_boxed ("nadd-val-symbols",
                                                       NULL,
                                                       "Additional value symbols",
                                                       G_TYPE_STRV,
                                                       G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));

  /**
   * NcmMSetCatalogDesc:m2lnp-var:
   *
   * The index of the additional value holding $-2\ln(p)$, or -1 when there is none.
   */
  g_object_class_install_property (object_class,
                                   PROP_M2LNP_VAR,
                                   g_param_spec_int ("m2lnp-var",
                                                     NULL,
                                                     "Index of the variable representing m2lnp",
                                                     -1, G_MAXINT, -1,
                                                     G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RUN_TYPE_STR,
                                   g_param_spec_string ("run-type-string",
                                                        NULL,
                                                        "Run type string",
                                                        NCM_MSET_CATALOG_RTYPE_UNDEFINED,
                                                        G_PARAM_READWRITE | G_PARAM_CONSTRUCT_ONLY | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
}

/* Checks the inputs that the parameter specifications cannot: the mset is given and its
 * free-parameter map is prepared, names and symbols agree in number, and m2lnp_var is -1 or
 * one of the user's additional values. */
static gboolean
_ncm_mset_catalog_desc_check_input (NcmMSet *mset, guint names_len, guint symbols_len, gint m2lnp_var, GError **error)
{
  if (mset == NULL)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INVALID_MSET,
                                "ncm_mset_catalog_desc: no NcmMSet given.");

    return FALSE;
  }

  if (!ncm_mset_fparam_map_valid (mset))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INVALID_MSET,
                                "ncm_mset_catalog_desc: the free-parameter map of the NcmMSet is not prepared.");

    return FALSE;
  }

  if (names_len != symbols_len)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INVALID_ADD_VALS,
                                "ncm_mset_catalog_desc: %u additional-value names and %u symbols.",
                                names_len, symbols_len);

    return FALSE;
  }

  if ((m2lnp_var < -1) || (m2lnp_var >= (gint) names_len))
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INVALID_M2LNP_VAR,
                                "ncm_mset_catalog_desc: m2lnp-var %d is not -1 or one of the %u additional values.",
                                m2lnp_var, names_len);

    return FALSE;
  }

  return TRUE;
}

/**
 * ncm_mset_catalog_desc_new:
 * @mset: a #NcmMSet
 * @nchains: number of chains
 * @weighted: whether the rows carry a weight column
 * @names: (in) (array zero-terminated=1) (element-type utf8) (nullable): names of the user's additional values
 * @symbols: (in) (array zero-terminated=1) (element-type utf8) (nullable): symbols of the user's additional values
 * @m2lnp_var: index of the additional value holding $-2\ln(p)$, or -1
 * @run_type: (nullable): run type string
 * @error: a #GError
 *
 * Creates a description from its inputs; see #NcmMSetCatalogDesc for the column layout. A
 * %NULL @run_type gives %NCM_MSET_CATALOG_RTYPE_UNDEFINED. The free-parameter map of @mset
 * must be prepared. On an invalid input it sets @error and returns %NULL.
 *
 * Returns: (transfer full) (nullable): a new #NcmMSetCatalogDesc
 */
NcmMSetCatalogDesc *
ncm_mset_catalog_desc_new (NcmMSet *mset, guint nchains, gboolean weighted, const gchar * const *names, const gchar * const *symbols, gint m2lnp_var, const gchar *run_type, GError **error)
{
  const guint names_len   = (names != NULL) ? g_strv_length ((gchar **) names) : 0;
  const guint symbols_len = (symbols != NULL) ? g_strv_length ((gchar **) symbols) : 0;

  if (!_ncm_mset_catalog_desc_check_input (mset, names_len, symbols_len, m2lnp_var, error))
    return NULL;

  return g_object_new (NCM_TYPE_MSET_CATALOG_DESC,
                       "mset",             mset,
                       "nchains",          nchains,
                       "weighted",         weighted,
                       "nadd-val-names",   names,
                       "nadd-val-symbols", symbols,
                       "m2lnp-var",        m2lnp_var,
                       "run-type-string",  (run_type != NULL) ? run_type : NCM_MSET_CATALOG_RTYPE_UNDEFINED,
                       NULL);
}

/**
 * ncm_mset_catalog_desc_ref:
 * @desc: a #NcmMSetCatalogDesc
 *
 * Increases the reference count of @desc by one.
 *
 * Returns: (transfer full): @desc
 */
NcmMSetCatalogDesc *
ncm_mset_catalog_desc_ref (NcmMSetCatalogDesc *desc)
{
  return g_object_ref (desc);
}

/**
 * ncm_mset_catalog_desc_free:
 * @desc: a #NcmMSetCatalogDesc
 *
 * Decreases the reference count of @desc by one.
 */
void
ncm_mset_catalog_desc_free (NcmMSetCatalogDesc *desc)
{
  g_object_unref (desc);
}

/**
 * ncm_mset_catalog_desc_clear:
 * @desc: a #NcmMSetCatalogDesc
 *
 * If *@desc is not %NULL, decreases its reference count and sets *@desc to %NULL.
 */
void
ncm_mset_catalog_desc_clear (NcmMSetCatalogDesc **desc)
{
  g_clear_object (desc);
}

/**
 * ncm_mset_catalog_desc_peek_mset:
 * @desc: a #NcmMSetCatalogDesc
 *
 * Returns the copy of the #NcmMSet held by @desc, the same object at every call.
 *
 * Returns: (transfer none): the #NcmMSet of @desc
 */
NcmMSet *
ncm_mset_catalog_desc_peek_mset (NcmMSetCatalogDesc *desc)
{
  return desc->mset;
}

/**
 * ncm_mset_catalog_desc_nchains:
 * @desc: a #NcmMSetCatalogDesc
 *
 * Returns: the number of chains
 */
guint
ncm_mset_catalog_desc_nchains (NcmMSetCatalogDesc *desc)
{
  return desc->nchains;
}

/**
 * ncm_mset_catalog_desc_weighted:
 * @desc: a #NcmMSetCatalogDesc
 *
 * Returns: whether the rows carry a weight column
 */
gboolean
ncm_mset_catalog_desc_weighted (NcmMSetCatalogDesc *desc)
{
  return desc->weighted;
}

/**
 * ncm_mset_catalog_desc_get_m2lnp_var:
 * @desc: a #NcmMSetCatalogDesc
 *
 * Returns: the index of the additional value holding $-2\ln(p)$, or -1
 */
gint
ncm_mset_catalog_desc_get_m2lnp_var (NcmMSetCatalogDesc *desc)
{
  return desc->m2lnp_var;
}

/**
 * ncm_mset_catalog_desc_get_run_type:
 * @desc: a #NcmMSetCatalogDesc
 *
 * Returns: (transfer none): the run type string
 */
const gchar *
ncm_mset_catalog_desc_get_run_type (NcmMSetCatalogDesc *desc)
{
  return desc->run_type;
}

/**
 * ncm_mset_catalog_desc_nadd_vals:
 * @desc: a #NcmMSetCatalogDesc
 *
 * Returns the number of additional values, the weight column included. It is the index of
 * the first free-parameter column, with the same meaning as ncm_mset_catalog_nadd_vals().
 *
 * Returns: the number of additional values
 */
guint
ncm_mset_catalog_desc_nadd_vals (NcmMSetCatalogDesc *desc)
{
  return desc->add_val_names->len;
}

/**
 * ncm_mset_catalog_desc_ncols:
 * @desc: a #NcmMSetCatalogDesc
 *
 * Returns: the number of columns, additional values and free parameters
 */
guint
ncm_mset_catalog_desc_ncols (NcmMSetCatalogDesc *desc)
{
  return desc->add_val_names->len + ncm_mset_fparams_len (desc->mset);
}

/**
 * ncm_mset_catalog_desc_col_name:
 * @desc: a #NcmMSetCatalogDesc
 * @i: column index
 *
 * Returns: (transfer none): the name of column @i
 */
const gchar *
ncm_mset_catalog_desc_col_name (NcmMSetCatalogDesc *desc, guint i)
{
  if (i < desc->add_val_names->len)
    return g_ptr_array_index (desc->add_val_names, i);
  else
    return ncm_mset_fparam_name (desc->mset, i - desc->add_val_names->len);
}

/**
 * ncm_mset_catalog_desc_col_full_name:
 * @desc: a #NcmMSetCatalogDesc
 * @i: column index
 *
 * Returns: (transfer none): the full name of column @i, model-qualified for a free parameter
 */
const gchar *
ncm_mset_catalog_desc_col_full_name (NcmMSetCatalogDesc *desc, guint i)
{
  if (i < desc->add_val_names->len)
    return g_ptr_array_index (desc->add_val_names, i);
  else
    return ncm_mset_fparam_full_name (desc->mset, i - desc->add_val_names->len);
}

/**
 * ncm_mset_catalog_desc_col_symb:
 * @desc: a #NcmMSetCatalogDesc
 * @i: column index
 *
 * Returns: (transfer none): the symbol of column @i
 */
const gchar *
ncm_mset_catalog_desc_col_symb (NcmMSetCatalogDesc *desc, guint i)
{
  if (i < desc->add_val_symbols->len)
    return g_ptr_array_index (desc->add_val_symbols, i);
  else
    return ncm_mset_fparam_symbol (desc->mset, i - desc->add_val_symbols->len);
}

/**
 * ncm_mset_catalog_desc_col_by_name:
 * @desc: a #NcmMSetCatalogDesc
 * @name: a free-parameter name, an additional-value name or a column index
 * @col_index: (out): the column index
 *
 * Finds the column of @name, tried in order as a free-parameter name, an additional-value
 * name, then a decimal column index; the same lookup as ncm_mset_catalog_col_by_name().
 *
 * Returns: whether @name was found
 */
gboolean
ncm_mset_catalog_desc_col_by_name (NcmMSetCatalogDesc *desc, const gchar *name, guint *col_index)
{
  const NcmMSetPIndex *pi = ncm_mset_fparam_get_pi_by_name (desc->mset, name, NULL);

  if (pi != NULL)
  {
    col_index[0] = ncm_mset_fparam_get_fpi (desc->mset, pi->mid, pi->pid) + desc->add_val_names->len;

    return TRUE;
  }
  else if (g_ptr_array_find_with_equal_func (desc->add_val_names, name, g_str_equal, col_index))
  {
    return TRUE;
  }
  else
  {
    gchar *end_ptr  = NULL;
    glong add_param = strtol (name, &end_ptr, 10);

    if (name != end_ptr)
    {
      col_index[0] = add_param;

      return TRUE;
    }

    return FALSE;
  }
}

static gboolean _ncm_mset_catalog_desc_check_add_vals (NcmMSetCatalogDesc *desc, NcmMSetCatalogDesc *other, GError **error);
static gboolean _ncm_mset_catalog_desc_check_fparam_set (NcmMSetCatalogDesc *desc, NcmMSetCatalogDesc *other, GError **error);

/**
 * ncm_mset_catalog_desc_check_compatible:
 * @desc: a #NcmMSetCatalogDesc
 * @other: a #NcmMSetCatalogDesc
 * @error: a #GError
 *
 * Checks whether rows described by @other can be added to rows described by @desc: the
 * same number of chains, the same weighting, the same m2lnp-var and run type, the same
 * additional values (names and symbols, in order) and the same set of free-parameter full
 * names, in any order. The #NcmMSet models and parameter values are not compared. When they are not compatible it sets @error with
 * %NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE, naming the first difference.
 *
 * Returns: whether @desc and @other are compatible
 */
gboolean
ncm_mset_catalog_desc_check_compatible (NcmMSetCatalogDesc *desc, NcmMSetCatalogDesc *other, GError **error)
{
  if (desc->nchains != other->nchains)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE,
                                "ncm_mset_catalog_desc_check_compatible: %u chains and %u chains.",
                                desc->nchains, other->nchains);

    return FALSE;
  }

  if (desc->weighted != other->weighted)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE,
                                "ncm_mset_catalog_desc_check_compatible: one description is weighted and the other is not.");

    return FALSE;
  }

  if (desc->m2lnp_var != other->m2lnp_var)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE,
                                "ncm_mset_catalog_desc_check_compatible: m2lnp-var %d and %d.",
                                desc->m2lnp_var, other->m2lnp_var);

    return FALSE;
  }

  if (strcmp (desc->run_type, other->run_type) != 0)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE,
                                "ncm_mset_catalog_desc_check_compatible: run types `%s' and `%s'.",
                                desc->run_type, other->run_type);

    return FALSE;
  }

  if (!_ncm_mset_catalog_desc_check_add_vals (desc, other, error))
    return FALSE;

  return _ncm_mset_catalog_desc_check_fparam_set (desc, other, error);
}

/* Compares the additional values of desc and other, names and symbols, position by position. */
static gboolean
_ncm_mset_catalog_desc_check_add_vals (NcmMSetCatalogDesc *desc, NcmMSetCatalogDesc *other, GError **error)
{
  guint i;

  if (desc->add_val_names->len != other->add_val_names->len)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE,
                                "ncm_mset_catalog_desc_check_compatible: %u additional values and %u additional values.",
                                desc->add_val_names->len, other->add_val_names->len);

    return FALSE;
  }

  for (i = 0; i < desc->add_val_names->len; i++)
  {
    const gchar *name       = g_ptr_array_index (desc->add_val_names, i);
    const gchar *other_name = g_ptr_array_index (other->add_val_names, i);
    const gchar *symb       = g_ptr_array_index (desc->add_val_symbols, i);
    const gchar *other_symb = g_ptr_array_index (other->add_val_symbols, i);

    if ((strcmp (name, other_name) != 0) || (strcmp (symb, other_symb) != 0))
    {
      ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE,
                                  "ncm_mset_catalog_desc_check_compatible: additional value %u is `%s' (`%s') and `%s' (`%s').",
                                  i, name, symb, other_name, other_symb);

      return FALSE;
    }
  }

  return TRUE;
}

/* Compares the free-parameter full names of desc and other as sets. */
static gboolean
_ncm_mset_catalog_desc_check_fparam_set (NcmMSetCatalogDesc *desc, NcmMSetCatalogDesc *other, GError **error)
{
  const guint len       = ncm_mset_fparams_len (desc->mset);
  const guint other_len = ncm_mset_fparams_len (other->mset);
  GHashTable *names;
  gboolean compatible = TRUE;
  guint i;

  if (len != other_len)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE,
                                "ncm_mset_catalog_desc_check_compatible: %u free parameters and %u free parameters.",
                                len, other_len);

    return FALSE;
  }

  names = g_hash_table_new (g_str_hash, g_str_equal);

  for (i = 0; i < other_len; i++)
    g_hash_table_add (names, (gpointer) ncm_mset_fparam_full_name (other->mset, i));

  for (i = 0; i < len; i++)
  {
    const gchar *name = ncm_mset_fparam_full_name (desc->mset, i);

    if (!g_hash_table_contains (names, name))
    {
      ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_INCOMPATIBLE,
                                  "ncm_mset_catalog_desc_check_compatible: free parameter `%s' is not free in the other description.",
                                  name);
      compatible = FALSE;
      break;
    }
  }

  g_hash_table_unref (names);

  return compatible;
}

/**
 * ncm_mset_catalog_desc_check_same_order:
 * @desc: a #NcmMSetCatalogDesc
 * @other: a #NcmMSetCatalogDesc
 * @error: a #GError
 *
 * Checks whether @desc and @other list the same free parameters in the same order, compared
 * by full name. A row copied by position from a catalog described by @other into one
 * described by @desc keeps each parameter in its column only when this holds. Otherwise it
 * sets @error with %NCM_MSET_CATALOG_DESC_ERROR_ORDER, naming the first position that
 * differs.
 *
 * Returns: whether the free parameters are the same and in the same order
 */
gboolean
ncm_mset_catalog_desc_check_same_order (NcmMSetCatalogDesc *desc, NcmMSetCatalogDesc *other, GError **error)
{
  const guint len       = ncm_mset_fparams_len (desc->mset);
  const guint other_len = ncm_mset_fparams_len (other->mset);
  guint i;

  if (len != other_len)
  {
    ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_ORDER,
                                "ncm_mset_catalog_desc_check_same_order: %u free parameters and %u free parameters.",
                                len, other_len);

    return FALSE;
  }

  for (i = 0; i < len; i++)
  {
    const gchar *name       = ncm_mset_fparam_full_name (desc->mset, i);
    const gchar *other_name = ncm_mset_fparam_full_name (other->mset, i);

    if (strcmp (name, other_name) != 0)
    {
      ncm_util_set_or_call_error (error, NCM_MSET_CATALOG_DESC_ERROR, NCM_MSET_CATALOG_DESC_ERROR_ORDER,
                                  "ncm_mset_catalog_desc_check_same_order: free parameter %u is `%s' and `%s'.",
                                  i, name, other_name);

      return FALSE;
    }
  }

  return TRUE;
}

