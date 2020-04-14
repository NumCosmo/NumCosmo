/* -*- Mode: C; indent-tabs-mode: nil; c-basic-offset: 2; tab-width: 2 -*-  */
/***************************************************************************
 *            nc_hipert_bg_var.c
 *
 *  Fri October 13 15:57:34 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * nc_hipert_bg_var.c
 * Copyright (C) 2017 Sandro Dias Pinto Vitenti <sandro@isoftware.com.br>
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
 * SECTION:nc_hipert_bg_var
 * @title: NcHIPertComp
 * @short_description: Perturbation background variables transport object
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "perturbations/nc_hipert_bg_var.h"

struct _NcHIPertBGVarPrivate
{
  gdouble zf;
};

enum
{
  PROP_0,
  PROP_DIST,
  PROP_RECOMB,
  PROP_SCALEFACTOR,
  PROP_ZF,
  PROP_LEN,
};

G_DEFINE_BOXED_TYPE (NcHIPertBGVarYDY, nc_hipert_bg_var_ydy, nc_hipert_bg_var_ydy_dup, nc_hipert_bg_var_ydy_free);
G_DEFINE_TYPE_WITH_PRIVATE (NcHIPertBGVar, nc_hipert_bg_var, G_TYPE_OBJECT);

static void
nc_hipert_bg_var_init (NcHIPertBGVar *bg_var)
{
  bg_var->priv = nc_hipert_bg_var_get_instance_private (bg_var);

  bg_var->cstructs = g_ptr_array_new ();

  bg_var->priv->zf = 0.0;
  
  g_ptr_array_set_size (bg_var->cstructs, nc_hipert_bg_var_len (bg_var));
  g_ptr_array_set_free_func (bg_var->cstructs, g_free);

  bg_var->recomb = NULL;
  bg_var->dist   = NULL;
  bg_var->a      = NULL;
  bg_var->t      = 0.0;
  bg_var->eta    = 0.0;
  bg_var->k      = 0.0;
  bg_var->x      = 0.0;
  bg_var->E      = 0.0;
}

static void
_nc_hipert_bg_var_set_property (GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec)
{
  NcHIPertBGVar *bg_var = NC_HIPERT_BG_VAR (object);
  g_return_if_fail (NC_IS_HIPERT_BG_VAR (object));

  switch (prop_id)
  {
    case PROP_DIST:
      nc_hipert_bg_var_set_dist (bg_var, g_value_get_object (value));
      break;
    case PROP_RECOMB:
      nc_hipert_bg_var_set_recomb (bg_var, g_value_get_object (value));
      break;
    case PROP_SCALEFACTOR:
      nc_hipert_bg_var_set_scalefactor (bg_var, g_value_get_object (value));
      break;
    case PROP_ZF:
      nc_hipert_bg_var_set_zf (bg_var, g_value_get_double (value));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_bg_var_get_property (GObject *object, guint prop_id, GValue *value, GParamSpec *pspec)
{
  NcHIPertBGVar *bg_var = NC_HIPERT_BG_VAR (object);
  g_return_if_fail (NC_IS_HIPERT_BG_VAR (object));

  switch (prop_id)
  {
    case PROP_DIST:
      g_value_take_object (value, nc_hipert_bg_var_get_dist (bg_var));
      break;
    case PROP_RECOMB:
      g_value_take_object (value, nc_hipert_bg_var_get_recomb (bg_var));
      break;
    case PROP_SCALEFACTOR:
      g_value_take_object (value, nc_hipert_bg_var_get_scalefactor (bg_var));
      break;
    case PROP_ZF:
      g_value_set_double (value, nc_hipert_bg_var_get_zf (bg_var));
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
_nc_hipert_bg_var_dispose (GObject *object)
{
  NcHIPertBGVar *bg_var = NC_HIPERT_BG_VAR (object);
  
  g_clear_pointer (&bg_var->cstructs, g_ptr_array_unref);

  nc_distance_clear (&bg_var->dist);
  nc_scalefactor_clear (&bg_var->a);
  nc_recomb_clear (&bg_var->recomb);

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_bg_var_parent_class)->dispose (object);
}

static void
_nc_hipert_bg_var_finalize (GObject *object)
{

  /* Chain up : end */
  G_OBJECT_CLASS (nc_hipert_bg_var_parent_class)->finalize (object);
}

static void
nc_hipert_bg_var_class_init (NcHIPertBGVarClass *klass)
{
  GObjectClass* object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = &_nc_hipert_bg_var_set_property;
  object_class->get_property = &_nc_hipert_bg_var_get_property;
  object_class->dispose      = &_nc_hipert_bg_var_dispose;
  object_class->finalize     = &_nc_hipert_bg_var_finalize;

  g_object_class_install_property (object_class,
                                   PROP_DIST,
                                   g_param_spec_object ("distance",
                                                        NULL,
                                                        "Distance object",
                                                        NC_TYPE_DISTANCE,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_RECOMB,
                                   g_param_spec_object ("recomb",
                                                        NULL,
                                                        "Recombination object",
                                                        NC_TYPE_RECOMB,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_SCALEFACTOR,
                                   g_param_spec_object ("scalefactor",
                                                        NULL,
                                                        "Scalefactor object",
                                                        NC_TYPE_SCALEFACTOR,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  g_object_class_install_property (object_class,
                                   PROP_ZF,
                                   g_param_spec_double ("zf",
                                                        NULL,
                                                        "Maximum redshift",
                                                        1.0, G_MAXDOUBLE, NC_HIPERT_BG_VAR_DEFAULT_ZF,
                                                        G_PARAM_READWRITE | G_PARAM_STATIC_NAME | G_PARAM_STATIC_BLURB));
  klass->bg_var_id_len     = 0;
  klass->ns_table          = g_hash_table_new (g_str_hash, g_str_equal);
  klass->bg_var_desc_array = g_array_new (FALSE, TRUE, sizeof (NcHIPertBGVarDesc));
}

G_LOCK_DEFINE_STATIC (last_bg_var_id);

/**
 * nc_hipert_bg_var_class_register_id: (skip)
 * @ns: object namespace
 * @desc: short description
 * @long_desc: long description
 * @cstruct_size: component struct size
 *
 * FIXME
 *
 */
void
nc_hipert_bg_var_class_register_id (const gchar *ns, const gchar *desc, const gchar *long_desc, guint cstruct_size)
{
  NcHIPertBGVarClass *bg_var_class = g_type_class_ref (NC_TYPE_HIPERT_BG_VAR);
  NcHIPertBGVarDesc *bg_var_desc   = NULL;
  NcHIPertBGVarID id;

  G_LOCK (last_bg_var_id);

  id = bg_var_class->bg_var_id_len;

  bg_var_class->bg_var_id_len++;
  
  g_array_set_size (bg_var_class->bg_var_desc_array, bg_var_class->bg_var_id_len);

  bg_var_desc       = &g_array_index (bg_var_class->bg_var_desc_array, NcHIPertBGVarDesc, id);
  bg_var_desc->init = TRUE;

  if (ns == NULL)
    g_error ("Cannot register background variables without a namespace.");
  if (desc == NULL)
    g_error ("Cannot register background variables without a description.");

  bg_var_desc->ns           = g_strdup (ns);
  bg_var_desc->desc         = g_strdup (desc);
  bg_var_desc->cstruct_size = cstruct_size;

  if (long_desc != NULL)
    bg_var_desc->long_desc = g_strdup (long_desc);
  else
    bg_var_desc->long_desc = NULL;

  if (g_hash_table_lookup (bg_var_class->ns_table, ns) != NULL)
    g_error ("Background variable namespace <%s> already registered.", ns);

  g_hash_table_insert (bg_var_class->ns_table, bg_var_desc->ns, GINT_TO_POINTER (id));

  G_UNLOCK (last_bg_var_id);
  g_type_class_unref (bg_var_class);
  
  return;
}

/**
 * nc_hipert_bg_var_class_get_id_by_gtype: (skip)
 * @gt: a #GType
 *
 * Gets the id associated with the GType @gt.
 * 
 * Returns: the id of @gt.
 */
NcHIPertBGVarID 
nc_hipert_bg_var_class_get_id_by_gtype (GType gt)
{
  NcHIPertBGVarClass *bg_var_class = g_type_class_ref (NC_TYPE_HIPERT_BG_VAR);
  gpointer id_ptr = NULL;
  gboolean has_it = FALSE;
  
  G_LOCK (last_bg_var_id);

  has_it = g_hash_table_lookup_extended (bg_var_class->ns_table, g_type_name (gt), NULL, &id_ptr);

  G_UNLOCK (last_bg_var_id);

  g_type_class_unref (bg_var_class);
  return has_it ? GPOINTER_TO_INT (id_ptr) : -1;
}

/**
 * nc_hipert_bg_var_class_get_id_by_ns: (skip)
 * @ns: an object namespace
 *
 * Gets the id associated with the namespace @ns.
 * 
 * Returns: the id of @ns.
 */
NcHIPertBGVarID 
nc_hipert_bg_var_class_get_id_by_ns (const gchar *ns)
{
  NcHIPertBGVarClass *bg_var_class = g_type_class_ref (NC_TYPE_HIPERT_BG_VAR);
  gpointer id_ptr = NULL;
  gboolean has_it = FALSE;
  
  G_LOCK (last_bg_var_id);

  has_it = g_hash_table_lookup_extended (bg_var_class->ns_table, ns, NULL, &id_ptr);

  G_UNLOCK (last_bg_var_id);

  g_type_class_unref (bg_var_class);
  return has_it ? GPOINTER_TO_INT (id_ptr) : -1;
}

/**
 * nc_hipert_bg_var_ydy_new:
 * 
 * Creates a new #NcHIPertBGVarYDY with all
 * entries set to zero.
 * 
 * Returns: (transfer full): a new #NcHIPertBGVarYDY.
 */
/**
 * nc_hipert_bg_var_ydy_dup:
 * @ydy: a #NcHIPertBGVarYDY
 * 
 * Duplicates @ydy.
 * 
 * Returns: (transfer full): a copy of @ydy.
 */
/**
 * nc_hipert_bg_var_ydy_free:
 * @ydy: a #NcHIPertBGVarYDY
 * 
 * Frees @ydy.
 * 
 */
/**
 * nc_hipert_bg_var_ydy_get_y_i:
 * @ydy: a #NcHIPertBGVarYDY
 * @i: variable index
 * 
 * Gets the @i-th variable.
 * 
 * Returns: the value of the @i-th variable.
 */
/**
 * nc_hipert_bg_var_ydy_set_dy_i:
 * @ydy: a #NcHIPertBGVarYDY
 * @i: variable index
 * @dy_i: the value of the @i-th variable derivative
 * 
 * Sets the @i-th variable derivative to @dy_i.
 */
/**
 * nc_hipert_bg_var_ydy_get_dy_i:
 * @ydy: a #NcHIPertBGVarYDY
 * @i: variable index
 * 
 * Gets the @i-th variable derivative.
 * 
 * Returns: the value of the @i-th variable derivative.
 */

/**
 * nc_hipert_bg_var_new:
 * 
 * Creates a new #NcHIPertBGVar.
 * 
 * Returns: (transfer full): the newly instantiated #NcHIPertBGVar.
 */
NcHIPertBGVar *
nc_hipert_bg_var_new (void)
{
  NcHIPertBGVar *bg_var = g_object_new (NC_TYPE_HIPERT_BG_VAR,
                                        NULL);

  return bg_var;
}

/**
 * nc_hipert_bg_var_ref:
 * @bg_var: a #NcHIPertBGVar
 *
 * Increases the reference count of @bg_var.
 *
 * Returns: (transfer full): @bg_var.
 */
NcHIPertBGVar *
nc_hipert_bg_var_ref (NcHIPertBGVar *bg_var)
{
  return g_object_ref (bg_var);
}

/**
 * nc_hipert_bg_var_free:
 * @bg_var: a #NcHIPertBGVar
 *
 * Decreases the reference count of @bg_var.
 *
 */
void 
nc_hipert_bg_var_free (NcHIPertBGVar *bg_var)
{
  g_object_unref (bg_var);
}

/**
 * nc_hipert_bg_var_clear:
 * @bg_var: a #NcHIPertBGVar
 *
 * Decreases the reference count of *@bg_var and sets the pointer *@bg_var to NULL.
 *
 */
void 
nc_hipert_bg_var_clear (NcHIPertBGVar **bg_var)
{
  g_clear_object (bg_var);
}

/**
 * nc_hipert_bg_var_prepare:
 * @bg_var: a #NcHIPertBGVar
 * @cosmo: a #NcHICosmo
 *
 * Prepares all computation objects inside @bg_var.
 *
 */
void 
nc_hipert_bg_var_prepare (NcHIPertBGVar *bg_var, NcHICosmo *cosmo)
{
  NcDistance *dist = nc_hipert_bg_var_peek_dist (bg_var);
  NcRecomb *recomb = nc_hipert_bg_var_peek_recomb (bg_var);
  NcScalefactor *a = nc_hipert_bg_var_peek_scalefactor (bg_var);
  
  if (dist != NULL)
    nc_distance_prepare (dist, cosmo);

  if (recomb != NULL)
    nc_recomb_prepare (recomb, cosmo);

  if (a != NULL)
    nc_scalefactor_prepare (a, cosmo);
}

/**
 * nc_hipert_bg_var_prepare_if_needed:
 * @bg_var: a #NcHIPertBGVar
 * @cosmo: a #NcHICosmo
 *
 * Prepares all computation objects inside @bg_var if necessary.
 *
 */
void 
nc_hipert_bg_var_prepare_if_needed (NcHIPertBGVar *bg_var, NcHICosmo *cosmo)
{
  NcDistance *dist = nc_hipert_bg_var_peek_dist (bg_var);
  NcRecomb *recomb = nc_hipert_bg_var_peek_recomb (bg_var);
  NcScalefactor *a = nc_hipert_bg_var_peek_scalefactor (bg_var);
  
  if (dist != NULL)
    nc_distance_prepare_if_needed (dist, cosmo);

  if (recomb != NULL)
    nc_recomb_prepare_if_needed (recomb, cosmo);

  if (a != NULL)
    nc_scalefactor_prepare_if_needed (a, cosmo);
}

/**
 * nc_hipert_bg_var_set_dist:
 * @bg_var: a #NcHIPertBGVar
 * @dist: a #NcDistance
 * 
 * Sets the #NcDistance object.
 * 
 */
/**
 * nc_hipert_bg_var_set_recomb:
 * @bg_var: a #NcHIPertBGVar
 * @recomb: a #NcRecomb
 * 
 * Sets the #NcRecomb object.
 * 
 */
/**
 * nc_hipert_bg_var_set_scalefactor:
 * @bg_var: a #NcHIPertBGVar
 * @a: a #NcScalefactor
 * 
 * Sets the #NcScalefactor object.
 * 
 */
/**
 * nc_hipert_bg_var_get_dist:
 * @bg_var: a #NcHIPertBGVar
 * 
 * Gets the #NcDistance object.
 * 
 * Returns: (transfer full) (nullable): the #NcDistance object used by @bg_var.
 */
/**
 * nc_hipert_bg_var_get_recomb:
 * @bg_var: a #NcHIPertBGVar
 * 
 * Gets the #NcRecomb object.
 * 
 * Returns: (transfer full) (nullable): the #NcRecomb object used by @bg_var.
 */
/**
 * nc_hipert_bg_var_get_scalefactor:
 * @bg_var: a #NcHIPertBGVar
 * 
 * Gets the #NcScalefactor object.
 * 
 * Returns: (transfer full) (nullable): the #NcScalefactor object used by @bg_var.
 */
/**
 * nc_hipert_bg_var_peek_dist:
 * @bg_var: a #NcHIPertBGVar
 * 
 * Peeks the #NcDistance object.
 * 
 * Returns: (transfer none) (nullable): the #NcDistance object used by @bg_var.
 */
/**
 * nc_hipert_bg_var_peek_recomb:
 * @bg_var: a #NcHIPertBGVar
 * 
 * Peeks the #NcRecomb object.
 * 
 * Returns: (transfer none) (nullable): the #NcRecomb object used by @bg_var.
 */
/**
 * nc_hipert_bg_var_peek_scalefactor:
 * @bg_var: a #NcHIPertBGVar
 * 
 * Peeks the #NcScalefactor object.
 * 
 * Returns: (transfer none) (nullable): the #NcScalefactor object used by @bg_var.
 */

/**
 * nc_hipert_bg_var_set_zf:
 * @bg_var: a #NcHIPertBGVar
 * @zf: the maximum redshift where calculations take place $z_f$
 * 
 * Requires the maximum redshift for computations to be @zf.
 * 
 */
void 
nc_hipert_bg_var_set_zf (NcHIPertBGVar *bg_var, const gdouble zf)
{
  g_assert_cmpfloat (zf, >, 0.0);
  if (zf != bg_var->priv->zf)
  {
    NcDistance *dist = nc_hipert_bg_var_peek_dist (bg_var);
    NcRecomb *recomb = nc_hipert_bg_var_peek_recomb (bg_var);
    NcScalefactor *a = nc_hipert_bg_var_peek_scalefactor (bg_var);

    if (dist != NULL)
      nc_distance_require_zf (dist, zf);

    if (recomb != NULL)
      nc_recomb_require_zi (recomb, zf);

    if (a != NULL)
      nc_scalefactor_require_zf (a, zf);

    bg_var->priv->zf = zf;
  }
}

/**
 * nc_hipert_bg_var_get_zf:
 * @bg_var: a #NcHIPertBGVar
 * 
 * Returns: the maximum redshift for computations @zf.
 */
gdouble 
nc_hipert_bg_var_get_zf (NcHIPertBGVar *bg_var)
{
  return bg_var->priv->zf;
}

/**
 * nc_hipert_bg_var_len:
 * @bg_var: a #NcHIPertBGVar
 *
 * Returns: the number of possible components installed.
 */
guint 
nc_hipert_bg_var_len (NcHIPertBGVar *bg_var)
{
  return NC_HIPERT_BG_VAR_GET_CLASS (bg_var)->bg_var_desc_array->len;
}

/**
 * nc_hipert_bg_var_cstruct_len:
 * @bg_var: a #NcHIPertBGVar
 * @id: a #NcHIPertBGVarID
 *
 * Returns: the number of possible components installed.
 */
guint 
nc_hipert_bg_var_cstruct_len (NcHIPertBGVar *bg_var, NcHIPertBGVarID id)
{
  const guint len = nc_hipert_bg_var_len (bg_var);
  g_assert_cmpuint (id, <, len);

  {
    NcHIPertBGVarDesc *desc = &g_array_index (NC_HIPERT_BG_VAR_GET_CLASS (bg_var)->bg_var_desc_array,
                                              NcHIPertBGVarDesc,
                                              id);
    return desc->cstruct_size;
  }
}

static void
_nc_hipert_bg_var_activate_id (NcHIPertBGVar *bg_var, NcHIPertBGVarID id)
{
  const guint len = nc_hipert_bg_var_len (bg_var);

  g_assert_cmpuint (id, <, len);

  g_ptr_array_set_size (bg_var->cstructs, len);

  if (g_ptr_array_index (bg_var->cstructs, id) != NULL)
    g_warning ("_nc_hipert_bg_var_activate_id: component %d already activated, ignoring...", id);
  else
    g_ptr_array_index (bg_var->cstructs, id) = g_malloc0 (nc_hipert_bg_var_cstruct_len (bg_var, id));
}

/**
 * nc_hipert_bg_var_activate_id: (skip)
 * @bg_var: a #NcHIPertBGVar
 * @...: a list of background ids
 *
 * Register which components that will use the object @bg_var.
 * the list ... must end with a -1 signaling the end of the list.
 *
 */
void
nc_hipert_bg_var_activate_id (NcHIPertBGVar *bg_var, ...)
{
  NcHIPertBGVarID id; 
  va_list ap;

  va_start (ap, bg_var);

  while ((id = va_arg (ap, NcHIPertBGVarID)) >= 0)
  {
    _nc_hipert_bg_var_activate_id (bg_var, id);
  }
  
  va_end (ap);
}

/**
 * nc_hipert_bg_var_activate_id_array:
 * @bg_var: a #NcHIPertBGVar
 * @ids: (array) (element-type NcHIPertBGVarID): an array of background ids
 *
 * Register which components that will use the object @bg_var.
 * the list ... must end with a -1 signaling the end of the list.
 *
 */
void
nc_hipert_bg_var_activate_id_array (NcHIPertBGVar *bg_var, GArray *ids)
{
  NcHIPertBGVarID id;
  guint i;

  for (i = 0; i < ids->len; i++)
  {
    id = g_array_index (ids, NcHIPertBGVarID, i);
    _nc_hipert_bg_var_activate_id (bg_var, id);
  }
}
